# Produces: Figure 6 -- number of strains detected per prevalent species by dRep (salmon),
# inStrain (green) and StrainPhlAn4 (blue), faceted by kefir type.
#
# Each bar counts all datasets: this study's samples plus the published comparison
# datasets (Walsh et al for milk kefir; Breselge et al, recorded as "Samuel et al" or
# "Mortensen et al 2023", for water kefir). There is one bar per species x method.
#
# Species names are harmonised across methods via the prevalence crosswalk
# (data/milk_taxonomic_profile_prevalence.csv, water_taxonomic_profile_prevalence.csv), and
# each method's output is restricted to species prevalent in that kefir type.

library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style()
DATA_DIR    <- here::here("data")
TREES_DIR   <- here::here("data", "06_trees")
OUT_DIR     <- here::here("output_reviewed", "06_strain_profiling")
FIGURES_DIR <- here::here("figures", "Figure_6")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

library(readr); library(dplyr); library(tidyr); library(ape); library(ggplot2)

## ---- Prevalent-species crosswalk ----
# WHY: the prevalence tables split some species into StrainPhlAn-only "subcluster" /
# "subspecies" entries that inStrain and dRep report under one name, and use pre-2020
# Lactobacillus genus names. Relabelling puts all three methods on the same names.
merge_subclusters <- function(x) {
  x <- gsub("Lactococcus_lactis subcluster 1", "Lactococcus_lactis", x, fixed = TRUE)
  x <- gsub("Lactococcus_lactis subcluster 2", "Lactococcus_cremoris", x, fixed = TRUE)
  x <- gsub("Pseudomonas_fragi_subspecies 1", "Pseudomonas_fragi", x, fixed = TRUE)
  x <- gsub("Zymomonas_mobilis_subcluster 1", "Zymomonas_mobilis", x, fixed = TRUE)
  x <- gsub("Lactobacillus_ghanensis", "Liquorilactobacillus_ghanensis", x, fixed = TRUE)
  x <- gsub("Lactobacillus_satsumensis", "Liquorilactobacillus_satsumensis", x, fixed = TRUE)
  x <- gsub("Lactobacillus_perolens", "Schleiferilactobacillus_perolens", x, fixed = TRUE)
  x
}
milk_prevalence  <- read_csv(file.path(DATA_DIR, "milk_taxonomic_profile_prevalence.csv"),  show_col_types = FALSE)
water_prevalence <- read_csv(file.path(DATA_DIR, "water_taxonomic_profile_prevalence.csv"), show_col_types = FALSE)
species_col <- setdiff(names(milk_prevalence), c("...1", "specie", "sum_specie", "kefir_type"))[1]
stopifnot("expected prevalence columns 'specie' and a species-name column not found" =
            "specie" %in% names(milk_prevalence) && !is.na(species_col))
milk_prevalence[[species_col]]  <- merge_subclusters(milk_prevalence[[species_col]])
water_prevalence[[species_col]] <- merge_subclusters(water_prevalence[[species_col]])
clade_to_species <- bind_rows(milk_prevalence, water_prevalence) |>
  transmute(clade = specie, species = gsub("_", " ", .data[[species_col]])) |>
  distinct(clade, .keep_all = TRUE)
prevalent_milk  <- unique(gsub("_", " ", milk_prevalence[[species_col]]))
prevalent_water <- unique(gsub("_", " ", water_prevalence[[species_col]]))

## ---- Sample metadata: this study + published datasets ----
kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"), show_col_types = FALSE)
kefir4all_metadata$merge_column <- gsub("_host_removed_R..fastq.gz", "", kefir4all_metadata$merge_column)
kefir4all_metadata <- kefir4all_metadata[!duplicated(kefir4all_metadata$merge_column), ]
global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv"), show_col_types = FALSE)
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv"), show_col_types = FALSE)
total_metadata <- bind_rows(
  dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
  dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`),
  dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`)) |>
  mutate(category = case_when(
    `kefir type` %in% c("ML", "MG") ~ "Milk.kefir",
    `kefir type` %in% c("WL", "WG") ~ "Water.kefir",
    TRUE ~ NA_character_))

# Maps each dataset label to "Kefir4All" or "Published"; anything else (controls, unknown
# sources) becomes NA and is dropped.
dataset_group <- function(data_source) {
  case_when(
    data_source == "This study" ~ "Kefir4All",
    data_source %in% c("Walsh et al", "Walsh et al 2023",
                       "Mortensen et al 2023", "Samuel et al") ~ "Published",
    TRUE ~ NA_character_)
}
kefir_type_label <- c(Milk.kefir = "Milk kefir", Water.kefir = "Water kefir")

## ---- 1. StrainPhlAn4: samples in each clade's tree ----
tree_files <- list.files(TREES_DIR, pattern = "\\.tre$", full.names = TRUE)
if (length(tree_files) == 0) stop("No StrainPhlAn4 tree files found in ", TREES_DIR)
strainphlan_counts <- lapply(tree_files, function(f) {
  clade <- sub("^RAxML_bestTree\\.(t__.+)\\.StrainPhlAn4\\.tre$", "\\1", basename(f))
  species <- clade_to_species$species[clade_to_species$clade == clade]
  if (!length(species)) return(NULL)
  tips <- gsub("-", "_", ape::read.tree(f)$tip.label)
  total_metadata |>
    filter(merge_column %in% tips, !is.na(category)) |>
    transmute(species = species[1], category, dataset = dataset_group(data_source))
}) |> bind_rows() |> mutate(tool = "StrainPhlAn4")

## ---- 2. inStrain: one row per genome-sample detection ----
instrain_counts <- read_csv(file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"), show_col_types = FALSE) |>
  filter(!is.na(classification), !is.na(category.y)) |>
  transmute(species = classification, category = category.y,
            dataset = dataset_group(data_source.y), tool = "inStrain")

## ---- 3. dRep: QC-passing MAGs per species ----
# WHY: mag_metadata_v3 covers every genome in Cdb.csv, but has no CheckM values for the
# published (Walsh, Samuel) MAGs. The HQ representative-MAG file covers fewer genomes but
# carries those QC values, so they are taken from there where v3 lacks them.
mag_pro <- read_csv(file.path(DATA_DIR, "mag_metadata_v3.csv.gz"), show_col_types = FALSE) |>
  transmute(user_genome, classification = sub(".*;s__", "", classification),
            data_source, category, Completeness, Contamination) |>
  distinct()
# HOW: genomes repeat in mag_metadata_v3; after distinct() each must have one record.
stopifnot("a genome has conflicting records in mag_metadata_v3" = !anyDuplicated(mag_pro$user_genome))
rep_qc <- read_csv(file.path(DATA_DIR, "HQ_prokaryotic_representatives_MAGs_per_sample.csv"), show_col_types = FALSE) |>
  distinct(user_genome, rep_completeness = Completeness, rep_contamination = Contamination)
stopifnot("a genome has conflicting QC values in the representative-MAG file" = !anyDuplicated(rep_qc$user_genome))
mag_pro <- mag_pro |>
  left_join(rep_qc, by = "user_genome", relationship = "one-to-one") |>
  mutate(Completeness  = coalesce(Completeness, rep_completeness),
         Contamination = coalesce(Contamination, rep_contamination)) |>
  select(-rep_completeness, -rep_contamination)
no_qc <- mag_pro |> filter(user_genome %in% sub("\\.fa$", "", read_csv(file.path(DATA_DIR, "Cdb.csv"), show_col_types = FALSE)$genome),
                           is.na(Completeness))
message("dRep genomes with no QC values (excluded): ", nrow(no_qc))
print(count(no_qc, data_source))
mag_euk <- read_csv(file.path(DATA_DIR, "HQ_eukaryotic_MAGs.csv"), show_col_types = FALSE) |>
  transmute(user_genome = `3`, classification = sub(".*;", "", X2),
            data_source = ifelse(grepl("_cs_bin", type), "This study", "Samuel et al"),
            category = "Water.kefir", Completeness = complete, Contamination = duplicated)
cdb <- read_csv(file.path(DATA_DIR, "Cdb.csv"), show_col_types = FALSE) |>
  mutate(genome = sub("\\.fa$", "", genome))
drep_counts <- cdb |>
  inner_join(bind_rows(mag_euk, mag_pro), by = c("genome" = "user_genome"), relationship = "one-to-one") |>
  filter(Completeness >= 80, Contamination <= 5) |>
  transmute(species = classification, category, dataset = dataset_group(data_source), tool = "dRep")

## ---- Combine, restrict to prevalent species, total across datasets ----
detections <- bind_rows(strainphlan_counts, instrain_counts, drep_counts) |>
  filter(!is.na(dataset), category %in% names(kefir_type_label)) |>
  mutate(kefir_type = kefir_type_label[category]) |>
  filter((kefir_type == "Milk kefir"  & species %in% prevalent_milk) |
         (kefir_type == "Water kefir" & species %in% prevalent_water))
strain_counts <- detections |>
  count(kefir_type, species, tool, dataset) |>
  pivot_wider(names_from = dataset, values_from = n, values_fill = 0) |>
  mutate(total = Kefir4All + Published) |>
  arrange(kefir_type, desc(total))
stopifnot("no detections left after restricting to prevalent species" = nrow(strain_counts) > 0)
write_csv(strain_counts, file.path(OUT_DIR, "Figure_6_method_comparison_data.csv"))
print(count(strain_counts, kefir_type, tool), n = Inf)

## ---- Panel B data: StrainPhlAn4 marker-sequence polymorphism per sample ----
# One value per sample x clade: the percentage of polymorphic sites across that clade's
# marker-gene alignment. Faceted by dataset, so published samples are kept separate here.
dataset_label <- function(data_source, category) {
  case_when(
    data_source == "This study" & category == "Milk.kefir"  ~ "Milk kefir - Kefir4All",
    data_source == "This study" & category == "Water.kefir" ~ "Water kefir - Kefir4All",
    data_source %in% c("Walsh et al", "Walsh et al 2023") ~ "Milk kefir - Walsh et al",
    data_source %in% c("Mortensen et al 2023", "Samuel et al") ~ "Water kefir - Breselge et al",
    TRUE ~ NA_character_)
}
sample_metadata <- total_metadata |>
  filter(!is.na(category)) |>
  distinct(merge_column, .keep_all = TRUE)
polymorphic_files <- list.files(file.path(DATA_DIR, "06_polymorphic"), pattern = "\\.polymorphic$", full.names = TRUE)
if (length(polymorphic_files) == 0) stop("No .polymorphic files found in ", file.path(DATA_DIR, "06_polymorphic"))
polymorphism_raw <- lapply(polymorphic_files, function(f) {
  read_tsv(f, show_col_types = FALSE) |>
    transmute(clade = sub("\\.polymorphic$", "", basename(f)),
              sample = gsub("-", "_", sample),
              polymorphic_sites_percent = percentage_of_polymorphic_sites)
}) |> bind_rows()
# HOW: prevalence is decided per kefir type from the clade IDs, as for panel A.
prevalent_clades <- bind_rows(
  tibble(clade = milk_prevalence$specie,  category = "Milk.kefir"),
  tibble(clade = water_prevalence$specie, category = "Water.kefir"))
polymorphism <- polymorphism_raw |>
  inner_join(sample_metadata, by = c("sample" = "merge_column"), relationship = "many-to-one") |>
  semi_join(prevalent_clades, by = c("clade", "category")) |>
  left_join(clade_to_species, by = "clade", relationship = "many-to-one") |>
  mutate(dataset = dataset_label(data_source, category)) |>
  filter(!is.na(dataset), !is.na(polymorphic_sites_percent)) |>
  select(dataset, species, clade, sample, polymorphic_sites_percent)
message("Polymorphism rows: ", nrow(polymorphism_raw), " read; ", nrow(polymorphism), " kept (prevalent species, samples with metadata)")
stopifnot("no polymorphism rows kept" = nrow(polymorphism) > 0)
write_csv(polymorphism, file.path(OUT_DIR, "Figure_6_marker_polymorphism_data.csv"))
print(count(polymorphism, dataset), n = Inf)

## ---- Plot ----
library(patchwork)
figure_theme_6 <- theme_bw(base_size = 12) +
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        plot.tag = element_text(size = 12, face = "bold"))

panel_a <- strain_counts |>
  mutate(species = reorder(species, total, FUN = max)) |>
  ggplot(aes(x = species, y = total, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  facet_wrap(~ kefir_type, scales = "free") +
  scale_fill_manual(values = c(dRep = "salmon", inStrain = "#2CA02C", StrainPhlAn4 = "#1F77B4")) +
  labs(x = "Species", y = "Number of strains detected across strain-profiling methods", fill = "Method") +
  coord_flip() +
  figure_theme_6 +
  theme(legend.position = "top")

panel_b <- ggplot(polymorphism, aes(x = polymorphic_sites_percent, y = species)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.6) +
  facet_wrap(~ dataset, scales = "free") +
  labs(x = "Marker-sequence polymorphism (%)", y = "Species") +
  figure_theme_6

# WHY: panel B holds two rows of facets with up to 18 species each, so it gets more height.
figure_6 <- panel_a / panel_b + plot_layout(heights = c(1, 1.6)) + plot_annotation(tag_levels = "A")
save_figure(figure_6, file.path(FIGURES_DIR, "Figure_6"),
            width_mm = 406.4, height_mm = 380,
            source = file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"),
            overwrite = TRUE, tiff = TRUE)
