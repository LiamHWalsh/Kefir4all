# Produces: Figure 6 -- "Method-specific detection counts and marker-sequence polymorphism
# across prevalent kefir species. Three approaches were used to examine aspects of strain-level
# variation: diversity recovered through genome assembly, reference-genome detection and
# assembly-free marker-sequence variation. (A) Species-level summaries obtained using dRep
# (salmon), inStrain (green) and StrainPhlAn4 (blue)..." (Kefir4All_Manuscript_TRACKED_SB.docx).
#
# Direct port of the original personal-workstation script "number of strains detected.R"
# (found at Kefir4all/temp_003/Temp_4Liam/R scripts/evolutionary/06_strain_profiling/ --
# never committed to this repo, so not visible to earlier provenance checks). Reproduces that
# script's Panel A exactly, including its data-source pooling: counts are NOT restricted to
# this study's own samples -- each method's bars also include matched published comparison
# datasets (Walsh et al 2023 for milk, Mortensen et al 2023 / Samuel et al for water), with
# data points labelled by both kefir type and data source (e.g. "Milk kefir - Kefir4All" vs
# "Milk kefir - Walsh et al"). This is a deliberate reversion to the original script's scope.
#
# Panel B (StrainPhlAn4 marker-polymorphism boxplots) is NOT reproduced here: the original
# script computes it from per-clade .polymorphic files (percentage_of_polymorphic_sites),
# which do not exist anywhere in this repo or on the local machine the original script was
# found on -- confirmed by search, not assumed. Panel A alone is produced.
#
# Species names are harmonised across all three methods via the prevalence crosswalk
# (data/milk_taxonomic_profile_prevalence.csv / water_taxonomic_profile_prevalence.csv), and
# each method's output is restricted to species prevalent in that kefir type.

# WHY: here::here() resolves relative to its own root-detection markers (.Rproj, .git),
# which are not present in the Docker image's working copy -- this sets the working
# directory to wherever the repo was actually mounted before here::here() is called.
for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {
  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("ape",  quietly = TRUE)) install.packages("ape")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style()
DATA_DIR    <- here::here("data")
TREES_DIR   <- here::here("data", "06_trees")
OUT_DIR     <- here::here("output_reviewed", "06_strain_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "Figure_6")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

library(readr); library(dplyr); library(tidyr); library(stringr); library(ape); library(ggplot2)

## ---- Prevalent-species crosswalk: t__SGBxxxx <-> species name, and the per-type prevalence lists ----
# The prevalence tables split a few species into dRep/StrainPhlAn-only
# "subcluster"/"subspecies" entries that inStrain and dRep report under a
# single merged name (confirmed in data/instrain_genome_species_
# primary_data_v4.csv) -- this relabelling pools them under one species name
# instead of leaving them as extra, StrainPhlAn4-only fragments.
merge_subclusters <- function(x) {
  x <- gsub("Lactococcus_lactis subcluster 1", "Lactococcus_lactis", x, fixed = TRUE)
  x <- gsub("Lactococcus_lactis subcluster 2", "Lactococcus_cremoris", x, fixed = TRUE)
  x <- gsub("Pseudomonas_fragi_subspecies 1", "Pseudomonas_fragi", x, fixed = TRUE)
  x <- gsub("Zymomonas_mobilis_subcluster 1", "Zymomonas_mobilis", x, fixed = TRUE)
  # Same post-2020 genus reclassifications the original script applied directly to these
  # three species labels (Lactobacillus -> Liquorilactobacillus / Schleiferilactobacillus).
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

## ---- Combined study + published-comparison metadata (legacy pools across all of these) ----
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

# Legacy's data-source -> pretty-label mapping, applied identically to all three methods so
# the same kefir-type/data-source combination reads the same way in every bar chart.
pretty_source <- function(data_source, category) {
  dplyr::case_when(
    data_source == "This study" & category == "Milk.kefir"  ~ "Milk kefir - Kefir4All",
    data_source == "This study" & category == "Water.kefir" ~ "Water kefir - Kefir4All",
    data_source %in% c("Walsh et al", "Walsh et al 2023")    ~ "Milk kefir - Walsh et al",
    data_source %in% c("Mortensen et al 2023", "Samuel et al") ~ "Water kefir - Breselge et al",
    TRUE ~ NA_character_)
}

## ---- 1. StrainPhlAn4: tree tip labels, categorised by data source (legacy used per-clade
## .mutation files for this count; those files do not exist in this repo or on the machine the
## original script was found on, so tree tip labels -- which encode the same "was this clade
## profiled in this sample" information -- are used instead, same substitution 06_figure6's
## prior version already made for the same reason). ----
tree_files <- list.files(TREES_DIR, pattern = "\\.tre$", full.names = TRUE)
if (length(tree_files) == 0) {
  stop("No StrainPhlAn4 tree files found in ", TREES_DIR,
       " -- Figure 6 needs all three methods; check TREES_DIR and the .tre file naming")
}
strainphlan_counts <- lapply(tree_files, function(f) {
  clade <- sub("^RAxML_bestTree\\.(t__.+)\\.StrainPhlAn4\\.tre$", "\\1", basename(f))
  species <- clade_to_species$species[clade_to_species$clade == clade]
  if (!length(species)) return(NULL)
  tips <- gsub("-", "_", ape::read.tree(f)$tip.label)
  hit <- total_metadata[total_metadata$merge_column %in% tips & !is.na(total_metadata$category), ]
  if (!nrow(hit)) return(NULL)
  hit$data_source.sample_id <- pretty_source(hit$data_source, hit$category)
  hit |> filter(!is.na(data_source.sample_id)) |>
    count(data_source.sample_id, name = "Freq") |>
    transmute(species = species[1], data_source.sample_id, Freq, tool = "StrainPhlAn4")
}) |> bind_rows()
stopifnot("no StrainPhlAn4 detections found across any tree -- check TREES_DIR tip-label parsing" =
            nrow(strainphlan_counts) > 0)

## ---- 2. inStrain: detection rows per species x data source, same pooling as StrainPhlAn4 ----
instrain <- read_csv(file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"), show_col_types = FALSE)
instrain_counts <- instrain |>
  filter(!is.na(classification), !is.na(category.y)) |>
  mutate(data_source.sample_id = pretty_source(data_source.y, category.y)) |>
  filter(!is.na(data_source.sample_id)) |>
  count(classification, data_source.sample_id, name = "Freq") |>
  transmute(species = classification, data_source.sample_id, Freq, tool = "inStrain")
stopifnot("no inStrain detections found -- check data_source.y/category.y values" =
            nrow(instrain_counts) > 0)

## ---- 3. dRep: distinct genomes (MAGs) per species x data source, same pooling, same prokaryotic
## + eukaryotic MAG source files and QC filter (Completeness >= 80, Contamination <= 5) the
## original script used. ----
mag_pro <- read_csv(file.path(DATA_DIR, "HQ_prokaryotic_representatives_MAGs_per_sample.csv"), show_col_types = FALSE) |>
  transmute(user_genome, classification = sub(".*;s__", "", classification),
            data_source, category, Completeness, Contamination)
mag_euk <- read_csv(file.path(DATA_DIR, "HQ_eukaryotic_MAGs.csv"), show_col_types = FALSE)
names(mag_euk)[names(mag_euk) == "3"] <- "user_genome"
names(mag_euk)[names(mag_euk) == "X2"] <- "classification_full"
names(mag_euk)[names(mag_euk) == "complete"] <- "Completeness"
names(mag_euk)[names(mag_euk) == "duplicated"] <- "Contamination"
mag_euk <- mag_euk |>
  mutate(classification = sub(".*;", "", classification_full),
         data_source = ifelse(grepl("_cs_bin", type), "This study", "Samuel et al"),
         category = "Water.kefir") |>
  transmute(user_genome, classification, data_source, category, Completeness, Contamination)
mag_metadata <- bind_rows(mag_euk, mag_pro)

cdb <- read_csv(file.path(DATA_DIR, "Cdb.csv"), show_col_types = FALSE) |>
  mutate(genome = sub("\\.fa$", "", genome))
drep_joined <- cdb |>
  inner_join(mag_metadata, by = c("genome" = "user_genome"), relationship = "many-to-one") |>
  filter(Completeness >= 80, Contamination <= 5) |>
  mutate(data_source.sample_id = pretty_source(data_source, category)) |>
  filter(!is.na(data_source.sample_id))
drep_counts <- drep_joined |>
  count(classification, data_source.sample_id, name = "Freq") |>
  transmute(species = classification, data_source.sample_id, Freq, tool = "dRep")
stopifnot("no dRep genomes found -- check data_source/category values in the MAG source files" =
            nrow(drep_counts) > 0)

## ---- Combine, restrict to species prevalent in that kefir type ----
test_db <- bind_rows(strainphlan_counts, instrain_counts, drep_counts) |>
  mutate(species = ifelse(species == "Pseudomonas fragi subspecies 1", "Pseudomonas fragi", species),
         new_type = sub(" -.*", "", data_source.sample_id)) |>
  filter(Freq > 0) |>
  filter((new_type == "Milk kefir"  & species %in% prevalent_milk) |
         (new_type == "Water kefir" & species %in% prevalent_water))
stopifnot("no rows left after restricting to prevalent species -- check prevalence crosswalk" =
            nrow(test_db) > 0)
write_csv(test_db, file.path(OUT_DIR, "Figure_6_method_comparison_data.csv"))
message("Species x data-source x tool combinations: ", nrow(test_db))
print(count(test_db, tool))

p <- test_db |>
  mutate(species = reorder(species, Freq, FUN = max)) |>
  ggplot(aes(x = species, y = Freq, fill = tool)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single"), width = 0.7) +
  facet_wrap(~ new_type, scales = "free") +
  scale_fill_manual(values = c("dRep" = "salmon", "inStrain" = "#2CA02C", "StrainPhlAn4" = "#1F77B4")) +
  labs(x = "Species", y = "Number of strains detected across strain-profiling methods", fill = "Method") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top",
        strip.text = element_text(size = 14, face = "bold", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"))

ggsave(file.path(FIGURES_DIR, "Figure_6_method_comparison.tiff"), p,
       width = 16, height = 10, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Figure_6_method_comparison.pdf"), p,
       width = 16, height = 10, bg = "white")

# NOTE: this figure combines three roughly-equal-importance sources (dRep, inStrain,
# StrainPhlAn4 tree tips); there is no single "primary" input. instrain_genome_species_
# primary_data_v4.csv is used here as the representative source for provenance purposes.
save_figure(p, file.path(FIGURES_DIR, "Figure_6_method_comparison_provenance"),
            width_mm = 16 * 25.4, height_mm = 10 * 25.4,
            source = file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"),
            seed = NULL)

message("Wrote Figure_6_method_comparison.{tiff,pdf} to ", FIGURES_DIR)
