# Produces: Figure 6 -- "Method-specific detection counts and marker-sequence polymorphism
# across prevalent kefir species. Three approaches were used to examine aspects of strain-level
# variation: diversity recovered through genome assembly, reference-genome detection and
# assembly-free marker-sequence variation. (A) Species-level summaries obtained using dRep
# (salmon), inStrain (green) and StrainPhlAn4 (blue)..." (Kefir4All_Manuscript_TRACKED_SB.docx).
#
# No script in this repo combines all three methods. Ported/adapted from a 2024 pre-repo
# personal-workstation script ("number of strains detected.R") that referenced an external
# Q:/ network drive, using:
#   - dRep: Cdb.csv + mag_metadata_v3.csv.gz, joined exactly as 06_supp_strain_profiling.R
#     already does for Figure 8's dRep panel -- distinct genomes (MAGs) per species x type.
#   - inStrain: instrain_genome_species_primary_data_v4.csv -- total detection rows per
#     species x type (same "Total Detections" definition used in 06_strain_profiling.R).
#   - StrainPhlAn4: detection counts use per-clade tree tip labels, and marker polymorphism
#     uses the downloaded per-clade .mutation matrices.
#
# Species names are harmonised across all three methods via the existing prevalence crosswalk
# (data/milk_taxonomic_profile_prevalence.csv / water_taxonomic_profile_prevalence.csv, the
# same "t__SGBxxxx" -> species-name table already used for Supplementary Figures 4 and 5),
# and each method's output is restricted to species prevalent in that kefir type, the
# convention used throughout this pipeline.

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
source(here::here("scripts_reviewed_and_figure_improved", "_lab_toolkit", "plotting_report.R"))
kefir_use_style()
DATA_DIR    <- here::here("data")
TREES_DIR   <- here::here("data", "06_trees")
OUT_DIR     <- here::here("output_reviewed", "06_strain_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures_ai_skill", "Figure_6")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

library(readr); library(dplyr); library(tidyr); library(ape); library(ggplot2)
library(patchwork)

## ---- Prevalent-species crosswalk: t__SGBxxxx <-> species name, and the per-type prevalence lists ----
# The prevalence tables split a few species into dRep/StrainPhlAn-only
# "subcluster"/"subspecies" entries that inStrain and dRep already report
# under a single merged name (confirmed in data/instrain_genome_species_
# primary_data_v4.csv) -- without this relabelling (same mapping the legacy
# script used), these show up as extra, StrainPhlAn4-only fragments instead
# of being pooled with the other two methods under one species name.
merge_subclusters <- function(x) {
  x <- gsub("Lactococcus_lactis subcluster 1", "Lactococcus_lactis", x, fixed = TRUE)
  x <- gsub("Lactococcus_lactis subcluster 2", "Lactococcus_cremoris", x, fixed = TRUE)
  x <- gsub("Pseudomonas_fragi_subspecies 1", "Pseudomonas_fragi", x, fixed = TRUE)
  x <- gsub("Zymomonas_mobilis_subcluster 1", "Zymomonas_mobilis", x, fixed = TRUE)
  x
}
milk_prevalence  <- read_csv(file.path(DATA_DIR, "milk_taxonomic_profile_prevalence.csv"),  show_col_types = FALSE)
water_prevalence <- read_csv(file.path(DATA_DIR, "water_taxonomic_profile_prevalence.csv"), show_col_types = FALSE)
stopifnot("expected prevalence columns 'specie'/'...5' not found -- readr's name repair may have changed" =
            all(c("specie", "...5") %in% names(milk_prevalence)) &&
            all(c("specie", "...5") %in% names(water_prevalence)))
milk_prevalence$`...5`  <- merge_subclusters(milk_prevalence$`...5`)
water_prevalence$`...5` <- merge_subclusters(water_prevalence$`...5`)
clade_to_species <- bind_rows(milk_prevalence, water_prevalence) |>
  transmute(clade = specie, species = gsub("_", " ", `...5`)) |>
  distinct(clade, .keep_all = TRUE)
prevalent_milk  <- unique(gsub("_", " ", milk_prevalence$`...5`))
prevalent_water <- unique(gsub("_", " ", water_prevalence$`...5`))

## ---- 1. StrainPhlAn4: tree tip labels present in this study, per clade x type ----
kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"), show_col_types = FALSE)
kefir4all_metadata$merge_column <- gsub("_host_removed_R..fastq.gz", "", kefir4all_metadata$merge_column)
kefir4all_metadata <- kefir4all_metadata[!duplicated(kefir4all_metadata$merge_column), ]
kefir4all_metadata$category.y <- case_when(
  kefir4all_metadata$`kefir type` %in% c("ML", "MG") ~ "Milk.kefir",
  kefir4all_metadata$`kefir type` %in% c("WL", "WG") ~ "Water.kefir",
  TRUE ~ NA_character_)

tree_files <- list.files(TREES_DIR, pattern = "\\.tre$", full.names = TRUE)
if (length(tree_files) == 0) {
  stop("No StrainPhlAn4 tree files found in ", TREES_DIR,
       " -- Figure 6 needs all three methods; check TREES_DIR and the .tre file naming")
}
strainphlan_counts <- lapply(tree_files, function(f) {
  clade <- sub("^RAxML_bestTree\\.(t__.+)\\.StrainPhlAn4\\.tre$", "\\1", basename(f))
  species <- clade_to_species$species[clade_to_species$clade == clade]
  if (!length(species)) return(NULL)
  tips <- ape::read.tree(f)$tip.label
  tips <- gsub("-", "_", tips)
  hit <- kefir4all_metadata[kefir4all_metadata$merge_column %in% tips & !is.na(kefir4all_metadata$category.y), ]
  if (!nrow(hit)) return(NULL)
  as.data.frame(table(category.y = hit$category.y)) |>
    transmute(species = species[1], category.y = as.character(category.y), Freq, tool = "StrainPhlAn4")
}) |> bind_rows()
stopifnot("no StrainPhlAn4 detections found across any tree -- check TREES_DIR tip-label parsing" =
            nrow(strainphlan_counts) > 0)

## ---- 2. inStrain: total detection rows per species x type ----
# Restricted to this study's own samples, same as dRep (data_source == "This
# study") and StrainPhlAn4 (tips filtered to kefir4all_metadata) below --
# without this, inStrain pooled in the published comparison datasets
# (Walsh et al / Mortensen et al 2023) while the other two methods didn't,
# inflating its bars relative to dRep/StrainPhlAn4 for the same species.
instrain <- read_csv(file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"), show_col_types = FALSE)
instrain_counts <- instrain |>
  filter(!is.na(classification), !is.na(category.y), data_source.y == "This study") |>
  count(classification, category.y, name = "Freq") |>
  transmute(species = classification, category.y, Freq, tool = "inStrain")
stopifnot("no inStrain 'This study' detections found -- check data_source.y values" =
            nrow(instrain_counts) > 0)

## ---- 3. dRep: distinct genomes (MAGs) per species x type (same join as 06_supp_strain_profiling.R) ----
cdb <- read_csv(file.path(DATA_DIR, "Cdb.csv"), show_col_types = FALSE)
md  <- read_csv(file.path(DATA_DIR, "mag_metadata_v3.csv.gz"), show_col_types = FALSE, guess_max = 50000)
strip_fa <- function(x) sub("\\.fa$", "", as.character(x))
cdb <- cdb |> mutate(g = strip_fa(genome))
md  <- md  |> mutate(g = strip_fa(user_genome))
species_short <- function(s) {
  m <- regmatches(s, regexpr("s__[^;]+$", s))
  ifelse(length(m), gsub("_", " ", trimws(sub("^s__", "", m))), s)
}
# WHY: mag_metadata_v3.csv.gz carries exact-duplicate metadata rows for ~42% of genomes
# (identical classification/data_source/kefir type, differing only in unrelated columns such
# as fastani stats) -- confirmed empirically: distinct() on these four columns alone collapses
# to exactly one row per genome with zero conflicting values. Deduplicating here, rather than
# just joining, avoids silently double-counting those genomes in the dRep bars below.
md_genome_lookup <- md |>
  select(g, classification, data_source, `kefir type`) |>
  distinct()
stopifnot("duplicate (g) after distinct -- some genomes have genuinely conflicting metadata" =
            !anyDuplicated(md_genome_lookup$g))

# relationship = "many-to-one": md_genome_lookup now has exactly one row per genome.
drep_joined <- cdb |>
  inner_join(md_genome_lookup, by = "g", relationship = "many-to-one") |>
  filter(data_source == "This study") |>
  mutate(species = vapply(classification, species_short, character(1)),
         category.y = case_when(
           `kefir type` %in% c("ML", "MG") ~ "Milk.kefir",
           `kefir type` %in% c("WL", "WG") ~ "Water.kefir",
           TRUE ~ NA_character_)) |>
  filter(!is.na(category.y))
drep_counts <- drep_joined |>
  count(species, category.y, name = "Freq") |>
  transmute(species, category.y, Freq, tool = "dRep")
stopifnot("no dRep 'This study' genomes found -- check data_source values in mag_metadata_v3.csv.gz" =
            nrow(drep_counts) > 0)

## ---- Panel B: StrainPhlAn4 marker-sequence polymorphism ----
global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv"),
                               show_col_types = FALSE)
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv"),
                               show_col_types = FALSE)
total_metadata <- bind_rows(
  select(kefir4all_metadata, merge_column, `kefir type`),
  select(global_mk_metadata, merge_column, `kefir type`),
  select(global_wk_metadata, merge_column, `kefir type`)
) |>
  distinct(merge_column, .keep_all = TRUE) |>
  mutate(category = case_when(
    `kefir type` %in% c("ML", "MG") ~ "Milk.kefir",
    `kefir type` %in% c("WL", "WG") ~ "Water.kefir",
    TRUE ~ NA_character_
  ))

marker_names <- read_csv(file.path(DATA_DIR, "marker_names.csv"), show_col_types = FALSE)
stopifnot("marker_names.csv must contain x and ...2 columns" =
            all(c("x", "...2") %in% names(marker_names)))
marker_names <- marker_names |>
  transmute(clade = paste0("t__", x), species = gsub("_", " ", ...2))
mutation_files <- list.files(file.path(DATA_DIR, "06_mutation"),
                             pattern = "\\.mutation$", full.names = TRUE)
stopifnot("no StrainPhlAn4 mutation files found in data/06_mutation" =
            length(mutation_files) > 0)

parse_mutation_ratio <- function(values) {
  vapply(strsplit(as.character(values), "/", fixed = TRUE), function(parts) {
    if (length(parts) == 1) return(suppressWarnings(as.numeric(parts[[1]])))
    numerator <- suppressWarnings(as.numeric(parts[[1]]))
    denominator <- suppressWarnings(as.numeric(parts[[2]]))
    if (is.na(numerator) || is.na(denominator) || denominator == 0) return(NA_real_)
    numerator / denominator
  }, numeric(1))
}

mutation_data <- lapply(mutation_files, function(path) {
  clade <- sub("\\.mutation$", "", basename(path))
  species <- marker_names$species[marker_names$clade == clade]
  if (length(species) != 1) return(NULL)
  matrix_data <- read_tsv(
    path,
    col_types = cols(.default = col_character()),
    na = c("", "NA"),
    show_col_types = FALSE
  )
  names(matrix_data)[1] <- "ids"
  matrix_data |>
    mutate(ids = gsub("-", "_", ids)) |>
    pivot_longer(-ids, names_to = "sample_id", values_to = "snps") |>
    mutate(
      sample_id = gsub("-", "_", sample_id),
      snps_updated = parse_mutation_ratio(snps),
      species = species
    )
}) |>
  bind_rows()
stopifnot("mutation matrices produced no rows" = nrow(mutation_data) > 0)

mutation_plot_data <- mutation_data |>
  left_join(select(total_metadata, merge_column, category),
            by = c("ids" = "merge_column"), relationship = "many-to-one") |>
  left_join(select(total_metadata, merge_column, category),
            by = c("sample_id" = "merge_column"),
            relationship = "many-to-one", suffix = c(".x", ".y")) |>
  filter(!is.na(category.x), !is.na(snps_updated)) |>
  filter(
    (category.x == "Milk.kefir" & species %in% prevalent_milk) |
      (category.x == "Water.kefir" & species %in% prevalent_water)
  ) |>
  mutate(
    species = gsub("Lactobacillus ghanensis", "Liquorilactobacillus ghanensis",
                   species, fixed = TRUE),
    species = gsub("Lactobacillus satsumensis", "Liquorilactobacillus satsumensis",
                   species, fixed = TRUE),
    species = gsub("Lactobacillus perolens", "Schleiferilactobacillus perolens",
                   species, fixed = TRUE),
    category = gsub("\\.", " ", category.x)
  )
stopifnot("no prevalent StrainPhlAn4 polymorphism values found" =
            nrow(mutation_plot_data) > 0)

panel_b <- ggplot(mutation_plot_data, aes(x = species, y = 100 * snps_updated)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ category, scales = "free") +
  labs(x = "Species", y = "Marker-sequence polymorphism (%)") +
  coord_flip() +
  theme_bw() +
  theme(
    strip.text = element_text(size = 18, face = "bold", colour = "black"),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size = 10, face = "italic"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 14, face = "bold")
  )

## ---- Combine, restrict to species prevalent in that kefir type ----
test_db <- bind_rows(strainphlan_counts, instrain_counts, drep_counts) |>
  filter(Freq > 0) |>
  filter((category.y == "Milk.kefir"  & species %in% prevalent_milk) |
         (category.y == "Water.kefir" & species %in% prevalent_water)) |>
  mutate(category.y = gsub("\\.", " ", category.y))
stopifnot("no rows left after restricting to prevalent species -- check prevalence crosswalk" =
            nrow(test_db) > 0)
write_csv(test_db, file.path(OUT_DIR, "Figure_6_method_comparison_data.csv"))
message("Species x type x tool combinations: ", nrow(test_db))
print(count(test_db, tool))

panel_a <- test_db |>
  mutate(species = reorder(species, Freq, FUN = max)) |>
  ggplot(aes(x = species, y = Freq, fill = tool)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single"), width = 0.7) +
  facet_wrap(~ category.y, scales = "free") +
  scale_fill_manual(values = c("dRep" = "salmon", "inStrain" = "#2CA02C", "StrainPhlAn4" = "#1F77B4")) +
  labs(x = "Species", y = "Number of strains detected across strain-profiling methods", fill = "Method") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top",
        strip.text = element_text(size = 22, face = "bold", colour = "black"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"))

p <- panel_a / panel_b +
  plot_layout(heights = c(1, 1))

# NOTE: this figure combines three roughly-equal-importance sources (dRep, inStrain,
# StrainPhlAn4 tree tips); there is no single "primary" input. instrain_genome_species_
# primary_data_v4.csv is used here as the representative source for provenance purposes.
save_figure(p, file.path(FIGURES_DIR, "Figure_6_method_comparison_provenance"),
            width_mm = 16 * 25.4, height_mm = 10 * 25.4,
            source = file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"),
            seed = NULL)

message("Wrote Figure_6_method_comparison.{tiff,pdf} to ", FIGURES_DIR)
