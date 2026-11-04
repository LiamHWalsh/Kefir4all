# Ported from a pre-repo personal-workstation script (04_mentel_compositonal_parameters.R,
# 2024-06-27) that referenced an external Q:/H2020 Master/... network drive. Reuses the
# metadata-loading pattern already proven in 04_community_stability_v2.R (merge_column <->
# survey-response crosswalk), adapted to build a Gower distance matrix over fermentation
# metadata and compare it against compositional Bray-Curtis dissimilarity via a Mantel test.
#
# Produces: Supplementary Figure 4 -- "Hexagonal heatmap visualizing the correlation between
# the distance in the compositional Bray-Curtis dissimilarity distance matrix and the metadata
# Gower distance matrix" (Kefir4All_Manuscript_TRACKED_SB.docx).

for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {
  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}

if (!requireNamespace("here",   quietly = TRUE)) install.packages("here")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
kefir_use_style()
DATA_DIR    <- here::here("data")
PRIVATE_DIR <- here::here("data", "private")
OUT_DIR     <- here::here("output", "04_taxonomic_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "04_taxonomic_profiling")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

SURVEY_MK <- file.path(PRIVATE_DIR, "fermentation_survey_mk_edited_dated.txt")
SURVEY_WK <- file.path(PRIVATE_DIR, "fermentation_survey_wk_edited_dated.txt")
if (!file.exists(SURVEY_MK) || !file.exists(SURVEY_WK)) {
  stop("This script requires the private fermentation survey files.\n",
       "Missing: ", SURVEY_MK, "\n", "Missing: ", SURVEY_WK, "\n",
       "See data/private/README.md for details.")
}

library(readr); library(readxl); library(dplyr); library(tidyr); library(tibble)
library(vegan); library(ggplot2)
pacman::p_load(StatMatch)
if (!requireNamespace("hexbin", quietly = TRUE)) install.packages("hexbin")

## ---- Per-metagenome metadata, keyed by merge_column (same crosswalk as 04_community_stability_v2.R) ----
kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"), show_col_types = FALSE)
kefir4all_metadata$merge_column <- gsub("_host_removed_R..fastq.gz", "", kefir4all_metadata$merge_column)
kefir4all_metadata <- kefir4all_metadata[-which(duplicated(kefir4all_metadata$merge_column)), ]
kefir4all_metadata <- kefir4all_metadata[kefir4all_metadata$`kefir type` != "Medium control", ]

myfiles <- list(
  mk = read_excel(file.path(DATA_DIR, "file_names_survey_responses_mk.xlsx")),
  wk = read_excel(file.path(DATA_DIR, "file_names_survey_responses_wk.xlsx"))
)
kefir4all_metadata <- merge(
  kefir4all_metadata,
  rbind(dplyr::select(myfiles[["mk"]], merge_column, observations, category_confirmed),
        dplyr::select(myfiles[["wk"]], merge_column, observations, category_confirmed)),
  by = "merge_column", all.x = TRUE)

## ---- Add fermentation-condition fields from the private survey files, joined by participant + stage ----
stage_to_week <- c(T0 = "T0", T1 = "wk01", T2 = "wk05", T3 = "wk09", T4 = "wk13", T5 = "wk17", T6 = "wk21")
kefir4all_metadata$survey_stage <- stage_to_week[kefir4all_metadata$Stage]

survey_mk <- read_delim(SURVEY_MK, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
survey_wk <- read_delim(SURVEY_WK, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
survey_mk$previous_weightofgrains <- as.numeric(gsub("[^0-9.-]", "", iconv(survey_mk$previous_weightofgrains, to = "UTF-8", sub = "")))
survey_wk$previous_weightofgrains <- as.numeric(gsub("[^0-9.-]", "", iconv(survey_wk$previous_weightofgrains, to = "UTF-8", sub = "")))
survey_mk$stage[is.na(survey_mk$stage)] <- "wk21"
survey_wk$stage[is.na(survey_wk$stage)] <- "wk21"
survey_fields <- c("Participant_ID", "stage", "observations", "Room_temperature", "pH",
                   "grain_wash", "previous_weightofgrains", "grains_weight", "Liquid_volume")
survey_all <- rbind(
  dplyr::select(survey_mk, dplyr::all_of(survey_fields)) %>% dplyr::rename(survey_observations = observations),
  dplyr::select(survey_wk, dplyr::all_of(survey_fields)) %>% dplyr::rename(survey_observations = observations)
)
# One fermentation can have multiple logged entries; keep the first per participant+stage.
survey_all <- survey_all[!duplicated(paste(survey_all$Participant_ID, survey_all$stage)), ]

kefir4all_metadata <- merge(
  kefir4all_metadata, survey_all,
  by.x = c("Sample", "survey_stage"), by.y = c("Participant_ID", "stage"), all.x = TRUE)

## ---- Build the metadata matrix for Gower distance (one row per metagenome, by merge_column) ----
meta_vars <- c("kefir type", "Stage", "category_confirmed", "Room_temperature", "pH",
               "grain_wash", "previous_weightofgrains", "grains_weight", "Liquid_volume")
gower_input <- kefir4all_metadata %>%
  dplyr::select(merge_column, dplyr::all_of(meta_vars)) %>%
  dplyr::filter(!duplicated(merge_column)) %>%
  tibble::column_to_rownames("merge_column")
cat("Gower distance computed over", ncol(gower_input), "metadata variables:",
    paste(colnames(gower_input), collapse = ", "), "\n")

gower_distance_matrix <- StatMatch::gower.dist(as.data.frame(unclass(gower_input), stringsAsFactors = TRUE))
rownames(gower_distance_matrix) <- rownames(gower_input)
colnames(gower_distance_matrix) <- rownames(gower_input)

## ---- Compositional (Bray-Curtis) distance over the same samples ----
metacache <- read_csv(file.path(DATA_DIR, "04_metacache_total_full_profile.csv"), show_col_types = FALSE)
metacache[is.na(metacache)] <- 0
metacache <- metacache %>% remove_rownames() %>% column_to_rownames("sample_id")

shared <- intersect(rownames(metacache), rownames(gower_distance_matrix))
cat("Samples with both compositional and metadata data:", length(shared), "\n")
dist_abund <- vegan::vegdist(metacache[shared, , drop = FALSE], method = "bray", na.rm = TRUE)
gower_distance_matrix <- gower_distance_matrix[shared, shared]
stopifnot(identical(rownames(as.matrix(dist_abund)), rownames(gower_distance_matrix)))

mantel_result <- vegan::mantel(as.matrix(dist_abund), as.matrix(gower_distance_matrix),
                               method = "spearman", permutations = 9999, na.rm = TRUE)
cat("Mantel test: r =", mantel_result$statistic, " p =", mantel_result$signif, "\n")

plot_data <- data.frame(
  dist_abund = as.vector(as.dist(dist_abund)),
  gower_distance = as.vector(as.dist(gower_distance_matrix))
)

p <- ggplot(plot_data, aes(x = dist_abund, y = gower_distance)) +
  geom_hex() +
  scale_fill_viridis_c(name = "Count") +
  geom_smooth(method = "lm", colour = "blue", se = FALSE, linewidth = 1) +
  labs(title = "Mantel test: correlation between distance matrices",
       subtitle = sprintf("Spearman r = %.3f, p = %.4g (9999 permutations, n = %d samples)",
                           mantel_result$statistic, mantel_result$signif, length(shared)),
       x = "Bray-Curtis dissimilarity (compositional)",
       y = "Gower distance (fermentation metadata)") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_4_mantel_gower.tiff"), p,
       width = 10, height = 8, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_4_mantel_gower.pdf"), p,
       width = 10, height = 8, bg = "white")
write_csv(plot_data, file.path(OUT_DIR, "Supplementary_Figure_4_mantel_gower_data.csv"))

message("Wrote Supplementary_Figure_4_mantel_gower.{tiff,pdf} to ", FIGURES_DIR)
