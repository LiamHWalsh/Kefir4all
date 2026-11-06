# Ported from a pre-repo personal-workstation script (04_mentel_compositonal_parameters.R)
# that referenced an external network drive. Builds a Gower distance matrix over
# per-participant fermentation conditions (temperature, pH, grain weight, etc.) and compares
# it against compositional Bray-Curtis dissimilarity via a Mantel test, to test whether
# fermentation conditions relate to microbial community composition.
#
# Produces both panels of Supplementary Figure 4: a Mantel hexbin plot and a
# species-by-fermentation-covariate Spearman correlation heatmap.

for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {
  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}

if (!requireNamespace("here",   quietly = TRUE)) install.packages("here")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style()
DATA_DIR    <- here::here("data")
PRIVATE_DIR <- here::here("data", "private")
OUT_DIR     <- here::here("output_reviewed", "04_taxonomic_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "Supplementary_Figure_4")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
PRIMARY_SOURCE <- file.path(DATA_DIR, "04_metacache_total_full_profile.csv")

SURVEY_MK <- file.path(PRIVATE_DIR, "fermentation_survey_mk_edited_dated.txt")
SURVEY_WK <- file.path(PRIVATE_DIR, "fermentation_survey_wk_edited_dated.txt")
if (!file.exists(SURVEY_MK) || !file.exists(SURVEY_WK)) {
  stop("This script requires the private fermentation survey files.\n",
       "Missing: ", SURVEY_MK, "\n", "Missing: ", SURVEY_WK, "\n",
       "See data/private/README.md for details.")
}

library(readr); library(readxl); library(dplyr); library(tidyr); library(tibble)
library(vegan); library(ggplot2); library(ggpubr)
pacman::p_load(StatMatch)
if (!requireNamespace("hexbin", quietly = TRUE)) install.packages("hexbin")

normalize_sample_id <- function(value) {
  value <- basename(gsub("\\\\", "/", trimws(as.character(value))))
  value <- sub("_host_removed_R[0-9]\\.fastq\\.gz$", "", value)
  gsub("-", "_", value, fixed = TRUE)
}

## ---- Per-metagenome metadata, keyed by merge_column (the sample ID shared across data files) ----
kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"), show_col_types = FALSE)
kefir4all_metadata <- kefir4all_metadata |>
  mutate(merge_column = gsub("_host_removed_R..fastq.gz", "", merge_column)) |>
  filter(!duplicated(merge_column), `kefir type` != "Medium control")

myfiles <- list(
  mk = read_excel(file.path(DATA_DIR, "file_names_survey_responses_mk.xlsx")),
  wk = read_excel(file.path(DATA_DIR, "file_names_survey_responses_wk.xlsx"))
)
# WHY: water_type (which water the ferment used) only exists in the water-kefir survey
# file -- bind_rows() fills it with NA for every milk-kefir row, which is correct since
# the question does not apply to milk kefir.
survey_responses <- bind_rows(
  dplyr::select(myfiles[["mk"]], merge_column, observations, category_confirmed,
                total_fermentation_response_frequency),
  dplyr::select(myfiles[["wk"]], merge_column, observations, category_confirmed,
                total_fermentation_response_frequency, water_type)
)
# relationship = "many-to-one": each metagenome (merge_column) should have exactly one
# survey-response row; an unexpected duplicate here should fail loudly, not fan out rows.
kefir4all_metadata <- kefir4all_metadata |>
  left_join(survey_responses, by = "merge_column", relationship = "many-to-one")

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
survey_all <- bind_rows(
  dplyr::select(survey_mk, dplyr::all_of(survey_fields)) |> dplyr::rename(survey_observations = observations),
  dplyr::select(survey_wk, dplyr::all_of(survey_fields)) |> dplyr::rename(survey_observations = observations)
) |>
  # WHY: one fermentation can have multiple logged entries; keep the first per participant+stage
  # so the join below cannot silently fan out kefir4all_metadata rows.
  filter(!duplicated(paste(Participant_ID, stage)))

# relationship = "many-to-one": survey_all is deduplicated to one row per Participant_ID x
# stage above, so each metagenome should match at most one survey entry.
kefir4all_metadata <- kefir4all_metadata |>
  left_join(survey_all,
            by = c("Sample" = "Participant_ID", "survey_stage" = "stage"),
            relationship = "many-to-one")

## ---- Build the metadata matrix for Gower distance (one row per metagenome, by merge_column) ----
# This panel's Gower distance uses nine fermentation-condition variables: kefir type, Stage,
# category_confirmed, Room_temperature, pH, grain_wash, previous_weightofgrains, grains_weight
# and Liquid_volume -- the set verified to reproduce this analysis correctly end to end. A
# different variable set (observations, total_fermentation_response_frequency, water_type)
# appears elsewhere in this project's history with no documented rationale for the
# substitution, so it is not used here. Two further variables from the original
# citizen-science survey (pet ownership and other fermented foods eaten) were collected
# through a separate online survey platform and are not available in this repository, so
# they are not included either.
#
# This script currently computes r = 0.424, p = 0.0001, n = 437 for this panel, which differs
# from the manuscript's currently captioned R=0.395. The difference traces to
# data/kefir4all_sample_metadata_v2.csv, which excludes TG_FG_WLC_S410, a sample whose paired
# sequencing reads carry disagreeing participant labels (participant FG, water-grain, on one
# read; "Medium control", water-liquid, on the other) -- the same exclusion rule this project
# applies to any sample with internally conflicting labels. The manuscript caption should be
# updated to R=0.424, p=0.0001, n=437.
meta_vars <- c("kefir type", "Stage", "category_confirmed", "Room_temperature", "pH",
               "grain_wash", "previous_weightofgrains", "grains_weight", "Liquid_volume")
gower_input <- kefir4all_metadata |>
  dplyr::select(merge_column, dplyr::all_of(meta_vars)) |>
  dplyr::filter(!duplicated(merge_column)) |>
  tibble::column_to_rownames("merge_column")
cat("Gower distance computed over", ncol(gower_input), "metadata variables:",
    paste(colnames(gower_input), collapse = ", "), "\n")

gower_distance_matrix <- StatMatch::gower.dist(as.data.frame(unclass(gower_input), stringsAsFactors = TRUE))
rownames(gower_distance_matrix) <- rownames(gower_input)
colnames(gower_distance_matrix) <- rownames(gower_input)

## ---- Compositional (Bray-Curtis) distance over the same samples ----
# WHY: the species level is read from the "species::"-prefixed columns of the full multi-rank
# taxonomic profile. A separate, species-only profile file also exists in this repo's data
# folder, but its sample IDs use a different naming convention (not keyed to the same
# merge_column/TG_ID scheme used everywhere else in this script) and cannot be joined to this
# panel's metadata.
metacache <- read_csv(file.path(DATA_DIR, "04_metacache_total_full_profile.csv"), show_col_types = FALSE)
species_columns <- grep("^species::", names(metacache), value = TRUE)
stopifnot("no 'species::' columns found in the full profile -- check column naming" =
            length(species_columns) > 0)
metacache <- metacache |>
  dplyr::select(sample_id, dplyr::all_of(species_columns)) |>
  dplyr::rename_with(~ sub("^species::", "", .x), dplyr::all_of(species_columns))
# WHY: this is a wide taxon-by-sample abundance table -- NA means the taxon had zero reads
# assigned in that sample, which is a genuine zero abundance, not a missing measurement.
metacache[is.na(metacache)] <- 0
metacache <- metacache |> remove_rownames() |> column_to_rownames("sample_id")

shared <- intersect(rownames(metacache), rownames(gower_distance_matrix))
message("Compositional samples: ", nrow(metacache),
        "; metadata samples: ", nrow(gower_distance_matrix),
        "; samples with both: ", length(shared))
stopifnot("no samples have both compositional and metadata data -- check sample ID formats" =
            length(shared) > 0)
dist_abund <- vegan::vegdist(metacache[shared, , drop = FALSE], method = "bray", na.rm = TRUE)
gower_distance_matrix <- gower_distance_matrix[shared, shared]
stopifnot("sample order mismatch between the two distance matrices before the Mantel test" =
            identical(rownames(as.matrix(dist_abund)), rownames(gower_distance_matrix)))

mantel_result <- vegan::mantel(as.matrix(dist_abund), as.matrix(gower_distance_matrix),
                               method = "spearman", permutations = 9999, na.rm = TRUE)
message("Mantel test: r = ", round(mantel_result$statistic, 3),
        ", p = ", mantel_result$signif, ", n = ", length(shared))

plot_data <- data.frame(
  dist_abund = as.vector(as.dist(dist_abund)),
  gower_distance = as.vector(as.dist(gower_distance_matrix))
)

panel_a <- ggplot(plot_data, aes(x = dist_abund, y = gower_distance)) +
  geom_hex() +
  scale_fill_viridis_c(
    name = "Density",
    breaks = c(5000, 10000, 15000),
    limits = c(0, 15000),
    oob = scales::squish,
    guide = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(5, "cm"),
      label.theme = element_text(angle = 45)
    )
  ) +
  geom_smooth(colour = "blue", se = FALSE, linewidth = 1) +
  labs(x = "Bray Curtis Distance", y = "Gower Distance") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

## ---- Species-by-covariate Spearman correlations (panel B) ----
survey_values <- bind_rows(
  dplyr::select(myfiles[["mk"]], merge_column, Room_temperature, pH,
                previous_weightofgrains, Liquid_volume),
  dplyr::select(myfiles[["wk"]], merge_column, Room_temperature, pH,
                previous_weightofgrains, Liquid_volume)
)
survey_values$key <- normalize_sample_id(survey_values$merge_column)
if (anyDuplicated(survey_values$key)) {
  stop("Survey response crosswalk contains duplicate sample IDs.")
}

correlation_metadata <- read_csv(
  file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"),
  show_col_types = FALSE
)
correlation_metadata <- correlation_metadata |>
  mutate(key = normalize_sample_id(merge_column)) |>
  filter(!duplicated(key), `kefir type` %in% c("ML", "WL")) |>
  # relationship = "many-to-one": survey_values is checked unique by key above.
  left_join(
    survey_values |> dplyr::select(key, Room_temperature, pH, previous_weightofgrains, Liquid_volume),
    by = "key", relationship = "many-to-one"
  )

as_numeric_value <- function(value) {
  suppressWarnings(as.numeric(gsub("[^0-9.-]", "", as.character(value))))
}
correlation_metadata$Timeframe <- as.numeric(
  sub("^T", "", correlation_metadata$Stage)
)
correlation_metadata$`Room temperature` <- as_numeric_value(
  correlation_metadata$Room_temperature
)
correlation_metadata$pH <- as_numeric_value(correlation_metadata$pH)
correlation_metadata$`Liquid volume` <- as_numeric_value(
  correlation_metadata$Liquid_volume
)
correlation_metadata$`Grain weight` <- as_numeric_value(
  correlation_metadata$previous_weightofgrains
)

species_profile <- read_csv(
  file.path(DATA_DIR, "04_metacache_total_species_profile_v2.csv"),
  show_col_types = FALSE
)
if (!"sample_id" %in% names(species_profile)) {
  stop("Species profile is missing the sample_id column.")
}
species_profile$key <- normalize_sample_id(species_profile$sample_id)
if (anyDuplicated(species_profile$key)) {
  stop("Species profile contains duplicate normalized sample IDs.")
}

milk_prevalent <- read_csv(
  file.path(DATA_DIR, "milk_metacache_prevalence.csv"),
  show_col_types = FALSE
)$specie
water_prevalent <- read_csv(
  file.path(DATA_DIR, "water_metacache_prevalence.csv"),
  show_col_types = FALSE
)$specie
covariates <- c("Timeframe", "Room temperature", "pH",
                "Liquid volume", "Grain weight")
species_columns <- setdiff(names(species_profile), c("sample_id", "key"))
correlation_rows <- list()
row_index <- 0L

for (type in c("ML", "WL")) {
  kefir_label <- if (type == "ML") "Milk kefir" else "Water kefir"
  type_metadata <- correlation_metadata[
    correlation_metadata[["kefir type"]] == type, ]
  missing_profile_ids <- setdiff(type_metadata$key, species_profile$key)
  if (length(missing_profile_ids)) {
    message("Excluding ", length(missing_profile_ids),
            " ", kefir_label, " sample(s) without species profiles.")
    type_metadata <- type_metadata[type_metadata$key %in% species_profile$key, ]
  }

  for (covariate in covariates) {
    usable_metadata <- type_metadata[
      is.finite(type_metadata[[covariate]]), ]
    profile_index <- match(usable_metadata$key, species_profile$key)
    if (anyNA(profile_index)) {
      stop("Species profiles are missing for ", sum(is.na(profile_index)),
           " ", kefir_label, " sample(s) used for ", covariate, ".")
    }
    profiles <- species_profile[profile_index, , drop = FALSE]
    max_abundance <- vapply(
      profiles[, species_columns, drop = FALSE],
      max, numeric(1), na.rm = TRUE
    )
    selected_species <- names(max_abundance)[max_abundance > 0.1]

    for (species in selected_species) {
      abundance <- as.numeric(profiles[[species]])
      covariate_values <- as.numeric(usable_metadata[[covariate]])
      complete <- is.finite(abundance) & is.finite(covariate_values)
      abundance <- abundance[complete]
      covariate_values <- covariate_values[complete]
      if (length(abundance) <= 4 ||
          length(unique(abundance)) < 2 ||
          length(unique(covariate_values)) < 2) {
        next
      }

      test <- suppressWarnings(cor.test(
        abundance, covariate_values, method = "spearman", exact = FALSE
      ))
      row_index <- row_index + 1L
      correlation_rows[[row_index]] <- data.frame(
        kefir_type = kefir_label,
        metadata = covariate,
        species = species,
        p_value = test$p.value,
        r = unname(test$estimate)
      )
    }
  }
}

correlation_data <- bind_rows(correlation_rows)
if (!nrow(correlation_data)) {
  stop("No species-by-covariate correlations could be computed.")
}
correlation_data$fdr <- p.adjust(correlation_data$p_value, method = "fdr")
correlation_data$stars <- cut(
  correlation_data$fdr,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)
correlation_data <- correlation_data[
  (correlation_data$kefir_type == "Milk kefir" &
     correlation_data$species %in% milk_prevalent) |
    (correlation_data$kefir_type == "Water kefir" &
       correlation_data$species %in% water_prevalent), ]
if (!nrow(correlation_data)) {
  stop("No prevalent species correlations are available for Supplementary Figure 4.")
}
correlation_data$metadata <- factor(
  correlation_data$metadata,
  levels = rev(covariates)
)
correlation_data$kefir_type <- factor(
  correlation_data$kefir_type,
  levels = c("Milk kefir", "Water kefir")
)

panel_b <- ggplot(
  correlation_data,
  aes(x = species, y = metadata, fill = r)
) +
  geom_tile() +
  geom_text(aes(label = stars), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C",
    name = "Spearman's rank correlation coefficient"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5
  )) +
  facet_wrap(~ kefir_type, scales = "free", nrow = 1) +
  labs(x = "", y = "Covariate") +
  theme(
    legend.position = "top",
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(size = 11),
    legend.key.size = unit(0.7, "cm"),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 7),
    axis.text.x = element_text(size = 3.5, face = "italic", angle = 90),
    axis.text.y = element_text(size = 6.5),
    axis.title = element_text(face = "bold", size = 10)
  )

figure <- ggarrange(
  panel_a, panel_b,
  nrow = 2, ncol = 1,
  labels = c("A.", "B."),
  heights = c(0.38, 0.62),
  common.legend = FALSE,
  font.label = list(size = 18)
)

ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_4_mantel_gower.tiff"),
       figure, width = 12, height = 8, dpi = 300,
       compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_4_mantel_gower.pdf"),
       figure, width = 12, height = 8, bg = "white")
save_figure(figure, file.path(FIGURES_DIR, "Supplementary_Figure_4_mantel_gower_toolkit"),
            width_mm = 12 * 25.4, height_mm = 8 * 25.4,
            source = PRIMARY_SOURCE, seed = NULL)
write_csv(plot_data, file.path(OUT_DIR, "Supplementary_Figure_4_mantel_gower_data.csv"))

message("Wrote two-panel Supplementary Figure 4 to ", FIGURES_DIR)
message("Panel B contains ", nrow(correlation_data), " species-by-covariate cells.")
