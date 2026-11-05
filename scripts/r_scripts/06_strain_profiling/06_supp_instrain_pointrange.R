# Ported from a pre-repo personal-workstation script (instrain_metrics_genome_point_range_r2_mean_etc.R,
# 2024-09-06) that referenced an external Q:/H2020 Master/... network drive. All inputs below
# already exist in this repo's data/ folder under their current names; only the plotting logic
# (and the path bootstrap, figure styling, and output naming) were adapted.
#
# Produces: Supplementary Figure 5 -- "The point range plot illustrates the distribution of
# various InStrain metrics, single nucleotide substitutions (SNS) count, and single nucleotide
# variants (SNV) count and r2 mean - across microbial species in Milk kefir and Water kefir
# samples" (Kefir4All_Manuscript_TRACKED_SB.docx).

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style()
DATA_DIR <- here::here("data")
OUT_DIR <- here::here("output_reviewed", "06_strain_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "Supplementary_Figure_5")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

METRIC_COLUMNS <- c("r2_mean", "SNS_count", "SNV_count")
N_SECONDARY_CLUSTERS <- 18L

instrain <- read_csv(file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"), show_col_types = FALSE)
stopifnot(
  "instrain data is missing expected columns" =
    all(c("classification", "category.y", "cluster", "sample_id", METRIC_COLUMNS) %in% names(instrain))
)

# ---- Prevalent-species lists (same files/convention used elsewhere in this pipeline) ----
milk_prevalence <- read_csv(file.path(DATA_DIR, "milk_taxonomic_profile_prevalence.csv"), show_col_types = FALSE)
water_prevalence <- read_csv(file.path(DATA_DIR, "water_taxonomic_profile_prevalence.csv"), show_col_types = FALSE)

# WHY: both files' literal header "...2" is renamed by readr's name-repair to "...5"
# (its own column position), not kept as "...2" -- confirmed empirically for this file.
stopifnot("expected prevalence column '...5' not found -- readr's name repair may have changed" =
            "...5" %in% names(milk_prevalence) && "...5" %in% names(water_prevalence))
total_prevalence <- bind_rows(milk_prevalence, water_prevalence) |>
  rename(species_clean = `...5`) |>
  mutate(species_clean = gsub("_", " ", species_clean))

milk_species <- total_prevalence$species_clean[total_prevalence$kefir_type == "milk"]
water_species <- total_prevalence$species_clean[total_prevalence$kefir_type == "water"]

# ---- For each metric x species x kefir type x secondary cluster: min/max/mean/n ----
metric_data <- instrain |>
  pivot_longer(cols = all_of(METRIC_COLUMNS), names_to = "unit", values_to = "value") |>
  filter(!is.na(value), category.y %in% c("Milk.kefir", "Water.kefir")) |>
  group_by(unit, species = classification, kefir_type = category.y, cluster) |>
  summarise(min = min(value), max = max(value), mean = mean(value),
            n = n_distinct(sample_id), .groups = "drop")

# ---- Restrict to the same prevalent-species list used throughout the rest of the pipeline ----
plot_data <- bind_rows(
  metric_data |> filter(kefir_type == "Milk.kefir", species %in% milk_species),
  metric_data |> filter(kefir_type == "Water.kefir", species %in% water_species)
) |>
  mutate(cluster = factor(gsub(".*_", "", cluster), levels = as.character(seq_len(N_SECONDARY_CLUSTERS))),
         unit = gsub("_", " ", unit),
         kefir_type = gsub("\\.", " ", kefir_type))
stopifnot("no rows left after restricting to prevalent species -- check prevalence files and ID formats" =
            nrow(plot_data) > 0)
if (anyNA(plot_data$cluster)) {
  stop("cluster labels outside 1:", N_SECONDARY_CLUSTERS,
       " found after parsing -- update N_SECONDARY_CLUSTERS or check the cluster column format")
}

p <- ggplot(plot_data, aes(x = species, y = mean, ymin = min, ymax = max, colour = cluster)) +
  geom_pointrange() +
  facet_wrap(~ kefir_type + unit, scales = "free") +
  theme_bw() +
  xlab("Species") + ylab("Value") +
  labs(colour = "Secondary\ncluster") +
  theme(strip.text = element_text(size = 22, face = "bold", colour = "black"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "italic"),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "right")

ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_5_instrain_metrics.tiff"), p,
       width = 20, height = 12, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_5_instrain_metrics.pdf"), p,
       width = 20, height = 12, bg = "white")

# Named "_toolkit", not matching the ggsave() names above, because those already wrote
# the .pdf at this path -- save_figure() would otherwise refuse to overwrite it.
save_figure(p, file.path(FIGURES_DIR, "Supplementary_Figure_5_instrain_metrics_toolkit"),
            width_mm = 20 * 25.4, height_mm = 12 * 25.4,
            source = file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"),
            seed = NULL)

write_csv(plot_data, file.path(OUT_DIR, "Supplementary_Figure_5_instrain_metrics_data.csv"))

message("Wrote Supplementary_Figure_5_instrain_metrics.{tiff,pdf} to ", FIGURES_DIR)
