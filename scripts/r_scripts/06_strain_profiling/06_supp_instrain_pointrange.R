# Ported from a pre-repo personal-workstation script (instrain_metrics_genome_point_range_r2_mean_etc.R,
# 2024-09-06) that referenced an external Q:/H2020 Master/... network drive. All inputs below
# already exist in this repo's data/ folder under their current names; only the plotting logic
# (and the path bootstrap, figure styling, and output naming) were adapted.
#
# Produces: Supplementary Figure 5 -- "The point range plot illustrates the distribution of
# various InStrain metrics, single nucleotide substitutions (SNS) count, and single nucleotide
# variants (SNV) count and r2 mean - across microbial species in Milk kefir and Water kefir
# samples" (Kefir4All_Manuscript_TRACKED_SB.docx).

if (!requireNamespace("here",   quietly = TRUE)) install.packages("here")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
kefir_use_style()
DATA_DIR    <- here::here("data")
OUT_DIR     <- here::here("output", "06_strain_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "06_strain_profiling")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

library(readr)
library(dplyr)
library(ggplot2)

instrain <- read_csv(file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"), show_col_types = FALSE)

# ---- Prevalent-species lists (same files/convention used elsewhere in this pipeline) ----
milk_prevalence  <- read_csv(file.path(DATA_DIR, "milk_taxonomic_profile_prevalence.csv"),  show_col_types = FALSE)
water_prevalence <- read_csv(file.path(DATA_DIR, "water_taxonomic_profile_prevalence.csv"), show_col_types = FALSE)
total_prevalence <- rbind(milk_prevalence, water_prevalence)
# readr's name-repair renames the literal header "...2" to "...5" (its own column
# position) to avoid colliding with its auto-naming scheme -- confirmed empirically.
total_prevalence$...5 <- gsub("_", " ", total_prevalence$...5)

# ---- For each metric x species x kefir type x strain cluster: min/max/mean/n ----
metric_data <- data.frame()
for (unit in c("r2_mean", "SNS_count", "SNV_count")) {
  for (species in levels(as.factor(instrain$classification))) {
    for (type in c("Milk.kefir", "Water.kefir")) {
      t1 <- instrain[which(instrain$classification == species &
                              !is.na(getElement(instrain, unit)) &
                              instrain$category.y == type), ]
      if (nrow(t1) == 0) next
      for (clust in levels(as.factor(t1$cluster))) {
        sub <- t1[which(t1$cluster == clust), ]
        metric_data <- rbind(metric_data, data.frame(
          species = species,
          cluster = clust,
          unit = as.character(unit),
          min = range(getElement(sub, unit))[1],
          max = range(getElement(sub, unit))[2],
          mean = mean(getElement(sub, unit)),
          n = as.numeric(length(unique(sub$sample_id))),
          kefir_type = type
        ))
      }
    }
  }
}

# ---- Restrict to the same prevalent-species list used throughout the rest of the pipeline ----
prevalence_metric_data <- rbind(
  metric_data[which(metric_data$kefir_type == "Milk.kefir" &
                       metric_data$species %in% total_prevalence$...5[which(total_prevalence$kefir_type == "milk")]), ],
  metric_data[which(metric_data$kefir_type == "Water.kefir" &
                       metric_data$species %in% total_prevalence$...5[which(total_prevalence$kefir_type == "water")]), ]
)

plot_data <- prevalence_metric_data %>%
  mutate(cluster = factor(gsub(".*_", "", cluster), levels = as.character(1:18)),
         unit = gsub("_", " ", unit),
         kefir_type = gsub("\\.", " ", kefir_type))

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
write_csv(prevalence_metric_data, file.path(OUT_DIR, "Supplementary_Figure_5_instrain_metrics_data.csv"))

message("Wrote Supplementary_Figure_5_instrain_metrics.{tiff,pdf} to ", FIGURES_DIR)
