# 06_instrain_temporal_alluvial.R -----------------------------------------
# Pipeline step 6 — Strain profiling
# Results sections: 3.7 (Figure 8, revised manuscript)
#
# inStrain temporal alluvial: per prevalent species, the count of
# inStrain-defined within-species cluster detections at each Stage (T0..T6,
# mapped to wk0..wk21), with alluvial flows connecting clusters across
# successive timepoints.
#
# This corresponds to panel B of the legacy Figure 10 (combined dRep +
# inStrain), now used as the single-panel revised Figure 8. The dRep
# equivalent is in 06_strain_profiling/06_supp_strain_profiling.R.
#
# Inputs (all in flat data/):
#   instrain_genome_species_primary_data_v4.csv
#   kefir4all_sample_metadata_v2.csv
#   milk_taxonomic_profile_prevalence.csv
#   water_taxonomic_profile_prevalence.csv
#
# Outputs:
#   figures/Figure_8.png
#   figures/Figure_8.pdf
#   figures/Figure_8_data.tsv

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

data_dir <- here::here("data")
OUT_DIR  <- here::here("output", "06_strain_profiling")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

instrain <- read_csv(file.path(data_dir, "instrain_genome_species_primary_data_v4.csv"),
                     show_col_types = FALSE)
kefir4all_md <- read_csv(file.path(data_dir, "kefir4all_sample_metadata_v2.csv"),
                         show_col_types = FALSE)
kefir4all_md$merge_column <- gsub("_host_removed_R..fastq.gz", "",
                                  kefir4all_md$merge_column)
kefir4all_md <- kefir4all_md[!duplicated(kefir4all_md$merge_column), ]

milk_prev  <- read_csv(file.path(data_dir, "milk_taxonomic_profile_prevalence.csv"),
                       show_col_types = FALSE)
water_prev <- read_csv(file.path(data_dir, "water_taxonomic_profile_prevalence.csv"),
                       show_col_types = FALSE)
total_prev <- bind_rows(
  milk_prev  %>% mutate(kefir_type = "milk"),
  water_prev %>% mutate(kefir_type = "water")
)

# Normalise a few historical species labels to current taxonomy
sp_col <- names(total_prev)[5]
recode_sp <- c(
  "Lactobacillus_ghanensis"             = "Liquorilactobacillus ghanensis",
  "Lactobacillus_satsumensis"           = "Liquorilactobacillus satsumensis",
  "Lactococcus_lactis subcluster 1"     = "Lactococcus lactis",
  "Lactococcus_lactis subcluster 2"     = "Lactococcus cremoris",
  "Pseudomonas_fragi_subspecies 1"      = "Pseudomonas fragi",
  "Zymomonas_mobilis_subcluster 1"      = "Zymomonas mobilis",
  "Lactobacillus_perolens"              = "Schleiferilactobacillus perolens"
)
total_prev[[sp_col]] <- ifelse(total_prev[[sp_col]] %in% names(recode_sp),
                               recode_sp[total_prev[[sp_col]]],
                               total_prev[[sp_col]])
total_prev[[sp_col]] <- gsub("_", " ", total_prev[[sp_col]])

milk_species  <- total_prev[[sp_col]][total_prev$kefir_type == "milk"]
water_species <- total_prev[[sp_col]][total_prev$kefir_type == "water"]

# Filter inStrain to high-confidence detections in the Kefir4All cohort.
# popANI_reference > 0.98 is the threshold used in the original Figure 10.R.
instrain_cs <- instrain %>%
  filter(popANI_reference > 0.98,
         sample_id %in% kefir4all_md$merge_column)

# Restrict to species that are prevalent in the relevant kefir type.
instrain_prevelant_cs_clust <- bind_rows(
  instrain_cs %>%
    filter(`category.y` == "Milk.kefir",
           classification %in% milk_species),
  instrain_cs %>%
    filter(`category.y` == "Water.kefir",
           classification %in% water_species)
)

# Extract the cluster index (e.g. "Lla.lactis_5" -> "5") and map T-stages to
# study weeks for axis labels.
instrain_prevelant_cs_clust <- instrain_prevelant_cs_clust %>%
  mutate(clusters = gsub(".*_", "", cluster),
         Stage = recode(Stage,
                        T1 = "wk01", T2 = "wk05", T3 = "wk09",
                        T4 = "wk13", T5 = "wk17", T6 = "wk21",
                        T0 = "T0"),
         Stage = factor(Stage,
                        levels = c("T0", "wk01", "wk05", "wk09",
                                   "wk13", "wk17", "wk21")))

# Build the count table for the alluvial plot.
plot_df <- as.data.frame(xtabs(~ Stage + clusters + classification + category.y,
                               data = instrain_prevelant_cs_clust)) %>%
  filter(Freq != 0) %>%
  mutate(category.y = gsub("\\.", " ", as.character(category.y)),
         clusters = factor(clusters, levels = as.character(1:20)))

p <- ggplot(plot_df,
            aes(x = Stage, stratum = clusters, alluvium = clusters,
                y = Freq, fill = clusters)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", colour = "darkgray") +
  geom_stratum() +
  facet_wrap(~ category.y + classification, scales = "free_y") +
  labs(title = "",
       x = "Timeframe",
       y = "Number of strains detected",
       fill = "Cluster") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 8, hjust = 0.5),
    axis.title = element_text(size = 10),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 8, face = "italic")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

n_panels <- nrow(unique(plot_df[, c("category.y", "classification")]))
fig_w <- 14
fig_h <- max(7, 2.0 * ceiling(n_panels / 4))

ggsave(file.path(OUT_DIR, "Figure_8.png"), p,
       width = fig_w, height = fig_h, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Figure_8.pdf"), p,
       width = fig_w, height = fig_h, bg = "white")

write_tsv(plot_df, file.path(OUT_DIR, "Figure_8_data.tsv"))
message("Wrote Figure_8.{png,pdf} + Figure_8_data.tsv to ", OUT_DIR)
