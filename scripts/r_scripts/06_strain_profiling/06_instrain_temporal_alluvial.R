# Anchor to the project root before here() loads. RStudio Server resets the session cwd
# to $HOME after sourcing .Rprofile, and here() computes its root AT PACKAGE LOAD
# TIME, so force the project root before any here load. This block is inert on non-portable
for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {



  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}

# 06_instrain_temporal_alluvial.R -----------------------------------------
# Pipeline step 6 — Strain profiling
# Results sections: 3.7 (Figure 8, revised manuscript)
#
# inStrain temporal alluvial: per prevalent species, the count of
# inStrain-defined within-species cluster detections at each Stage (T0..T6,
# mapped to wk0..wk21), with alluvial flows connecting clusters across
# successive timepoints.
#
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
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style(base_size = 8)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
})

data_dir <- here::here("data")
OUT_DIR  <- here::here("output_reviewed", "06_strain_profiling")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIGURES_DIR <- here::here("figures", "Figure_8")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

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
  milk_prev  |> mutate(kefir_type = "milk"),
  water_prev |> mutate(kefir_type = "water")
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
# popANI_reference > 0.98 is the threshold for a high-confidence within-species
# cluster assignment.
instrain_cs <- instrain |>
  filter(popANI_reference > 0.98,
         sample_id %in% kefir4all_md$merge_column)
stopifnot("no inStrain detections pass the popANI_reference > 0.98 / Kefir4All-cohort filter" =
            nrow(instrain_cs) > 0)

# Restrict to species that are prevalent in the relevant kefir type.
instrain_prevelant_cs_clust <- bind_rows(
  instrain_cs |>
    filter(`category.y` == "Milk.kefir",
           classification %in% milk_species),
  instrain_cs |>
    filter(`category.y` == "Water.kefir",
           classification %in% water_species)
)

# Extract the cluster index (e.g. "Lla.lactis_5" -> "5") and map T-stages to
# study weeks for axis labels.
instrain_prevelant_cs_clust <- instrain_prevelant_cs_clust |>
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
                               data = instrain_prevelant_cs_clust)) |>
  filter(Freq != 0) |>
  mutate(category.y = gsub("\\.", " ", as.character(category.y)),
         clusters = factor(clusters, levels = as.character(1:20)))
stopifnot("no rows left for the alluvial plot -- check prevalence/species-recode filtering" =
            nrow(plot_df) > 0)

# WHY: species names never repeat across the two kefir types (confirmed -- no
# classification value appears under both category.y groups), so faceting by
# classification alone cannot create any ambiguity between them.
# HOW: explicitly ordering the factor levels (Milk kefir species first, Water
# kefir second, alphabetical within each) keeps species grouped by kefir type
# in the panel layout, without needing "Milk kefir"/"Water kefir" repeated as
# a second strip-text line on every panel.
species_order <- plot_df |>
  dplyr::distinct(category.y, classification) |>
  dplyr::arrange(category.y != "Milk kefir", classification) |>
  dplyr::pull(classification)
plot_df$classification <- factor(plot_df$classification, levels = species_order)

p <- ggplot(plot_df,
            aes(x = Stage, stratum = clusters, alluvium = clusters,
                y = Freq, fill = clusters)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", colour = "darkgray") +
  geom_stratum() +
  facet_wrap(~ classification, scales = "free_y") +
  labs(title = "",
       x = "Timeframe",
       y = "Number of strains detected",
       fill = "Cluster") +
  kefir_theme(base_size = 8) +
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
    strip.text = element_text(size = 10, face = "bold.italic")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

n_panels <- nrow(unique(plot_df[, c("category.y", "classification")]))
fig_w <- 16
fig_h <- 10.5

ggsave(file.path(FIGURES_DIR, "Figure_8.tiff"), p,
       width = fig_w, height = fig_h, dpi = 300, bg = "white", compression = "lzw")
ggsave(file.path(FIGURES_DIR, "Figure_8.pdf"), p,
       width = fig_w, height = fig_h, bg = "white")

# Named "_toolkit", not "Figure_8", because the ggsave() calls above already wrote
# Figure_8.tiff/.pdf -- save_figure() would otherwise refuse to overwrite that .pdf.
save_figure(p, file.path(FIGURES_DIR, "Figure_8_toolkit"),
            width_mm = fig_w * 25.4, height_mm = fig_h * 25.4,
            source = file.path(data_dir, "instrain_genome_species_primary_data_v4.csv"),
            seed = NULL)

write_tsv(plot_df, file.path(OUT_DIR, "Figure_8_data.tsv"))
message("Wrote Figure_8.{tiff,pdf} + Figure_8_data.tsv to ", OUT_DIR)