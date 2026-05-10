# Figure_8.R --------------------------------------------------------------
# Revised-manuscript Figure 8.
#
# Single-panel inStrain temporal small-multiples: per prevalent species, a
# stacked bar at each Stage (T0..T6, mapped to wk0..wk21) coloured by
# inStrain-defined within-species cluster.
#
# Inputs:
#   data/instrain_genome_species_primary_data_v4.csv
#
# Outputs:
#   figures/Figure_8.png
#   figures/Figure_8.pdf
#   figures/Figure_8_data.tsv
#
# Filters mirror the Methods: per-genome breadth >= 0.35 and breadth /
# breadth_expected >= 0.75, restricted to Stage T0..T6 and the Kefir4All
# (this study) data source.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."))
data_path <- file.path(repo_root, "data",
                      "instrain_genome_species_primary_data_v4.csv")
out_dir   <- file.path(repo_root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(data_path, show_col_types = FALSE) %>%
  filter(!is.na(breadth), !is.na(breadth_expected), !is.na(Stage)) %>%
  filter(breadth >= 0.35, (breadth / breadth_expected) >= 0.75) %>%
  filter(str_detect(as.character(Stage), "^T[0-9]$"))

if ("data_source.x" %in% names(df)) {
  df <- df %>% filter(`data_source.x` == "This study")
}

df <- df %>%
  mutate(kg = case_when(
    str_detect(genome, "_M[GL]_") ~ "Milk kefir",
    str_detect(genome, "_W[GL]_") ~ "Water kefir",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(kg))

stage_levels <- c("T0", "T1", "T2", "T3", "T4", "T5", "T6")
stage_labels <- c("T0", "wk01", "wk05", "wk09", "wk13", "wk17", "wk21")

# Prevalent species (>= 10 detections within their kefir type)
prev_sp <- df %>%
  count(kg, classification) %>%
  filter(n >= 10) %>%
  arrange(kg, desc(n))

df_plot <- df %>%
  filter(classification %in% prev_sp$classification) %>%
  mutate(
    Stage = factor(Stage, levels = stage_levels, labels = stage_labels),
    species_label = paste0(classification, "  (", tolower(sub(" kefir", "", kg)), " kefir)"),
    species_label = factor(species_label,
                           levels = unique(species_label[order(kg, classification)]))
  )

p <- ggplot(df_plot, aes(x = Stage, fill = factor(cluster))) +
  geom_bar(position = "stack", width = 0.9, colour = "white", linewidth = 0.15) +
  facet_wrap(~ species_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = unname(grDevices::hcl.colors(20, "Set 3")), guide = "none") +
  labs(
    title = "Figure 8. inStrain-defined within-species cluster detections across the Kefir4All study, by species and timepoint.",
    x = NULL, y = "metagenomes"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, hjust = 0),
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid.minor = element_blank()
  )

n_panels <- length(unique(df_plot$species_label))
fig_w <- 12
fig_h <- max(6, 2.0 * ceiling(n_panels / 4))

ggsave(file.path(out_dir, "Figure_8.png"), p,
       width = fig_w, height = fig_h, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "Figure_8.pdf"), p,
       width = fig_w, height = fig_h, bg = "white")

df_plot %>%
  transmute(species = classification,
            kefir_type = kg,
            stage = as.character(Stage),
            cluster = cluster,
            sample_id = sample_id) %>%
  write_tsv(file.path(out_dir, "Figure_8_data.tsv"))

message("Wrote Figure_8.png + .pdf + Figure_8_data.tsv to ", out_dir)
