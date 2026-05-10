# Supplementary_Fig_8.R --------------------------------------------------
# Within-sample polymorphism (per-genome SNV count and nucleotide diversity)
# across prevalent kefir species, indexed by inStrain. Four panels:
#   (A) milk kefir  - per-genome SNV sites (log y-axis)
#   (B) water kefir - per-genome SNV sites (log y-axis)
#   (C) milk kefir  - nucleotide diversity
#   (D) water kefir - nucleotide diversity
#
# Inputs:
#   data/instrain_genome_species_primary_data_v4.csv
#
# Outputs:
#   figures/Supplementary_Fig_8.png
#   figures/Supplementary_Fig_8.pdf
#   figures/Supplementary_Fig_8_data.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."))
data_path <- file.path(repo_root, "data",
                      "instrain_genome_species_primary_data_v4.csv")
out_dir   <- file.path(repo_root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(data_path, show_col_types = FALSE) %>%
  filter(!is.na(breadth), !is.na(breadth_expected),
         !is.na(SNV_count), !is.na(nucl_diversity)) %>%
  filter(breadth >= 0.35, (breadth / breadth_expected) >= 0.75)

milk  <- df %>% filter(`kefir type.x` %in% c("ML", "MG"))
water <- df %>% filter(`kefir type.x` %in% c("WL", "WG"))

prev_species <- function(d, n = 10) {
  d %>% count(classification) %>% filter(n >= !!n) %>% pull(classification)
}

mp <- milk  %>% filter(classification %in% prev_species(milk))
wp <- water %>% filter(classification %in% prev_species(water))

order_by_median_snv <- function(d) {
  d %>%
    group_by(classification) %>%
    summarise(med = median(SNV_count, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>% pull(classification)
}
mp <- mp %>% mutate(classification = factor(classification,
                                            levels = order_by_median_snv(.)))
wp <- wp %>% mutate(classification = factor(classification,
                                            levels = order_by_median_snv(.)))

box_theme <- theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, hjust = 0),
    axis.text.x = element_text(angle = 70, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid.minor = element_blank()
  )
fill_col <- "#4C9F70"

p_a <- ggplot(mp, aes(x = classification, y = SNV_count)) +
  geom_boxplot(fill = fill_col, alpha = 0.7, outlier.shape = NA) +
  scale_y_log10(labels = label_comma()) +
  labs(title = "(A) Milk kefir - per-genome SNV counts",
       x = NULL, y = "SNV sites per genome (log scale)") +
  box_theme

p_b <- ggplot(wp, aes(x = classification, y = SNV_count)) +
  geom_boxplot(fill = fill_col, alpha = 0.7, outlier.shape = NA) +
  scale_y_log10(labels = label_comma()) +
  labs(title = "(B) Water kefir - per-genome SNV counts",
       x = NULL, y = "SNV sites per genome (log scale)") +
  box_theme

p_c <- ggplot(mp, aes(x = classification, y = nucl_diversity)) +
  geom_boxplot(fill = fill_col, alpha = 0.7, outlier.shape = NA) +
  labs(title = "(C) Milk kefir - nucleotide diversity",
       x = NULL, y = "Nucleotide diversity") +
  box_theme

p_d <- ggplot(wp, aes(x = classification, y = nucl_diversity)) +
  geom_boxplot(fill = fill_col, alpha = 0.7, outlier.shape = NA) +
  labs(title = "(D) Water kefir - nucleotide diversity",
       x = NULL, y = "Nucleotide diversity") +
  box_theme

pp <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(title = "Supplementary Fig. 8. Within-sample polymorphism (per-genome SNV sites and nucleotide diversity) across prevalent kefir species, indexed by inStrain.",
                  theme = theme(plot.title = element_text(size = 10)))

ggsave(file.path(out_dir, "Supplementary_Fig_8.png"), pp,
       width = 11, height = 10, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Fig_8.pdf"), pp,
       width = 11, height = 10, bg = "white")

# Per-species summary
summ <- df %>%
  filter(classification %in% c(levels(mp$classification),
                               levels(wp$classification))) %>%
  group_by(`kefir type.x`, classification) %>%
  summarise(n_detections = n(),
            median_snv = median(SNV_count),
            median_nucl_div = median(nucl_diversity),
            .groups = "drop")
write_tsv(summ, file.path(out_dir, "Supplementary_Fig_8_data.tsv"))
message("Wrote Supplementary_Fig_8 to ", out_dir)
