# Supplementary_Fig_9.R --------------------------------------------------
# dRep secondary-cluster assignments across the Kefir4All study, by
# species and timepoint. Per-prevalent-species alluvial flows; one panel
# per species, with stacks at each Stage (T0..T6 mapped to wk0..wk21)
# coloured by dRep secondary cluster.
#
# Mirrors the Figure 8 layout but uses dRep secondary clusters rather
# than inStrain within-species clusters; the two together correspond to
# the two halves of the legacy Figure 10.
#
# Inputs (all in flat data/):
#   Cdb.csv                         (dRep cluster assignments)
#   mag_metadata_v3.csv.gz          (MAG -> species, sample, kefir type)
#
# Outputs:
#   figures/Supplementary_Fig_9.png
#   figures/Supplementary_Fig_9.pdf
#   figures/Supplementary_Fig_9_data.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
repo_root <- here::here()
data_dir  <- file.path(repo_root, "data")
out_dir   <- file.path(repo_root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cdb <- read_csv(file.path(data_dir, "Cdb.csv"), show_col_types = FALSE)
md  <- read_csv(file.path(data_dir, "mag_metadata_v3.csv.gz"),
                show_col_types = FALSE, guess_max = 50000)

strip_fa <- function(x) sub("\\.fa$", "", as.character(x))
cdb <- cdb %>% mutate(g = strip_fa(genome))
md  <- md  %>% mutate(g = strip_fa(user_genome))

joined <- cdb %>%
  inner_join(md %>%
               select(g, classification, Sample, data_source, `kefir type`),
             by = "g") %>%
  mutate(Stage = str_match(genome, "_T(\\d+)_")[, 2]) %>%
  filter(!is.na(Stage), !is.na(classification)) %>%
  mutate(Stage = paste0("T", Stage)) %>%
  filter(data_source == "This study")

species_short <- function(s) {
  m <- regmatches(s, regexpr("s__[^;]+$", s))
  ifelse(length(m), trimws(sub("^s__", "", m)), s)
}
joined <- joined %>%
  mutate(species = vapply(classification, species_short, character(1)),
         category.y = case_when(
           `kefir type` %in% c("ML", "MG") ~ "Milk.kefir",
           `kefir type` %in% c("WL", "WG") ~ "Water.kefir",
           TRUE ~ NA_character_
         )) %>%
  filter(!is.na(category.y))

stage_levels <- c("T0", "T1", "T2", "T3", "T4", "T5", "T6")
stage_labels <- c("T0", "wk01", "wk05", "wk09", "wk13", "wk17", "wk21")
joined <- joined %>%
  filter(Stage %in% stage_levels) %>%
  mutate(Stage = factor(Stage, levels = stage_levels, labels = stage_labels))

# Restrict to prevalent species (>= 10 MAGs in the kefir type)
prev <- joined %>% count(category.y, species) %>% filter(n >= 10) %>% select(category.y, species)
joined_p <- joined %>% inner_join(prev, by = c("category.y", "species"))

joined_p <- joined_p %>% mutate(clusters = gsub(".*_", "", as.character(secondary_cluster)))

plot_df <- as.data.frame(xtabs(~ Stage + clusters + species + category.y, data = joined_p)) %>%
  filter(Freq != 0) %>%
  mutate(category.y = gsub("\\.", " ", as.character(category.y)),
         clusters   = factor(clusters, levels = sort(unique(clusters), method = "radix")))

p <- ggplot(plot_df,
            aes(x = Stage, stratum = clusters, alluvium = clusters,
                y = Freq, fill = clusters)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", colour = "darkgray") +
  geom_stratum() +
  facet_wrap(~ category.y + species, scales = "free_y") +
  labs(title = "",
       x = "Timeframe",
       y = "Number of MAGs assigned",
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

n_panels <- nrow(unique(plot_df[, c("category.y", "species")]))
fig_w <- 14
fig_h <- max(7, 2.0 * ceiling(n_panels / 4))

ggsave(file.path(out_dir, "Supplementary_Fig_9.png"), p,
       width = fig_w, height = fig_h, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Fig_9.pdf"), p,
       width = fig_w, height = fig_h, bg = "white")

write_tsv(plot_df, file.path(out_dir, "Supplementary_Fig_9_data.tsv"))
message("Wrote Supplementary_Fig_9.{png,pdf} + Supplementary_Fig_9_data.tsv to ", out_dir)
