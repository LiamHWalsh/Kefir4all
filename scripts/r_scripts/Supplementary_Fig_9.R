# Supplementary_Fig_9.R --------------------------------------------------
# dRep secondary-cluster assignments across the Kefir4All study, by
# species and timepoint. Per-prevalent-species small multiples; stacked
# bar at each Stage (T0..T6, mapped to wk0..wk21) coloured by dRep
# secondary cluster.
#
# Inputs:
#   data/Figure_8_data/Cdb.csv                  (dRep cluster assignments)
#   data/Figure_10_data/mag_metadata_v3.csv.gz  (MAG -> species, sample)
#
# Outputs:
#   figures/Supplementary_Fig_9.png
#   figures/Supplementary_Fig_9.pdf
#   figures/Supplementary_Fig_9_data.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."))
cdb_path  <- file.path(repo_root, "data", "Figure_8_data",  "Cdb.csv")
md_path   <- file.path(repo_root, "data", "Figure_10_data", "mag_metadata_v3.csv.gz")
out_dir   <- file.path(repo_root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cdb <- read_csv(cdb_path, show_col_types = FALSE)
md  <- read_csv(md_path,  show_col_types = FALSE, guess_max = 50000)

strip_fa <- function(x) sub("\\.fa$", "", as.character(x))
cdb <- cdb %>% mutate(g = strip_fa(genome))
md  <- md  %>% mutate(g = strip_fa(user_genome))

m <- cdb %>%
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
m <- m %>% mutate(species = vapply(classification, species_short, character(1)))

m <- m %>% mutate(kg = case_when(
  `kefir type` %in% c("ML", "MG") ~ "Milk kefir",
  `kefir type` %in% c("WL", "WG") ~ "Water kefir",
  TRUE ~ NA_character_
)) %>% filter(!is.na(kg))

stage_levels <- c("T0", "T1", "T2", "T3", "T4", "T5", "T6")
stage_labels <- c("T0", "wk01", "wk05", "wk09", "wk13", "wk17", "wk21")
m <- m %>%
  filter(Stage %in% stage_levels) %>%
  mutate(Stage = factor(Stage, levels = stage_levels, labels = stage_labels))

prev <- m %>% count(species) %>% filter(n >= 10) %>% pull(species)
m_p <- m %>% filter(species %in% prev) %>%
  mutate(species_label = paste0(species, "  (", tolower(sub(" kefir", "", kg)), " kefir)"),
         species_label = factor(species_label,
                                levels = unique(species_label[order(kg, species)])))

p <- ggplot(m_p, aes(x = Stage, fill = factor(secondary_cluster))) +
  geom_bar(position = "stack", width = 0.9, colour = "white", linewidth = 0.15) +
  facet_wrap(~ species_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = unname(grDevices::hcl.colors(20, "Set 3")), guide = "none") +
  labs(
    title = "Supplementary Fig. 9. dRep secondary-cluster assignments across the Kefir4All study, by species and timepoint.",
    x = NULL, y = "MAGs"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, hjust = 0),
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid.minor = element_blank()
  )

n_panels <- length(unique(m_p$species_label))
fig_w <- 12
fig_h <- max(6, 2.0 * ceiling(n_panels / 4))

ggsave(file.path(out_dir, "Supplementary_Fig_9.png"), p,
       width = fig_w, height = fig_h, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Fig_9.pdf"), p,
       width = fig_w, height = fig_h, bg = "white")

m_p %>%
  transmute(species,
            kefir_type = kg,
            stage = as.character(Stage),
            secondary_cluster, primary_cluster, Sample) %>%
  write_tsv(file.path(out_dir, "Supplementary_Fig_9_data.tsv"))

message("Wrote Supplementary_Fig_9 to ", out_dir)
