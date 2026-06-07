# Supplementary_Fig_5.R ---------------------------------------------------
# Per-species proportion of within-species MAG pairs reaching >= 99% ANI,
# matching the legacy Extended Data Figure 3 style: stacked horizontal
# bars showing the percentage of "highly related" (>=99% ANI) versus
# "less related" pairwise comparisons within each species, faceted by
# kefir type, with labels for both segments.
#
# Inputs (all in flat data/):
#   ndb.csv.gz                      (dRep pairwise ANI)
#   Cdb.csv                         (dRep cluster assignments)
#   mag_metadata_v3.csv.gz          (MAG -> species, source, kefir type)
#
# Outputs:
#   figures/Supplementary_Fig_5.png
#   figures/Supplementary_Fig_5.pdf
#   figures/Supplementary_Fig_5_data.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  m <- grep(file_arg, args)
  if (length(m) > 0) return(dirname(normalizePath(sub(file_arg, "", args[m]))))
  for (fr in sys.frames()) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
  getwd()
}
repo_root <- normalizePath(file.path(get_script_dir(), "..", ".."))
data_dir  <- file.path(repo_root, "data")
out_dir   <- file.path(repo_root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ndb <- read_csv(file.path(data_dir, "ndb.csv.gz"), show_col_types = FALSE)
cdb <- read_csv(file.path(data_dir, "Cdb.csv"),    show_col_types = FALSE)
md  <- read_csv(file.path(data_dir, "mag_metadata_v3.csv.gz"),
                show_col_types = FALSE, guess_max = 50000)

strip_fa <- function(x) sub("\\.fa$", "", as.character(x))
ndb <- ndb %>% mutate(q = strip_fa(querry), r = strip_fa(reference))
cdb <- cdb %>% mutate(g = strip_fa(genome))
md  <- md  %>% mutate(g = strip_fa(user_genome))

g2pc <- setNames(cdb$primary_cluster, cdb$g)
g2sp <- setNames(md$classification,    md$g)
g2kt <- setNames(md$`kefir type`,      md$g)
this_study <- if ("data_source" %in% names(md)) md$g[md$data_source == "This study"] else md$g

within <- ndb %>%
  mutate(q_pc = g2pc[q], r_pc = g2pc[r],
         q_sp = g2sp[q], q_kt = g2kt[q]) %>%
  filter(!is.na(q_pc), q_pc == r_pc, q != r,
         q %in% this_study, r %in% this_study)

within <- within %>% mutate(type = case_when(
  q_kt %in% c("ML", "MG") ~ "Milk kefir",
  q_kt %in% c("WL", "WG") ~ "Water kefir",
  TRUE ~ NA_character_
)) %>% filter(!is.na(type))

species_short <- function(s) {
  m <- regmatches(s, regexpr("s__[^;]+$", s))
  ifelse(length(m), trimws(sub("^s__", "", m)), s)
}
within <- within %>% mutate(species = vapply(q_sp, species_short, character(1)))

agg <- within %>%
  group_by(type, species) %>%
  summarise(
    n_pairs   = n(),
    n_ge99    = sum(ani >= 0.99),
    pct_ge99  = 100 * n_ge99 / n_pairs,
    .groups   = "drop"
  ) %>%
  filter(n_pairs >= 3) %>%
  arrange(type, desc(pct_ge99))

write_tsv(agg, file.path(out_dir, "Supplementary_Fig_5_data.tsv"))

stacked <- bind_rows(
  agg %>% mutate(strain_type = "% highly related strains (>=99% ANI)",
                 frequency   = pct_ge99),
  agg %>% mutate(strain_type = "% less related strains (<99% ANI)",
                 frequency   = 100 - pct_ge99)
)

# Order species within each kefir type by descending pct_ge99
species_order <- agg %>% group_by(type) %>%
  arrange(type, desc(pct_ge99)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>% select(type, species, rank)

stacked <- stacked %>% left_join(species_order, by = c("type", "species")) %>%
  mutate(species = factor(species, levels = unique(species[order(type, rank)])))

p <- ggplot(stacked,
            aes(x = species, y = frequency,
                fill = factor(strain_type,
                              levels = c("% highly related strains (>=99% ANI)",
                                         "% less related strains (<99% ANI)")))) +
  geom_col(colour = "white", width = 0.85) +
  geom_text(aes(label = round(frequency, 1)),
            position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold", colour = "black", show.legend = FALSE) +
  facet_wrap(~ type, scales = "free_y") +
  scale_fill_manual(values = c(
    "% highly related strains (>=99% ANI)" = "#4C9F70",
    "% less related strains (<99% ANI)"    = "#D9D9D9"
  )) +
  labs(title = "",
       x = NULL,
       y = "Pairwise comparisons (%)",
       fill = "Strain comparisons") +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9, face = "italic"),
    axis.title = element_text(size = 10),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 10)
  ) +
  guides(fill = guide_legend(nrow = 1))

n_sp <- max(table(agg$type))
fig_w <- 12
fig_h <- max(7, 0.32 * n_sp + 2)

ggsave(file.path(out_dir, "Supplementary_Fig_5.png"), p,
       width = fig_w, height = fig_h, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Fig_5.pdf"), p,
       width = fig_w, height = fig_h, bg = "white")
message("Wrote Supplementary_Fig_5.{png,pdf} + Supplementary_Fig_5_data.tsv to ", out_dir)
