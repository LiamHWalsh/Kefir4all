# Anchor to the project root before here() loads. RStudio Server resets the session cwd
# to $HOME after sourcing .Rprofile, and here() computes its root AT PACKAGE LOAD
# TIME, so force the project root before any here load. This block is inert on non-portable
# machines (none of the candidate dirs exist), so local/native R runs are unaffected.
for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {



  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}

# Extended_Data_Figure_3.R -------------------------------------------------
# Per-species proportion of within-species MAG pairs reaching >= 99% ANI:
# stacked horizontal bars showing the percentage of "highly related"
# (>=99% ANI) versus "less related" pairwise comparisons within each
# species, faceted by kefir type, with labels for both segments.
#
# Renamed from "Supplementary_Fig_5" -- that was a stale name left over from
# an earlier manuscript numbering; this is Extended Data Figure 3 in the
# current manuscript (see scripts/README.md's script-to-results table).
#
# Inputs (all in flat data/):
#   ndb.csv.gz                      (dRep pairwise ANI)
#   Cdb.csv                         (dRep cluster assignments)
#   mag_metadata_v3.csv.gz          (MAG -> species, source, kefir type)
#
# Outputs:
#   figures/Extended_Data_Figure_3/Extended_Data_Figure_3.tiff
#   figures/Extended_Data_Figure_3/Extended_Data_Figure_3.pdf
#   output_reviewed/04_taxonomic_profiling/Extended_Data_Figure_3_data.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
})

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style(base_size = 8)
repo_root <- here::here()
data_dir  <- file.path(repo_root, "data")
OUT_DIR   <- file.path(repo_root, "output_reviewed", "04_taxonomic_profiling")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIGURES_DIR <- file.path(repo_root, "figures", "Extended_Data_Figure_3")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

ndb <- read_csv(file.path(data_dir, "ndb.csv.gz"), show_col_types = FALSE)
cdb <- read_csv(file.path(data_dir, "Cdb.csv"),    show_col_types = FALSE)
md  <- read_csv(file.path(data_dir, "mag_metadata_v3.csv.gz"),
                show_col_types = FALSE, guess_max = 50000)

strip_fa <- function(x) sub("\\.fa$", "", as.character(x))
ndb <- ndb |> mutate(q = strip_fa(querry), r = strip_fa(reference))
cdb <- cdb |> mutate(g = strip_fa(genome))
md  <- md  |> mutate(g = strip_fa(user_genome))

g2pc <- setNames(cdb$primary_cluster, cdb$g)
g2sp <- setNames(md$classification,    md$g)
g2kt <- setNames(md$`kefir type`,      md$g)
this_study <- if ("data_source" %in% names(md)) md$g[md$data_source == "This study"] else md$g

within <- ndb |>
  mutate(q_pc = g2pc[q], r_pc = g2pc[r],
         q_sp = g2sp[q], q_kt = g2kt[q]) |>
  filter(!is.na(q_pc), q_pc == r_pc, q != r,
         q %in% this_study, r %in% this_study)

within <- within |> mutate(type = case_when(
  q_kt %in% c("ML", "MG") ~ "Milk kefir",
  q_kt %in% c("WL", "WG") ~ "Water kefir",
  TRUE ~ NA_character_
)) |> filter(!is.na(type))

species_short <- function(s) {
  m <- regmatches(s, regexpr("s__[^;]+$", s))
  ifelse(length(m), trimws(sub("^s__", "", m)), s)
}
within <- within |> mutate(species = vapply(q_sp, species_short, character(1)))

agg <- within |>
  group_by(type, species) |>
  summarise(
    n_pairs   = n(),
    n_ge99    = sum(ani >= 0.99),
    pct_ge99  = 100 * n_ge99 / n_pairs,
    .groups   = "drop"
  ) |>
  filter(n_pairs >= 3) |>
  arrange(type, desc(pct_ge99))
stopifnot("no species have >=3 within-cluster ANI pairs -- check ndb.csv.gz / cluster assignment" =
            nrow(agg) > 0)

write_tsv(agg, file.path(OUT_DIR, "Extended_Data_Figure_3_data.tsv"))

stacked <- bind_rows(
  agg |> mutate(strain_type = "% highly related strains (>=99% ANI)",
                 frequency   = pct_ge99),
  agg |> mutate(strain_type = "% less related strains (<99% ANI)",
                 frequency   = 100 - pct_ge99)
)

# Order species within each kefir type by descending pct_ge99
species_order <- agg |> group_by(type) |>
  arrange(type, desc(pct_ge99)) |>
  mutate(rank = row_number()) |>
  ungroup() |> select(type, species, rank)

# relationship = "many-to-one": stacked has 2 rows per (type, species) (one per strain_type
# segment); species_order has exactly 1, built from the same grouping as agg.
stacked <- stacked |> left_join(species_order, by = c("type", "species"), relationship = "many-to-one") |>
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
  kefir_theme(base_size = 8) +
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
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 1))

n_sp <- max(table(agg$type))
fig_w <- 26.2
fig_h <- 17.3

ggsave(file.path(FIGURES_DIR, "Extended_Data_Figure_3.tiff"), p,
       width = fig_w, height = fig_h, dpi = 300, bg = "white", compression = "lzw")
ggsave(file.path(FIGURES_DIR, "Extended_Data_Figure_3.pdf"), p,
       width = fig_w, height = fig_h, bg = "white")

save_figure(p, file.path(FIGURES_DIR, "Extended_Data_Figure_3_provenance"),
            width_mm = fig_w * 25.4, height_mm = fig_h * 25.4,
            source = file.path(data_dir, "ndb.csv.gz"),
            seed = NULL)

message("Wrote Extended_Data_Figure_3.{tiff,pdf} + Extended_Data_Figure_3_data.tsv to ", OUT_DIR)