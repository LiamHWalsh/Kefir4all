# Supplementary_Fig_5.R ---------------------------------------------------
# Per-species proportion of within-species MAG pairs reaching >= 99% ANI,
# computed from the dRep ANI table on within-primary-cluster pairs in the
# Kefir4All cohort.
#
# Inputs:
#   data/ndb.csv.gz                (dRep pairwise ANI)
#   data/Cdb.csv                    (dRep cluster assignments)
#   data/mag_metadata_v3.csv.gz    (MAG -> species, source)
#
# Outputs:
#   figures/Supplementary_Fig_5.png
#   figures/Supplementary_Fig_5.pdf
#   figures/Supplementary_Fig_5_data.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
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
ndb_path  <- file.path(repo_root, "data", "ndb.csv.gz")
cdb_path  <- file.path(repo_root, "data",  "Cdb.csv")
md_path   <- file.path(repo_root, "data", "mag_metadata_v3.csv.gz")
out_dir   <- file.path(repo_root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ndb <- read_csv(ndb_path, show_col_types = FALSE)
cdb <- read_csv(cdb_path, show_col_types = FALSE)
md  <- read_csv(md_path,  show_col_types = FALSE, guess_max = 50000)

strip_fa <- function(x) sub("\\.fa$", "", as.character(x))
ndb <- ndb %>% mutate(q = strip_fa(querry), r = strip_fa(reference))
cdb <- cdb %>% mutate(g = strip_fa(genome))
md  <- md  %>% mutate(g = strip_fa(user_genome))

g2pc <- setNames(cdb$primary_cluster, cdb$g)
g2sp <- setNames(md$classification, md$g)
this_study <- if ("data_source" %in% names(md)) {
  md$g[md$data_source == "This study"]
} else md$g

within <- ndb %>%
  mutate(q_pc = g2pc[q], r_pc = g2pc[r], q_sp = g2sp[q]) %>%
  filter(!is.na(q_pc), q_pc == r_pc, q != r) %>%
  filter(q %in% this_study, r %in% this_study)

species_short <- function(s) {
  m <- regmatches(s, regexpr("s__[^;]+$", s))
  ifelse(length(m), trimws(sub("^s__", "", m)), s)
}
within <- within %>%
  mutate(species = vapply(q_sp, species_short, character(1)))

agg <- within %>%
  group_by(species) %>%
  summarise(
    n_pairs = n(),
    n_ge99  = sum(ani >= 0.99),
    pct_ge99 = 100 * n_ge99 / n_pairs,
    .groups = "drop"
  ) %>%
  filter(n_pairs >= 3) %>%
  arrange(desc(pct_ge99))

write_tsv(agg, file.path(out_dir, "Supplementary_Fig_5_data.tsv"))

top <- agg %>% slice_head(n = 30) %>% arrange(pct_ge99) %>%
  mutate(species = factor(species, levels = species))

p <- ggplot(top, aes(x = pct_ge99, y = species)) +
  geom_col(fill = "#4C9F70", colour = "white", width = 0.85) +
  geom_text(aes(label = paste0("n=", n_pairs)),
            hjust = -0.05, size = 2.5) +
  coord_cartesian(xlim = c(0, 105)) +
  labs(
    title = "Supplementary Fig. 5. Per-species proportion of within-species MAG pairs reaching >=99% ANI.",
    x = "Pairwise comparisons within species >= 99% ANI (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, hjust = 0),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(out_dir, "Supplementary_Fig_5.png"), p,
       width = 8, height = 9, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Fig_5.pdf"), p,
       width = 8, height = 9, bg = "white")
message("Wrote Supplementary_Fig_5 to ", out_dir)
