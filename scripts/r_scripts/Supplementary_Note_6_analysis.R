# Supplementary_Note_6_analysis.R ---------------------------------------
# Per-species co-detection counts of within-species secondary clusters
# across the Kefir4All metagenomes, backing the per-pair frequencies in
# Supplementary Note 6. For every prevalent multi-cluster species, every
# pair of secondary clusters is enumerated, and the number of metagenomes
# in which both members of the pair are detected by inStrain is counted.
#
# Reproduces the n-values reported in Supplementary Note 6 (and provides
# the full pair table per species, which the Note summarises selectively).
#
# Inputs (in flat data/):
#   instrain_genome_species_primary_data_v4.csv
#
# Output:
#   figures/Supplementary_Note_6_data.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
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

instrain <- read_csv(file.path(data_dir, "instrain_genome_species_primary_data_v4.csv"),
                     show_col_types = FALSE)

# Apply the same per-genome detection filter stated in the Methods
instrain <- instrain %>%
  filter(!is.na(breadth), !is.na(breadth_expected),
         breadth >= 0.35, (breadth / breadth_expected) >= 0.75)

# For each species, enumerate every pair of secondary clusters and count
# the number of metagenomes in which both are detected.
species_list <- instrain %>% pull(classification) %>% unique() %>% sort()

results <- lapply(species_list, function(sp) {
  sub <- instrain %>% filter(classification == sp)
  if (length(unique(sub$cluster)) < 2) return(NULL)
  pres <- sub %>% distinct(sample_id, cluster) %>% arrange(sample_id, cluster)
  per_sample <- pres %>% group_by(sample_id) %>%
                summarise(clusters = list(sort(unique(cluster))), .groups = "drop") %>%
                filter(lengths(clusters) >= 2)
  if (nrow(per_sample) == 0) return(NULL)
  pair_rows <- do.call(rbind, lapply(per_sample$clusters, function(cs) {
    if (length(cs) < 2) return(NULL)
    p <- t(combn(cs, 2)); data.frame(c1 = p[, 1], c2 = p[, 2])
  }))
  pair_rows %>% group_by(c1, c2) %>%
    summarise(n_co_detected = n(), .groups = "drop") %>%
    arrange(desc(n_co_detected)) %>%
    mutate(species = sp, .before = 1)
})

results_df <- bind_rows(results) %>%
  filter(n_co_detected >= 5)   # robustness: matches the >=5 threshold used in the Note

write_tsv(results_df, file.path(out_dir, "Supplementary_Note_6_data.tsv"))
message("Wrote Supplementary_Note_6_data.tsv with ", nrow(results_df),
        " species-pair co-detection rows (>=5 metagenomes).")
