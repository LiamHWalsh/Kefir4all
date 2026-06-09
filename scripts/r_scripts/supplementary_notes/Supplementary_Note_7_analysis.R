# Supplementary_Note_7_analysis.R ---------------------------------------
# Per-species median per-genome SNV counts and nucleotide diversity from
# inStrain, backing the prose in Supplementary Note 7.
#
# Reproduces the species-level summary statistics quoted in Note 7
# (e.g. milk-kefir SNV-count medians ranging from <200 to >4,000;
# per-genome nucleotide diversity in the 10^-3 range).
#
# Inputs (in flat data/):
#   instrain_genome_species_primary_data_v4.csv
#
# Output:
#   figures/supplementary_notes/Supplementary_Note_7_data.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
repo_root <- here::here()
data_dir  <- file.path(repo_root, "data")
OUT_DIR   <- file.path(repo_root, "output", "supplementary_notes")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

instrain <- read_csv(file.path(data_dir, "instrain_genome_species_primary_data_v4.csv"),
                     show_col_types = FALSE)

# Apply the per-genome detection filter stated in the Methods
instrain <- instrain %>%
  filter(!is.na(breadth), !is.na(breadth_expected),
         !is.na(SNV_count), !is.na(nucl_diversity),
         breadth >= 0.35, (breadth / breadth_expected) >= 0.75)

# Tag kefir type from the embedded `kefir type.x` column
kefir_group <- function(k) {
  if (is.na(k)) return(NA_character_)
  if (k %in% c("ML", "MG")) "Milk kefir"
  else if (k %in% c("WL", "WG")) "Water kefir"
  else NA_character_
}

instrain <- instrain %>%
  mutate(kefir_type = vapply(`kefir type.x`, kefir_group, character(1))) %>%
  filter(!is.na(kefir_type))

# Restrict to species detected in at least 10 metagenomes within their kefir type
prev <- instrain %>%
  count(kefir_type, classification) %>%
  filter(n >= 10) %>% select(kefir_type, classification)

instrain_p <- instrain %>% inner_join(prev, by = c("kefir_type", "classification"))

summ <- instrain_p %>%
  group_by(kefir_type, classification) %>%
  summarise(
    n_detections    = n(),
    median_snv      = median(SNV_count, na.rm = TRUE),
    p25_snv         = quantile(SNV_count, 0.25, na.rm = TRUE),
    p75_snv         = quantile(SNV_count, 0.75, na.rm = TRUE),
    median_nucl_div = median(nucl_diversity, na.rm = TRUE),
    p25_nucl_div    = quantile(nucl_diversity, 0.25, na.rm = TRUE),
    p75_nucl_div    = quantile(nucl_diversity, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(kefir_type, desc(median_snv))

write_tsv(summ, file.path(OUT_DIR, "Supplementary_Note_7_data.tsv"))
message("Wrote Supplementary_Note_7_data.tsv with ", nrow(summ),
        " species-level rows (>=10 detections per kefir type).")
