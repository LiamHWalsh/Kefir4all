# ==============================================================================
# resistome_diversity_bayesian_v1.R -- Bayesian test for whether resistome
# (antimicrobial-resistance gene) Shannon diversity differs between
# household-fermented and laboratory-controlled samples
#
# 05_resistome.R already tests this with a frequentist Wilcoxon rank-sum
# test, separately for milk liquid (ML) and water liquid (WL) samples. This
# script re-derives the same per-sample Shannon diversity values from the
# same inputs and asks the same question with a Bayesian model instead, so
# the result can be reported as a posterior difference and P(direction)
# rather than (or alongside) a p-value.
#
# Does not source 05_resistome.R: that script is a flat, ~3000-line file
# with no functions and no run-guard, so sourcing it would also execute
# everything else in it (its own plots, additional analyses). The data
# pipeline up to the per-sample Shannon value is reproduced here instead --
# it is the same merge/normalise/diversity() steps 05_resistome.R itself
# uses, just pulled out on their own.
#
# Requires a Stan backend (cmdstanr + compiled CmdStan) and brms; the
# project's slim Docker image does not include one -- build and use
# kefir4all:4.4.2-brms instead (docker/Dockerfile.brms in this repo).
# ==============================================================================

find_root <- function() {
  args <- commandArgs(FALSE)
  script <- sub("^--file=", "", args[grepl("^--file=", args)])
  starts <- c(Sys.getenv("KEFIR4ALL_ROOT"),
              if (requireNamespace("here", quietly = TRUE))
                suppressMessages(tryCatch(here::here(), error = function(e) NA_character_)) else NA_character_,
              if (length(script)) dirname(script[1]) else NA_character_, getwd())
  for (start in starts[nzchar(starts) & !is.na(starts)]) {
    for (p in c(start, file.path(start, "kefir4all_review_history"))) {
      if (file.exists(file.path(p, "data", "kefir4all_sample_metadata_v2.csv")))
        return(normalizePath(p, winslash = "/"))
    }
  }
  stop("Cannot locate data root. Set KEFIR4ALL_ROOT to the review-history root.")
}

if (!requireNamespace("brms", quietly = TRUE))
  stop("Install brms (and a Stan backend, e.g. cmdstanr) to run this script -- ",
       "see docker/Dockerfile.brms, which builds kefir4all:4.4.2-brms with both.")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(vegan)
})

ROOT <- find_root()
DATA_DIR    <- file.path(ROOT, "data")
PRIVATE_DIR <- file.path(ROOT, "data", "private")
CS_METADATA_PRIVATE <- file.path(PRIVATE_DIR, "Citizen Scientist metadata_v8.csv")
if (!file.exists(CS_METADATA_PRIVATE)) {
  stop("Private metadata not found: ", CS_METADATA_PRIVATE,
       "\nCopy to data/private/ to run this script.")
}

RESISTOME_SEED <- 20261107L

#' Rebuild the per-sample resistome Shannon diversity table that
#' 05_resistome.R computes, restricted to ML/WL (liquid) samples with a known
#' kefir type -- the same scope its Wilcoxon test uses.
#' Reports the row count at each filtering step.
load_resistome_diversity <- function(root) {
  global_mk_metadata <- read_csv(file.path(root, "data", "global_milk_kefir_metadata_v1.csv"), show_col_types = FALSE)
  global_wk_metadata <- read_csv(file.path(root, "data", "global_water_kefir_metadata_v1.csv"), show_col_types = FALSE)
  global_mk_metadata$Stage <- NA
  global_wk_metadata$Stage <- NA

  kefir4all_metadata <- read_csv(file.path(root, "data", "kefir4all_sample_metadata_v2.csv"), show_col_types = FALSE)
  kefir4all_metadata$merge_column <- gsub("_host_removed_R..fastq.gz", "", kefir4all_metadata$merge_column)
  kefir4all_metadata <- kefir4all_metadata[!duplicated(kefir4all_metadata$merge_column), ]

  total_metadata <- bind_rows(
    dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`, Stage, Sample),
    dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`, Stage, Sample),
    dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`, Stage, Sample))

  gene_total_resistome <- read_csv(file.path(root, "data", "gene_total_resistome.csv"), show_col_types = FALSE)
  reads_per_sample <- read_csv(file.path(root, "data", "reads_per_sample_interleaved_total_data.csv"), show_col_types = FALSE)
  reads_per_sample$sample_id[grep("WK0", reads_per_sample$sample_id)] <-
    gsub("-", "_", reads_per_sample$sample_id[grep("WK0", reads_per_sample$sample_id)])
  reads_per_sample$X1 <- NULL

  gene_total_resistome$Sample[grep("WK0", gene_total_resistome$Sample)] <-
    gsub("-", "_", gene_total_resistome$Sample[grep("WK0", gene_total_resistome$Sample)])
  gene_total_resistome$Sample <- gsub("_antimicrobial_resistance_mapped_and_unmapped", "", gene_total_resistome$Sample)

  n_before <- nrow(gene_total_resistome)
  gene_total_resistome <- merge(gene_total_resistome, reads_per_sample,
                                 by.x = "Sample", by.y = "sample_id", all.x = TRUE)
  stopifnot("gene_total_resistome row count changed after merging reads_per_sample -- check for duplicate sample_id" =
              nrow(gene_total_resistome) == n_before)

  # CPM-style normalisation, same as 05_resistome.R: hits per million reads (X2).
  gene_total_resistome$nor_hits <- gene_total_resistome$Hits * 1 / gene_total_resistome$X2 * 10^6

  gene_total_resistome_wide <- gene_total_resistome |>
    dplyr::select(nor_hits, Gene, Sample) |>
    pivot_wider(names_from = "Gene", values_from = "nor_hits", values_fill = 0)

  message("Resistome samples before metadata merge: ", nrow(gene_total_resistome_wide))

  shannon <- diversity(gene_total_resistome_wide |> column_to_rownames("Sample"))
  diversity_df <- data.frame(site = names(shannon), shannon = as.numeric(shannon))

  n_before <- nrow(diversity_df)
  diversity_df <- merge(diversity_df, total_metadata, by.x = "site", by.y = "merge_column", all.x = TRUE)
  stopifnot("diversity_df row count changed after merging total_metadata -- check for duplicate merge_column" =
              nrow(diversity_df) == n_before)

  diversity_df <- diversity_df[!is.na(diversity_df$`kefir type`), ]
  message("Resistome samples with a known kefir type: ", nrow(diversity_df))

  diversity_df$conditions <- "Household conditions"
  diversity_df$conditions[diversity_df$data_source %in% c("Walsh et al 2023", "Mortensen et al 2023")] <- "Laboratory controlled"

  diversity_df <- diversity_df[diversity_df$`kefir type` %in% c("ML", "WL"), ]
  message("Resistome samples restricted to ML/WL (matching 05_resistome.R's Wilcoxon test scope): ", nrow(diversity_df))

  diversity_df
}

#' Bayesian equivalent of wilcox_test(shannon ~ conditions, ...) for one
#' kefir type: an independent-samples comparison, not a mixed-effects model
#' (no repeated measures here, unlike the Figure 3 turnover analysis) --
#' brms is used for the posterior difference and P(direction) it gives,
#' not because the model itself needs to be hierarchical.
#' No manually-specified prior: unlike Bray-Curtis distance (bounded to
#' [0, 1], where community_turnover_mixed_v1.R uses a prior centred at 0.5),
#' Shannon diversity here has no natural bound to centre a prior on, so
#' brms's default weakly-informative priors are used as-is.
fit_bayesian_diversity <- function(d) {
  d$conditions <- factor(d$conditions, levels = c("Laboratory controlled", "Household conditions"))
  fit <- brms::brm(shannon ~ conditions, data = d,
                    chains = 4, iter = 4000, seed = RESISTOME_SEED, refresh = 0, silent = 2)
  draws <- brms::as_draws_df(fit)
  delta <- draws[["b_conditionsHouseholdconditions"]]
  list(estimate = mean(delta), ci_low = unname(stats::quantile(delta, 0.025)),
       ci_high = unname(stats::quantile(delta, 0.975)),
       p_increase = mean(delta > 0), p_decrease = mean(delta < 0))
}

run_resistome_diversity <- function() {
  root <- ROOT
  parent <- file.path(root, "output", "05_functional_profiling", "resistome_diversity_bayesian_v1")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  message("Output: ", out)
  set.seed(RESISTOME_SEED)

  d <- load_resistome_diversity(root)
  write.csv(d, file.path(out, "resistome_diversity_data.csv"), row.names = FALSE)

  results <- lapply(c("ML", "WL"), function(type) {
    dt <- d[d$`kefir type` == type, ]
    bayes <- fit_bayesian_diversity(dt)
    data.frame(type = type, n = nrow(dt), estimate = bayes$estimate,
               ci_low = bayes$ci_low, ci_high = bayes$ci_high,
               p_increase = bayes$p_increase, p_decrease = bayes$p_decrease)
  })
  results <- do.call(rbind, results)
  write.csv(results, file.path(out, "resistome_diversity_bayesian_tests.csv"), row.names = FALSE)

  print(results, row.names = FALSE)
  invisible(list(output = out, data = d, results = results))
}

if (isTRUE(getOption("resistome_diversity_bayesian_v1.run", TRUE))) run_resistome_diversity()
