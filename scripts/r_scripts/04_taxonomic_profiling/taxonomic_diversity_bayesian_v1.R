# ==============================================================================
# taxonomic_diversity_bayesian_v1.R -- Bayesian test for whether taxonomic
# (MetaCache species) Shannon diversity differs between household-fermented
# and laboratory-controlled samples (Figure 5 Panel B statistics)
#
# Figure 5 Panel B (04_environmental_microbes_v2.R) already tests this with a
# frequentist test (ggbetweenstats, Kruskal-Wallis/pairwise), for milk liquid
# (ML) samples. This script re-derives the same per-sample Shannon diversity
# values from the same inputs and asks the same question with a Bayesian
# model instead, reporting a posterior difference and P(direction).
#
# ML only, matching Panel B's own scope: that panel's total_metadata never
# merges in the water-kefir global comparison dataset, so there is no
# laboratory-controlled water liquid (WL) group to compare against -- all WL
# samples there are household conditions. Confirmed empirically (table of
# kefir type x conditions) before writing this script, not assumed.
#
# Does not source 04_environmental_microbes_v2.R: that script also builds
# Figure 5 Panels A and C (persistence detection, plotting), which this
# script does not need. The data pipeline up to the per-sample Shannon value
# is reproduced here instead -- the same read/merge/diversity() steps that
# script itself uses for Panel B, pulled out on their own.
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
  library(vegan)
})

ROOT <- find_root()
TAXONOMY_SEED <- 20261107L

#' Rebuild Figure 5 Panel B's per-sample taxonomic Shannon diversity table,
#' restricted to ML with a known kefir type -- the same scope that
#' panel's statistical test uses. Reports the row count at each step.
load_taxonomic_diversity <- function(root) {
  mc_raw <- read_csv(file.path(root, "data", "04_metacache_total_species_profile_v2.csv"), show_col_types = FALSE)
  mc_raw$sample_id <- sub("^.*/", "", mc_raw$sample_id)
  mc_raw <- mc_raw[mc_raw$sample_id != "TG_EC_S72", ]

  kefir4all_metadata <- read_csv(file.path(root, "data", "kefir4all_sample_metadata_v2.csv"), show_col_types = FALSE) |>
    mutate(merge_column = gsub("_host_removed_R..fastq.gz", "", merge_column)) |>
    filter(!duplicated(merge_column))

  num_block <- suppressWarnings(mc_raw[, -1])
  num_block[] <- lapply(num_block, function(x) suppressWarnings(as.numeric(x)))
  num_block[is.na(num_block)] <- 0

  global_mk_metadata <- read_csv(file.path(root, "data", "global_milk_kefir_metadata_v1.csv"), show_col_types = FALSE)
  total_metadata <- bind_rows(
    dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
    dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`))
  stopifnot("total_metadata has duplicate merge_column values -- Kefir4All and global datasets overlap" =
              !anyDuplicated(total_metadata$merge_column))

  alpha_mat <- as.matrix(num_block)
  rownames(alpha_mat) <- mc_raw$sample_id
  keep <- rowSums(alpha_mat) > 0
  alpha_mat <- alpha_mat[keep, , drop = FALSE]
  message("Taxonomic samples before metadata merge: ", nrow(alpha_mat))

  shannon <- diversity(alpha_mat)
  alpha_df <- data.frame(site = names(shannon), shannon = as.numeric(shannon))
  alpha_df <- left_join(alpha_df, total_metadata, by = c("site" = "merge_column"), relationship = "many-to-one")

  # Same fallback Figure 5 Panel B itself uses: an unmatched sample's kefir type and
  # data_source are assumed ML / Walsh et al 2023 (no distinct water-kefir global
  # comparison cohort is merged in here, matching that panel's own scope).
  alpha_df$`kefir type`[is.na(alpha_df$`kefir type`)] <- "ML"
  alpha_df$data_source[is.na(alpha_df$`kefir type`)] <- "Walsh et al 2023"
  alpha_df <- alpha_df[!(alpha_df$site %in%
    total_metadata$merge_column[total_metadata$`kefir type` %in% c("Extraction control", "Medium control")]), ]
  alpha_df$data_source[is.na(alpha_df$data_source)] <- "Walsh et al 2023"
  message("Taxonomic samples with a known kefir type: ", nrow(alpha_df))

  alpha_df$conditions <- "Household conditions"
  alpha_df$conditions[alpha_df$data_source %in% c("Walsh et al 2023", "Mortensen et al 2023")] <- "Laboratory controlled"

  # WL is excluded entirely here, not just left untested: Figure 5 Panel B's own
  # total_metadata never merges in the water-kefir global comparison dataset, so
  # every WL sample it has is "Household conditions" -- no WL group comparison
  # exists in the analysis this script is reproducing.
  alpha_df <- alpha_df[alpha_df$`kefir type` == "ML", ]
  message("Taxonomic samples restricted to ML (matching Figure 5 Panel B's scope -- no WL ",
          "laboratory-controlled group exists in that panel's own data): ", nrow(alpha_df))

  alpha_df
}

#' Bayesian equivalent of Figure 5 Panel B's between-groups test for one
#' kefir type: independent samples, not a mixed-effects model (no repeated
#' measures here). No manually-specified prior: Shannon diversity has no
#' natural bound like Bray-Curtis distance does, so brms's default
#' weakly-informative priors are used as-is (same choice as
#' resistome_diversity_bayesian_v1.R, for the same reason).
fit_bayesian_diversity <- function(d) {
  d$conditions <- factor(d$conditions, levels = c("Laboratory controlled", "Household conditions"))
  fit <- brms::brm(shannon ~ conditions, data = d,
                    chains = 4, iter = 4000, seed = TAXONOMY_SEED, refresh = 0, silent = 2)
  draws <- brms::as_draws_df(fit)
  delta <- draws[["b_conditionsHouseholdconditions"]]
  list(estimate = mean(delta), ci_low = unname(stats::quantile(delta, 0.025)),
       ci_high = unname(stats::quantile(delta, 0.975)),
       p_increase = mean(delta > 0), p_decrease = mean(delta < 0))
}

run_taxonomic_diversity <- function() {
  root <- ROOT
  parent <- file.path(root, "output", "04_taxonomic_profiling", "taxonomic_diversity_bayesian_v1")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  message("Output: ", out)
  set.seed(TAXONOMY_SEED)

  d <- load_taxonomic_diversity(root)
  write.csv(d, file.path(out, "taxonomic_diversity_data.csv"), row.names = FALSE)

  bayes <- fit_bayesian_diversity(d)
  results <- data.frame(type = "ML", n = nrow(d), estimate = bayes$estimate,
                         ci_low = bayes$ci_low, ci_high = bayes$ci_high,
                         p_increase = bayes$p_increase, p_decrease = bayes$p_decrease)
  write.csv(results, file.path(out, "taxonomic_diversity_bayesian_tests.csv"), row.names = FALSE)

  print(results, row.names = FALSE)
  invisible(list(output = out, data = d, results = results))
}

if (isTRUE(getOption("taxonomic_diversity_bayesian_v1.run", TRUE))) run_taxonomic_diversity()
