# ==============================================================================
# superfocus_diversity_bayesian_v1.R -- Bayesian test for whether SUPER-FOCUS
# (broad functional subsystem) Shannon diversity differs between
# household-fermented and laboratory-controlled samples
#
# No existing script in this repo runs this comparison on SUPER-FOCUS data --
# 05_resistome.R tests it for antimicrobial-resistance genes specifically
# (a narrower functional profile), and this script asks the same question of
# SUPER-FOCUS's broader subsystem-level functional annotation, using the
# same model shape as resistome_diversity_bayesian_v1.R and
# taxonomic_diversity_bayesian_v1.R (Figure 5 Panel B's taxonomic version).
#
# data/functional/superfocus_p3.csv's "WK0..."-prefixed sample_ids use dashes
# (e.g. "WK069-2-L-48h-05D2_S81"), while global_water_kefir_metadata_v1.csv's
# merge_column uses underscores for the same samples -- normalised here the
# same way 05_resistome.R normalises this file's "WK0..." IDs before merging.
# Without this, those rows (the Mortensen et al 2023 water-kefir comparison
# cohort) silently fail to match and WL ends up with no laboratory-controlled
# group at all.
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
  library(vegan)
})

ROOT <- find_root()
SUPERFOCUS_SEED <- 20261107L

#' Build the per-sample SUPER-FOCUS Shannon diversity table, restricted to
#' ML/WL samples with a known kefir type and a household/laboratory
#' condition. Reports the row count at each filtering step.
load_superfocus_diversity <- function(root) {
  global_mk_metadata <- read_csv(file.path(root, "data", "global_milk_kefir_metadata_v1.csv"), show_col_types = FALSE)
  global_wk_metadata <- read_csv(file.path(root, "data", "global_water_kefir_metadata_v1.csv"), show_col_types = FALSE)

  kefir4all_metadata <- read_csv(file.path(root, "data", "kefir4all_sample_metadata_v2.csv"), show_col_types = FALSE)
  kefir4all_metadata$merge_column <- gsub("_host_removed_R..fastq.gz", "", kefir4all_metadata$merge_column)
  kefir4all_metadata <- kefir4all_metadata[!duplicated(kefir4all_metadata$merge_column), ]

  total_metadata <- bind_rows(
    dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
    dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`),
    dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`))

  sf <- read_csv(file.path(root, "data", "functional", "superfocus_p3.csv"), show_col_types = FALSE)
  sf$sample_id[grep("WK0", sf$sample_id)] <- gsub("-", "_", sf$sample_id[grep("WK0", sf$sample_id)])
  sf <- sf[sf$sample_id %in% total_metadata$merge_column, , drop = FALSE]
  message("SUPER-FOCUS samples matching project metadata: ", length(unique(sf$sample_id)))

  sf_wide <- sf |>
    dplyr::select(sample_id, `Subsystem Level 3`, `reads_%`) |>
    dplyr::group_by(sample_id, `Subsystem Level 3`) |>
    dplyr::summarise(`reads_%` = sum(`reads_%`), .groups = "drop") |>
    pivot_wider(names_from = `Subsystem Level 3`, values_from = `reads_%`, values_fill = 0)

  shannon <- diversity(sf_wide |> tibble::column_to_rownames("sample_id"))
  diversity_df <- data.frame(site = names(shannon), shannon = as.numeric(shannon))

  n_before <- nrow(diversity_df)
  diversity_df <- merge(diversity_df, total_metadata, by.x = "site", by.y = "merge_column", all.x = TRUE)
  stopifnot("diversity_df row count changed after merging total_metadata -- check for duplicate merge_column" =
              nrow(diversity_df) == n_before)

  diversity_df <- diversity_df[!is.na(diversity_df$`kefir type`), ]
  message("SUPER-FOCUS samples with a known kefir type: ", nrow(diversity_df))

  diversity_df$conditions <- "Household conditions"
  diversity_df$conditions[diversity_df$data_source %in% c("Walsh et al 2023", "Mortensen et al 2023")] <- "Laboratory controlled"

  diversity_df <- diversity_df[diversity_df$`kefir type` %in% c("ML", "WL"), ]
  message("SUPER-FOCUS samples restricted to ML/WL (matching the resistome/taxonomic scripts' scope): ", nrow(diversity_df))

  diversity_df
}

#' Bayesian equivalent of an independent-samples diversity-by-condition test
#' for one kefir type -- same model shape as resistome_diversity_bayesian_v1.R
#' and taxonomic_diversity_bayesian_v1.R. No manually-specified prior, for
#' the same reason given there: Shannon diversity has no natural bound to
#' centre one on.
fit_bayesian_diversity <- function(d) {
  d$conditions <- factor(d$conditions, levels = c("Laboratory controlled", "Household conditions"))
  fit <- brms::brm(shannon ~ conditions, data = d,
                    chains = 4, iter = 4000, seed = SUPERFOCUS_SEED, refresh = 0, silent = 2)
  draws <- brms::as_draws_df(fit)
  delta <- draws[["b_conditionsHouseholdconditions"]]
  list(estimate = mean(delta), ci_low = unname(stats::quantile(delta, 0.025)),
       ci_high = unname(stats::quantile(delta, 0.975)),
       p_increase = mean(delta > 0), p_decrease = mean(delta < 0))
}

run_superfocus_diversity <- function() {
  root <- ROOT
  parent <- file.path(root, "output", "05_functional_profiling", "superfocus_diversity_bayesian_v1")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  message("Output: ", out)
  set.seed(SUPERFOCUS_SEED)

  d <- load_superfocus_diversity(root)
  write.csv(d, file.path(out, "superfocus_diversity_data.csv"), row.names = FALSE)

  results <- lapply(c("ML", "WL"), function(type) {
    dt <- d[d$`kefir type` == type, ]
    bayes <- fit_bayesian_diversity(dt)
    data.frame(type = type, n = nrow(dt), estimate = bayes$estimate,
               ci_low = bayes$ci_low, ci_high = bayes$ci_high,
               p_increase = bayes$p_increase, p_decrease = bayes$p_decrease)
  })
  results <- do.call(rbind, results)
  write.csv(results, file.path(out, "superfocus_diversity_bayesian_tests.csv"), row.names = FALSE)

  print(results, row.names = FALSE)
  invisible(list(output = out, data = d, results = results))
}

if (isTRUE(getOption("superfocus_diversity_bayesian_v1.run", TRUE))) run_superfocus_diversity()
