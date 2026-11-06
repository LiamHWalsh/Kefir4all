# ==============================================================================
# community_turnover_mixed_v1.R -- mixed-effects test for within-participant
# taxonomic turnover across consecutive visit intervals (Figure 3 statistics)
#
# WHY this script exists: the manuscript Results/Methods describe a segmented
# linear model with an a priori breakpoint at week 9 (two slopes, a test for
# whether they differ). That model lives only in community_changes_v2.R, an
# older script no longer linked to Figure 3 -- the current Figure 3 script
# (community_changes_v3.R) computes no such test ("Descriptive only -- no
# significance test is computed in this script"). Week 9 was also not found
# in the data; it was an a priori guess, not an inflection point the turnover
# curve actually shows.
#
# This script instead tests the claim the manuscript text actually makes
# ("turnover was highest in the first interval and lower thereafter") without
# assuming where any change happens: distance ~ interval, interval as a plain
# 5-level factor (T1->T2 .. T5->T6), no breakpoint. A planned contrast then
# asks whether the first interval differs from the mean of the other four.
#
# Reuses find_root(), normalize_id(), load_metadata(), load_taxonomy(),
# prepare_panel() and build_intervals() from 04_community_stability_v3.R
# verbatim, sourced with community_stability_v3.run = FALSE so sourcing it
# does not also run its own (different) analysis. Only MetaCache taxonomic
# profiles are used -- SUPER-FOCUS and volatiles are out of scope for this
# specific test (the manuscript sentence concerns community composition).
#
# Primary model per sample type (MG/WG/ML/WL):
#   m1: distance ~ interval + (1 | participant), ML (for the LRT)
#   m0: distance ~ 1        + (1 | participant), ML
#   Overall interval effect = anova(m0, m1), a likelihood-ratio test, df = 4.
#   Planned contrast (REML refit, emmeans): mean of the four later intervals
#   minus the first interval, weights c(-1, 0.25, 0.25, 0.25, 0.25) against
#   interval levels in order T1, T2, T3, T4, T5.
#   BH correction is applied across the four sample types' contrast p-values.
#
# Sensitivity: a Bayesian version of the same model (brms) would report the
# posterior probability of a decrease per type, as a robustness check -- not
# the primary analysis, and only run if brms is installed. The project's slim
# Docker image does not include a Stan backend, so this step is skipped with
# an explicit message rather than silently omitted; it is not required to
# answer the manuscript's question, which the primary model already answers.
#
# EXPECTED RESULTS (from an independent Python check on the same data/rule,
# for validation -- not a target to match, a discrepancy >~0.01 should be
# investigated, not silently accepted):
#   Contrast (later - first): WG -0.176 (-0.285, -0.067) p_BH=0.003
#                              WL -0.247 (-0.369, -0.126) p_BH=0.0003
#                              MG -0.088 (-0.194,  0.018) p_BH=0.14
#                              ML -0.145 (-0.342,  0.053) p_BH=0.15
#   Overall LRT p: MG 0.45, WG 0.029, ML 0.37, WL 0.0009
#   Interval means (T1->T2 .. T5->T6):
#     MG 0.27 0.23 0.17 0.17 0.20 | WG 0.43 0.23 0.28 0.27 0.28
#     ML 0.71 0.59 0.49 0.53 0.63 | WL 0.68 0.43 0.54 0.43 0.26
#   Total intervals across all four types: 251, from 40 distinct participants.
# ==============================================================================

options(community_stability_v3.run = FALSE)

# Minimal, duplicated root-finder: needed to locate the sibling script to
# source below, before any of ITS helpers (including its own find_root()) are
# available. Kept identical to 04_community_stability_v3.R's own find_root()
# so both scripts agree on the project root from the same invocation.
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

ROOT <- find_root()
source(file.path(ROOT, "scripts", "r_scripts", "04_taxonomic_profiling", "04_community_stability_v3.R"))

if (!requireNamespace("lme4", quietly = TRUE)) stop("Install lme4 for the mixed-effects model")
if (!requireNamespace("lmerTest", quietly = TRUE)) stop("Install lmerTest for the mixed-effects model")
if (!requireNamespace("emmeans", quietly = TRUE)) stop("Install emmeans for the planned contrast")
library(lmerTest)  # masks stats::lmer with lme4's; loads lme4

TURNOVER_SEED <- 20260916L
INTERVAL_LEVELS <- paste0("T", 1:5)  # the FROM stage of each of the 5 consecutive intervals
# Later-minus-first contrast weights, in INTERVAL_LEVELS order (T1 is the
# first interval T1->T2; T2..T5 are the four later intervals, equally weighted).
CONTRAST_WEIGHTS <- c(-1, 0.25, 0.25, 0.25, 0.25)

#' Build the per-type turnover table from build_intervals()'s output, keeping
#' only participants with >=2 intervals within this type (an interval-level
#' mixed model needs within-participant replication across at least two
#' interval categories to estimate the random intercept meaningfully).
turnover_table <- function(bc, md) {
  iv <- build_intervals(bc, md)  # from 04_community_stability_v3.R; columns include from/to/distance
  if (!nrow(iv)) return(iv)
  iv$interval <- factor(iv$from, levels = INTERVAL_LEVELS)
  n_per_participant <- table(iv$participant)
  keep <- iv$participant %in% names(n_per_participant)[n_per_participant >= 2]
  iv[keep, , drop = FALSE]
}

#' Fit the primary mixed model for one sample type's turnover table.
#' Returns the LRT, the planned contrast, and descriptive interval means.
fit_turnover_model <- function(d) {
  d$participant <- factor(d$participant)
  means <- stats::aggregate(distance ~ interval, data = d, FUN = mean)
  means <- means[match(INTERVAL_LEVELS, means$interval), ]
  m1_ml <- lmerTest::lmer(distance ~ interval + (1 | participant), data = d, REML = FALSE)
  m0_ml <- lmerTest::lmer(distance ~ 1 + (1 | participant), data = d, REML = FALSE)
  lrt <- anova(m0_ml, m1_ml)
  lrt_p <- lrt[["Pr(>Chisq)"]][2]
  m1_reml <- lmerTest::lmer(distance ~ interval + (1 | participant), data = d, REML = TRUE)
  emm <- emmeans::emmeans(m1_reml, "interval")
  contrast_spec <- stats::setNames(list(CONTRAST_WEIGHTS), "later_vs_first")
  ctr <- as.data.frame(emmeans::contrast(emm, contrast_spec, infer = c(TRUE, TRUE)))
  list(n_intervals = nrow(d), n_participants = length(unique(d$participant)),
       lrt_chisq = unname(lrt[["Chisq"]][2]), lrt_p = lrt_p,
       contrast_estimate = ctr$estimate[1], ci_low = ctr$lower.CL[1], ci_high = ctr$upper.CL[1],
       contrast_p = ctr$p.value[1], interval_means = setNames(means$distance, means$interval))
}

#' Bayesian sensitivity check (brms): same model structure, vague-but-centred
#' prior on the interval means. Only run if brms is installed -- the project's
#' slim Docker image does not include a Stan backend, so this is expected to
#' be skipped there; that is reported explicitly, not treated as a pass.
fit_bayesian_sensitivity <- function(d) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    message("brms not installed -- skipping Bayesian sensitivity check (primary model above is unaffected)")
    return(NULL)
  }
  d$participant <- factor(d$participant)
  fit <- brms::brm(distance ~ 0 + interval + (1 | participant), data = d,
                    prior = brms::prior(normal(0.5, 1), class = "b"),
                    chains = 4, iter = 4000, seed = TURNOVER_SEED, refresh = 0, silent = 2)
  draws <- brms::as_draws_df(fit)
  coefs <- paste0("b_interval", INTERVAL_LEVELS)
  if (!all(coefs %in% names(draws))) stop("Unexpected brms coefficient names: ", paste(names(draws), collapse = ", "))
  later <- rowMeans(draws[, coefs[-1]])
  first <- draws[[coefs[1]]]
  delta <- later - first
  list(estimate = mean(delta), ci_low = unname(stats::quantile(delta, 0.025)),
       ci_high = unname(stats::quantile(delta, 0.975)), p_drop = mean(delta < 0))
}

run_turnover_mixed <- function() {
  root <- ROOT
  parent <- file.path(root, "output", "04_taxonomic_profiling", "community_turnover_mixed_v1")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  message("Output: ", out)
  set.seed(TURNOVER_SEED)

  md <- load_metadata(root)
  taxonomy <- load_taxonomy(root, md)

  all_intervals <- list()
  model_rows <- list()
  bayes_rows <- list()
  for (type in names(TYPE_LABELS)) {
    pp <- prepare_panel(taxonomy$mat, md, type, relative = TRUE)
    d <- turnover_table(pp$bc, pp$md)
    if (!nrow(d)) { message(type, ": no eligible intervals, skipping"); next }
    all_intervals[[type]] <- data.frame(type = type, d[, c("participant", "from", "to", "interval", "distance")])
    fit <- fit_turnover_model(d)
    model_rows[[type]] <- data.frame(type = type, n_intervals = fit$n_intervals,
      n_participants = fit$n_participants, lrt_chisq = fit$lrt_chisq, lrt_p = fit$lrt_p,
      contrast_estimate = fit$contrast_estimate, ci_low = fit$ci_low, ci_high = fit$ci_high,
      contrast_p = fit$contrast_p)
    bayes <- tryCatch(fit_bayesian_sensitivity(d), error = function(e) {
      message(type, ": Bayesian sensitivity check failed: ", conditionMessage(e)); NULL
    })
    if (!is.null(bayes)) bayes_rows[[type]] <- data.frame(type = type, estimate = bayes$estimate,
      ci_low = bayes$ci_low, ci_high = bayes$ci_high, p_drop = bayes$p_drop)
  }

  intervals <- do.call(rbind, all_intervals)
  stopifnot("expected 251 total turnover intervals across MG/WG/ML/WL -- check prepare_panel/build_intervals filtering upstream" =
              nrow(intervals) == 251)
  stopifnot("expected 40 distinct participants across MG/WG/ML/WL -- check prepare_panel/build_intervals filtering upstream" =
              length(unique(intervals$participant)) == 40)
  write.csv(intervals, file.path(out, "intervals.csv"), row.names = FALSE)

  interval_means <- do.call(rbind, lapply(names(TYPE_LABELS), function(type) {
    d <- all_intervals[[type]]
    if (is.null(d)) return(NULL)
    means <- stats::aggregate(distance ~ interval, data = d, FUN = mean)
    data.frame(type = type, means[match(INTERVAL_LEVELS, means$interval), ])
  }))
  write.csv(interval_means, file.path(out, "interval_means.csv"), row.names = FALSE)

  tests <- do.call(rbind, model_rows)
  tests$p_BH <- p.adjust(tests$contrast_p, "BH")
  write.csv(tests, file.path(out, "mixed_model_tests.csv"), row.names = FALSE)

  if (length(bayes_rows)) {
    write.csv(do.call(rbind, bayes_rows), file.path(out, "bayesian_sensitivity.csv"), row.names = FALSE)
  }

  writeLines(c(
    "Tests whether within-participant taxonomic (MetaCache) turnover differs by",
    "consecutive-visit interval, per sample type (MG/WG/ML/WL), without assuming",
    "where any change happens (no breakpoint).",
    "Primary model: distance ~ interval + (1 | participant), interval a 5-level",
    "factor (T1->T2 .. T5->T6). Overall effect: likelihood-ratio test (ML fit),",
    "df = 4. Planned contrast (REML refit, emmeans): mean of the four later",
    "intervals minus the first interval.",
    "BH correction applied across the four sample types' contrast p-values.",
    "Participants need >=2 intervals within a sample type to be included.",
    "Bayesian sensitivity (brms), when available: same model structure, vague",
    "prior centred at 0.5 on the interval means; reports the posterior",
    "probability that turnover decreased (later < first). Not the primary",
    "analysis -- a robustness check only.",
    paste("Run directory:", out)), file.path(out, "METHODS.txt"))
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

  print(tests, row.names = FALSE)
  invisible(list(output = out, intervals = intervals, tests = tests))
}

if (isTRUE(getOption("community_turnover_mixed_v1.run", TRUE))) run_turnover_mixed()
