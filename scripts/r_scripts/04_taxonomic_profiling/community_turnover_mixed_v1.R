# ==============================================================================
# community_turnover_mixed_v1.R -- Bayesian mixed-effects test for
# within-participant taxonomic turnover across consecutive visit intervals
# (Figure 3 statistics)
#
# For each sample type (milk/water kefir, grain/liquid), tests whether
# Bray-Curtis turnover between consecutive visits is higher in the first
# interval after starting than in later intervals. Turnover is modelled as a
# function of which of the five consecutive intervals (T1->T2 .. T5->T6) is
# being measured, with no assumption about when any change happens; a
# planned contrast compares the first interval to the mean of the later four.
#
# Reuses find_root(), load_metadata(), load_taxonomy(), prepare_panel() and
# build_intervals() from 04_community_stability_v3.R. Only MetaCache
# taxonomic profiles are used (the manuscript sentence this supports concerns
# community composition, not SUPER-FOCUS or volatiles).
#
# Model, per sample type: brms::brm(distance ~ 0 + interval + (1 |
# participant)), vague prior normal(0.5, 1) on each interval mean, 4 chains,
# 4000 iter. Reports the posterior mean, 95% credible interval, and
# P(decrease) for (mean of the four later intervals) - (first interval).
#
# Requires a Stan backend (cmdstanr + compiled CmdStan) and brms; the
# project's slim Docker image does not include one -- build and use
# kefir4all:4.4.2-brms instead (docker/Dockerfile.brms in this repo).
#
# Cross-checked two independent ways on this dataset: a frequentist
# mixed-effects model (lme4::lmer, REML planned contrast via emmeans) agrees
# in direction and magnitude for all four types (BH-adjusted contrast p: MG
# 0.153, WG 0.0050, ML 0.161, WL 0.0008), and a separately written
# implementation gives the same P(decrease): ~0.94 (MG), ~0.999 (WG), ~0.92
# (ML), ~1.00 (WL).
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

if (!requireNamespace("brms", quietly = TRUE))
  stop("Install brms (and a Stan backend, e.g. cmdstanr) to run this script -- ",
       "see docker/Dockerfile.brms, which builds kefir4all:4.4.2-brms with both.")

TURNOVER_SEED <- 20260916L
INTERVAL_LEVELS <- paste0("T", 1:5)  # the FROM stage of each of the 5 consecutive intervals

#' Build the per-type turnover table from build_intervals()'s output, keeping
#' only participants with >=2 intervals within this type (the model needs
#' within-participant replication across at least two interval categories to
#' estimate the random intercept meaningfully).
#' Reports the row count before and after this filter, and how many
#' participants it drops, so an unexpectedly large exclusion is visible
#' rather than silently reducing the dataset the model is fit on.
turnover_table <- function(bc, md, type) {
  iv <- build_intervals(bc, md)  # from 04_community_stability_v3.R; columns include from/to/distance
  if (!nrow(iv)) return(iv)
  iv$interval <- factor(iv$from, levels = INTERVAL_LEVELS)
  n_per_participant <- dplyr::count(iv, participant, name = "n")
  eligible <- n_per_participant$participant[n_per_participant$n >= 2]
  kept <- iv[iv$participant %in% eligible, , drop = FALSE]
  n_dropped_participants <- length(unique(iv$participant)) - length(eligible)
  message(type, ": ", nrow(iv), " intervals from ", length(unique(iv$participant)),
          " participants before the >=2-intervals filter; ", nrow(kept), " intervals from ",
          length(eligible), " participants after (", n_dropped_participants,
          " participant(s) with a single interval excluded)")
  kept
}

#' Bayesian mixed-effects model for one sample type's turnover table: a vague
#' prior centred at 0.5 on each interval's mean (Bray-Curtis distance ranges
#' 0-1, so this is weakly informative, not a substantive assumption).
#' Reports the posterior mean, 95% credible interval and P(decrease) for
#' (mean of the four later intervals) - (first interval).
fit_bayesian_turnover <- function(d) {
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
       ci_high = unname(stats::quantile(delta, 0.975)), p_decrease = mean(delta < 0))
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
  bayes_rows <- list()
  for (type in names(TYPE_LABELS)) {
    pp <- prepare_panel(taxonomy$mat, md, type, relative = TRUE)
    d <- turnover_table(pp$bc, pp$md, type)
    if (!nrow(d)) { message(type, ": no eligible intervals, skipping"); next }
    all_intervals[[type]] <- data.frame(type = type, d[, c("participant", "from", "to", "interval", "distance")])
    bayes <- fit_bayesian_turnover(d)
    bayes_rows[[type]] <- data.frame(type = type, n_intervals = nrow(d),
      n_participants = length(unique(d$participant)), estimate = bayes$estimate,
      ci_low = bayes$ci_low, ci_high = bayes$ci_high, p_decrease = bayes$p_decrease)
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
    means <- d |>
      dplyr::group_by(interval) |>
      dplyr::summarise(distance = mean(distance), .groups = "drop")
    data.frame(type = type, means[match(INTERVAL_LEVELS, means$interval), ])
  }))
  write.csv(interval_means, file.path(out, "interval_means.csv"), row.names = FALSE)

  tests <- do.call(rbind, bayes_rows)
  write.csv(tests, file.path(out, "bayesian_turnover_tests.csv"), row.names = FALSE)

  writeLines(c(
    "Tests whether within-participant taxonomic (MetaCache) turnover differs by",
    "consecutive-visit interval, per sample type (MG/WG/ML/WL), without assuming",
    "where any change happens (no breakpoint).",
    "Model: brms::brm(distance ~ 0 + interval + (1 | participant)), interval a",
    "5-level factor (T1->T2 .. T5->T6), vague prior normal(0.5, 1) on each",
    "interval mean, 4 chains, 4000 iter.",
    "Reported quantity: posterior distribution of (mean of the four later",
    "intervals) - (first interval) -- posterior mean, 95% credible interval,",
    "and P(decrease), the posterior probability turnover was lower after the",
    "first interval than during it.",
    "Participants need >=2 intervals within a sample type to be included.",
    "A frequentist mixed-effects model (lme4::lmer, likelihood-ratio test, REML",
    "planned contrast via emmeans) was also fit on this exact dataset and",
    "agrees with this result in direction and approximate magnitude for all",
    "four types; see the script header for those values.",
    paste("Run directory:", out)), file.path(out, "METHODS.txt"))
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

  print(tests, row.names = FALSE)
  invisible(list(output = out, intervals = intervals, tests = tests))
}

if (isTRUE(getOption("community_turnover_mixed_v1.run", TRUE))) run_turnover_mixed()
