# ==============================================================================
# functional_volatile_turnover_bayesian_v1.R -- Bayesian test for whether
# within-participant SUPER-FOCUS (functional) and volatile (metabolomic)
# turnover across consecutive visit intervals follows the same pattern
# already found for taxonomic composition in Figure 3
# (community_turnover_mixed_v1.R)
#
# Same question, same model, two more data types: does turnover (Bray-Curtis
# distance between a participant's consecutive visits) differ by WHICH of
# the five consecutive intervals (T1->T2 .. T5->T6) it is, with no assumption
# about where any change happens. community_turnover_mixed_v1.R answers this
# for MetaCache taxonomic profiles only, by design (the manuscript sentence
# it supports concerns community composition specifically). This script asks
# the identical question of SUPER-FOCUS (broad functional subsystems, all
# four sample types) and the volatile metabolomics profiles (ML/WL only, no
# grain volatile data exists), as a check on whether the same "more turnover
# early, less later" pattern (or its absence) generalises beyond taxonomy.
#
# A new script, not an edit to community_turnover_mixed_v1.R: that script
# already has a reported, cross-checked result for Figure 3 and its own
# tightly scoped purpose; extending it in place would blur that scope and
# risk the Figure 3 statistics output changing as a side effect.
#
# Reuses find_root(), load_metadata(), load_superfocus(), load_volatiles(),
# prepare_panel() and build_intervals() from 04_community_stability_v3.R,
# sourced with community_stability_v3.run = FALSE, exactly as
# community_turnover_mixed_v1.R already does for load_taxonomy(). The
# turnover_table()/fit_bayesian_turnover() model functions are duplicated
# from that script rather than shared, following the same convention it
# uses for find_root() -- each script stays self-contained and readable on
# its own.
#
# Requires a Stan backend (cmdstanr + compiled CmdStan) and brms; the
# project's slim Docker image does not include one -- build and use
# kefir4all:4.4.2-brms instead (docker/Dockerfile.brms in this repo).
# ==============================================================================

options(community_stability_v3.run = FALSE)

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
INTERVAL_LEVELS <- paste0("T", 1:5)

#' Identical to community_turnover_mixed_v1.R's own turnover_table(): keeps
#' only participants with >=2 intervals within this type, and reports the
#' row count before/after.
turnover_table <- function(bc, md, label) {
  iv <- build_intervals(bc, md)
  if (!nrow(iv)) return(iv)
  iv$interval <- factor(iv$from, levels = INTERVAL_LEVELS)
  n_per_participant <- dplyr::count(iv, participant, name = "n")
  eligible <- n_per_participant$participant[n_per_participant$n >= 2]
  kept <- iv[iv$participant %in% eligible, , drop = FALSE]
  n_dropped_participants <- length(unique(iv$participant)) - length(eligible)
  message(label, ": ", nrow(iv), " intervals from ", length(unique(iv$participant)),
          " participants before the >=2-intervals filter; ", nrow(kept), " intervals from ",
          length(eligible), " participants after (", n_dropped_participants,
          " participant(s) with a single interval excluded)")
  kept
}

#' Identical model to community_turnover_mixed_v1.R's fit_bayesian_turnover():
#' Bray-Curtis distance is bounded [0, 1] regardless of whether the profile it
#' was computed from is relative (SUPER-FOCUS) or raw-concentration
#' (volatiles) -- the vague normal(0.5, 1) prior on each interval mean is
#' equally weakly informative for both.
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

#' Run the turnover test for every type in a single pre-loaded matrix panel
#' (SUPER-FOCUS: one matrix shared across all four types, like MetaCache).
run_panel_matrix <- function(mat, md, types, relative, label) {
  rows <- list()
  bayes_rows <- list()
  for (type in types) {
    pp <- tryCatch(prepare_panel(mat, md, type, relative = relative),
                    error = function(e) { message(label, "/", type, ": ", conditionMessage(e)); NULL })
    if (is.null(pp)) next
    d <- turnover_table(pp$bc, pp$md, paste0(label, "/", type))
    if (!nrow(d)) { message(label, "/", type, ": no eligible intervals, skipping"); next }
    # The model needs every one of the 5 interval categories observed at least
    # once (it estimates a separate mean per interval) -- too few remaining
    # participants/visits after the >=2-intervals filter can leave one
    # interval with zero rows, which brms would otherwise fit silently with
    # one fewer coefficient than expected. Checked explicitly and skipped,
    # rather than caught only as a cryptic "unexpected coefficient names" error.
    if (!setequal(levels(droplevels(d$interval)), INTERVAL_LEVELS)) {
      message(label, "/", type, ": only ", nlevels(droplevels(d$interval)), " of 5 interval",
              " categories observed (", nrow(d), " intervals, ", length(unique(d$participant)),
              " participants) -- too sparse to fit this model, skipping")
      next
    }
    rows[[type]] <- data.frame(panel = label, type = type, d[, c("participant", "from", "to", "interval", "distance")])
    bayes <- fit_bayesian_turnover(d)
    bayes_rows[[type]] <- data.frame(panel = label, type = type, n_intervals = nrow(d),
      n_participants = length(unique(d$participant)), estimate = bayes$estimate,
      ci_low = bayes$ci_low, ci_high = bayes$ci_high, p_decrease = bayes$p_decrease)
  }
  list(intervals = do.call(rbind, rows), tests = do.call(rbind, bayes_rows))
}

run_functional_volatile_turnover <- function() {
  root <- ROOT
  parent <- file.path(root, "output", "04_taxonomic_profiling", "functional_volatile_turnover_bayesian_v1")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out)
  message("Output: ", out)
  set.seed(TURNOVER_SEED)

  md <- load_metadata(root)

  # SUPER-FOCUS: one matrix, all four sample types, relative profiles --
  # exactly parallel to community_turnover_mixed_v1.R's MetaCache handling.
  superfocus <- load_superfocus(root, md)
  sf_result <- run_panel_matrix(superfocus$mat, md, names(TYPE_LABELS), relative = TRUE, label = "SUPER_FOCUS")

  # Volatiles: no grain data exists, and ML/WL each come from a separate
  # workbook loaded per type -- not a single shared matrix like SUPER-FOCUS.
  vol_rows <- list(); vol_tests <- list()
  for (type in c("ML", "WL")) {
    vol <- tryCatch(load_volatiles(root, md, type),
                      error = function(e) { message("volatiles/", type, ": ", conditionMessage(e)); NULL })
    if (is.null(vol)) next
    res <- run_panel_matrix(vol$mat, md, type, relative = FALSE, label = "volatiles")
    vol_rows[[type]] <- res$intervals
    vol_tests[[type]] <- res$tests
  }
  vol_result <- list(intervals = do.call(rbind, vol_rows), tests = do.call(rbind, vol_tests))

  intervals <- rbind(sf_result$intervals, vol_result$intervals)
  write.csv(intervals, file.path(out, "intervals.csv"), row.names = FALSE)

  tests <- rbind(sf_result$tests, vol_result$tests)
  write.csv(tests, file.path(out, "bayesian_turnover_tests.csv"), row.names = FALSE)

  writeLines(c(
    "Tests whether within-participant SUPER-FOCUS (functional) and volatile",
    "(metabolomic) turnover differs by consecutive-visit interval, the same",
    "question and model community_turnover_mixed_v1.R answers for MetaCache",
    "taxonomic profiles.",
    "SUPER-FOCUS: all four sample types (MG/WG/ML/WL), relative profiles.",
    "Volatiles: ML/WL only (no grain volatile data exists), raw concentrations.",
    "Model: brms::brm(distance ~ 0 + interval + (1 | participant)), interval a",
    "5-level factor (T1->T2 .. T5->T6), vague prior normal(0.5, 1) on each",
    "interval mean, 4 chains, 4000 iter.",
    "Reported quantity: posterior distribution of (mean of the four later",
    "intervals) - (first interval) -- posterior mean, 95% credible interval,",
    "and P(decrease), the posterior probability turnover was lower after the",
    "first interval than during it.",
    "Participants need >=2 intervals within a panel/type to be included.",
    paste("Run directory:", out)), file.path(out, "METHODS.txt"))
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

  print(tests, row.names = FALSE)
  invisible(list(output = out, intervals = intervals, tests = tests))
}

if (isTRUE(getOption("functional_volatile_turnover_bayesian_v1.run", TRUE))) run_functional_volatile_turnover()
