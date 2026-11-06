# ==============================================================================
# 04_community_stability_v3.R -- MetaCache, SUPER-FOCUS and volatile metabolomics
# Temporal stability is within-participant change, NOT between-sample convergence.
# Primary endpoint: participant mean Bray-Curtis over T1->T2 and T2->T3 (early)
# versus T3->T4, T4->T5 and T5->T6 (late). All intervals span four weeks.
# T0 and gaps are excluded. Participants need >=1 interval in each window;
# complete five-interval trajectories provide a sensitivity analysis.
# Paired Wilcoxon tests operate on ONE early/late pair per participant, two-sided.
# The stage cadence is an ANALYST ASSUMPTION, not a measured variable: metadata
# records Stage only, so T1..T6 are taken to be weeks 1,5,9,13,17,21 and every
# interval spans four assumed weeks. Interval tables label those columns
# assumed_week_*; the stage labels remain the only measured time variable.
# BH correction is across panels with a computable p-value, separately for
# primary/sensitivity; the family rule is fixed a priori and its realised size
# is reported in METHODS.txt.
# CIs bootstrap participants, preserving their early/late pairing. Independent
# participants are assumed; school/household clustering is not modelled.
# Missing visits may be informative. No equivalence margin is specified, so a
# non-significant test is NOT evidence of stability or a plateau.
# Within-stage between-participant convergence is exported descriptively only.
# Distances come from vegan::vegdist and CIs from boot::boot.ci (the project
# image provides both); readxl is needed only for the volatile workbooks.
# Volatile panels use raw concentrations while MetaCache/SUPER-FOCUS use row
# proportions, so their distance values are not directly comparable.
# Source with options(community_stability_v3.run = FALSE) for function tests.
# ==============================================================================

# Honors the documented contract (REPRODUCIBILITY.md: scripts are standalone and
# run from the repository root). KEFIR4ALL_ROOT overrides; here::here(), the
# script location and the working directory are tried. The data file itself is
# the test, so a wrong root fails loudly instead of silently reading nothing.
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
STAGES <- paste0("T", 1:6)
# Assumed stage cadence, not a measured variable (see header).
WEEKS_ASSUMED <- c(1, 5, 9, 13, 17, 21)
BOOT_R <- 5000L    # bootstrap replicates for the participant-mean CI
MIN_CI_N <- 5L     # below this, a bootstrap interval is not interpretable
SEED <- 20260916L
TYPE_LABELS <- c(MG = "Milk kefir (grain)", WG = "Water kefir (grain)",
                 ML = "Milk kefir (liquid)", WL = "Water kefir (liquid)")
normalize_id <- function(x) {
  x <- basename(gsub("\\\\", "/", trimws(as.character(x))))
  x <- sub("_host_removed_R[0-9]\\.fastq\\.gz$", "", x)
  gsub("-", "_", x, fixed = TRUE)
}
read_table <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
require_columns <- function(x, cols) {
  if (!all(cols %in% names(x))) stop("Missing columns: ", paste(setdiff(cols, names(x)), collapse = ", "))
}
load_metadata <- function(root) {
  md <- read_table(file.path(root, "data", "kefir4all_sample_metadata_v2.csv"))
  require_columns(md, c("merge_column", "Sample", "kefir type", "Stage"))
  md$key <- normalize_id(md$merge_column)
  md$Stage <- toupper(trimws(md$Stage))
  md$type <- trimws(md[["kefir type"]])
  md$participant <- toupper(trimws(as.character(md$Sample)))
  numeric_id <- grepl("^(ID)?0*[0-9]+$", md$participant)
  md$participant[numeric_id] <- paste0("ID", as.integer(sub("^ID", "", md$participant[numeric_id])))
  # Baselines and unrelated cohorts cannot contribute to temporal inference.
  # Filter them before join validation (one T0 medium control has conflicting labels).
  md <- md[!is.na(md$type) & md$type %in% names(TYPE_LABELS) &
             !is.na(md$Stage) & md$Stage %in% STAGES &
             !is.na(md$participant) & nzchar(md$participant), , drop = FALSE]
  signature <- paste(md$participant, md$type, md$Stage, sep = "|")
  conflicts <- vapply(split(signature, md$key), function(z) length(unique(z)) > 1, logical(1))
  if (any(conflicts)) stop("Conflicting metadata for ", sum(conflicts), " normalized sample keys")
  md <- md[!duplicated(md$key), , drop = FALSE]
  md[!is.na(md$type) & md$type %in% names(TYPE_LABELS) &
       !is.na(md$Stage) & md$Stage %in% c("T0", STAGES) &
       !is.na(md$participant) & nzchar(md$participant), , drop = FALSE]
}

# Strict numeric conversion: workbook blanks are non-detections (zero), but
# unexpected text, negative values and infinities are errors, never clipped.
numeric_matrix <- function(x, blanks_zero = FALSE) {
  raw <- as.matrix(x)
  m <- matrix(suppressWarnings(as.numeric(raw)), nrow(raw), ncol(raw), dimnames = dimnames(raw))
  blank <- is.na(raw) | trimws(raw) == ""
  if (any(is.na(m) & !blank)) stop("Non-numeric abundance values")
  if (blanks_zero) m[blank] <- 0
  if (any(!is.finite(m)) || any(m < 0)) stop("Missing, negative or infinite abundance values")
  m
}
load_taxonomy <- function(root, md) {
  x <- read_table(file.path(root, "data", "04_metacache_total_full_profile.csv"))
  require_columns(x, "sample_id")
  keys <- normalize_id(x$sample_id)
  keep <- keys %in% md$key
  if (anyDuplicated(keys[keep])) stop("Duplicate normalized MetaCache sample IDs")
  m <- numeric_matrix(x[keep, setdiff(names(x), "sample_id"), drop = FALSE], blanks_zero = TRUE)
  rownames(m) <- keys[keep]
  list(mat = m, notes = "")
}
load_superfocus <- function(root, md) {
  folder <- Sys.getenv("KEFIR4ALL_FUNC_DIR", file.path(root, "data", "private", "functional"))
  paths <- c(file.path(folder, "Temp_4Liam", "05_short_read_functional_profiling", "05_superfocus", "superfocus_p3.csv"),
             file.path(folder, "superfocus_p3.csv"),
             file.path(root, "data", "functional", "superfocus_p3.csv"),
             file.path(root, "data", "metabolomics", "superfocus_p3.csv"))
  path <- paths[file.exists(paths)][1]
  if (is.na(path)) stop("SUPER-FOCUS input missing under ", folder)
  x <- read_table(path)
  require_columns(x, c("sample_id", "Subsystem Level 3", "reads_%"))
  x$key <- normalize_id(x$sample_id)
  x <- x[x$key %in% md$key, , drop = FALSE]
  if (!nrow(x)) stop("No SUPER-FOCUS sample IDs match metadata")
  x$value <- as.vector(numeric_matrix(x["reads_%"]))
  x$feature <- x[["Subsystem Level 3"]]
  if (anyNA(x$feature) || any(!nzchar(trimws(x$feature)))) stop("Missing SUPER-FOCUS feature labels")
  # Duplicate sample/pathway entries are summed; absent pathways are zero.
  m <- as.matrix(xtabs(value ~ key + feature, data = x))
  list(mat = m, notes = "")
}

load_volatiles <- function(root, md, type) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read volatile workbooks")
  file <- if (type == "ML") "Final results table all milk kefirs 0913-001.xlsx" else
    "Final results table all water kefirs 0913-002.xlsx"
  x <- as.data.frame(readxl::read_excel(file.path(root, "data", file), .name_repair = "minimal"))
  if ("Name" %in% names(x)) names(x)[names(x) == "Name"] <- "Compound"
  require_columns(x, "Compound")
  cn <- toupper(gsub("__", "_", names(x), fixed = TRUE))
  keep <- grepl("^0*[0-9]{1,3}_T[1-6]$", cn)
  if (!any(keep)) stop("No participant/stage columns in ", file)
  num <- as.integer(sub("_T[1-6]$", "", cn[keep]))
  stage <- sub("^.*_", "", cn[keep])
  keys <- paste(type, paste0("ID", num), stage, sep = "|")
  if (anyDuplicated(keys)) stop("Duplicate volatile participant/stage columns")
  mdx <- md[md$type == type, , drop = FALSE]
  mk <- paste(mdx$type, mdx$participant, mdx$Stage, sep = "|")
  # Workbook identifies visits, not sequencing libraries. Equivalent metadata
  # records for a visit are collapsed; the selected key is only a join label.
  mdx <- mdx[!duplicated(mk), , drop = FALSE]
  mk <- mk[!duplicated(mk)]
  pos <- match(keys, mk)
  matched <- !is.na(pos)
  if (!any(matched)) stop("No volatile columns match metadata for ", type)
  # pos and keys are aligned to the kept columns only; keep is aligned to all
  # columns, so map back through which(keep) rather than mixing the two indices.
  sel <- which(keep)[matched]
  dropped <- keys[!matched]
  notes <- if (length(dropped))
    paste0(length(dropped), " volatile column(s) without metadata: ", paste(dropped, collapse = ", ")) else ""
  message(type, " volatile columns: ", sum(keep), "; unmatched: ", length(dropped),
          "; non-participant/baseline/annotation columns excluded: ", sum(!keep))
  if (length(dropped))
    warning("Unmatched volatile columns excluded: ", paste(dropped, collapse = ", "))
  # Original workbook is features x samples. Select COLUMNS, then transpose.
  m <- t(numeric_matrix(x[, sel, drop = FALSE], blanks_zero = TRUE))
  compound <- as.character(x$Compound)
  if (anyNA(compound) || any(!nzchar(trimws(compound)))) stop("Missing compound labels")
  colnames(m) <- make.unique(compound)
  rownames(m) <- mdx$key[pos[matched]]
  list(mat = m, notes = notes)
}

# Bray-Curtis via vegan::vegdist, keeping the input checks this pipeline relies
# on (finite, nonnegative, nonempty profiles with positive row totals).
bray_matrix <- function(m) {
  if (!nrow(m) || any(!is.finite(m)) || any(m < 0) || any(rowSums(m) <= 0))
    stop("Bray-Curtis requires finite nonnegative, nonempty sample profiles")
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Install vegan for Bray-Curtis distances")
  bc <- as.matrix(vegan::vegdist(m, method = "bray"))
  diag(bc) <- 0
  bc
}

prepare_panel <- function(m, md, type, relative) {
  pos <- match(rownames(m), md$key)
  if (anyNA(pos)) stop("Matrix contains unmatched sample IDs")
  meta <- md[pos, , drop = FALSE]
  keep <- meta$type == type & meta$Stage %in% STAGES
  m <- m[keep, , drop = FALSE]; meta <- meta[keep, , drop = FALSE]
  if (!nrow(m)) stop("No T1-T6 samples for ", type)
  if (anyDuplicated(rownames(m))) stop("Duplicate sample rows")
  if (any(!is.finite(m)) || any(m < 0)) stop("Invalid abundance matrix")
  zero <- rowSums(m) == 0
  exclusions <- meta[zero, c("key", "participant", "type", "Stage"), drop = FALSE]
  m <- m[!zero, , drop = FALSE]; meta <- meta[!zero, , drop = FALSE]
  if (!nrow(m)) stop("All sample profiles have zero total")
  visit <- paste(meta$participant, meta$Stage)
  dup <- duplicated(visit) | duplicated(visit, fromLast = TRUE)
  exclusions$reason <- rep("zero_profile", nrow(exclusions))
  if (any(dup)) {
    extra <- meta[dup, c("key", "participant", "type", "Stage"), drop = FALSE]
    extra$reason <- "ambiguous_duplicate_visit"
    exclusions <- rbind(exclusions, extra)
    m <- m[!dup, , drop = FALSE]; meta <- meta[!dup, , drop = FALSE]
  }
  if (!nrow(m)) stop("No unambiguous positive profiles")
  m <- m[, colSums(m) > 0, drop = FALSE]
  # MetaCache and SUPER-FOCUS: relative profiles, controlling total abundance.
  # Volatiles: unscaled measured abundances, retaining concentration differences.
  if (relative) m <- m / rowSums(m)
  stopifnot(identical(rownames(m), meta$key))
  list(mat = m, md = meta, excluded = exclusions, bc = bray_matrix(m))
}

build_intervals <- function(bc, md) {
  rows <- list()
  for (id in unique(md$participant)) {
    idx <- which(md$participant == id)
    at <- idx[match(STAGES, md$Stage[idx])]
    for (s in 1:5) {
      if (anyNA(at[c(s, s + 1L)])) next  # never bridge missing visits
      i <- at[s]; j <- at[s + 1L]
      rows[[length(rows) + 1L]] <- data.frame(participant = id,
        from = STAGES[s], to = STAGES[s + 1L], assumed_week_from = WEEKS_ASSUMED[s],
        assumed_week_to = WEEKS_ASSUMED[s + 1L], window = if (s <= 2) "early" else "late",
        distance = bc[cbind(i, j)], stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(participant = character(), from = character(),
    to = character(), assumed_week_from = numeric(), assumed_week_to = numeric(),
    window = character(), distance = numeric()))
  do.call(rbind, rows)
}
participant_summary <- function(intervals) {
  rows <- lapply(unique(intervals$participant), function(id) {
    d <- intervals[intervals$participant == id, ]
    early <- d$distance[d$window == "early"]; late <- d$distance[d$window == "late"]
    data.frame(participant = id, n_early = length(early), n_late = length(late),
      early = if (length(early)) mean(early) else NA_real_,
      late = if (length(late)) mean(late) else NA_real_)
  })
  if (!length(rows)) return(data.frame(participant = character(), n_early = integer(),
    n_late = integer(), early = numeric(), late = numeric(), delta = numeric()))
  out <- do.call(rbind, rows)
  out$delta <- out$late - out$early
  out
}

#' Descriptive between-participant convergence: within-stage pairwise
#' Bray-Curtis among different participants sharing a stage. Pairs are NOT
#' independent and these summaries are exported without any significance test.
within_stage_convergence <- function(bc, md) {
  empty <- data.frame(Stage = character(), n_samples = integer(),
    n_pairs = integer(), mean = numeric(), median = numeric(),
    q25 = numeric(), q75 = numeric())
  if (nrow(md) < 2) return(empty)
  cp <- utils::combn(seq_len(nrow(md)), 2)
  sel <- md$participant[cp[1, ]] != md$participant[cp[2, ]] &
         md$Stage[cp[1, ]] == md$Stage[cp[2, ]]
  if (!any(sel)) return(empty)
  i <- cp[1, sel]; j <- cp[2, sel]
  d <- data.frame(Stage = md$Stage[i], distance = bc[cbind(i, j)], stringsAsFactors = FALSE)
  do.call(rbind, lapply(sort(unique(d$Stage)), function(s) {
    v <- d$distance[d$Stage == s]
    data.frame(Stage = s, n_samples = sum(md$Stage == s), n_pairs = length(v),
      mean = mean(v), median = median(v), q25 = unname(stats::quantile(v, 0.25)),
      q75 = unname(stats::quantile(v, 0.75)))
  }))
}

#' Two-sided paired Wilcoxon signed-rank on participant deltas. The exact null
#' distribution needs untied, non-zero differences; otherwise the normal
#' approximation is used, and which one was applied is reported (p_exact) rather
#' than left implicit. All-zero deltas leave the statistic undefined, so p is NA
#' -- never 1, which would wrongly take up weight in the BH family.
wilcox_test <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) return(list(p = NA_real_, exact = NA))
  if (!any(x != 0)) {
    warning("Signed-rank test undefined: all ", n, " participant deltas are zero; p reported as NA")
    return(list(p = NA_real_, exact = FALSE))
  }
  exact <- !anyDuplicated(x) && n < 50L
  list(p = stats::wilcox.test(x, alternative = "two.sided", exact = exact)$p.value,
       exact = exact)
}

#' Participant bootstrap CI for the mean delta. BCa is preferred because it
#' corrects for bias and skew at the small n of these panels; percentile is the
#' fallback when BCa is unavailable. Below MIN_CI_N, or when the resampled mean
#' is degenerate, the interval is NA rather than a misleading number.
boot_ci <- function(x, R = BOOT_R, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < MIN_CI_N) return(c(lo = NA_real_, hi = NA_real_))
  if (!requireNamespace("boot", quietly = TRUE)) stop("Install boot for bootstrap CIs")
  b <- boot::boot(x, function(d, i) mean(d[i]), R = R)
  ci <- tryCatch(boot::boot.ci(b, conf = conf, type = "bca")$bca, error = function(e) NULL)
  if (is.null(ci)) ci <- tryCatch(boot::boot.ci(b, conf = conf, type = "perc")$perc, error = function(e) NULL)
  ends <- if (is.null(ci)) c(NA_real_, NA_real_) else tail(as.numeric(ci), 2)
  if (any(!is.finite(ends))) return(c(lo = NA_real_, hi = NA_real_))
  c(lo = ends[1], hi = ends[2])
}

#' One paired early/late observation per participant: the mean over available
#' intervals in each window (>= 1 required per window). Primary test: two-sided
#' Wilcoxon signed-rank on participant deltas. Sensitivity: complete
#' five-interval trajectories only. BH correction across panels is applied
#' separately for primary and sensitivity families.
analyse_panel <- function(bc, md) {
  iv <- build_intervals(bc, md)
  ps <- participant_summary(iv)
  elig <- ps[ps$n_early >= 1 & ps$n_late >= 1, , drop = FALSE]
  # Empty paired cohorts still produce descriptive summaries and an NA test.
  primary_d <- elig$late - elig$early
  full <- ps[ps$n_early == 2 & ps$n_late == 3, , drop = FALSE]
  full_d <- full$delta
  ci <- boot_ci(primary_d)
  wt <- wilcox_test(primary_d)
  ws <- wilcox_test(full_d)
  list(participants = elig, intervals = iv, complete = full,
       primary = list(mean_early = mean(elig$early), mean_late = mean(elig$late),
                      mean_delta = mean(primary_d), p = wt$p, exact = wt$exact,
                      ci_lo = ci[["lo"]], ci_hi = ci[["hi"]], n = length(primary_d)),
       sensitivity = list(mean_delta = if (length(full_d)) mean(full_d) else NA_real_,
                          p = ws$p, exact = ws$exact,
                          n = length(full_d)),
       convergence = within_stage_convergence(bc, md))
}

write_outputs <- function(res, dir) {
  named <- Filter(function(z) !is.null(z), res)
  part <- do.call(rbind, lapply(names(named), function(k) {
    z <- named[[k]]
    data.frame(layer = rep(k, nrow(z$participants)), z$participants)
  }))
  conv <- do.call(rbind, lapply(names(named), function(k)
    data.frame(layer = rep(k, nrow(named[[k]]$convergence)), named[[k]]$convergence)))
  rows <- lapply(names(named), function(k) {
    z <- named[[k]]; p1 <- z$primary; p2 <- z$sensitivity
    data.frame(layer = k, n_participants = p1$n, mean_early = p1$mean_early,
      mean_late = p1$mean_late, mean_delta = p1$mean_delta,
      ci_low = p1$ci_lo, ci_high = p1$ci_hi, p_value = p1$p, p_exact = p1$exact,
      sensitivity_n = p2$n, sensitivity_mean_delta = p2$mean_delta,
      sensitivity_p_value = p2$p, sensitivity_p_exact = p2$exact)
  })
  write.csv(part, file.path(dir, "participant_deltas.csv"), row.names = FALSE, na = "")
  write.csv(conv, file.path(dir, "between_participant_convergence.csv"), row.names = FALSE, na = "")
  tests <- do.call(rbind, rows)
  tests$p_BH <- p.adjust(tests$p_value, "BH")
  tests$sensitivity_p_BH <- p.adjust(tests$sensitivity_p_value, "BH")
  write.csv(tests, file.path(dir, "primary_and_sensitivity_tests.csv"), row.names = FALSE, na = "")
  invisible(tests)
}

plot_panel <- function(z, title, path) {
  draw <- function() {
    par(mfrow = c(1, 2), mar = c(5, 4, 4, 1))
    d <- z$participants
    plot(c(1, 2), c(0, 1), type = "n", xaxt = "n", xlab = "Window",
         ylab = "Within-participant Bray-Curtis", main = title)
    axis(1, at = c(1, 2), labels = c("Early (weeks 1-9)", "Late (weeks 9-21)"))
    for (i in seq_len(nrow(d))) lines(c(1, 2), c(d$early[i], d$late[i]), col = "grey75", type = "b", pch = 16, cex = 0.5)
    lines(c(1, 2), c(mean(d$early), mean(d$late)), col = "#1F78B4", lwd = 3, type = "b", pch = 16)
    mtext(paste("Paired participants:", nrow(d)), side = 3, cex = 0.8)
    cv <- z$convergence
    plot(c(1, 21), c(0, 1), type = "n", xlab = "Week (assumed cadence)",
         ylab = "Between-participant Bray-Curtis",
         main = "Convergence (descriptive only)")
    if (nrow(cv)) {
      week <- WEEKS_ASSUMED[match(cv$Stage, STAGES)]
      segments(week, cv$q25, week, cv$q75, col = "grey50")
      lines(week, cv$median, type = "b", pch = 16, col = "#33A02C")
    }
    mtext("Median and IQR; pairs are dependent", side = 3, cex = 0.8)
  }
  pdf(paste0(path, ".pdf"), width = 6, height = 4)
  tryCatch(draw(), finally = dev.off())
  tiff(paste0(path, ".tiff"), width = 6, height = 4, units = "in", res = 300, compression = "lzw")
  tryCatch(draw(), finally = dev.off())
}

run_stability <- function() {
  root <- find_root()
  # Fresh directory prevents historical results being mistaken for this run.
  parent <- file.path(root, "output", "04_taxonomic_profiling", "community_stability_v3")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  dir.create(out); dir.create(file.path(out, "figures"))
  message("Output: ", out)
  set.seed(SEED)
  md <- load_metadata(root)
  # Single source of truth for each panel: which kefir types it contributes and
  # whether profiles become proportions before distances are computed are
  # declared together here, rather than each being inferred separately from
  # the panel name elsewhere (which could silently drift out of sync).
  panels <- list(
    list(key = "MetaCache", relative = TRUE, types = names(TYPE_LABELS),
         load = function() load_taxonomy(root, md)),
    list(key = "SUPER_FOCUS", relative = TRUE, types = names(TYPE_LABELS),
         load = function() load_superfocus(root, md)),
    list(key = "volatiles_ML", relative = FALSE, types = "ML",
         load = function() load_volatiles(root, md, "ML")),
    list(key = "volatiles_WL", relative = FALSE, types = "WL",
         load = function() load_volatiles(root, md, "WL")))
  res <- list(); manifest <- list()
  record <- function(panel, status, reason, n_samples = 0, n_participants = 0) {
    manifest[[length(manifest) + 1L]] <<- data.frame(panel = panel, status = status,
      reason = reason, n_samples = n_samples, n_participants = n_participants)
  }
  for (cfg in panels) {
    loaded <- tryCatch(cfg$load(), error = function(e) e)
    if (inherits(loaded, "error")) {
      record(cfg$key, "SKIPPED", conditionMessage(loaded)); next
    }
    for (type in cfg$types) {
      # Single-type panels keep their own name, so volatile keys read
      # "volatiles_ML" rather than the duplicated "volatiles_ML_ML".
      key <- if (length(cfg$types) == 1L) cfg$key else paste(cfg$key, type, sep = "_")
      result <- tryCatch({
        pp <- prepare_panel(loaded$mat, md, type, relative = cfg$relative)
        write.csv(pp$md[, c("key", "participant", "type", "Stage")], file.path(out, paste0(key, "_sample_alignment.csv")), row.names = FALSE)
        write.csv(pp$excluded, file.path(out, paste0(key, "_zero_profile_exclusions.csv")), row.names = FALSE)
        counts <- data.frame(Stage = STAGES, n_samples = as.integer(table(factor(pp$md$Stage, levels = STAGES))))
        write.csv(counts, file.path(out, paste0(key, "_sample_counts.csv")), row.names = FALSE)
        z <- analyse_panel(pp$bc, pp$md)
        write.csv(z$intervals, file.path(out, paste0(key, "_intervals.csv")), row.names = FALSE)
        write.csv(participant_summary(z$intervals), file.path(out, paste0(key, "_all_participants.csv")), row.names = FALSE)
        # This per-panel plot is a descriptive diagnostic, not a numbered
        # manuscript figure, so it carries no figure number in its filename.
        # Stem kept short on purpose: output lives under a deep OneDrive path and
        # "community_stability_" pushed the SUPER_FOCUS TIFFs past Windows' 260
        # character MAX_PATH, silently losing them.
        plot_panel(z, paste(cfg$key, TYPE_LABELS[type]), file.path(out, "figures", paste0("stability_", key)))
        list(z = z, n_samples = nrow(pp$md))
      }, error = function(e) e)
      if (inherits(result, "error")) {
        record(key, "ERROR", conditionMessage(result))
      } else {
        res[[key]] <- result$z
        status <- if (result$z$primary$n > 0) "OK" else "DESCRIPTIVE_ONLY"
        # Loader notes (e.g. volatile columns with no metadata) reach the manifest,
        # so anything warned about on the console is also visible in the outputs.
        reason <- paste(c(if (status == "OK") NULL else
                            "No participant has adjacent intervals in both windows",
                          if (nzchar(loaded$notes)) loaded$notes), collapse = "; ")
        record(key, status, reason, result$n_samples, result$z$primary$n)
      }
    }
  }
  man <- do.call(rbind, manifest)
  write.csv(man, file.path(out, "community_stability_v3_manifest.csv"), row.names = FALSE)
  tests <- if (length(res)) write_outputs(res, out) else NULL
  bh_n <- function(v) if (is.null(v)) 0L else sum(!is.na(v))
  writeLines(c("MetaCache / SUPER-FOCUS / volatile metabolomics only.",
    "MetaCache and SUPER-FOCUS: row proportions. Volatiles: raw measured abundances; blanks = zero.",
    "Volatile and taxonomic distances are therefore not directly comparable.",
    "ASSUMED time axis: metadata records Stage only; T1-T6 are taken as weeks 1,5,9,13,17,21.",
    "Early: assumed weeks 1->5 and 5->9. Late: 9->13, 13->17, 17->21. No gap bridging; T0 excluded.",
    "Participant mean over available four-week intervals; >=1 per window. Sensitivity: all five intervals.",
    "MetaCache missing feature cells are zero (legacy convention). Negative/infinite values rejected.",
    "Sequencing visits with multiple libraries are excluded entirely; reasons exported in exclusions tables.",
    "Volatile visits map to distinct participant/type/stage metadata; library keys are join labels only.",
    "Independent participants assumed; school/household clustering not modelled.",
    "Two-sided paired Wilcoxon signed-rank (stats::wilcox.test); exact when differences are untied and non-zero.",
    "Mean delta = late minus early; negative indicates reduced turnover, not proof of stability.",
    paste0("95% bootstrap CI on the participant mean delta: ", BOOT_R, " replicates, BCa (percentile fallback), seed ",
           SEED, "; computed only when n >= ", MIN_CI_N, "."),
    "Signed-rank test concerns paired-difference distribution (symmetry assumed), not specifically its mean.",
    paste0("BH family rule fixed a priori as every panel with a computable p-value: primary N = ",
           bh_n(tests$p_value), ", sensitivity N = ", bh_n(tests$sensitivity_p_value),
           ". Undefined (all-zero) p-values are NA and excluded from the family."),
    "p_exact records whether the exact (TRUE) or normal-approximation (FALSE) null was used.",
    "No equivalence test or plateau claim. Between-participant convergence is descriptive only.",
    "Missingness may bias available-interval comparisons; review complete-trajectory sensitivity.",
    "Participant IDs in outputs: treat alignment/interval/participant tables as potentially private.",
    paste("Run directory:", out)), file.path(out, "METHODS.txt"))
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  print(man, row.names = FALSE)
  if (any(man$status %in% c("ERROR", "SKIPPED")))
    warning("Incomplete run: consult manifest for ERROR/SKIPPED panels")
  invisible(list(output = out, manifest = man, results = res, tests = tests))
}
if (isTRUE(getOption("community_stability_v3.run", TRUE))) run_stability()

