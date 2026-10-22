# ===========================================================================
# community_changes_v2.R
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Kefir4All is a citizen-science trial: participants were given a shared milk
# or water kefir starter culture and fermented it themselves for 21 weeks,
# sending in a sample at each of six visits (T1-T6, assumed to fall at weeks
# 1, 5, 9, 13, 17, 21 - the metadata records only the visit label, never a
# measured date). This script describes how much each participant's kefir
# community drifted from the starting point, and how much it moved between
# consecutive visits, across the trial.
#
# It analyses ten "panels" - one per combination of data source and kefir
# type:
#   data sources: MetaCache (rank-qualified taxonomic features), SUPER-FOCUS (functional
#                 subsystems, level 3), volatiles (measured VOC concentrations)
#   kefir types:  MG/WG (milk/water kefir GRAIN), ML/WL (milk/water kefir
#                 LIQUID) - volatiles were only measured for the two liquid
#                 types, so there are 4 + 4 + 2 = 10 panels in total.
# Everything - data loading, descriptive summaries, and the figures - lives
# in this one file; there is no separate script to run afterwards.
#
# WHY (ALMOST) NO SIGNIFICANCE TESTS
# -----------------------------------
# This is a deliberately descriptive counterpart to community_changes.R.
# It reports distances and their patterns without p-values or plateau claims,
# with one specific exception: a single, pre-planned pooled taxonomic test of
# whether the TREND in consecutive-visit distance changes after week 9 (see
# "THE ONE STATISTICAL TEST" below). Everything else - the stage-distance
# summary, the peak-stage tally, the figures - stays purely descriptive. The
# tested companion (community_changes.R) runs the broader battery of
# omnibus/pairwise/PERMANOVA tests; both scripts use the same full-rank
# MetaCache input.
#
# THE ONE STATISTICAL TEST
# -------------------------
# S_5_segmented() pools the four taxonomic panels (MetaCache MG/WG/ML/WL) and
# fits consecutive-visit Bray-Curtis distance against a segmented (piecewise)
# function of the interval's midpoint week, with a knot at week 9 - the
# project's own pre-existing hypothesised water-kefir peak week (not chosen
# by inspecting this dataset's medians). Every participant x type trajectory
# gets its own fixed-effect intercept; uncertainty is clustered by
# PARTICIPANT (not participant x type), since one participant can contribute
# more than one kefir type and those are not independent. Trajectories with
# fewer than two usable intervals are dropped (a one-point trajectory has no
# within-trajectory information once it has its own saturated intercept).
# Inference is a null-imposed, participant-clustered wild-cluster bootstrap
# (Cameron, Gelbach & Miller 2008) with CR1 standard errors, because ordinary
# lm() p-values assume independent residuals and CR0 does not correct for a
# small number of clusters. The test asks whether the LATER slope differs
# from the EARLY slope - it does NOT test whether the later slope is zero,
# so it cannot by itself demonstrate a plateau; a shallower but still
# nonzero later trend is not ruled out. Result:
# community_changes_v2_segmented_slope_test.csv. This is one single planned
# contrast, not the broader battery community_changes.R runs, and multiple-
# testing correction across any wider family of tests is not attempted here.
#
# THE TWO DISTANCE MEASURES
# -------------------------
# For every panel, two things are measured (Bray-Curtis distance for the
# compositional data sources; Euclidean distance for volatiles, since row
# proportions would discard the concentration scale there):
#   1. reference - how far has this visit's community drifted from the
#      starting point? The distance from the visit's sample back to a T0
#      ("time zero") starter sample of the same kefir type. See "THE T0
#      REFERENCE" below for how that starter sample is chosen.
#   2. interval  - how much did the community move between two consecutive
#      visits (T1-T2, T2-T3, ..., T5-T6)? The same distance metric, computed
#      between each participant's own pair of consecutive samples.
#
# THE T0 REFERENCE
# -----------------
# T0 is not a per-participant baseline - every participant started from the
# same shared, laboratory-prepared starter culture, so there is no "each
# participant's own T0" to compare against. That starter culture was itself
# prepared in one of several ways, though, encoded directly in the lab
# sample's ID (e.g. "TG_FU_T3_MG_..."). Verified against the actual kefir-type
# field in the metadata (not just guessed from the letters):
#   milk (MG/ML):  FU = whole/full-fat, GO = goat, SK = skim, LO = low-fat
#   water (WG/WL): AP = white sugar + dried apricot, BR = brown sugar,
#                  FG = white sugar + dried fig, FL = white sugar + dried fig
#                  and fresh lemon
#
# For LIQUID panels (ML/WL) this matters biologically: the liquid is the
# product of fermenting the grain in one of these specific media, so a
# participant who fermented skim milk should be compared to the SK starter,
# not an arbitrary one. Every liquid visit is therefore matched to the T0 of
# its own participant's fermentation substrate, looked up from the private
# citizen-scientist metadata (see "DATA THIS SCRIPT NEEDS" below).
#
# For GRAIN panels (MG/WG), medium doesn't apply the same way - participants
# share one grain culture regardless of what they later ferment it in - so
# every grain visit uses one single, consistent T0 reference (the same
# physical starter sample across MetaCache/SUPER-FOCUS/volatiles).
#
# Either way, if a participant's substrate is unknown or their matched
# starter has no profile in a given data source, the visit falls back to one
# default T0 per type (the first T0 row for that type in metadata order).
# Which reference was actually used, for every visit, is recorded in
# community_changes_v2_t0_reference.csv (the per-panel default) and
# <panel>_reference_distances.csv (the actual reference per visit, and
# whether it was substrate-matched).
#
# DESCRIPTIVE SUMMARIES (no statistical tests)
# ----------------------------------------------
# Per panel: a mean/median distance-to-starter per stage (every available
# profile, no complete-trajectory filter), and a peak-stage tally (which
# visit holds each participant's maximum distance to the starter, among
# participants with at least four visits). Both are plain available-case
# reads of the data, not tests - a difference between stages or participants
# can partly reflect attrition (later visits have fewer participants) rather
# than within-participant change, and that is not something a descriptive
# summary alone can rule out.
#
# DATA THIS SCRIPT NEEDS
# -----------------------
# Public (distributed with the repository, under data/):
#   kefir4all_sample_metadata_v2.csv                  - sample/participant/stage map
#   04_metacache_total_full_profile.csv               - classified MetaCache rank::name read counts
#   functional/superfocus_p3.csv                       - SUPER-FOCUS level-3 profiles
#   "Final results table all milk/water kefirs *.xlsx" - volatile workbooks
# Private (NOT distributed with the public repository - see
# data/private/README.md for how to obtain it):
#   data/private/Citizen Scientist metadata_v8.csv - each participant's own
#     fermentation substrate (category_confirmed), needed for the liquid-panel
#     T0 substrate matching described above. If this file is missing, the
#     script stops immediately with a clear error rather than silently
#     falling back to some public-data-only mode.
#
# OUTPUTS
# -------
# Tables land under output/04_taxonomic_profiling/community_changes_v2/:
#   <panel>_reference_distances.csv, <panel>_interval_distances.csv - raw distances
#   community_changes_v2_t0_reference.csv     - which T0 sample was used, per panel
#   community_changes_v2_peak_stage.csv       - per-participant peak-visit tally
#   community_changes_v2_stage_distance_summary.csv - descriptive stage-distance summary
#   community_changes_v2_manifest.csv         - which panels ran OK / were skipped, and why
#   community_changes_v2_segmented_slope_test.csv - the one statistical test (see header)
#   METHODS.txt                                - full prose methods for every number above
# Figures land under figures/:
#   community_changes_v2_figure.{tiff,pdf}        - reference vs interval distance, all panels
#   community_changes_v2_peak_stage.{tiff,pdf}    - per-panel peak-visit bar charts
#
# USAGE
# -----
# Run from the kefir4all_review_history root (S_0_find_root() locates it by
# searching for data/kefir4all_sample_metadata_v2.csv):
#   Rscript scripts/r_scripts/04_taxonomic_profiling/community_changes_v2.R
# Or via the pinned, reproducible Docker image (see docker/README.md) - this
# is the environment that should be treated as authoritative, since it is the
# one guaranteed to have every package this script actually needs:
#   docker run --rm -v "$(pwd):/home/rstudio/kefir4all" \
#     -w /home/rstudio/kefir4all ghcr.io/liamhwalsh/kefir4all:4.4.2-slim \
#     /opt/R/4.4.2/bin/Rscript scripts/r_scripts/04_taxonomic_profiling/community_changes_v2.R
#
# LAYOUT OF THIS FILE
# --------------------
# Function definitions and the script that applies them are kept in two
# separate, non-interleaved blocks below: every function is defined once, in
# one contiguous "FUNCTION DEFINITIONS" section, numbered by the order its
# *category* first runs (0_ = generic utility, 1_ = data loading, 2_ = data
# preparation, 3_ = descriptive summary, 4_ = per-panel result assembly, i.e.
# the steps of the main panel loop, 5_ = the one pooled statistical test, run
# once after the main loop). The "MAIN SCRIPT" section after it only calls
# these functions; it defines no functions of its own.
# ===========================================================================

# ===========================================================================
# CONSTANTS
# ===========================================================================
STAGES <- paste0("T", 1:6)
# Assumed stage cadence, not a measured variable (see header).
WEEKS_ASSUMED <- c(1, 5, 9, 13, 17, 21)
TYPE_LABELS <- c(MG = "Milk kefir (grain)", WG = "Water kefir (grain)",
                 ML = "Milk kefir (liquid)", WL = "Water kefir (liquid)")
TYPE_SHORT <- c(MG = "Milk grain", WG = "Water grain", ML = "Milk liquid", WL = "Water liquid")
SRC_SHORT  <- c(MetaCache = "MetaCache", SUPER_FOCUS = "SUPER-FOCUS", volatiles = "Volatiles")

# ===========================================================================
# FUNCTION DEFINITIONS
# ===========================================================================

# --- 0_x: generic utilities --------------------------------------------------
# --- Small shared helpers ----------------------------------------------------
S_0_find_root <- function() {
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
S_0_normalize_id <- function(x) {
  x <- basename(gsub("\\\\", "/", trimws(as.character(x))))
  x <- sub("_host_removed_R[0-9]\\.fastq\\.gz$", "", x)
  gsub("-", "_", x, fixed = TRUE)
}
S_0_read_table <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
S_0_require_columns <- function(x, cols) {
  if (!all(cols %in% names(x))) stop("Missing columns: ", paste(setdiff(cols, names(x)), collapse = ", "))
}
# Strict numeric conversion: workbook blanks are non-detections (zero), but
# unexpected text, negative values and infinities are errors, never clipped.
S_0_numeric_matrix <- function(x, blanks_zero = FALSE) {
  raw <- as.matrix(x)
  m <- matrix(suppressWarnings(as.numeric(raw)), nrow(raw), ncol(raw), dimnames = dimnames(raw))
  blank <- is.na(raw) | trimws(raw) == ""
  if (any(is.na(m) & !blank)) stop("Non-numeric abundance values")
  if (blanks_zero) m[blank] <- 0
  if (any(!is.finite(m)) || any(m < 0)) stop("Missing, negative or infinite abundance values")
  m
}

# --- 1_x: data loading --------------------------------------------------------
S_1_load_taxonomy <- function(root, md) {
  x <- S_0_read_table(file.path(root, "data", "04_metacache_total_full_profile.csv"))
  S_0_require_columns(x, "sample_id")
  keys <- S_0_normalize_id(x$sample_id)
  keep <- keys %in% md$key
  if (anyDuplicated(keys[keep])) stop("Duplicate normalized MetaCache sample IDs")
  m <- S_0_numeric_matrix(x[keep, setdiff(names(x), "sample_id"), drop = FALSE], blanks_zero = TRUE)
  rownames(m) <- keys[keep]
  list(mat = m, notes = "")
}
S_1_load_superfocus <- function(root, md) {
  # Distributed copy first. A researcher holding the full internal tree can
  # point KEFIR4ALL_FUNC_DIR at it instead; no private path is assumed.
  folder <- Sys.getenv("KEFIR4ALL_FUNC_DIR", "")
  paths <- c(file.path(root, "data", "functional", "superfocus_p3.csv"),
             if (nzchar(folder))
               c(file.path(folder, "superfocus_p3.csv"),
                 file.path(folder, "Temp_4Liam", "05_short_read_functional_profiling",
                           "05_superfocus", "superfocus_p3.csv")) else character(0),
             file.path(root, "data", "metabolomics", "superfocus_p3.csv"))
  path <- paths[file.exists(paths)][1]
  if (is.na(path)) stop("SUPER-FOCUS input missing: expected ",
                        file.path(root, "data", "functional", "superfocus_p3.csv"),
                        " or KEFIR4ALL_FUNC_DIR")
  x <- S_0_read_table(path)
  S_0_require_columns(x, c("sample_id", "Subsystem Level 3", "reads_%"))
  x$key <- S_0_normalize_id(x$sample_id)
  x <- x[x$key %in% md$key, , drop = FALSE]
  if (!nrow(x)) stop("No SUPER-FOCUS sample IDs match metadata")
  x$value <- as.vector(S_0_numeric_matrix(x["reads_%"]))
  x$feature <- x[["Subsystem Level 3"]]
  if (anyNA(x$feature) || any(!nzchar(trimws(x$feature)))) stop("Missing SUPER-FOCUS feature labels")
  # Duplicate sample/pathway entries are summed; absent pathways are zero.
  m <- as.matrix(xtabs(value ~ key + feature, data = x))
  list(mat = m, notes = "")
}
# --- Volatile workbooks -----------------------------------------------------
# The volatile workbooks hold one column per compound measurement, with two
# kinds of column name that we use:
#   participant visit      "09_T2"      -> participant ID9, visit T2
#   laboratory reference   "FU_T1_0AN"  -> starter sample FU, used as T0
# The reference names also carry a _T1_/_T2_ token (an artefact of how the
# starters were sequenced); every such column is pooled into that type's T0
# reference, which is why the loader returns its own metadata.
# Returns list(mat, md); md plays the role of the study metadata for this panel
# because the reference rows are not listed in the sample metadata.
S_1_load_volatile_with_references <- function(root, meta, type) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("install readxl for the volatile workbooks")
  file <- if (type == "ML") "Final results table all milk kefirs 0913-001.xlsx" else
    "Final results table all water kefirs 0913-002.xlsx"
  sheet <- as.data.frame(readxl::read_excel(file.path(root, "data", file),
                                            .name_repair = "minimal"))
  if ("Name" %in% names(sheet)) names(sheet)[names(sheet) == "Name"] <- "Compound"
  S_0_require_columns(sheet, "Compound")

  column_names <- toupper(gsub("__", "_", names(sheet), fixed = TRUE))
  is_visit <- grepl("^0*[0-9]{1,3}_T[1-6]$", column_names)
  is_reference <- grepl("^[A-Z]+_T[1-6]_{1,2}0A[EN]T?$", column_names)
  if (!any(is_visit | is_reference)) stop("no visit or reference columns in ", file)

  # compounds are columns in the workbook, so transpose to samples in rows
  sel <- which(is_visit | is_reference)     # sheet columns we use, in sheet order
  values <- t(S_0_numeric_matrix(sheet[, sel, drop = FALSE], blanks_zero = TRUE))
  colnames(values) <- make.unique(as.character(sheet$Compound))

  # keep the columns that match a participant or a study starter of this kefir
  # type; every kept row needs a unique key and a participant label
  rows <- list(); keys <- character(0)
  participants <- character(0); stages <- character(0)
  for (k in seq_along(sel)) {               # k = row of values, j = column of sheet
    j <- sel[k]
    label <- column_names[j]
    if (is_visit[j]) {
      participant <- paste0("ID", as.integer(sub("_T[1-6]$", "", label)))
      stage <- sub("^.*_", "", label)
      hit <- meta$type == type & meta$participant == participant & meta$Stage == stage
      if (!any(hit)) next
      key <- meta$key[which(hit)[1]]
    } else {
      starter <- sub("_T[1-6].*$", "", label)
      hit <- meta$type == type & meta$participant == starter & meta$Stage == "T0"
      if (!any(hit)) next
      # reference rows are keyed REF_<type>_<starter>_<n> so that they stay unique
      key <- paste0("REF_", type, "_", starter, "_", length(keys) + 1L)
      participant <- paste0(starter, "_", length(keys) + 1L)
      stage <- "T0"
    }
    rows[[length(rows) + 1L]] <- values[k, , drop = FALSE]
    keys <- c(keys, key)
    participants <- c(participants, participant)
    stages <- c(stages, stage)
  }
  if (!length(rows)) stop("no volatile columns matched the metadata for ", type)

  mat <- do.call(rbind, rows)
  rownames(mat) <- keys
  md <- data.frame(key = keys, participant = participants, type = type,
                   Stage = stages, stringsAsFactors = FALSE)
  list(mat = mat, md = md)
}

# --- 2_x: data preparation -----------------------------------------------------
# Distance matrix builder. Bray-Curtis for compositional data (relative=TRUE:
# MetaCache, SUPER-FOCUS); Euclidean for measured concentrations (relative=FALSE:
# volatiles). In both cases profiles must be finite, nonnegative, nonempty.
S_2_dist_matrix <- function(m, method = c("bray", "euclidean")) {
  method <- match.arg(method)
  if (!nrow(m) || any(!is.finite(m)) || any(m < 0) || any(rowSums(m) <= 0))
    stop("Distance requires finite nonnegative, nonempty sample profiles")
  if (method == "bray") {
    if (!requireNamespace("vegan", quietly = TRUE)) stop("Install vegan for Bray-Curtis distances")
    d <- as.matrix(vegan::vegdist(m, method = "bray"))
    diag(d) <- 0
  } else {
    # Euclidean on the raw concentrations; sqrt(2) is the max distance for
    # 0/1 profiles, but volatile concentrations are unbounded.
    d <- as.matrix(dist(m, method = "euclidean"))
    diag(d) <- 0
  }
  d
}

S_2_prepare_panel <- function(m, md, type, relative) {
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
  method <- if (relative) "bray" else "euclidean"
  list(mat = m, md = meta, excluded = exclusions, dist = S_2_dist_matrix(m, method))
}

S_2_build_intervals <- function(dist, md) {
  rows <- list()
  for (id in unique(md$participant)) {
    idx <- which(md$participant == id)
    at <- idx[match(STAGES, md$Stage[idx])]
    for (s in 1:5) {
      if (anyNA(at[c(s, s + 1L)])) next  # never bridge missing visits
      i <- at[s]; j <- at[s + 1L]
      rows[[length(rows) + 1L]] <- data.frame(participant = id,
        from = STAGES[s], to = STAGES[s + 1L], assumed_week_from = WEEKS_ASSUMED[s],
        assumed_week_to = WEEKS_ASSUMED[s + 1L],
        distance = dist[cbind(i, j)], stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(participant = character(), from = character(),
    to = character(), assumed_week_from = numeric(), assumed_week_to = numeric(),
    distance = numeric()))
  do.call(rbind, rows)
}

# --- Measure 1: distance of every visit to the laboratory T0 reference ------
# mat      : abundance or compound matrix, samples in rows (names = sample keys)
# meta     : key / participant / type / Stage table covering rownames(mat)
# type     : kefir type to analyse ("MG", "WG", "ML", "WL")
# relative : TRUE converts each sample to proportions before the distance
#            (MetaCache, SUPER-FOCUS); FALSE keeps the measured values, which
#            are already comparable across samples (volatiles)
#
# Each visit is compared to the T0 library of ITS OWN participant's
# fermentation substrate (PARTICIPANT_SUBSTRATE, e.g. a participant who
# fermented skim milk is compared to the SK T0, not an arbitrary one) - see the
# "Substrate-matched T0 reference" section above for why this matters and how
# the codes are derived. A visit falls back to preferred_participant's T0 (or,
# failing that, the first T0 row for this type in metadata order - the
# original analysis's selection) when its participant's substrate is unknown
# or that substrate's T0 has no profile in this data source's matrix (volatile
# reference rows are keyed "<starter>_<n>", so the code is matched by prefix,
# not exact string). MetaCache and SUPER-FOCUS visit profiles are row
# proportions, so everything (visits and every candidate reference alike) is
# rescaled together before any distance is computed. Returns one row per
# visit: panel, participant, Stage, distance, reference_participant_used,
# matched_substrate; the default (fallback) reference key/participant/whether
# it was itself the preferred one are attached as attributes for the caller to
# log (see reference_rows below) - per-visit substrate-match detail lives in
# the returned data frame itself, written to <panel>_reference_distances.csv.
S_2_reference_distances <- function(mat, meta, type, relative, preferred_participant = NA_character_,
                                participant_substrate = NULL) {
  pos <- match(rownames(mat), meta$key)
  hits <- !is.na(pos) & meta$type[pos] == type
  if (!any(hits)) stop("no samples of kefir type ", type)
  x <- mat[hits, , drop = FALSE]
  info <- meta[pos[hits], , drop = FALSE]

  # Bray-Curtis cannot place an empty profile on a distance scale
  nonempty <- rowSums(x) > 0
  x <- x[nonempty, , drop = FALSE]
  info <- info[nonempty, , drop = FALSE]
  if (relative) x <- x / rowSums(x)

  is_reference <- info$Stage == "T0"
  if (!any(is_reference)) stop("no laboratory T0 reference for ", type)
  visits <- !is_reference
  if (!any(visits)) stop("no longitudinal visits for ", type)
  visit_idx <- which(visits)
  ref_idx_all <- which(is_reference)

  # default/fallback reference: preferred_participant's T0 if it has a profile
  # here, else the first T0 row for this type in metadata order
  used_preferred <- FALSE
  default_idx <- NA_integer_
  if (!is.na(preferred_participant)) {
    is_preferred <- info$participant == preferred_participant |
      grepl(paste0("^", preferred_participant, "_[0-9]+$"), info$participant)
    if (any(is_reference & is_preferred)) {
      default_idx <- which(is_reference & is_preferred)[1]
      used_preferred <- TRUE
    }
  }
  if (is.na(default_idx)) {
    default_idx <- ref_idx_all[1]
    if (!is.na(preferred_participant))
      message("No T0 profile for preferred reference participant ", preferred_participant,
              " (", type, "); falling back to ", info$participant[default_idx])
  }

  # substrate code each candidate T0 row itself represents (its own
  # participant field IS the code, e.g. "SK"; volatile reference rows are
  # keyed "SK_<n>", hence the suffix strip)
  t0_code <- sub("_[0-9]+$", "", info$participant[ref_idx_all])
  t0_by_code <- split(ref_idx_all, t0_code)

  ref_for_visit <- rep(default_idx, length(visit_idx))
  matched_substrate <- rep(FALSE, length(visit_idx))
  if (!is.null(participant_substrate)) {
    codes <- unname(participant_substrate[info$participant[visit_idx]])
    hit <- !is.na(codes) & codes %in% names(t0_by_code)
    ref_for_visit[hit] <- vapply(t0_by_code[codes[hit]], `[`, integer(1), 1)
    matched_substrate[hit] <- TRUE
  }

  # distance every visit profile to its own assigned reference, computed once
  # over the combined matrix. Bray-Curtis for compositional profiles;
  # Euclidean for measured volatile concentrations.
  combined <- rbind(x[visit_idx, , drop = FALSE], x[ref_idx_all, , drop = FALSE])
  method <- if (relative) "bray" else "euclidean"
  d <- S_2_dist_matrix(combined, method)
  n_visits <- length(visit_idx)
  ref_pos_in_combined <- n_visits + match(ref_for_visit, ref_idx_all)
  dist_to_own_ref <- d[cbind(seq_len(n_visits), ref_pos_in_combined)]

  out <- data.frame(panel = type,
             participant = info$participant[visit_idx],
             Stage = info$Stage[visit_idx],
             distance = dist_to_own_ref,
             reference_participant_used = info$participant[ref_for_visit],
             matched_substrate = matched_substrate,
             stringsAsFactors = FALSE)
  attr(out, "reference_key") <- info$key[default_idx]
  attr(out, "reference_participant") <- info$participant[default_idx]
  attr(out, "used_preferred") <- used_preferred
  out
}

# --- 3_x: descriptive summary -------------------------------------------------
# Plain available-case read of the raw reference-distance signal: mean/median
# per stage using every profile observed at that stage, no complete-
# trajectory filter.
S_3_stage_summary <- function(references) {
  do.call(rbind, lapply(STAGES, function(s) {
    v <- references$distance[references$Stage == s]
    v <- v[is.finite(v)]
    week <- WEEKS_ASSUMED[match(s, STAGES)]
    if (!length(v)) return(data.frame(Stage = s, week = week,
      n = 0L, mean = NA_real_, median = NA_real_))
    data.frame(Stage = s, week = week, n = length(v),
      mean = mean(v), median = stats::median(v))
  }))
}

# --- 4_x: per-panel result assembly ------------------------------------------
# One step of the main panel loop per function below, so the loop itself
# reads as a short list of steps instead of one long inlined block.

# Liquids (ML/WL): distance to the T0 reference matching the participant's
# own fermentation substrate - the liquid is fermented from the grain in a
# specific medium, so a different medium is a genuinely different starting
# point. Grains (MG/WG): a single consistent reference regardless of the
# participant's later medium choice - the grain culture itself does not
# depend on medium (MG has only one lab replicate anyway; WG's two, FL and
# AP, are independent grain preparations rather than medium-specific
# variants). Falls back to the single default reference when substrate
# matching is off or a participant's substrate is unknown/unavailable.
# Writes <panel>_reference_distances.csv; returns the visit-level distances
# and the one-row default-reference summary for community_changes_v2_t0_reference.csv.
S_4_reference_step <- function(loaded, panel_meta, spec, preferred_participant,
                                participant_substrate, out_dir) {
  references <- S_2_reference_distances(loaded$mat, panel_meta, spec$type,
                                    spec$relative, preferred_participant, participant_substrate)
  row <- data.frame(
    panel = spec$key, source = spec$source, type = spec$type,
    default_reference_key = attr(references, "reference_key"),
    default_reference_participant = attr(references, "reference_participant"),
    used_preferred_default = attr(references, "used_preferred"),
    n_visits = nrow(references),
    n_substrate_matched = sum(references$matched_substrate),
    stringsAsFactors = FALSE)
  write.csv(references, file.path(out_dir, paste0(spec$key, "_reference_distances.csv")),
            row.names = FALSE)
  list(references = references, row = row)
}

# Long-format participant distances for the figure (reference per stage,
# interval at the interval midpoint week), and the peak-stage tally: which
# visit sits furthest from the starter, per participant (>=4 visits; ties -
# none observed at the participant level - resolve to the earliest).
S_4_figures_and_peak <- function(references, intervals, spec) {
  fig_reference <- data.frame(
    source = spec$source, type = spec$type, measure = "reference",
    week = WEEKS_ASSUMED[match(references$Stage, STAGES)],
    distance = references$distance, participant = references$participant,
    stringsAsFactors = FALSE)
  fig_interval <- data.frame(
    source = spec$source, type = spec$type, measure = "interval",
    week = (intervals$assumed_week_from + intervals$assumed_week_to) / 2,
    distance = intervals$distance, participant = intervals$participant,
    stringsAsFactors = FALSE)

  ref_mean <- stats::aggregate(distance ~ participant + Stage,
                               data = references, FUN = mean)
  peak <- unlist(lapply(split(ref_mean, ref_mean$participant), function(z) {
    if (nrow(z) < 4L) return(NULL)
    z$Stage[which.max(z$distance)]
  }), use.names = FALSE)
  peak_row <- NULL
  if (length(peak)) {
    tb <- table(factor(peak, levels = STAGES))
    peak_row <- data.frame(
      source = spec$source, type = spec$type, panel = spec$key,
      stage = names(tb), n_peak = as.integer(tb), n_eligible = length(peak),
      stringsAsFactors = FALSE)
  }
  list(figures = list(fig_reference, fig_interval), peak = peak_row)
}

# --- 5_x: the one pooled statistical test ------------------------------------
# See header "THE ONE STATISTICAL TEST" for the full rationale. d must have
# columns participant, type, assumed_week_from, assumed_week_to, distance,
# already restricted to the taxonomic (MetaCache) panels and pooled across
# MG/WG/ML/WL.
S_5_segmented <- function(d, B = 9999L, seed = 20261022L + 201L, knot = 9) {
  d$trajectory <- interaction(d$participant, d$type, drop = TRUE)
  keep_traj <- names(which(table(d$trajectory) >= 2))
  d <- d[d$trajectory %in% keep_traj, , drop = FALSE]
  d$trajectory <- droplevels(d$trajectory)

  d$midpoint_week <- (d$assumed_week_from + d$assumed_week_to) / 2
  d$early_segment <- pmin(d$midpoint_week, knot) - knot
  d$later_segment <- pmax(d$midpoint_week - knot, 0)

  X <- model.matrix(~ early_segment + later_segment + trajectory, data = d)
  y <- d$distance
  cluster <- d$participant
  n <- nrow(X); k <- ncol(X)
  clusters <- unique(cluster)
  G <- length(clusters)
  XtX_inv <- solve(crossprod(X))

  # participant-clustered CR1 variance of a fitted beta, given its residuals
  cr1_vcov <- function(u) {
    meat <- matrix(0, k, k)
    for (g in clusters) {
      idx <- cluster == g
      s <- crossprod(X[idx, , drop = FALSE], u[idx])
      meat <- meat + tcrossprod(s)
    }
    adj <- (G / (G - 1)) * ((n - 1) / (n - k))
    adj * (XtX_inv %*% meat %*% XtX_inv)
  }

  beta <- XtX_inv %*% crossprod(X, y)
  resid_u <- as.numeric(y - X %*% beta)
  V <- cr1_vcov(resid_u)

  i_early <- which(colnames(X) == "early_segment")
  i_later <- which(colnames(X) == "later_segment")
  contrast <- rep(0, k); contrast[i_later] <- 1; contrast[i_early] <- -1
  slope_change <- as.numeric(contrast %*% beta)
  se_change <- sqrt(as.numeric(t(contrast) %*% V %*% contrast))
  t_obs <- slope_change / se_change
  tcrit <- stats::qt(0.975, df = G - 1)

  # null-restricted fit (one shared slope) supplies the wild-bootstrap DGP
  combined <- X[, i_early] + X[, i_later]
  Xr <- cbind(X[, -c(i_early, i_later), drop = FALSE], combined_segment = combined)
  beta_r <- solve(crossprod(Xr), crossprod(Xr, y))
  fitted_r <- as.numeric(Xr %*% beta_r)
  resid_r <- y - fitted_r

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  weights <- matrix(sample(c(-1, 1), G * B, replace = TRUE), nrow = G, ncol = B,
                     dimnames = list(clusters, NULL))
  w_per_obs <- weights[match(cluster, clusters), , drop = FALSE]
  Ystar <- fitted_r + resid_r * w_per_obs
  Beta_star <- XtX_inv %*% crossprod(X, Ystar)
  slope_change_star <- as.numeric(contrast %*% Beta_star)
  Resid_star <- Ystar - X %*% Beta_star
  se_star <- vapply(seq_len(B), function(b)
    sqrt(as.numeric(t(contrast) %*% cr1_vcov(Resid_star[, b]) %*% contrast)), numeric(1))
  t_star <- slope_change_star / se_star
  p_wild <- (sum(abs(t_star) >= abs(t_obs), na.rm = TRUE) + 1) / (B + 1)

  se_early <- sqrt(V[i_early, i_early]); se_later <- sqrt(V[i_later, i_later])
  data.frame(
    n_intervals = n, n_participants = G, n_trajectories = nlevels(d$trajectory),
    early_slope = unname(beta[i_early]),
    early_ci_low = unname(beta[i_early] - tcrit * se_early),
    early_ci_high = unname(beta[i_early] + tcrit * se_early),
    later_slope = unname(beta[i_later]),
    later_ci_low = unname(beta[i_later] - tcrit * se_later),
    later_ci_high = unname(beta[i_later] + tcrit * se_later),
    slope_change = slope_change, t_obs = t_obs, change_p_wild = p_wild, B = B)
}

# ===========================================================================
# MAIN SCRIPT
# ===========================================================================

# --- Setup --------------------------------------------------
root <- S_0_find_root()
out_dir <- file.path(root, "output", "04_taxonomic_profiling", "community_changes_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Trial order of the visits and the consecutive pairs between them
VISITS <- STAGES                                                    # T1..T6
INTERVALS <- paste(VISITS[-length(VISITS)], VISITS[-1], sep = "-")  # T1-T2..T5-T6

# --- Metadata --------------------------------------------------------------
# The metadata has one row per sequenced library, so the same sample can appear
# several times; collapse it to one row per sample key. Only study samples
# (kefir types MG/WG/ML/WL, stages T0-T6) are relevant to this script.
metadata <- S_0_read_table(file.path(root, "data", "kefir4all_sample_metadata_v2.csv"))
S_0_require_columns(metadata, c("merge_column", "Sample", "kefir type", "Stage"))
metadata$key <- S_0_normalize_id(metadata$merge_column)
metadata$type <- trimws(metadata[["kefir type"]])
metadata$Stage <- toupper(trimws(metadata$Stage))
metadata$participant <- toupper(trimws(as.character(metadata$Sample)))

# Participant labels are spelled inconsistently ("5", "ID5", "ID05") -> "ID5"
numeric_label <- grepl("^(ID)?0*[0-9]+$", metadata$participant)
metadata$participant[numeric_label] <- paste0(
  "ID", as.integer(sub("^ID", "", metadata$participant[numeric_label])))

candidate <- metadata$type %in% names(TYPE_LABELS) &
             metadata$Stage %in% c("T0", VISITS) &
             !is.na(metadata$participant) & nzchar(metadata$participant)

# A handful of merge_column keys have mate reads (R1/R2 of the same physical
# library) whose metadata disagrees on participant and/or kefir type - e.g. one
# T0 library is labelled participant FG / type WG on one mate and Medium
# control / type WL on the other. Which label is correct cannot be resolved
# from the metadata alone, so any such key is dropped entirely rather than
# silently keeping whichever mate wins the row-order dedup below.
signature <- paste(metadata$participant, metadata$type, metadata$Stage, sep = "|")
conflicted <- vapply(split(signature[candidate], metadata$key[candidate]),
                     function(z) length(unique(z)) > 1, logical(1))
conflicted_keys <- names(conflicted)[conflicted]
if (length(conflicted_keys))
  message("Excluding ", length(conflicted_keys),
          " sample key(s) with conflicting metadata: ", paste(conflicted_keys, collapse = ", "))
metadata <- metadata[!(metadata$key %in% conflicted_keys), , drop = FALSE]

study_sample <- metadata$type %in% names(TYPE_LABELS) &
                metadata$Stage %in% c("T0", VISITS) &
                !is.na(metadata$participant) & nzchar(metadata$participant)
metadata <- metadata[study_sample & !duplicated(metadata$key), , drop = FALSE]

# Fallback T0 reference participant per kefir type, used when a visit's own
# participant has no fermentation-substrate match available (see below).
# Picked from the raw metadata, independent of which data source's matrix
# happens to carry it, so MetaCache and SUPER-FOCUS fall back to the same
# physical starter for a given type instead of each independently picking
# whichever T0 replicate is first in *their own* sample matrix.
T0_REFERENCE_PARTICIPANT <- c(MG = NA_character_, WG = NA_character_, ML = NA_character_, WL = NA_character_)
for (ty in names(T0_REFERENCE_PARTICIPANT)) {
  rows <- metadata[metadata$Stage == "T0" & metadata$type == ty, , drop = FALSE]
  if (nrow(rows)) T0_REFERENCE_PARTICIPANT[[ty]] <- rows$participant[1]
}

# --- Substrate-matched T0 reference (see header for the full rationale) -----
# Lab code -> fermentation substrate, verified against the metadata's own
# kefir-type field for each code (note: 04_substrate_comparison.R's own
# lab_labels table has FL and FU swapped - trust this mapping, not that one).
SUBSTRATE_CODE <- c(
  "Cow milk (whole)" = "FU", "Cow milk (low fat)" = "LO", "Cow milk (skim)" = "SK", "Goat" = "GO",
  "Brown sugar-None" = "BR", "White sugar-dried apricot" = "AP", "White sugar-dried fig" = "FG",
  "White sugar-dried fig and fresh lemon" = "FL")
cs_private_path <- file.path(root, "data", "private", "Citizen Scientist metadata_v8.csv")
if (!file.exists(cs_private_path))
  stop("This script requires the private citizen-scientist metadata file (participant fermentation\n",
       "substrate) which is not distributed with the public repository.\n",
       "Missing: ", cs_private_path, "\nSee data/private/README.md for details.")
cs_meta <- S_0_read_table(cs_private_path)
S_0_require_columns(cs_meta, c("ID", "category_confirmed"))
cs_participant <- toupper(trimws(as.character(cs_meta$ID)))
numeric_cs <- grepl("^(ID)?0*[0-9]+$", cs_participant)
cs_participant[numeric_cs] <- paste0("ID", as.integer(sub("^ID", "", cs_participant[numeric_cs])))
cs_substrate_code <- unname(SUBSTRATE_CODE[trimws(cs_meta$category_confirmed)])
unmapped <- unique(trimws(cs_meta$category_confirmed)[is.na(cs_substrate_code) & nzchar(trimws(cs_meta$category_confirmed))])
if (length(unmapped))
  message("Substrate categories with no known T0 code (falling back to the default reference): ",
          paste(unmapped, collapse = ", "))
PARTICIPANT_SUBSTRATE <- setNames(cs_substrate_code, cs_participant)
PARTICIPANT_SUBSTRATE <- PARTICIPANT_SUBSTRATE[!duplicated(names(PARTICIPANT_SUBSTRATE)) & !is.na(PARTICIPANT_SUBSTRATE)]

# --- Panels -----------------------------------------------------------------
# One panel is one data source x one kefir type. relative = TRUE converts
# samples to proportions before the distances (MetaCache, SUPER-FOCUS); volatile
# concentrations are already comparable and stay as measured.
panels <- rbind(
  data.frame(key = paste("MetaCache", names(TYPE_LABELS), sep = "_"),
             source = "MetaCache", type = names(TYPE_LABELS), relative = TRUE,
             stringsAsFactors = FALSE),
  data.frame(key = paste("SUPER_FOCUS", names(TYPE_LABELS), sep = "_"),
             source = "SUPER_FOCUS", type = names(TYPE_LABELS), relative = TRUE,
             stringsAsFactors = FALSE),
  data.frame(key = c("volatiles_ML", "volatiles_WL"), source = "volatiles",
             type = c("ML", "WL"), relative = FALSE, stringsAsFactors = FALSE))

# Each data source is read once. A source that cannot be read only skips its own
# panels; those are listed with the error message in the manifest.
sources <- list(
  MetaCache = try(S_1_load_taxonomy(root, metadata), silent = TRUE),
  SUPER_FOCUS = try(S_1_load_superfocus(root, metadata), silent = TRUE),
  volatiles_ML = try(S_1_load_volatile_with_references(root, metadata, "ML"), silent = TRUE),
  volatiles_WL = try(S_1_load_volatile_with_references(root, metadata, "WL"), silent = TRUE))

manifest <- data.frame(panel = character(0), status = character(0),
                       reason = character(0), stringsAsFactors = FALSE)
peak_rows <- list()      # peak-stage (argmax distance to starter) tallies
stage_rows <- list()     # descriptive stage-distance summaries
figure_rows <- list()    # long-format participant distances, for the figure
reference_rows <- list() # which T0 library was used as the reference, per panel
metacache_interval_rows <- list() # the 4 taxonomic panels' intervals, for S_5_segmented()

for (i in seq_len(nrow(panels))) {
  spec <- panels[i, , drop = FALSE]
  loaded <- sources[[if (spec$source == "volatiles") spec$key else spec$source]]
  if (inherits(loaded, "try-error")) {
    manifest <- rbind(manifest, data.frame(
      panel = spec$key, status = "SKIPPED",
      reason = conditionMessage(attr(loaded, "condition")),
      stringsAsFactors = FALSE))
    next
  }
  # the matrix panels share the study metadata; a volatile panel brings its own
  # because its reference rows are not listed in the metadata
  panel_meta <- if (is.null(loaded$md)) metadata else loaded$md

  outcome <- try({
    # distance between each participant's own consecutive visits
    prepared <- S_2_prepare_panel(loaded$mat, panel_meta, spec$type,
                              relative = spec$relative)
    intervals <- S_2_build_intervals(prepared$dist, prepared$md)
    intervals$panel <- spec$key
    intervals$type <- spec$type
    write.csv(intervals, file.path(out_dir, paste0(spec$key, "_interval_distances.csv")),
              row.names = FALSE)
    if (spec$source == "MetaCache") metacache_interval_rows[[spec$key]] <- intervals

    # distance of every visit to its T0 reference (substrate-matched for
    # liquids - see S_4_reference_step())
    ref <- S_4_reference_step(loaded, panel_meta, spec, T0_REFERENCE_PARTICIPANT[[spec$type]],
                          if (spec$type %in% c("ML", "WL")) PARTICIPANT_SUBSTRATE else NULL,
                          out_dir)
    references <- ref$references
    reference_rows[[length(reference_rows) + 1L]] <- ref$row

    # figure data and the peak-stage tally
    fp <- S_4_figures_and_peak(references, intervals, spec)
    figure_rows[[length(figure_rows) + 1L]] <- fp$figures[[1]]
    figure_rows[[length(figure_rows) + 1L]] <- fp$figures[[2]]
    if (!is.null(fp$peak)) peak_rows[[length(peak_rows) + 1L]] <- fp$peak

    # descriptive stage-distance summary (every available profile, no
    # complete-trajectory filter)
    ss <- S_3_stage_summary(references)
    ss$panel <- spec$key
    stage_rows[[length(stage_rows) + 1L]] <- ss[, c("panel", setdiff(names(ss), "panel"))]
  }, silent = TRUE)

  failed <- inherits(outcome, "try-error")
  manifest <- rbind(manifest, data.frame(
    panel = spec$key,
    status = if (failed) "SKIPPED" else "OK",
    reason = if (failed) conditionMessage(attr(outcome, "condition")) else "",
    stringsAsFactors = FALSE))
}

# --- The one pooled statistical test (see header) ----------------------------
# Requires all four taxonomic panels to have loaded; a partial pool (e.g. one
# type skipped) would silently understate the trajectory count and is not run.
if (length(metacache_interval_rows) == 4L) {
  segmented <- S_5_segmented(do.call(rbind, metacache_interval_rows))
  write.csv(segmented, file.path(out_dir, "community_changes_v2_segmented_slope_test.csv"),
            row.names = FALSE)
  print(segmented, row.names = FALSE)
} else {
  message("Skipping the pooled taxonomic slope-change test: only ",
          length(metacache_interval_rows), " of 4 MetaCache panels loaded.")
}

# --- T0 reference audit trail (which library was used as "the" starter) -----
if (length(reference_rows)) {
  write.csv(do.call(rbind, reference_rows),
            file.path(out_dir, "community_changes_v2_t0_reference.csv"), row.names = FALSE)
}

# --- Peak-stage table (when is each participant furthest from the starter?) ---
if (length(peak_rows)) {
  write.csv(do.call(rbind, peak_rows),
            file.path(out_dir, "community_changes_v2_peak_stage.csv"), row.names = FALSE)
}

# --- Descriptive stage-distance summary --------------------------------------
if (length(stage_rows)) {
  write.csv(do.call(rbind, stage_rows),
            file.path(out_dir, "community_changes_v2_stage_distance_summary.csv"), row.names = FALSE)
}

# --- Manifest and methods wording -------------------------------------------
write.csv(manifest, file.path(out_dir, "community_changes_v2_manifest.csv"),
          row.names = FALSE)

writeLines(c(
  "Shared laboratory T0 reference distances and all five consecutive participant intervals.",
  "T0 is not a participant-specific baseline.",
  "This is a deliberately descriptive counterpart to community_changes.R, with one specific exception (see below): no other significance test is computed anywhere in this script. Both scripts now use the full rank-qualified MetaCache profile; formal-test results must be regenerated from that input before interpretation.",
  "POOLED TAXONOMIC SLOPE-CHANGE TEST (the one exception): community_changes_v2_segmented_slope_test.csv pools the four taxonomic panels (MetaCache MG/WG/ML/WL) and fits consecutive-visit Bray-Curtis distance against a segmented function of the interval midpoint week, with a knot at week 9 (the project's own pre-existing hypothesised water-kefir peak week, not chosen from this dataset's medians). Every participant x type trajectory gets its own fixed-effect intercept; trajectories with fewer than two usable intervals are dropped. Uncertainty is clustered by participant (not participant x type, since one participant can contribute more than one kefir type), using CR1 standard errors and a null-imposed participant-clustered wild-cluster bootstrap (Cameron, Gelbach and Miller 2008) for change_p_wild. This tests whether the later slope differs from the early slope; it does NOT test whether the later slope is zero, so it cannot by itself demonstrate a plateau - a shallower but still nonzero later trend is not ruled out, and this result says nothing about whether any individual panel, or the functional/volatile data, behaves the same way. No multiple-testing correction across a wider family is attempted here.",
  "T0 REFERENCE (SUBSTRATE-MATCHED FOR LIQUIDS): ML/WL visits are compared to the T0 library of their own participant's fermentation substrate, since the liquid is fermented from the grain in a specific medium; MG/WG visits use a single consistent reference regardless of medium, since the grain culture itself does not depend on it (MG has only one lab replicate; WG's two, FL and AP, are independent grain preparations, not medium-specific variants). T0 lab codes are substrate-specific: milk FU=whole/full-fat, GO=goat, SK=skim, LO=low-fat; water AP=white sugar+dried apricot, BR=brown sugar, FG=white sugar+dried fig, FL=white sugar+dried fig and fresh lemon (confirmed against the actual kefir-type tagging in the metadata, correcting a swapped FL/FU pair documented elsewhere in this repo). Participant substrate comes from data/private/Citizen Scientist metadata_v8.csv (category_confirmed), a private file not distributed with the public repository - this script shares that dependency with 04_community_stability.R and 04_substrate_comparison.R rather than running from public data alone. A visit falls back to a single default T0 per type (the first T0 row in metadata order, matching Figure_3.R/04_community_stability.R, not a pooled average) when substrate matching is off (grains) or a participant's substrate is unknown or that substrate's T0 has no profile in the relevant data source's matrix (liquids); MetaCache/SUPER-FOCUS/volatiles share the same default participant per type rather than each independently picking whichever T0 replicate their own matrix carries first, wherever possible. community_changes_v2_t0_reference.csv records the default reference and substrate-match coverage (n_substrate_matched / n_visits) per panel; <panel>_reference_distances.csv records the actual reference used and match status per visit (reference_participant_used, matched_substrate).",
  "EXCLUDED CONFLICTING SAMPLE: one merge_column key (TG_FG_WLC_S410) has mate reads with disagreeing metadata (participant FG/type WG on one mate, Medium control/type WL on the other) and is dropped entirely rather than guessing which label is correct; any future such conflicts are detected the same way and logged to the console.",
  "MetaCache and SUPER-FOCUS profiles are converted to row proportions before Bray-Curtis distances; volatile concentrations use Euclidean distances.",
  "MetaCache uses rank-qualified rank::name features; unclassified reads are excluded and samples below 10,000 total reads are omitted during profile construction. Rows are normalized over classified reads. The build audit records total/classified reads, classified fraction, raw-report checksum, and exclusion reason; differing rank resolution can contribute to dissimilarity.",
  "Volatile distances use raw measured concentrations and are not numerically comparable to taxonomic or functional distances.",
  "Distances are averaged per participant per stage/interval so partial coverage still contributes.",
  "ATTRITION: later visits have fewer participants than earlier ones (see n in community_changes_v2_stage_distance_summary.csv and the per-week point counts in the figure); an apparent front-loaded pattern of change can partly reflect this rather than purely within-participant dynamics.",
  "PEAK STAGE: community_changes_v2_peak_stage.csv tallies which visit has each participant's largest observed distance to the starter; it is a descriptive summary for participants with at least four visits, not a test of a peak or plateau.",
  "STAGE-DISTANCE SUMMARY: community_changes_v2_stage_distance_summary.csv gives the mean/median distance-to-starter per stage using every available profile (no complete-trajectory filter). This is a descriptive available-case read, not a paired test: each stage can draw on different participants, so a difference between stages can partly reflect attrition rather than within-participant change.",
  paste("Run directory:", out_dir)), file.path(out_dir, "METHODS.txt"))

print(manifest, row.names = FALSE)

# --- Figure: reference vs interval distance per panel, and peak-stage bars --
# Builds directly from the figure_rows/peak_rows already collected above (not
# by re-reading the CSVs just written), so this single file goes from raw
# data straight to the finished figures. Panels with no data (e.g. volatiles
# has no MG/WG panels) are simply omitted rather than drawn as empty
# placeholder cells. No significance annotation is drawn - this figure is
# descriptive only (see header).
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 for the figures")
figures_dir <- file.path(root, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

type_order <- c("MG", "WG", "ML", "WL")
panel_order <- rbind(
  data.frame(source = "MetaCache", type = type_order, stringsAsFactors = FALSE),
  data.frame(source = "SUPER_FOCUS", type = type_order, stringsAsFactors = FALSE),
  data.frame(source = "volatiles", type = c("ML", "WL"), stringsAsFactors = FALSE))
panel_order$panel <- ifelse(panel_order$source == "volatiles",
                            paste0("volatiles_", panel_order$type),
                            paste(panel_order$source, panel_order$type, sep = "_"))
panel_order$facet <- paste(TYPE_SHORT[panel_order$type], SRC_SHORT[panel_order$source], sep = " · ")

if (length(figure_rows)) {
  d <- do.call(rbind, figure_rows)
  d$measure <- factor(d$measure, levels = c("reference", "interval"))
  d$panel <- ifelse(d$source == "volatiles", paste0("volatiles_", d$type), paste(d$source, d$type, sep = "_"))
  d$facet <- factor(panel_order$facet[match(d$panel, panel_order$panel)], levels = panel_order$facet)

  # participant-mean +/- standard error per panel/measure/week
  mean_d <- aggregate(distance ~ panel + measure + week, data = d, FUN = mean, na.rm = TRUE)
  se_d   <- aggregate(distance ~ panel + measure + week, data = d,
                      FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  names(mean_d)[4] <- "mean"; names(se_d)[4] <- "se"
  m <- merge(mean_d, se_d, by = c("panel", "measure", "week"))
  m$lo <- m$mean - m$se; m$hi <- m$mean + m$se
  m$facet <- factor(panel_order$facet[match(m$panel, panel_order$panel)], levels = panel_order$facet)

  p <- ggplot2::ggplot(d, ggplot2::aes(week, distance)) +
    ggplot2::geom_point(ggplot2::aes(color = measure), alpha = 0.18, size = 0.7, show.legend = FALSE) +
    ggplot2::geom_errorbar(data = m, ggplot2::aes(y = mean, ymin = lo, ymax = hi, color = measure),
                           width = 0.6, show.legend = FALSE) +
    ggplot2::geom_line(data = m, ggplot2::aes(y = mean, color = measure, linetype = measure)) +
    ggplot2::geom_point(data = m, ggplot2::aes(y = mean, color = measure), size = 1.3) +
    ggplot2::facet_wrap(~ facet, scales = "free_y", ncol = 4) +
    ggplot2::scale_color_manual(values = c(reference = "#2c7bb6", interval = "#d7191c"),
                                labels = c("starter (T0)", "previous visit"), name = "Distance to") +
    ggplot2::scale_linetype_manual(values = c(reference = "solid", interval = "dashed"), guide = "none") +
    ggplot2::scale_x_continuous(breaks = unname(WEEKS_ASSUMED),
                                labels = paste0("wk", sprintf("%02d", WEEKS_ASSUMED))) +
    ggplot2::labs(x = "Time since start (week)", y = "Distance",
                 caption = "Points: individual participants; lines: mean +/- SE. Descriptive only - no significance test is computed in this script.") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(size = 6))

  ggplot2::ggsave(file.path(figures_dir, "community_changes_v2_figure.tiff"), p, width = 9, height = 6, dpi = 300, compression = "lzw")
  ggplot2::ggsave(file.path(figures_dir, "community_changes_v2_figure.pdf"), p, width = 9, height = 6)
  write.table(d, file.path(figures_dir, "community_changes_v2_figure_data.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- Peak-stage bar chart: when is each participant furthest from T0? -------
if (length(peak_rows)) {
  peak <- do.call(rbind, peak_rows)
  peak$panel <- paste(peak$source, peak$type, sep = "_")
  peak$stage <- factor(peak$stage, levels = STAGES)
  elig <- stats::aggregate(n_eligible ~ panel, data = peak, FUN = function(x) x[1])
  peak$facet <- paste0(peak$panel, "  (n = ", elig$n_eligible[match(peak$panel, elig$panel)], ")")

  p2 <- ggplot2::ggplot(peak, ggplot2::aes(stage, n_peak)) +
    ggplot2::geom_col(fill = "#2c7bb6") +
    ggplot2::facet_wrap(~ facet, scales = "free_y") +
    ggplot2::labs(x = "Peak visit (max distance to T0)", y = "Participants",
                 caption = "Peak-stage distribution: descriptive only - no single visit dominates; maxima scatter across wk01-wk21.") +
    ggplot2::theme_bw(base_size = 9)

  ggplot2::ggsave(file.path(figures_dir, "community_changes_v2_peak_stage.tiff"), p2, width = 9, height = 6, dpi = 300, compression = "lzw")
  ggplot2::ggsave(file.path(figures_dir, "community_changes_v2_peak_stage.pdf"), p2, width = 6, height = 4)
  write.table(peak, file.path(figures_dir, "community_changes_v2_peak_stage_data.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

message("Wrote figures to: ", figures_dir)
