# ===========================================================================
# community_changes.R
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Kefir4All is a citizen-science trial: participants were given a shared milk
# or water kefir starter culture and fermented it themselves for 21 weeks,
# sending in a sample at each of six visits (T1-T6, assumed to fall at weeks
# 1, 5, 9, 13, 17, 21 - the metadata records only the visit label, never a
# measured date). This script asks: how much does each participant's kefir
# community change over that trial, and does it change more at some points
# than others?
#
# It analyses ten "panels" - one per combination of data source and kefir
# type:
#   data sources: MetaCache (rank-qualified taxonomic features), SUPER-FOCUS (functional
#                 subsystems, level 3), volatiles (measured VOC concentrations)
#   kefir types:  MG/WG (milk/water kefir GRAIN), ML/WL (milk/water kefir
#                 LIQUID) - volatiles were only measured for the two liquid
#                 types, so there are 4 + 4 + 2 = 10 panels in total.
# Everything - data loading, statistics, and the summary figures - lives in
# this one file; there is no separate script to run afterwards.
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
# community_changes_t0_reference.csv (the per-panel default) and
# <panel>_reference_distances.csv (the actual reference per visit, and
# whether it was substrate-matched).
#
# STATISTICAL TESTS
# ------------------
# Both distance measures (reference and interval) get, per panel:
#   - an OMNIBUS test asking "do any of the sampling timepoints differ from
#     each other?" (perm_p) - a participant-blocked permutation test: values
#     are reshuffled only among a given participant's OWN groups, because
#     ordinary Kruskal-Wallis assumes independent groups and here the same
#     participant contributes to every group (see the comment above S_3_kw_H()
#     for the full reasoning). perm_p_balanced repeats the same test
#     restricted to participants with a complete trajectory, as a robustness
#     check.
#   - PAIRWISE Wilcoxon tests between every pair of timepoints, BH-corrected
#     within the panel (paired when the same participants appear in both
#     timepoints, unpaired otherwise).
# Separately, PERMANOVA (adonis2, participant-blocked permutation) runs
# directly on the underlying T1-T6 sample-by-sample distance matrix, testing
# whether composition differs by stage at all - a different, more complete
# question than the scalar reference-distance test above, since two visits
# equally far from T0 but different from each other would be invisible to
# that test. A betadisper/permutest check runs alongside it, because a
# significant PERMANOVA can reflect a difference in within-stage spread
# rather than a genuine shift in typical composition.
# Finally, a purely descriptive stage-distance summary (mean/median distance
# per stage, using every available profile) is reported - this is descriptive
# only, not a statistical test.
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
# Tables land under output/04_taxonomic_profiling/community_changes/:
#   <panel>_reference_distances.csv, <panel>_interval_distances.csv - raw distances
#   <panel>_<measure>_pairwise.csv         - pairwise Wilcoxon results
#   community_changes_omnibus.csv          - omnibus test result per panel/measure
#   community_changes_permanova.csv        - PERMANOVA + betadisper per panel
#   community_changes_t0_reference.csv     - which T0 sample was used, per panel
#   community_changes_peak_stage.csv       - per-participant peak-visit tally
#   community_changes_stage_distance_summary.csv - descriptive stage-distance summary
#   community_changes_manifest.csv         - which panels ran OK / were skipped, and why
#   METHODS.txt                             - full prose methods for every number above
# Figures land under figures/:
#   community_changes_figure.{tiff,pdf}        - reference vs interval distance, all panels
#   community_changes_peak_stage.{tiff,pdf}    - per-panel peak-visit bar charts
#
# USAGE
# -----
# Run from the kefir4all_review_history root (S_0_find_root() locates it by
# searching for data/kefir4all_sample_metadata_v2.csv):
#   Rscript scripts/r_scripts/04_taxonomic_profiling/community_changes.R
# Or via the pinned, reproducible Docker image (see docker/README.md) - this
# is the environment that should be treated as authoritative, since it is the
# one guaranteed to have every package this script actually needs:
#   docker run --rm -v "$(pwd):/home/rstudio/kefir4all" \
#     -w /home/rstudio/kefir4all ghcr.io/liamhwalsh/kefir4all:4.4.2-slim \
#     /opt/R/4.4.2/bin/Rscript scripts/r_scripts/04_taxonomic_profiling/community_changes.R
#
# LAYOUT OF THIS FILE
# --------------------
# Function definitions and the script that applies them are kept in two
# separate, non-interleaved blocks below: every function (21 total) is
# defined once, in one contiguous "FUNCTION DEFINITIONS" section, numbered by
# the order its *category* first runs (0_ = generic utility, 1_ = data
# loading, 2_ = data preparation, 3_ = statistical analysis, 4_ = per-panel
# result assembly, i.e. the steps of the main panel loop; e.g. MetaCache
# is loaded before SUPER-FOCUS, so S_1_load_taxonomy before S_1_load_superfocus).
# The "MAIN SCRIPT" section after it only calls these functions; it defines
# no functions of its own.
# ===========================================================================

# ===========================================================================
# CONSTANTS
# ===========================================================================
STAGES <- paste0("T", 1:6)
# Assumed stage cadence, not a measured variable (see header).
WEEKS_ASSUMED <- c(1, 5, 9, 13, 17, 21)
PERM_B <- 2000L    # participant-blocked permutations for the omnibus robustness check
NPERM <- 4999L     # PERMANOVA permutations (adonis2), also participant-blocked
SEED <- 20260916L  # permutation seed (reproducibility)
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
        assumed_week_to = WEEKS_ASSUMED[s + 1L], window = if (s <= 2) "early" else "late",
        distance = dist[cbind(i, j)], stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(participant = character(), from = character(),
    to = character(), assumed_week_from = numeric(), assumed_week_to = numeric(),
    window = character(), distance = numeric()))
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

# --- 3_x: statistical analysis -------------------------------------------------
# --- Participant-blocked permutation omnibus ---------------------------------
# Ordinary Kruskal-Wallis assumes independent observations between groups, but
# here the same participant contributes a value to every group (every stage,
# or every interval), so the groups are not independent: a participant who is
# simply more (or less) divergent than average pulls their value up or down in
# every group at once, which can bias a plain KW p-value in either direction,
# not just conservatively.
#
# The Monte Carlo test reshuffles each participant's observed values among
# their observed groups and compares KW H statistics. Interpretation requires
# within-participant exchangeability. This is especially uncertain for adjacent
# interval distances, which share visits; treat those tests as exploratory.
S_3_kw_H <- function(m) {
  # m is participants (rows) x groups (columns). R flattens column-major, so
  # the group labels run per column: rep(each = nrow), not times - times would
  # interleave 1..G,1..G,... and mislabel every group.
  x <- as.vector(m)
  g <- factor(rep(seq_len(ncol(m)), each = nrow(m)))
  keep <- is.finite(x)
  x <- x[keep]; g <- droplevels(g[keep])
  N <- length(x)
  if (N < 3L || nlevels(g) < 2L) return(NA_real_)
  r <- rank(x)
  ng <- as.integer(table(g))
  rg <- as.numeric(tapply(r, g, mean))
  H <- 12 / (N * (N + 1)) * sum(ng * (rg - (N + 1) / 2)^2)
  ties <- table(x); ties <- ties[ties > 1L]
  if (length(ties)) H <- H / (1 - sum(ties^3 - ties) / (N^3 - N))
  H
}

S_3_blocked_kruskal_p <- function(long, B = PERM_B) {
  wide <- tidyr::pivot_wider(long[, c("participant", "group", "distance")],
                             names_from = group, values_from = distance)
  groups <- setdiff(names(wide), "participant")
  m <- as.matrix(wide[, groups, drop = FALSE])
  obs <- S_3_kw_H(m)
  if (!is.finite(obs)) return(NA_real_)
  # One reproducible stream for the whole script run: the seed is set once by
  # the caller (before the panel loop), never re-seeded here, so successive
  # panels/measures draw fresh permutations instead of replaying one sequence.
  # NA cells are structural missingness (participant lacks that visit), not
  # values to shuffle: only the observed cells move within each participant.
  cnt <- 0L
  for (b in seq_len(B)) {
    mp <- m
    for (i in seq_len(nrow(m))) {
      idx <- which(!is.na(m[i, ]))
      if (length(idx) > 1L) mp[i, idx] <- m[i, sample(idx)]
    }
    if (isTRUE(S_3_kw_H(mp) >= obs)) cnt <- cnt + 1L
  }
  (1 + cnt) / (B + 1)
}

# --- Participant-blocked permutation omnibus plus all pairwise comparisons ---
# long    : participant, group (stage or interval), distance
# panel   : panel name, used for file names and inside the result tables
# measure : "reference" (groups = T1..T6) or "interval" (groups = the five
#           consecutive visit pairs)
#
# Distances are averaged per participant and group, so participants with
# partial visit coverage still contribute. The omnibus is the
# participant-blocked permutation test (perm_p, with perm_p_balanced on the
# complete-trajectory subset - see S_3_blocked_kruskal_p()). n_min/n_max expose
# the between-group attrition (later visits typically have fewer
# participants than earlier ones). Each pair of groups is tested with a
# paired Wilcoxon when at least two participants overlap, or unpaired Wilcoxon
# only when groups are disjoint. A single overlapping participant is not tested.
# Pairwise results are BH-corrected within the panel. Returns omnibus and
# pairwise tables; the caller writes them so the omnibus p-values can be
# BH-adjusted across panels within each measure.
S_3_analyse_measure <- function(long, panel, measure) {
  unit <- if (measure == "reference") "stage" else "interval"
  groups <- if (measure == "reference") VISITS else INTERVALS

  # one distance per participant and group: average repeat samples, no drop
  long <- stats::aggregate(distance ~ participant + group, data = long,
                           FUN = mean, na.rm = TRUE)
  long$group <- factor(long$group, levels = groups)
  long$participant <- factor(long$participant)
  wide <- tidyr::pivot_wider(long[, c("participant", "group", "distance")],
                             names_from = group, values_from = distance)
  groups <- intersect(groups, names(wide))

  # group sizes expose the attrition that leaves the later visits unbalanced
  sizes <- vapply(groups, function(g) sum(is.finite(wide[[g]])), integer(1))
  balanced <- wide[stats::complete.cases(wide[, groups, drop = FALSE]), , drop = FALSE]

  omnibus <- data.frame(panel = panel, measure = measure,
                        n_participants = sum(rowSums(!is.na(wide[, groups, drop = FALSE])) > 0),
                        stringsAsFactors = FALSE)
  omnibus[[paste0("n_", unit, "s")]] <- length(groups)
  omnibus$n_min <- if (length(sizes)) min(sizes) else NA_integer_
  omnibus$n_max <- if (length(sizes)) max(sizes) else NA_integer_
  omnibus$n_balanced <- nrow(balanced)
  omnibus$perm_p <- NA_real_
  omnibus$perm_p_balanced <- NA_real_
  pairwise <- data.frame(panel = character(0), a = character(0), b = character(0),
                         n_a = integer(0), n_b = integer(0), n_paired = integer(0),
                         test_method = character(0), test_note = character(0),
                         p_value = numeric(0), p_BH = numeric(0),
                         stringsAsFactors = FALSE)

  # headline omnibus (see S_3_blocked_kruskal_p): the null is within-participant
  # exchangeability rather than independent samples
  if (length(groups) >= 2L && nrow(long) > 0L) {
    omnibus$perm_p <- S_3_blocked_kruskal_p(long)
    balanced_long <- long[long$participant %in% balanced$participant, , drop = FALSE]
    if (nrow(balanced) >= 2L) omnibus$perm_p_balanced <- S_3_blocked_kruskal_p(balanced_long)
  }

  # every group pair, BH-corrected within this panel
  if (length(groups) >= 2L) {
    pairwise <- do.call(rbind, lapply(utils::combn(groups, 2L, simplify = FALSE),
      function(pair) {
        a <- wide[[pair[1]]]; b <- wide[[pair[2]]]
        both <- is.finite(a) & is.finite(b)
        n_paired <- sum(both)
        test_method <- "not_testable"
        test_note <- "insufficient observations"
        if (n_paired >= 2L) {
          p <- stats::wilcox.test(a[both], b[both], paired = TRUE, exact = FALSE)$p.value
          test_method <- "paired_wilcoxon"
          test_note <- ""
        } else {
          one <- is.finite(a); two <- is.finite(b)
          if (n_paired == 0L && sum(one) >= 2L && sum(two) >= 2L) {
            p <- stats::wilcox.test(a[one], b[two], exact = FALSE)$p.value
            test_method <- "unpaired_wilcoxon"
            test_note <- ""
          } else {
            p <- NA_real_
            if (n_paired == 1L) test_note <- "one_shared_participant; independence not met"
          }
        }
        data.frame(panel = panel, a = pair[1], b = pair[2],
                   n_a = sum(is.finite(a)), n_b = sum(is.finite(b)),
                   n_paired = n_paired, test_method = test_method, test_note = test_note,
                   p_value = p, stringsAsFactors = FALSE)
      }))
    pairwise$p_BH <- stats::p.adjust(pairwise$p_value, method = "BH", n = nrow(pairwise))
  }
  names(pairwise)[2:3] <- paste0(unit, c("_1", "_2"))

  # the caller writes the files so that the omnibus p-values can first be
  # BH-adjusted across panels within each measure family
  list(omnibus = omnibus, pairwise = pairwise, measure = measure, panel = panel)
}

# --- PERMANOVA on the T1-T6 distance matrix ----------------------------------
# S_2_prepare_panel() already restricts m/dist to T1-T6 (STAGES excludes T0), so no
# extra filtering is needed. Tests whether composition differs by Stage using
# every available profile (not just complete trajectories), with permutations
# restricted within participant blocks so a participant appearing at several
# stages does not let the test treat those as independent draws.
#
# One follow-up alongside the categorical omnibus, using the same
# participant-blocked permutation scheme: betadisper/permutest, because
# PERMANOVA is sensitive to multivariate dispersion, not just location, so a
# significant result cannot be read as "typical composition differs by
# stage" without checking whether it is instead a spread difference.
# Stage is categorical so the test does not impose a linear temporal trend.
S_3_permanova_test <- function(prepared) {
  if (!requireNamespace("permute", quietly = TRUE)) stop("Install permute for PERMANOVA")
  n <- nrow(prepared$mat)
  np <- length(unique(prepared$md$participant))
  empty <- list(n = n, n_participants = np, R2 = NA_real_, F = NA_real_, p = NA_real_,
                betadisper_p = NA_real_)
  if (n < 10L || np < 3L) return(empty)
  d <- stats::as.dist(prepared$dist)
  dd <- data.frame(Stage = factor(prepared$md$Stage, levels = STAGES),
                   participant = prepared$md$participant)
  ctl <- permute::how(nperm = NPERM, blocks = dd$participant)
  res <- vegan::adonis2(d ~ Stage, data = dd, permutations = ctl)
  bd_p <- tryCatch({
    bd <- vegan::betadisper(d, dd$Stage)
    vegan::permutest(bd, permutations = ctl)$tab$`Pr(>F)`[1]
  }, error = function(e) NA_real_)
  list(n = n, n_participants = np, R2 = unname(res$R2[1]),
       F = unname(res$F[1]), p = unname(res$`Pr(>F)`[1]), betadisper_p = bd_p)
}

# --- Descriptive stage-distance summary --------------------------------------
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
# and the one-row default-reference summary for community_changes_t0_reference.csv.
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

# S_3_analyse_measure() returns the omnibus/pairwise tables for both distance
# measures; this writes the pairwise CSVs and hands back the two omnibus rows
# for the caller to collect (so p-values can later be BH-adjusted across panels).
S_4_omnibus_step <- function(references, intervals, spec, out_dir) {
  results <- list(
    S_3_analyse_measure(data.frame(participant = references$participant,
                               group = references$Stage,
                               distance = references$distance,
                               stringsAsFactors = FALSE),
                    spec$key, "reference"),
    S_3_analyse_measure(data.frame(participant = intervals$participant,
                               group = paste(intervals$from, intervals$to, sep = "-"),
                               distance = intervals$distance,
                               stringsAsFactors = FALSE),
                    spec$key, "interval"))
  lapply(results, function(res) {
    write.csv(res$pairwise,
              file.path(out_dir, paste0(res$panel, "_", res$measure, "_pairwise.csv")),
              row.names = FALSE)
    res$omnibus
  })
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

# PERMANOVA on the T1-T6 distance matrix: does composition differ by stage at
# all, using every available profile (participant-blocked permutation), plus
# the betadisper dispersion check (see S_3_permanova_test()). This uses the
# full distance matrix directly, rather than reducing each sample to a single
# distance-to-starter number, so two visits equally far from T0 but different
# from each other are caught here even though the reference-distance test
# above would treat them as identical.
S_4_permanova_step <- function(prepared, spec) {
  prm <- tryCatch(S_3_permanova_test(prepared), error = function(e)
    list(n = nrow(prepared$mat), n_participants = NA_integer_,
         R2 = NA_real_, F = NA_real_, p = NA_real_, betadisper_p = NA_real_))
  data.frame(
    panel = spec$key, n = prm$n, n_participants = prm$n_participants,
    R2 = prm$R2, F = prm$F, p_value = prm$p, betadisper_p = prm$betadisper_p,
    stringsAsFactors = FALSE)
}

# ===========================================================================
# MAIN SCRIPT
# ===========================================================================

# --- Setup --------------------------------------------------
root <- S_0_find_root()
source(file.path(root, "scripts", "r_scripts", "figure_style.R"))
kefir_use_style(base_size = 12)
out_dir <- file.path(root, "output", "04_taxonomic_profiling", "community_changes")
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
omnibus_rows <- list()   # collected so p-values can be BH-adjusted across panels
peak_rows <- list()      # peak-stage (argmax distance to starter) tallies
stage_rows <- list()     # descriptive stage-distance summaries
permanova_rows <- list() # PERMANOVA + betadisper
figure_rows <- list()    # long-format participant distances, for the figure
reference_rows <- list() # which T0 library was used as the reference, per panel

# Single reproducible RNG stream for every permutation test in this run; each
# S_3_analyse_measure() call continues the sequence rather than re-seeding.
set.seed(SEED)

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
    write.csv(intervals, file.path(out_dir, paste0(spec$key, "_interval_distances.csv")),
              row.names = FALSE)

    # distance of every visit to its T0 reference (substrate-matched for
    # liquids - see S_4_reference_step())
    ref <- S_4_reference_step(loaded, panel_meta, spec, T0_REFERENCE_PARTICIPANT[[spec$type]],
                          if (spec$type %in% c("ML", "WL")) PARTICIPANT_SUBSTRATE else NULL,
                          out_dir)
    references <- ref$references
    reference_rows[[length(reference_rows) + 1L]] <- ref$row

    # omnibus + pairwise tests for both distance measures
    for (om in S_4_omnibus_step(references, intervals, spec, out_dir))
      omnibus_rows[[length(omnibus_rows) + 1L]] <- om

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

    # PERMANOVA + betadisper on the T1-T6 distance matrix
    permanova_rows[[length(permanova_rows) + 1L]] <- S_4_permanova_step(prepared, spec)
  }, silent = TRUE)

  failed <- inherits(outcome, "try-error")
  manifest <- rbind(manifest, data.frame(
    panel = spec$key,
    status = if (failed) "SKIPPED" else "OK",
    reason = if (failed) conditionMessage(attr(outcome, "condition")) else "",
    stringsAsFactors = FALSE))
}

# --- Omnibus tables ----------------------------------------------------------
# BH is applied across the panels within each measure family (ten tests each),
# which the per-panel tests do not do on their own.
if (length(omnibus_rows)) {
  # Reference rows carry n_stages and interval rows n_intervals, so the collected
  # frames do not share an identical column set and base rbind() refuses to
  # combine them ("names do not match previous names"). Align the columns first,
  # padding with NA where a column does not apply to that measure.
  cols <- unique(unlist(lapply(omnibus_rows, names), use.names = FALSE))
  omnibus_rows <- lapply(omnibus_rows, function(d) {
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
    d[, cols, drop = FALSE]
  })
  omnibus_all <- do.call(rbind, omnibus_rows)
  omnibus_all$perm_p_BH <- stats::ave(omnibus_all$perm_p, omnibus_all$measure,
                                      FUN = function(p) stats::p.adjust(p, method = "BH", n = nrow(panels)))
  write.csv(omnibus_all, file.path(out_dir, "community_changes_omnibus.csv"),
            row.names = FALSE)
  for (i in seq_len(nrow(omnibus_all))) {
    one <- omnibus_all[i, , drop = FALSE]
    write.csv(one, file.path(out_dir, paste0(one$panel, "_", one$measure, "_omnibus.csv")),
              row.names = FALSE)
  }
}

# --- T0 reference audit trail (which library was used as "the" starter) -----
if (length(reference_rows)) {
  write.csv(do.call(rbind, reference_rows),
            file.path(out_dir, "community_changes_t0_reference.csv"), row.names = FALSE)
}

# --- Peak-stage table (when is each participant furthest from the starter?) ---
if (length(peak_rows)) {
  write.csv(do.call(rbind, peak_rows),
            file.path(out_dir, "community_changes_peak_stage.csv"), row.names = FALSE)
}

# --- Descriptive stage-distance summary --------------------------------------
if (length(stage_rows)) {
  write.csv(do.call(rbind, stage_rows),
            file.path(out_dir, "community_changes_stage_distance_summary.csv"), row.names = FALSE)
}

# --- PERMANOVA (primary community-composition test; BH across the ten panels) -
if (length(permanova_rows)) {
  permanova_all <- do.call(rbind, permanova_rows)
  permanova_all$p_BH <- stats::p.adjust(permanova_all$p_value, method = "BH", n = nrow(panels))
  write.csv(permanova_all, file.path(out_dir, "community_changes_permanova.csv"), row.names = FALSE)
}

# --- Manifest and methods wording -------------------------------------------
write.csv(manifest, file.path(out_dir, "community_changes_manifest.csv"),
          row.names = FALSE)

writeLines(c(
  "Shared laboratory T0 reference distances and all five consecutive participant intervals.",
  "T0 is not a participant-specific baseline.",
  "T0 REFERENCE (SUBSTRATE-MATCHED FOR LIQUIDS): ML/WL visits are compared to the T0 library of their own participant's fermentation substrate, since the liquid is fermented from the grain in a specific medium; MG/WG visits use a single consistent reference regardless of medium, since the grain culture itself does not depend on it (MG has only one lab replicate; WG's two, FL and AP, are independent grain preparations, not medium-specific variants). T0 lab codes are substrate-specific: milk FU=whole/full-fat, GO=goat, SK=skim, LO=low-fat; water AP=white sugar+dried apricot, BR=brown sugar, FG=white sugar+dried fig, FL=white sugar+dried fig and fresh lemon (confirmed against the actual kefir-type tagging in the metadata, correcting a swapped FL/FU pair documented elsewhere in this repo). Participant substrate comes from data/private/Citizen Scientist metadata_v8.csv (category_confirmed), a private file not distributed with the public repository - this script now shares that dependency with 04_community_stability.R and 04_substrate_comparison.R rather than running from public data alone. A visit falls back to a single default T0 per type (the first T0 row in metadata order, matching Figure_3.R/04_community_stability.R, not a pooled average) when substrate matching is off (grains) or a participant's substrate is unknown or that substrate's T0 has no profile in the relevant data source's matrix (liquids); MetaCache/SUPER-FOCUS/volatiles share the same default participant per type rather than each independently picking whichever T0 replicate their own matrix carries first, wherever possible. community_changes_t0_reference.csv records the default reference and substrate-match coverage (n_substrate_matched / n_visits) per panel; <panel>_reference_distances.csv records the actual reference used and match status per visit (reference_participant_used, matched_substrate).",
  "EXCLUDED CONFLICTING SAMPLE: one merge_column key (TG_FG_WLC_S410) has mate reads with disagreeing metadata (participant FG/type WG on one mate, Medium control/type WL on the other) and is dropped entirely rather than guessing which label is correct; any future such conflicts are detected the same way and logged to the console.",
  "MetaCache and SUPER-FOCUS profiles are converted to row proportions before Bray-Curtis distances; volatile concentrations use Euclidean distances.",
  "MetaCache uses rank-qualified rank::name features from the full profile; unclassified reads are excluded, and rows are normalized over classified reads. Rank-resolution differences between samples can contribute to profile dissimilarity. The build audit reports total and classified reads, classified fraction, raw-report filename checksum, and QC reason.",
  "Volatile distances use raw measured concentrations and are not numerically comparable to taxonomic or functional distances.",
  "OMNIBUS: perm_p is a participant-blocked Monte Carlo permutation test using the Kruskal-Wallis H statistic (2000 permutations; seed 20260916). It reshuffles each participant's observed values among groups and requires within-participant exchangeability. It allows incomplete coverage but does not remove distributional or time-dependence assumptions; interval results are exploratory because adjacent intervals share visits. perm_p_balanced repeats the test among participants covering every group; the complete-case sample can be too small to be informative.",
  "Reference tests compare distance-to-T0 distributions across wk01-wk21 using available participant means. Pairwise Wilcoxon is paired with at least two shared participants, unpaired only when groups are disjoint, and not run when exactly one participant overlaps. BH correction covers all planned pairs within each panel.",
  "Interval tests compare T1-T2, T2-T3, T3-T4, T4-T5, and T5-T6 distance distributions using available participant means. Pairwise tests use the same overlap rule. Because adjacent intervals share samples, their permutation omnibus is exploratory and does not establish a general early-versus-late effect.",
  "Distances are averaged per participant per stage/interval so partial coverage still contributes; n_participants reports distinct participants with at least one value; n_min/n_max and n_balanced report group sizes and the complete-trajectory subset.",
  "ATTRITION: group sizes vary by stage and interval; n_min/n_max and per-pair n report the available coverage. Available-case summaries can reflect both temporal change and changes in which participants were sampled.",
  "PAIRWISE TESTS: p_BH is corrected across all planned group pairs within each panel. Paired Wilcoxon is used with at least two shared participants; unpaired Wilcoxon only when groups are disjoint. A pair with exactly one shared participant is not tested because neither independence nor a paired estimate is supported.",
  "Omnibus p-values are uncorrected within the per-panel files; community_changes_omnibus.csv adds perm_p_BH across the ten panels of each measure family.",
  "PEAK STAGE: community_changes_peak_stage.csv tallies which visit has each participant's largest observed distance to the starter; it is a descriptive summary for participants with at least four visits, not a test of a peak or plateau.",
  "STAGE-DISTANCE SUMMARY: community_changes_stage_distance_summary.csv gives the mean/median distance-to-starter per stage using every available profile (no complete-trajectory filter). This is a descriptive available-case read, not a paired test: each stage can draw on different participants, so a difference between stages can partly reflect attrition rather than within-participant change.",
  paste0("PERMANOVA (community_changes_permanova.csv): adonis2 tests whether T1-T6 composition is associated with Stage, using all available profiles and permutations restricted within participant blocks. This is distinct from testing distance-to-T0, and is exploratory because inference requires within-participant exchangeability. p_BH adjusts across all ten planned panels. betadisper_p is an unadjusted dispersion diagnostic, not a definitive separation of location and spread effects."),
  "Non-significant PERMANOVA results are not evidence of stability (no equivalence margin).",
  paste("Run directory:", out_dir)), file.path(out_dir, "METHODS.txt"))

print(manifest, row.names = FALSE)

# --- Figure: reference vs interval distance per panel, and peak-stage bars --
# Builds directly from the figure_rows/omnibus_all/peak_rows already
# collected above (not by re-reading the CSVs just written), so this single
# file goes from raw data straight to the finished figures. Panels with no
# data (e.g. volatiles has no MG/WG panels) are simply omitted rather than
# drawn as empty placeholder cells.
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

  # annotation text: headline permutation p per panel/measure, from omnibus_all
  ann <- do.call(rbind, lapply(panel_order$panel, function(key) {
    ro <- if (exists("omnibus_all")) omnibus_all[omnibus_all$panel == key & omnibus_all$measure == "reference", ] else NULL
    io <- if (exists("omnibus_all")) omnibus_all[omnibus_all$panel == key & omnibus_all$measure == "interval", ] else NULL
    rp <- if (!is.null(ro) && nrow(ro)) ro$perm_p[1] else NA_real_
    ip <- if (!is.null(io) && nrow(io)) io$perm_p[1] else NA_real_
    data.frame(panel = key,
               label = sprintf("ref perm p = %.3g\nint perm p = %.3g%s", rp, ip, if (isTRUE(ip < 0.05)) " *" else ""),
               sig = isTRUE(ip < 0.05), stringsAsFactors = FALSE)
  }))
  ann <- merge(ann, aggregate(hi ~ panel, data = m, FUN = max, na.rm = TRUE), by = "panel")
  ann$facet <- factor(panel_order$facet[match(ann$panel, panel_order$panel)], levels = panel_order$facet)
  ann <- ann[!is.na(ann$facet), , drop = FALSE]

  p <- ggplot2::ggplot(d, ggplot2::aes(week, distance)) +
    ggplot2::geom_point(ggplot2::aes(color = measure), alpha = 0.18, size = 0.7, show.legend = FALSE) +
    ggplot2::geom_errorbar(data = m, ggplot2::aes(y = mean, ymin = lo, ymax = hi, color = measure),
                           width = 0.6, show.legend = FALSE) +
    ggplot2::geom_line(data = m, ggplot2::aes(y = mean, color = measure, linetype = measure)) +
    ggplot2::geom_point(data = m, ggplot2::aes(y = mean, color = measure), size = 1.3) +
    ggplot2::geom_text(data = ann, ggplot2::aes(x = 1, y = hi * 1.05, label = label,
                                                fontface = ifelse(sig, "bold", "plain")),
                       hjust = 0, vjust = 1, size = 3.2, color = "grey20", inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ facet, scales = "free_y", ncol = 4) +
    ggplot2::scale_color_manual(values = c(reference = "#2c7bb6", interval = "#d7191c"),
                                labels = c("starter (T0)", "previous visit"), name = "Distance to") +
    ggplot2::scale_linetype_manual(values = c(reference = "solid", interval = "dashed"), guide = "none") +
    ggplot2::scale_x_continuous(breaks = unname(WEEKS_ASSUMED),
                                labels = paste0("wk", sprintf("%02d", WEEKS_ASSUMED))) +
    ggplot2::labs(x = "Time since start (week)", y = "Distance",
                 caption = "Points: individual participants; lines: mean +/- SE. * participant-blocked permutation p < 0.05 (interval).") +
    kefir_theme(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 16, face = "bold"))

  ggplot2::ggsave(file.path(figures_dir, "community_changes_figure.tiff"), p, width = 15, height = 11, dpi = 300, compression = "lzw")
  ggplot2::ggsave(file.path(figures_dir, "community_changes_figure.pdf"), p, width = 15, height = 11)
  write.table(d, file.path(figures_dir, "community_changes_figure_data.tsv"),
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
                 caption = "Peak-stage distribution: no single visit dominates; maxima scatter across wk01-wk21.") +
    kefir_theme(base_size = 9)

  ggplot2::ggsave(file.path(figures_dir, "community_changes_peak_stage.tiff"), p2, width = 9, height = 6, dpi = 300, compression = "lzw")
  ggplot2::ggsave(file.path(figures_dir, "community_changes_peak_stage.pdf"), p2, width = 6, height = 4)
  write.table(peak, file.path(figures_dir, "community_changes_peak_stage_data.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

message("Wrote figures to: ", figures_dir)
