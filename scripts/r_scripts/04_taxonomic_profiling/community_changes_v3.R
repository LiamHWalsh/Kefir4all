# ===========================================================================
# community_changes_v3.R
# Describes participant trajectories; tests common stage effects separately.
# Run from the repository root in ghcr.io/liamhwalsh/kefir4all:4.4.2-slim:
# /opt/R/4.4.2/bin/Rscript scripts/r_scripts/04_taxonomic_profiling/community_changes_v3.R
# Optional KEFIR4ALL_COMMUNITY_TYPES: CSV with sample_id, community_type.
# Assignments must come from the existing community-type analysis, not this file.
# Output: output/04_taxonomic_profiling/community_changes_v3/ and figures/.
# ===========================================================================

STAGES <- paste0("T", 1:6)
# T0 is week 0; T1 is 7 days later. Later visit weeks are assumed from the study schedule.
WEEKS <- c(1, 5, 9, 13, 17, 21)
TYPES <- c("MG", "WG", "ML", "WL")
SEED <- 20260916L
NPERM <- 4999L
NBOOT <- 2000L
DEPTH_SENSITIVITY_CUTOFFS <- c(30000, 1000000)
READ_DEPTH_SENSITIVITY <- c(30000, 1000000)

find_root <- function() {
  starts <- c(Sys.getenv("KEFIR4ALL_ROOT"), getwd(), file.path(getwd(), "kefir4all_review_history"))
  for (start in starts[nzchar(starts)]) {
    if (file.exists(file.path(start, "data", "kefir4all_sample_metadata_v2.csv")))
      return(normalizePath(start, winslash = "/"))
  }
  stop("Run from the repository root or set KEFIR4ALL_ROOT.")
}

normalize_id <- function(value) {
  value <- basename(gsub("\\\\", "/", trimws(as.character(value))))
  value <- sub("_host_removed_R[0-9]\\.fastq\\.gz$", "", value)
  gsub("-", "_", value, fixed = TRUE)
}

participant_id <- function(value) {
  value <- toupper(trimws(as.character(value)))
  numeric <- !is.na(value) & grepl("^(ID)?0*[0-9]+$", value)
  value[numeric] <- paste0("ID", as.integer(sub("^ID", "", value[numeric])))
  value
}

read_table <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

require_columns <- function(table, columns) {
  if (!all(columns %in% names(table)))
    stop("Missing columns: ", paste(setdiff(columns, names(table)), collapse = ", "))
}

numeric_matrix <- function(table, blanks_zero = FALSE) {
  raw <- as.matrix(table)
  blank <- is.na(raw) | trimws(raw) == ""
  values <- matrix(suppressWarnings(as.numeric(raw)), nrow(raw), ncol(raw), dimnames = dimnames(raw))
  if (any(is.na(values) & !blank)) stop("Unexpected nonnumeric abundance")
  if (blanks_zero) values[blank] <- 0
  if (any(!is.finite(values)) || any(values < 0)) stop("Missing, negative, or infinite abundance")
  values
}

load_metadata <- function(path) {
  metadata <- read_table(path)
  require_columns(metadata, c("merge_column", "Sample", "kefir type", "Stage"))
  metadata$key <- normalize_id(metadata$merge_column)
  metadata$participant <- participant_id(metadata$Sample)
  metadata$type <- trimws(metadata[["kefir type"]])
  metadata$Stage <- toupper(trimws(metadata$Stage))
  candidate <- metadata$type %in% TYPES & metadata$Stage %in% c("T0", STAGES)
  # Check mate-read disagreements before dropping controls or other non-study rows.
  relevant <- metadata$key %in% metadata$key[candidate]
  signature <- paste(metadata$participant, metadata$type, metadata$Stage, sep = "|")
  conflict <- vapply(split(signature[relevant], metadata$key[relevant]),
                     function(values) length(unique(values)) > 1L, logical(1))
  bad_keys <- names(conflict)[conflict]
  invalid <- is.na(metadata$key) | !nzchar(metadata$key) |
    is.na(metadata$participant) | !nzchar(metadata$participant)
  excluded <- metadata[(metadata$key %in% bad_keys | invalid) & relevant,
                       c("key", "participant", "type", "Stage"), drop = FALSE]
  excluded$reason <- ifelse(excluded$key %in% bad_keys, "conflicting_metadata", "missing_identifier")
  metadata <- metadata[candidate & !invalid & !(metadata$key %in% bad_keys), , drop = FALSE]
  metadata <- metadata[!duplicated(metadata$key), , drop = FALSE]
  list(md = metadata, excluded = excluded)
}

load_substrates <- function(path) {
  if (!file.exists(path)) stop("Private substrate metadata required: ", path)
  table <- read_table(path)
  require_columns(table, c("ID", "category_confirmed"))
  # Map the recorded substrate labels explicitly; do not infer codes from names.
  codes <- c("Cow milk (whole)" = "FU", "Cow milk (low fat)" = "LO",
             "Cow milk (skim)" = "SK", "Goat" = "GO", "Brown sugar-None" = "BR",
             "White sugar-dried apricot" = "AP", "White sugar-dried fig" = "FG",
             "White sugar-dried fig and fresh lemon" = "FL")
  table$participant <- participant_id(table$ID)
  table$category <- trimws(table$category_confirmed)
  keep <- !is.na(table$participant) & nzchar(table$participant)
  table <- table[keep, , drop = FALSE]
  conflict <- vapply(split(table$category, table$participant),
                     function(values) length(unique(values)) > 1L, logical(1))
  if (any(conflict)) stop("Conflicting private substrate categories for a participant")
  table <- table[!duplicated(table$participant), , drop = FALSE]
  table$code <- unname(codes[table$category])
  list(codes = setNames(table$code, table$participant),
       unmapped = table[is.na(table$code), c("participant", "category"), drop = FALSE])
}

load_taxonomy <- function(path, md) {
  table <- read_table(path)
  require_columns(table, "sample_id")
  keys <- normalize_id(table$sample_id)
  keep <- keys %in% md$key
  if (anyDuplicated(keys[keep])) stop("Duplicate normalized MetaCache keys")
  mat <- numeric_matrix(table[keep, setdiff(names(table), "sample_id"), drop = FALSE], TRUE)
  rownames(mat) <- keys[keep]
  list(mat = mat, md = md)
}

load_superfocus <- function(path, md) {
  table <- read_table(path)
  require_columns(table, c("sample_id", "Subsystem Level 3", "reads_%"))
  table$key <- normalize_id(table$sample_id)
  table <- table[table$key %in% md$key, , drop = FALSE]
  if (!nrow(table)) stop("No matched SUPER-FOCUS profiles")
  table$feature <- table[["Subsystem Level 3"]]
  if (anyNA(table$feature) || any(!nzchar(trimws(table$feature)))) stop("Missing pathway labels")
  table$value <- as.vector(numeric_matrix(table["reads_%"]))
  list(mat = as.matrix(stats::xtabs(value ~ key + feature, table)), md = md)
}

find_superfocus_file <- function(root, functional_dir = Sys.getenv("KEFIR4ALL_FUNC_DIR", "")) {
  candidates <- c(file.path(root, "data", "functional", "superfocus_p3.csv"),
    if (nzchar(functional_dir)) c(file.path(functional_dir, "superfocus_p3.csv"),
      file.path(functional_dir, "Temp_4Liam", "05_short_read_functional_profiling",
        "05_superfocus", "superfocus_p3.csv")),
    file.path(root, "data", "metabolomics", "superfocus_p3.csv"))
  available <- candidates[file.exists(candidates)]
  if (length(available)) available[1] else candidates[1]
}

load_volatiles <- function(path, md, type) {
  table <- as.data.frame(readxl::read_excel(path, .name_repair = "minimal"))
  if ("Name" %in% names(table)) names(table)[names(table) == "Name"] <- "Compound"
  require_columns(table, "Compound")
  labels <- toupper(gsub("__", "_", names(table), fixed = TRUE))
  # Workbooks store participant visits and lab starters as sample columns.
  visit <- grepl("^0*[0-9]{1,3}_T[1-6]$", labels)
  reference <- grepl("^[A-Z]+_T[1-6]_{1,2}0A[EN]T?$", labels)
  selected <- which(visit | reference)
  if (!length(selected)) stop("No recognized volatile columns")
  values <- t(numeric_matrix(table[, selected, drop = FALSE], TRUE))
  compounds <- trimws(as.character(table$Compound))
  if (anyNA(compounds) || any(!nzchar(compounds)) || anyDuplicated(compounds))
    stop("Missing or duplicate volatile compound labels; resolve before analysis")
  colnames(values) <- compounds
  rows <- list()
  info <- list()
  excluded <- list()
  for (index in seq_along(selected)) {
    column <- selected[index]
    label <- labels[column]
    if (visit[column]) {
      participant <- participant_id(sub("_T[1-6]$", "", label))
      stage <- sub("^.*_", "", label)
      hit <- which(md$type == type & md$participant == participant & md$Stage == stage)
      if (length(hit) != 1L) {
        excluded[[length(excluded) + 1L]] <- data.frame(column = names(table)[column],
          reason = if (length(hit)) "ambiguous_metadata_visit" else "no_matching_study_visit")
        next
      }
      key <- md$key[hit]
    } else {
      starter <- sub("_T[1-6].*$", "", label)
      if (!any(md$type == type & md$participant == starter & md$Stage == "T0")) {
        excluded[[length(excluded) + 1L]] <- data.frame(column = names(table)[column], reason = "no_matching_study_T0")
        next
      }
      key <- paste("REF", type, label, column, sep = "_")
      participant <- paste(starter, column, sep = "_")
      stage <- "T0"
    }
    rows[[length(rows) + 1L]] <- values[index, , drop = FALSE]
    info[[length(info) + 1L]] <- data.frame(key = key, participant = participant,
      type = type, Stage = stage, workbook_column = names(table)[column])
  }
  if (!length(rows)) stop("No volatile profiles matched study metadata")
  metadata <- do.call(rbind, info)
  mat <- do.call(rbind, rows)
  rownames(mat) <- metadata$key
  list(mat = mat, md = metadata, unmatched_columns = if (length(excluded)) do.call(rbind, excluded)
    else data.frame(column = character(), reason = character()))
}

prepare_panel <- function(loaded, type, relative) {
  if (anyNA(match(rownames(loaded$mat), loaded$md$key))) stop("Unmatched matrix keys")
  # Filter metadata in its original order (used for the default T0) and align
  # the abundance matrix to the selected keys without reshaping the matrix.
  md <- loaded$md |>
    dplyr::filter(.data$key %in% rownames(loaded$mat),
                  .data$type == .env$type,
                  .data$Stage %in% c("T0", STAGES))
  mat <- loaded$mat[match(md$key, rownames(loaded$mat)), , drop = FALSE]
  if (!nrow(mat) || any(!is.finite(mat)) || any(mat < 0)) stop("Invalid panel matrix")
  # Exclude every copy of a duplicate participant visit rather than choose one.
  ambiguous_keys <- md |>
    dplyr::filter(.data$Stage %in% STAGES) |>
    dplyr::add_count(.data$participant, .data$Stage, name = "n_visits") |>
    dplyr::filter(.data$n_visits > 1L) |>
    dplyr::pull("key")
  zero_keys <- md$key[rowSums(mat) == 0]
  # Empty profiles have undefined compositional distances; retain why each row drops.
  excluded <- md |>
    dplyr::mutate(reason = dplyr::case_when(
      .data$key %in% ambiguous_keys ~ "ambiguous_duplicate_visit",
      .data$key %in% zero_keys ~ "zero_profile",
      TRUE ~ NA_character_)) |>
    dplyr::filter(!is.na(.data$reason)) |>
    dplyr::select("key", "participant", "type", "Stage", "reason")
  keep <- !(md$key %in% excluded$key)
  md <- md[keep, , drop = FALSE]
  mat <- mat[keep, , drop = FALSE]
  if (!nrow(mat) || !any(md$Stage %in% STAGES)) stop("No usable visits")
  # Bray-Curtis compares composition; Euclidean volatiles retain concentration scale.
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (relative) mat <- mat / rowSums(mat)
  stopifnot(identical(rownames(mat), md$key))
  list(mat = mat, md = md, excluded = excluded, relative = relative,
       distance = as.matrix(if (relative) vegan::vegdist(mat, method = "bray") else stats::dist(mat)))
}

reference_distances <- function(panel, default_key, default_participant, substrates) {
  md <- panel$md
  starters <- which(md$Stage == "T0")
  if (!length(starters)) stop("No usable T0 reference")
  preferred <- starters[md$key[starters] == default_key]
  if (!length(preferred)) preferred <- starters[sub("_[0-9]+$", "", md$participant[starters]) == default_participant]
  # Use the preferred starter if available; otherwise use the first usable T0.
  default <- if (length(preferred)) preferred[1] else starters[1]
  codes <- sub("_[0-9]+$", "", md$participant[starters])
  visits <- which(md$Stage %in% STAGES)
  rows <- lapply(visits, function(index) {
    code <- if (is.null(substrates)) NA_character_ else unname(substrates[md$participant[index]])
    candidates <- if (is.na(code)) integer() else starters[codes == code]
    reference <- if (length(candidates)) candidates[1] else default
    data.frame(key = md$key[index], participant = md$participant[index], Stage = md$Stage[index],
      distance = panel$distance[index, reference], reference_key = md$key[reference],
      reference_participant = md$participant[reference], matched_substrate = length(candidates) > 0L,
      n_reference_candidates = if (length(candidates)) length(candidates) else length(preferred),
      fallback_reason = if (length(candidates)) "" else if (is.null(substrates)) "grain_consistent_reference"
        else if (is.na(code)) "unknown_substrate" else "matched_T0_unavailable",
      default_reference_key = md$key[default], preferred_default_available = length(preferred) > 0L)
  })
  do.call(rbind, rows)
}

build_intervals <- function(panel) {
  md <- panel$md
  rows <- list()
  for (participant in unique(md$participant[md$Stage %in% STAGES])) {
    indices <- which(md$participant == participant & md$Stage %in% STAGES)
    at <- indices[match(STAGES, md$Stage[indices])]
    for (stage in seq_len(5L)) {
      pair <- at[c(stage, stage + 1L)]
      # A missing visit breaks the pair; never bridge across it.
      if (anyNA(pair)) next
      rows[[length(rows) + 1L]] <- data.frame(participant = participant,
        from_key = md$key[pair[1]], to_key = md$key[pair[2]], from = STAGES[stage], to = STAGES[stage + 1L],
        week = mean(WEEKS[c(stage, stage + 1L)]),
        distance = panel$distance[pair[1], pair[2]])
    }
  }
  if (!length(rows)) return(data.frame(participant = character(), from_key = character(),
    to_key = character(), from = character(), to = character(), week = numeric(), distance = numeric()))
  do.call(rbind, rows)
}

participant_summary <- function(panel, intervals) {
  md <- panel$md
  rows <- lapply(unique(md$participant[md$Stage %in% STAGES]), function(participant) {
    indices <- which(md$participant == participant & md$Stage %in% STAGES)
    indices <- indices[order(match(md$Stage[indices], STAGES))]
    # Net endpoint distance is not a path length; the interval sum uses only observed adjacent pairs.
    first <- indices[1]
    last <- tail(indices, 1)
    distances <- intervals$distance[intervals$participant == participant]
    data.frame(participant = participant, n_visits = length(indices), n_observed_intervals = length(distances),
      first_key = md$key[first], last_key = md$key[last], first_stage = md$Stage[first], last_stage = md$Stage[last],
      net_first_last_distance = if (length(indices) >= 2L) panel$distance[first, last] else NA_real_,
      mean_interval_distance = if (length(distances)) mean(distances) else NA_real_,
      median_interval_distance = if (length(distances)) stats::median(distances) else NA_real_,
      sum_observed_interval_distances = if (length(distances)) sum(distances) else NA_real_,
      complete_trajectory = all(STAGES %in% md$Stage[indices]))
  })
  do.call(rbind, rows)
}

mean_interval <- function(values) {
  values <- values[is.finite(values)]
  result <- c(n = length(values), mean = if (length(values)) mean(values) else NA_real_,
              lower = NA_real_, upper = NA_real_)
  if (length(values) >= 5L && length(unique(values)) > 1L) {
    # Resample participants (one value each), not individual repeated visits.
    fit <- boot::boot(values, function(data, indices) mean(data[indices]), R = NBOOT)
    interval <- boot::boot.ci(fit, type = "perc")$percent[4:5]
    result[c("lower", "upper")] <- interval
  }
  result
}

stage_test <- function(panel) {
  keep <- panel$md$Stage %in% STAGES
  md <- panel$md[keep, , drop = FALSE]
  repeated <- names(which(table(md$participant) >= 2L))
  keep <- keep & panel$md$participant %in% repeated
  md <- panel$md[keep, , drop = FALSE]
  result <- data.frame(n = nrow(md), n_participants = length(repeated), n_stages = length(unique(md$Stage)),
    R2_total = NA_real_, partial_R2 = NA_real_, F = NA_real_, p_value = NA_real_,
    dispersion_p = NA_real_, n_permutations = 0L, reason = "", dispersion_reason = "")
  if (length(repeated) < 3L || length(unique(md$Stage)) < 2L) {
    result$reason <- "insufficient_repeated_participants_or_stages"
    return(result)
  }
  dd <- data.frame(participant = factor(md$participant), Stage = factor(md$Stage))
  design <- stats::model.matrix(~ participant + Stage, dd)
  if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) {
    result$reason <- "aliased_stage_or_no_residual_degrees_of_freedom"
    return(result)
  }
  distance <- stats::as.dist(panel$distance[keep, keep, drop = FALSE])
  # Control for participant-level differences and restrict permutations within participants.
  # Inference still assumes exchangeability across each participant's observed visits.
  control <- permute::how(nperm = NPERM, blocks = dd$participant)
  permutations <- permute::shuffleSet(nrow(dd), nset = NPERM, control = control)
  result$n_permutations <- nrow(permutations)
  fit <- vegan::adonis2(distance ~ participant + Stage, data = dd,
    permutations = permutations, by = "margin", sqrt.dist = panel$relative)
  result$R2_total <- fit["Stage", "R2"]
  result$partial_R2 <- fit["Stage", "SumOfSqs"] / (fit["Stage", "SumOfSqs"] + fit["Residual", "SumOfSqs"])
  result$F <- fit["Stage", "F"]
  result$p_value <- fit["Stage", "Pr(>F)"]
  dispersion_distance <- if (panel$relative) sqrt(distance) else distance
  if (any(table(dd$Stage) < 2L)) {
    result$dispersion_reason <- "stage_has_fewer_than_two_samples"
  } else {
    diagnostic <- tryCatch({
      dispersion <- vegan::betadisper(dispersion_distance, dd$Stage, type = "median", bias.adjust = TRUE)
      list(p = vegan::permutest(dispersion, permutations = permutations)$tab[1, "Pr(>F)"], reason = "")
    }, error = function(error) list(p = NA_real_, reason = conditionMessage(error)))
    result$dispersion_p <- diagnostic$p
    result$dispersion_reason <- diagnostic$reason
  }
  result
}

summarize_rank_sensitivity <- function(loaded, type, panel_key, default_key,
                                       default_participant, substrates) {
  feature_ranks <- sub("::.*$", "", colnames(loaded$mat))
  ranks <- sort(unique(feature_ranks))
  summaries <- list()
  manifest <- list()

  for (rank in ranks) {
    rank_loaded <- list(mat = loaded$mat[, feature_ranks == rank, drop = FALSE], md = loaded$md)
    result <- tryCatch({
      panel <- prepare_panel(rank_loaded, type, relative = TRUE)
      references <- reference_distances(panel, default_key, default_participant, substrates)
      intervals <- build_intervals(panel)
      reference_summary <- references |>
        dplyr::group_by(.data$Stage) |>
        dplyr::summarise(n = dplyr::n(), mean = mean(.data$distance),
          median = stats::median(.data$distance), .groups = "drop") |>
        dplyr::transmute(panel = panel_key, rank = rank, measure = "reference",
          group = .data$Stage, n = .data$n, mean = .data$mean, median = .data$median)
      interval_summary <- intervals |>
        dplyr::group_by(.data$from, .data$to) |>
        dplyr::summarise(n = dplyr::n(), mean = mean(.data$distance),
          median = stats::median(.data$distance), .groups = "drop") |>
        dplyr::transmute(panel = panel_key, rank = rank, measure = "interval",
          group = paste(.data$from, .data$to, sep = "-"), n = .data$n,
          mean = .data$mean, median = .data$median)
      list(status = "OK", reason = "", n_samples = nrow(panel$mat),
        summaries = dplyr::bind_rows(reference_summary, interval_summary))
    }, error = function(error) list(status = "SKIPPED", reason = conditionMessage(error),
      n_samples = 0L, summaries = NULL))
    manifest[[length(manifest) + 1L]] <- data.frame(panel = panel_key, rank = rank,
      status = result$status, reason = result$reason, n_samples = result$n_samples,
      n_features = if (result$status == "OK") sum(feature_ranks == rank) else 0L)
    if (!is.null(result$summaries))
      summaries[[length(summaries) + 1L]] <- result$summaries
  }
  list(summary = if (length(summaries)) dplyr::bind_rows(summaries) else data.frame(),
       manifest = dplyr::bind_rows(manifest))
}

test_depth_sensitivity <- function(loaded, type, read_depth, cutoff) {
  keep <- is.finite(read_depth) & read_depth >= cutoff
  result <- tryCatch({
    panel <- prepare_panel(list(mat = loaded$mat[keep, , drop = FALSE], md = loaded$md),
                           type, relative = TRUE)
    model <- stage_test(panel)
    data.frame(minimum_read_depth = cutoff, n_matrix_samples = sum(keep),
      n_samples = model$n, n_participants = model$n_participants,
      partial_R2 = model$partial_R2, F = model$F, p_value = model$p_value,
      dispersion_p = model$dispersion_p, reason = model$reason)
  }, error = function(error) data.frame(minimum_read_depth = cutoff,
    n_matrix_samples = sum(keep), n_samples = 0L, n_participants = NA_integer_,
    partial_R2 = NA_real_, F = NA_real_, p_value = NA_real_, dispersion_p = NA_real_,
    reason = conditionMessage(error)))
  result
}

load_assignments <- function(path) {
  if (!nzchar(path)) return(NULL)
  table <- read_table(path)
  require_columns(table, c("sample_id", "community_type"))
  table$key <- normalize_id(table$sample_id)
  table$community_type <- trimws(as.character(table$community_type))
  if (anyNA(table$key) || any(!nzchar(table$key)) || anyDuplicated(table$key) ||
      anyNA(table$community_type) || any(!nzchar(table$community_type)))
    stop("Community assignments require unique normalized keys and nonempty labels")
  table[, c("key", "community_type")]
}

run_analysis <- function() {
  packages <- c("dplyr", "vegan", "permute", "readxl", "ggplot2", "boot")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))
  root <- find_root()
  # This script has no figure_style.R/here() usage, so the lab-toolkit helper is sourced here,
  # right after `root` becomes available, as the closest equivalent anchor point.
  source(file.path(root, "scripts", "_lab_toolkit", "plotting_report.R"))
  out <- file.path(root, "output_reviewed", "04_taxonomic_profiling", "community_changes_v3")
  figures <- file.path(root, "figures", "Figure_3")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  save_csv <- function(table, filename) write.csv(table, file.path(out, filename), row.names = FALSE)
  inputs <- c(metadata = file.path(root, "data", "kefir4all_sample_metadata_v2.csv"),
    private = file.path(root, "data", "private", "Citizen Scientist metadata_v8.csv"),
    MetaCache = file.path(root, "data", "04_metacache_total_full_profile.csv"),
    MetaCache_audit = file.path(root, "output", "04_taxonomic_profiling", "metacache_full_profile_build_audit.csv"),
    SUPER_FOCUS = find_superfocus_file(root),
    volatiles_ML = file.path(root, "data", "Final results table all milk kefirs 0913-001.xlsx"),
    volatiles_WL = file.path(root, "data", "Final results table all water kefirs 0913-002.xlsx"))
  assignment_path <- Sys.getenv("KEFIR4ALL_COMMUNITY_TYPES")
  assignments <- load_assignments(assignment_path)
  if (nzchar(assignment_path)) inputs <- c(inputs, community_types = assignment_path)
  metadata <- load_metadata(inputs["metadata"])
  substrates <- load_substrates(inputs["private"])
  taxonomic_audit <- read_table(inputs["MetaCache_audit"])
  require_columns(taxonomic_audit, c("sample_id", "read_depth", "classified_fraction"))
  # Save input and runtime provenance before panel processing can fail.
  save_csv(metadata$excluded, "metadata_exclusions.csv")
  save_csv(substrates$unmapped, "unmapped_substrates.csv")
  save_csv(data.frame(input = names(inputs), path = unname(inputs), exists = file.exists(inputs),
    md5 = unname(tools::md5sum(inputs))), "input_manifest.csv")
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  save_csv(data.frame(package = packages, version = vapply(packages,
    function(package) as.character(utils::packageVersion(package)), character(1))), "package_versions.csv")
  sources <- list(MetaCache = tryCatch(load_taxonomy(inputs["MetaCache"], metadata$md), error = identity),
    SUPER_FOCUS = tryCatch(load_superfocus(inputs["SUPER_FOCUS"], metadata$md), error = identity),
    volatiles_ML = tryCatch(load_volatiles(inputs["volatiles_ML"], metadata$md, "ML"), error = identity),
    volatiles_WL = tryCatch(load_volatiles(inputs["volatiles_WL"], metadata$md, "WL"), error = identity))
  if (!inherits(sources$MetaCache, "error")) {
    taxonomy_metadata <- sources$MetaCache$md
    audit_index <- match(taxonomy_metadata$key, taxonomic_audit$sample_id)
    taxonomy_metadata$classified_fraction <- taxonomic_audit$classified_fraction[audit_index]
    sources$MetaCache$md <- taxonomy_metadata
    matrix_metadata <- taxonomy_metadata[
      match(rownames(sources$MetaCache$mat), taxonomy_metadata$key), , drop = FALSE]
    classified_fraction_summary <- matrix_metadata |>
      dplyr::filter(.data$Stage %in% c("T0", STAGES)) |>
      dplyr::group_by(.data$type, .data$Stage) |>
      dplyr::summarise(n = sum(is.finite(.data$classified_fraction)),
        mean = mean(.data$classified_fraction, na.rm = TRUE),
        median = stats::median(.data$classified_fraction, na.rm = TRUE), .groups = "drop")
  } else {
    classified_fraction_summary <- data.frame(type = character(), Stage = character(),
      n = integer(), mean = numeric(), median = numeric())
  }
  specs <- rbind(expand.grid(type = TYPES, source = c("MetaCache", "SUPER_FOCUS"), stringsAsFactors = FALSE),
                 data.frame(type = c("ML", "WL"), source = "volatiles"))
  # Results are accumulated by family, then written together for consistent adjustment and provenance.
  stages <- list()
  interval_summaries <- list()
  summaries <- list()
  magnitudes <- list()
  trajectories <- list()
  compositions <- list()
  rank_summaries <- list()
  rank_manifests <- list()
  depth_sensitivity_rows <- list()
  manifest <- list()
  for (index in seq_len(nrow(specs))) {
    spec <- specs[index, ]
    key <- paste(spec$source, spec$type, sep = "_")
    set.seed(SEED + index)
    messages <- character()
    status <- tryCatch(withCallingHandlers({
      loaded <- sources[[if (spec$source == "volatiles") key else spec$source]]
      if (inherits(loaded, "error")) stop(conditionMessage(loaded))
      if (!is.null(loaded$unmatched_columns))
        save_csv(loaded$unmatched_columns, paste0(key, "_unmatched_workbook_columns.csv"))
      relative <- spec$source != "volatiles"
      panel <- prepare_panel(loaded, spec$type, relative)
      save_csv(panel$excluded, paste0(key, "_exclusions.csv"))
      save_csv(panel$md, paste0(key, "_sample_metadata.csv"))
      visits <- panel$md$Stage %in% STAGES
      starter <- metadata$md[metadata$md$type == spec$type & metadata$md$Stage == "T0", , drop = FALSE]
      references <- reference_distances(panel, if (nrow(starter)) starter$key[1] else NA_character_,
        if (nrow(starter)) starter$participant[1] else NA_character_,
        if (spec$type %in% c("ML", "WL")) substrates$codes else NULL)
      intervals <- build_intervals(panel)
      save_csv(references, paste0(key, "_reference_distances.csv"))
      save_csv(intervals, paste0(key, "_interval_distances.csv"))
      for (interval_index in seq_len(length(STAGES) - 1L)) {
        from <- STAGES[interval_index]
        to <- STAGES[interval_index + 1L]
        values <- intervals$distance[intervals$from == from & intervals$to == to]
        estimate <- mean_interval(values)
        interval_summaries[[length(interval_summaries) + 1L]] <- data.frame(
          panel = key, from = from, to = to,
          assumed_week_from = WEEKS[interval_index], assumed_week_to = WEEKS[interval_index + 1L],
          n_participants = estimate["n"], mean = estimate["mean"],
          median = if (length(values)) stats::median(values) else NA_real_,
          lower = estimate["lower"], upper = estimate["upper"])
      }
      participants <- participant_summary(panel, intervals)
      save_csv(participants, paste0(key, "_participant_change.csv"))
      estimate <- mean_interval(participants$mean_interval_distance)
      magnitudes[[length(magnitudes) + 1L]] <- data.frame(panel = key,
        n_participants = nrow(participants), n_with_intervals = estimate["n"],
        n_complete = sum(participants$complete_trajectory),
        participant_weighted_mean_interval = estimate["mean"], lower = estimate["lower"], upper = estimate["upper"])
      for (stage in STAGES) {
        values <- references$distance[references$Stage == stage]
        estimate <- mean_interval(values)
        summaries[[length(summaries) + 1L]] <- data.frame(panel = key, Stage = stage,
          week = WEEKS[match(stage, STAGES)], n = estimate["n"], mean = estimate["mean"],
          median = if (length(values)) stats::median(values) else NA_real_,
          lower = estimate["lower"], upper = estimate["upper"])
      }
      reference_plot <- data.frame(panel = key, participant = references$participant, measure = "reference",
        week = WEEKS[match(references$Stage, STAGES)], distance = references$distance)
      interval_plot <- data.frame(panel = rep(key, nrow(intervals)), participant = intervals$participant,
        measure = rep("interval", nrow(intervals)), week = intervals$week, distance = intervals$distance)
      trajectories[[length(trajectories) + 1L]] <- rbind(reference_plot, interval_plot)
      model <- tryCatch(stage_test(panel), error = function(error) data.frame(n = sum(visits),
        n_participants = NA_integer_, n_stages = NA_integer_, R2_total = NA_real_, partial_R2 = NA_real_,
        F = NA_real_, p_value = NA_real_, dispersion_p = NA_real_, n_permutations = 0L,
        reason = conditionMessage(error), dispersion_reason = "stage_model_failed"))
      if (nzchar(model$reason)) messages <- c(messages, paste("Stage model:", model$reason))
      if (nzchar(model$dispersion_reason)) messages <- c(messages, paste("Dispersion:", model$dispersion_reason))
      stages[[length(stages) + 1L]] <- cbind(panel = key, model)
      if (spec$source == "MetaCache") {
        rank_result <- summarize_rank_sensitivity(loaded, spec$type, key,
          if (nrow(starter)) starter$key[1] else NA_character_,
          if (nrow(starter)) starter$participant[1] else NA_character_,
          if (spec$type %in% c("ML", "WL")) substrates$codes else NULL)
        rank_summaries[[length(rank_summaries) + 1L]] <- rank_result$summary
        rank_manifests[[length(rank_manifests) + 1L]] <- rank_result$manifest
        raw_depth <- taxonomic_audit$read_depth[
          match(rownames(loaded$mat), taxonomic_audit$sample_id)]
        depth_sensitivity_rows[[length(depth_sensitivity_rows) + 1L]] <- lapply(
          DEPTH_SENSITIVITY_CUTOFFS,
          function(cutoff) cbind(panel = key,
            test_depth_sensitivity(loaded, spec$type, raw_depth, cutoff))) |>
          dplyr::bind_rows()

        # These plots show abundant rank-qualified features, not only species or community types.
        visit_mat <- panel$mat[visits, , drop = FALSE]
        top <- head(names(sort(colMeans(visit_mat), decreasing = TRUE)), 8L)
        reduced <- cbind(visit_mat[, top, drop = FALSE], Other = rowSums(visit_mat[, setdiff(colnames(visit_mat), top), drop = FALSE]))
        composition <- data.frame(panel = key, key = rep(rownames(reduced), times = ncol(reduced)),
          feature = rep(colnames(reduced), each = nrow(reduced)), proportion = as.vector(reduced))
        position <- match(composition$key, panel$md$key)
        composition$participant <- panel$md$participant[position]
        composition$Stage <- panel$md$Stage[position]
        compositions[[length(compositions) + 1L]] <- composition
        transitions <- intervals
        transitions$from_type <- if (is.null(assignments)) rep(NA_character_, nrow(intervals))
          else assignments$community_type[match(transitions$from_key, assignments$key)]
        transitions$to_type <- if (is.null(assignments)) rep(NA_character_, nrow(intervals))
          else assignments$community_type[match(transitions$to_key, assignments$key)]
        transitions$changed <- ifelse(is.na(transitions$from_type) | is.na(transitions$to_type),
          NA, transitions$from_type != transitions$to_type)
        save_csv(transitions, paste0(key, "_community_type_transitions.csv"))
      }
      "OK"
    }, warning = function(warning) {
      messages <<- c(messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }), error = function(error) {
      messages <<- c(messages, conditionMessage(error))
      "FAILED"
    })
    manifest[[index]] <- data.frame(panel = key, status = status, messages = paste(unique(messages), collapse = " | "),
      community_types = if (spec$source != "MetaCache") "not_applicable" else if (is.null(assignments)) "not_supplied" else "supplied")
  }
  manifest_table <- do.call(rbind, manifest)
  # Preserve missing panel tests in the planned 10-test BH family.
  if (nrow(manifest_table) != nrow(specs)) stop("Manifest does not cover every planned panel")
  save_csv(manifest_table, "community_changes_v3_manifest.csv")
  if (length(stages)) {
    stage_results <- do.call(rbind, stages)
    stage_results$p_BH <- stats::p.adjust(stage_results$p_value, method = "BH", n = nrow(specs))
    save_csv(stage_results, "community_changes_v3_stage_tests.csv")
  }
  if (length(interval_summaries))
    save_csv(do.call(rbind, interval_summaries), "community_changes_v3_interval_summary.csv")
  if (length(summaries)) save_csv(do.call(rbind, summaries), "community_changes_v3_stage_summary.csv")
  if (length(magnitudes)) save_csv(do.call(rbind, magnitudes), "community_changes_v3_change_magnitude.csv")
  save_csv(classified_fraction_summary, "community_changes_v3_classified_fraction_by_stage.csv")
  if (length(rank_summaries))
    save_csv(dplyr::bind_rows(rank_summaries), "community_changes_v3_rank_sensitivity.csv")
  if (length(rank_manifests))
    save_csv(dplyr::bind_rows(rank_manifests), "community_changes_v3_rank_sensitivity_manifest.csv")
  if (length(depth_sensitivity_rows)) {
    depth_results <- dplyr::bind_rows(depth_sensitivity_rows)
    depth_results$p_BH <- NA_real_
    for (cutoff in unique(depth_results$minimum_read_depth)) {
      selected <- depth_results$minimum_read_depth == cutoff
      depth_results$p_BH[selected] <- stats::p.adjust(depth_results$p_value[selected],
        method = "BH", n = length(TYPES))
    }
    save_csv(depth_results, "community_changes_v3_depth_sensitivity.csv")
  }
  if (length(trajectories)) {
    plot_data <- do.call(rbind, trajectories)
    save_csv(plot_data, "community_changes_v3_trajectory_data.csv")
    # Figure 3: one panel per source x type, with distance-to-T0 ("reference") and
    # distance-to-previous-visit ("interval") overlaid as two mean +/- SE series,
    # individual participant values shown as faint points underneath. Panel order/
    # labels and the blue/red colour convention match the manuscript caption exactly.
    panel_label_map <- c(
      MetaCache_MG = "Milk grain · MetaCache", MetaCache_WG = "Water grain · MetaCache",
      MetaCache_ML = "Milk liquid · MetaCache", MetaCache_WL = "Water liquid · MetaCache",
      SUPER_FOCUS_MG = "Milk grain · SUPER-FOCUS", SUPER_FOCUS_WG = "Water grain · SUPER-FOCUS",
      SUPER_FOCUS_ML = "Milk liquid · SUPER-FOCUS", SUPER_FOCUS_WL = "Water liquid · SUPER-FOCUS",
      volatiles_ML = "Milk liquid · Volatiles", volatiles_WL = "Water liquid · Volatiles")
    plot_data$panel <- factor(plot_data$panel, levels = names(panel_label_map)[names(panel_label_map) %in% unique(plot_data$panel)])

    measure_labels <- c(reference = "starter (T0)", interval = "previous visit")
    plot_data$measure_label <- factor(measure_labels[plot_data$measure], levels = measure_labels)

    summary_data <- plot_data |>
      dplyr::group_by(panel, measure_label, week) |>
      dplyr::summarise(n = dplyr::n(), mean = mean(distance),
                        se = if (n > 1) stats::sd(distance) / sqrt(n) else NA_real_,
                        .groups = "drop")

    plot <- ggplot2::ggplot(summary_data, ggplot2::aes(week, mean, colour = measure_label)) +
      ggplot2::geom_point(data = plot_data, ggplot2::aes(week, distance, colour = measure_label),
                           alpha = 0.15, size = 0.6, inherit.aes = FALSE) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(linetype = measure_label)) +
      ggplot2::geom_point(size = 1.4) +
      ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = 4,
                          labeller = ggplot2::as_labeller(panel_label_map)) +
      ggplot2::scale_x_continuous(breaks = WEEKS, labels = sprintf("wk%02d", WEEKS)) +
      ggplot2::scale_colour_manual(values = c("starter (T0)" = "blue", "previous visit" = "red")) +
      ggplot2::scale_linetype_manual(values = c("starter (T0)" = "solid", "previous visit" = "dashed")) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(x = "Time since start (week)", y = "Distance", colour = "Distance to", linetype = "Distance to",
                    caption = "Points: individual participants; lines: mean ± SE. Descriptive only -- no significance test is computed in this script.")
    ggplot2::ggsave(file.path(figures, "community_changes_v3_trajectories.pdf"), plot, width = 13, height = 9)
    ggplot2::ggsave(file.path(figures, "community_changes_v3_trajectories.tiff"), plot, width = 13, height = 9, dpi = 300, compression = "lzw")
    # save_figure() provenance wrapper: additional output alongside the files above.
    # "_toolkit" suffix: the existing ggsave() calls above already write
    # "community_changes_v3_trajectories.pdf"; save_figure() would collide with that exact
    # filename (its overwrite guard correctly refused it) if given the same base name.
    # source = NULL: this panel combines MetaCache, SUPER-FOCUS and two volatiles workbooks
    # (see `inputs`), so no single file is "the" source without misrepresenting the others.
    # seed = SEED (the base constant): each panel's own randomness used SEED + index inside the
    # per-panel loop above; the plot itself (pure aggregation of already-computed `trajectories`)
    # involves no further randomness, so SEED is the most accurate single value to record.
    save_figure(plot, file.path(figures, "community_changes_v3_trajectories_toolkit"),
                width_mm = 13 * 25.4, height_mm = 9 * 25.4, source = NULL, seed = SEED)
  }
  if (length(compositions)) {
    # Only the per-panel composition data table is kept as provenance (CSV). The per-panel
    # rank-qualified-feature bar charts this data used to drive (one PDF/TIFF per MetaCache
    # panel -- MG/WG/ML/WL) were never a manuscript figure, just a diagnostic breakdown; this
    # script's only manuscript output is the Figure 3 trajectories plot saved above.
    composition <- do.call(rbind, compositions)
    save_csv(composition, "community_changes_v3_composition.csv")
  }
  writeLines(c(
    "V3 separates magnitude of observed change from common stage-associated change. No zero-distance null or stability claim is tested.",
    "participant_change files describe each participant's observed consecutive-interval distances and net first-to-last distance. Net distances span the recorded endpoints even with intervening missing visits; sums include ONLY observed consecutive pairs and are not full path lengths for incomplete trajectories. change_magnitude averages participant mean interval distances, giving each participant with at least one interval equal weight; participants without intervals remain in the counts but not that mean.",
    "T0 is one laboratory sample, not a participant baseline or pooled reference. Grains use the metadata-first default key; liquids use the first usable substrate-matched starter in metadata/workbook order. Every choice and fallback is recorded. T0 uncertainty is not quantified by these analyses.",
    "MetaCache uses the full rank-qualified rank::name profile; unclassified reads and reports below 10,000 total reads were excluded during profile construction. SUPER-FOCUS and MetaCache use row proportions and Bray-Curtis for descriptive distances; volatiles use raw concentrations and Euclidean distances. Blanks are assumed nondetections, encoded as zero. Distances from different sources are not interchangeable.",
    "community_changes_v3_classified_fraction_by_stage.csv summarizes the fraction of total raw-report reads classified by MetaCache at each type and stage; it is a coverage diagnostic, not a covariate adjustment.",
    "community_changes_v3_rank_sensitivity.csv recomputes descriptive reference and adjacent-interval distances separately within each available taxonomic rank, renormalizing within that rank. It reports rank-level sample counts and is not a set of significance tests. The manifest lists ranks with no usable T0 or visit profiles.",
    paste0("community_changes_v3_depth_sensitivity.csv reruns the participant-adjusted MetaCache Stage model after 30,000-read and 1,000,000-read filters. The primary profile remains the 10,000-read build; the stricter cutoffs are sensitivity scenarios, not recommended QC thresholds. P-values are BH-adjusted across the four MetaCache types separately at each cutoff, retaining unavailable tests in the family."),
    "Mate conflicts are checked across ALL rows sharing a study key, including nonstudy mate labels. Zero profiles and all ambiguous duplicate participant visits are excluded consistently from references and longitudinal analyses.",
    "Stage analysis is secondary: vegan::adonis2 distance ~ participant + Stage, marginal Stage test, permutations restricted within participant. Singleton participants are omitted from this model, not descriptive outputs. Bray-Curtis is square-root transformed ONLY for the multivariate tests to reduce non-Euclidean geometry; raw descriptive distances are unchanged. R2_total uses total inertia; partial_R2 uses Stage/(Stage+residual) sums of squares.",
    "Stage inference requires exchangeable residuals within participants. Blocking alone does not establish this: serial correlation, heterogeneous variances, and informative missingness can invalidate inference. This is an exploratory conditional analysis, not proof of a universal temporal trajectory or a causal effect. No time-series permutation scheme is imposed on irregular, incomplete trajectories.",
    "betadisper uses stage spatial medians, small-sample bias adjustment, the same model subset and permutation matrix. Its p-value is an unadjusted exploratory diagnostic and does not prove that a PERMANOVA result is solely location or solely dispersion. Participant adjustment in PERMANOVA is not reproduced by betadisper.",
    "community_changes_v3_interval_summary.csv reports each adjacent interval separately, using participant-level mean, median, n, and participant-bootstrap 95% interval. It does not pool intervals into early and late windows or test for a plateau. The bootstrap resamples independent participants; estimates are available-case, pointwise, and do not correct informative attrition or T0-reference uncertainty.",
    "Stage p-values are BH-adjusted across the ten planned panels, retaining failed or unavailable panel tests in the family size. The stage analysis is exploratory and is not a test of equivalence or stability.",
    "Community types are imported only if KEFIR4ALL_COMMUNITY_TYPES names a CSV with unique normalized sample_id and community_type. Only adjacent observed taxonomic visits are linked; unmatched labels remain NA. No labels are derived here. Label switching is descriptive and depends on assignment robustness; it does not itself prove distinct ecological states. Transition files are always written with NA labels when unavailable to prevent stale results from a previous labelled run.",
    "Outputs use a separate v3 directory and figure prefix. Manifest records panel/model failures and warnings; partial panel files can exist for failed panels and must not be interpreted as completed analyses. Input hashes, package versions, per-panel deterministic seeds and sessionInfo record provenance. Private input hashes are recorded, not distributed private data.",
    paste("Seed:", SEED, "plus panel index; permutations:", NPERM, "bootstrap replicates:", NBOOT)),
    file.path(out, "METHODS.txt"))
  print(manifest_table, row.names = FALSE)
  if (any(manifest_table$status != "OK"))
    stop("Some panels failed; see community_changes_v3_manifest.csv")
  cat("V3 outputs: ", out, "\n", sep = "")
}

run_self_tests <- function() {
  stopifnot(identical(participant_id(c("5", "ID05", "FU")), c("ID5", "ID5", "FU")))
  metadata <- data.frame(key = c("starter", "first", "third"),
    participant = c("FU", "ID1", "ID1"), type = "MG", Stage = c("T0", "T1", "T3"))
  mat <- rbind(c(0.5, 0.5), c(1, 0), c(0, 1))
  rownames(mat) <- metadata$key
  panel <- prepare_panel(list(mat = mat[c(3, 2, 1), ], md = metadata), "MG", TRUE)
  stopifnot(identical(panel$md$key, metadata$key), nrow(build_intervals(panel)) == 0L)
  references <- reference_distances(panel, "starter", "FU", NULL)
  stopifnot(all(references$distance == 0.5), all(references$reference_key == "starter"))
  metadata$Stage[3] <- "T2"
  panel <- prepare_panel(list(mat = mat, md = metadata), "MG", TRUE)
  stopifnot(build_intervals(panel)$distance == 1)
  metadata <- rbind(metadata, metadata[2, ])
  metadata$key[4] <- "duplicate"
  mat <- rbind(mat, duplicate = c(1, 0))
  panel <- prepare_panel(list(mat = mat, md = metadata), "MG", TRUE)
  stopifnot(nrow(panel$excluded) == 2L, !any(panel$md$Stage == "T1"))
  intervals <- build_intervals(panel)
  stopifnot(all(intervals$to == paste0("T", as.integer(sub("T", "", intervals$from)) + 1L)))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  write.csv(data.frame(merge_column = c("key_host_removed_R1.fastq.gz", "key_host_removed_R2.fastq.gz"),
    Sample = c("ID1", "control"), `kefir type` = c("MG", "Medium control"),
    Stage = c("T1", "T0"), check.names = FALSE), path, row.names = FALSE)
  result <- load_metadata(path)
  stopifnot(nrow(result$md) == 0L, nrow(result$excluded) == 2L)
  set.seed(SEED)
  metadata <- expand.grid(Stage = STAGES, participant = paste0("ID", 1:5), stringsAsFactors = FALSE)
  metadata$key <- paste(metadata$participant, metadata$Stage, sep = "_")
  metadata$type <- "MG"
  mat <- matrix(stats::runif(nrow(metadata) * 4), nrow(metadata), 4)
  rownames(mat) <- metadata$key
  panel <- prepare_panel(list(mat = mat, md = metadata), "MG", TRUE)
  result <- stage_test(panel)
  stopifnot(is.finite(result$p_value), result$n == 30L, result$n_participants == 5L)
  selected <- metadata$Stage != "T6" | metadata$participant == "ID1"
  panel <- prepare_panel(list(mat = mat[selected, ], md = metadata[selected, ]), "MG", TRUE)
  result <- stage_test(panel)
  stopifnot(is.finite(result$p_value), is.na(result$dispersion_p),
    result$dispersion_reason == "stage_has_fewer_than_two_samples")
  cat("V3 focused self-tests PASS\n")
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if ("--self-test" %in% arguments) {
    run_self_tests()
  } else if ("--verify-reproducibility" %in% arguments) {
    run_analysis()
    folder <- file.path(find_root(), "output", "04_taxonomic_profiling", "community_changes_v3")
    files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
    before <- tools::md5sum(files)
    tables <- lapply(files, read_table)
    run_analysis()
    after <- tools::md5sum(files)
    changed <- which(before != after)
    for (index in changed) {
      cat("Reproducibility mismatch:", basename(files[index]), "\n")
      print(all.equal(tables[[index]], read_table(files[index]), tolerance = 0))
    }
    stopifnot(identical(before, after))
    cat("V3 Docker CSV reproducibility PASS\n")
  } else {
    run_analysis()
  }
}