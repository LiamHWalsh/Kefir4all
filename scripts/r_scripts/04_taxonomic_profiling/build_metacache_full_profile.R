# Build the mixed-rank MetaCache matrix used by community_changes.R, v2, and v3.
# Set KEFIR4ALL_METACACHE_RAW_DIRS to the colon-separated raw report directories.
# Values are read counts; community_changes scripts convert classified reads to proportions.

MIN_READ_DEPTH <- 10000

find_root <- function() {
  starts <- c(Sys.getenv("KEFIR4ALL_ROOT"), getwd(),
              file.path(getwd(), "kefir4all_review_history"))
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

study_metadata <- function(path) {
  metadata <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("merge_column", "Sample", "kefir type", "Stage")
  if (!all(required %in% names(metadata)))
    stop("Metadata missing columns: ", paste(setdiff(required, names(metadata)), collapse = ", "))

  metadata$key <- normalize_id(metadata$merge_column)
  metadata$participant <- toupper(trimws(as.character(metadata$Sample)))
  numeric_id <- !is.na(metadata$participant) & grepl("^(ID)?0*[0-9]+$", metadata$participant)
  metadata$participant[numeric_id] <- paste0(
    "ID", as.integer(sub("^ID", "", metadata$participant[numeric_id])))
  metadata$type <- trimws(metadata[["kefir type"]])
  metadata$Stage <- toupper(trimws(metadata$Stage))
  candidate <- metadata$type %in% c("MG", "WG", "ML", "WL") &
    metadata$Stage %in% c("T0", paste0("T", 1:6)) &
    !is.na(metadata$participant) & nzchar(metadata$participant)
  candidate_keys <- unique(metadata$key[candidate])

  # Compare all mate rows for keys that have at least one eligible study row.
  relevant <- metadata$key %in% candidate_keys
  signature <- paste(metadata$participant, metadata$type, metadata$Stage, sep = "|")
  conflicts <- vapply(split(signature[relevant], metadata$key[relevant]),
    function(values) length(unique(values)) > 1L, logical(1))
  conflict_keys <- names(conflicts)[conflicts]
  list(keys = candidate_keys[!candidate_keys %in% conflict_keys],
       conflict_keys = conflict_keys)
}

parse_report_lines <- function(lines) {
  data_lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*#", lines)]
  if (!length(data_lines)) stop("Report contains no classification rows")

  parsed <- lapply(seq_along(data_lines), function(index) {
    fields <- trimws(strsplit(data_lines[index], "\\s*\\|\\s*", perl = TRUE)[[1]])
    if (length(fields) < 4L) stop("Malformed classification row: ", data_lines[index])
    rank <- tolower(fields[1])
    reads <- suppressWarnings(as.numeric(fields[length(fields) - 1L]))
    if (is.na(reads) || !is.finite(reads) || reads < 0)
      stop("Invalid read count in classification row: ", data_lines[index])

    feature <- NA_character_
    if (tolower(rank) != "unclassified") {
      if (length(fields) < 5L || !nzchar(fields[2]))
        stop("Classified row has no taxon name: ", data_lines[index])
      feature <- paste(rank, fields[2], sep = "::")
    }
    data.frame(rank = rank, feature = feature, reads = reads, stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, parsed)
  total_reads <- sum(rows$reads)
  unclassified_reads <- sum(rows$reads[rows$rank == "unclassified"])
  assigned <- rows[!is.na(rows$feature), , drop = FALSE]
  classified_reads <- sum(assigned$reads)
  counts <- if (nrow(assigned)) tapply(assigned$reads, assigned$feature, sum) else numeric()
  counts <- counts[counts > 0]
  list(total_reads = total_reads, classified_reads = classified_reads,
       unclassified_reads = unclassified_reads,
       classified_fraction = if (total_reads > 0) classified_reads / total_reads else NA_real_,
       counts = counts)
}

parse_report <- function(path) {
  parse_report_lines(readLines(path, warn = FALSE, encoding = "UTF-8"))
}

build_full_profile <- function(root = find_root(), raw_dirs = NULL,
                                min_read_depth = MIN_READ_DEPTH) {
  if (is.null(raw_dirs)) {
    raw_dirs <- strsplit(Sys.getenv("KEFIR4ALL_METACACHE_RAW_DIRS"),
                         .Platform$path.sep, fixed = TRUE)[[1]]
  }
  raw_dirs <- raw_dirs[nzchar(raw_dirs)]
  if (!length(raw_dirs) || any(!dir.exists(raw_dirs)))
    stop("Set KEFIR4ALL_METACACHE_RAW_DIRS to existing raw MetaCache report directories.")
  if (length(min_read_depth) != 1L || is.na(min_read_depth) ||
      !is.finite(min_read_depth) || min_read_depth < 0)
    stop("min_read_depth must be one finite, nonnegative number.")

  reports <- sort(unlist(lapply(raw_dirs, list.files, pattern = "^TG_.*\\.txt$",
                                full.names = TRUE, recursive = FALSE), use.names = FALSE))
  if (!length(reports)) stop("No TG_*.txt reports found in raw directories")
  filenames <- basename(reports)
  if (anyDuplicated(filenames))
    stop("Duplicate raw report filenames across snapshots: ",
         paste(unique(filenames[duplicated(filenames)]), collapse = ", "))
  report_ids <- normalize_id(sub("\\.txt$", "", filenames))
  if (anyDuplicated(report_ids))
    stop("Different report filenames normalize to the same sample ID")
  report_hashes <- unname(tools::md5sum(reports))

  metadata <- study_metadata(file.path(root, "data", "kefir4all_sample_metadata_v2.csv"))
  profiles <- list()
  audit <- vector("list", length(reports) + length(metadata$keys))
  audit_index <- 0L
  audit_row <- function(sample_id, raw_file = NA_character_, raw_file_md5 = NA_character_,
                        parsed = NULL, included = FALSE, reason) {
    data.frame(sample_id = sample_id, raw_file = raw_file, raw_file_md5 = raw_file_md5,
      read_depth = if (is.null(parsed)) NA_real_ else parsed$total_reads,
      classified_reads = if (is.null(parsed)) NA_real_ else parsed$classified_reads,
      unclassified_reads = if (is.null(parsed)) NA_real_ else parsed$unclassified_reads,
      classified_fraction = if (is.null(parsed)) NA_real_ else parsed$classified_fraction,
      minimum_read_depth = min_read_depth, included = included, reason = reason,
      stringsAsFactors = FALSE)
  }

  for (index in seq_along(reports)) {
    sample_id <- report_ids[index]
    if (sample_id %in% metadata$conflict_keys) {
      audit_index <- audit_index + 1L
      audit[[audit_index]] <- audit_row(sample_id, filenames[index], report_hashes[index],
        reason = "conflicting_metadata")
      next
    }
    if (!sample_id %in% metadata$keys) {
      audit_index <- audit_index + 1L
      audit[[audit_index]] <- audit_row(sample_id, filenames[index], report_hashes[index],
        reason = "not_in_study_metadata")
      next
    }

    parsed <- tryCatch(parse_report(reports[index]), error = function(error)
      stop("Failed to parse ", reports[index], ": ", conditionMessage(error)))
    reason <- if (parsed$total_reads < min_read_depth) "below_minimum_read_depth" else
      if (!length(parsed$counts)) "no_classified_reads" else "included"
    audit_index <- audit_index + 1L
    audit[[audit_index]] <- audit_row(sample_id, filenames[index], report_hashes[index],
      parsed, reason == "included", reason)
    if (reason == "included") profiles[[sample_id]] <- parsed$counts
  }

  missing_reports <- setdiff(metadata$keys, report_ids)
  if (length(missing_reports)) {
    for (sample_id in missing_reports) {
      audit_index <- audit_index + 1L
      audit[[audit_index]] <- audit_row(sample_id, reason = "missing_raw_report")
    }
  }
  audit <- do.call(rbind, audit[seq_len(audit_index)])

  included_ids <- metadata$keys[metadata$keys %in% names(profiles)]
  if (!length(included_ids)) stop("No metadata-matched reports passed read-depth QC")
  features <- sort(unique(unlist(lapply(profiles[included_ids], names), use.names = FALSE)))
  abundance <- matrix(0, nrow = length(included_ids), ncol = length(features),
    dimnames = list(included_ids, features))
  for (sample_id in included_ids) {
    counts <- profiles[[sample_id]]
    abundance[sample_id, names(counts)] <- counts
  }

  data_dir <- file.path(root, "data")
  output_path <- file.path(data_dir, "04_metacache_total_full_profile.csv")
  audit_path <- file.path(root, "output", "04_taxonomic_profiling",
                          "metacache_full_profile_build_audit.csv")
  dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
  profile_tmp <- tempfile("metacache_full_profile_", tmpdir = data_dir)
  audit_tmp <- tempfile("metacache_profile_audit_", tmpdir = dirname(audit_path))
  on.exit(unlink(c(profile_tmp, audit_tmp)), add = TRUE)
  write.csv(data.frame(sample_id = included_ids, abundance, check.names = FALSE),
            profile_tmp, row.names = FALSE, na = "")
  write.csv(audit, audit_tmp, row.names = FALSE, na = "")
  if (!file.copy(profile_tmp, output_path, overwrite = TRUE) ||
      !file.copy(audit_tmp, audit_path, overwrite = TRUE))
    stop("Failed to publish full profile or build audit")

  cat("Raw reports:", length(reports), "\n")
  cat("Matched, QC-passing samples:", length(included_ids), "\n")
  cat("Rank-qualified features:", length(features), "\n")
  cat("Minimum total reads:", min_read_depth, "\n")
  cat("Profile:", output_path, "\n")
  cat("Audit:", audit_path, "\n")
  invisible(list(profile = output_path, audit = audit_path,
                 n_samples = length(included_ids), n_features = length(features)))
}

run_self_tests <- function() {
  parsed <- parse_report_lines(c(
    "species | Lactococcus lactis | 1358 | 5.12719e+07 | 58.7%",
    "species | Lactococcus lactis | 1358 | 2 | 0.01%",
    "genus | Lactococcus | 1358 | 500 | 0.5%",
    "unclassified | 0 | 700 | 40.79%"))
  stopifnot(parsed$total_reads == 51273102,
    parsed$classified_reads == 51272402,
    parsed$unclassified_reads == 700,
    isTRUE(all.equal(parsed$classified_fraction, 51272402 / 51273102)),
    parsed$counts[["species::Lactococcus lactis"]] == 51271902,
    parsed$counts[["genus::Lactococcus"]] == 500,
    !any(grepl("unclassified", names(parsed$counts))))

  root <- tempfile("metacache_profile_test_")
  raw_dir <- file.path(root, "raw")
  dir.create(file.path(root, "data"), recursive = TRUE)
  dir.create(raw_dir)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  metadata <- data.frame(
    merge_column = c("TG_ID001_MG_T1", "TG_ID002_MG_T1", "TG_ID003_MG_T1",
                     "TG_FG_WLC_S410", "TG_FG_WLC_S410"),
    Sample = c("ID1", "ID2", "ID3", "FG", "Medium control"),
    `kefir type` = c("MG", "MG", "MG", "WG", "WL"),
    Stage = c("T1", "T1", "T1", "T0", "T0"), check.names = FALSE)
  write.csv(metadata, file.path(root, "data", "kefir4all_sample_metadata_v2.csv"), row.names = FALSE)
  writeLines(c("species | Lactococcus lactis | 1358 | 12000 | 80.0%",
               "unclassified | 0 | 3000 | 20.0%"),
    file.path(raw_dir, "TG_ID001_MG_T1.txt"))
  writeLines(c("genus | Lactococcus | 1358 | 7000 | 70.0%",
               "unclassified | 0 | 3000 | 30.0%"),
    file.path(raw_dir, "TG_ID002_MG_T1.txt"))
  writeLines("unclassified | 0 | 1 | 100.0%", file.path(raw_dir, "TG_FG_WLC_S410.txt"))
  outputs <- build_full_profile(root, raw_dir, min_read_depth = 10001)
  profile <- read.csv(outputs$profile, check.names = FALSE)
  audit <- read.csv(outputs$audit, stringsAsFactors = FALSE)
  stopifnot(nrow(profile) == 1L, profile$sample_id == "TG_ID001_MG_T1",
    profile[["species::Lactococcus lactis"]] == 12000,
    audit$reason[audit$sample_id == "TG_ID002_MG_T1"] == "below_minimum_read_depth",
    audit$reason[audit$sample_id == "TG_FG_WLC_S410"] == "conflicting_metadata",
    audit$reason[audit$sample_id == "TG_ID003_MG_T1"] == "missing_raw_report",
    audit$minimum_read_depth[audit$sample_id == "TG_ID001_MG_T1"] == 10001,
    audit$classified_fraction[audit$sample_id == "TG_ID001_MG_T1"] == 0.8,
    nzchar(audit$raw_file_md5[audit$sample_id == "TG_ID001_MG_T1"]))
  cat("Full-profile parser self-tests PASS\n")
}

if (sys.nframe() == 0L) {
  if ("--self-test" %in% commandArgs(trailingOnly = TRUE)) run_self_tests()
  else build_full_profile()
}
