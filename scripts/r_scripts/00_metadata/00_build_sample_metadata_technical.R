# Build an enriched, technical sample-metadata file for reviewer/internal use.
# Starts from data/kefir4all_sample_metadata_v2.csv (not rebuilt from
# filenames -- its existing corrections are kept, not redone) and adds ENA
# sequencing depth (from the cached PRJEB77409 run report) and MetaCache
# full-profile inclusion status (from metacache_full_profile_build_audit.csv,
# when available).
#
# Standalone by design (no shared helper file to source), so find_root() and
# normalize_id() are duplicated from build_metacache_full_profile.R rather
# than factored out -- a deliberate, small trade-off, not an oversight.
#
# Outputs:
#   data/sample_metadata_technical.csv
#   output/04_taxonomic_profiling/sample_counts_by_type_timepoint.csv
#
# Does not overwrite kefir4all_sample_metadata_v2.csv, does not modify any
# other script, and is not wired into any analysis pipeline. Sourcing this
# file loads the functions only; running it as a script also executes main().

find_root <- function() {
  starts <- c(Sys.getenv("KEFIR4ALL_ROOT"), getwd(), file.path(getwd(), "kefir4all_review_history"))
  for (start in starts[nzchar(starts)]) {
    if (file.exists(file.path(start, "data", "kefir4all_sample_metadata_v2.csv")))
      return(normalizePath(start, winslash = "/"))
  }
  stop("Run from the repository root or set KEFIR4ALL_ROOT.")
}

normalize_id <- function(x) {
  x <- basename(gsub("\\\\", "/", trimws(as.character(x))))
  x <- sub("_host_removed_R[0-9]\\.fastq\\.gz$", "", x)
  gsub("-", "_", x, fixed = TRUE)
}

# A participant label matching "ID###" is a citizen-scientist visit
# (biological_sample). Anything else with a real kefir type (MG/WG/ML/WL) at
# stage T0 is a lab-prepared starter reference for that type (T0_reference) --
# distinct from a generic lab control, and must keep its type visible rather
# than being lumped in with substrate-code controls. Everything remaining
# (2-letter lab substrate codes, "Medium control", etc.) is a control.
classify_sample_role <- function(participant, type_cleaned, stage_cleaned) {
  known_types <- c("MG", "WG", "ML", "WL")
  is_participant <- grepl("^ID[0-9]+$", participant)
  is_t0_reference <- !is_participant & type_cleaned %in% known_types & stage_cleaned == "T0"
  ifelse(is_participant, "biological_sample", ifelse(is_t0_reference, "T0_reference", "control"))
}

# Loads the base metadata as-is, normalises just enough to join and compare
# (participant label, key), and -- critically -- checks conflicts only among
# "candidate" rows (a real type, a real stage, a non-empty participant),
# matching the scoping community_changes.R/v2.R already use, rather than
# comparing every raw row regardless of whether it looks like a study sample.
# Conflicting keys are flagged AND have the actual disagreeing labels
# recorded, not just a boolean -- the point of this file is to let a reader
# see the disagreement without going back to the raw CSV themselves.
load_and_dedupe_metadata <- function(base_path) {
  base <- read.csv(base_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("merge_column", "Sample", "kefir type", "Stage")
  missing <- setdiff(required, names(base))
  if (length(missing)) stop("Base metadata missing columns: ", paste(missing, collapse = ", "))

  base$sample_id_normalized <- normalize_id(base$merge_column)
  base$participant_cleaned <- toupper(trimws(as.character(base$Sample)))
  numeric_id <- grepl("^(ID)?0*[0-9]+$", base$participant_cleaned)
  base$participant_cleaned[numeric_id] <- paste0("ID", as.integer(sub("^ID", "", base$participant_cleaned[numeric_id])))
  base$type_cleaned <- trimws(base[["kefir type"]])
  base$stage_cleaned <- toupper(trimws(base$Stage))
  base$sample_role <- classify_sample_role(base$participant_cleaned, base$type_cleaned, base$stage_cleaned)

  known_types <- c("MG", "WG", "ML", "WL")
  known_stages <- c("T0", paste0("T", 1:6))
  candidate <- base$type_cleaned %in% known_types & base$stage_cleaned %in% known_stages &
    !is.na(base$participant_cleaned) & nzchar(base$participant_cleaned)

  signature <- paste(base$participant_cleaned, base$type_cleaned, base$stage_cleaned, sep = "|")
  by_key <- split(signature[candidate], base$sample_id_normalized[candidate])
  conflict_labels <- vapply(by_key, function(z) paste(unique(z), collapse = " vs "), character(1))
  is_conflicting <- vapply(by_key, function(z) length(unique(z)) > 1L, logical(1))

  base$conflicting_metadata <- is_conflicting[base$sample_id_normalized]
  base$conflicting_metadata[is.na(base$conflicting_metadata)] <- FALSE
  base$conflicting_labels <- ifelse(base$conflicting_metadata,
                                     conflict_labels[base$sample_id_normalized], NA_character_)

  dup_counts <- table(base$sample_id_normalized)
  diagnostics <- list(
    n_raw_rows = nrow(base),
    n_unique_keys = length(dup_counts),
    n_duplicate_keys = sum(dup_counts > 1),
    n_conflicting_keys = sum(is_conflicting))

  enriched <- base[!duplicated(base$sample_id_normalized), , drop = FALSE]
  list(enriched = enriched, diagnostics = diagnostics)
}

# Fetches and caches the ENA run report when it's missing locally (e.g. a
# fresh clone, since data/ena_PRJEB77409_run_report.tsv is not committed).
# Never errors out -- a failed fetch (no network, API change, etc.) just
# leaves the file absent, and the caller falls back to NA ena_* fields.
fetch_ena_report <- function(ena_path) {
  url <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/filereport",
    "?accession=PRJEB77409&result=read_run",
    "&fields=run_accession,sample_accession,sample_alias,library_name,library_layout,read_count,base_count",
    "&format=tsv")
  ok <- tryCatch({
    utils::download.file(url, destfile = ena_path, quiet = TRUE, mode = "wb")
    file.exists(ena_path) && file.info(ena_path)$size > 0
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!isTRUE(ok)) { unlink(ena_path); FALSE } else TRUE
}

# NULL if the ENA report isn't present -- callers must handle that, not
# invent fallback values. Refuses to aggregate over a duplicated
# run_accession (would silently inflate summed counts), and keeps a sample's
# deposited reads/bases as NA -- not 0 -- if any of its runs are missing a
# count, rather than letting na.rm=TRUE hide the gap.
load_ena_report <- function(ena_path) {
  if (!file.exists(ena_path)) return(NULL)
  ena <- read.delim(ena_path, stringsAsFactors = FALSE)
  if (anyDuplicated(ena$run_accession))
    stop("Duplicate run_accession values in ", ena_path, " -- refusing to aggregate silently.")
  ena$sample_id_normalized <- normalize_id(ena$sample_alias)
  do.call(rbind, lapply(split(ena, ena$sample_id_normalized), function(z) data.frame(
    sample_id_normalized = z$sample_id_normalized[1],
    ena_run_accessions = paste(z$run_accession, collapse = ";"),
    ena_n_runs = nrow(z),
    ena_library_layout = paste(unique(z$library_layout), collapse = ";"),
    ena_deposited_reads = if (anyNA(z$read_count)) NA_real_ else sum(z$read_count),
    ena_deposited_bases = if (anyNA(z$base_count)) NA_real_ else sum(z$base_count),
    stringsAsFactors = FALSE)))
}

# NULL if the MetaCache build audit isn't present. Refuses to merge on a
# duplicated sample_id (would silently fan out enriched into duplicate rows).
load_metacache_audit <- function(audit_path) {
  if (!file.exists(audit_path)) return(NULL)
  audit <- read.csv(audit_path, stringsAsFactors = FALSE)
  if (anyDuplicated(audit$sample_id))
    stop("Duplicate sample_id values in ", audit_path, " -- refusing to merge silently.")
  data.frame(
    sample_id_normalized = audit$sample_id,
    metacache_included = audit$included,
    metacache_exclusion_reason = ifelse(audit$included, NA_character_, audit$reason),
    metacache_read_depth = audit$read_depth,
    metacache_classified_fraction = audit$classified_fraction,
    stringsAsFactors = FALSE)
}

# NA (unknown), not 0, when inclusion status is unavailable for every row in
# a cell -- sum(..., na.rm = TRUE) over an all-NA vector silently returns 0,
# which would read as "confirmed excluded" rather than "no information".
safe_inclusion_sum <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_integer_)
  sum(x, na.rm = TRUE)
}

build_summary <- function(enriched) {
  bio <- enriched[enriched$sample_role == "biological_sample", , drop = FALSE]
  t0_ref <- enriched[enriched$sample_role == "T0_reference", , drop = FALSE]
  ctrl <- enriched[enriched$sample_role == "control", , drop = FALSE]
  stages <- c("T0", paste0("T", 1:6))
  types <- c("MG", "WG", "ML", "WL")

  row_for <- function(category, type, sub) data.frame(
    category = category, type = type, n_in_metadata = nrow(sub),
    n_in_metacache_profile = safe_inclusion_sum(sub$metacache_included),
    n_conflicting = sum(sub$conflicting_metadata, na.rm = TRUE),
    stringsAsFactors = FALSE)

  bio_rows <- do.call(rbind, lapply(types, function(ty) do.call(rbind, lapply(stages, function(st)
    cbind(stage = st, row_for("biological_sample", ty, bio[bio$type_cleaned == ty & bio$stage_cleaned == st, ]))))))

  t0_ref_rows <- do.call(rbind, lapply(types, function(ty)
    cbind(stage = "T0", row_for("T0_reference", ty, t0_ref[t0_ref$type_cleaned == ty, ]))))

  ctrl_rows <- do.call(rbind, lapply(stages, function(st)
    cbind(stage = st, row_for("control", NA_character_, ctrl[ctrl$stage_cleaned == st, ]))))

  rbind(bio_rows, t0_ref_rows, ctrl_rows)
}

# Orchestrates the whole build; returns the pieces rather than writing them.
build_technical_metadata <- function(root) {
  data_dir <- file.path(root, "data")
  out_dir <- file.path(root, "output", "04_taxonomic_profiling")

  base <- load_and_dedupe_metadata(file.path(data_dir, "kefir4all_sample_metadata_v2.csv"))
  enriched <- base$enriched
  diagnostics <- base$diagnostics

  ena_path <- file.path(data_dir, "ena_PRJEB77409_run_report.tsv")
  if (!file.exists(ena_path)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    diagnostics$ena_fetched <- fetch_ena_report(ena_path)
  }
  ena_agg <- load_ena_report(ena_path)
  if (is.null(ena_agg)) {
    for (col in c("ena_run_accessions", "ena_n_runs", "ena_library_layout",
                  "ena_deposited_reads", "ena_deposited_bases")) enriched[[col]] <- NA
    diagnostics$ena_status <- paste("MISSING (not cached, and fetch did not succeed):", ena_path)
  } else {
    enriched <- merge(enriched, ena_agg, by = "sample_id_normalized", all.x = TRUE)
    diagnostics$ena_matched <- sum(!is.na(enriched$ena_n_runs))
    diagnostics$ena_unmatched <- sum(is.na(enriched$ena_n_runs))
    diagnostics$ena_multi_run <- sum(enriched$ena_n_runs > 1, na.rm = TRUE)
  }

  audit_path <- file.path(out_dir, "metacache_full_profile_build_audit.csv")
  audit <- load_metacache_audit(audit_path)
  if (is.null(audit)) {
    enriched$metacache_included <- NA
    enriched$metacache_exclusion_reason <- NA_character_
    enriched$metacache_read_depth <- NA
    enriched$metacache_classified_fraction <- NA
    diagnostics$metacache_status <- paste("MISSING:", audit_path)
  } else {
    enriched <- merge(enriched, audit, by = "sample_id_normalized", all.x = TRUE)
    diagnostics$metacache_matched <- sum(!is.na(enriched$metacache_included))
    diagnostics$metacache_unmatched <- sum(is.na(enriched$metacache_included))
  }

  if (nrow(enriched) != diagnostics$n_unique_keys)
    stop("Row count after joins (", nrow(enriched), ") does not match unique sample keys (",
         diagnostics$n_unique_keys, ") -- a join key was not unique.")

  list(enriched = enriched, summary = build_summary(enriched),
       diagnostics = diagnostics, out_dir = out_dir, data_dir = data_dir)
}

main <- function() {
  root <- find_root()
  result <- build_technical_metadata(root)
  d <- result$diagnostics

  cat("=== Base metadata ===\n")
  cat("Raw rows (one per mate read):", d$n_raw_rows, "\n")
  cat("Unique normalised sample keys:", d$n_unique_keys, "\n")
  cat("Keys with >1 raw row (mate pairs or more):", d$n_duplicate_keys, "\n")
  cat("Keys with conflicting participant/type/stage across their raw rows:", d$n_conflicting_keys, "\n\n")

  if (!is.null(d$ena_fetched))
    message(if (d$ena_fetched) "Fetched and cached the ENA (PRJEB77409) run report." else "ENA fetch failed (no network?).")

  if (!is.null(d$ena_status)) {
    message(d$ena_status, " -- ena_* fields left NA. Re-fetch via the ENA Portal API if needed.")
  } else {
    cat("=== ENA (PRJEB77409) join ===\n")
    cat("Samples matched to >=1 ENA run:", d$ena_matched, "\n")
    cat("Samples with NO ENA match:", d$ena_unmatched, "\n")
    cat("Samples with >1 ENA run (counts summed, accessions concatenated):", d$ena_multi_run, "\n\n")
  }

  if (!is.null(d$metacache_status)) {
    message(d$metacache_status, " -- metacache_* fields left NA. Run build_metacache_full_profile.R first if needed.")
  } else {
    cat("=== MetaCache audit join ===\n")
    cat("Samples matched to an audit record:", d$metacache_matched, "\n")
    cat("Samples with no audit record:", d$metacache_unmatched, "\n\n")
  }

  dir.create(result$data_dir, recursive = TRUE, showWarnings = FALSE)
  enriched_path <- file.path(result$data_dir, "sample_metadata_technical.csv")
  write.csv(result$enriched, enriched_path, row.names = FALSE, na = "")
  cat("Wrote:", enriched_path, "(", nrow(result$enriched), "rows )\n\n")

  dir.create(result$out_dir, recursive = TRUE, showWarnings = FALSE)
  summary_path <- file.path(result$out_dir, "sample_counts_by_type_timepoint.csv")
  write.csv(result$summary, summary_path, row.names = FALSE)
  cat("Wrote:", summary_path, "\n\n")

  bio <- result$enriched[result$enriched$sample_role == "biological_sample", ]
  t0_ref <- result$enriched[result$enriched$sample_role == "T0_reference", ]
  cat("=== Totals (for comparison against any manuscript-reported counts -- not forced to match) ===\n")
  cat("Biological samples (MG+WG+ML+WL, participant visits, stages T1-T6):", nrow(bio), "\n")
  for (ty in c("MG", "WG", "ML", "WL")) cat("  ", ty, ":", sum(bio$type_cleaned == ty), "\n")
  cat("T0 lab-prepared starter references (by type, not tied to a participant):", nrow(t0_ref), "\n")
  for (ty in c("MG", "WG", "ML", "WL")) cat("  ", ty, ":", sum(t0_ref$type_cleaned == ty), "\n")
  cat("Control rows (no real kefir type):", sum(result$enriched$sample_role == "control"), "\n")
  invisible(result)
}

if (sys.nframe() == 0L) main()
