# slim_packages.R — strip unused R packages from the library
# Reduces the 551-package library to only its strong-dependency closure
# for the 75 packages actually required by the analysis scripts.
#
# Usage in Dockerfile AFTER library tar extraction:
#   Rscript slim_packages.R <keep_list.txt>
#
# keep_list.txt: one package name per line (the scripts' library() targets)

args <- commandArgs(trailingOnly = TRUE)
target_file <- args[1]
if (!file.exists(target_file)) {
  cat("ERROR: keep-list file not found:", target_file, "\n")
  quit(status = 1)
}

targets <- scan(target_file, what = character(), quiet = TRUE)
lib <- "/opt/R/4.4.2/lib/R/library"

# All currently installed packages (excluding base/recommended via priority filter)
ip <- installed.packages(lib.loc = lib, priority = "NA")
installed_pkgs <- rownames(ip)

# Compute recursive STRONG dependencies (Imports/LinkingTo/Depends only —
# the true runtime minimum; Suggests are excluded to allow maximal pruning)
strong_deps <- tools::package_dependencies(
  targets, db = ip, which = "strong", recursive = TRUE
)

all_needed <- unique(c(targets, unlist(strong_deps, use.names = FALSE)))
all_needed <- intersect(all_needed, installed_pkgs)

to_remove <- setdiff(installed_pkgs, all_needed)

cat(sprintf("Keeping:   %d packages (%.1f MB)\n",
            length(all_needed),
            sum(sapply(all_needed, file.size,  # quick approximation
                       simplify = TRUE) | 0) / 1048576))

cat(sprintf("Removing:  %d packages\n", length(to_remove)))

if (length(to_remove) > 0) {
  removed_size <- 0
  for (p in to_remove) {
    pkg_dir <- file.path(lib, p)
    if (dir.exists(pkg_dir)) {
      fs <- list.files(pkg_dir, recursive = TRUE, full.names = TRUE,
                       include.dirs = TRUE)
      info <- file.info(fs)
      removed_size <- removed_size +
        sum(info$size[!info$isdir], na.rm = TRUE)
      unlink(pkg_dir, recursive = TRUE)
    }
  }
  cat(sprintf("Freed:     %.1f MB\n", removed_size / 1048576))
}

cat(sprintf("Final:     %d packages in library\n",
            length(list.dirs(lib, recursive = FALSE))))
