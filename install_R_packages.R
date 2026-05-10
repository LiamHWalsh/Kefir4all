# install_R_packages.R ---------------------------------------------------
# Helper script that installs every R package the figure scripts need.
# Run once on a fresh R installation:
#
#   Rscript install_R_packages.R
#
# Tested with R >= 4.2. For exact-version reproducibility consider using
# renv (see REPRODUCIBILITY.md).

required <- c(
  "readr",      # CSV / TSV I/O
  "dplyr",      # data manipulation
  "tidyr",      # pivoting / reshaping
  "stringr",    # text utilities
  "forcats",    # factor utilities
  "ggplot2",    # plotting
  "patchwork",  # multi-panel layout (Supplementary Fig 8)
  "scales",     # log labels
  "vegan",      # community ecology (used by some R figure scripts)
  "RColorBrewer"
)

missing <- setdiff(required, rownames(installed.packages()))
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
message("All required R packages are available.")
sessionInfo()
