# Publication-style ggplot2 and ggstatsplot figures that can be regenerated and checked.
#
# Every figure is saved as a vector PDF plus a PNG preview, with a sidecar JSON
# file recording the source data hash, seed, package versions and sizes. The
# theme uses one text size and no plot title. Statistics are drawn on the figure,
# not in a caption.

WONG_PALETTE <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7",
                  "#D55E00", "#F0E442", "#0072B2", "#000000")

#' Return a ggplot2 theme with one text size and no plot title
#'
#' @param text_size_pt Numeric. Size used for every text element: axis titles,
#'   tick labels, legend text, strip labels and the subtitle. The title is blank.
#' @param base_family Character. Font family. Times New Roman is the default for
#'   manuscripts; it needs to be installed on the machine that renders the PDF.
#'
#' @return A ggplot2 theme. Add it with `plot + figure_theme()`.
#'
#' @details Statistics belong on the figure, in the subtitle or the plotted
#'   annotations, not in a caption. Mapping colours is the caller's job. Use
#'   `WONG_PALETTE` for up to eight discrete levels and a viridis scale for
#'   continuous values.
figure_theme <- function(text_size_pt = 8, base_family = "Times New Roman") {
  ggplot2::theme_bw(base_size = text_size_pt, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = text_size_pt),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_text(size = text_size_pt, hjust = 0),
      legend.title = ggplot2::element_text(size = text_size_pt),
      legend.text = ggplot2::element_text(size = text_size_pt),
      axis.title = ggplot2::element_text(size = text_size_pt),
      axis.text = ggplot2::element_text(size = text_size_pt),
      strip.text = ggplot2::element_text(size = text_size_pt),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Commit and dirty flag of the repository containing the working directory.
# Outside a git repository, or without git installed, both fields are NA. That
# value is written to the provenance file, so the missing state is recorded.
git_state <- function() {
  not_available <- list(commit = NA_character_, dirty = NA)
  if (!nzchar(Sys.which("git"))) return(not_available)
  run <- function(args) {
    out <- suppressWarnings(system2("git", args, stdout = TRUE, stderr = FALSE))
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) return(NULL)
    out
  }
  commit <- run(c("rev-parse", "HEAD"))
  if (is.null(commit)) return(not_available)
  status <- run(c("status", "--porcelain"))
  # A failed status is recorded as unknown (NA), never as a clean tree.
  dirty <- if (is.null(status)) NA else length(status) > 0
  list(commit = commit[1], dirty = dirty)
}

#' Save a plot as PDF and PNG and write a provenance JSON file
#'
#' @param plot A ggplot object with `figure_theme()` already added.
#' @param path Character. Output path without extension, for example
#'   `"outputs/fig1"`. Writes `fig1.pdf`, `fig1.png` and `fig1.json`.
#' @param width_mm,height_mm Numeric. Physical size. Use 85 mm for a single
#'   column and 180 mm for a double column. Figures with statistics on them
#'   need 180 mm to stay readable.
#' @param source Character or NULL. Input data file. Its SHA-256 goes in the
#'   provenance file.
#' @param seed Integer or NULL. Seed used for any randomness in the data
#'   preparation.
#' @param overwrite Logical. If FALSE (the default), refuse to replace existing
#'   outputs.
#' @param extra Named list or NULL. Extra provenance fields, such as the
#'   transformations applied.
#'
#' @return The provenance list that was written, invisibly.
#'
#' @details The PNG is a preview at 300 dpi. The PDF is the file to submit.
save_figure <- function(plot, path, width_mm = 85, height_mm = 60,
                        source = NULL, seed = NULL, overwrite = FALSE,
                        extra = NULL) {
  base <- sub("\\.[^./]*$", "", path)
  outputs <- paste0(base, c(".pdf", ".png", ".json"))
  if (!overwrite && any(file.exists(outputs))) {
    stop("Refusing to overwrite: ",
         paste(outputs[file.exists(outputs)], collapse = ", "),
         ". Pass overwrite = TRUE.", call. = FALSE)
  }
  # Hash the source before writing anything, so a bad path leaves no partial output.
  source_sha256 <- if (is.null(source)) {
    NULL
  } else {
    if (!file.exists(source)) {
      stop("Source file not found: ", source, call. = FALSE)
    }
    digest::digest(file = source, algo = "sha256")
  }
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)

  # Write every file under a temporary name first, then move them into place with
  # the JSON last. A failure part-way leaves no file that looks complete.
  temporary <- paste0(dirname(base), "/.", basename(base),
                      c(".pdf", ".png", ".json"), ".tmp")
  finished <- FALSE
  on.exit(if (!finished) unlink(temporary[file.exists(temporary)]), add = TRUE)

  ggplot2::ggsave(temporary[1], plot, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::cairo_pdf)
  # The PNG uses the same cairo engine as the PDF, so the two match.
  cairo_png <- function(filename, width, height, units = "in", dpi = 300, ...) {
    grDevices::png(filename, width = width, height = height, units = units,
                   res = dpi, type = "cairo")
  }
  ggplot2::ggsave(temporary[2], plot, width = width_mm, height = height_mm,
                  units = "mm", dpi = 300, device = cairo_png)

  state <- git_state()
  record <- list(
    figure = basename(outputs[1]),
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    width_mm = width_mm,
    height_mm = height_mm,
    seed = seed,
    source_file = source,
    source_sha256 = source_sha256,
    lab_toolkit_commit = state$commit,
    lab_toolkit_uncommitted_changes = state$dirty,
    versions = list(
      R = as.character(getRversion()),
      ggplot2 = as.character(utils::packageVersion("ggplot2")),
      ggstatsplot = if (requireNamespace("ggstatsplot", quietly = TRUE)) {
        as.character(utils::packageVersion("ggstatsplot"))
      } else {
        NA_character_
      }
    )
  )
  if (!is.null(extra)) record <- utils::modifyList(record, extra)
  jsonlite::write_json(record, temporary[3], auto_unbox = TRUE, pretty = TRUE,
                       null = "null")
  for (i in seq_along(outputs)) file.rename(temporary[i], outputs[i])
  finished <- TRUE
  invisible(record)
}

#' Between-groups comparison with ggstatsplot, styled for the manuscript
#'
#' @param data Data frame.
#' @param x,y Unquoted column names: grouping factor and numeric outcome.
#' @param seed Integer. Set before the call, because ggstatsplot bootstraps
#'   some statistics.
#' @param text_size_pt Numeric. Passed to `figure_theme()`.
#' @param x_label,y_label Character or NULL. Axis titles. NULL uses the column
#'   names, which are often not fit for a figure.
#' @param level_labels Named character vector or NULL. Names are the data values
#'   as stored, values are the labels to draw, in the order to plot them, for
#'   example `c(control = "Control", treated_high = "Treated high")`. Every
#'   value in the data must be named here.
#'
#' @return A list with:
#'   `plot`: the ggplot, with ggstatsplot's statistics on it as the subtitle,
#'   pairwise brackets and Bayesian line.
#'   `stats`: the ANOVA results (`extract_stats()$subtitle_data`).
#'   `bayes`: the Bayes factor and Bayesian R-squared table.
#'   `pairwise`: the pairwise comparisons.
#'   The `expression` columns are dropped from each table, because they hold
#'   plotmath calls that cannot be written to CSV or converted to Python.
#'
#' @details One ggstatsplot call fits the model, runs the tests and draws the
#'   figure, so the statistics and the figure come from the same fit. Titles
#'   are removed. The legend is removed too, because the x-axis already names
#'   each group. Use a 180 mm width when saving.
ggstat_between <- function(data, x, y, seed, text_size_pt = 8,
                           x_label = NULL, y_label = NULL, level_labels = NULL) {
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  x_name <- rlang::as_name(x)
  y_name <- rlang::as_name(y)

  if (!is.null(level_labels)) {
    values <- as.character(data[[x_name]])
    unlabelled <- setdiff(unique(values), names(level_labels))
    if (length(unlabelled) > 0) {
      stop("level_labels has no label for: ",
           paste(unlabelled, collapse = ", "), call. = FALSE)
    }
    data[[x_name]] <- factor(values, levels = names(level_labels),
                             labels = unname(level_labels))
  }

  present <- data[!is.na(data[[y_name]]), ]
  counts <- table(as.character(present[[x_name]]))
  if (length(counts) < 2) {
    stop("ggstat_between() needs at least two groups with data.", call. = FALSE)
  }
  if (any(counts < 2)) {
    stop("Each group needs at least two non-missing values; too few in: ",
         paste(names(counts)[counts < 2], collapse = ", "), call. = FALSE)
  }
  if (stats::var(present[[y_name]]) == 0) {
    stop("The outcome has zero variance, so the tests cannot be computed.",
         call. = FALSE)
  }

  xlab <- if (is.null(x_label)) x_name else x_label
  ylab <- if (is.null(y_label)) y_name else y_label

  set.seed(seed)
  built <- ggstatsplot::ggbetweenstats(
    data = data, x = !!x, y = !!y, xlab = xlab, ylab = ylab,
    title = NULL, caption = NULL,
    results.subtitle = TRUE, bf.message = TRUE,
    pairwise.display = "all"
  )

  tables <- ggstatsplot::extract_stats(built)
  drop_expression <- function(table) {
    table <- as.data.frame(table)
    table$expression <- NULL
    table
  }

  figure <- built + ggplot2::labs(title = NULL) + figure_theme(text_size_pt) +
    ggplot2::theme(legend.position = "none")
  list(
    plot = figure,
    stats = drop_expression(tables$subtitle_data),
    bayes = drop_expression(tables$caption_data),
    pairwise = drop_expression(tables$pairwise_comparisons_data)
  )
}
