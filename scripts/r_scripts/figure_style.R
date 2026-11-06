# Shared figure style for manuscript, extended-data and supplementary figures.

KEFIR_PALETTE <- c(
  vermillion = "#D55E00",
  sky_blue = "#56B4E9",
  bluish_green = "#009E73",
  reddish_purple = "#CC79A7",
  orange = "#E69F00",
  blue = "#0072B2",
  yellow = "#F0E442",
  black = "#000000"
)

kefir_palette <- function(n) {
  if (n <= length(KEFIR_PALETTE)) {
    return(unname(KEFIR_PALETTE[seq_len(n)]))
  }
  grDevices::hcl.colors(n, palette = "Dark 3")
}

kefir_theme <- function(base_size = 14, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = base_size + 4, face = "bold", hjust = 0.5, colour = "black"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size, hjust = 0.5, colour = "black"
      ),
      plot.tag = ggplot2::element_text(
        size = base_size + 6, face = "bold", colour = "black"
      ),
      strip.text = ggplot2::element_text(
        size = base_size + 4, face = "bold", colour = "black", hjust = 0.5
      ),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(
        size = base_size + 2, face = "plain", colour = "black"
      ),
      axis.text = ggplot2::element_text(size = base_size, colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(
        size = base_size + 1, face = "bold", colour = "black"
      ),
      legend.text = ggplot2::element_text(size = base_size, colour = "black"),
      legend.key = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

kefir_use_style <- function(base_size = 14) {
  ggplot2::theme_set(kefir_theme(base_size = base_size))
  invisible(TRUE)
}
