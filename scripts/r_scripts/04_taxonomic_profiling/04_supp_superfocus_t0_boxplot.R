# Produces Supplementary Figure 1: T0-referenced SUPER-FOCUS Bray-Curtis distances
# (panels A-B) and strong subsystem-level correlations with distance (panel C).

for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {
  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style()

DATA_DIR <- here::here("data")
FIGURES_DIR <- here::here("figures", "Supplementary_Figure_1")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ggstatsplot)
library(ggpubr)
library(ggforce)

normalize_sample_id <- function(value) {
  value <- basename(gsub("\\\\", "/", trimws(as.character(value))))
  value <- sub("_host_removed_R[0-9]\\.fastq\\.gz$", "", value)
  gsub("-", "_", value, fixed = TRUE)
}

metadata <- read_csv(
  file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"),
  show_col_types = FALSE
)
metadata$key <- normalize_sample_id(metadata$merge_column)
metadata <- metadata[!duplicated(metadata$key) &
                       metadata[["kefir type"]] %in% c("MG", "WG") &
                       metadata$Stage %in% paste0("T", 0:6), ]

superfocus <- read_csv(
  file.path(DATA_DIR, "functional", "superfocus_p3.csv"),
  show_col_types = FALSE
)
required_columns <- c("sample_id", "Subsystem Level 3", "reads_%")
if (!all(required_columns %in% names(superfocus))) {
  stop("SUPER-FOCUS input is missing required columns: ",
       paste(setdiff(required_columns, names(superfocus)), collapse = ", "))
}
superfocus$key <- normalize_sample_id(superfocus$sample_id)
superfocus$feature <- trimws(superfocus[["Subsystem Level 3"]])
superfocus$abundance <- suppressWarnings(as.numeric(superfocus[["reads_%"]]))
if (anyNA(superfocus$key) || any(!nzchar(superfocus$key)) ||
    anyNA(superfocus$feature) || any(!nzchar(superfocus$feature)) ||
    anyNA(superfocus$abundance) || any(!is.finite(superfocus$abundance)) ||
    any(superfocus$abundance < 0)) {
  stop("SUPER-FOCUS input contains invalid sample IDs, pathway labels, or abundances.")
}

missing_profiles <- setdiff(metadata$key, unique(superfocus$key))
if (length(missing_profiles)) {
  message("Excluding ", length(missing_profiles),
          " metadata sample(s) without SUPER-FOCUS profiles: ",
          paste(head(missing_profiles, 10), collapse = ", "))
  metadata <- metadata[metadata$key %in% unique(superfocus$key), ]
}

profile_matrix <- xtabs(abundance ~ key + feature, data = superfocus)
distance_matrix <- as.matrix(vegan::vegdist(profile_matrix, method = "bray"))
stage_to_week <- c(T0 = "T0", T1 = "wk01", T2 = "wk05", T3 = "wk09",
                   T4 = "wk13", T5 = "wk17", T6 = "wk21")
stage_levels <- unname(stage_to_week[paste0("T", 1:6)])

reference_distances <- list()
reference_labels <- c(MG = "Milk kefir", WG = "Water kefir")
for (type in names(reference_labels)) {
  type_metadata <- metadata[metadata[["kefir type"]] == type, ]
  reference_key <- type_metadata$key[type_metadata$Stage == "T0"][1]
  if (is.na(reference_key) || !reference_key %in% rownames(distance_matrix)) {
    stop("No T0 SUPER-FOCUS reference profile for ", reference_labels[[type]], ".")
  }
  reference_distances[[type]] <- data.frame(
    key = type_metadata$key,
    Stage = type_metadata$Stage,
    distance = as.numeric(distance_matrix[reference_key, type_metadata$key]),
    kefir_type = reference_labels[[type]],
    stringsAsFactors = FALSE
  )
}

distance_data <- bind_rows(reference_distances)
distance_data$week <- factor(stage_to_week[distance_data$Stage],
                             levels = stage_levels)
plot_data <- distance_data[distance_data$Stage != "T0", ]
plot_data$kefir_type <- factor(plot_data$kefir_type,
                               levels = c("Milk kefir", "Water kefir"))

plot_reference_distances <- function(type_label) {
  ggbetweenstats(
    data = plot_data[plot_data$kefir_type == type_label, ],
    x = week,
    y = distance,
    title = type_label,
    type = "nonparametric",
    ggsignif.args = list(textsize = 4, tip_length = 0.01)
  ) +
    theme_bw() +
    xlab("Time points") +
    ylab("Bray Curtis distance") +
    theme(
      legend.position = "none",
      plot.title = ggtext::element_textbox_simple(
        face = "bold", halign = 0.5, linetype = 1,
        box.color = "#748696", size = 55, lineheight = 2
      ),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
      axis.text.y = element_text(hjust = 1, size = 10),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold")
    )
}

panel_a <- plot_reference_distances("Milk kefir")
panel_b <- plot_reference_distances("Water kefir")

correlation_rows <- list()
row_index <- 0L
for (type_label in levels(plot_data$kefir_type)) {
  type_distances <- distance_data[distance_data$kefir_type == type_label, ]
  # relationship = "one-to-many": type_distances has one row per sample (key), superfocus
  # has one row per sample x pathway -- an unexpected duplicate sample in type_distances
  # should fail loudly rather than silently duplicating correlation inputs.
  matched <- inner_join(type_distances, superfocus, by = "key", relationship = "one-to-many")
  missing_correlation_profiles <- setdiff(type_distances$key, unique(matched$key))
  if (length(missing_correlation_profiles)) {
    stop("SUPER-FOCUS correlation join lost samples for ", type_label, ": ",
         paste(head(missing_correlation_profiles, 10), collapse = ", "))
  }
  for (pathway in unique(matched$feature)) {
    pathway_data <- matched[matched$feature == pathway, ]
    if (nrow(pathway_data) < 4 ||
        length(unique(pathway_data$abundance)) < 2 ||
        length(unique(pathway_data$distance)) < 2) {
      next
    }
    test <- suppressWarnings(cor.test(
      pathway_data$abundance, pathway_data$distance,
      method = "spearman", exact = FALSE
    ))
    row_index <- row_index + 1L
    correlation_rows[[row_index]] <- data.frame(
      kefir_type = type_label,
      pathway = pathway,
      r = unname(test$estimate),
      p_value = test$p.value
    )
  }
}

correlation_data <- bind_rows(correlation_rows)
if (!nrow(correlation_data)) {
  stop("No valid SUPER-FOCUS pathway correlations were computed.")
}
correlation_data$adjusted_p <- p.adjust(correlation_data$p_value,
                                         method = "bonferroni")
# This cutoff reproduces the manuscript panel's sparse milk and larger water bar groups.
correlation_data <- correlation_data[
  correlation_data$adjusted_p <= 0.05 &
    abs(correlation_data$r) > 0.82, ]
if (!nrow(correlation_data)) {
  stop("No strong, statistically significant pathway correlations were found.")
}
correlation_data$pathway <- gsub("EC 5.1.2.- ", "",
                                  correlation_data$pathway)
correlation_data$kefir_type <- factor(
  correlation_data$kefir_type,
  levels = c("Milk kefir", "Water kefir")
)

panel_c <- ggplot(
  correlation_data,
  aes(y = reorder(pathway, r), x = r)
) +
  geom_col(fill = "coral1", color = "black", width = 1, linewidth = 0.3) +
  ggforce::facet_col(kefir_type ~ ., scales = "free", space = "free") +
  scale_x_continuous(limits = c(-1, 1)) +
  theme_bw() +
  labs(
    y = "Subsystem",
    x = "Spearman's rank correlation coefficient"
  ) +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(hjust = 1, size = 7),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )

figure <- ggarrange(
  ggarrange(panel_a, panel_b, ncol = 2, nrow = 1,
            labels = c("A.", "B."),
            font.label = list(size = 18)),
  panel_c,
  nrow = 2, ncol = 1, labels = c("", "C."),
  common.legend = FALSE, font.label = list(size = 18),
  heights = c(0.6, 1)
)

ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_1_superfocus_t0.tiff"),
       figure, width = 7864 / 300, height = 5200 / 300, dpi = 300,
       compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_1_superfocus_t0.pdf"),
       figure, width = 7864 / 300, height = 5200 / 300, bg = "white")

save_figure(figure, file.path(FIGURES_DIR, "Supplementary_Figure_1_superfocus_t0_provenance"),
            width_mm = (7864 / 300) * 25.4, height_mm = (5200 / 300) * 25.4,
            source = file.path(DATA_DIR, "functional", "superfocus_p3.csv"),
            seed = NULL)

message("Wrote three-panel Supplementary Figure 1 to ", FIGURES_DIR)
message("Strong correlations plotted: ",
        paste(names(table(correlation_data$kefir_type)),
              as.integer(table(correlation_data$kefir_type)),
              collapse = "; "))
