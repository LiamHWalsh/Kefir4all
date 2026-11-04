# Produces: Supplementary Figure 1 -- "(A-B) Beta-diversity boxplots displaying the Bray
# Curtis distance values as determined based on functional composition, across (A) milk
# grain metagenomes and (B) water grain metagenomes, representing different timeframe of
# the Kefir4All study (wk01, wk05, wk09, wk13, wk17 and wk21) compared with the
# corresponding T0 milk or water kefir grain profile" (Kefir4All_Manuscript_TRACKED_SB.docx).
#
# No script in this repo plots this (community_changes_v3.R only renders a trajectory-line
# version combining all types/sources). The underlying per-sample distance-to-T0 values it
# already computes and writes to output/04_taxonomic_profiling/community_changes_v3/ are
# reused here as a boxplot, grain types only, one panel per kefir type.

for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {
  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
kefir_use_style()
CC3_DIR     <- here::here("output", "04_taxonomic_profiling", "community_changes_v3")
FIGURES_DIR <- here::here("figures", "04_taxonomic_profiling")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

library(readr); library(dplyr); library(ggplot2)

stage_to_week <- c(T0 = "T0", T1 = "wk01", T2 = "wk05", T3 = "wk09", T4 = "wk13", T5 = "wk17", T6 = "wk21")
stage_levels  <- c("wk01", "wk05", "wk09", "wk13", "wk17", "wk21")

read_reference <- function(type, label) {
  path <- file.path(CC3_DIR, paste0("SUPER_FOCUS_", type, "_reference_distances.csv"))
  if (!file.exists(path)) stop("Missing: ", path, " -- run community_changes_v3.R first.")
  d <- read_csv(path, show_col_types = FALSE)
  d$week <- factor(stage_to_week[d$Stage], levels = stage_levels)
  d$type_label <- label
  d
}

plot_data <- rbind(
  read_reference("MG", "Milk grain"),
  read_reference("WG", "Water grain")
)
plot_data$type_label <- factor(plot_data$type_label, levels = c("Milk grain", "Water grain"))

p <- ggplot(plot_data, aes(x = week, y = distance, fill = type_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  facet_wrap(~ type_label, nrow = 1) +
  scale_fill_manual(values = c("Milk grain" = "#D55E00", "Water grain" = "#009E73"), guide = "none") +
  labs(x = "Fermentation stage (compared to T0 of the same type)",
       y = "Bray-Curtis dissimilarity\n(SUPER-FOCUS functional composition, vs. T0)") +
  theme_bw() +
  theme(strip.text = element_text(size = 22, face = "bold", colour = "black"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"))

ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_1_superfocus_t0.tiff"), p,
       width = 12, height = 6, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Supplementary_Figure_1_superfocus_t0.pdf"), p,
       width = 12, height = 6, bg = "white")

message("Wrote Supplementary_Figure_1_superfocus_t0.{tiff,pdf} to ", FIGURES_DIR)
