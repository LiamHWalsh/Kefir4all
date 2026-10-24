# Figure 6 — full (A + B + C), portable, repo-data only.
# A: Shannon alpha diversity across kefir types (ggbetweenstats with pairwise).
# B: Shannon by fermentation conditions (Laboratory controlled vs Household).
# C: persistence bubble plot of putative environmental species with
#    same-participant (Sample x species) connecting lines.
# + C-alt1: participant-coloured lines (multiple participants separable by colour).
# + C-alt2: faceted by participant (one small panel per participant).
# Inputs:  data/04_metacache_total_species_profile_v2.csv,
#          data/kefir4all_sample_metadata_v2.csv,
#          data/global_milk_kefir_metadata_v1.csv
# Outputs: output/04_taxonomic_profiling/Figure_6A_alpha_kefirtype.pdf,
#          Figure_6B_alpha_condition.pdf,
#          Figure_6C_persistence.pdf (+ _alt1_participant_colour.pdf,
#          _alt2_facet_participant.pdf), CSVs of underlying data.

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(here)
suppressPackageStartupMessages(library(pacman))

DATA_DIR <- here::here("data")
OUT_DIR  <- here::here("output", "04_taxonomic_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "04_taxonomic_profiling")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

pacman::p_load(readr, dplyr, tidyr, tibble, ggplot2, scales, vegan, ggstatsplot)

message("Reading inputs from ", DATA_DIR)
mc_raw <- read_csv(file.path(DATA_DIR, "04_metacache_total_species_profile_v2.csv"),
                   show_col_types = FALSE)
# Strip path prefix: metacache IDs are like '04_Metacache_individual/TG_ID..._S12'
mc_raw$sample_id <- sub("^.*/", "", mc_raw$sample_id)
mc_raw <- mc_raw[mc_raw$sample_id != "TG_EC_S72", ]

kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"),
                               show_col_types = FALSE)
kefir4all_metadata$merge_column <-
  gsub("_host_removed_R..fastq.gz", "", kefir4all_metadata$merge_column)
kefir4all_metadata <-
  kefir4all_metadata[!duplicated(kefir4all_metadata$merge_column), ]

# Numeric compositional block (coerce; NAs -> 0)
num_block <- suppressWarnings(mc_raw[, -1])
num_block[] <- lapply(num_block, function(x) suppressWarnings(as.numeric(x)))
num_block[is.na(num_block)] <- 0

####################################################################
# PANEL A+B: alpha diversity (Shannon) via vegan
####################################################################
global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv"),
                               show_col_types = FALSE)
total_metadata <- rbind(
  dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
  dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`))

# Shannon diversity across all metagenomes
alpha_mat <- as.matrix(num_block)
rownames(alpha_mat) <- mc_raw$sample_id
# drop all-zero rows
keep <- rowSums(alpha_mat) > 0
alpha_mat <- alpha_mat[keep, , drop = FALSE]

my.files_summary_shannon <- diversity(alpha_mat)
my.files_richness <- specnumber(alpha_mat)
my.files_evenness <- my.files_summary_shannon / log(my.files_richness)

alpha_df <- data.frame(
  site = names(my.files_summary_shannon),
  shannon = as.numeric(my.files_summary_shannon),
  richness = as.numeric(my.files_richness),
  pielou = as.numeric(my.files_evenness),
  stringsAsFactors = FALSE)
alpha_df <- merge(alpha_df, total_metadata,
                  by.x = "site", by.y = "merge_column", all.x = TRUE)

# Fill in any unmatched Kefir4All samples
alpha_df$`kefir type`[which(is.na(alpha_df$`kefir type`))] <- "ML"
alpha_df$data_source[which(is.na(alpha_df$`kefir type`))] <- "Walsh et al 2023"
alpha_df <- alpha_df[!(alpha_df$site %in%
  total_metadata$merge_column[which(total_metadata$`kefir type` %in%
    c("Extraction control", "Medium control"))]), ]
alpha_df$data_source[which(is.na(alpha_df$data_source))] <- "Walsh et al 2023"

# Conditions: Laboratory controlled (Walsh/Mortensen) vs Household conditions
alpha_df$conditions <- "Household conditions"
alpha_df$conditions[alpha_df$data_source %in%
  c("Walsh et al 2023", "Mortensen et al 2023")] <- "Laboratory controlled"

# Factor levels for kefir type
alpha_df$`kefir type` <- factor(alpha_df$`kefir type`,
                                levels = c("MG", "ML", "WG", "WL"))

####################################################################
# PANEL A: Shannon by kefir type (ggbetweenstats, all data sources)
####################################################################
panelA <- ggbetweenstats(
  data = alpha_df,
  x = `kefir type`,
  y = shannon,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  ggsignif.args = list(textsize = 4, tip_length = 0.01),
  bf.message = FALSE,
  palette = "Dark2"
) +
  theme_bw() +
  xlab("Kefir type") +
  ylab("Alpha diversity values (Shannon)") +
  theme(plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
        legend.title = element_text(size = 25, face = "bold"),
        axis.text.x = element_text(hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text(size = 15, face = "bold",
                                    hjust = 0.5, vjust = -1),
        axis.title.y = element_text(size = 15, face = "bold",
                                    hjust = 0.5, vjust = 2),
        legend.position = "none") +
  scale_x_discrete(labels = c("Milk Grain", "Milk Liquid",
                              "Water Grain", "Water Liquid"))

####################################################################
# PANEL B: Shannon by condition (Laboratory vs Household), ML/WL only
####################################################################
alpha_df$merge <- paste0(gsub("ML", "Milk Liquid",
                               gsub("WL", "Water Liquid",
                                    alpha_df$`kefir type`)),
                         "\n                   ",
                         alpha_df$conditions)

panelB <- ggbetweenstats(
  data = alpha_df %>% filter(`kefir type` %in% c("ML", "WL")),
  x = merge, y = shannon,
  grouping.var = merge,
  type = "nonparametric",
  plot.type = "box",
  ggsignif.args = list(textsize = 4, tip_length = 0.01),
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
) +
  theme_bw() +
  xlab("Fermentation conditions") +
  ylab("Alpha diversity values (Shannon)") +
  theme(plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
        legend.title = element_text(size = 25, face = "bold"),
        axis.text.x = element_text(hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text(size = 15, face = "bold",
                                    hjust = 0.5, vjust = -1),
        axis.title.y = element_text(size = 15, face = "bold",
                                    hjust = 0.5, vjust = 2),
        legend.position = "none")

# ---- PANEL C data: species table + persistence detection ----
# Species-level table at >0.1% max abundance (as in run_fig6c_persistence.R)
maxab_species <- apply(num_block, 2, max, na.rm = TRUE)
metacache_num <- cbind(sample_id = mc_raw$sample_id, num_block)
species_data <- metacache_num[, c(TRUE, maxab_species >= 0.1)] %>%
  column_to_rownames("sample_id") %>%
  rownames_to_column("sample") %>%
  pivot_longer(!sample, names_to = "clade_name",
               values_to = "relative_abundance")
species_data <- merge(species_data, kefir4all_metadata,
                      by.x = "sample", by.y = "merge_column")
species_data$type <- NA
is_wg_wl <- species_data$sample %in% kefir4all_metadata$merge_column[
  kefir4all_metadata$`kefir type` %in% c("WG", "WL") &
    kefir4all_metadata$Stage != "T0"]
is_mg_ml <- species_data$sample %in% kefir4all_metadata$merge_column[
  kefir4all_metadata$`kefir type` %in% c("MG", "ML") &
    kefir4all_metadata$Stage != "T0"]
is_wg_wl_t0 <- species_data$sample %in% kefir4all_metadata$merge_column[
  kefir4all_metadata$`kefir type` %in% c("WG", "WL") &
    kefir4all_metadata$Stage == "T0"]
is_mg_ml_t0 <- species_data$sample %in% kefir4all_metadata$merge_column[
  kefir4all_metadata$`kefir type` %in% c("MG", "ML") &
    kefir4all_metadata$Stage == "T0"]
species_data$type[is_wg_wl] <- "Water.kefir"
species_data$type[is_mg_ml] <- "Milk.kefir"
species_data$type[is_wg_wl_t0] <- "Water.kefir.T0"
species_data$type[is_mg_ml_t0] <- "Milk.kefir.T0"
species_data <- species_data[!is.na(species_data$type), ]
message("species_data rows: ", nrow(species_data))

# Candidate environmental microbes: absent at T0, present later
comp_df <- data.frame(
  species = as.character(levels(as.factor(species_data$clade_name))),
  check.names = FALSE)
for (cn in c("Milk.kefir", "Milk.kefir.T0", "Water.kefir", "Water.kefir.T0"))
  comp_df[[cn]] <- NA_real_
for (sp in comp_df$species) {
  for (t2 in c("Milk.kefir", "Milk.kefir.T0",
               "Water.kefir", "Water.kefir.T0")) {
    comp_df[comp_df$species == sp, t2] <- nrow(subset(
      species_data, clade_name == sp & type == t2 & relative_abundance > 0.1))
  }
}
soi <- subset(comp_df, Milk.kefir.T0 == 0 & Water.kefir.T0 == 0 &
                Milk.kefir > 0 & Water.kefir > 0) %>%
  dplyr::select(species) %>% mutate(type = "Occurs in both")
t_m <- comp_df[comp_df$Milk.kefir.T0 == 0 & comp_df$Milk.kefir > 0, ] %>%
  dplyr::select(species) %>% mutate(type = "Milk.kefir")
t_m <- t_m[!(t_m$species %in% soi$species), ]
t_w <- comp_df[comp_df$Water.kefir.T0 == 0 & comp_df$Water.kefir > 0, ] %>%
  dplyr::select(species) %>% mutate(type = "Water.kefir")
t_w <- t_w[!(t_w$species %in% soi$species), ]
species_of_interest <- rbind(soi, t_m, t_w)
message("Candidate environmental microbes: ", nrow(species_of_interest))

# Persistence: same participant (Sample) at consecutive stages in grain
persist_df <- data.frame(sample = character(), clade_name = character(),
  relative_abundance = numeric(), Sample = character(),
  `kefir type` = character(), Stage = character(),
  data_source = character(), type = character(), check.names = FALSE)
for (sp in species_of_interest$species) {
  t2 <- species_of_interest$type[species_of_interest$species == sp]
  if (t2 == "Occurs in both") {
    t2 <- c("Milk.kefir", "Water.kefir"); grain_type <- c("MG", "WG")
  } else if (t2 == "Water.kefir") {
    grain_type <- "WG"
  } else {
    grain_type <- "MG"
  }
  for (grains in levels(as.factor(grain_type))) {
    if (nrow(subset(species_data, clade_name == sp & type == t2 &
                    relative_abundance > .1 &
                    `kefir type` == grains)) == 0) next
    soi_rows <- subset(species_data, clade_name == sp & type == t2 &
                         relative_abundance > .1) %>%
      dplyr::select(sample, clade_name, relative_abundance, Sample,
                    `kefir type`, Stage, data_source, type)
    if (TRUE %in% duplicated(soi_rows$Sample)) {
      soi_rows <- soi_rows[soi_rows$Sample %in%
                             soi_rows$Sample[duplicated(soi_rows$Sample)], ]
      soi_red <- soi_rows[soi_rows$`kefir type` == grains, ]
      dup_s <- soi_red$Sample[duplicated(soi_red$Sample)]
      if (length(dup_s) == 0) next
      st <- soi_red$Stage[soi_red$Sample == dup_s[1]]
      tn <- suppressWarnings(as.numeric(gsub("T", "", st)))
      if (length(tn) >= 2 && TRUE %in% ((tn[1] + 1) == tn[2]))
        persist_df <- rbind(soi_rows, persist_df)
    }
  }
}
message("Persistence rows: ", nrow(persist_df))
# ---- END PANEL C data ----

# ---- PANEL C plots: reuse persist_df (built above) ----
# Labels: MG -> Milk Grain, WG -> Water Grain; T1..T6 -> wk1..wk21
persist_df[["kefir type"]][persist_df[["kefir type"]] == "MG"] <- "Milk Grain"
persist_df[["kefir type"]][persist_df[["kefir type"]] == "WG"] <- "Water Grain"
persist_df$stage_label <- persist_df$Stage
persist_df$stage_label[persist_df$Stage == "T1"] <- "wk1"
persist_df$stage_label[persist_df$Stage == "T2"] <- "wk5"
persist_df$stage_label[persist_df$Stage == "T3"] <- "wk9"
persist_df$stage_label[persist_df$Stage == "T4"] <- "wk13"
persist_df$stage_label[persist_df$Stage == "T5"] <- "wk17"
persist_df$stage_label[persist_df$Stage == "T6"] <- "wk21"
persist_df$stage_label <- factor(persist_df$stage_label,
  levels = c("wk0", "wk1", "wk5", "wk9", "wk13", "wk17", "wk21"))
persist_df$kefir_type_lab <- persist_df[["kefir type"]]
persist_df$n_tp <- ave(seq_len(nrow(persist_df)),
                       persist_df$Sample, persist_df$clade_name,
                       FUN = length)
persist_df$traj_id <- paste(persist_df$Sample, persist_df$clade_name, sep = " | ")
# n_participants: how many unique Samples share this clade_name (sharing breadth)
clade_sharing <- tapply(persist_df$Sample, persist_df$clade_name,
                        function(x) length(unique(x)))
persist_df$n_participants <- as.integer(clade_sharing[persist_df$clade_name])
# ---- END PANEL C PREP ----

panelC_theme <- theme_bw() +
  theme(legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(size = 17.5),
    axis.text.y = element_text(size = 13.5, face = "italic"),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(size = 20),
    axis.title.y = element_text(size = 20, face = "bold",
                                hjust = 0.5, vjust = 1.5),
    legend.text = element_text(size = 20, face = "italic"),
    legend.title = element_text(size = 20, vjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))

# Panel C: same-participant lines coloured by number of participants
# sharing the clade_name (sharing breadth), with a sparse integer legend
n_max <- max(persist_df$n_participants, na.rm = TRUE)
n_breaks <- unique(round(seq(1, n_max, length.out = min(n_max, 5))))
panelC <- ggplot(persist_df[order(persist_df$traj_id), ],
                 aes(x = stage_label, y = clade_name)) +
  geom_line(aes(group = traj_id, colour = n_participants),
            alpha = 0.85, linewidth = 0.9) +
  geom_point(aes(size = relative_abundance, fill = relative_abundance),
             shape = 21, colour = "grey20", alpha = 0.9) +
  scale_size(range = c(2, 20), name = "Relative abundance (%)") +
  scale_fill_gradientn(
    colours = c("Dark blue", "yellow", "yellow2", "orange", "red", "darkred"),
    name = "Relative abundance (%)") +
  scale_colour_gradientn(
    colours = c("#30123B", "#1AE4B6", "#F9C74F", "#F3722C", "#F94144"),
    breaks = n_breaks,
    labels = n_breaks,
    name = "Number of\nparticipants") +
  labs(x = "", y = "Species") +
  facet_wrap(~kefir_type_lab, dir = "v", scales = "free") +
  panelC_theme +
  guides(fill = guide_colourbar(label.position = "bottom"),
         colour = guide_colourbar(label.position = "bottom",
                                  barwidth = 25, barheight = 1.6,
                                  ticks.linewidth = 1),
         size = guide_legend(label.position = "bottom"))
# ---- END PANEL C ----

# ---- OUTPUT ----
ggsave(file.path(FIGURES_DIR, "Figure_6A_alpha_kefirtype.tiff"),
       panelA, width = 6, height = 4, dpi = 300, compression = "lzw")
ggsave(file.path(FIGURES_DIR, "Figure_6B_alpha_condition.tiff"),
       panelB, width = 6, height = 4, dpi = 300, compression = "lzw")
ggsave(file.path(FIGURES_DIR, "Figure_6C_persistence.tiff"),
       panelC, width = 6, height = 4, dpi = 300, compression = "lzw")
# ---- END PANEL C OUTPUT ----





# Combined Figure 6 (A + B + C)
library(ggpubr)
fig6_combined <- ggarrange( ggarrange(panelA, panelB,
                                      nrow=1,ncol=2, labels=c("A.","B."),common.legend = FALSE),
panelC,
                           labels = c("", "C"),
                           ncol = 1, nrow = 2,
                           heights = c(.4, 1))
ggexport(fig6_combined, filename = file.path(FIGURES_DIR, "Figure_environmental_microbes.pdf"),
         width = 18, height = 28)

# Underlying data CSVs
write.csv(alpha_df, file.path(OUT_DIR, "Figure_environmental_microbes_B_alpha_data.csv"),
          row.names = FALSE)
write.csv(persist_df, file.path(OUT_DIR, "Figure_environmental_microbes_persistence_data.csv"),
          row.names = FALSE)

message("All done :) outputs written to ", OUT_DIR)


