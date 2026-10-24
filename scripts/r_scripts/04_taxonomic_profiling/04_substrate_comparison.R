# ==============================================================================
# 04_substrate_comparison.R
#
# Effect of fermentation substrate on the kefir microbiome, function and
# volatiles. Addresses the reviewer comment: citizen scientists chose different
# substrates (milk: whole / low-fat / skim / goat; water: sugar type +/- dried
# fruit). How many arms were studied, and does substrate shape the microbiome?
#
# Two complementary datasets are analysed:
#   (i)  Citizen-science cohort (TG_ID###_...): substrate per participant from
#        the private citizen-scientist metadata (category_confirmed).
#   (ii) Controlled laboratory time-course (TG_<substrate code>_...): substrate
#        is encoded in the sample ID (milk: SK skim, LO low-fat, GO goat,
#        FL/FU whole-milk variants; water: BR brown sugar, AP white+apricot,
#        FG white+fig ...), including anaerobic (0AN) vs aerobic (0AE) arms.
#
# Analyses
#   1. Substrate arm census + substrate stability (did participants switch?)
#   2. Beta diversity (Bray-Curtis): adonis2 ~ substrate (+ covariates),
#      permutations restricted within participant/lab arm; betadisper check;
#      BH-corrected pairwise contrasts
#   3. Alpha diversity (Shannon, richness) across arms
#   4. Volatile profiles (published water-kefir datasets): PERMANOVA of
#      fermentation time vs vessel
#   5. Functional/resistome: computed only if the resistome tables exist
#      (they are not part of the distributed repository)
#
# Output: output/04_substrate_comparison/ (tables, PDF + TIFF figure)
# ==============================================================================

# ---- Path configuration (portable, anchors at repo root) ---------------------
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
here::i_am("scripts/r_scripts/04_taxonomic_profiling/04_substrate_comparison.R")

DATA_DIR    <- here::here("data")
PRIVATE_DIR <- here::here("data", "private")
OUT_DIR     <- here::here("output", "04_substrate_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "04_substrate_comparison")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Libraries ---------------------------------------------------------------
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ggpubr)

cat("=== SUBSTRATE COMPARISON ANALYSIS ===\n")

# ============================================================================
# 1. LOAD DATA + PARSE SAMPLE KEYS
# ============================================================================
cat("Loading data...\n")

metacache <- read_csv(file.path(DATA_DIR, "04_metacache_total_species_profile_v2.csv"),
                      show_col_types = FALSE)
cs_meta <- read.csv(file.path(PRIVATE_DIR, "Citizen Scientist metadata_v8.csv"),
                    check.names = FALSE)

# metacache sample IDs: strip folder prefix, normalise dashes (same as Figure 2)
metacache$sample_key <- gsub(".*/", "", gsub("-", "_", metacache$sample_id))

key_info <- tibble(sample_key = metacache$sample_key) %>%
  mutate(
    tg          = grepl("^TG_", sample_key),
    token1      = ifelse(tg, sub("_.*$", "", sub("^TG_", "", sample_key)), NA_character_),
    kind        = case_when(
      grepl("_MG_", sample_key) ~ "MG",
      grepl("_ML_|_MLC", sample_key) ~ "ML",
      grepl("_WG_", sample_key) ~ "WG",
      grepl("_WL_|_WLC", sample_key) ~ "WL",
      TRUE ~ NA_character_),
    compartment = case_when(kind %in% c("MG", "WG") ~ "Grain",
                            kind %in% c("ML", "WL") ~ "Liquid",
                            TRUE ~ NA_character_),
    kefir_kind  = case_when(kind %in% c("MG", "ML") ~ "Milk",
                            kind %in% c("WG", "WL") ~ "Water",
                            TRUE ~ NA_character_),
    is_citizen  = tg & grepl("^ID", token1),
    o2          = case_when(grepl("0AN", sample_key) ~ "Anaerobic",
                            grepl("0AE", sample_key) ~ "Aerobic",
                            TRUE ~ NA_character_),
    stage       = {
      st <- str_extract(sample_key, "_T[0-6]_")
      ifelse(!is.na(st), gsub("_", "", st), NA_character_)})

# participant ID padding (same convention as Figure 2: ID -> ID0/ID00)
pad_id <- function(x) {
  x <- as.character(x); out <- x
  out[nchar(x) == 3] <- gsub("ID", "ID00", x[nchar(x) == 3])
  out[nchar(x) == 4] <- gsub("ID", "ID0",  x[nchar(x) == 4])
  out
}
cs_meta$ID_pad <- pad_id(cs_meta$ID)

substrate_map <- cs_meta %>%
  transmute(participant = ID_pad,
            substrate   = category_confirmed,
            ferment_type = `Fermentation type`,
            ferment_freq = category,
            part_num    = as.numeric(sub("^ID0*", "", ID_pad)))

lab_labels <- c(SK = "Skim milk", LO = "Low-fat milk", GO = "Goat milk",
                FL = "Full-fat milk", FU = "Whole-milk variant",
                BR = "Brown sugar", AP = "White sugar + dried apricot",
                FG = "White sugar + dried fig", WA = "White sugar")

key_info <- key_info %>%
  mutate(part_num = ifelse(is_citizen, as.numeric(sub("^ID0*", "", token1)), NA_real_)) %>%
  left_join(substrate_map, by = "part_num") %>%
  mutate(source    = case_when(is_citizen ~ "Citizen", tg ~ "Lab", TRUE ~ "Other"),
         substrate = case_when(is_citizen ~ substrate,
                               tg & token1 %in% names(lab_labels) ~
                                 unname(lab_labels[token1]),
                               tg ~ token1,
                               TRUE ~ NA_character_),
         participant = ifelse(is_citizen, token1, paste0("lab_", token1)))

metacache$sample_id <- NULL
species_wide <- metacache %>% column_to_rownames("sample_key")
species_wide[is.na(species_wide)] <- 0
# ============================================================================
# 2. SUBSTRATE ARM CENSUS + STABILITY
# ============================================================================
cat("\n--- 2. Substrate arms ---\n")

cit <- key_info %>% filter(source == "Citizen", !is.na(substrate))
lab <- key_info %>% filter(source == "Lab", !is.na(substrate))

cat("Citizen participants per milk arm:\n")
print(table(cit$substrate[cit$kefir_kind == "Milk"]))
cat("Citizen participants per water arm:\n")
print(table(cit$substrate[cit$kefir_kind == "Water"]))

arm_census <- bind_rows(
  cit %>% count(source = "Citizen", kefir_kind, compartment, substrate),
  lab %>% count(source = "Lab", kefir_kind, compartment, substrate))
write_csv(arm_census, file.path(OUT_DIR, "sample_counts_by_substrate.csv"))
cat("\nSample counts by arm:\n"); print(arm_census)

# Substrate stability: did any participant change substrate during the study?
mk_survey <- read.delim(file.path(PRIVATE_DIR, "fermentation_survey_mk_edited_dated.txt"),
                        sep = "\t", check.names = FALSE)
wk_survey <- read.delim(file.path(PRIVATE_DIR, "fermentation_survey_wk_edited_dated.txt"),
                        sep = "\t", check.names = FALSE)
mk_changes <- mk_survey %>% group_by(Participant_ID) %>%
  summarise(n_arms = n_distinct(Category, na.rm = TRUE)) %>% filter(n_arms > 1)
wk_changes <- wk_survey %>% group_by(Participant_ID) %>%
  summarise(n_arms = n_distinct(sugar_type, na.rm = TRUE)) %>% filter(n_arms > 1)
arm_stability <- tibble(kefir = c("Milk", "Water"),
                        n_participants = c(n_distinct(mk_survey$Participant_ID),
                                           n_distinct(wk_survey$Participant_ID)),
                        n_changed_substrate = c(nrow(mk_changes), nrow(wk_changes)))
cat("\nSubstrate stability (survey rows per participant):\n")
print(arm_stability)
write_csv(arm_stability, file.path(OUT_DIR, "substrate_stability.csv"))

# ============================================================================
# 3. PERMANOVA HELPER
# ============================================================================
cat("\n--- 3. Beta diversity (Bray-Curtis, adonis2) ---\n")

run_permanova <- function(df, label, rhs) {
  df <- df %>% filter(sample_key %in% rownames(species_wide))
  rhs_terms <- strsplit(rhs, " \\+")[[1]]
  req <- unique(c("substrate", "participant", rhs_terms))
  req <- req[req %in% names(df)]
  df <- df[complete.cases(df[, req, drop = FALSE]), ]
  if (nrow(df) < 8 || length(unique(df$substrate)) < 2) {
    cat("Not enough samples for", label, "- skipping\n"); return(NULL)
  }
  df$participant <- as.factor(df$participant)
  X <- species_wide[df$sample_key, , drop = FALSE]
  d <- vegdist(X, method = "bray")

  ctrl <- how(nperm = 999)
  form <- as.formula(paste("d ~", rhs))
  res <- adonis2(form, data = df, permutations = ctrl)
  cat("\n### PERMANOVA:", label, "###\n"); print(res)

  bd <- betadisper(d, group = df$substrate)
  pd <- permutest(bd, permutations = 999)
  cat("\n### Dispersion (betadisper) —", label, "###\n"); print(pd)
  write_csv(tibble(dataset = label, F_statistic = unname(pd$tab$F[1])),
            file.path(OUT_DIR, paste0("betadisper_", label, ".csv")))

  # pairwise contrasts (BH-corrected), restricted permutations
  arms <- sort(unique(df$substrate))
  pw <- list()
  if (length(arms) > 1) {
    for (i in seq_len(length(arms) - 1)) {
      for (j in (i + 1):length(arms)) {
        sel <- which(df$substrate %in% c(arms[i], arms[j]))
        if (length(sel) > 4 && length(unique(df$participant[sel])) > 2) {
          mm <- df[sel, ]
          dd <- vegdist(X[sel, , drop = FALSE], method = "bray")
          rr <- adonis2(dd ~ substrate, data = mm,
                        permutations = how(nperm = 499))
          pw[[length(pw) + 1]] <-
            tibble(contrast = paste(arms[i], "vs", arms[j]),
                   R2 = round(as.numeric(rr$R2[1]), 4),
                   p = as.numeric(rr$`Pr(>F)`[1]))
        }
      }
    }
  }
  pw_df <- bind_rows(pw)
  if (nrow(pw_df) > 0) {
    pw_df <- pw_df %>% mutate(p_BH = p.adjust(p, "BH"))
    cat("\n### Pairwise contrasts —", label, "###\n"); print(pw_df)
    write_csv(pw_df, file.path(OUT_DIR, paste0("pairwise_", label, ".csv")))
  }

  out <- tibble(dataset = label, term = rownames(res), R2 = round(res$R2, 4),
                F = round(res$F, 3), p = res$`Pr(>F)`)
  write_csv(out, file.path(OUT_DIR, paste0("adonis2_", label, ".csv")))

  pcoa <- cmdscale(d, k = 2, eig = TRUE, add = TRUE)
  ord <- as_tibble(pcoa$points, .name_repair = "unique") %>%
    setNames(c("Axis1", "Axis2"))
  ord$substrate   <- df$substrate
  ord$Stage       <- df$Stage
  ord$o2          <- df$o2
  ord$kefir_kind  <- df$kefir_kind
  tot <- sum(pcoa$eig[pcoa$eig > 0])
  attr(ord, "var1") <- pcoa$eig[1] / tot
  attr(ord, "var2") <- pcoa$eig[2] / tot
  list(adonis = res, ord = ord)
}
# ============================================================================
# 4. RUN PERMANOVAs
# ============================================================================
res_cit_milk <- run_permanova(
  cit %>% filter(kefir_kind == "Milk", compartment == "Liquid"),
  "citizen_milk_liquid", "substrate")
res_cit_water <- run_permanova(
  cit %>% filter(kefir_kind == "Water", compartment == "Liquid"),
  "citizen_water_liquid", "substrate")
res_lab_milk <- run_permanova(
  lab %>% filter(kefir_kind == "Milk", kind %in% c("ML", "MG")),
  "lab_milk", "substrate")
res_lab_water <- run_permanova(
  lab %>% filter(kefir_kind == "Water", kind %in% c("WL", "WG")),
  "lab_water", "substrate")

# ============================================================================
# 5. ALPHA DIVERSITY ACROSS SUBSTRATE ARMS
# ============================================================================
cat("\n--- 5. Alpha diversity ---\n")

alpha_df <- key_info %>% filter(!is.na(kind)) %>%
  mutate(shannon = NA_real_, richness = NA_real_)
for (i in seq_len(nrow(alpha_df))) {
  v <- as.numeric(species_wide[alpha_df$sample_key[i], ])
  s <- sum(v)
  if (s > 0) {
    p <- v / s
    alpha_df$shannon[i]  <- -sum(p[p > 0] * log(p[p > 0]))
    alpha_df$richness[i] <- sum(v > 0)
  }
}
write_csv(alpha_df, file.path(OUT_DIR, "alpha_diversity_per_sample.csv"))

alpha_test <- function(d, label) {
  dd <- d %>% filter(!is.na(substrate), !is.na(shannon))
  if (length(unique(dd$substrate)) < 2) return(NULL)
  kw_s <- kruskal.test(shannon ~ substrate, data = dd)
  kw_r <- kruskal.test(richness ~ substrate, data = dd)
  out <- tibble(dataset = label,
                test = c("Shannon ~ substrate (KW)", "Richness ~ substrate (KW)"),
                statistic = c(round(unname(kw_s$statistic), 2),
                              round(unname(kw_r$statistic), 2)),
                df = c(unname(kw_s$parameter), unname(kw_r$parameter)),
                p = c(kw_s$p.value, kw_r$p.value))
  cat("\n### Alpha diversity —", label, "###\n"); print(out)
  write_csv(out, file.path(OUT_DIR, paste0("alpha_", label, ".csv")))
  out
}
invisible(alpha_test(alpha_df %>% filter(source == "Citizen", kefir_kind == "Milk", compartment == "Liquid"),
                     "citizen_milk_liquid"))
invisible(alpha_test(alpha_df %>% filter(source == "Citizen", kefir_kind == "Water", compartment == "Liquid"),
                     "citizen_water_liquid"))

# ============================================================================
# 6. FUNCTIONAL / RESISTOME (guarded - not in the distributed data)
# ============================================================================
cat("\n--- 6. Functional comparison ---\n")
resist_files <- file.path(DATA_DIR, c("gene_total_resistome.csv",
                                      "mech_total_resistome.csv"))
if (all(file.exists(resist_files))) {
  gene_res <- read_csv(resist_files[1], show_col_types = FALSE)
  gkey <- names(gene_res)[1]
  gw <- gene_res %>% column_to_rownames(gkey)
  gX <- gw[rownames(gw) %in% key_info$sample_key, , drop = FALSE]
  gi <- key_info %>% filter(sample_key %in% rownames(gX), !is.na(substrate))
  gX <- gX[gi$sample_key, , drop = FALSE]
  gd <- vegdist(gX, method = "bray")
  gres <- adonis2(gd ~ substrate + Stage, data = gi,
                  permutations = how(nperm = 999, plots = Plots(strata = gi$participant)))
  print(gres)
  write_csv(tibble(term = rownames(gres), R2 = round(gres$R2, 4),
                   F = round(gres$F, 3), p = gres$`Pr(>F)`),
            file.path(OUT_DIR, "adonis2_resistome.csv"))
} else {
  msg <- paste("Resistome/SuperFocus tables (gene_total_resistome.csv,",
               "mech_total_resistome.csv) are not part of the distributed",
               "repository, so the functional substrate comparison cannot be",
               "re-run here. Pipeline for context: scripts/unix_scripts/",
               "05_short_read_functional_profiling/05_superfocus.txt")
  cat(msg, "\n")
  write_lines(msg, file.path(OUT_DIR, "functional_comparison_status.txt"))
}
# ============================================================================
# 7. VOLATILE PROFILES (published water-kefir datasets)
# ============================================================================
cat("\n--- 7. Volatile profiles ---\n")

run_vol <- function(fpath, label) {
  if (!file.exists(fpath)) return(NULL)
  raw <- read.csv(fpath, check.names = FALSE)
  samp_cols <- names(raw)[-(1:6)]
  mat <- t(as.matrix(raw[, samp_cols]))
  colnames(mat) <- raw$Compound
  mat[is.na(mat)] <- 0
  mat_log <- log10(mat + 1)

  info <- tibble(sample = samp_cols) %>%
    mutate(time_h = as.factor(str_extract(sample, "(\\d+)h")),
           vessel = str_remove(sample, "(\\d+)h"))
  info <- info %>% filter(!is.na(time_h), time_h != "")
  if (nrow(info) < 6 || n_distinct(info$time_h) < 2) return(NULL)

  tryCatch({
    d <- vegdist(mat_log, method = "bray")
    r <- adonis2(d ~ time_h + vessel, data = info, permutations = 999)
    cat("\n### Volatiles —", label, "###\n"); print(r)
    tibble(dataset = label, term = rownames(r), R2 = round(r$R2, 4),
           F = round(r$F, 3), p = r$`Pr(>F)`)
  }, error = function(e) {
    cat("volatile analysis failed for", label, ":", conditionMessage(e), "\n")
    NULL
  })
}

vol_res <- bind_rows(
  run_vol(file.path(DATA_DIR, "metabolomics", "other_datasets",
                    "global_wk_raw_metabolomics.csv"), "walsh_wk"),
  run_vol(file.path(DATA_DIR, "metabolomics", "other_datasets",
                    "wk_Breselge_et_al.csv"), "breselge_wk"))
if (!is.null(vol_res) && nrow(vol_res) > 0)
  write_csv(vol_res, file.path(OUT_DIR, "adonis2_volatiles.csv"))
# ============================================================================
# 8. FIGURE (PDF + TIFF, 180 x 130 mm)
# ============================================================================
cat("\n--- 8. Figure ---\n")

plot_pcoa <- function(ord, title) {
  ggplot(ord, aes(Axis1, Axis2, colour = substrate)) +
    geom_point(size = 1.8, alpha = 0.85) +
    stat_ellipse(linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
    labs(title = title,
         x = paste0("PCoA 1 (", round(100 * attr(ord, "var1"), 1), "%)"),
         y = paste0("PCoA 2 (", round(100 * attr(ord, "var2"), 1), "%)"),
         colour = "Substrate") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7))
}

panels <- list(); labs <- c()
if (!is.null(res_cit_milk)) {
  panels <- c(panels, list(plot_pcoa(res_cit_milk$ord, "Milk kefir, citizen cohort")))
  labs <- c(labs, "A")
}
if (!is.null(res_cit_water)) {
  panels <- c(panels, list(plot_pcoa(res_cit_water$ord, "Water kefir, citizen cohort")))
  labs <- c(labs, "B")
}

alpha_lab <- alpha_df %>% filter(!is.na(substrate), !is.na(shannon),
                                 source == "Citizen", compartment == "Liquid")
if (nrow(alpha_lab) > 0) {
  pC <- ggplot(alpha_lab, aes(substrate, shannon, fill = substrate)) +
    geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
    facet_wrap(~kefir_kind, scales = "free_x") +
    labs(title = "Alpha diversity (Shannon)", x = NULL, y = "Shannon index",
         fill = "Substrate") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
          legend.position = "none",
          strip.background = element_rect(fill = "grey95", colour = "black"),
          strip.text = element_text(size = 9, face = "bold"))
  panels <- c(panels, list(pC)); labs <- c(labs, "C")
}
lab_all <- bind_rows(
  res_lab_milk$ord %>% mutate(dataset = "Milk, lab"),
  res_lab_water$ord %>% mutate(dataset = "Water, lab"))
if (nrow(lab_all) > 0) {
  v1 <- mean(c(attr(res_lab_milk$ord, "var1"), attr(res_lab_water$ord, "var1")))
  v2 <- mean(c(attr(res_lab_milk$ord, "var2"), attr(res_lab_water$ord, "var2")))
  pD <- ggplot(lab_all, aes(Axis1, Axis2, colour = substrate, shape = o2)) +
    geom_point(size = 2) +
    facet_wrap(~dataset) +
    labs(title = "Controlled lab time-course",
         x = paste0("PCoA 1 (", round(100 * v1, 1), "%)"),
         y = paste0("PCoA 2 (", round(100 * v2, 1), "%)"),
         colour = "Substrate", shape = "O2") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))
  panels <- c(panels, list(pD)); labs <- c(labs, "D")
}

fig <- ggarrange(plotlist = panels, nrow = 2, ncol = 2, labels = labs,
                 font.label = list(size = 12), common.legend = FALSE)

pdf(file.path(FIGURES_DIR, "Figure_substrate_comparison.pdf"),
    width = 6, height = 4, useDingbats = FALSE)
print(fig)
dev.off()

tiff(filename = file.path(FIGURES_DIR, "Figure_substrate_comparison.tiff"),
     width = 1800, height = 1200, units = "px", res = 300,
     compression = "lzw", type = "cairo")
print(fig)
dev.off()

cat("\n=== DONE: results + PDF/TIFF written to", OUT_DIR, "===\n")




