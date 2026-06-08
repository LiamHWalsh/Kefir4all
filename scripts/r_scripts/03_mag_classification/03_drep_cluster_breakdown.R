# 03_drep_cluster_breakdown.R ---------------------------------------------
# Pipeline step 3 -- MAG classification
# Results sections: 3.5 (strain sub-clusters), Supplementary Result S6
#
# Loads dRep secondary-cluster assignments (Cdb.csv) and MAG metadata,
# filters to prevalent species, and produces a bar chart of MAG counts
# per species coloured by secondary cluster.
#
# Inputs (all in data/):
#   Cdb.csv                              dRep cluster assignments
#   mag_metadata_v3.csv.gz               GTDB taxonomy + CheckM quality
#   kefir4all_sample_metadata_v2.csv
#   global_milk_kefir_metadata_v1.csv
#   global_water_kefir_metadata_v1.csv
#   milk_taxonomic_profile_prevalence.csv
#   water_taxonomic_profile_prevalence.csv
#
# Outputs:
#   figures/Figure_S6_drep_cluster_breakdown.png

if (!requireNamespace("here",   quietly = TRUE)) install.packages("here")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(here)
pacman::p_load(readr, dplyr, tidyr, ggplot2, stringr, tidytext)

DATA_DIR    <- here::here("data")
FIGURES_DIR <- here::here("figures")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# Private metadata -- not committed to repo
CS_METADATA_PRIVATE <- file.path(here::here("data", "private"),
                                  "Citizen Scientist metadata_v8.csv")
if (!file.exists(CS_METADATA_PRIVATE)) {
  stop(
    "Private metadata file not found: ", CS_METADATA_PRIVATE,
    "\nThis file contains citizen-scientist identifiers and is not committed ",
    "to the repository. Copy it to data/private/ to run this script."
  )
}

# -------------------------------------------------------------------------
# Import metadata
# -------------------------------------------------------------------------
global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv"),
                                show_col_types = FALSE)
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv"),
                                show_col_types = FALSE)
global_mk_metadata$Stage <- NA
global_wk_metadata$Stage <- NA

kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"),
                                show_col_types = FALSE)
kefir4all_metadata$merge_column <- gsub("_host_removed_R..fastq.gz", "",
                                         kefir4all_metadata$merge_column)
kefir4all_metadata <- kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$merge_column))), ]

total_metadata <- rbind(
  dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`, Stage, Sample),
  dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`, Stage, Sample),
  dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`, Stage, Sample)
)
total_metadata$category <- NA
total_metadata$category[total_metadata$`kefir type` %in% c("WL", "WG")] <- "Water.kefir"
total_metadata$category[total_metadata$`kefir type` %in% c("ML", "MG")] <- "Milk.kefir"

# -------------------------------------------------------------------------
# Load dRep cluster assignments and MAG metadata
# -------------------------------------------------------------------------
Cdb <- read_csv(file.path(DATA_DIR, "Cdb.csv"), show_col_types = FALSE)

mag_metadata <- read_csv(file.path(DATA_DIR, "mag_metadata_v3.csv.gz"),
                          show_col_types = FALSE)

mag_metadata$bin <- gsub(".*_bins_", "", mag_metadata$user_genome)
mag_metadata$bin <- gsub(".orig|.permissive|.strict", "", mag_metadata$bin)
mag_metadata$bin <- gsub(".*_", "", mag_metadata$bin)
mag_metadata$bin_id <- paste(mag_metadata$base_name, "_", mag_metadata$bin, sep = "")

mag_metadata <- mag_metadata[-c(which(duplicated(mag_metadata$user_genome))), ]
mag_metadata$classification_full <- mag_metadata$classification
mag_metadata$classification <- gsub(".*;s__", "", mag_metadata$classification)

mag_metadata <- dplyr::select(
  mag_metadata,
  "user_genome", "base_name", "bin_id", "classification", "classification_full",
  "Sample.y", "kefir type", "Stage", "data_source", "type.x",
  "Timepoint", "category", "Completeness", "Contamination", "linker"
)
colnames(mag_metadata)[colnames(mag_metadata) == "Sample.y"] <- "Sample"
colnames(mag_metadata)[colnames(mag_metadata) == "type.x"]   <- "type"

# -------------------------------------------------------------------------
# Merge dRep clusters with MAG metadata
# -------------------------------------------------------------------------
Cdb$genome <- gsub(".fa", "", Cdb$genome)
Cdb <- merge(mag_metadata, Cdb, by.x = "user_genome", by.y = "genome", all.y = TRUE)

# Fill blank species names with genus + ".species" placeholder
Cdb$classification[Cdb$classification == ""] <- paste(
  gsub(".*g__|;s__", "", Cdb$classification_full[Cdb$classification == ""]),
  ".species", sep = ""
)

# -------------------------------------------------------------------------
# Load prevalence tables and filter to prevalent species only
# -------------------------------------------------------------------------
milk_prev  <- read_csv(file.path(DATA_DIR, "milk_taxonomic_profile_prevalence.csv"),
                        show_col_types = FALSE)
water_prev <- read_csv(file.path(DATA_DIR, "water_taxonomic_profile_prevalence.csv"),
                        show_col_types = FALSE)
total_prevalence <- rbind(milk_prev, water_prev)

# Harmonise legacy species labels
sp_col <- names(total_prevalence)[5]
recode_sp <- c(
  "Lactobacillus_ghanensis"         = "Liquorilactobacillus ghanensis",
  "Lactococcus_lactis subcluster 1" = "Lactococcus lactis",
  "Lactococcus_lactis subcluster 2" = "Lactococcus cremoris",
  "Pseudomonas_fragi_subspecies 1"  = "Pseudomonas fragi",
  "Zymomonas mobilis_subcluster 1"  = "Zymomonas mobilis"
)
total_prevalence[[sp_col]] <- ifelse(
  total_prevalence[[sp_col]] %in% names(recode_sp),
  recode_sp[total_prevalence[[sp_col]]],
  total_prevalence[[sp_col]]
)
total_prevalence[[sp_col]] <- gsub("_", " ", total_prevalence[[sp_col]])

cdb_prevalent <- rbind(
  Cdb[Cdb$classification %in%
        total_prevalence[[sp_col]][total_prevalence$kefir_type == "milk"] &
        Cdb$category == "Milk.kefir", ],
  Cdb[Cdb$classification %in%
        total_prevalence[[sp_col]][total_prevalence$kefir_type == "water"] &
        Cdb$category == "Water.kefir", ]
)

# -------------------------------------------------------------------------
# Summarise and plot
# -------------------------------------------------------------------------
mag_prevalent_breakdown <- as.data.frame(
  xtabs(~ classification + data_source + secondary_cluster + category,
        cdb_prevalent)
)
mag_prevalent_breakdown$db_merge <- paste(
  mag_prevalent_breakdown$data_source,
  mag_prevalent_breakdown$category, sep = "_"
)

p <- mag_prevalent_breakdown %>%
  mutate(detection_category = gsub(
    "This study_Milk.kefir",        "Milk kefir - Kefir4all",
    gsub("This study_Water.kefir",  "Water kefir - Kefir4all",
    gsub("Walsh et al_Milk.kefir",  "Milk kefir - Walsh et al 2023",
    gsub("Walsh et al 2023_Milk.kefir", "Milk kefir - Walsh et al 2023",
    gsub("Samuel et al_Water.kefir","Water kefir - Mortensen et al 2023",
         db_merge)))))) %>%
  mutate(new_type = gsub(" -.*", "", detection_category),
         secondary_cluster_v2 = gsub(".*_", "", secondary_cluster)) %>%
  filter(Freq != 0) %>%
  ggplot(aes(x = reorder(classification, -Freq), y = Freq,
             fill = secondary_cluster_v2)) +
  geom_col() +
  facet_wrap(~ new_type, scales = "free") +
  labs(x = "Species", y = "Number of MAGs", fill = "Secondary cluster") +
  theme_bw() +
  theme(
    legend.position  = "top",
    axis.ticks       = element_blank(),
    axis.text.x      = element_text(size = 15),
    axis.title       = element_text(size = 20),
    axis.text.y      = element_text(size = 12.5, face = "italic"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background  = element_blank(),
    legend.text      = element_text(size = 20),
    legend.key.size  = unit(1.5, "cm"),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text       = element_text(size = 15.5)
  ) +
  coord_flip()

ggsave(file.path(FIGURES_DIR, "Figure_S6_drep_cluster_breakdown.png"),
       plot = p, width = 16, height = 12, dpi = 300)

# -------------------------------------------------------------------------
# Summary statistics (MAG quality, cluster counts)
# -------------------------------------------------------------------------
message("Primary cluster count:  ", nrow(as.data.frame(table(Cdb$primary_cluster))))
message("Secondary cluster count: ", nrow(as.data.frame(table(Cdb$secondary_cluster))))
message("Completeness range: ",
        paste(round(range(Cdb$Completeness, na.rm = TRUE), 1), collapse = " - "))
message("Contamination range: ",
        paste(round(range(Cdb$Contamination, na.rm = TRUE), 1), collapse = " - "))
