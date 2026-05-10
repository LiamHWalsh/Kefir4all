# Set library paths for R packages
.libPaths("E:/STORE N GO/R/R-4.4.2/win-library/4.0")

# Load necessary libraries using pacman
pacman::p_load(readxl, readr, reshape2, dplyr, gplots, Heatplus, vegan, RColorBrewer, tidyr, gtools, 
               stringr, tidyverse, ComplexHeatmap, magick, viridis, Hotelling, stringdist, dplyr, data.table)

# Set working directory to where the Excel files are stored
setwd("Q:/H2020 Master/Citizen Science Project_Kefir4All/Results/00/")

# Read all Excel files in the directory and store them in a list
temp = list.files(pattern=".xlsx", recursive = FALSE)
myfiles = lapply(temp, read_excel)
names(myfiles) <- gsub("file_names_survey_responses_|.xlsx", "", temp)

# Load kefir metadata from a CSV file
kefir4all_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project_Kefir4All/Citizen science metadata/sample_metadata/kefir4all_sample_metadata_v2.csv")

# Recode kefir type for specific sample names
kefir4all_metadata$`kefir type`[grep("MLC|WLC", kefir4all_metadata$merge_column)] <- "Liquid control"
kefir4all_metadata$merge_column_v2 <- gsub("_host_removed_R..fastq.gz", "", kefir4all_metadata$merge_column)

# Merge metadata with observations from survey response files
kefir4all_metadata <- merge(
  kefir4all_metadata, 
  rbind(dplyr::select(myfiles[["mk"]], merge_column, Sample, observations, category_confirmed),
        dplyr::select(myfiles[["wk"]], merge_column, Sample, observations, category_confirmed)),
  by.y = "merge_column",
  by.x = "merge_column_v2",
  all.x = TRUE
)

# Recode kefir types for clarity
kefir4all_metadata <- kefir4all_metadata %>%
  mutate(`kefir type` = recode(`kefir type`,
                               "MG" = "Milk kefir grain",
                               "ML" = "Milk kefir liquid",
                               "WG" = "Water kefir grain",
                               "WL" = "Water kefir liquid"))

# Rename timepoints from T1-T6 to weeks
kefir4all_metadata$Stage <- gsub("T1", "wk01", kefir4all_metadata$Stage)
kefir4all_metadata$Stage <- gsub("T2", "wk05", kefir4all_metadata$Stage)
kefir4all_metadata$Stage <- gsub("T3", "wk09", kefir4all_metadata$Stage)
kefir4all_metadata$Stage <- gsub("T4", "wk13", kefir4all_metadata$Stage)
kefir4all_metadata$Stage <- gsub("T5", "wk17", kefir4all_metadata$Stage)
kefir4all_metadata$Stage <- gsub("T6", "wk21", kefir4all_metadata$Stage)

# Create a new column describing the sample
kefir4all_metadata$description = paste(kefir4all_metadata$Stage, kefir4all_metadata$`kefir type`, sep="- ")

# Generate filenames for forward and reverse sequencing reads
kefir4all_metadata$forward_file_name = paste(kefir4all_metadata$merge_column_v2, "_host_removed_R1.fastq.gz", sep="")
kefir4all_metadata$reverse_file_name = paste(kefir4all_metadata$merge_column_v2, "_host_removed_R2.fastq.gz", sep="")

# Remove duplicated rows based on merge_column_v2
kefir4all_metadata = kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$merge_column_v2))),]

# Define column names for output files
column_names <- c("tax_id", "scientific_name", "sample_alias", "sample_title", "sample_description", 
                  "metagenomic source", "project name", "collection date", 
                  "geographic location (country and/or sea)")

# Load MD5 checksum files for integrity check
md5sums_R1 <- read_table("Q:/H2020 Master/Citizen Science Project_Kefir4All/ena submission/md5sums_R1.txt", col_names = FALSE)
md5sums_R2 <- read_table("Q:/H2020 Master/Citizen Science Project_Kefir4All/ena submission/md5sums_R2.txt", col_names = FALSE)

# Ensure forward and reverse file names match
identical(
  gsub("_host_removed_R1.fastq.gz", "", md5sums_R1$X2),
  gsub("_host_removed_R2.fastq.gz", "", md5sums_R2$X2)
)

# Merge MD5 checksum data
md5sums_total = cbind(md5sums_R1, md5sums_R2)
colnames(md5sums_total) = c("forward_file_md5", "forward_file_name", "reverse_file_md5", "reverse_file_name")
md5sums_total$merge_column = gsub("_host_removed_R1.fastq.gz", "", md5sums_total$forward_file_name)

# Merge MD5 checksum information with metadata
kefir4all_metadata <- merge(kefir4all_metadata, md5sums_total, by.y="merge_column", by.x="merge_column_v2", all.x=TRUE)

# Create submission file for ENA metadata
output_dir = "Q:/H2020 Master/Citizen Science Project_Kefir4All/ena submission"
data_file_for_use_in_ena_submission <- data.frame(
  NCBI_Taxonomy = as.character("1326787"),
  scientific_name = as.character("fermentation metagenome"),
  sample_alias = as.character(kefir4all_metadata$merge_column_v2),
  sample_title = as.character("Kefir metagenomic reads"),
  sample_description = kefir4all_metadata$description,
  `metagenomic source` = as.character("Kefir metagenome"),
  `project name` = "PRJEB77409",
  `collection date` = as.character("not provided"),
  `geographic location (country and/or sea)` = as.character("Ireland")
)
colnames(data_file_for_use_in_ena_submission) = column_names
fwrite(data_file_for_use_in_ena_submission, file.path(output_dir, "data_file_for_use_in_ena_submission.tsv"), sep = "\t")

# Create submission file for sequencing reads
data_file_for_use_in_ena_submission_for_files <- data.frame(
  sample = as.character(kefir4all_metadata$merge_column_v2),
  study = "PRJEB77409",
  instrument_model = "Illumina NovaSeq 6000",
  library_name = as.character(kefir4all_metadata$merge_column_v2),
  library_source = "METAGENOMIC",
  library_selection = "PCR",
  library_strategy = "WGS",
  library_layout = "PAIRED",
  forward_file_name = as.character(kefir4all_metadata$forward_file_name.x),
  forward_file_md5 = as.character(kefir4all_metadata$forward_file_md5),
  reverse_file_name = as.character(kefir4all_metadata$reverse_file_name.x),
  reverse_file_md5 = as.character(kefir4all_metadata$reverse_file_md5)
)
fwrite(data_file_for_use_in_ena_submission_for_files, file.path(output_dir, "data_file_for_use_in_ena_submission_for_files.tsv"), sep = "\t")
