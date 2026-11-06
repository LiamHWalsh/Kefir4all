# Anchor to the project root before here() loads. RStudio Server resets the session cwd
# to $HOME after sourcing .Rprofile, and here() computes its root AT PACKAGE LOAD
# TIME, so force the project root before any here load. This block is inert on non-portable
for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {



  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}


# ==============================================================================
# PATH CONFIGURATION â€” portable relative paths via the here package
# For files listed in data/private/README.md, provide your own copies.
# NOTE: no install.packages() calls are needed here -- the Docker image has all packages pre-installed.
# ==============================================================================


# The anchor block at lines 1-10 ensures we are in the project root already.
# Just confirm here() will resolve against the correct root.
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style()
here::i_am("scripts/r_scripts/04_taxonomic_profiling/04_taxonomic_profiling.R")
DATA_DIR    <- here::here("data")                         # public data files
PRIVATE_DIR <- here::here("data", "private")              # private/sensitive files
OUT_DIR <- here::here("output_reviewed", "04_taxonomic_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "Figure_2")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# Private metadata â€” not distributed with repo (participant identifiers)
CS_METADATA_PRIVATE <- file.path(PRIVATE_DIR, "Citizen Scientist metadata_v8.csv")
if (!file.exists(CS_METADATA_PRIVATE)) {
  stop(
    "This script requires the private citizen-scientist metadata file which is not\n",
    "distributed with the public repository (contains participant identifiers).\n",
    "Missing: ", CS_METADATA_PRIVATE, "\n",
    "See data/private/README.md for details."
  )
}
# ==============================================================================
########################################################################################################################
#Research objectives in this script
#Describe the compositional profile of the CS dataset only
# Identify prevalent species
#identify dominating species
# track persistence of microbes in the grain that may be environmentally sourced
# Identify pathogenic microbes
########################################################################################################################


########################################################################################################################
#Libraries used 
########################################################################################################################


library(tibble)
library(tidyr)
library(scales)

library(readxl)
library(readr)
library(dplyr)
library(ggplot2)
library(grid)
library(ggpubr)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ggh4x)


########################################################################################################################
#Import metacache data
########################################################################################################################
metacache <- read_csv(file.path(DATA_DIR, "04_metacache_total_species_profile_v2.csv"))
metacache_strain <-  read_csv(file.path(DATA_DIR, "04_metacache_total_strain_profile_v2.csv"))

metacache$sample_id <- gsub("-","_",metacache$sample_id)

metacache$sample_id <- gsub(".*/","",metacache$sample_id)


########################################################################################################################
#Import metadata
########################################################################################################################

global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv"))
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv"))
Citizen_Scientist_metadata_v8 <- read_csv(CS_METADATA_PRIVATE)

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )

########################################################################################################################
#Import and modify survey_responses 
########################################################################################################################
# Survey response xlsx files â€” already in data/ (non-identifying aggregate responses)
myfiles <- list(
  mk = read_excel(file.path(DATA_DIR, "file_names_survey_responses_mk.xlsx")),
  wk = read_excel(file.path(DATA_DIR, "file_names_survey_responses_wk.xlsx"))
)




kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"))

#kefir4all_metadata$`kefir type`[grep("MLC|WLC",kefir4all_metadata$merge_column)] <- "Liquid control"
kefir4all_metadata$merge_column <-  gsub("_host_removed_R..fastq.gz","",kefir4all_metadata$merge_column)
kefir4all_metadata <- kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$merge_column))),]

print(paste("The citizen science project is composed of ",nrow(kefir4all_metadata))) 

print(paste("number of sample according to types in this study"))

print(table(kefir4all_metadata$`kefir type`))
########################################################################################################################
# Merge metadata into one
########################################################################################################################
total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`))


########################################################################################################################
#Get data from the global_MK study
########################################################################################################################

#metacache <- metacache[which(gsub("_S.*","", metacache$sample_id) %in% global_mk_metadata$Sample),]

########################################################################################################################
#Modufy dataset into metacache citizen scientist and baseline used
########################################################################################################################


metacache [is.na(metacache )] <- 0
kefir4all_metadata <- 
kefir4all_metadata[-c(which(
kefir4all_metadata$`kefir type`=="Medium control")),]

metacache <- metacache[which(metacache$sample_id %in% 
kefir4all_metadata$merge_column), ]

#n1_species <-which(maxab_species < 0.1)


metacache <- 
metacache |> column_to_rownames("sample_id")



print(paste("Number of samples in metacache is",length(rownames(metacache))))

# samples  not accounted for by metacache- base_line typically did not have a compostional profile 



kefir4all_metadata$merge_column[-c(which(kefir4all_metadata$merge_column %in% rownames(metacache)))]


###################
#Prevalence info
specie <- c()
sum_specie <- c()
species_prevalence <- c()


species_profile <- c()
species_profile[["water"]] <- metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL"))]),]#& 
                                                                                                               #kefir4all_metadata$Stage!="T0")]),]

species_profile[["milk"]] <- metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML"))]),] 
          
                                                                                                                                                                                                                

 #species_profile[["water"]] <- metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL")& 
                                                                                                              #kefir4all_metadata$Stage!="T0")]),]
                                                                                                                                                                                                                          
                                                                                                                                                                                                                           
#species_profile[["water_baseline"]] <- 
#metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                                                                                                                                                                       
        
#species_profile[["milk_baseline"]] <- 
  #metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                   #kefir4all_metadata$Stage=="T0")]),]                             
  
  


#taxize_name <- classification(tax_info$name, db = "ncbi")

tax_phylo<- c()
pa <- c()
metadata <- c()
type <-  names(species_profile)[4]
maxab_species <- c()
total_compositional_data <- data.frame(sample=character(), clade_name=character(), relative_abundance=numeric(), type=character(), stringsAsFactors = FALSE)
type <- "water"

for (type in names(species_profile)){
  
  maxab_species <- apply(species_profile[[type]],2, max, na.rm=TRUE)
  
  n1_species <-names(which(maxab_species > .1))
  
  species_profile[[type]] <- species_profile[[type]][,which(colnames(species_profile[[type]]) %in%  n1_species )]
  
  species_prevalence[[type]] <- data.frame(specie=as.character(),
                                       sum_specie=as.numeric())
  
  #i <- "Lentilactobacillus hilgardii" 
for (i in colnames(species_profile[[type]])){
  specie <- i
  
  sum_specie <-  length(which(species_profile[[type]][,i] >.1))
  
  data.prevalence <- cbind(specie,sum_specie)
  species_prevalence[[type]] <- merge(data.prevalence,species_prevalence[[type]], all=TRUE)
  
}
  species_prevalence[[type]] <- species_prevalence[[type]][which(as.numeric(species_prevalence[[type]]$sum_specie)> nrow(species_profile[[type]])*.1),]
  species_prevalence[[type]]$kefir_type<- type
  species_profile[[type]] <- species_profile[[type]][,which(colnames(species_profile[[type]]) %in%    species_prevalence[[type]]$specie)]
  
  
  
species_profile[[type]] <-
  species_profile[[type]] |> rownames_to_column("sample") |> pivot_longer(!sample, names_to = "clade_name",values_to = "relative_abundance")

species_profile[[type]]$type <- type
total_compositional_data <- rbind(species_profile[[type]],total_compositional_data )

print(paste(type, "has", length(levels(as.factor( species_profile[[type]]$sample))),"number of sample"))



print(paste(type, "has", length(levels(as.factor( species_profile[[type]]$clade_name))),"prevalent microbes"))
  #.02
  
}

#write.csv(species_prevalence[["milk"]],"Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/milk_metacache_prevalence.csv" )
#write.csv(species_prevalence[["water"]],"Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/water_metacache_prevalence.csv" )



#ggheatmap(
#metacache[,which(colnames(metacache) %in% levels(as.factor(total_compositional_data$clade_name)))]) 


total_compositional_data <- metacache[,which(colnames(metacache) %in% levels(as.factor(total_compositional_data$clade_name)))] |>
  rownames_to_column("sample") |> pivot_longer(!sample  , names_to  ="clade_name" , values_to ="relative_abundance")
                                     



.n_kefir4all_metadata_before_survey_merge <- nrow(kefir4all_metadata)
kefir4all_metadata <-
  merge(
    kefir4all_metadata,
    rbind(dplyr::select(myfiles[["mk"]], merge_column,Sample, observations, category_confirmed),
          dplyr::select(myfiles[["wk"]], merge_column,Sample, observations, category_confirmed)),
    by="merge_column",
    all.x=TRUE)
# WHY: all.x = TRUE should preserve exactly one output row per kefir4all_metadata row;
# a row-count increase here would mean the mk/wk survey tables have a duplicate
# merge_column, silently duplicating every downstream analysis of this sample.
stopifnot("survey-response merge changed kefir4all_metadata's row count -- check mk/wk survey files for duplicate merge_column values" =
            nrow(kefir4all_metadata) == .n_kefir4all_metadata_before_survey_merge)




stopifnot("kefir4all_metadata has duplicate merge_column values -- the next merge (inner join) would silently duplicate rows" =
            !anyDuplicated(kefir4all_metadata$merge_column))
total_compositional_data <-
  merge(total_compositional_data,kefir4all_metadata,by.x="sample",by.y="merge_column")



total_compositional_data$type <- NA
#total_compositional_data$type <- gsub("milk","Milk kefir",total_compositional_data$type)
#total_compositional_data$type <- gsub("water","Water kefir",total_compositional_data$type)


#table(dplyr::select(total_compositional_data,clade_name,type))

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                                                 kefir4all_metadata$Stage!="T0")])] <- "Water kefir"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                                                 kefir4all_metadata$Stage!="T0")])] <- "Milk kefir"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                                     kefir4all_metadata$Stage=="T0")])] <- "Water\n kefir\n T0"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                                                 kefir4all_metadata$Stage=="T0")])] <- "Milk\n kefir\n T0"


# remove extraction controls
#total_compositional_data <- total_compositional_data[-c(which(is.na(
    #total_compositional_data$type))),]

#total_compositional_data$type <- gsub("Milk kefir_baseline","Milk kefir - T0",total_compositional_data$type)
#total_compositional_data$type <- gsub("Water kefir_baseline","Water kefir - T0",total_compositional_data$type)


#print(paste(table(dplyr::select(total_compositional_data,`kefir type`))))

###########################################################################################################################################################

# Heatmap of all prevalent microbes

###########################################################################################################################################################





total_compositional_data$kefirdataset <- NA
total_compositional_data$kefirdataset[which(total_compositional_data$`kefir type` %in% c("WL","WG"))] <- "Water kefir"

total_compositional_data$kefirdataset[which(total_compositional_data$`kefir type` %in% c("ML","MG"))] <- "Milk kefir"




#total_compositional_data$`kefir type`[which(total_compositional_data$`kefir type` %in% c("MG","WG"))] <- "Grain"

#total_compositional_data$`kefir type`[which(total_compositional_data$`kefir type` %in% c("ML","WL"))] <- "Liquid"



#which(is.na(total_compositional_data$kefirdataset))
if (!requireNamespace("ggh4x", quietly=TRUE)) install.packages("ggh4x")
library(ggh4x)

# Map Stage codes to readable timepoint labels for x-axis strips
total_compositional_data$timepoint <- recode(total_compositional_data$Stage,
  "T0"="T0", "T1"="W1", "T2"="W5", "T3"="W9", "T4"="W13", "T5"="W17", "T6"="W21")

p_heatmap <- c()
i <- "Milk kefir"


total_compositional_data <- total_compositional_data[-c(which(total_compositional_data$Sample.x %in% c("EC",   "Medium control") )),]



names(
species_prevalence)[1] <- "Water kefir"


names(
  species_prevalence)[2] <-"Milk kefir"


# ---- Shared relative-abundance colour scale for both Figure 2 heatmap panels ----
# A simple five-stop gradient (dark blue to dark red) with fixed breaks at
# 0/10/30/70/100 and a single-line "Relative abundance" legend title, shared
# by both panels so they read on the same scale.
shared_ra_scale <- scale_fill_gradientn(
  colours = c("Dark blue","yellow","yellow2","orange","red","darkred"),
  breaks  = c(0, 10, 30, 70, 100),
  limits  = c(0, 100),
  guide   = "colorbar",
  name    = "Relative abundance",
  labels  = c(0, 10, 30, 70, 100)
)


fig2_species_levels <- list()

for (i in levels(as.factor((total_compositional_data$kefirdataset)))){



  # ------------------------------------------------------------------------------

  # Panel A (milk) / B (water): ONE single continuous heatmap (no grain/liquid

  # facet split), with the columns DIVIDED INTO THE WEEK BLOCKS OF THE STUDY

  # (wk01, wk05, wk09, wk13, wk17, wk21, T0). Each week block gets a white

  # header cell on top with the week label centered inside it, separated from

  # the next block by a bold vertical divider, and the grain/liquid composition

  # of every column is shown by the colour bar underneath the heatmap

  # (Grain = salmon, Liquid = teal).

  # NOTE: the heatmap axes are NUMERIC (x_pos, y_pos) rather than discrete
  # factors, and columns are drawn with geom_rect-style numeric positioning
  # rather than geom_tile(width=n). This renders every week block at its true
  # width, so even a narrow block (e.g. wk21, with few samples) stays wide
  # enough for its label to stay readable.

  # ------------------------------------------------------------------------------

  fig2_hm_df <- total_compositional_data[which(total_compositional_data$kefirdataset==i &

                                               total_compositional_data$clade_name  %in% species_prevalence[[i]]$specie

                                               ),]|>

    mutate(`kefir type`=gsub("ML|WL","Liquid",

                             gsub("MG|WG", "Grain",`kefir type` ))) |>

    # Stage codes -> study weeks: T1=wk01, T2=wk05, T3=wk09, T4=wk13, T5=wk17, T6=wk21

    mutate(week = factor(c("T0"="T0", "T1"="wk01", "T2"="wk05", "T3"="wk09",

                           "T4"="wk13", "T5"="wk17", "T6"="wk21")[as.character(Stage)],

                         levels = c("T0","wk01","wk05","wk09","wk13","wk17","wk21")))



  # Panel title carried in `type` so facet_wrap renders a single-strip
  # boxed title (white fill, black border), consistent with panels A/B.
  fig2_hm_df$type <- i

  n_dropped_week <- sum(is.na(fig2_hm_df$week))

  if (n_dropped_week > 0) {

    fig2_hm_df <- fig2_hm_df[-c(which(is.na(fig2_hm_df$week))),]

  }



  # WHY: alphabetical by species name, matching the "explicit taxonomic axis ordering"
  # promised to reviewers in the rebuttal (REBUTTAL_total.docx, responses to Reviewer 1
  # Comment 8 and Reviewer 2 Comment 10) -- the prevalence pipeline's delivery order was
  # arbitrary (inherited from the source CSV's column order), not a deliberate ordering.
  # HOW: species_levels[1] renders at the BOTTOM of this y-axis (header " " is appended
  # last, placing it at the top) -- descending order puts the alphabetically-first species
  # at the top (just under the header) and reads A-to-Z top-to-bottom.
  species_levels <- sort(unique(as.character(fig2_hm_df$clade_name)), decreasing = TRUE)

  fig2_species_levels[[i]] <- species_levels

  fig2_hm_df <- fig2_hm_df |>

    mutate(clade_name = factor(clade_name, levels = c(species_levels, " ")))

  # y_pos: numeric position of each row on the y-axis. The header level " "

  # (appended as the last factor level) gets the highest y_pos, placing it at
  # the top of the heatmap. as.numeric(factor) returns the level index.
  fig2_hm_df <- fig2_hm_df |>
    mutate(y_pos = as.numeric(clade_name))

  # header_y: the y-position of the week-label header row (top of heatmap).
  header_y <- length(species_levels) + 1

  # One row per sample, ordered: week block -> grain before liquid -> sample id

  sample_order_df <- fig2_hm_df |>

    dplyr::distinct(sample, `kefir type`, week) |>

    dplyr::arrange(week, `kefir type`, sample) |>

    mutate(sample_combined = paste(`kefir type`, sample, sep = "_"),
           x_pos = dplyr::row_number())



  fig2_hm_df <- fig2_hm_df |>

    mutate(sample_combined = factor(paste(`kefir type`, sample, sep = "_"),
                                    levels = sample_order_df$sample_combined)) |>
    dplyr::left_join(sample_order_df |> dplyr::select(sample, x_pos), by = "sample")



  # Week-block geometry on the NUMERIC x scale: block width (n) and the numeric

  # of each block (for the header cell/label) and the block end (for dividers)

  block_info <- sample_order_df |>

    dplyr::mutate(row_idx = dplyr::row_number()) |>

    dplyr::group_by(week) |>

    dplyr::summarise(n = dplyr::n(),

                     x_mid_pos = mean(row_idx),

                     x_min_pos = min(row_idx) - 0.5,
                     x_max_pos = max(row_idx) + 0.5,

                     .groups = "drop")



  print(paste("Figure 2 (", i, ") columns divided by week: ",


              sep = ""))




  p_heatmap[[i]] <- ggplot(
    fig2_hm_df,
    aes(x = x_pos, y = y_pos)
  ) +
    
    # ------------------------------------------------------------
  # HEATMAP
  # ------------------------------------------------------------
  geom_tile(
    aes(fill = relative_abundance),
    linewidth = 0.2
  ) +
    
    # ------------------------------------------------------------
  # KEFIR TYPE BAR
  # ------------------------------------------------------------
  geom_point(
    aes(
      x = x_pos,
      y = 0,
      col = `kefir type`
    ),
    size = 4.6,
    shape = 15,
    inherit.aes = FALSE
  ) +
    
    # ------------------------------------------------------------
  # FACET BY WEEK
  # ------------------------------------------------------------
  facet_wrap(
    ~ week,
    scales = "free_x",
    nrow = 1
  ) +
    
    # ------------------------------------------------------------
  # LABELS / SCALES
  # ------------------------------------------------------------
  labs(
    x = "",
    y = "",
    title = ""
  ) +
    
    shared_ra_scale +
    
    scale_y_continuous(
      breaks = seq_along(species_levels),
      labels = species_levels,
      expand = expansion(mult = c(0, 0))
    ) +
    
    # ------------------------------------------------------------
  # THEME
  # ------------------------------------------------------------
  theme_bw() +
    
    theme(
      
      legend.position = "top",

      # White facet-strip box, not theme_bw()'s grey85 default.
      strip.background = element_rect(fill = "white", colour = "black"),

      # Keep the week labels, since these now identify each panel
      strip.text = element_text(
        size = 22,
        face = "bold"
      ),
      
      legend.key.size = unit(3, "cm"),
      legend.key.height = unit(3, "cm"),
      legend.key.width = unit(3, "cm"),
      
      legend.title = element_text(
        size = 17.5,
        face = "bold",
        margin = margin(r = 4)
      ),
      
      legend.text = element_text(
        size = 17.5
      ),
      
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      
      axis.text.y = element_text(
        size = 13.5,
        face = "italic",
        hjust = 1,
        margin = margin(r = 1)
      ),
      
      axis.title = element_text(
        face = "bold",
        size = 10.5
      ),
      
      panel.grid = element_blank(),
      
      # Make the facet panels visually distinct
      panel.spacing.x = unit(0.15, "cm"),
      
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    ) +
    
    guides(
      color = guide_legend(
        override.aes = list(size = 16),
        title = "kefir type"
      )
    )
}







 #graphics.off()
 

 
 ###########################################################################################################################################################
 
 
 # Identify unique prevalent microbes in either milk or water kefir
 
 ###########################################################################################################################################################
 
  # 
   #species <-  "Leclercia adecarboxylata"  
   #species <- "Acetobacter aceti"
   
   #Look for microbes that are not present during the project but are present at T0

   total_compositional_data$type <- gsub("\n", "", total_compositional_data$type)
   total_compositional_data$type <- gsub(" ", ".", total_compositional_data$type)
   
   compositional_differences <- data.frame(species=as.character( levels(as.factor(total_compositional_data$clade_name))),
                                              `Milk kefir`=as.character(NA),
                                              `Milk kefir T0`=as.character(NA),
                                              `Water kefir`=as.character(NA),
                                              `Water kefir T0`=as.character(NA)
                                              )
   
print(paste("This dataset contains a total of", length(levels(as.factor(total_compositional_data$clade_name))), "prevalent microbes"))
   
 
  species <-  "Lactobacillus kefiranofaciens"
  type_2 <- "Water.kefir"
  for (species in compositional_differences$species){

    for (type_2 in levels(as.factor(total_compositional_data$type))){
      compositional_differences[which(   compositional_differences$species==species),which(colnames(   compositional_differences )==type_2)] <-
      nrow(subset(total_compositional_data, subset=clade_name==species &
               type==type_2 &
               relative_abundance >0.1))
    }
  }
   


# identify prevalent microbes in milk kefir but not water
t <- c()
  t[["milk"]] <- 
subset( compositional_differences, subset=
          compositional_differences$Milk.kefir.T0 >=0 &
          compositional_differences$Water.kefir.T0 ==0 &
          compositional_differences$Milk.kefir>=0&
          compositional_differences$Water.kefir==0
)$species




t[["water"]] <- subset( compositional_differences, subset=
          compositional_differences$Milk.kefir.T0 ==0 &
          compositional_differences$Water.kefir.T0 >=0 &
          compositional_differences$Milk.kefir==0&
          compositional_differences$Water.kefir>=0
)$species




#           compositional_differences$Milk.kefir.T0 >=0 &
#           compositional_differences$Water.kefir.T0 ==0 &
#           compositional_differences$Milk.kefir>=0&
#           compositional_differences$Water.kefir==0
# )$species

###########################################################################################################################################################


# Of these prevlent microbes, what is the range of RA they occur at and whats their number of samples tht they dominate in 

###########################################################################################################################################################

sample_name <- "TG_AP_T1_WL_0AE_S350"


dominanting_species=data.frame(sample=as.character(),
                               clade_name=as.character(),
                               relative_abundance=as.numeric(),
                               Sample=as.character(),
                               `kefir type`=as.character(),
                               Stage =as.character(),
                               data_source =as.character(),
                               type=as.character()
  
  
)               

# Get the dominant species per sample
sample_name <- "TG_ID012_MG_T3_S189"




for (sample_name in levels(as.factor(total_compositional_data$sample))){
  if(sample_name %in% c("TG_ID012_MG_T3_S189", "TG_ID115_MG_T4_S96")){
    next
  }
  dominanting_species <- rbind(dominanting_species,
                               
  total_compositional_data |> 
    filter(sample== sample_name) |> 
    filter(relative_abundance==max(relative_abundance))
  
  )
}


dominanting_species[which(duplicated(dominanting_species$sample)),]

# breakdown of dominate species per kefir type



dominate_breakdown <- as.data.frame(
table(dplyr::select(  dominanting_species, clade_name,`kefir type`))) 

dominate_breakdown$Freq <- as.numeric(dominate_breakdown$Freq)

dominate_breakdown <- 
dominate_breakdown|> pivot_wider(names_from="kefir.type",values_from = "Freq")



# range of values per dominate type

dominate_breakdown$min <- NA
dominate_breakdown$max <- NA


for (sample_name in levels(as.factor(dominate_breakdown$clade_name))){
  
  

  
  dominate_breakdown$min[which(dominate_breakdown$clade_name==sample_name)] <-  getElement(   dominanting_species |> 
                                           filter(clade_name== sample_name) |> 
                                           filter(relative_abundance==min(relative_abundance)),"relative_abundance")
  
  
  dominate_breakdown$max [which(dominate_breakdown$clade_name==sample_name)]<-  getElement(   dominanting_species |> 
                                           filter(clade_name== sample_name) |> 
                                           filter(relative_abundance==max(relative_abundance)),"relative_abundance")
  
  
  
}


# get the percent of samples with these microbes 
samples_mk<- 
length(rownames(metacache)[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG"))])])

samples_ml<- 
  length(rownames(metacache)[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("ML"))])])

samples_wg<- 
  length(rownames(metacache)[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG"))])])

samples_wl<- 
  length(rownames(metacache)[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WL"))])])


dominate_breakdown$MG_percent <- (dominate_breakdown$MG / samples_mk)*100

dominate_breakdown$ML_percent <- (dominate_breakdown$ML / samples_ml)*100
dominate_breakdown$WG_percent <- (dominate_breakdown$WG / samples_wg)*100
dominate_breakdown$WL_percent <- (dominate_breakdown$WL / samples_wl)*100

length(rownames(metacache)[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL"))])])



sum(dominate_breakdown$MG_percent)

# = ggplot(pcoa, aes(x=PC1, y=PC2)) + geom_point(size=3) + theme_bw()#


#pcoa$kefir_dataset <- "Milk.kefir"
#pcoa$kefir_dataset[which(pcoa$`kefir type.x` %in% c("WL","WG"))] <-"Water kefir"

#pcoa <- pcoa[-c(which(pcoa$`kefir type.x`=="Liquid control")),]


#total_prevalence <- rbind(species_prevalence[["milk"]],species_prevalence[["water"]])





########################################################################################################################
# get grain samples from metacache 
########################################################################################################################




grain_microbes <- c()


grain_microbes[["Milk.kefir"]] <- 
  metacache[
    which(rownames(metacache) %in%
            kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                    kefir4all_metadata$`kefir type` %in% c("MG"))]
    ),]



grain_microbes[["Water.kefir"]]<- 
  
  metacache[
    which(rownames(metacache) %in%
            kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                    kefir4all_metadata$`kefir type` %in% c("WG"))]
    ),]


grain_microbes[["Water.kefir"]] <- 
  grain_microbes[["Water.kefir"]][rownames(grain_microbes[["Water.kefir"]])=="TG_FL_T2_WG_0AN_S328",]


for (i in names(grain_microbes)){
  
 unclassified_rate=as.numeric(100-
                                                    rowSums(grain_microbes[[i]]))
  
  maxab_species <- apply( grain_microbes[[i]],2, max, na.rm=TRUE)
  
  n1_species <-names(which(maxab_species > .1))
  
  rare_biosphere_rate= rowSums(grain_microbes[[i]] [,-c(which(colnames( grain_microbes[[i]] ) %in%  n1_species ))])
    
  grain_microbes[[i]] <-  grain_microbes[[i]] [,which(colnames( grain_microbes[[i]] ) %in%  n1_species )]
  
  grain_microbes[[i]] <- grain_microbes[[i]] |> rownames_to_column("sample_id") |> pivot_longer(!sample_id,names_to="species",values_to = "relative_abundance")
  grain_microbes[[i]]$type <- i
  
  
  grain_microbes[[i]]=rbind(grain_microbes[[i]],
  data.frame(sample_id=as.character(unique(grain_microbes[[i]]$sample_id)),
             species="Unclassified",
             relative_abundance=unclassified_rate,
             type=i),
  
  
  data.frame(sample_id=as.character(unique(grain_microbes[[i]]$sample_id)),
             species="<0.1%" ,
             relative_abundance= rare_biosphere_rate,
             type=i)
  
  
  )
 
}



########################################################################################################################
# Identify microbes present in kefir grans >.1 RA
########################################################################################################################

dplyr::select(grain_microbes[[1]],species, relative_abundance)
dplyr::select(grain_microbes[[2]],species, relative_abundance)

########################################################################################################################
# Identify prevalent micorbes not in grains >.1 RA
########################################################################################################################


species_prevalence$`Milk kefir`[-c(which(species_prevalence$`Milk kefir`$specie %in% grain_microbes[["Milk.kefir"]]$species )),]

species_prevalence$`Water kefir`[-c(which(species_prevalence$`Water kefir`$specie %in% grain_microbes[["Water.kefir"]]$species )),]


print(paste(nrow(species_prevalence$`Milk kefir`[-c(which(species_prevalence$`Milk kefir`$specie %in% grain_microbes[["Milk.kefir"]]$species )),]),"prevelant species not detected at T0"))

print(paste(nrow(species_prevalence$`Water kefir`[-c(which(species_prevalence$`Water kefir`$specie %in% grain_microbes[["Water.kefir"]]$species )),]),"prevelant species not detected at T0"))



########################################################################################################################
#add this prevalence info to the p_heatmap data 
########################################################################################################################

# star labels for species absent from the T0 grain. Aligned to the loop's

# species levels; the empty top " " header level is excluded via `breaks`.

mk_levels <- fig2_species_levels[["Milk kefir"]]

mk_lab <- mk_levels

mk_t0 <- grain_microbes[["Milk.kefir"]]$species[grain_microbes[["Milk.kefir"]]$relative_abundance > 0.1]

mk_lab[!mk_lab %in% mk_t0] <- paste0("*", mk_lab[!mk_lab %in% mk_t0])

# y is NUMERIC (y_pos), so these must be scale_y_continuous -- a discrete
# scale here would silently blank the whole y axis of panels A/B.
p_heatmap$`Milk kefir` <- p_heatmap$`Milk kefir` +

  scale_y_continuous(breaks = seq_along(mk_levels), labels = mk_lab,

                     expand = expansion(mult = 0.01, add = c(0.6, 0.6)))



wk_levels <- fig2_species_levels[["Water kefir"]]

wk_lab <- wk_levels

wk_t0 <- grain_microbes[["Water.kefir"]]$species[grain_microbes[["Water.kefir"]]$relative_abundance > 0.1]

wk_lab[!wk_lab %in% wk_t0] <- paste0("*", wk_lab[!wk_lab %in% wk_t0])

p_heatmap$`Water kefir` <- p_heatmap$`Water kefir` +

  scale_y_continuous(breaks = seq_along(wk_levels), labels = wk_lab,

                     expand = expansion(mult = 0.01, add = c(0.6, 0.6)))








# ============================================================================
# Figure 2 assembly (deterministic grid composition).
#
# Each heatmap panel is drawn directly into an explicit grid viewport sized
# as a fraction of the page, so every panel STRETCHES to fill its box: the
# two panels tile the page with no wasted white space, panel titles stay in
# their white facet-strip boxes, the shared legend sits in a clean band below
# both panels, and the week blocks (wk01..wk21, T0) span the full width of
# each heatmap so even the narrow wk21 block is readable.
# ============================================================================
compose_figure2 <- function() {
  grid.newpage()

  # Panel grobs -- legends suppressed (drawn once as the shared legend at the
  # bottom of the page). Panel titles ARE the facet strips of each heatmap
  # (white box, black border) rendered in their default position at the very
  # TOP of every panel. Figure 2 has two panels (milk kefir, water kefir),
  # each one continuous heatmap spanning T0 through wk21 -- no separate
  # T0-only grain panel is needed, since every heatmap already includes its
  # T0 grain samples alongside the later liquid/grain timepoints.
  gC <- ggplotGrob(p_heatmap[[1]] + theme(legend.position = "none"))
  gD <- ggplotGrob(p_heatmap[[2]] + theme(legend.position = "none"))

  # Shared legend: Relative-abundance colourbar + kefir-type key.  Both C/D
  # heatmaps carry an identical legend, so take it from the milk-kefir panel.
  leg <- NULL
  try(leg <- ggpubr::get_legend(p_heatmap[[1]]), silent = TRUE)
  if (is.null(leg)) try(leg <- ggpubr::get_legend(p_heatmap[[2]]), silent = TRUE)

  # --- shared legend band, at the very bottom of the page (well below the
  #     heatmaps) so it can never be mistaken for a panel title ---
  # npc (page-fraction) coordinates: see the EXPORT FIGURE 2 note below for why
  # these don't need to change with canvas size.
  if (!is.null(leg)) {
    tryCatch(
      {
        pushViewport(viewport(x = unit(0.5, "npc"), y = unit(0.078, "npc"),
                              width = unit(0.97, "npc"), height = unit(0.156, "npc"),
                              just = c("centre", "centre")))
        grid.draw(leg)
        popViewport()
      },
      error = function(e) {
        try(popViewport(), silent = TRUE)
        message("shared legend skipped: ", conditionMessage(e))
      }
    )
  }

  # --- panel viewport helper ---
  vp_panel <- function(x0, y0, x1, y1) {
    viewport(x = unit(x0, "npc"), y = unit(y0, "npc"),
             width = unit(x1 - x0, "npc"), height = unit(y1 - y0, "npc"),
             just = c("left", "bottom"))
  }

  # A: milk kefir heatmap; B: water kefir heatmap
  pushViewport(vp_panel(0.015, 0.106, 0.490, 0.929)); grid.draw(gC); popViewport()
  pushViewport(vp_panel(0.510, 0.106, 0.985, 0.929)); grid.draw(gD); popViewport()

  # --- panel labels (A. B.) ---
  panel_label <- function(txt, x, y) {
    grid.text(txt, x = unit(x, "npc"), y = unit(y, "npc"),
              gp = gpar(fontsize = 30, fontface = "bold"),
              just = c("left", "top"))
  }
  panel_label("A.", 0.015, 0.940)
  panel_label("B.", 0.510, 0.940)
  # Milk / water kefir title on top of each panel
  grid.text("Milk kefir", x = unit(0.2525, "npc"), y = unit(0.965, "npc"),
            gp = gpar(fontsize = 22, fontface = "bold"),
            just = c("center", "top"))
  grid.text("Water kefir", x = unit(0.7475, "npc"), y = unit(0.965, "npc"),
            gp = gpar(fontsize = 22, fontface = "bold"),
            just = c("center", "top"))

  invisible(NULL)
}

# ============================================================================
# EXPORT FIGURE 2
# ============================================================================
# NOTE: compose_figure2() positions panels with npc (page-fraction) viewports
# and uses ABSOLUTE point font sizes that were tuned for a 26.2 x 17.3 in
# canvas (7864 x 5200 px @ 300 dpi, pointsize = 15) -- matching the canvas
# other full-page figures in this project use (e.g. Figure 7). The heatmap
# panels' row heights come from npc fractions of this canvas, so the full
# height gives each species row enough room to read clearly; a shorter canvas
# compresses every row by the same proportion. Vector PDF text scales cleanly,
# so it can be re-placed at any print size (e.g. 178 mm wide) downstream
# without re-rendering.

pdf(
  file = file.path(FIGURES_DIR, "Figure_2.pdf"),
  width  = 7864 / 300,          # 26.213 in
  height = 5200 / 300,          # 17.333 in
  pointsize = 15,               # match the jpeg() pointsize
  useDingbats = FALSE
)

compose_figure2()

dev.off()

# High-resolution raster version on the SAME canvas as the PDF above,
# LZW-compressed TIFF for publication.
tiff(
  filename   = file.path(FIGURES_DIR, "Figure_2.tiff"),
  width      = 7864,
  height     = 5200,
  units      = "px",
  res        = 300,
  pointsize  = 15,
  compression = "lzw",
  type       = "cairo"
)

compose_figure2()

dev.off()



########################################################################################################################
# Explore the unclassified component
########################################################################################################################
cat("\n=== UNCLASSIFIED COMPONENT ANALYSIS ===\n")

# Compute unclassified rates for ALL T0 grain samples from metacache
grain_samples_mk <- kefir4all_metadata$merge_column[kefir4all_metadata$Stage == "T0" & kefir4all_metadata$`kefir type` %in% c("MG")]
grain_samples_wk <- kefir4all_metadata$merge_column[kefir4all_metadata$Stage == "T0" & kefir4all_metadata$`kefir type` %in% c("WG")]

get_unclassified_rate <- function(samples, metacache) {
  if (length(samples) == 0) return(NULL)
  mat <- metacache[rownames(metacache) %in% samples, , drop = FALSE]
  mat[is.na(mat)] <- 0
  data.frame(
    sample_id = rownames(mat),
    unclassified_pct = 100 - rowSums(mat),
    n_species_detected = rowSums(mat > 0),
    stringsAsFactors = FALSE
  )
}

mk_unclass <- get_unclassified_rate(grain_samples_mk, metacache)
wk_unclass <- get_unclassified_rate(grain_samples_wk, metacache)
if (!is.null(mk_unclass)) mk_unclass$kefir_type <- "Milk kefir"
if (!is.null(wk_unclass)) wk_unclass$kefir_type <- "Water kefir"
unclassified_all <- bind_rows(mk_unclass, wk_unclass)

if (is.null(unclassified_all) || nrow(unclassified_all) == 0) {
  cat("No grain samples found for unclassified analysis.\n")
} else {
  cat("\n1. SUMMARY OF UNCLASSIFIED READS IN T0 GRAIN\n")
  cat("   Total grain samples:", nrow(unclassified_all), "\n")
  cat(sprintf("   Range: %.1f%% - %.1f%%\n", min(unclassified_all$unclassified_pct), max(unclassified_all$unclassified_pct)))
  cat(sprintf("   Mean: %.1f%% (SD: %.1f)\n", mean(unclassified_all$unclassified_pct), sd(unclassified_all$unclassified_pct)))

  cat("\n2. BY KEFIR TYPE:\n")
  type_summary <- unclassified_all |>
    group_by(kefir_type) |>
    summarise(n = n(), mean = mean(unclassified_pct), median = median(unclassified_pct),
              sd = sd(unclassified_pct), min = min(unclassified_pct), max = max(unclassified_pct), .groups = "drop")
  print(as.data.frame(type_summary))

  cat("\n3. STATISTICAL COMPARISON:\n")
  mk_rates <- unclassified_all |> filter(kefir_type == "Milk kefir") |> pull(unclassified_pct)
  wk_rates <- unclassified_all |> filter(kefir_type == "Water kefir") |> pull(unclassified_pct)
  if (length(mk_rates) > 1 & length(wk_rates) > 1) {
    wt <- wilcox.test(mk_rates, wk_rates)
    cat(sprintf("   Wilcoxon test: p = %.4f\n", wt$p.value))
    cat(sprintf("   Milk kefir median: %.1f%% (n=%d)\n", median(mk_rates), length(mk_rates)))
    cat(sprintf("   Water kefir median: %.1f%% (n=%d)\n", median(wk_rates), length(wk_rates)))
  }

  cat("\n4. CORRELATION: UNCLASSIFIED RATE vs. SPECIES DETECTED:\n")
  cor_test <- cor.test(unclassified_all$unclassified_pct, unclassified_all$n_species_detected, method = "spearman")
  cat(sprintf("   Spearman rho: %.3f (p = %.4f)\n", cor_test$estimate, cor_test$p.value))

  cat("\n5. GRAIN VS. LIQUID COMPARISON:\n")
  liquid_idx <- kefir4all_metadata$Stage != "T0"
  liquid_samples <- kefir4all_metadata$merge_column[liquid_idx]
  liquid_kefir_types <- kefir4all_metadata$`kefir type`[liquid_idx]
  names(liquid_kefir_types) <- liquid_samples
  get_liquid_unclass <- function(samples, kefir_types_named, metacache) {
    if (length(samples) == 0) return(NULL)
    mat <- metacache[rownames(metacache) %in% samples, , drop = FALSE]
    mat[is.na(mat)] <- 0
    matched_types <- kefir_types_named[rownames(mat)]
    data.frame(sample_id = rownames(mat), unclassified_pct = 100 - rowSums(mat),
               n_species_detected = rowSums(mat > 0),
               kefir_type = ifelse(matched_types %in% c("MG", "ML"), "Milk kefir", "Water kefir"),
               stringsAsFactors = FALSE)
  }
  liquid_unclass <- get_liquid_unclass(liquid_samples, liquid_kefir_types, metacache)
  if (!is.null(liquid_unclass) && nrow(liquid_unclass) > 0) {
    combined_compare <- bind_rows(
      unclassified_all |> select(sample_id, unclassified_pct, n_species_detected, kefir_type) |> mutate(compartment = "Grain (T0)"),
      liquid_unclass |> select(sample_id, unclassified_pct, n_species_detected, kefir_type) |> mutate(compartment = "Liquid")
    )
    compare_summary <- combined_compare |>
      group_by(compartment, kefir_type) |>
      summarise(n = n(), mean = mean(unclassified_pct), median = median(unclassified_pct),
                sd = sd(unclassified_pct), .groups = "drop")
            print(as.data.frame(compare_summary))
  }
}



cat("\n=== END UNCLASSIFIED ANALYSIS ===\n")
