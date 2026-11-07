# Anchor to the project root before here() loads. RStudio Server resets the session cwd
# to $HOME after sourcing .Rprofile, and here() computes its root AT PACKAGE LOAD
# TIME, so force the project root before any here load. This block is inert on non-portable
for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {



  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}


# ==============================================================================
# PATH CONFIGURATION — portable relative paths via the here package
# Run install.packages("here") if not installed.
# For files listed in data/private/README.md, provide your own copies.
# ==============================================================================
if (!requireNamespace("here",   quietly = TRUE)) install.packages("here")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
source(here::here("scripts", "_lab_toolkit", "plotting_report.R"))
kefir_use_style()
DATA_DIR    <- here::here("data")                         # public data files
PRIVATE_DIR <- here::here("data", "private")              # private/sensitive files
OUT_DIR <- here::here("output_reviewed", "06_strain_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "Figure_7")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# Private metadata — not distributed with repo (participant identifiers)
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






###############################################################################################################
#look into the genetic distance between clusters
###############################################################################################################
library(ggalluvial)   # geom_alluvium/geom_stratum/geom_flow
library(ggpubr)       # ggarrange
library(gtools)       # mixedsort
library(stringr)      # str_split_fixed
library(tibble)       # column_to_rownames/rownames_to_column
library(tidyr)        # spread/gather
library(widyr)        # pairwise_count
library(readr)
library(dplyr)
library(ggplot2)

instrain <- read_csv(
  file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv")
)

########################################################################################################################
#Import sample metadata
########################################################################################################################

global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv"))
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv"))
global_mk_metadata$Stage <- NA
global_wk_metadata$Stage <- NA

Citizen_Scientist_metadata_v8 <- read_csv(CS_METADATA_PRIVATE)

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"))
kefir4all_metadata$merge_column <-  gsub("_host_removed_R..fastq.gz","",kefir4all_metadata$merge_column)
kefir4all_metadata <- kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$merge_column))),]
########################################################################################################################
# Merge metadata into one
########################################################################################################################
total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`,Stage,Sample),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`,Stage,Sample),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`,Stage,Sample))


total_metadata $category <- NA
total_metadata $category[which(total_metadata $`kefir type` %in% c("WL","WG"))] <- "Water.kefir"
total_metadata $category[which(total_metadata $`kefir type` %in% c("ML","MG"))] <- "Milk.kefir"




########################################################################################################################
# Import prevalence metadata
########################################################################################################################

milk_taxonomic_profile_prevalence <- read_csv(file.path(DATA_DIR, "milk_taxonomic_profile_prevalence.csv"))
water_taxonomic_profile_prevalence <- read_csv(file.path(DATA_DIR, "water_taxonomic_profile_prevalence.csv"))
total_prevalence <- rbind(milk_taxonomic_profile_prevalence, water_taxonomic_profile_prevalence)

# WHY: the prevalence tables split a few species into dRep/StrainPhlAn-only
# "subcluster"/"subspecies" variants that inStrain reports under one merged
# name -- without this relabelling, these show up as extra fragments instead
# of pooling with the rest of that species' detections.
total_prevalence$...5[which(total_prevalence$...5 == "Lactobacillus_ghanensis")] <- "Liquorilactobacillus ghanensis"
total_prevalence$...5[which(total_prevalence$...5 == "Lactobacillus_satsumensis")] <- "Liquorilactobacillus satsumensis"
total_prevalence$...5[which(total_prevalence$...5 == "Lactococcus_lactis subcluster 1")] <- "Lactococcus lactis"
total_prevalence$...5[which(total_prevalence$...5 == "Lactococcus_lactis subcluster 2")] <- "Lactococcus cremoris"
total_prevalence$...5[which(total_prevalence$...5 == "Pseudomonas_fragi_subspecies 1")] <- "Pseudomonas fragi"
total_prevalence$...5[which(total_prevalence$...5 == "Zymomonas_mobilis_subcluster 1")] <- "Zymomonas mobilis"
total_prevalence$...5[which(total_prevalence$...5 == "Lactobacillus_perolens")] <- "Schleiferilactobacillus perolens"
# WHY left un-gsubbed here: every downstream comparison against total_prevalence$...5 in
# this file applies its own gsub("_", " ", ...) inline at the point of use (see the
# occurence_data_v2 / multi_strain_detection_prevelant filters below). Applying the
# substitution globally here instead changes which species pass the
# instrain_prevelant_cs_clust filter further down this script -- confirmed empirically --
# so the underscore form is kept here as the canonical value this column holds throughout.

# relationship: a left merge (all.x=TRUE) must preserve instrain's row count exactly.
instrain_rows_before_merge <- nrow(instrain)
instrain <- merge(instrain, total_metadata, by.x = "sample_id", by.y = "merge_column", all.x = TRUE)
stopifnot("merging total_metadata into instrain changed the row count -- check merge_column for duplicates" =
            nrow(instrain) == instrain_rows_before_merge)

#we applied the authors' recommended thresholds of 99.999% identity and 50% coverage.

colnames(global_mk_metadata)[which(colnames(global_mk_metadata) == "country")] <- "Sample_grain"
global_wk_metadata$Sample_grain <- global_wk_metadata$Sample

global_metadata <- list(Milk.kefir = global_mk_metadata, Water.kefir = global_wk_metadata)





multi_strain_detection <- c()
occurence_data <- c()
cluster_size_total <- c()

for (species in levels(as.factor(instrain$classification))) {

  for (type in c("Milk.kefir", "Water.kefir")) {

    # BUG FIX: this script previously applied only the popANI_reference > 0.98 cut,
    # missing the breadth >= 0.35 and observed:expected breadth ratio >= 0.75 criteria
    # that the Figure 7 caption states and that 04_supp_community_types.R and
    # Supplementary_Note_6_analysis.R already apply for the same per-genome detection
    # decision -- the three filters together are the Methods-stated detection criterion.
    t1 <- instrain[which(instrain$classification == species &
      instrain$category.y == type &
        instrain$breadth >= 0.35 &
        (instrain$breadth / instrain$breadth_expected) >= 0.75 &
        instrain$popANI_reference > .98), ]

    if (nrow(t1) == 0) {
      next
    }
    t1 <- dplyr::select(t1, "sample_id", "genome", "coverage", "breadth",
                        "Stage.y", "Sample.y", "category.y", "data_source.y",
                        "cluster", "classification")

    # relationship: all.y=TRUE preserves every row of t1; global_metadata[[type]] is
    # one row per merge_column, so the result should have exactly nrow(t1) rows.
    t1_rows_before_merge <- nrow(t1)
    t1 <- merge(global_metadata[[type]], t1, by.y = "sample_id", by.x = "merge_column", all.y = TRUE)
    stopifnot("merging global_metadata into t1 changed the row count -- check merge_column for duplicates" =
                nrow(t1) == t1_rows_before_merge)

    t1$Sample_grain[grep("ID", t1$Sample.y)] <- t1$Sample.y[grep("ID", t1$Sample.y)]

   t2 <- as.data.frame(table(t1$merge_column))

    
    multi_strain_detection=rbind(multi_strain_detection,
                                 data.frame(species=species,
                                            kefir=type,
                                            n_total_detections=    nrow(t2), 
                                            n_multi_detections= length(which(t2$Freq>1))
                                 
                                            
                                            ))
    
    
    if(length(which(t2$Freq>1))==0){
      
      print(paste( " We did not detect more than one strain cluster in more than one sample for the",species, "within the", type, "microbiome"))
      next     }
    
    
    
    
    
    
    
    
    cooccurrence_matrix <-
      dplyr::select(t1, merge_column, classification, cluster) |>
      mutate(cluster = paste(classification, cluster)) |>
      pairwise_count(cluster, merge_column) |>
      spread(item2, n, fill = 0) |>
      column_to_rownames(var = "item1")

    cooccurrence_long <- as.data.frame(cooccurrence_matrix) |>
      rownames_to_column("Var2") |>
      pivot_longer(!Var2, names_to = "Var1", values_to = "Freq")

    cluster_size <- dplyr::select(t1, merge_column, classification, cluster) |>
      dplyr::count(cluster)
    
    
    
    
#     
#     
#     
#     # Plot the network with vertex labels as clusters
# 
#     # Plot the network
#     
#     
#     
#     
#     
#     
#     
#     
#     # Set the diagonal to 0 (no self-co-occurrence)
#     
#     # Create a graph object from the co-occurrence matrix
#     
#     # Convert igraph object to tidygraph
#     
#     # Step 2: Plot the Co-occurrence Network
#     set.seed(42)  # For reproducibility of the layout
#         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#         legend.position = "none"
#     
#     
#     
#     
#     graphics.off()
    
    cooccurrence_long$species=species
    cooccurrence_long$kefir=type
    
    
    cluster_size$species=species
    cluster_size$kefir=type
    
    cluster_size$var= paste(species,cluster_size$cluster)
    
    
    
occurence_data=rbind(               cooccurrence_long,
                                    occurence_data)
    
    cluster_size_total <- rbind( cluster_size,
                                 cluster_size_total)



    
    
  }
  
  
}



occurence_data_v2=
  occurence_data[which(occurence_data$Freq>=5),]

occurence_data_v2=
rbind(
occurence_data_v2[which(occurence_data_v2$species %in% gsub("_"," ",total_prevalence$...5[which(total_prevalence$kefir_type=="water")])&
                       
                       occurence_data_v2$kefir=="Water.kefir"
                       ),],




occurence_data_v2[which(occurence_data_v2$species %in% gsub("_"," ",total_prevalence$...5[which(total_prevalence$kefir_type=="milk")])&
                       
                       occurence_data_v2$kefir=="Milk.kefir"
),] ) 
  
  
 plot_data=  
   occurence_data_v2
 
 
  plot_data$Var1_c=
  
  str_split_fixed(plot_data$Var1,"_",2)[,2]


plot_data$Var2_c=
  
  str_split_fixed(plot_data$Var2,"_",2)[,2]

plot_data$Var1_primary=
  
  gsub(".* ","",plot_data$Var1)
 
plot_data$Var2_primary=
  
  gsub(".* ","",plot_data$Var2)

plot_data$kefir <- gsub("\\."," ",plot_data$kefir )

  
  
b=

ggplot(plot_data, aes(axis1 = Var1_primary, axis2 = Var2_primary, y = Freq)) +
  geom_alluvium(aes(fill = factor(Var1_c, levels = mixedsort(unique(Var1_c))))) +
  geom_stratum() +
  labs(title = "",
       x = "",
       y = "Number of strain Detections",
       fill = "Secondary cluster") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  facet_wrap(~species + kefir, scales = "free") +
  theme_bw() +
  theme(
    legend.position = c(0.78, 0.1),  # Adjust position to bottom left corner
    legend.justification = "left",  # Align legend to the left
    legend.direction = "horizontal", # Make the legend horizontal
    legend.box = "horizontal",  # Arrange in a horizontal box
    legend.box.just = "left",  # Align legend box to the left
    legend.spacing.x = unit(0.5, 'cm'), # Adjust spacing between items
    legend.spacing.y = unit(0.5, 'cm'), # Adjust spacing between rows
    legend.text = element_text(face = "italic", size = 20),
    legend.key.size = unit(1.5, 'cm'),  # Adjust legend key size
    legend.key.height = unit(1.5, 'cm'),  # Adjust legend key height
    legend.key.width = unit(1.5, 'cm'),  # Adjust legend key width
    legend.title = element_text(size = 20,hjust = 0.5),
    #legend.title.align = 0.5,  # Center-align title
    axis.text.y = element_text(size = 15, face = "italic"),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid labels
    panel.grid.minor = element_blank(),  # Remove minor grid labels
    plot.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "white", colour = "black")
  )+
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5))


# 
# graphics.off()
# #



  #theme_void()




multi_strain_detection_prevelant=
  rbind(
    multi_strain_detection[which(multi_strain_detection$species %in% gsub("_"," ",total_prevalence$...5[which(total_prevalence$kefir_type=="water")])&
                              
                                   multi_strain_detection$kefir=="Water.kefir"
    ),],
    
    
    multi_strain_detection[which(multi_strain_detection$species %in% gsub("_"," ",total_prevalence$...5[which(total_prevalence$kefir_type=="milk")])&
                                   
                                   multi_strain_detection$kefir=="Milk.kefir"
    ),]
     ) 




a=
multi_strain_detection_prevelant |> 
  mutate(kefir=gsub("\\."," ",kefir)) |> 
  filter(n_multi_detections !=0) |> 
tidyr::gather(key = "type", value = "count", n_total_detections, n_multi_detections) |> 
  mutate(type=
                gsub( "n_multi_detections", "Multi-Strain Detections",
                      gsub("n_total_detections" ,"Total Detections",type))) |> 



  ggplot(aes(x = reorder(species, -count), y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Species",
       y = "Number of Detections",
       fill = "Detection Type") +
  coord_flip() +
  facet_wrap(~ kefir, scales = "free_y") +
  theme_bw()+
  theme( legend.position = c(0.35, .7), 
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.y = element_text(size=15,face="italic"),
        axis.title = element_text(size = 20),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(size=20,hjust = 1), # remove x-axis labels
        #axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.text=element_text(size = 20),
        legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(1.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'), #change legend key width
        legend.title = element_text(size=20),
        strip.text = element_text(size = 18, face = "bold"),
        strip.background = element_rect(fill = "white", colour = "black"))




# 
# graphics.off()
# #


library(ggpubr)

# This script's sole output is Figure 7, a two-panel strain-detection
# summary saved below (the variable name fig9_combined reflects the
# panel's internal figure number, not its final manuscript numbering).
tiff(filename=file.path(FIGURES_DIR, 'Figure_7.tiff'), width = 7864, height = 5200, units = "px", res = 300, compression = "lzw", type = "cairo", pointsize = 15)


fig9_combined <- ggarrange(a, b,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(size = 30),heights = c(.4,1), common.legend = FALSE)
print(fig9_combined)
graphics.off()

save_figure(fig9_combined, file.path(FIGURES_DIR, "Figure_7"),
            width_mm = 7864 / 300 * 25.4, height_mm = 5200 / 300 * 25.4,
            source = file.path(DATA_DIR, "instrain_genome_species_primary_data_v4.csv"),
            seed = NULL)

# This script produces only Figure 7. Figure 8's alluvial plot (a different
# figure, built from the same xtabs(~Stage+clusters+classification+category.y)
# style of analysis) is produced by its own dedicated script
# (06_instrain_temporal_alluvial.R), not duplicated here.
