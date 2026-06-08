# ---- Portable path bootstrap (added during repo migration) ----
if (!requireNamespace("here",   quietly = TRUE)) install.packages("here")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(here)
DATA_DIR    <- here::here("data")
PRIVATE_DIR <- here::here("data", "private")
FIGURES_DIR <- here::here("figures")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
CS_METADATA_PRIVATE <- file.path(PRIVATE_DIR, "Citizen Scientist metadata_v8.csv")
if (!file.exists(CS_METADATA_PRIVATE)) {
  stop("Private metadata not found: ", CS_METADATA_PRIVATE,
       "\nCopy to data/private/ to run this script.")
}
# ---- End bootstrap ----

﻿



###############################################################################################################
#look into the genetic distance between clusters
###############################################################################################################

pacman::p_load(rlang,tibble,ape,colorspace,concaveman,ggnewscale,readxl,hrbrthemes,Biostrings,ggtree,flextable,devtools,R4RNA,taxize,rotl,ape,treeio,DECIPHER,ggdendro,ggplot2,tidyr,RSQLite,optmatch,rentrez,dplyr,seqinr,RColorBrewer,ggtext)
pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis)
pacman::p_load(cutpointr,readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)
#remotes::install_github("YuLab-SMU/ggtree")
# if (!requireNamespace("devtools", quietly=TRUE))
#   install.packages("devtools")
# devtools::install_github("YuLab-SMU/ggmsa")
library(ggmsa)
library(readr)
library(dplyr)
library(ggplot2)



########################################################################################################################
#Import metadata
########################################################################################################################

global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv")
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv")
global_mk_metadata$Stage <- NA
global_wk_metadata$Stage <- NA

Citizen_Scientist_metadata_v8 <- read_csv(CS_METADATA_PRIVATE)

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv")
kefir4all_metadata$merge_column <-  gsub("_host_removed_R..fastq.gz","",kefir4all_metadata$merge_column)
kefir4all_metadata <- kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$merge_column))),]

marker_names <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/06_strain_profiling//marker_names.csv")

marker_names$x <- paste("t__",marker_names$x,sep="")
########################################################################################################################
# Merge metadata into one
########################################################################################################################
total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`,Stage,Sample),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`,Stage,Sample),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`,Stage,Sample))


total_metadata $category <- NA
total_metadata $category[which(total_metadata $`kefir type` %in% c("WL","WG"))] <- "Water.kefir"
total_metadata $category[which(total_metadata $`kefir type` %in% c("ML","MG"))] <- "Milk.kefir"


###############################################################################################################
#For loop to add metadata to a dataframe list
###############################################################################################################

my.list_metadata <- vector(mode = "list", length = length(names(myfiles)))
names(my.list_metadata) <- names(myfiles)

#i <- "Lactococcus_lactis"
pos <- c()
pos.1 <- c()

for (i in names(myfiles)){
  myfiles[[i]]$tip.label <- gsub("-","_",myfiles[[i]]$tip.label)
  pos <- which(total_metadata$merge_column %in%  myfiles[[i]]$tip.label)
  
  my.list_metadata[[i]] <- total_metadata[pos,]
  my.list_metadata[[i]] <-   my.list_metadata[[i]] %>% relocate(merge_column, .before = data_source)
  
  
}


###############################################################################################################
#look into the SNPs between strains
###############################################################################################################

setwd("Q:/H2020 Master/Citizen Science Project/Results/06_strain_profiling/06_strainphlan/06_strainphlan_sensitive/06_mutation")
temp_phylo_mut = list.files(pattern=".mutation", recursive = FALSE, full.names = TRUE,include.dirs = TRUE)
myfiles_mut = lapply(temp_phylo_mut,read_table)

names(myfiles_mut) <- gsub("./|.tsv|.mutation","",temp_phylo_mut)

unit_values=c()

strainphlan_mut <-c()
for (species in names(myfiles_mut)){
  if(length(which(myfiles_mut[[species]]$ids %in% kefir4all_metadata$merge_column ))==0){
    print(paste(species, "only found in global dataset"))
    next
    
  }
  # 
  # npte subset will remove strain comparisions between strain fron the global studies
  
 myfiles_mut[[species]]$ids <-  gsub("-","_",myfiles_mut[[species]]$ids)
 myfiles_mut[[species]]$ids <- gsub("-","_",myfiles_mut[[species]]$ids)
  
  
  
  
  t <- 
   myfiles_mut[[species]] %>% 
    mutate_all(as.character) %>% 
   
    pivot_longer(cols = -ids, names_to = "sample_id", values_to = "snps") %>% 
   # mutate(X3=as.numeric(X3)) %>% 
    #myfiles_mut[[species]][which(myfiles_dist[[species]]$X1 %in% kefir4all_metadata$merge_column ),] %>% 
    #myfiles_mut[[species]]$X2 %in% kefir4all_metadata$merge_column),] %>% 
    #     mutate(across(everything(), as.character)) %>% 
    #     pivot_longer(!ids,names_to = "sample_id",values_to = "snps") %>% 
    merge(.,total_metadata,by.x="ids",by.y="merge_column",all.x=TRUE)%>% 
    merge(.,total_metadata,by.x="sample_id",by.y="merge_column",all.x=TRUE) %>% 
    mutate(comparision_type="Unrelated")
  
  
  t$snps_1 <- NA
  
  t$snps_2 <- NA
  t$snps_updated <-   t$snps
  
  
  
  t$snps_1[grep("/",t$snps)] <- str_split_fixed(t$snps[grep("/",t$snps)],"/",2)[,1]
  t$snps_2[grep("/",t$snps)] <- str_split_fixed(t$snps[grep("/",t$snps)],"/",2)[,2]
  
  t$snps_updated[grep("/",t$snps)] <- as.numeric(t$snps_1[grep("/",t$snps)]) /  as.numeric(t$snps_2[grep("/",t$snps)])
  
  t$snps_updated <- as.numeric(t$snps_updated)
  
  
  #for a comparision to be defined as unrelated it can  A. from a different datasets B. different ID in the global dataset C. different categories in the CS dataset the same sample ID (note this includes the global dataset)
  
  # t$comparision_type[which(t$`kefir type.x`==t$`kefir type.y`)] <- "Intra"
  
  
  
  #for a comparision to be defined as inter it needs to be in the cs dataset and be of the same category e.g. milk or water kefir
  # note this command will override the intra command below so run it first
  t$comparision_type[which(t$data_source.x=="This study" &
                             t$data_source.y=="This study" &
                             t$category.x==t$category.y)] <- "Inter"
  
  
  #for a comparision to be defined as intra it needs to be from the same sample ID (note this includes the global dataset)
  t$comparision_type[which(t$Sample.x==t$Sample.y)] <- "Intra"
  
  
  

  
  for (kefir in levels(as.factor( t$category.x))){
    
    
    t1=t[which(t$category.x==kefir),]
    
    if(TRUE %in% is.na(t1$snps_updated)){
    t1= t1[-c(which(is.na(t1$snps_updated))),]
    
    }
    if(TRUE %in% is.nan(t1$snps_updated)){
      t1= t1[-c(which(is.nan(t1$snps_updated))),]
      
    }
      
    t1$species=as.character(    marker_names$...2[which(marker_names$x==species)])
    t1$type=kefir
    unit_values <- rbind(unit_values,
                         data.frame(species=as.character(    marker_names$...2[which(marker_names$x==species)]),
                                    type=as.character(kefir),
                                    mean=as.numeric(mean(t1$snps_updated) ), 
                                    min=as.numeric(min (t1$snps_updated) ),
                                    max=as.numeric(max(t1$snps_updated) ),
                                    min=as.numeric(min(t1$snps_updated)),
                                    n_strains=as.numeric(length(unique(levels(as.factor(t1$sample_id)),levels(as.factor(t1$ids)))))
                                    
                         ))
    
    
    strainphlan_mut <- 
    rbind(strainphlan_mut,
          t1)
    
  }
  
  
  
}




milk_taxonomic_profile_prevalence <- read_csv(file.path(DATA_DIR, "milk_taxonomic_profile_prevalence.csv"))

water_taxonomic_profile_prevalence <- read_csv(file.path(DATA_DIR, "water_taxonomic_profile_prevalence.csv"))


total_prevalence <- rbind(milk_taxonomic_profile_prevalence,water_taxonomic_profile_prevalence)



total_prevalence$...5[which(total_prevalence$...5=="Lactobacillus_ghanensis")] ="Liquorilactobacillus ghanensis"
#total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 1" )] ="Lactococcus lactis" 
#total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 2")] ="Lactococcus cremoris"

total_prevalence$...5[which(total_prevalence$...5=="Zymomonas mobilis subcluster 1")] ="Zymomonas mobilis"



# total_prevalence$...5 <- 
#   gsub("_"," ",total_prevalence$...5)


unit_values_prevalence <- 
  
  rbind(
    
    unit_values[which(unit_values$type=="Milk.kefir" &
                        unit_values$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]),],
    
    unit_values[which(unit_values$type=="Water.kefir" &
                        unit_values$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]),]
    
    
    
    
  )

ggplot(unit_values_prevalence %>% 
         
         mutate(species=gsub("Lactococcus_lactis subcluster 1", "Lactococcus lactis",species )) %>% 
         mutate(species=gsub("Lactococcus_lactis subcluster 2", "Lactococcus cremoris",species )) %>% 
         mutate(species=gsub("Zymomonas_mobilis_subcluster 1", "Zymomonas mobilis",species )) %>% 
         mutate(species=gsub("Pseudomonas_fragi_subspecies 1", "Pseudomonas fragi",species )) %>% 
         mutate(species=gsub("_", " ",species )
                
         ),
       aes(x = species, y = mean, ymin = min, ymax = max)) + 
  geom_pointrange(colour="red") +
  #geom_label(aes(label = n_strains), vjust = -1, size = 5, label.padding = unit(0.2, "lines")) +  # Add text labels in boxes
  geom_label(aes(y = max + 0.05 * (max - mean), label = n_strains), vjust = 0, size = 5, label.padding = unit(0.2, "lines")) +  # Adjust the position of n_strains label
  #geom_text(aes(label = n_strains), vjust = -1, size = 5) +  # Add text labels
  theme_minimal() +
  facet_wrap(~type, scales = "free") +
  theme_bw() +
  xlab("Species") +
  ylab("Phylogenetic dist") +
  labs(color = "Number of Strains") +
  theme(plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
        legend.title = element_text(size = 25, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text(size = 15, face = "bold", hjust = 0.5, vjust = -1),
        axis.title.y = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 2)
  )




strainphlan_mut_prevalence <- 
  
  rbind(
    
    strainphlan_mut[which(strainphlan_mut$type=="Milk.kefir" &
                        strainphlan_mut$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]),],
    
    strainphlan_mut[which(strainphlan_mut$type=="Water.kefir" &
                        strainphlan_mut$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]),]
    
    
    
    
  )


strainphlan_plot <- 
ggplot(strainphlan_mut_prevalence  %>% 
         
         mutate(species=gsub("Lactococcus_lactis subcluster 1", "Lactococcus lactis",species )) %>% 
         mutate(species=gsub("Lactococcus_lactis subcluster 2", "Lactococcus cremoris",species )) %>% 
         mutate(species=gsub("Zymomonas_mobilis_subcluster 1", "Zymomonas mobilis",species )) %>% 
         mutate(species=gsub("Pseudomonas_fragi_subspecies 1", "Pseudomonas fragi",species )) %>% 
         mutate(species=gsub("_", " ",species )) %>% 
                mutate(type=gsub("\\.", " ",type )
         ),
       aes(x = species, y = snps_updated)) + 
 # geom_jitter(color="black", size=0.4, alpha=0.9) +
  geom_boxplot()+
  #geom_label(aes(label = n_strains), vjust = -1, size = 5, label.padding = unit(0.2, "lines")) +  # Add text labels in boxes
  #geom_label(aes(y = max + 0.05 * (max - mean), label = n_strains), vjust = 0, size = 5, label.padding = unit(0.2, "lines")) +  # Adjust the position of n_strains label
  #geom_text(aes(label = n_strains), vjust = -1, size = 5) +  # Add text labels
  theme_minimal() +
  facet_wrap(~type, scales = "free") +
  theme_bw() +
  xlab("Species") +
  ylab("Rate of nucleotide change") +
  labs(color = "Number of Strains") +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15,face="italic"),
    strip.background = element_rect(
      color="black", fill="white"),
    strip.text = element_text(size=15.5),
    #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
    axis.title.y = element_text( size=20, face="bold",hjust = 0.5, vjust = 1.5),
    legend.key.size = unit(2, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm'), #change legend key width
    legend.direction = "horizontal", 
    legend.position = "top",
    legend.box = "horizontal",
    # legend.text = element_text(size=30),
    legend.title = element_text(size=17.5),
    legend.text=element_text(size=17.5))







###############################################################################################################
#look into the ANI between MAGs
###############################################################################################################






########################################################################################################################
#Import MAG metadata 
########################################################################################################################

wdb <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/drep/data_tables/Wdb.csv")
wdb$genome <- substr(wdb$genome,1,nchar(wdb$genome)-3)
cdb <- read_csv(file.path(DATA_DIR, "Cdb.csv"))
cdb$genome <- substr(cdb$genome,1,nchar(cdb$genome)-3)

Ndb <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/drep/data_tables/Ndb.csv")
Ndb$querry <- substr(Ndb$querry,1,nchar(Ndb$querry)-3)

Ndb$reference <- substr(Ndb$reference,1,nchar(Ndb$reference)-3)


#ndb <- ndb[-c(which(ndb$querry=="Lactococcus_lactis_GCA_015476255.1_ASM1547625v1_genomic.fa")),]
mag_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/mag_metadata_processed_v1.csv")


mag_metadata$bin <- 
  gsub(".*_bins_","",mag_metadata$user_genome)

mag_metadata$bin <- 
  gsub(".orig|.permissive|.strict","",mag_metadata$bin)


mag_metadata$bin <- 
  gsub(".*_","",mag_metadata$bin)


mag_metadata$bin_id <- paste(mag_metadata$base_name,"_",mag_metadata$bin,sep="")




mag_metadata$classification_full <- mag_metadata$classification
mag_metadata$classification <- gsub(".*;s__","",mag_metadata$classification )



# mag_metadata <- 
#   mag_metadata [-c(which(mag_metadata $Completeness <=80 |
#                            mag_metadata $Contamination >=5)),]



wdb <- merge(wdb,mag_metadata, 
             by.x="genome",
             by.y="user_genome",
             all.x=TRUE)



wdb$primary <- gsub("_.","",wdb$cluster)

wdb$species_rep <- NA

for (i in levels(as.factor(wdb$primary))){
  
  wdb$species_rep [which(wdb$primary==i &   wdb$score ==
                           max(
                             wdb$score[which(wdb$primary==i)]))] <- "*"
  
}








Cdb <- read_csv(file.path(DATA_DIR, "Cdb.csv"))


# import relevant data
Cdb$genome <- gsub(".fa","",Cdb$genome)

nrow( mag_metadata)
Cdb <-  merge(mag_metadata,
              Cdb,
              by.x="user_genome",
              by.y="genome",
              all.x=TRUE)



nrow( Cdb)


Cdb[
  which(is.na(Cdb$secondary_cluster)),]


Cdb$classification[which(Cdb$classification=="")] <- paste( gsub(".*g__|;s__","",Cdb$classification_full[which(Cdb$classification=="")]),".species",sep="")





milk_taxonomic_profile_prevalence <- read_csv(file.path(DATA_DIR, "milk_taxonomic_profile_prevalence.csv"))

water_taxonomic_profile_prevalence <- read_csv(file.path(DATA_DIR, "water_taxonomic_profile_prevalence.csv"))


total_prevalence <- rbind(milk_taxonomic_profile_prevalence,water_taxonomic_profile_prevalence)



total_prevalence$...5[which(total_prevalence$...5=="Lactobacillus_ghanensis")] ="Liquorilactobacillus ghanensis"
total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 1" )] ="Lactococcus lactis" 
total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 2")] ="Lactococcus cremoris"

total_prevalence$...5[which(total_prevalence$...5=="Pseudomonas_fragi_subspecies 1")] ="Pseudomonas fragi"



total_prevalence$...5[which(total_prevalence$...5=="Zymomonas mobilis_subcluster 1")] ="Zymomonas mobilis"


total_prevalence$...5  = gsub("_", " ",total_prevalence$...5)

# Cdb <- 
#   Cdb[-c(which(Cdb$Completeness <80 |
#                  Cdb$Contamination >5)),]




"Saccharomyces cerevisiae"  
#[which(Cdb$category=="Milk.kefir")] 


cdb_prevalent <- rbind(
  
  
  
  
  Cdb[
    which(Cdb$classification%in%  total_prevalence$...5[which(total_prevalence$kefir_type=="milk")] & Cdb$category=="Milk.kefir"),],
  
  
  
  
  Cdb[
    which(Cdb$classification%in%  total_prevalence$...5[which(total_prevalence$kefir_type=="water")] & Cdb$category=="Water.kefir"),]
)


# this will get get all ani values to a refernce genome aligning with instrain and strainphlan comparisions
# ani_value_prevalent <- c()
# for (species in levels(as.factor(cdb_prevalent$classification))){
#   
#   
#   for (kefir in c("Milk.kefir","Water.kefir")){
#     
#     
#     t <- 
#     cdb_prevalent[which(cdb_prevalent$classification==species &
#                           cdb_prevalent$category==kefir),]
# 
#   

#   
# 
#     
#   if(length(levels(as.factor(t$primary_cluster)))!=1){
#     print(paste("something has gone wrong with", species,kefir))
#     next
#     
#   }
#     t1 <- Ndb[
#     which(
#       Ndb$reference==   wdb$genome[which(wdb$primary==levels(as.factor(t$primary_cluster)) &
#                                                    wdb$species_rep=="*")]),]
#     
#     
#     t1 <- t1[which(t1$querry %in% 
#     cdb_prevalent$user_genome[which(cdb_prevalent$primary_cluster==levels(as.factor(t$primary_cluster)) &
#                           cdb_prevalent$category==kefir)]),]
#     
#   
#   
#   ani_value_prevalent <- rbind(ani_value_prevalent,
#   data.frame(species=as.character(species),
#              type=as.character(kefir),
#              min=as.numeric(min(t1$ani)),
#              mean=as.numeric(mean(t1$ani)),
#              max=as.numeric(max(t1$ani)),
#                n_strains=as.numeric(length(unique(c(levels(as.factor(t1$querry)),levels(as.factor(t1$reference))))))
#   )
#   
#   )
#  
#   
#   
#   }
# }


#this will just compare all mags of the same species to each other

ani_value_prevalent <- c()

drep_prevelant_data <- c()
for (species in levels(as.factor(cdb_prevalent$classification))){
  
  
  for (kefir in c("Milk.kefir","Water.kefir")){
    
    
    t <- 
      cdb_prevalent[which(cdb_prevalent$classification==species &
                            cdb_prevalent$category==kefir),]
    
    
    
    
    
    
    if(length(levels(as.factor(t$primary_cluster)))!=1){
      print(paste("something has gone wrong with", species,kefir))
      next
      
    }
    
    t1 <-    Ndb[ which(Ndb$querry %in% t$user_genome &
                          Ndb$reference %in% t$user_genome
    ),]
    
    
    
    
    t1 <- t1[which(t1$querry %in% 
                     cdb_prevalent$user_genome[which(cdb_prevalent$primary_cluster==levels(as.factor(t$primary_cluster)) &
                                                       cdb_prevalent$category==kefir)]),]
    
    
    
    
    
    t1$species=species
    t1$type=kefir
    
    
    drep_prevelant_data <- rbind(drep_prevelant_data , t1)
    
    ani_value_prevalent <- rbind(ani_value_prevalent,
                                 data.frame(species=as.character(species),
                                            type=as.character(kefir),
                                            min=as.numeric(min(t1$ani)),
                                            mean=as.numeric(mean(t1$ani)),
                                            max=as.numeric(max(t1$ani)),
                                            n_strains=as.numeric(length(unique(c(levels(as.factor(t1$querry)),levels(as.factor(t1$reference))))))
                                 )
                                 
    )
    
    
    
  }
}


#ani_value_prevalent 

library(rstatix)



#jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/test_assembly.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)




# Create the plot
ggplot( ani_value_prevalent,
        aes(x = species, y = mean, ymin = min, ymax = max)) + #
  geom_pointrange() +
  geom_label(aes(y = max + 0.05 * (max - mean), label = n_strains), vjust = 0, size = 5, label.padding = unit(0.2, "lines")) +  # Adjust the position of n_strains label
  theme_minimal() +
  facet_wrap(~type,scales = "free")+
  theme_bw()+
  xlab("Species")+
  ylab("iRep Value")+
  labs(fill = "Timepoint")+
  #guides(colour = guide_legend(override.aes = list(size=25)))+
  #ylim(0,250)+
  theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
        legend.title = element_text( size=25, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
        axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2)#,
        # legend.position = "none"
  )



drep_prevelant_data <- 
merge( dplyr::select(ani_value_prevalent, species,n_strains),
       drep_prevelant_data,
       by="species",
       all.y=TRUE)



drep_plot <- 
  ggplot(drep_prevelant_data %>% 
           
           mutate(species=gsub("Lactococcus_lactis subcluster 1", "Lactococcus lactis",species )) %>% 
           mutate(species=gsub("Lactococcus_lactis subcluster 2", "Lactococcus cremoris",species )) %>% 
           mutate(species=gsub("Zymomonas_mobilis_subcluster 1", "Zymomonas mobilis",species )) %>% 
           mutate(species=gsub("Pseudomonas_fragi_subspecies 1", "Pseudomonas fragi",species )) %>% 
           mutate(species=gsub("_", " ",species )) %>% 
                  mutate(type=gsub("\\.", " ",type )
                  
           ),
         aes(x = species, y = ani)) + 
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  geom_boxplot()+
  #geom_label(aes(label = n_strains), vjust = -1, size = 5, label.padding = unit(0.2, "lines")) +  # Add text labels in boxes
  #geom_label(aes(y = max + 0.05 * (max - mean), label = n_strains), vjust = 0, size = 5, label.padding = unit(0.2, "lines")) +  # Adjust the position of n_strains label
  #geom_text(aes(label = n_strains), vjust = -1, size = 5) +  # Add text labels
  theme_minimal() +
  facet_wrap(~type, scales = "free") +
  theme_bw() +
  xlab("Species") +
  ylab("ANI") +
  labs(color = "Number of Strains") +
  
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15,face="italic"),
    strip.background = element_rect(
      color="black", fill="white"),
    strip.text = element_text(size=15.5),
    #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
    axis.title.y = element_text( size=20, face="bold",hjust = 0.5, vjust = 1.5),
    legend.key.size = unit(2, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm'), #change legend key width
    legend.direction = "horizontal", 
    legend.position = "top",
    legend.box = "horizontal",
    # legend.text = element_text(size=30),
    legend.title = element_text(size=17.5),
    legend.text=element_text(size=17.5))



  






library(ggpubr)

jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/supplementary figure 2.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)

ggarrange(strainphlan_plot,drep_plot,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(size = 30))
graphics.off()
