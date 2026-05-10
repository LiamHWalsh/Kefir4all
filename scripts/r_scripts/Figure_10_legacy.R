
###############################################################################################################
#Libraries used 
###############################################################################################################
#"E:/R/win-library/4.0")
.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")

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



###############################################################################################################
#Data importation trees
###############################################################################################################


setwd("Q:/H2020 Master/Citizen Science Project/Results/06_strain_profiling/06_strainphlan/06_strainphlan_sensitive/06_trees")
temp = list.files(pattern=".StrainPhlAn4.tre", recursive = FALSE)
myfiles = lapply(temp,read.tree)
###############################################################################################################
#Data modification
###############################################################################################################
# Data renaming 
names(myfiles) <- gsub("RAxML_bestTree.t__|.StrainPhlAn4.tre","",temp)

i <- "EUK4932"
###############################################################################################################
#Removing contaiminated samples

for (i in names(myfiles)){
  
  
  myfiles[[i]] <- drop.tip(myfiles[[i]], tip= grep("H3_S87_", myfiles[[i]]$tip.label))
  
  #myfiles[[i]]$tip.label <-gsub("_S.*","",myfiles[[i]]$tip.label)
}
###############################################################################################################
#Removing control samples

controls_to_remove <- c("A1","A2","A3","A4","A5","T6","T7","T8","T9","T10","U9","U10","U11","U12", "W11","W12","Y1","Y2","Y3")



control_1.1 <-c()
control_2.2 <- c()


y <- 1
#i <- "Lactococcus_lactis"
#i <- "A2"

for (i in controls_to_remove){ 
  
  control_1.1 <-  which(myfiles[[y]]$tip.label %in% controls_to_remove)
  
  myfiles[[y]] <- drop.tip(myfiles[[y]], tip= myfiles[[y]]$tip.label[c(control_1.1)])
  y <- y+1
  
  if (y==length(myfiles)){
    break
  }
  
}

########################################################################################################################
#Import metadata
########################################################################################################################

global_mk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_milk_kefir_metadata_v1.csv")
global_wk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_water_kefir_metadata_v1.csv")
global_mk_metadata$Stage <- NA
global_wk_metadata$Stage <- NA

Citizen_Scientist_metadata_v8 <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/Citizen Scientist metadata_v8.csv")

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/kefir4all_sample_metadata_v2.csv")
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


###############################################################################################################################
# drep data import 
###############################################################################################################################

# (C), which is generated using multidimensional scaling (MDS) of the ANI distance matrix.

# distance values drep
ndb <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/drep/data_tables/ndb.csv")
#ndb <- ndb[-c(which(ndb$querry=="Lactococcus_lactis_GCA_015476255.1_ASM1547625v1_genomic.fa")),]
mag_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/mag_metadata_v3.csv")


mag_metadata$bin <- 
  gsub(".*_bins_","",mag_metadata$user_genome)

mag_metadata$bin <- 
  gsub(".orig|.permissive|.strict","",mag_metadata$bin)


mag_metadata$bin <- 
  gsub(".*_","",mag_metadata$bin)


mag_metadata$bin_id <- paste(mag_metadata$base_name,"_",mag_metadata$bin,sep="")


#levels(as.factor(mag_metadata$bin ))


mag_metadata <- mag_metadata[-c(which(duplicated(mag_metadata$user_genome))),]

mag_metadata$classification_full <- mag_metadata$classification
mag_metadata$classification <- gsub(".*;s__","",mag_metadata$classification )

#mag_metadata <- mag_metadata[-c(which(mag_metadata$classification=="")) ,]





# &
# as.numeric(mag_metadata$fastani_ani)   <95),])


#mag_metadata[c(which(as.numeric(mag_metadata$)   <95),]))



mag_metadata <- dplyr::select(mag_metadata,"user_genome",	"base_name","bin_id",	"classification","classification_full","Sample.y","kefir type",	"Stage",	"data_source","type.x","Timepoint","category","Completeness",	"Contamination", "linker")

colnames(mag_metadata)[which(colnames(mag_metadata)=="Sample.y")] <- "Sample"
colnames(mag_metadata)[which(colnames(mag_metadata)=="type.x")] <- "type"






# leave only distance values
#closest_placement_ani
###############################################################################################################################
# drep data manipulation
###############################################################################################################################

ndb_plot <- dplyr::select(ndb,c("querry","reference","ani","alignment_coverage","querry_coverage","ref_coverage"))




#make into wide format making a distance matrix
ndb_plot$querry <- substr(ndb_plot$querry,1,nchar(ndb_plot$querry)-3)
ndb_plot$reference <- substr(ndb_plot$reference,1,nchar(ndb_plot$reference)-3)

#ndb_plot$querry <- gsub("-","_",ndb_plot$querry)
#ndb_plot$reference <- gsub("-","_",ndb_plot$reference)

#total_metadata$merge_column <- gsub("-","_",total_metadata$merge_column )

ndb_plot_v2 <-
  merge(ndb_plot,mag_metadata, by.x="querry",by.y="user_genome",all.x=TRUE)
ndb_plot_v2 <-
  merge(ndb_plot_v2,mag_metadata, by.x="reference",by.y="user_genome",all.x=TRUE)





###############################################################################################################################
# find same strains in the dataset
# 
#gonna to change this to highlight related strains 
###############################################################################################################################



# same_strains_drep <- 
# ndb_plot_v2[which(ndb_plot_v2$ani>=0.99 &
#                     ndb_plot_v2$alignment_coverage>=.90),]

related_strains_drep <- 
  ndb_plot_v2[which(ndb_plot_v2$ani>=0.99),]


###############################################################################################################################
# look at the putative new species ensure no mistake
###############################################################################################################################

#note .x will be querry and .y will be reference 
bins_of_interest <- 
  dplyr::select(related_strains_drep[which(related_strains_drep$classification.x==""),],querry,bin_id.x)

bins_of_interest <- 
  
  bins_of_interest[-c(which(duplicated(bins_of_interest$querry))),]


t <- c()
i="J1_S109_bin.4"

problem_bins <- c()
for (i in levels(as.factor(bins_of_interest$bin_id.x))){
  
  t <- mag_metadata[which(mag_metadata$bin_id %in%  i),] 
  
  if(nrow(t)==length(which(t$classification==""))){
    print(paste(i," is actually a putative spcies"))
  }else{
    print(paste(i," is not a putative spcies"))
    problem_bins <- c(problem_bins,bins_of_interest$querry[which(bins_of_interest$bin_id.x==i)])
  }
  
  
}
# REMOVE bins we are not sure of taxonomic origin
related_strains_drep <- 
  related_strains_drep[-c(which(related_strains_drep$reference %in% problem_bins)),]



related_strains_drep <- 
  related_strains_drep[-c(which(related_strains_drep$querry %in% problem_bins)),]


related_strains_drep$classification.x[which(related_strains_drep$classification.x=="")] <- paste("putative novel",gsub(".*;g__|;s__","",related_strains_drep$classification_full.x[which(related_strains_drep$classification.x=="")]),"species")


related_strains_drep$classification.y[which(related_strains_drep$classification.y=="")] <- paste("putative novel",gsub(".*;g__|;s__","",related_strains_drep$classification_full.y[which(related_strains_drep$classification.y=="")]),"species")

#of thoose putative species which seem to cluster with mags of a defined species

samples_with_species_to_change <- 
  dplyr::select(related_strains_drep[which(related_strains_drep$classification.x != related_strains_drep$classification.y &
                                             gsub("putative novel.*","",related_strains_drep$classification.x) == ""),],querry,classification.y)

# okay there were a couple update metadata accordingly leaving the full classifcation path alone to remind myself of update


#note .x will be querry and .y will be reference 

for (i in samples_with_species_to_change$querry){
  
  related_strains_drep$classification.x[which(related_strains_drep$querry==i)] <- samples_with_species_to_change$classification.y[which(samples_with_species_to_change$querry==i)]
  
  related_strains_drep$classification.y[which(related_strains_drep$reference==i)] <- samples_with_species_to_change$classification.y[which(samples_with_species_to_change$querry==i)]
  
  
}




#correlations with mutation rate and same strain recovery


###############################################################################################################################
# lOOK AT SAME VS DIFFERENT STRAINS frequecny 
###############################################################################################################################



if(nrow(
  related_strains_drep[which(related_strains_drep$classification.x != related_strains_drep$classification.y),])==0){
  print(paste("proceed no unusal results detected all mag group into the approprite species bin and are labelled to reflect this"))
}

# FIRST get number of comparsions per species in drep need to also compare  to total comparisions and remove self comparisions


related_strains_drep <- 
  related_strains_drep[-c(which(related_strains_drep$reference==related_strains_drep$querry)
  ),]

ndb_plot_v2 <- ndb_plot_v2[-c(which(ndb_plot_v2$reference==ndb_plot_v2$querry)),]



drep_species_sample_strain <- rbind(data.frame(species=as.character(levels(as.factor(related_strains_drep$classification.x))),
                                               number_of_comparisions=as.numeric(NA),
                                               total_comparisions=as.numeric(NA),
                                               type=as.character("Milk.kefir")), 
                                    data.frame(species=as.character(levels(as.factor(related_strains_drep$classification.x))),
                                               number_of_comparisions=as.numeric(NA),
                                               total_comparisions=as.numeric(NA),
                                               type=as.character("Water.kefir")))

t <- c()


ani_range <- c()
for (i in levels(as.factor(related_strains_drep$classification.x))){ # get the species detected 
  
  
  t <- related_strains_drep[which(related_strains_drep$classification.x==i &
                                    #related_strains_drep$classification.y==i &
                                    related_strains_drep$data_source.x=="This study" &
                                    related_strains_drep$data_source.y=="This study"  ),]
  
  
  
  for (type in  levels(as.factor(t$category.x ))){ # see if they are in milk or water kefir
    
    
    
    drep_species_sample_strain$number_of_comparisions[which(drep_species_sample_strain$species==i &
                                                              drep_species_sample_strain$type==type)]  <- 
      
      # this is a subset of the dataframe with the number of related strains comparisions within kefir4all
      nrow(t[which(
        t$category.x==type &
          t$category.y==type
      ),])
    
    
    
    # this is a subset of the dataframe with the number of total strains comparisions within kefir4all
    
    
    
    
    
    
    
    t1=  ndb_plot_v2[which((ndb_plot_v2$classification.x == i | ndb_plot_v2$classification.y == i) &
                             
                             
                             #ndb_plot_v2$classification.x==i &
                             #ndb_plot_v2$classification.y==i & # missing some classification y don't know why maybe a merge error anyway they are species specific comparisons so not needed here 
                             ndb_plot_v2$data_source.x=="This study" &
                             ndb_plot_v2$data_source.y=="This study" &
                             ndb_plot_v2$category.x==type &
                             ndb_plot_v2$category.y==type
    ),]
    
    
    
    
    
    
    
    drep_species_sample_strain$total_comparisions[which(drep_species_sample_strain$species==i &
                                                          drep_species_sample_strain$type==type)] <- 
      nrow(t1)
    
    
    
    ani_range <- rbind(ani_range,
                       data.frame(species=i,
                                  category=type,
                                  range_min=min(t1$ani),
                                  range_max=max(t1$ani)))
    
    
    # can remove leaving for trouble shooting only 
    # 
    #  ndb_plot_v2$merge <- paste(ndb_plot_v2$reference,ndb_plot_v2$querry,sep="_" )
    # 
    #  
    #  
    t$merge <- paste(t$reference,t$querry,sep="_")
    t1$merge <- paste(t1$reference,t1$querry,sep="_")
    # #   
    # # 
    #    t$merge[-c(which(t$merge %in% t1$merge))]
    t1[-c(which(t1$merge %in% t$merge)),]  
    # 
    #    
    #    
    #      ndb_plot_v2$classification.y[ which(ndb_plot_v2$merge==
    #    "TG_ID048_ML_T6_S55_metawrap_bin_reassembly_original_bins_bin.9_MK_TG_ID020_merged_T1_metawrap_bin_reassembly_original_bins_bin.8")]
    # 
    #    
    #    
    #   
    #   View(
    #     t[-c(which(t$merge %in% ndb_plot_v2$merge )),])
    #   
    
    
  }
}

drep_species_sample_strain$frequency <-   (drep_species_sample_strain$number_of_comparisions/  drep_species_sample_strain$total_comparisions)*100



drep_species_sample_strain$difference_frequency <- 100-drep_species_sample_strain$frequency

drep_species_sample_strain <-  drep_species_sample_strain[-c(which(is.na(drep_species_sample_strain$frequency=="Na"))),]



milk_taxonomic_profile_prevalence <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_metaphlan/prevalence/milk_taxonomic_profile_prevalence.csv")

water_taxonomic_profile_prevalence <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_metaphlan/prevalence/water_taxonomic_profile_prevalence.csv")


total_prevalence <- rbind(milk_taxonomic_profile_prevalence,water_taxonomic_profile_prevalence)





total_prevalence$...5[which(total_prevalence$...5=="Lactobacillus_ghanensis")] ="Liquorilactobacillus ghanensis"

total_prevalence$...5[which(total_prevalence$...5=="Lactobacillus_satsumensis")] ="Liquorilactobacillus satsumensis"

total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 1" )] ="Lactococcus lactis" 
total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 2")] ="Lactococcus cremoris"

total_prevalence$...5[which(total_prevalence$...5=="Pseudomonas_fragi_subspecies 1")] ="Pseudomonas fragi"



total_prevalence$...5[which(total_prevalence$...5=="Zymomonas_mobilis_subcluster 1")] ="Zymomonas mobilis"


total_prevalence$...5[which(total_prevalence$...5=="Lactobacillus_perolens")]="Schleiferilactobacillus perolens"



drep_species_sample_strain$species[which(drep_species_sample_strain$species=="Lactococcus lactis")] <- "Lactococcus lactis"
drep_species_sample_strain$species[which(drep_species_sample_strain$species=="Lactococcus lactis_E")] <- "Lactococcus cremoris"

drep_species_sample_strain$species[which(drep_species_sample_strain$species=="Gluconobacter oxydans_B")] <- "Gluconobacter oxydans"


total_prevalence$...5 <- 
  gsub("_"," ",total_prevalence$...5)


prevalent_cs_related_strain_frequency <- rbind(
  
  drep_species_sample_strain[which(drep_species_sample_strain$type=="Milk.kefir" &
                                     drep_species_sample_strain$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]),],
  
  drep_species_sample_strain[which(drep_species_sample_strain$type=="Water.kefir" &
                                     drep_species_sample_strain$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]),]
)

prevalent_cs_related_strain_frequency$ranking <- NA

prevalent_cs_related_strain_frequency$ranking[which(prevalent_cs_related_strain_frequency$type=="Milk.kefir")] <- 
  rank(-prevalent_cs_related_strain_frequency$frequency[which(prevalent_cs_related_strain_frequency$type=="Milk.kefir")])

prevalent_cs_related_strain_frequency$ranking[which(prevalent_cs_related_strain_frequency$type=="Water.kefir")] <- 
  rank(-prevalent_cs_related_strain_frequency$frequency[which(prevalent_cs_related_strain_frequency$type=="Water.kefir")])





prevalent_cs_related_strain_frequency$difference_frequency <- 100-prevalent_cs_related_strain_frequency$frequency



pa <- 
  rbind(
    
    prevalent_cs_related_strain_frequency %>% dplyr::select(species,frequency,type,ranking) %>% mutate(strain_type="% highly related strain") ,
    
    
    prevalent_cs_related_strain_frequency %>% dplyr::select(species,difference_frequency,type,ranking) %>% mutate(strain_type="% less related strains:") %>% dplyr::rename("frequency"="difference_frequency")
  )  %>% #pivot_wider(names_from="species", values_from="frequency", values_fill = 0) %>% 
  
  #
  
  
  
  
  ggplot(aes(x=reorder(species,ranking),y=frequency,fill=factor(strain_type),order=factor(strain_type)))+
  geom_col()+
  geom_text(aes(y=frequency,color=factor(strain_type), label = round(frequency,1)),vjust=-0.25,hjust=.5,size=3,fontface='bold',show.legend=FALSE, position = position_stack(vjust = 0.5))+
  labs(x="Species", y="Strain comparisons", title="",fill="Strain comparisons") +
  #coord_equal() +
  theme_bw()+
  facet_wrap(~type,scales="free")+
  theme(legend.position = "top",#axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.x = element_text(size=15),
        axis.title = element_text(size = 20),
        # axis.text.x = element_text(size=15,face = "italic",angle=90,hjust = .5,vjust = .6),
        axis.text.y = element_text(size=15,face = "italic"),
        #axis.text.x = element_text(size=10,angle = 45,hjust = 1), # remove x-axis labels
        #axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        # position = position_stack(vjust = 0.5),
        legend.text=element_text(face="italic",size = 20),
        legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(1.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'), #change legend key width
        legend.title = element_text(size=20),
        strip.background = element_rect(
          color="black", fill="white"),
        strip.text = element_text(size=15.5))+
  scale_colour_discrete(l = 40)+
  coord_flip()+
  guides(fill = guide_legend(nrow = 1))



########################################################################################################################
#Import metadata
########################################################################################################################

global_mk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_milk_kefir_metadata_v1.csv")
global_wk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_water_kefir_metadata_v1.csv")
global_mk_metadata$Stage <- NA
global_wk_metadata$Stage <- NA

Citizen_Scientist_metadata_v8 <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/Citizen Scientist metadata_v8.csv")

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/kefir4all_sample_metadata_v2.csv")
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




# distance values drep
ndb <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/drep/data_tables/ndb.csv")


Cdb <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/drep/data_tables/Cdb.csv")


#ndb <- ndb[-c(which(ndb$querry=="Lactococcus_lactis_GCA_015476255.1_ASM1547625v1_genomic.fa")),]
mag_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/mag_metadata_v3.csv")





mag_metadata$bin <- 
  gsub(".*_bins_","",mag_metadata$user_genome)

mag_metadata$bin <- 
  gsub(".orig|.permissive|.strict","",mag_metadata$bin)


mag_metadata$bin <- 
  gsub(".*_","",mag_metadata$bin)


mag_metadata$bin_id <- paste(mag_metadata$base_name,"_",mag_metadata$bin,sep="")


#levels(as.factor(mag_metadata$bin ))


mag_metadata <- mag_metadata[-c(which(duplicated(mag_metadata$user_genome))),]

mag_metadata$classification_full <- mag_metadata$classification
mag_metadata$classification <- gsub(".*;s__","",mag_metadata$classification )

#mag_metadata <- mag_metadata[-c(which(mag_metadata$classification=="")) ,]





# &
# as.numeric(mag_metadata$fastani_ani)   <95),])


#mag_metadata[c(which(as.numeric(mag_metadata$)   <95),]))



mag_metadata <- dplyr::select(mag_metadata,"user_genome",	"base_name","bin_id",	"classification","classification_full","Sample.y","kefir type",	"Stage",	"data_source","type.x","Timepoint","category","Completeness",	"Contamination", "linker")

colnames(mag_metadata)[which(colnames(mag_metadata)=="Sample.y")] <- "Sample"
colnames(mag_metadata)[which(colnames(mag_metadata)=="type.x")] <- "type"






# leave only distance values
#closest_placement_ani
###############################################################################################################################
# drep data manipulation
###############################################################################################################################

# import relevant data
Cdb$genome <- gsub(".fa","",Cdb$genome)

nrow( Cdb)
Cdb <-  merge(mag_metadata,
              Cdb,
              by.x="user_genome",
              by.y="genome",
              all.y=TRUE)




Cdb[which(Cdb$classification==""),] <- paste( gsub(".*g__|;s__","",Cdb$classification_full[which(Cdb$classification=="")]),".species",sep="")


# identify problem samples
Cdb[which(is.na(Cdb$Completeness)),]

# get the range of completness and contamination of MAGs
range(Cdb$Completeness[-c(which(is.na(Cdb$Completeness)))])

range(Cdb$Contamination[-c(which(is.na(Cdb$Contamination)))])



# get the frequency of primary clusters
drep_n_primary <- 
  as.data.frame(table(Cdb$primary_cluster))


# get the frequency of secondary clusters
drep_n_secondary <- 
  as.data.frame(table(Cdb$secondary_cluster))

nrow(drep_n_primary)

nrow(
  drep_n_primary[which(drep_n_primary$Freq>=10),])


nrow(drep_n_secondary)
nrow(
  drep_n_secondary[which(drep_n_secondary$Freq>=10),])



# filter to just look at >=10 primary clusters 

Cdb_top=
  
  Cdb[which(Cdb$primary_cluster %in%    drep_n_primary$Var1[which(drep_n_primary$Freq>=10)]),]




i="d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lentilactobacillus;s__Lentilactobacillus hilgardii"


i="d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus lactis"

i="d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus raffinolactis" 
i="d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus helveticus" 

T0_mag_rep <- c()

t_time <- c()

kefir4all_secondary_clusters <-c()
species_mag_divering <- c()


secondary_cluster <- c()

kefir4all_secondary_cluster <- c()
problem_samples_cluster <- c()


for (category in c("Milk.kefir","Water.kefir")){
  for (i in levels(as.factor(Cdb_top$classification_full))){
    
    # subset based on variables of interest
    t= Cdb_top[which(Cdb_top$classification_full==i &
                       Cdb_top$category==category),]
    
    
    
    
    # skip features without variables of interest
    
    if(nrow(t)==0){
      print(paste("No entries for ",  i,"-",category, " secondary_cluster analysis, skipping"))
      next
      
    }
    
    t_cs=t[which(t$data_source=="This study"),]
    
    
    
    if(nrow(t_cs)==0){
      print(paste("No entries for kefir4all -  ",  i,"-",category, " secondary_cluster analysis, skipping"))
      next
      
    }
    
    
    
    t1=(as.data.frame(
      xtabs(~data_source+secondary_cluster,t)))
    
    
    t1$category=category
    
    t1$classification=i
    
    
    
    secondary_cluster <- rbind(secondary_cluster,t1)
    
    
    t_time <-as.data.frame(
      xtabs(~data_source+secondary_cluster+Stage,t[which(t$data_source=="This study"),] ))
    
    
    t_time$category=category
    
    t_time$classification=i
    
    
    
    kefir4all_secondary_cluster <- rbind(kefir4all_secondary_cluster, t_time)
    
    
    # ggplot(t_time, aes(x = Stage, y = Freq, fill = secondary_cluster)) +
    #   geom_bar(stat = "identity", position = "dodge") +
    #   labs(title = "Frequency of Clusters Over Different Stages",
    #        x = "Stage", y = "Frequency",
    #        fill = "Secondary Cluster") +
    #   theme_minimal()
    # 
    
    
    print(paste("for classification", i, "we found", length(levels(as.factor(t1$secondary_cluster))),"secondary clusters across all datasets for",category ))
    
    
    
    if(length(
      
      levels(as.factor(t_cs$secondary_cluster[which(t_cs$Stage==levels(as.factor(t_cs$Stage))[1])])))==1){
      
      
      print(paste(
        "the representative genome at baseline for", i, "is", levels(as.factor(
          t$secondary_cluster[which(t$Stage==  levels(as.factor(t$Stage))[1]  
          )]))))
      
      
      T0_mag_rep <- rbind( T0_mag_rep,
                           data.frame(classification=as.character(i),
                                      category=as.character(category),
                                      rep_method=as.character(levels(as.factor(t$Stage))[1]),
                                      rep=as.character(levels(as.factor(t$secondary_cluster[which(t$Stage==levels(as.factor(t$Stage))[1])]))))
                           
                           
      )
      
      
      
      
      
      
    }else{
      
      print(paste("error occured for ",  i,"-",category, " secondary_cluster for T0 analysis, saving to problem_samples_cluster"))
      problem_samples_cluster <- rbind(  problem_samples_cluster ,
                                         data.frame(classification=as.character(i),
                                                    category=as.character(category),
                                                    rep_method=as.character(levels(as.factor(t$Stage))[1]))
      )
      
    }
    
    
    
    
    
    
    
    
    
    
  }#end of i
  
} # end of category












kefir4all_secondary_cluster$clusters <- gsub(".*_","",kefir4all_secondary_cluster$secondary_cluster)


#number of MAGs within secondary clusters per classifcation and category
kefir4all_secondary_cluster_classification_sum <- c()


#number of secondary clusters per classifcation and category
kefir4all_secondary_cluster_per_classification <- c()





i="d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lacticaseibacillus;s__Lacticaseibacillus paracasei" 

for (category in c("Milk.kefir","Water.kefir")){
  for (i in levels(as.factor(kefir4all_secondary_cluster$classification))){
    
    
    t <-  as.data.frame(kefir4all_secondary_cluster[which(kefir4all_secondary_cluster$classification==i &
                                                            kefir4all_secondary_cluster$category==category
    ),])
    
    
    kefir4all_secondary_cluster_classification_sum <- rbind(kefir4all_secondary_cluster_classification_sum,
                                                            data.frame(classification=as.character(i),
                                                                       category=as.character(category),
                                                                       sum_freq=as.numeric( sum(t$Freq)))
    )
    
    kefir4all_secondary_cluster_per_classification <- rbind(kefir4all_secondary_cluster_per_classification , 
                                                            data.frame(classification=as.character(i),
                                                                       category=as.character(category),
                                                                       n_secondary_clusters=as.numeric(length(levels(as.factor(as.character(t$secondary_cluster))))))
                                                            
                                                            #
    )
    
    
    
  }
  
  
}










kefir4all_secondary_cluster_per_classification$merge <-  paste(kefir4all_secondary_cluster_per_classification$classification, kefir4all_secondary_cluster_per_classification$category,sep="_")








input_plot_data <- 
  
  
  rbind(
    kefir4all_secondary_cluster_per_classification[which(
      kefir4all_secondary_cluster_per_classification$n_secondary_clusters>1 &
        
        
        kefir4all_secondary_cluster_per_classification$classification %in% kefir4all_secondary_cluster_classification_sum$classification[which(kefir4all_secondary_cluster_classification_sum$sum_freq>20 &
                                                                                                                                                 kefir4all_secondary_cluster_classification_sum$category=="Milk.kefir"   
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
        )] &
        kefir4all_secondary_cluster_per_classification$category=="Milk.kefir"                               
      
    ),]
    , 
    kefir4all_secondary_cluster_per_classification[which(
      kefir4all_secondary_cluster_per_classification$n_secondary_clusters>1 &
        
        
        kefir4all_secondary_cluster_per_classification$classification %in% kefir4all_secondary_cluster_classification_sum$classification[which(kefir4all_secondary_cluster_classification_sum$sum_freq>20 &
                                                                                                                                                 kefir4all_secondary_cluster_classification_sum$category=="Water.kefir"   
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
        )] &
        kefir4all_secondary_cluster_per_classification$category=="Water.kefir"                               
      
    ),]
  )



library(ggsankey)

library(networkD3)


library(ggplot2)
library(ggforce)
library(dplyr)
pacman::p_load(ggalluvial)

kefir4all_secondary_cluster$classification[
  which(kefir4all_secondary_cluster$classification==
          "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus lactis_E")] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus cremoris"



kefir4all_secondary_cluster$classification[
  which(kefir4all_secondary_cluster$classification==
          "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Pseudomonas_E fragi_B")] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Pseudomonas fragi"



kefir4all_secondary_cluster$classification[
  which(kefir4all_secondary_cluster$classification==
          "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Gluconobacter oxydans_B")] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Gluconobacter oxydans"




kefir4all_secondary_cluster$merge <- paste(kefir4all_secondary_cluster$classification, kefir4all_secondary_cluster$category,sep="_")


kefir4all_secondary_cluster$classification_v2=sub(".*s__","",kefir4all_secondary_cluster$classification )






kefir4all_secondary_cluster_prevalent<- 
  rbind(
    kefir4all_secondary_cluster[which(kefir4all_secondary_cluster$category=="Milk.kefir" &
                                        kefir4all_secondary_cluster$classification_v2 %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]),],
    
    kefir4all_secondary_cluster[which(kefir4all_secondary_cluster$category=="Water.kefir" &
                                        kefir4all_secondary_cluster$classification_v2 %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]),]
    
  )


total_prevalence$...5[-c(
which(total_prevalence$...5 %in% kefir4all_secondary_cluster_prevalent$classification_v2))]



kefir4all_secondary_cluster$classification_v2 %in% total_prevalence$...5



pb <- 
  kefir4all_secondary_cluster_prevalent %>% 
  # kefir4all_secondary_cluster[which(
  #   kefir4all_secondary_cluster$merge %in% 
  #     input_plot_data$merge),] %>%   
  mutate(classification= gsub(".*s__","",classification )) %>% 
  mutate(category= gsub("\\."," ",category )) %>% 
  dplyr::select(Stage, Freq, clusters,classification,category) %>%
  distinct() %>%
  
  
  
  mutate(  Stage =  gsub("T1", "wk01",
                         gsub( "T2", "wk05",
                               gsub( "T3", "wk09",
                                     gsub("T4", "wk13",
                                          gsub("T5", "wk17",
                                               gsub( "T6", "wk21", Stage)))))),
           Stage = factor(Stage, levels = c("T0", "wk01", "wk05", "wk09", "wk13", "wk17","wk21")),
           clusters=factor(clusters,levels=c(1:15))
           
  )%>% 
  #mutate(classification= gsub(".*s__","",classification )) %>% 
  
  
  
  
  ggplot( aes(x = Stage, stratum = clusters, alluvium = clusters, y = Freq, fill = clusters)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum() +
  # scale_fill_brewer(type = "qual")+#, palette = "Set1") +
  theme_bw() +
  labs(title = "",
       x = "Timeframes",
       y = "Number of metagenome assembled genomes",
       fill="Clusters")+
  
  facet_wrap(~category+classification,scales = "free_y")+
  theme(legend.position = "top",#axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=15,,hjust = .5,vjust = .6),
        #axis.text.x = element_text(size=10,angle = 45,hjust = 1), # remove x-axis labels
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
        strip.background = element_rect(
          color="black", fill="white"),
        strip.text = element_text(size=15.5,face="italic"))+
  scale_colour_discrete(l = 40)+
    guides(fill = guide_legend(nrow = 1))




pb

###############################################################################################################
#look into the genetic distance between clusters
###############################################################################################################
.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")

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



instrain =read_csv(
  "Q:/H2020 Master/Citizen Science Project/Results/06_strain_profiling/06_instrain/combined_outputs/instrain_genome_species_primary_data_V4.csv"
)



library(readr)

all_genomes_strain_v2 <- read_delim("Q:/H2020 Master/Citizen Science Project/Results/06_strain_profiling/06_instrain/all.genomes_strain_v2.stb", 
                                    delim = "\t", escape_double = FALSE, 
                                    col_names = FALSE, trim_ws = TRUE)

all_genomes_strain_v2$X2 <- substr(all_genomes_strain_v2$X2 ,1,nchar(all_genomes_strain_v2$X2 )-3)



instrain_gene_total_data <- c()
scaffolds <- c()
t2 <- c()


########################################################################################################################
#
########################################################################################################################



########################################################################################################################
#Import sample metadata
########################################################################################################################

global_mk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_milk_kefir_metadata_v1.csv")
global_wk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_water_kefir_metadata_v1.csv")
global_mk_metadata$Stage <- NA
global_wk_metadata$Stage <- NA

Citizen_Scientist_metadata_v8 <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/Citizen Scientist metadata_v8.csv")

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/kefir4all_sample_metadata_v2.csv")
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
# MImport prevalence metadata
########################################################################################################################

milk_taxonomic_profile_prevalence <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_metaphlan/prevalence/milk_taxonomic_profile_prevalence.csv")

water_taxonomic_profile_prevalence <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_metaphlan/prevalence/water_taxonomic_profile_prevalence.csv")


total_prevalence <- rbind(milk_taxonomic_profile_prevalence,water_taxonomic_profile_prevalence)



total_prevalence$...5[which(total_prevalence$...5=="Lactobacillus_ghanensis")] ="Liquorilactobacillus ghanensis"

total_prevalence$...5[which(total_prevalence$...5=="Lactobacillus_satsumensis")] ="Liquorilactobacillus satsumensis"

total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 1" )] ="Lactococcus lactis" 
total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 2")] ="Lactococcus cremoris"

total_prevalence$...5[which(total_prevalence$...5=="Pseudomonas_fragi_subspecies 1")] ="Pseudomonas fragi"



total_prevalence$...5[which(total_prevalence$...5=="Zymomonas_mobilis_subcluster 1")] ="Zymomonas mobilis"


total_prevalence$...5[which(total_prevalence$...5=="Lactobacillus_perolens")]="Schleiferilactobacillus perolens"



total_prevalence$...5 <- 
  gsub("_"," ",total_prevalence$...5)


# 
# prevalence_metric_data <- c()
# for (unit in names(metric_data)){
#   
#   prevalence_metric_data[[unit]]<-  rbind(metric_data[[unit]][which( metric_data[[unit]]$kefir_type=="Milk.kefir" &
#                                                                        metric_data[[unit]]$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]
#                                                                      
#   ),],
#   metric_data[[unit]][which( metric_data[[unit]]$kefir_type=="Water.kefir" &
#                                metric_data[[unit]]$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]
#   ),]
#   )
#   
#   
#   
# }




########################################################################################################################

########################################################################################################################


species= "Lactobacillus helveticus"  



instrain<- merge(instrain, total_metadata,by.x="sample_id",by.y="merge_column",all.x=TRUE)




clus_breakdown <- c()






  species="Lactococcus cremoris"
  
  for (species in  levels(as.factor(instrain$classification))){
    
    for (type in c("Milk.kefir", "Water.kefir")){
    
    t1 <- instrain[which(instrain$classification==species &
                           instrain$category.y==type&
                           instrain$popANI_reference>.98),]
    
    # 
    # t2 <- 
    # 
    # t1[which(t1$sample_id %in% kefir4all_metadata$merge_column),]
    # 
    # #
    # if(nrow(t2)==0){
    #   
    #   print(paste("Did not identify ",species, " in any cs ", type," metagenomes",sep=""))
    #   next
    # }
    # 
    
    
    #dplyr::select( t2, sample_id,genome, Stage.y,Stage.x)
    
    
    
    if(nrow(t1)==0){next }else{
      
      t2 <- data.frame( table(t1$cluster))
      
      
      clus_breakdown <- 
        rbind(clus_breakdown, data.frame(species=species,
                                         type=type, 
                                         clust=t2$Var1, 
                                         Freq=t2$Freq
        ))
    }
  }
  
  }
  
  
  library(upstartr)

clus_breakdown_prevalent <- 
  rbind(
    clus_breakdown[which(clus_breakdown$type=="Milk.kefir" &
                           clus_breakdown$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]),],
    clus_breakdown[which(clus_breakdown$type=="Water.kefir" &
                           clus_breakdown$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]),]) 



clus_breakdown_prevalent %>% 
  mutate(clust=gsub(".*_","",clust)) %>% 
  
  
  #ggplot(aes(x=  reorder(...2, count),y=count,fill=detection_category))+
  ggplot( aes(x = fct_reorder(species, Freq, .fun = sum),y=Freq,fill=clust))+
  geom_col()+
  facet_wrap(~type,scales="free")+
  scale_x_reordered()+
  # geom_boxplot()+
  #geom_point()+
  labs(x="Species", y="Number of strains detected", title="",fill="Secondary_clusters") +
  #coord_equal() +
  theme_bw()+
  theme(legend.position = "top",#axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.x = element_text(size=15),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size=12.5,face = "italic",hjust = .5,vjust = .6),
        #axis.text.x = element_text(size=10,angle = 45,hjust = 1), # remove x-axis labels
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
        strip.background = element_rect(
          color="black", fill="white"),
        strip.text = element_text(size=15.5))+
  guides(fill = guide_legend(nrow = 2))+
  coord_flip()












instrain_prevelant_cs_clust <- 
instrain[which(
                 instrain$popANI_reference>.98 &
                   
                 instrain$sample_id %in% kefir4all_metadata$merge_column),]




instrain_prevelant_cs_clust <- 
  rbind(
    instrain_prevelant_cs_clust[which(instrain_prevelant_cs_clust$category.y=="Milk.kefir" &
                                        instrain_prevelant_cs_clust$classification %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]),],
    
    instrain_prevelant_cs_clust[which(instrain_prevelant_cs_clust$category.y=="Water.kefir" &
                                        instrain_prevelant_cs_clust$classification %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]),]
  )
    
    
    

library(ggalluvial)


  
instrain_prevelant_cs_clust <- 
  instrain_prevelant_cs_clust %>% 
  mutate(clusters=gsub(".*_","",cluster)) %>% 

  mutate(  Stage =  gsub("T1", "wk01",
                         gsub( "T2", "wk05",
                               gsub( "T3", "wk09",
                                     gsub("T4", "wk13",
                                          gsub("T5", "wk17",
                                               gsub( "T6", "wk21", Stage.x)))))),
           Stage = factor(Stage, levels = c("T0", "wk01", "wk05", "wk09", "wk13", "wk17","wk21")))
       
  
  
pc<- 
as.data.frame(
xtabs(~Stage+clusters+ classification+category.y,instrain_prevelant_cs_clust )) %>% 
  filter(Freq !=0) %>% 
  
  mutate(
           Stage = factor(Stage, levels = c("T0", "wk01", "wk05", "wk09", "wk13", "wk17","wk21")),
           clusters=factor(clusters,levels=c(1:20))
           
  )%>% 
  mutate(category.y=gsub("\\."," ",category.y)) %>% 
  
  ggplot( aes(x = Stage, stratum = clusters, alluvium = clusters, y = Freq, fill = clusters)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum() +
  # scale_fill_brewer(type = "qual")+#, palette = "Set1") +
  theme_bw() +
  labs(title = "",
       x = "Timeframes",
       y = "Number of strains detected",
  fill="Clusters")+
  facet_wrap(~category.y+classification,scales = "free_y")+
  theme(legend.position = "right",#axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=15,,hjust = .5,vjust = .6),
        #axis.text.x = element_text(size=10,angle = 45,hjust = 1), # remove x-axis labels
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
        strip.background = element_rect(
          color="black", fill="white"),
        strip.text = element_text(size=15.5,face="italic"))+
  scale_colour_discrete(l = 40)+
  guides(fill = guide_legend(nrow = 1))




#########################################


#########################################

jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 13_v2.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)



#ggarrange( pc,pa,nrow=2,labels=c("A.","B."),font.label = list(size = 30),heights=c(1.4), common.legend = FALSE)


#ggarrange( pb, pc,pa,nrow=3,labels=c("A.","B.","C."),font.label = list(size = 30),heights=c(1,1,.5), common.legend = FALSE)



#ggarrange(
ggarrange( pb, pc,labels=c("A.","B."),font.label = list(size = 30),heights=c(1,1), common.legend = TRUE,nrow=2)+
  guides(colour = guide_legend(nrow = 1))#,
           
  #          pa,nrow=2,labels=c("","C."),font.label = list(size = 30),heights=c(1,.3), common.legend = FALSE)+
  # theme(legend.position = "top")+
  # guides(fill = guide_legend(nrow = 1))



graphics.off()



jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/supplementary figure 3.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)


pa

graphics.off()




# jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/instrain_strain_clusters_cs.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
# 
# 
# plot(plot)
# graphics.off()


# track nucl_diversity or pop ani across one strain sub clust over the stages 







########################################################################################################################
# how releated our strain using instrain 
########################################################################################################################



total_prevalence$kefir_type=
gsub("milk","Milk.kefir",
     gsub("water","Water.kefir",total_prevalence$kefir_type))

for (species in levels(as.factor(total_prevalence$...5))){
  
  t1= instrain[ which(instrain$data_source.y=="This study" &
                        instrain$classification==species&
                   instrain$category.y==
                            total_prevalence$kefir_type[which(total_prevalence$...5==species)]), ]




for (species in  levels(as.factor(instrain$classification))){
  
  for (type in c("Milk.kefir", "Water.kefir")){
t <- 
instrain_prevelant_cs_clust[which(instrain_prevelant_cs_clust$classification==species &
                                    instrain_prevelant_cs_clust$category.x==type),]

t1 <- 
dplyr::select(t, sample_id,genome,Sample.x,clusters, Stage,category.y,iRep)


t2 <- 
as.data.frame(xtabs(~Sample.x+clusters,t1)) %>% filter(Freq!=0)


names(which(table(t2$Sample.x)>1))




t3 <- rbind( t1[which(t1$Stage=="T0"),],
  
  t1[which(t1$Sample.x %in% names(which(table(t2$Sample.x)>1))),])


#fwrite( t3,paste("Q:/H2020 Master/Citizen Science Project/Results/06_strain_profiling/06_instrain/instrain_compare/instrain_compare_selected_sample_",species,"_",type,".csv",sep="" ))



  }
  
}










for (species in levels(as.factor(prevelance_data[[type]]$...5))){
  
  
  
  t1 <- instrain[which(instrain$classification==species &
                         instrain$category.y==type&
                         instrain$popANI_reference>.98),]
  
  
  if(length(levels(as.factor(    t1$cluster)))>1){
    
    
    
    t2 <-   
      t1 %>%
      group_by(cluster) %>%
      top_n(5, wt = popANI_reference) %>%
      arrange(cluster, desc(popANI_reference))
    
    
    samples_with_multiple_clus[[type]] <- 
      rbind(samples_with_multiple_clus[[type]], 
            t2)
    
    
    
  }
}

samples_with_multiple_clus[[type]] <- rbind( samples_with_multiple_clus[[type]], 
                                             t_baseline
)
























species="Lactococcus lactis"

type="Milk.kefir"


samples_with_multiple_clus <- c()

prevelance_data <- c()
prevelance_data[["Milk.kefir"]] <- 
milk_taxonomic_profile_prevalence

prevelance_data[["Water.kefir"]] <- 
water_taxonomic_profile_prevalence

dir.create("Q:/H2020 Master/Citizen Science Project/Results/06_strain_profiling/06_instrain/instrain_compare")

library(data.table)
for (type in c("Milk.kefir", "Water.kefir")){
  
  
  
  prevelance_data[[type]]$...5[which(prevelance_data[[type]]$...5=="Lactobacillus_ghanensis")] ="Liquorilactobacillus ghanensis"
 prevelance_data[[type]]$...5[which(prevelance_data[[type]]$...5=="Lactococcus_lactis subcluster 1" )] ="Lactococcus lactis" 
 prevelance_data[[type]]$...5[which(prevelance_data[[type]]$...5=="Lactococcus_lactis subcluster 2")] ="Lactococcus cremoris"
  
 prevelance_data[[type]]$...5[which(prevelance_data[[type]]$...5=="Zymomonas_mobilis_subcluster 1")] ="Zymomonas mobilis"
  
  
 t_baseline=instrain[which(
                             instrain$category.y==type&
                           
                                            instrain$Stage.x=="T0"),]
   

 
  for (species in levels(as.factor(prevelance_data[[type]]$...5))){
    
    
    
    t1 <- instrain[which(instrain$classification==species &
                           instrain$category.y==type&
                           instrain$popANI_reference>.98),]
  
        
    if(length(levels(as.factor(    t1$cluster)))>1){
      
    
      
    t2 <-   
      t1 %>%
        group_by(cluster) %>%
        top_n(5, wt = popANI_reference) %>%
        arrange(cluster, desc(popANI_reference))
      
      
    samples_with_multiple_clus[[type]] <- 
    rbind(samples_with_multiple_clus[[type]], 
          t2)
          

      
    }
  }
  
 samples_with_multiple_clus[[type]] <- rbind( samples_with_multiple_clus[[type]], 
                                              t_baseline
                                              )
  
 fwrite( samples_with_multiple_clus[[type]],paste("Q:/H2020 Master/Citizen Science Project/Results/06_strain_profiling/06_instrain/instrain_compare/instrain_compare_selected_sample_",type,".csv",sep="" ))
  levels(as.factor(samples_with_multiple_clus[[type]]$sample_id))
}
        
   






     t2 <-
          
          t1[which(t1$sample_id %in% kefir4all_metadata$merge_column),]



        
        View(t2)

for (analysis in c("secondary", "primary")){
  
  
  if (analysis=="primary"){
    
    t1 <-   instrain_primary [
      which(
        instrain_primary$classification == species & 
          
          !is.na( getElement(instrain_primary, unit))&
          
          instrain_primary$category.y==type
        
        
      ),]
    
    
    
    
    if(nrow(t1)==0){
      next    }
    
    metric_data[["total_clust"]]= rbind(metric_data[["total_clust"]], 
                                        data.frame(species=species,
                                                   cluster="total_clust",
                                                   unit=as.character(unit),
                                                   min=range(getElement(t1,unit))[1], 
                                                   max=  range(getElement(t1,unit))[2],
                                                   mean= mean(getElement(t1,unit)), 
                                                   n=as.numeric(length(unique(t1$sample_id))), 
                                                   
                                                   kefir_type=type
                                                   
                                                   
                                                   
                                        )
    )
    
    
    for (clust in  levels(as.factor(t1$cluster))){
      
      
      metric_data[["clusters"]]= rbind( metric_data[["clusters"]], 
                                        data.frame(species=species,
                                                   cluster=clust,
                                                   
                                                   unit=as.character(unit),
                                                   min=range(getElement(t1[which(t1$cluster==clust),],unit))[1], 
                                                   max=  range(getElement(t1[which(t1$cluster==clust),],unit))[2],
                                                   mean= mean(getElement(t1[which(t1$cluster==clust),],unit)), 
                                                   n=as.numeric(length(unique(t1$sample_id[which(t1$cluster==clust)]))), 
                                                   
                                                   kefir_type=type
                                        )
                                        
      )
      
      
      
      
    }
    
  }else if (analysis=="secondary"){
    t1 <- 
      
      instrain_genome_total_data[
        which(
          instrain_genome_total_data$classification == species & 
            
            !is.na( getElement(instrain_genome_total_data, unit))&
            
            instrain_genome_total_data$category.y==type
          
          
        ),]
    
    if(nrow(t1)==0){
      next    }
    
    metric_data[["total"]]= rbind(metric_data[["total"]], 
                                  data.frame(species=species,
                                             cluster="total",
                                             unit=as.character(unit),
                                             min=range(getElement(t1,unit))[1], 
                                             max=  range(getElement(t1,unit))[2],
                                             mean= mean(getElement(t1,unit)), 
                                             n=as.numeric(length(unique(t1$sample_id))), 
                                             
                                             kefir_type=type
                                             
                                             
                                             
                                  )
    )
  }
  
  
  #  
  # which( instrain_genome_total_data$classification==species &
  #    instrain_genome_total_data$iRep!=NA)
  
}


}

}

}


prevalence_metric_data <- c()
for (unit in names(metric_data)){
  
  prevalence_metric_data[[unit]]<-  rbind(metric_data[[unit]][which( metric_data[[unit]]$kefir_type=="Milk.kefir" &
                                                                       metric_data[[unit]]$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]
                                                                     
  ),],
  metric_data[[unit]][which( metric_data[[unit]]$kefir_type=="Water.kefir" &
                               metric_data[[unit]]$species %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]
  ),]
  )
  
  
  
}




#"179_1" weird primary cluster of lc.raffinolactis seperate to the main primary cluster, I suggest removing
prevalence_metric_data[[ "clusters" ]] <- 
  prevalence_metric_data[[ "clusters" ]][
    -c(which(prevalence_metric_data[[ "clusters" ]]$cluster=="179_1")),]


library(rstatix)



plot <- 
  # Create the plot
  ggplot( prevalence_metric_data[["total"]]%>% mutate(cluster=gsub(".*_","",cluster)), aes(x = species, y = mean, ymin = min, ymax = max, color = cluster)) +
  geom_pointrange() +
  theme_minimal() +
  facet_wrap(~kefir_type+unit,scales = "free")+
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





