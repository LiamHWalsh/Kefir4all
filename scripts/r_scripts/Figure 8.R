

.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")
pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis)
pacman::p_load(readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)
library(vegan)
library(ggplot2)
library(grid)
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
#ndb <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/drep/data_tables/ndb.csv")




#ndb <- ndb[-c(which(ndb$querry=="Lactococcus_lactis_GCA_015476255.1_ASM1547625v1_genomic.fa")),]
mag_metadata_pro <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/HQ_prokaryotic_representatives_MAGs_per_sample.csv")


mag_metadata_euk <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/HQ_eukaryotic_MAGs.csv")


mag_metadata_pro$bin <- 
  gsub(".*_bins_","",mag_metadata_pro$user_genome)

mag_metadata_pro$bin <- 
  gsub(".orig|.permissive|.strict","",mag_metadata_pro$bin)


mag_metadata_pro$bin <- 
  gsub(".*_","",mag_metadata_pro$bin)


mag_metadata_pro$bin_id <- paste(mag_metadata_pro$base_name,"_",mag_metadata_pro$bin,sep="")


#levels(as.factor(mag_metadata$bin ))


#mag_metadata_pro <- mag_metadata_pro[-c(which(duplicated(mag_metadata_pro$user_genome))),]

mag_metadata_pro$classification_full <- mag_metadata_pro$classification
mag_metadata_pro$classification <- gsub(".*;s__","",mag_metadata_pro$classification )

#mag_metadata <- mag_metadata[-c(which(mag_metadata$classification=="")) ,]





# &
# as.numeric(mag_metadata$fastani_ani)   <95),])


#mag_metadata[c(which(as.numeric(mag_metadata$)   <95),]))



mag_metadata_pro <- dplyr::select(mag_metadata_pro,"user_genome",	"base_name",	"classification","classification_full","kefir type",	"data_source","type.x","Timepoint","category","Completeness",	"Contamination", "linker")


colnames(mag_metadata_pro )[which(colnames(mag_metadata_pro )=="type.x")] <- "type"

colnames(mag_metadata_euk)[which(colnames(mag_metadata_euk)=="3" )] <- "user_genome"

mag_metadata_euk$base_name <- gsub("_metawrap_bins.*","",mag_metadata_euk$user_genome)

mag_metadata_euk$base_id <- paste(mag_metadata_euk$base_name, gsub(".*_bins_","",mag_metadata_euk$user_genome),sep="_")

mag_metadata_euk$classification <- gsub(".*;","",mag_metadata_euk$X2)


colnames(mag_metadata_euk )[which(colnames(mag_metadata_euk )=="X2")] <- "classification_full"


mag_metadata_euk$`kefir type` <- "merge"

mag_metadata_euk$`kefir type`[grep("-L-",mag_metadata_euk$user_genome )]="WL"
mag_metadata_euk$`kefir type`[grep("-G-",mag_metadata_euk$user_genome )]="WG"


mag_metadata_euk$data_source <- "Samuel et al"

mag_metadata_euk$data_source[grep("_cs_bin", mag_metadata_euk$type)] <- "This study"

mag_metadata_euk$type <- sub(".*_cs_|.*_global_","",mag_metadata_euk$type)

mag_metadata_euk$Timepoint <- NA

mag_metadata_euk$Timepoint[which(mag_metadata_euk$data_source=="This study")] <- gsub("_","",unlist(str_extract_all(mag_metadata_euk$base_name[which(mag_metadata_euk$data_source=="This study")] , "_T\\d+")))


mag_metadata_euk$ category<- "Water.kefir"


grep("MK",mag_metadata_euk$base_name )



colnames(mag_metadata_euk)[which(colnames(mag_metadata_euk)=="complete" )] <- "Completeness"


colnames(mag_metadata_euk)[which(colnames(mag_metadata_euk)=="duplicated" )] <- "Contamination"

mag_metadata_euk$linker<- NA


mag_metadata_euk$linker[which(mag_metadata_euk$data_source=="This study")] <-      gsub(".K_TG_|_merged.*","",mag_metadata_euk$base_id[which(mag_metadata_euk$data_source=="This study")] )


mag_metadata_euk$linker[which(mag_metadata_euk$data_source== "Samuel et al")] <-      gsub("-.*","",mag_metadata_euk$base_id[which(mag_metadata_euk$data_source== "Samuel et al")] )



mag_metadata_euk <- dplyr::select(mag_metadata_euk,"user_genome",	"base_name",	"classification","classification_full","kefir type",	"data_source","type","Timepoint","category","Completeness",	"Contamination", "linker")




mag_metadata <- rbind(mag_metadata_euk, 
                      mag_metadata_pro)


#fwrite(mag_metadata, "Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/mag_metadata_processed_v1.csv")                


# leave only distance values
#closest_placement_ani
###############################################################################################################################
# drep data manipulation
###############################################################################################################################


Cdb <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/drep/data_tables/Cdb.csv")


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





total_prevalence$...5  = gsub("_", " ",total_prevalence$...5)

Cdb <- 
  Cdb[-c(which(Cdb$Completeness <80 |
                 Cdb$Contamination >5)),]




"Saccharomyces cerevisiae"  
#[which(Cdb$category=="Milk.kefir")] 


cdb_prevalent <- rbind(
  
  
  
  
  Cdb[
    which(Cdb$classification%in%  total_prevalence$...5[which(total_prevalence$kefir_type=="milk")] & Cdb$category=="Milk.kefir"),],
  
  
  
  
  Cdb[
    which(Cdb$classification%in%  total_prevalence$...5[which(total_prevalence$kefir_type=="water")] & Cdb$category=="Water.kefir"),]
)









cdb_prevalent$secondary_cluster[which(cdb_prevalent$classification=="Saccharomyces cerevisiae")] <- 1


mag_prevalent_breakdown <- 
  as.data.frame(
    xtabs(~classification+data_source+secondary_cluster+category, cdb_prevalent))


mag_prevalent_breakdown$db_merge <- paste(mag_prevalent_breakdown$data_source, mag_prevalent_breakdown$category,sep="_")



pacman::p_load(upstartr)

a<- 
  
  mag_prevalent_breakdown %>% 
  mutate(detection_category=gsub("This study_Milk.kefir" ,"Milk kefir - Kefir4all",
                                 gsub("This study_Water.kefir","Water kefir - Kefir4all",
                                      gsub("Walsh et al_Milk.kefir" ,"Milk kefir - Walsh et al 2023",
                                           gsub("Walsh et al 2023_Milk.kefir" ,"Milk kefir - Walsh et al 2023",
                                                gsub( "Samuel et al_Water.kefir","Water kefir - Mortensen et al 2023", db_merge)))))) %>% 
  mutate(new_type=gsub(" -.*","",detection_category)) %>% 
  mutate(secondary_cluster_v2=gsub(".*_","",secondary_cluster)) %>% 
  filter(Freq!=0) %>% 
  mutate(secondary_cluster_v2 = factor(secondary_cluster_v2, levels = mixedsort(unique(secondary_cluster_v2)))) %>% 
  #filter(new_type=="Milk kefir") %>% 
  
  
  #ggplot(aes(x=  reorder(...2, count),y=count,fill=detection_category))+
  ggplot(aes(x = fct_reorder(classification, Freq, .fun = sum),,y=Freq,fill=secondary_cluster_v2))+
  geom_col()+
  facet_wrap(~new_type,scales="free")+
  scale_x_reordered()+
  # geom_boxplot()+
  #geom_point()+
  labs(x="Species", y="Number of MAGs detected", title="",fill="Secondary clusters") +
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




print(paste("The kefir4all dataset has the following MAG breakdown details"))

print(paste("The kefir4all dataset has ", length(which(duplicated(Cdb$user_genome))), "duplicates" ))



print(paste("The kefir4all dataset has ", length(Cdb$user_genome), "MAGs" ))

print(paste("MAGs ranged from", paste(range(Cdb$Completeness),collapse="-"), "Completeness" ))

print(paste("MAGs ranged from", paste(range(Cdb$Contamination),collapse="-"), "Contamination" ))


print(paste("The kefir4all dataset has ", length(levels(as.factor( Cdb$primary_cluster))), "primary clusters" ))

print(paste("The kefir4all dataset has ", length(levels(as.factor( Cdb$secondary_cluster))), "secondary clusters" ))




primary_by_secondary <- as.data.frame(str_split_fixed(levels(as.factor(Cdb$secondary_cluster)),"_",2))

t <- 
  as.data.frame(
    table(primary_by_secondary$V1))


nrow(t[which(t$Freq>1),])
t[which(t$Freq>1),]

str_split_fixed(levels(as.factor(Cdb$secondary_cluster)),"_",2)
xtabs(~Cdb)

length(levels(as.factor(
  Cdb$secondary_cluster)))


secondary_cluster_kefir_type_breakdown <-
  unique(
    dplyr::select(cdb_prevalent, classification, secondary_cluster, category))


as.data.frame(table(secondary_cluster_kefir_type_breakdown$classification[which(secondary_cluster_kefir_type_breakdown$category=="Water.kefir")]))%>% filter(Freq>1) 


secondary_cluster_kefir_type_breakdown <- 
  as.data.frame(
    xtabs(~classification+secondary_cluster+category, cdb_prevalent)) %>% filter(Freq!=0) 









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


clus_breakdown_prevalent$species[
which(clus_breakdown_prevalent$species=="Lactococcus raffinolactis" &
        clus_breakdown_prevalent$clust=="178_1"
        
        )]="Lactococcus raffinolactis_1"
  


clus_breakdown_prevalent$species[
  which(clus_breakdown_prevalent$species=="Lactococcus raffinolactis" &
          clus_breakdown_prevalent$clust!="178_1"
        
  )]="Lactococcus raffinolactis_2"
  
b=
clus_breakdown_prevalent %>% 
  mutate(clust=gsub(".*_","",clust)) %>% 
  
  mutate(secondary_cluster_v2 = factor(clust, levels = mixedsort(unique(clust)))) %>% 
  #ggplot(aes(x=  reorder(...2, count),y=count,fill=detection_category))+
  ggplot( aes(x = fct_reorder(species, Freq, .fun = sum),y=Freq,fill=secondary_cluster_v2 ))+
  geom_col()+
  facet_wrap(~type,scales="free")+
  scale_x_reordered()+
  # geom_boxplot()+
  #geom_point()+
  labs(x="Species", y="Number of strains detected", title="",fill="Secondary clusters") +
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



#ggplot(aes(x=  reorder(...2, count),y=count,fill=detection_category))+
ggplot(aes(x = fct_reorder(classification, Freq, .fun = sum),,y=Freq,fill=secondary_cluster_v2))+
  geom_col()+
  facet_wrap(~new_type,scales="free")+
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


























library(ggpubr)




jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 9_v2.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)


ggarrange(a,b,ncol=1,nrow=2,labels=c("A.","B."),common.legend = TRUE,font.label = list(size = 30))
graphics.off()


