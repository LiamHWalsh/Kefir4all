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

﻿#############################################################################################################################################################################################################################################

#Libraries used
#############################################################################################################################################################################################################################################

pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis,tidytext)
pacman::p_load(readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales,ggrepel,ggside,DESeq2,ggnewscale )
library(vegan)
library(ggplot2)
library(grid)

#############################################################################################################################################################################################################################################

#import and merge all seperate gene files into one big dataframe
#############################################################################################################################################################################################################################################



# 
 setwd("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/")
# 
# 
# 
# temp <- c()
# temp[["mech"]] = list.files(pattern="*mech.tsv", recursive = TRUE, full.names = TRUE)
# # #temp[["class"]] = list.files(pattern="*class.tsv", recursive = TRUE, full.names = TRUE)
# # #temp[["group"]] = list.files(pattern="*group.tsv", recursive = TRUE, full.names = TRUE)
# # temp[["gene"]] = list.files(pattern="*gene.tsv", recursive = TRUE, full.names = TRUE)
# # #looks like gene contains all the necessary info so only processed gene files and mech
# # 
# # 
# # 
# myfiles <- c()
# restistome <- c()
# for (i in names(temp)){
# 
# 
# 
# myfiles[[i]] = lapply(temp[[i]],read_delim, delim = "\t", escape_double = FALSE,
#                                 trim_ws = TRUE)
# 
# 
# 
# 
# names(myfiles[[i]]) <- gsub("./resistome_output/|.tsv","",temp[[i]])
# 
# 
# for (sample in names(myfiles[[i]])){
# 
#   restistome[[i]] <- rbind(  restistome[[i]] ,
#                                 myfiles[[i]][[sample]])
# }
# 
# if(length(names(myfiles[[i]])[-c(
#   which(names(myfiles[[i]]) %in% restistome[[i]]$Sample))])==0){
#   print(paste("alll outputs for",i,"had some value", sep=""))
# }
# 
# write.csv(  restistome[[i]],paste("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/",i,"_total_resistome.csv",sep=""),row.names = FALSE,quote=FALSE)
# }


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
########################################################################################################################
# Merge metadata into one
########################################################################################################################
total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`,Stage,Sample),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`,Stage,Sample),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`,Stage,Sample))


total_metadata $category <- NA
total_metadata $category[which(total_metadata $`kefir type` %in% c("WL","WG"))] <- "Water.kefir"
total_metadata $category[which(total_metadata $`kefir type` %in% c("ML","MG"))] <- "Milk.kefir"



#############################################################################################################################################################################################################################################

#import TOTAL gene resistome data and other data 
#############################################################################################################################################################################################################################################


gene_total_resistome <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/gene_total_resistome.csv")

mech_total_resistome <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/mech_total_resistome.csv")

reads_per_sample <-  read_csv("Q:/H2020 Master/Citizen Science Project/Results/01_preprocessing/01_readlength/reads_per_sample_interleaved_total_data.csv")

reads_per_sample$sample_id[grep("WK0",reads_per_sample$sample_id)] <- gsub("-","_",reads_per_sample$sample_id[grep("WK0",reads_per_sample$sample_id)])
reads_per_sample$X1 <- NULL

gene_total_resistome$Sample[grep("WK0",gene_total_resistome$Sample)]<- gsub("-","_",gene_total_resistome$Sample[grep("WK0",gene_total_resistome$Sample)])

gene_total_resistome$Sample <- gsub("_antimicrobial_resistance_mapped_and_unmapped","",gene_total_resistome$Sample)


mech_total_resistome$Sample[grep("WK0",mech_total_resistome$Sample)] <- gsub("-","_",mech_total_resistome$Sample[grep("WK0",mech_total_resistome$Sample)])


mech_total_resistome$Sample <- gsub("U-2-","U_2_",mech_total_resistome$Sample)

mech_total_resistome$Sample <- gsub("_antimicrobial_resistance_mapped_and_unmapped","",mech_total_resistome$Sample)


nrow(gene_total_resistome)
gene_total_resistome <- merge(gene_total_resistome,
                              reads_per_sample,
                              by.x="Sample",
                              by.y="sample_id",
                              all.x=TRUE)
nrow(gene_total_resistome)




nrow(mech_total_resistome)
mech_total_resistome <- merge(mech_total_resistome,
                              reads_per_sample,
                              by.x="Sample",
                              by.y="sample_id",
                              all.x=TRUE)
nrow(mech_total_resistome)

#############################################################################################################################################################################################################################################

#normalize TOTAL gene resistome data according to the cpm method see https://www.metagenomics.wiki/pdf/qc/RPKM and merge with metadata
#############################################################################################################################################################################################################################################



gene_total_resistome$nor_hits <- gene_total_resistome$Hits*1/gene_total_resistome$X2 * 10^6

nrow(gene_total_resistome)
gene_total_resistome <- merge(gene_total_resistome,
                             total_metadata,
                              by.x="Sample",
                              by.y="merge_column",
                              all.x=TRUE)

nrow(gene_total_resistome)

  
  gene_total_resistome$nor_hits <- gene_total_resistome$Hits*1/gene_total_resistome$X2 * 10^6
  
  
  nrow(mech_total_resistome)
  
  mech_total_resistome <- merge(mech_total_resistome,
                                total_metadata,
                                by.x="Sample",
                                by.y="merge_column",
                                all.x=TRUE)
  
  nrow(mech_total_resistome)

  mech_total_resistome$nor_hits <- mech_total_resistome$Hits*1/mech_total_resistome$X2 * 10^6

  gene_total_resistome <- cbind(gene_total_resistome,str_split_fixed(gene_total_resistome$Gene,"\\|",7) )
  
  gene_total_resistome <- gene_total_resistome %>% dplyr::rename("class"="2",
                                        "mech"="3",
                                         "group"="4",
                                         "gene"="5")

  gene_total_resistome_wide <- dplyr::select(gene_total_resistome,  nor_hits, Gene,Sample) %>% pivot_wider(names_from = "Gene",values_from = "nor_hits",values_fill = 0)
  
  
  ########################################################################################################################
  #description of the kefir restistome at the class level
  ########################################################################################################################
  


 #library(ggpubr)
  # jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/CS_metagenomics/test.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
  # 
  # pa
  # 
  # graphics.off()
  
  
  
  mech_total_resistome$kefir_category<- "Milk.kefir"
  mech_total_resistome$kefir_category[which(mech_total_resistome$`kefir type` %in% c("WG","WL"))] <- "Water.kefir"
  
  
  
  prevalence_value <- c()
  
  

  amr_class_summary <- c()
    
    for (j in  levels(as.factor(mech_total_resistome$kefir_category))) {
      
      amr_class_summary[[j]] <- 
        
        data.frame(mechanisms=as.character(levels(as.factor(  mech_total_resistome$Mechanism))),
                   prevalence=as.character(NA),
                   prevalence_per=as.character(NA),
                   max=as.character(NA),
                   min=as.character(NA),
                   mean=as.character(NA))
      
      
      for (i in levels(as.factor(  mech_total_resistome$Mechanism))){
        
        
      
      prevalence_value[[j]] <- length(  levels(as.factor(mech_total_resistome$Sample[which(mech_total_resistome$kefir_category== j)])))
      
      
      amr_class_summary[[j]]$prevalence[
      which(    amr_class_summary[[j]]$mechanisms==i)] <-  length(levels(as.factor(mech_total_resistome$Sample[which(mech_total_resistome$Mechanism==i &
                                                                                                                       mech_total_resistome$kefir_category==j )])))
      
      amr_class_summary[[j]]$prevalence_per[
        which(    amr_class_summary[[j]]$mechanisms==i)] <- 
      
      length(levels(as.factor(mech_total_resistome$Sample[which(mech_total_resistome$Mechanism==i &
                                   mech_total_resistome$kefir_category==j )]))) /     prevalence_value[[j]]
      
      
      
      
      amr_class_summary[[j]]$max[
        which(    amr_class_summary[[j]]$mechanisms==i)] <-  max(mech_total_resistome$nor_hits[which(mech_total_resistome$Mechanism==i &
                                                                                                 mech_total_resistome$kefir_category==j )])
        
        amr_class_summary[[j]]$min[
          which(    amr_class_summary[[j]]$mechanisms==i)] <- min(mech_total_resistome$nor_hits[which(mech_total_resistome$Mechanism==i &
                                                                                                        mech_total_resistome$kefir_category==j )])
        
        amr_class_summary[[j]]$mean[
          which(    amr_class_summary[[j]]$mechanisms==i)] <- mean(mech_total_resistome$nor_hits[which(mech_total_resistome$Mechanism==i &
                                                                                                        mech_total_resistome$kefir_category==j )])
        
      
        

      }
      
      amr_class_summary[[j]]$kefir=j
      
      
      amr_class_summary[[j]]$prevalence <- as.numeric(amr_class_summary[[j]]$prevalence)
      amr_class_summary[[j]]$prevalence_per <- as.numeric(amr_class_summary[[j]]$prevalence_per)
      amr_class_summary[[j]]$max <- as.numeric(amr_class_summary[[j]]$max)
      amr_class_summary[[j]]$min <- as.numeric(amr_class_summary[[j]]$min)
      amr_class_summary[[j]]$mean <- as.numeric(amr_class_summary[[j]]$mean)
  }
  
  
  nrow(
  amr_class_summary[["Milk.kefir"]][
  which(amr_class_summary[["Milk.kefir"]]$prevalence_per>0.1),])
  
  

    amr_class_summary[["Milk.kefir"]]$mechanisms[
      which(amr_class_summary[["Milk.kefir"]]$prevalence_per>0.1)]
    
   #  
   #  levels(as.factor(
   #  gene_total_resistome$Gene[which( gene_total_resistome$mech=="MLS")]))
   #  
   #  
   #  
   #  
   #  View(as.data.frame(levels(as.factor(
   #    gene_total_resistome$Gene[which( gene_total_resistome$mech=="MLS" &
   #                                       gene_total_resistome$category=="Milk.kefir")]))))
   #  View( as.data.frame( levels(as.factor(
   #    gene_total_resistome$Gene[which( gene_total_resistome$mech=="MLS" &
   #                                       gene_total_resistome$category=="Water.kefir")]))))
   #  
   #  
   #  
   #  
   #  View(as.data.frame(levels(as.factor(
   #    gene_total_resistome$Gene[which( gene_total_resistome$mech=="Aminoglycosides" &
   #                                       gene_total_resistome$category=="Milk.kefir")]))))
   # View( as.data.frame( levels(as.factor(
   #    gene_total_resistome$Gene[which( gene_total_resistome$mech=="Aminoglycosides" &
   #                                       gene_total_resistome$category=="Water.kefir")]))))
   #  
   #  
   #  levels(as.factor(
   #    gene_total_resistome$Gene[which( gene_total_resistome$mech=="Oxazolidinone")]))
    
  ########################################################################################################################
  #number of cpms per sample for each resistome 
  ########################################################################################################################
  
  
total_metadata$amr <- NA
  
  for (j in  levels(as.factor(mech_total_resistome$kefir_category))) {
    
 
  
  
  for (sample in levels(as.factor( mech_total_resistome$Sample[which(
                                                      mech_total_resistome$kefir_category==j )]))){
    
    total_metadata$amr[which(total_metadata$merge_column ==sample)] <- 
    
   sum( mech_total_resistome$nor_hits[which(mech_total_resistome$kefir_category==j &
                                          mech_total_resistome$Sample==sample )])
    
    
  }
    
    
    
    
  }
  
 amr_sample <- total_metadata 
  
 amr_sample <-  amr_sample[-c(which(is.na( amr_sample$amr))),]
 
 
 
 for (j in  levels(as.factor(mech_total_resistome$kefir_category))) {
   
   
  print(paste(j, "has", mean(amr_sample$amr[which(amr_sample$category==j)]),"value",sep=" "))
  
 }
  
  
  ########################################################################################################################
  #Alpha diversity total dataset
  ########################################################################################################################
  
  my.files_summary_shannon <- diversity(  gene_total_resistome_wide %>% column_to_rownames("Sample") )
  my.files_richness <- specnumber(  gene_total_resistome_wide %>% column_to_rownames("Sample"))
  my.files_evenness<- my.files_summary_shannon/log(my.files_richness)
  my.files_beta<- vegdist(  gene_total_resistome_wide %>% column_to_rownames("Sample"), method = "bray")
  
  my.files_summary<- cbind(shannon = my.files_summary_shannon, richness = my.files_richness, pielou = my.files_evenness,site= gene_total_resistome_wide$Sample)
  
  
  my.files_summary<- as.data.frame(my.files_summary)
  my.files_summary <- merge(my.files_summary,total_metadata,by.x="site",by.y="merge_column",all.x=TRUE)
  
  
  my.files_summary<-
  my.files_summary[-c(which(is.na(my.files_summary$`kefir type`))),]
  #my.files_summary$data_source[which(is.na(my.files_summary$`kefir type`))] <- "Walsh et al 2022" 
  
  
  
  my.files_summary <- my.files_summary[-c(which(my.files_summary$site %in% total_metadata$merge_column[which(total_metadata$`kefir type` %in% c("Extraction control", "Medium control"))])),]
  
  #my.files_summary$data_source[which(is.na(my.files_summary$data_source))] <-  "Walsh et al 2023" 
  
  my.files_summary$conditions <-  "Household conditions"
  
  my.files_summary$conditions[which(my.files_summary$data_source %in% c("Walsh et al 2023","Mortensen et al 2023" ))] <- "Laboratory controlled"
  my.files_summary$shannon <- as.numeric(  my.files_summary$shannon)
  #setwd("Q:/H2020 Master/Citizen Science Project/Plots/Evolution")
  #setwd("Q:/H2020 Master/Citizen Science Project/Plots/Evolution")
  
  #jpeg(filename='Alpha diversity_kefir_type.jpeg', width = 35*700, height=30*700,res=1700,pointsize = 15) #, width=2000, height=1950)
  ########################################################################################################################
  #tstatistically test samples in sterile conditions (global studies) and kefir4all (non sterile conditions) uing wilcoxon 
  ########################################################################################################################
  
  library(rstatix)
  wilcoxon_mk <- wilcox_test(shannon ~ conditions, data = my.files_summary[which(my.files_summary$`kefir type`=="ML"),])#
  wilcoxon_mk <- wilcoxon_mk %>% add_xy_position(x = "conditions", fun = "mean_se", scales = "free", step.increase = 0)
  wilcoxon_mk$type <-"mk"
  
  
  wilcoxon_wk <- wilcox_test(shannon ~ conditions, data = my.files_summary[which(my.files_summary$`kefir type`=="WL"),])
  

  wilcoxon_wk <- wilcoxon_wk %>% add_xy_position(x = "conditions", fun = "mean_se", scales = "free", step.increase = 0)
  wilcoxon_wk$type <-"wk"
  
  

  
  wilcoxon <- rbind(wilcoxon_mk,
                    wilcoxon_wk)
  
  wilcoxon$p_mark <- NA
  wilcoxon$p_mark[which(wilcoxon$p <0.01)] <- "****"
  library(ggpubr) # for stat_pvalue_manual
  ########################################################################################################################
  #Plot alpha diversity between samples in sterile conditions (global studies) and kefir4all (non sterile conditions)
  ########################################################################################################################
  
  
  
  dominating_species <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/dominating species/04_metacache_dominating_species.csv")
  
  
  dominating_species <- dominating_species[-c(which(duplicated(dominating_species$sample_id))),]
  
  
  my.files_summary$merge <- paste(gsub("ML","Milk Liquid",
                                       gsub("WL","Water Liquid",my.files_summary$`kefir type`)),"\n                   ",
                                  my.files_summary$conditions)
  
  my.files_summary$shannon <- as.numeric(my.files_summary$shannon)
  
  nrow( my.files_summary)
  my.files_summary <- 
  merge(
  my.files_summary,
  dominating_species,by.x="site",by.y="sample_id",all.x=TRUE)
  nrow( my.files_summary)
  
  
  
  
  
library(ggstatsplot)
  

    pb <-
    
    ggbetweenstats(data=my.files_summary %>%
                     filter(`kefir type.x` %in% c("ML","WL")), x=merge, y=shannon, 
                   grouping.var     = merge,
                   type = "nonparametric", # ANOVA or Kruskal-Wallis
                   plot.type = "box",
                   ggsignif.args    = list(textsize = 4, tip_length = 0.01),
                   pairwise.comparisons = TRUE,
                   pairwise.display = "significant",
                   centrality.plotting = FALSE,
                   bf.message = FALSE
    )+
    theme_bw()+
    xlab("Fermentation conditions")+
    ylab("Alpha diversity values (Shannon)")+
    #labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          #legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text( hjust = 1, size = 15),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")
  
  alpha_species_milk_kefir <- 
  ggbetweenstats(data=my.files_summary %>%
                   #filter(`kefir type.x` %in% c("ML","WL")), 
                   filter(category == "Milk.kefir") %>%  
                   filter(species  %in% c("Lacticaseibacillus paracasei", "Lentilactobacillus hilgardii", "Zymomonas mobilis","Lactobacillus helveticus","Lactobacillus kefiranofaciens","Lactococcus cremoris", "Lactococcus lactis")),
                 x=species, y=shannon, 
                 grouping.var     = species,
                 type = "nonparametric", # ANOVA or Kruskal-Wallis
                 plot.type = "box",
                 ggsignif.args    = list(textsize = 4, tip_length = 0.01),
                 pairwise.comparisons = TRUE,
                 pairwise.display = "significant",
                 centrality.plotting = FALSE,
                 bf.message = FALSE
  )+
    theme_bw()+
    xlab("Fermentation conditions")+
    ylab("Alpha diversity values (Shannon)")+
    #labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          #legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text( hjust = 1, size = 15,angle = 45),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")
    
    
  alpha_species_water_kefir <- 
  ggbetweenstats(data=my.files_summary %>%
                   #filter(`kefir type.x` %in% c("ML","WL")), 
                   filter(category == "Water.kefir") %>%  
                   filter(species  %in% c("Lacticaseibacillus paracasei", "Lentilactobacillus hilgardii", "Zymomonas mobilis","Lactobacillus helveticus","Lactobacillus kefiranofaciens","Lactococcus cremoris", "Lactococcus lactis")),
                 x=species, y=shannon, 
                 grouping.var     = species,
                 type = "nonparametric", # ANOVA or Kruskal-Wallis
                 plot.type = "box",
                 ggsignif.args    = list(textsize = 4, tip_length = 0.01),
                 pairwise.comparisons = TRUE,
                 pairwise.display = "significant",
                 centrality.plotting = FALSE,
                 bf.message = FALSE
  )+
    theme_bw()+
    xlab("Fermentation conditions")+
    ylab("Alpha diversity values (Shannon)")+
    #labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          #legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text( hjust = 1, size = 15,angle = 45),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")
  
  
  
  # jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/CS_metagenomics/resistome_alpha_species_diversity.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
  # 
  # 
  # 
  #           ggarrange(     alpha_species_water_kefir,     alpha_species_milk_kefir, nrow=1,ncol=2, labels=c("A.","B."),common.legend = TRUE,  font.label = list(size = 30))
  # 
  # graphics.off()
  # 
  
  

  pa2 <- 
    
    ggbetweenstats(
      data = my.files_summary,
      x = `kefir type.x`,
      y = shannon,
      fill=`kefir type`,
      color=`kefir type`,
      type = "nonparametric", # ANOVA or Kruskal-Wallis
      plot.type = "box",
      pairwise.comparisons = TRUE,
      pairwise.display = "significant",
      centrality.plotting = FALSE,
      ggsignif.args    = list(textsize = 4, tip_length = 0.01),
      bf.message = FALSE
    )+
    theme_bw()+
    xlab("Time points")+
    ylab("Alpha diversity values (Shannon)")+
    labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text(hjust = 1, size = 15),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")+
    scale_x_discrete(labels=  c('Milk - Grains', 'Milk Liquid', 'Water - Grains', 'Water - Liquid'))




  # the below just confirms we see the same trend as above with the cs samples
  
  my.files_summary_cs <- my.files_summary[which(my.files_summary$data_source.x=="This study"),]
  
  ggplot(my.files_summary_cs, aes(x=`kefir type.x`, y=as.numeric(shannon),fill=`kefir type.x`)) +
    geom_boxplot() +
    #labs(title= 'Alpha diversity of timepoints') +
    geom_point()+
    theme_bw()+
    xlab("Time points")+
    ylab("Alpha diversity values (Shannon)")+
    labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")+
    scale_x_discrete(labels=c('Milk - Grains', 'Milk - Liquid', 'Water - Grains', 'Water - Liquid'))
  
  ########################################################################################################################
  #Plot alpha diversity overtime between cs samples
  ########################################################################################################################
  my.files_summary_cs <- 
  my.files_summary_cs %>% dplyr::rename("kefir_type"="kefir type.x")
  
  nrow(my.files_summary_cs)
  my.files_summary_cs <- 
  merge(  my.files_summary_cs, 
          Citizen_Scientist_metadata_v8,
          by.x="Sample",
          by.y="ID",
          all.x=TRUE)
  
  nrow(my.files_summary_cs)
  
  
  my.files_summary_cs[which(  my.files_summary_cs$Stage=="T0" &   my.files_summary_cs$`kefir type` != "Media control"),] <- 
    my.files_summary_cs[which(  my.files_summary_cs$Stage=="T0" &   my.files_summary_cs$`kefir type` != "Media control"),] %>% 
    mutate(category_confirmed=
             gsub("GO",	"Goat",
                  gsub( "FU",	"Cow milk (whole)",
                        gsub("LO", 	"Cow milk (low fat)",
                             gsub("SK",	"Cow milk (skim)",
                                  gsub("FG",	"White sugar-dried fig",
                                       gsub("FL" , "White sugar-dried fig and fresh lemon",
                                            gsub( "AP",	"White sugar-dried apricot",
                                                  gsub( "BR",	"Brown sugar-None",
                                                        Sample    
                                                        
                                                        
                                                  )))))))))
  
  
  
  # fermentation categories 
  grouped_ggbetweenstats(
    data = my.files_summary_cs,
    x = category_confirmed,
    y = shannon,
    fill=category_confirmed,
    color=category_confirmed,
    grouping.var =kefir_type,
    type = "nonparametric", # ANOVA or Kruskal-Wallis
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    ggsignif.args    = list(textsize = 4, tip_length = 0.01),
    bf.message = FALSE
  )+
    theme_bw()+
    xlab("Time points")+
    ylab("Alpha diversity values (Shannon)")+
    labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text(hjust = 1, size = 15),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")
  
  
  
  grouped_ggbetweenstats(
    data = my.files_summary_cs,
    x =  category.y,
    y = shannon,
    fill=category_confirmed,
    color=category_confirmed,
    grouping.var =group,
    type = "nonparametric", # ANOVA or Kruskal-Wallis
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    ggsignif.args    = list(textsize = 4, tip_length = 0.01),
    bf.message = FALSE
  )+
    theme_bw()+
    xlab("Time points")+
    ylab("Alpha diversity values (Shannon)")+
    labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text(hjust = 1, size = 15),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")
  
  
  
  
  
  
  
  
  
  
  nrow(my.files_summary_cs)
  
  #[which(my.files_summary_cs$`kefir type`=="WL"),])
  
  kruskal_stage=c()
  dunn_stage=c()
  library(FSA)
  for (i in levels(as.factor(my.files_summary_cs$kefir_type))){
    
    
    kruskal_stage[[i]] <-  kruskal.test(shannon ~ Stage , data = my.files_summary_cs[which(my.files_summary_cs$kefir_type==i),])
    
    kruskal_stage[[i]]$res$P.unadj[which(  kruskal_stage[[i]]$res$P.unadj <0.05)] <- "*"
    kruskal_stage[[i]]$res$P.unadj[which(  kruskal_stage[[i]]$res$P.unadj <0.01)] <- "****"
    
    dunn_stage[[i]] <-   dunnTest(shannon ~ Stage , data = my.files_summary_cs[which(my.files_summary_cs$kefir_type==i),],
             method="bonferroni")
    
    
    
    dunn_stage[[i]]$res$P.unadj[which( dunn_stage[[i]]$res$P.unadj <0.05)] <- "*"
    dunn_stage[[i]]$res$P.unadj[which( dunn_stage[[i]]$res$P.unadj <0.01)] <- "****"
    
    
  }

  
  ggplot(my.files_summary_cs, aes(x= Stage, y=as.numeric(shannon),fill=kefir_type)) +
    geom_boxplot() +
    facet_wrap(~kefir_type)+
    #labs(title= 'Alpha diversity of timepoints') +
    geom_point()+
    theme_bw()+
    xlab("Time points")+
    ylab("Alpha diversity values (Shannon)")+
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
  
  
  
  grouped_ggbetweenstats(
    data = my.files_summary,
    x = Stage,
    y = shannon,
    fill=`kefir type.x`,
    color=`kefir type.x`,
    grouping.var =`kefir type.x`,
    type = "nonparametric", # ANOVA or Kruskal-Wallis
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    ggsignif.args    = list(textsize = 4, tip_length = 0.01),
    bf.message = FALSE
  )+
    theme_bw()+
    xlab("Time points")+
    ylab("Alpha diversity values (Shannon)")+
    labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text(hjust = 1, size = 15),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")
  
  
  ########################################################################################################################
  #  Look at fermentation categories and fermentation type
  ######################################################################################################################## 
  
  grouped_ggbetweenstats(
    data = my.files_summary_cs,
    x = Stage,
    y = shannon,
    fill=`kefir_type`,
    color=`kefir_type`,
    grouping.var =`kefir_type`,
    type = "nonparametric", # ANOVA or Kruskal-Wallis
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    ggsignif.args    = list(textsize = 4, tip_length = 0.01),
    bf.message = FALSE
  )+
    theme_bw()+
    xlab("Time points")+
    ylab("Alpha diversity values (Shannon)")+
    labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
          legend.title = element_text( size=25, face="bold"),
          axis.text.x = element_text(hjust = 1, size = 15),
          axis.text.y = element_text(hjust = 1, size = 10),
          axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
          axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
          legend.position = "none")
  
  
  
  ########################################################################################################################
  # creat pcoa plot
  ########################################################################################################################
  
  data.b.pcoa=cmdscale(my.files_beta,k=2,eig=TRUE,add = TRUE) #ordination
  
  
  
  pcoa = data.frame(PC1 = data.b.pcoa$points[,1], PC2 = data.b.pcoa$points[,2])
  percent_explained <- 100* data.b.pcoa$eig/sum(data.b.pcoa$eig)
  
  
  #############################################################################################################################
  #Clustering by dataset
  #############################################################################################################################
  pcoa <- merge(pcoa,total_metadata,by.x=0,by.y="merge_column",all.x=TRUE)
  

  pcoa <-
  pcoa[-c(which(is.na(pcoa$`kefir type`))),] 
  
  #pcoa$data_source[which(is.na(pcoa$data_source))] <- "Walsh et al - Milk kefir"
  
  
  pcoa$data_source_specific <- pcoa$data_source
  pcoa$data_source_specific[which(pcoa$data_source=="This study" &
                                    pcoa$`kefir type` %in% c("WL","WG"))] <- "Kefir4All - Water kefir"
  
  
  pcoa$data_source_specific[which(pcoa$data_source=="This study" &
                                    pcoa$`kefir type` %in% c("ML","MG"))] <- "Kefir4All - Milk kefir"
  
  pcoa$data_source_specific[which(pcoa$data_source=="Mortensen et al 2023")] <- "Breselge et al - Water kefir"
  
  pcoa$data_source_specific[which(pcoa$data_source=="Walsh et al 2023")] <- "Walsh et al - Milk kefir"
  
  
 pcoa <- 
    pcoa[-c(which(pcoa$data_source_specific=="This study")),] 
 #pcoa <- 
   #pcoa[c(which(pcoa$data_source_specific=="Walsh et al 2023")),] 
 
  
  centroid <- data.frame(data_source=as.character(levels(as.factor(pcoa$data_source_specific))),
                         PC1=as.numeric(0),
                         PC2=as.numeric(0))
  
  for (i in centroid$data_source){
    centroid$PC1[which(centroid$data_source==i)] <- mean(pcoa$PC1[which(pcoa$data_source_specific ==i)])
    centroid$PC2[which(centroid$data_source==i)] <- mean(pcoa$PC2[which(pcoa$data_source_specific==i)])
    
  }
  
  
  p_total <- 
    
    ggplot(pcoa, aes(PC1, y=PC2,colour=data_source_specific))+
    geom_point(size=5) +
    geom_point(data=centroid,size=10,shape=21, color="black",aes(fill=data_source),  show.legend=FALSE)+ #aes(, colour=data_source)),size=10)+
    #geom_convexhull(alpha=.1)+
    stat_ellipse(geom = "polygon",
                 aes(fill=data_source_specific),
                 alpha = 0.25,
                 type = "norm",
                 show.legend=FALSE)+
    #geom_text(colour="blue", check_overlap = TRUE, size=2.5, 
    #hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.025)+
    # directlabels::geom_dl(data=labels, aes(label = species), method = "smart.grid")+
    # Filter data first
    #geom_segment(data=species.long3, 
    #aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
    #colour="red", size=0.7, arrow=arrow()) +
    
    ########
  labs(x=paste("PCoA1 - ", round(percent_explained[1]), "%", sep=""), y=paste("PCoA2 - ", round(percent_explained[2]), "%", sep=""), title="") +
    coord_equal() +
    theme_bw()+
    theme(legend.position = "right",#c(.85,.2),#axis.text.x = element_blank(),  # remove x-axis text
          #axis.text.y = element_blank(), # remove y-axis text
          axis.ticks = element_blank(),  # remove axis ticks
          axis.text = element_blank(),
          axis.title.x = element_text(size=20), # remove x-axis labels
          axis.title.y = element_text(size=20), # remove y-axis labels
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          plot.background = element_blank(),
          legend.text=element_text(size = 20),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          #legend.key = element_rect(fill = "white", color = NA),
          panel.border = element_rect(colour = "black", fill=NA, size = 2.5))+
    #legend.box.background = element_rect(colour = "red"))+
    scale_colour_discrete(name="Data source")+
    scale_shape_discrete(name="Data source")+
    guides(colour = guide_legend(override.aes = list(size=10)))
  
  #############################################################################################################################
  # look at pcoa pecific for kefir4alll
  #############################################################################################################################
  
  
  
  my.files_beta_cs<- vegdist(  gene_total_resistome_wide[which(gene_total_resistome_wide$Sample %in% kefir4all_metadata$merge_column[-c(which(kefir4all_metadata$`kefir type` %in% c("Extraction control", "Medium control")))]), ] %>% column_to_rownames("Sample"), method = "bray")
  
  
  
  
  data.b.pcoa_cs=cmdscale(my.files_beta_cs,k=2,eig=TRUE,add = TRUE) #ordination
  
  
  
  pcoa_cs = data.frame(PC1 = data.b.pcoa_cs$points[,1], PC2 = data.b.pcoa_cs$points[,2])
  percent_explained <- 100* data.b.pcoa_cs$eig/sum(data.b.pcoa_cs$eig)
  
  
  #############################################################################################################################
  #Clustering kefir4all  samples  by time/category
  #spolier nothing 
  #############################################################################################################################
  
  kefir4all_metadata$category_kefir <- NA
  
  kefir4all_metadata$category_kefir[
which(kefir4all_metadata$`kefir type` %in%  c("MG","ML"))] <- "Milk.kefir"

  
  kefir4all_metadata$category_kefir[
    which(kefir4all_metadata$`kefir type` %in%  c("WG","WL"))] <- "Water.kefir"
  
            
  
  
  pcoa_cs<- merge(pcoa_cs,kefir4all_metadata,by.x=0,by.y="merge_column",all.x=TRUE)
  
  pcoa_cs <- 
    merge(  pcoa_cs,Citizen_Scientist_metadata_v8,by.x="Sample",by.y="ID",all.x=TRUE)
  
  
  
  library(stringi)
  
  pcoa_cs[which(pcoa_cs$Stage=="T0" & pcoa_cs$`kefir type` != "Media control"),] <- 
    pcoa_cs[which(pcoa_cs$Stage=="T0" & pcoa_cs$`kefir type` != "Media control"),] %>% 
    mutate(category_confirmed=
             gsub("GO",	"Goat",
                  gsub( "FU",	"Cow milk (whole)",
                        gsub("LO", 	"Cow milk (low fat)",
                             gsub("SK",	"Cow milk (skim)",
                                  gsub("FG",	"White sugar-dried fig",
                                       gsub("FL" , "White sugar-dried fig and fresh lemon",
                                            gsub( "AP",	"White sugar-dried apricot",
                                                  gsub( "BR",	"Brown sugar-None",
                                                        Sample    
                                                        
                                                        
                                                  )))))))))
  
  
  
  ggplot(pcoa_cs, aes(PC1, y=PC2,colour= category_confirmed,fill= category_confirmed))+
    geom_point(size=5) +
    facet_wrap(~`kefir type`)+
    #geom_convexhull(alpha=.1)+
    #stat_ellipse(geom = "polygon",
    #alpha = 0.25,
    #type = "norm")+
    #geom_text(colour="blue", check_overlap = TRUE, size=2.5, 
    #hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.025)+
    # directlabels::geom_dl(data=labels, aes(label = species), method = "smart.grid")+
    # Filter data first
    #geom_segment(data=species.long3, 
    #aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
  #colour="red", size=0.7, arrow=arrow()) +
  
  ########
  labs(x=paste("PCoA1 - ", round(percent_explained[1]), "%", sep=""), y=paste("PCoA2 - ", round(percent_explained[2]), "%", sep=""), title="") +
    coord_equal() +
    theme_bw()+
    theme(#legend.position = "none",#axis.text.x = element_blank(),  # remove x-axis text
      #axis.text.y = element_blank(), # remove y-axis text
      axis.ticks = element_blank(),  # remove axis ticks
      axis.text = element_blank(),
      axis.title.x = element_text(size=18), # remove x-axis labels
      axis.title.y = element_text(size=18), # remove y-axis labels
      panel.background = element_blank(), 
      panel.grid.major = element_blank(),  #remove major-grid labels
      panel.grid.minor = element_blank(),  #remove minor-grid labels
      plot.background = element_blank(),
      legend.text=element_text(face="italic",size = 15),
      legend.key.size = unit(1, 'cm'), #change legend key size
      legend.key.height = unit(1, 'cm'), #change legend key height
      legend.key.width = unit(1, 'cm'), #change legend key width
      legend.title = element_text(size=20) #change legend title font size
    )
  
  
  
  
  #############################################################################################################################
  #permanova based on type
  #############################################################################################################################

  

  
  
 # metadata_adonis <- kefir4all_metadata[which(kefir4all_metadata$merge_column %in% rownames(as.data.frame(as.matrix(my.files_beta_cs)))),]
 
  metadata_adonis <-   pcoa_cs[which(  pcoa_cs$Row.names %in% rownames(as.data.frame(as.matrix(my.files_beta_cs)))),]
  
  colnames(metadata_adonis)[which(colnames(metadata_adonis)=="Row.names")] <- "merge_column"

  

  metadata_adonis <- 
    metadata_adonis[order(match(metadata_adonis$merge_column,rownames(as.data.frame(as.matrix(my.files_beta_cs))))),]
  
  identical( metadata_adonis$merge_column,rownames(as.data.frame(as.matrix(my.files_beta_cs))) )
  
  
tester <-  as.data.frame( cbind(metadata_adonis$merge_column,rownames(as.data.frame(as.matrix(my.files_beta_cs))) ))

length(
which(tester$V1 != tester$V2))

length(metadata_adonis$merge_column)


length(rownames(as.data.frame(as.matrix(my.files_beta_cs))) )


  
  
  adonis<-c()
  adonis[["Kefir_type"]] <- 
    adonis2( my.files_beta_cs~`kefir type`, data= metadata_adonis ,permutations =  10000)
  
  
  # adonis[["Kefir_type_stage"]] <- 
  #   adonis2( my.files_beta_cs~`kefir type`+Stage, data= metadata_adonis ,permutations =  10000)

  res <- c()
  library(pairwiseAdonis)
  #default is 999 permutations
  res[["Kefir_type"]] <-pairwiseAdonis::pairwise.adonis(my.files_beta_cs,metadata_adonis$`kefir type`)
  

 #  test <- 
 #  merge(
 # as.data.frame(as.matrix(my.files_beta_cs)),
 #  metadata_adonis,
 #  by.x=0,
 #  by.y="merge_column",
 #  all.x=TRUE)
  
  #############################################################################################################################
  #permanova  based on fermentation categories
  #############################################################################################################################
  
  metadata <- c()
  amr_profile <-c() 
  
  
  
  for (i in c("Milk.kefir", "Water.kefir")){
  
    metadata <-  pcoa_cs[which( pcoa_cs$category_kefir==i &
                                     pcoa_cs$Row.names %in% rownames(as.data.frame(as.matrix(my.files_beta_cs)))),]
    
      metadata
    
    amr_profile <- as.data.frame(as.matrix(my.files_beta_cs))
      
    amr_profile <-     amr_profile[which(     rownames(amr_profile ) %in%  metadata$Row.names), which(     colnames(amr_profile ) %in%  metadata$Row.names) ]
      
    
    
    

    
    
    metadata  <- 
      metadata [order(match(metadata$Row.names,rownames(    amr_profile))),]
    
    #identical( metadata$Row.names,rownames(  amr_profile) )
    
    
    tester <-  as.data.frame( cbind(metadata$Row.names, rownames(amr_profile))) 
    
    if(
    length(
      which(tester$V1 != tester$V2))!= 0){
      print(paste("category confirmed metadata does not match the AMR profile matrix"))
        break
      }
    
    metadata$category[
   which(is.na(metadata$category))] <- "Baseline"
    
  adonis[["category_confirmed"]][[i]] <- 
    adonis2( amr_profile~`category_confirmed`, data= metadata ,permutations =  10000)
  
  res[["category_confirmed"]][[i]] <-pairwiseAdonis::pairwise.adonis(amr_profile,metadata$category_confirmed)
  
  
  
  
  adonis[["fermenter_category"]][[i]] <- 
    adonis2( amr_profile~category, data= metadata ,permutations =  10000)
  
  res[["fermenter_category"]][[i]] <-pairwiseAdonis::pairwise.adonis(amr_profile,metadata$category)
  
  
  }
  
  #############################################################################################################################
  #permanova  based on stage 
  #############################################################################################################################
  

  t <- c()
  t1 <- c()
  
  p_pairwise_stage <- c()
  p_pairwise_stage_total <- c()
  for (i in levels(as.factor(metadata_adonis$`kefir type`))){
    
  t <- metadata_adonis[which(metadata_adonis$`kefir type`==i),]
  
  t1 <- as.data.frame(as.matrix(my.files_beta_cs))

  t1 <- t1[which(rownames(t1) %in% t$merge_column ),which(colnames(t1) %in% t$merge_column )]
    
    
    metadata_adonis <- 
      metadata_adonis[order(match(metadata_adonis$merge_column,rownames(as.data.frame(as.matrix(my.files_beta_cs))))),]
    
    if(
    identical( c(t$merge_column),c(rownames(t1)) )){
      print(paste("for", i, "dataframes match up"))
    }else{
      print(paste("for", i, "dataframes do match up, going to crash now"))
break
    }
    
 
    p_pairwise_stage[[i]] <-   pairwiseAdonis::pairwise.adonis2(t1 ~ Stage, data = t)
    
    
    
    
      
      
      for (j in names(p_pairwise_stage[[i]])){
        if (j=="parent_call"){
          next
        }else{
          p_pairwise_stage[[i]][[j]] <- as.data.frame(p_pairwise_stage[[i]][[j]])
          p_pairwise_stage[[i]][[j]]$kefir_type=i
          p_pairwise_stage[[i]][[j]]$comparision_type=j
          
          p_pairwise_stage_total <- rbind(p_pairwise_stage_total,
                                          p_pairwise_stage[[i]][[j]])
        }
      }
    
  }
  p_pairwise_stage_total <- 
  p_pairwise_stage_total[ grep("Stage",rownames(p_pairwise_stage_total) ),]

  p_pairwise_stage_total$fdr <- p.adjust(  p_pairwise_stage_total$`Pr(>F)`, method ="bonferroni")
  

  
  
  #############################################################################################################################
  #Cconfirm no clear clustering based on community type using permanova 
  #############################################################################################################################
  
  

  
  # Are p values are siginificant all that means is that one of the groups at least is significnatly different from the others that  raises the question of which pairs are different from each other
  
  
  
dominating_species <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/dominating species/04_metacache_dominating_species.csv")
  
  
  dominating_species <- dominating_species[-c(which(duplicated(dominating_species$sample_id))),]
  
  
  
  metadata_adonis <- 
    merge(metadata_adonis,
          dominating_species,
          by.x="merge_column",
          by.y="sample_id",
          all.x=TRUE)
    

identical( metadata_adonis$merge_column,rownames(as.data.frame(as.matrix(my.files_beta_cs))) )



pcoa <- merge(pcoa,
              dominating_species,
              by.x="Row.names",
              by.y="sample_id",
              all.x=TRUE)

species <- 
table(
  pcoa$species)

pcoa$species[which(  pcoa$species %in% names(which(species<80)))] <- "Other"


pc <- 
  ggplot(pcoa, aes(PC1, y=PC2,colour=species,shape=species))+
  geom_point(size=5) 
  
  
  
  
#############################################################################################################################
#envfit anlaysis on total dataset
#############################################################################################################################


  #adonis2( my.files_beta_cs~`kefir type.x`+species, data= metadata_adonis ,permutations =  10000)
  
  
  t <- c()
  t1 <- c()
  
  p_pairwise_species <- c()
  p_pairwise_species_total <- c()
  for (i in levels(as.factor(metadata_adonis$group))){
    
    t <- metadata_adonis[which(metadata_adonis$group==i),]
    
    
    species <- 
    table(t$species)
    t <-  t[which(t$species %in% names(which(species>10))),]
    
    t1 <- as.data.frame(as.matrix(my.files_beta_cs))
    
    t1 <- t1[which(rownames(t1) %in% t$merge_column ),which(colnames(t1) %in% t$merge_column )]
    
    
    
    if(
      identical( c(t$merge_column),c(rownames(t1)) )){
      print(paste("for", i, "dataframes match up"))
    }else{
      print(paste("for", i, "dataframes do match up, going to crash now"))
      break
    }
    
    
    p_pairwise_species[[i]] <-   pairwiseAdonis::pairwise.adonis2(t1 ~species, data = t)
    
    
    
    
    for (j in names(p_pairwise_species[[i]])){
      if (j=="parent_call"){
        next
      }else{
        p_pairwise_species[[i]][[j]] <- as.data.frame(p_pairwise_species[[i]][[j]])
        p_pairwise_species[[i]][[j]]$kefir_type=i
        p_pairwise_species[[i]][[j]]$comparision_type=j
        
        p_pairwise_species_total <- rbind(p_pairwise_species_total,
                                        p_pairwise_species[[i]][[j]])
      }
    }
    
  }
  
  
  p_pairwise_species_total <- 
    p_pairwise_species_total[ grep("species",rownames(p_pairwise_species_total) ),]
  
  p_pairwise_species_total$fdr <- p.adjust(  p_pairwise_species_total$`Pr(>F)`, method ="bonferroni")
  
  p_pairwise_species_total$comparision_type[which(p_pairwise_species_total$fdr<0.05)]
  
  
  
  
  
  metadata_adonis$conditions <-  "Household conditions"
  
  metadata_adonis$conditions[which(metadata_adonis$data_source %in% c("Walsh et al 2023","Mortensen et al 2023" ))] <- "Laboratory controlled"
  
  
  metadata_adonis$conditions[which(metadata_adonis$data_source.x == "This study" &
                                     metadata_adonis$Stage=="T0" )] <- "Laboratory controlled"
  
  
  t <- c()
  t1 <- c()
  
  p_conditions <- c()
  my.files_beta <- c()
  
  bray_distance <- c()
    
  functional_data=c()
  
  correlation_data <- c()
  
  for (i in levels(as.factor(metadata_adonis$group))){
    
    t <- metadata_adonis[which(metadata_adonis$group==i),]
    
    
    t1 <- as.data.frame(as.matrix(my.files_beta_cs))
    
    t1 <- t1[which(rownames(t1) %in% t$merge_column ),which(colnames(t1) %in% t$merge_column )]
    
    
    
    if(
      identical( c(t$merge_column),c(rownames(t1)) )){
      print(paste("for", i, "dataframes match up"))
    }else{
      print(paste("for", i, "dataframes do match up, going to crash now"))
      break
    }
    
    
    p_conditions[[i]] <-   adonis2(t1 ~conditions, data = t)
    

    for (grain_name in levels(as.factor( metadata_adonis$`kefir type.x`[which( metadata_adonis$group==i)]))[1]){
      print(grain_name)
    
    
    my.files_beta[[grain_name]] <- t1[which(rownames(t1)==metadata_adonis$merge_column[which(metadata_adonis$Stage=="T0"&
                                                                                         metadata_adonis$`kefir type.x` %in% grain_name)][1]),]                                       
    
    
    my.files_beta[[grain_name]]<- 
      as.data.frame(my.files_beta[[grain_name]] ) %>% rownames_to_column("reference") %>% pivot_longer(!reference, names_to = "merge_column",values_to = "distance")


    
    
    my.files_beta[[grain_name]] <- 
      merge(metadata_adonis,
            my.files_beta[[grain_name]] , by="merge_column",all.y=TRUE)
    
    
    # this will restrict everything to just grains
    my.files_beta[[grain_name]]<-
      my.files_beta[[grain_name]][which(my.files_beta[[grain_name]]$`kefir type.x`==grain_name),]
    
    
    
    #my.files_beta[[grain_name]]$type <- gsub("\\."," ",grain_name)
    
    bray_distance[[grain_name]] <- 
      ggbetweenstats(data = my.files_beta[[grain_name]] %>% 
                       filter(Stage!="T0"),
                     x=Stage, 
                     y=distance,
                    title=grain_name,
                     type = "nonparametric", # A
                     ggsignif.args    = list(textsize = 2, tip_length = 0.01)
                     
      )+
      theme_bw()+
      xlab("Time points")+
      ylab("Bray Curtis distance")+
      labs(fill = "Timepoint")+
      #guides(colour = guide_legend(override.aes = list(size=25)))+
      #ylim(0,250)+
      theme(#plot.title = element_text(hjust = 0.5,size=35,face="bold"),
        legend.title = element_text( size=25, face="bold"),
        axis.text.x = element_text( hjust = 1, size = 15),#angle = 45,
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
        axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
        legend.position = "none",
        plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
                                                    box.color = "#748696",size=35, lineheight = 2))                                                                                     
    
    


  
    
    for (mech in 
   levels(as.factor( mech_total_resistome$Mechanism))){
      
      functional_data <- 
      mech_total_resistome[which(mech_total_resistome$Mechanism==mech),]
    
   data <- 
    merge(dplyr::select(functional_data,Sample, nor_hits,Mechanism),
          dplyr::select(my.files_beta[[grain_name]],merge_column,reference,distance),#-c(`kefir type.y`, data_source.y,Sample )),
          by.x="Sample",
          by.y="merge_column",
          all.y="TRUE")
    
          data <-   data %>% 
      dplyr::rename("nor_hits_sample"="nor_hits")
    
    
    
    data <- 
      merge(dplyr::select(functional_data,Sample, nor_hits,Mechanism),
           data,
            by.y="reference",
            by.x="Sample",
           all.y="TRUE")
    
    data <-     data %>% 
      dplyr::rename("nor_hits_reference"="nor_hits")
    
    
    
    data$nor_hits_reference[which(is.na(data$nor_hits_reference))] <- 0
    data$nor_hits_sample[which(is.na(data$nor_hits_sample))] <- 0
    
    
    data$nor_hits_change <- data$nor_hits_sample- data$nor_hits_reference
    
    
    res<- cor.test(as.numeric(   data$nor_hits_change),
                   as.numeric(   data$distance),
                   method="kendall")
    
    
    correlation_data[[i]] <- rbind(data.frame(
      group_type=as.character(i),
      mech_type=as.character(mech),
      p_value=as.numeric(res$p.value),
      r=as.numeric(res$estimate)),
      correlation_data[[i]])
    
  
  
}
}

    correlation_data[[i]]$fdr <- p.adjust(  correlation_data[[i]]$p_value, method ="bonferroni")
    
    
    
  }
  correlation_data[["water"]] <- 
  correlation_data[["water"]][-c(which(is.na( correlation_data[["water"]]$fdr))),]
  correlation_data[["milk"]] <- 
  correlation_data[["milk"]][-c(which(is.na( correlation_data[["milk"]]$fdr))),]
  
  
  correlation_data[["water"]][  which(correlation_data[["water"]]$fdr<0.05),]
  correlation_data[["milk"]][  which(correlation_data[["milk"]]$fdr<0.05),]
  
  #more correlations with datasource and kefir type
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 #t <-  dplyr::select(mech_total_resistome, Sample, Mechanism, nor_hits) %>% pivot_wider(names_from = Mechanism, values_from = nor_hits) %>% column_to_rownames("Sample")
  
  mech_total_resistome$kefir_type_v2 <- NA
  
  mech_total_resistome$kefir_type_v2[
 grep("L", mech_total_resistome$`kefir type`)] <- "Liquid"
 
  
  mech_total_resistome$kefir_type_v2[
    grep("G", mech_total_resistome$`kefir type`)] <- "Grain"
  
  
  biplot_data <- c()
  
  
  mech_total_resistome <- mech_total_resistome[-c(which(mech_total_resistome$Sample %in%  c("TG_ID084_WL_T2_S369", "TG_ID086_WL_T2_S310"))),]
  
 
  
  
  mech_total_resistome$conditions <-  "Household conditions"
  
  mech_total_resistome$conditions[which(mech_total_resistome$data_source %in% c("Walsh et al 2023","Mortensen et al 2023" ))] <- "Laboratory controlled"
  
  
  mech_total_resistome$conditions[which(mech_total_resistome$conditions == "This study" &
                                          mech_total_resistome$Stage=="T0" )] <- "Laboratory controlled"
  
  
  
  
  pa <-  mech_total_resistome %>%  
    filter(`kefir type` %in% c("MG","ML","WG","WL")) %>% 
    mutate(`Fermentation parameter`=conditions) %>% 
    mutate(sample_v2=paste(conditions,Sample,sep="_")) %>% 
    mutate(Mechanism_v2=gsub(" ","\n", Mechanism)) %>% 
    mutate(Mechanism_v2=gsub("Drug\nand\nbiocide\nresistance","Drug and\nbiocide esistance", Mechanism_v2)) %>% 
    mutate(Mechanism_v2=gsub("Cationic\nantimicrobial\npeptides","Cationic\nantimicrobial peptides", Mechanism_v2)) %>% 
    mutate(Mechanism_v2=gsub("Drug\nand\nbiocide\nand\nmetal\nresistance","Drug and biocide and\nmetal resistance", Mechanism_v2)) %>% 
    
      #mutate(conditions=as.factor(conditions)) %>% 
    #  filter(`kefir type` %in% c("MG","ML","WG","WL")) %>% 
    # arrange(Sample,Mechanism) %>% 
    ggplot(aes(y = sample_v2,x =Mechanism_v2)) +
    geom_tile(aes(fill = nor_hits),#color = "white",
              #lwd = 1.5,
              linetype = 1) +
    labs(x="", y="", title="", fill="AMR-encoding reads (CPM)")+ #y="Feature"
    #geom_ysidetile(aes(x = "Fermentation conditions", yfill = conditions))+
    
    facet_wrap(~`kefir type` ,scales="free",
               labeller = labeller(`kefir type` = c(
                 'MG'=  'Milk - Grains',
                 'ML'= 'Milk - Liquid',
                 'WG'= 'Water - Grains',
                 'WL'= 'Water - Liquid')))+
    scale_y_reordered() +
    
    scale_fill_gradient(low = "green", high = "red") + 
    # scale_fill_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
    #                      breaks= c(-0.05,0,0.10,0.30,0.70,0.90,1),
    #                      values= rescale(as.numeric(c(-0.05,0,0.10,0.30,0.70,0.90,1))),
    #                      guide="colorbar",
    #                      name="Bray Curtis dissimilarity" ,
    #                      labels=c(-0.05,0,0.10,0.30,0.70,0.90,1))+
    theme_bw()+
    theme(#plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5),size = 55), #element_text(color="red", size=14, face="bold",hjust = 0.5),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=12.5,angle = 90,hjust = 1),
      strip.background = element_rect(
        color="black", fill="white"),
      strip.text = element_text(size=25),
      #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
      axis.title.y = element_text( size=10, face="bold",hjust = 0.5, vjust = 1.5),
      legend.key.size = unit(1, 'cm'), #change legend key size
      legend.key.height = unit(1, 'cm'), #change legend key height
      legend.key.width = unit(1, 'cm'), #change legend key width
      legend.text = element_text(size=10),
      legend.title = element_text(size=20, vjust = 0.5),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing.x = unit(-5, "lines"))+
    geom_ysidetile(aes(x = "Fermentation\nconditions", yfill = `Fermentation parameter`))#+
  #  labs(fill="Fermention parameter")
  #scale_y_discrete(breaks=test$Sample)
  
  
  library(ggbiplot)
  
  
  pcoa_data <- c()
  envfit_output <- c()
  i="Milk.kefir"
  
  a <- 
  nrow(  mech_total_resistome)
  
  mech_total_resistome <- 
    merge(  mech_total_resistome, 
            dplyr::select(dominating_species,sample_id,species),
            by.x="Sample",
            by.y="sample_id",
            all.x=TRUE)
  
  b <- 
    nrow(  mech_total_resistome)
  
  if(a==b){
    print("merge command above was correct")
  }else{
    print("merge command above needs more attention")
  }
  
  correlations <- c()
  env_correlations <- c()
  
  for (i in c("Milk.kefir","Water.kefir")){
    

    t <-  dplyr::select(mech_total_resistome, -c(Hits,X2)) %>% pivot_wider(names_from = Mechanism, values_from = nor_hits) %>% column_to_rownames("Sample")                         
    
    t[is.na(t)]=0
    
    
    t <-
    t[which(t$category==i),]
    
    t_names <- t
    
    t1 <- t
    
    
    metadata <-  dplyr::select( t1, c(data_source,category, kefir_type_v2,conditions,Stage,species)) # Stage needs to be included in cs only dataset
      
      
    t <- 
    dplyr::select( t, -c(data_source, `kefir type`, Stage, Sample.y,category, kefir_category, kefir_type_v2,conditions,species))
    
    t <- t[,-c(which(colnames(t) %in% names(which(colSums(t)==0)) ))]
    

    # pca <- prcomp(t,
    #               scale = TRUE)
    
    
    
# 
#    
#     biplot_data[[i]] <- 
#     ggbiplot(pca,
#             # groups = t_names$data_source,
#              groups = t_names$conditions ,
#             # groups =t_names$kefir_type_v2,
#              #labels = crime$st,
#              # labels.size = 4,
#              var.factor = 1.4,
#              ellipse = TRUE, ellipse.level = 0.5, ellipse.alpha = 0.1,
#              circle = TRUE,
#              varname.size = 4,
#              varname.color = "black") +
#       labs(fill = "Region", color = "Region")+
#       theme(legend.direction = 'horizontal', legend.position = 'top')
#     
#     
    
    
    
    ######################################
    # envfit
    x_matrix<-vegdist(t,method = "bray")
    pcoa <- cmdscale(x_matrix, k=2)

    #
    data.scores = as.data.frame(pcoa)
    
    data.scores =
    merge(data.scores,
          t1,
          by=0,
          all.x=TRUE)
    
    en_pcoa <- vegan::envfit(pcoa,   t1[,which(colnames(t1) %in%   c("kefir_type_v2","conditions","category_confirmed","Fermentation type","species",levels(as.factor(mech_total_resistome$Mechanism)) ) ) ], permutations = 1000)
 
    

   # en_pcoa <- vegan::envfit(pcoa, dplyr::select(t1,-c(Sample.y)), permutations = 1000)
    
    en_pcoa$vectors$pvals_fdr <-  p.adjust(en_pcoa$vectors$pvals,method = "bonferroni")
    

    envfit_output[[i]] <- en_pcoa

    #en_pcoa <- vegan::envfit(pcoa, dplyr::select(t1,-c(data_source, `kefir type`, Stage, Sample.y,category, kefir_category, kefir_type_v2)), permutations = 1000)
    
    #en_pcoa <- vegan::envfit(pcoa ~ conditions, t1, permutations = 1000)
    

#     library(MASS)
#     t1$conditions <- as.factor(    t1$conditions)
# lda(conditions~., data=t1)
#     
# str(t1
#     )
    

    
    
    A <- as.list(en_pcoa$vectors)
    pvals<-as.data.frame(A$pvals)
    pvals_fdr<-as.data.frame(A$pvals_fdr)
    arrows<-as.data.frame(A$arrows*sqrt(A$r))
    C<- cbind(arrows, pvals,pvals_fdr, as.data.frame(A$r))
    
    Cred<-subset(C,pvals_fdr<0.05)
    Cred <- cbind(Cred, Species = rownames(Cred))
    
    colnames(Cred)[which(colnames(Cred)=="A$r")] <- "cor"
    

    

    #https://www.mdpi.com/1420-3049/26/20/6254 - for interpretatting pca plots
    #https://cran.r-project.org/web/packages/ggbiplot/readme/README.html
   # https://medium.com/@RaharditoDP/principal-component-analysis-with-biplot-analysis-in-r-ee39d17096a1
 

    
    for (j in levels(as.factor( mech_total_resistome$Mechanism))){
      
      a <- cor.test(getElement(data.scores,"V1"), getElement(data.scores,j),method = "spearman", exact = FALSE)
    
      
     b <-  cor.test(getElement(data.scores,"V2"), getElement(data.scores,j),method = "spearman", exact = FALSE)
     
     correlations[[i]] <- rbind(data.frame(name=j, estimate.x=a$estimate, p.value.x= a$p.value,
                estimate.y=b$estimate, p.value.y= b$p.value),  correlations[[i]] )
     
      
    }
    
    
    
    
    for (target in c("kefir_type_v2","conditions","species")){
      
      data.scores_plot <- dplyr::select(data.scores,Row.names,V1,V2, target,levels(as.factor(mech_total_resistome$Mechanism))) 
      
    
      colnames(data.scores_plot)[which(      colnames(data.scores_plot)==target)] <- "target"
      
      
      if(target=="species"){
        if (i=="Milk.kefir"){
          
          data.scores_plot$target[-c(which(data.scores_plot$target %in% c("Lactobacillus helveticus","Lactobacillus kefiranofaciens","Lactococcus cremoris", "Lactococcus lactis")))] <- "Other"
          #work here 

                                 

        }else{
          data.scores_plot$target[-c(which(data.scores_plot$target %in%  c("Lacticaseibacillus paracasei", "Lentilactobacillus hilgardii", "Zymomonas mobilis")))] <- "Other"
          
        }
      }
      
    centroid <- data.frame(target=as.character(levels(as.factor(data.scores_plot$target))),
                           V1=as.numeric(0),
                           V2=as.numeric(0))
    
    for (e in centroid$target){
      centroid$V1[which(centroid$target==e)] <- mean(data.scores_plot$V1 [which(data.scores_plot$target ==e)])
      centroid$V2[which(centroid$target==e)] <- mean(data.scores_plot$V2[which(data.scores_plot$target==e)])
      
    }#end of e
    
    
    
    #test <- 
   #dplyr::select( data.scores, c( Row.names, V1, V2, levels(as.factor( mech_total_resistome$Mechanism)))) %>% pivot_longer(!c( Row.names, V1, V2),values_to = "RA", names_to = "mech") 

    
    library(ggnewscale)
    
    colnames(data.scores_plot)[which(colnames(data.scores_plot)=="target")] <- "Fermentation parameter"
   # colnames(data.scores_plot)[4] <- "Fermentation parameter"
    

    if(i=="Milk.kefir"){
    
   pcoa_data[[i]][[target]] <- 
      
    ggplot(data = data.scores_plot, aes(x = V1, y = V2)) + 
    # geom_point(data = data.scores, aes(colour = conditions), size = 3, alpha = 0.5) + 
     geom_point(data = data.scores_plot, aes(colour= `Fermentation parameter`), size = 3, alpha = 0.5) + 
     guides(color = guide_legend(order = 1))+ # this miracle peiece of code sets the fermentation paramter legend first
     geom_point(data=centroid,size=10,shape=21, color="black",aes(fill=target),  show.legend=FALSE)+ 
     stat_ellipse(geom = "polygon",
                  aes(fill= `Fermentation parameter`),
                  alpha = 0.25,
                  type = "norm",
                  show.legend=FALSE)+
      #scale_colour_manual(values = c("orange", "steelblue"))  + 
     #scale_fill_manual(values = c("orange", "steelblue"))  + 
     new_scale_color()+
    # geom_point(data = data.scores_plot, aes(fill= `Fermentation parameter`), size = 3, alpha = 0.5) + 
      # geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
      #              data =  Cred, size =1, alpha = 0.5, colour = "grey30") +
    
      geom_text_repel(data = Cred, aes(x = Dim1, y = Dim2, colour = cor), 
                label = row.names(Cred), fontface = "bold",min.segment.length=0.01) + 
     geom_point(data = Cred, aes(x = Dim1, y = Dim2, colour = cor), 
                shape = "diamond", size = 4, alpha = 0.6) +
      #geom_text(data = Cred, aes(x = Dim1, y = Dim2), colour = "grey30", 
                #fontface = "bold", label = row.names(Cred)) + 
      theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
            panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
            axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
            legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
            legend.text = element_text(size = 9, colour = "grey30"),
            legend.position = c(.6,.9),
            legend.direction="horizontal",
            legend.box="vertical",
            plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
                                                        box.color = "#748696",size=25, lineheight = 2)) + 
      labs(colour = "Variance (r2)",
           x="PC-1",
           y="PC-2")+
     scale_color_distiller(palette = "RdYlBu")
  
   
    }else{

      
      
      pcoa_data[[i]][[target]] <- 
      ggplot(data = data.scores_plot, aes(x = V1, y = V2)) + 
        # geom_point(data = data.scores, aes(colour = conditions), size = 3, alpha = 0.5) + 
        geom_point(data = data.scores_plot, aes(colour= `Fermentation parameter`), size = 3, alpha = 0.5) + 
        guides(color = guide_legend(order = 1))+ # this miracle peiece of code sets the fermentation paramter legend first
        geom_point(data=centroid,size=10,shape=21, color="black",aes(fill=target),  show.legend=FALSE)+ 
        stat_ellipse(geom = "polygon",
                     aes(fill= `Fermentation parameter`),
                     alpha = 0.25,
                     type = "norm",
                     show.legend=FALSE)+
        #scale_colour_manual(values = c("orange", "steelblue"))  + 
        #scale_fill_manual(values = c("orange", "steelblue"))  + 
        new_scale_color()+
        # geom_point(data = data.scores_plot, aes(fill= `Fermentation parameter`), size = 3, alpha = 0.5) + 
        # geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
        #              data =  Cred, size =1, alpha = 0.5, colour = "grey30") +
        
        geom_text_repel(data = Cred, aes(x = Dim1, y = Dim2, colour = cor), 
                        label = row.names(Cred), fontface = "bold",min.segment.length=0.01) + 
        geom_point(data = Cred, aes(x = Dim1, y = Dim2, colour = cor), 
                   shape = "diamond", size = 4, alpha = 0.6) +
        #geom_text(data = Cred, aes(x = Dim1, y = Dim2), colour = "grey30", 
        #fontface = "bold", label = row.names(Cred)) + 
        theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
              panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
              axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
              legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
              legend.text = element_text(size = 9, colour = "grey30"),
              legend.position ="none",
              legend.direction="horizontal",
              legend.justification = c(1, 1) )+
        #legend.box="Horizontal") + 
        labs(colour = "Variance (r2)",
             x="PC-1",
             y="PC-2")+
        scale_color_distiller(palette = "RdYlBu")
      
      
    }
  }
   
    
   ######################################
   # identify which vectors are closest to the centroid points  
   
   Cred$closest_unit <- NA
  

   Cred$dist_h <- NA
   Cred$dist_l <- NA
   
  h_c.x <-  centroid$V1[1]
   
  h_c.y <- centroid$V2[1]
  
  l_c.x <-  centroid$V1[2]
  
  l_c.y <- centroid$V2[2]
  
  
  
  for (j in rownames(Cred)){
    
    
  x2 <- Cred$Dim1[which(rownames(Cred)==j)]
  
  y2 <- Cred$Dim2[which(rownames(Cred)==j)]
  
  Cred$dist_h[which(rownames(Cred)==j)] <-   
  sqrt((h_c.x -x2)^2  + (h_c.y -y2)^2)
  
  Cred$dist_l[which(rownames(Cred)==j)] <-   
    sqrt((l_c.x -x2)^2  + (l_c.y -y2)^2)
  
  
  
  #which(min(c( Cred$dist_h[which(rownames(Cred)==j)],   Cred$dist_l[which(rownames(Cred)==j)] )))
  
  if(Cred$dist_h[which(rownames(Cred)==j)] <  Cred$dist_l[which(rownames(Cred)==j)]){
    
    Cred$closest_unit[which(rownames(Cred)==j)] <- "dist_h"
  }else{
    Cred$closest_unit[which(rownames(Cred)==j)] <- "dist_l"
  }
  
  }#end of j
  
  Cred$rank <- 
    
    rank(-Cred$cor)
  
  
  env_correlations[[i]] <- Cred
  

   }#end of i
    
  
  for (j in names(env_correlations[[i]])){
  
    nrow(  
      
      env_correlations[[i]])
    
  rownames( env_correlations[[i]])[which( env_correlations[[i]]$rank %in% 1:5)]
                 
                 }
  

  a <- 
  as.list(envfit_output[[1]]$vectors)
  
  b <- 
    as.list(envfit_output[[2]]$vectors)
  

  test <- 
 rbind(
  data.frame(p=a$pvals, 
             mech=  a$r,
             type=names( envfit_output)[1]),
  data.frame(p=b$pvals, 
             mech=  b$r,
             type=names( envfit_output)[2]))
  
  test$fdr <- p.adjust(test$p,method="bonferroni")
  
  nrow(test)
  length(which(
  test$fdr<=0.05))
  
  # title here 
  #############################################################################################################################
  #envfit analysis on kefir4all dataset
  #############################################################################################################################
  
  

    
  pcoa_data_cs <- c()
  envfit_output_cs <- c()
  i="Milk.kefir"
  for (i in c("Milk.kefir","Water.kefir")){
    
    
    t <-  dplyr::select(mech_total_resistome[which(mech_total_resistome$data_source=="This study"),], -c(Hits,X2)) %>% pivot_wider(names_from = Mechanism, values_from = nor_hits) %>% column_to_rownames("Sample")                         
    
    
    t <- merge(t,
               Citizen_Scientist_metadata_v8,
               by.x="Sample.y",
               by.y="ID",
               all.x=TRUE)
    t[is.na(t)]=0
    
    
    t[which(t$Stage=="T0"),] <- 
      t[which(t$Stage=="T0"),]  %>% 
      mutate(category_confirmed=
               # wildcards in r are .
               gsub(".WL.*","", 
                    gsub(".WG.*","",
                         gsub(".ML.*","",
                              gsub(".MG.*","",
                                   Sample.y))))) %>% 
      mutate(category_confirmed=
               gsub("GO",	"Goat",
                    gsub( "FU",	"Cow milk (whole)",
                          gsub("LO", 	"Cow milk (low fat)",
                               gsub("SK",	"Cow milk (skim)",
                                    gsub("FG",	"White sugar-dried fig",
                                         gsub("FL" , "White sugar-dried fig and fresh lemon",
                                              gsub( "AP",	"White sugar-dried apricot",
                                                    gsub( "BR",	"Brown sugar-None",
                                                          category_confirmed     
                                                          
                                                          
                                                    )))))))))
    
    
    
    t$category.y[which(t$Stage=="T0")] <- "Baseline"
    
    
  
    
    t <-
      t[which(t$category.x==i),]
    
    t_names <- t
    
    t1 <- t
    
    
    #metadata <-  dplyr::select( t1, c(category.x, kefir_type_v2,conditions,category_confirmed,`Fermentation type`    )) # Stage needs to be included in cs only dataset
    
    
    t <- t[,which(colnames(t) %in% levels(as.factor(mech_total_resistome$Mechanism)))]
   
    t <- t[,-c(which(colnames(t) %in% names(which(colSums(t)==0)) ))]
    
    
    pca <- prcomp(t,
                  scale = TRUE)
    
    
    
    biplot_data[[i]] <- 
      ggbiplot(pca,
               # groups = t_names$data_source,
               groups = t_names$conditions ,
               # groups =t_names$kefir_type_v2,
               #labels = crime$st,
               # labels.size = 4,
               var.factor = 1.4,
               ellipse = TRUE, ellipse.level = 0.5, ellipse.alpha = 0.1,
               circle = TRUE,
               varname.size = 4,
               varname.color = "black") +
      labs(fill = "Region", color = "Region")+
      theme(legend.direction = 'horizontal', legend.position = 'top')
    
    
    
    
    
    ######################################
    # envfit
    x_matrix<-vegdist(t,method = "bray")
    pcoa <- cmdscale(x_matrix, k=2)
    
    #
    data.scores = as.data.frame(pcoa)
    
    data.scores =
      merge(data.scores,
            t1,
            by=0,
            all.x=TRUE)
    
  
    en_pcoa <- vegan::envfit(pcoa, dplyr::select(t1,-c(Sample.y,data_source, kefir_category,`kefir type`, Kefir,School,type,category.x,category.y,total_fermentation_response_frequency,`Answered the getting started survey`,conditions, kefir_type_v2  )), permutations = 1000)
    
    
    en_pcoa$factors$padjust <- 
      p.adjust(en_pcoa$factors$pvals, method ="bonferroni")
    
    envfit_output_cs[[i]] <- en_pcoa
    
    #en_pcoa <- vegan::envfit(pcoa, dplyr::select(t1,-c(data_source, , Stage, Sample.y,category,, kefir_type_v2)), permutations = 1000)
    
    #en_pcoa <- vegan::envfit(pcoa ~ conditions, t1, permutations = 1000)
    
    
    #     library(MASS)
    #     t1$conditions <- as.factor(    t1$conditions)
    # lda(conditions~., data=t1)
    #     
    # str(t1
    #     )
    
    
    
    
    A <- as.list(en_pcoa$vectors)
    pvals<-as.data.frame(A$pvals)
    arrows<-as.data.frame(A$arrows*sqrt(A$r))
    C<-cbind(arrows, pvals)
    
    Cred<-subset(C,pvals<0.05)
    Cred <- cbind(Cred, Species = rownames(Cred))
    
    
    
    pcoa_data_cs[[i]] <- 
      
      ggplot(data = data.scores, aes(x = V1, y = V2)) + 
      geom_point(data = data.scores, aes(colour = conditions), size = 3, alpha = 0.5) + 
      scale_colour_manual(values = c("orange", "steelblue"))  + 
      geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
                   data =  Cred, size =1, alpha = 0.5, colour = "grey30") +
      geom_point(data = Cred, aes(x = Dim1, y = Dim2), 
                 shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
      geom_text(data = Cred, aes(x = Dim1, y = Dim2+0.04), 
                label = row.names(Cred), colour = "navy", fontface = "bold") + 
      #geom_text(data = Cred, aes(x = Dim1, y = Dim2), colour = "grey30", 
      #fontface = "bold", label = row.names(Cred)) + 
      theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
            panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
            axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
            legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
            legend.text = element_text(size = 9, colour = "grey30")) + 
      labs(colour = "Season")
    
    
    
    
    
    
  }
  
  
  #############################################################################################################################
  #lefse and deseq2 total data
  #############################################################################################################################

  # pacman::p_load(lefser,microbiomeMarker)
  # mm=c()
  # mm_deseq<-c()
  # # dir.create("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/DESeq2/")
  # i="Water.kefir"
  # deseq_output <- c()
  # for (i in c("Milk.kefir","Water.kefir")){
  # 
  # t <-  dplyr::select(mech_total_resistome, -c(Hits,X2)) %>% pivot_wider(names_from = Mechanism, values_from = nor_hits) %>% column_to_rownames("Sample")
  # t <- t[which(t$category==i),]
  # t1 <- t
  # 
  # 
  # t <- t[,which(colnames(t) %in% levels(as.factor(mech_total_resistome$Mechanism)))]
  # 
  # 
  # t[is.na(t)]=0
  # 
  # 
  # t<- t(t)
  # 
  # 
  # 
  # 
  # metadata <-  dplyr::select( t1, c( kefir_type_v2,conditions,species)) # Stage needs to be included in cs only dataset
  # 
  # 
  # 
  # s_tax_tab <- t %>%  as.data.frame()%>% rownames_to_column("species") %>%
  #   dplyr::select(species) %>%
  #   dplyr::mutate(Species = trimws(species),
  #                 spec_row = Species) %>%
  #   dplyr::select(-species) %>%
  #   tibble::column_to_rownames(var = "spec_row")
  # 
  # 
  # 
  # s_otu_tab <-  t%>%
  #   as.data.frame() %>% rownames_to_column("species") %>%
  #   #dplyr::rename("Pathways" = "#Classification") %>%
  #   tibble::column_to_rownames(var = "species")
  # 
  # 
  # s_meta_time <- data.frame(seq_id = names(s_otu_tab))
  # 
  # s_meta_time <- merge(s_meta_time, metadata, by.x="seq_id",by.y=0,all.x=TRUE, sort = F)
  # 
  # s_meta_time_rows <- s_meta_time$seq_id
  # rownames(s_meta_time) <- s_meta_time_rows
  # s_meta_time$seq_id <- NULL
  # #colnames(s_meta_time) <- "time"
  # 
  # 
  # s_meta_time$seq_id <- rownames(  s_meta_time)
  # 
  # s_meta <- data.frame(seq_id = names(s_otu_tab))
  # s_meta <- s_meta %>%
  #   dplyr::mutate(sampleNames_row = seq_id) %>%
  #   tibble::column_to_rownames(var = "sampleNames_row")
  # contrast_lists <- c()
  # 
  # pacman::p_load(phyloseq)
  # 
  # library("phyloseq")
  # # Create a phyloseq object to compare 0.8 and 24 hours
  # (ps_bracken_species <- phyloseq(sample_data(s_meta_time),
  #                                 otu_table(s_otu_tab, taxa_are_rows = TRUE),
  #                                 tax_table(as.matrix(s_tax_tab))))
  # 
  # 
  # 
  # for (j in  colnames( s_meta_time)){print(paste("begining anlaysis for", j,"-",i))
  #   if(j=="seq_id"){
  #     next
  #   }else{
  # 
  # 
  #     (ps_bracken_species_deseq <- phyloseq(sample_data(dplyr::select(s_meta_time, j, seq_id) %>%
  #                                                         dplyr::rename("group"=j)),
  #                                           otu_table(s_otu_tab+1, taxa_are_rows = TRUE),
  #                                           tax_table(as.matrix(s_tax_tab))))
  # 
  # 
  # 
  # 
  #    #  mm[[i]][[j]] <-run_lefse(ps_bracken_species, group = j)
  #    #
  #    # # mm[[i]][[j]] <-run_lefse(ps_bracken_species, group = j)
  #    #  mm[[i]][[j]]@marker_table$kefir_type=i
  #    #  mm[[i]][[j]]@marker_table$category=j
  #    #  write.csv( mm[[i]][[j]]@marker_table,paste("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/lefse/lefse_",i,"_",j,".csv",sep=""),quote = FALSE,row.names = FALSE)
  #    #
  # 
  # 
  #     diagdds = phyloseq_to_deseq2(ps_bracken_species_deseq, ~ group)
  #     pen= DESeq(diagdds, test="Wald", fitType="parametric")
  # 
  #     if(length(levels(as.factor(getElement(s_meta_time, j))))>2){
  #       
  #     groupings <- levels(as.factor(getElement(s_meta_time, j)))
  #     for (t2 in groupings){
  #       
  #       groupings <-   groupings[-c(which(groupings==t2))]
  #       
  #       for (t3 in groupings){
  #         
  #         contrast_lists <- c( contrast_lists, list(c( paste("group"),t2, t3))[[1]])
  #         
  #         
  #         res <- results(pen,contrast = c("group",t2, t3))
  #         res$kefir_type <- i
  #         res$category=j
  #         
  #         deseq_output <- capture.output(res)
  #         
  #         
  #         res$log_comparision <-  deseq_output[1]
  #         res$wald_comparision <-  deseq_output[2]
  #         mm_deseq[[i]][[j]] <- rbind(  mm_deseq[[i]][[j]], res)
  #         
  #        #here  
  #       }#end of t3
  #     }#end of t2
  #     
  # 
  #     }else{
  #       
  #       
  #       res <- results(pen) #contrast = c("group","Laboratory controlled", "Household conditions"))
  #       #note I do not need the contrast option above it will pick one at random and will be displayed in the message section at the start, furthermore the log2-fold changes represent the difference in expression between the "Control" group and the "Treatment" group, stated in the text on top of the deseq2 output.
  #       
  #       # res <- cbind(res, tax_table(ps_bracken_species_deseq))
  #       
  #       res$kefir_type <- i
  #       res$category=j
  #       
  #       mm_deseq[[i]][[j]] <- res
  #       
  #       deseq_output <- capture.output(res)
  #       
  #       
  #       mm_deseq[[i]][[j]]$log_comparision <-  deseq_output[1]
  #       mm_deseq[[i]][[j]]$wald_comparision <-  deseq_output[2]
  #     
  #     }# end of else
  #   
  # }# end of else
  #   write.csv(  mm_deseq[[i]][[j]],paste("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/DESeq2/DESeq2_",i,"_",j,".csv",sep=""),quote = FALSE,row.names = TRUE)
  #   
  #   
  # }#end of j
  # }#end of i
  # 
  # 
  # 
  
  
  
  
  #############################################################################################################################
  #lefse AND DESEQ2 cs 
  #############################################################################################################################

#   pacman::p_load(lefser,microbiomeMarker)
# 
#   mm_cs=c()
#   mm_deseq_cs=c()
#   for (i in c("Milk.kefir","Water.kefir")){
# 
# 
# 
# 
#     t <-  dplyr::select(mech_total_resistome[which(mech_total_resistome$data_source=="This study"),], -c(Hits,X2)) %>% pivot_wider(names_from = Mechanism, values_from = nor_hits) %>% column_to_rownames("Sample")
# 
# 
#     t <- merge(t %>% rownames_to_column("sample_id"),
#                Citizen_Scientist_metadata_v8,
#                by.x="Sample.y",
#                by.y="ID",
#                all.x=TRUE)
#     t[is.na(t)]=0
# 
# 
# 
#     t[which(t$Stage=="T0"),] <-
#       t[which(t$Stage=="T0"),]  %>%
#       mutate(category_confirmed=
#                # wildcards in r are .
#                gsub(".WL.*","",
#                     gsub(".WG.*","",
#                          gsub(".ML.*","",
#                               gsub(".MG.*","",
#                                    Sample.y))))) %>%
#       mutate(category_confirmed=
#                gsub("GO",	"Goat",
#                     gsub( "FU",	"Cow milk (whole)",
#                           gsub("LO", 	"Cow milk (low fat)",
#                                gsub("SK",	"Cow milk (skim)",
#                                     gsub("FG",	"White sugar-dried fig",
#                                          gsub("FL" , "White sugar-dried fig and fresh lemon",
#                                               gsub( "AP",	"White sugar-dried apricot",
#                                                     gsub( "BR",	"Brown sugar-None",
#                                                           category_confirmed
# 
# 
#                                                     )))))))))
# 
# 
# 
#     t$category.y[which(t$Stage=="T0")] <- "Baseline"
# 
# 
#     t <- t[which(t$category.x==i),]
#     t1 <- t
# 
#     t <- t[,which(colnames(t) %in% c("sample_id",levels(as.factor(mech_total_resistome$Mechanism))))]
# 
#     t[is.na(t)]=0
# 
# 
#     t<- t(t %>% remove_rownames()%>%  column_to_rownames("sample_id"))
# 
# 
#     metadata=dplyr::select(t1,-c(Sample.y,data_source, kefir_category,`kefir type`,conditions, Kefir,School,type,category.x,total_fermentation_response_frequency,`Answered the getting started survey`,`Fermentation type`, total_fermentation_response_rank  ))
# 
# 
#    # metadata <-  dplyr::select( t1, c( kefir_type_v2,conditions)) # Stage needs to be included in cs only dataset
# 
#     metadata <-   metadata  %>% remove_rownames()%>%  column_to_rownames("sample_id")
# 
# 
#     metadata <-    metadata [, -c(which(colnames(   metadata) %in% levels(as.factor(mech_total_resistome$Mechanism)))) ]
# 
#     s_tax_tab <- t %>%  as.data.frame()%>% rownames_to_column("species") %>%
#       dplyr::select(species) %>%
#       dplyr::mutate(Species = trimws(species),
#                     spec_row = Species) %>%
#       dplyr::select(-species) %>%
#       tibble::column_to_rownames(var = "spec_row")
# 
# 
# 
#     s_otu_tab <-  t%>%
#       as.data.frame() %>% rownames_to_column("species") %>%
#       #dplyr::rename("Pathways" = "#Classification") %>%
#       tibble::column_to_rownames(var = "species")
# 
# 
#     s_meta_time <- data.frame(seq_id = names(s_otu_tab))
# 
#     s_meta_time <- merge(s_meta_time, metadata, by.x="seq_id",by.y=0,all.x=TRUE, sort = F)
# 
#     s_meta_time_rows <- s_meta_time$seq_id
#     rownames(s_meta_time) <- s_meta_time_rows
#     s_meta_time$seq_id <- NULL
#     #colnames(s_meta_time) <- "time"
# 
# 
#     s_meta_time$seq_id <- rownames(  s_meta_time)
# 
#     s_meta <- data.frame(seq_id = names(s_otu_tab))
#     s_meta <- s_meta %>%
#       dplyr::mutate(sampleNames_row = seq_id) %>%
#       tibble::column_to_rownames(var = "sampleNames_row")
# 
# 
#     pacman::p_load(phyloseq)
# 
#     library("phyloseq")
#     # Create a phyloseq object to compare 0.8 and 24 hours
#     (ps_bracken_species <- phyloseq(sample_data(s_meta_time),
#                                     otu_table(s_otu_tab, taxa_are_rows = TRUE),
#                                     tax_table(as.matrix(s_tax_tab))))
# 
# 
# 
# 
# 
#     for (j in  colnames( s_meta_time)){
#       if(j=="seq_id"){
#         next
#       }else{
# 
# 
#         (ps_bracken_species_deseq <- phyloseq(sample_data(dplyr::select(s_meta_time, j, seq_id) %>%
#                                                                                                                     dplyr::rename("group"=j)),
#                                                                                                       otu_table(s_otu_tab+1, taxa_are_rows = TRUE),
#                                                                                                       tax_table(as.matrix(s_tax_tab))))
# 
# 
#         diagdds = phyloseq_to_deseq2(ps_bracken_species_deseq, ~ group)
#             pen= DESeq(diagdds, test="Wald", fitType="parametric")
# 
# 
# 
#             if(length(levels(as.factor(getElement(s_meta_time, j))))>2){
# 
#               res <- results(pen)
# 
#            # res <- results(pen, contrast = c("group",levels(as.factor(getElement(s_meta_time, j)))))
# 
#             contrast_lists <- c()
#           compar_vec <- levels(as.factor(getElement(s_meta_time, j)))
#             for (k in levels(as.factor(getElement(s_meta_time, j)))){
# 
#               for (k_1 in compar_vec[-c(which( compar_vec==k))]){
# 
#                 # if(k==k_1){
#                 #   next
#                 # }else{
# 
#                 #}
#                 #contrast_lists <- c( contrast_lists,paste("'group'",",","'", k,"'",",","'",k_1,"'",sep=""))
# 
# 
#                 #contrast_lists <- c( contrast_lists,paste("group", k,k_1,sep=" "))
#                 contrast_lists <- c( contrast_lists, list(c( paste("group"),k, k_1)))
# 
# 
#                res <-  results(pen,list(c( paste("group"),k, k_1))[[1]])
# 
# 
#                res$kefir_type <- i
#                res$category=j
# 
#                deseq_output <- capture.output(res)
# 
# 
#                res$log_comparision <-  deseq_output[1]
#                res$wald_comparision <-  deseq_output[2]
# 
# 
# 
#                mm_deseq_cs[[i]][[j]] <- rbind(mm_deseq_cs[[i]][[j]],
#                                               as.data.frame(res))
# 
#            # contrast_lists <- paste(levels(as.factor(getElement(s_meta_time, j)))
# 
#             }# end of k_1
#             }#end of k
# 
# 
# 
#           write.csv(  mm_deseq_cs[[i]][[j]],paste("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/DESeq2/DESeq2_cs_",i,"_",j,".csv",sep=""),quote = FALSE,row.names = TRUE)
#           #
# 
# 
# 
# 
#             }else{
# 
# 
#             res <- results(pen) #contrast = c("group","Laboratory controlled", "Household conditions"))
#             #note I do not need the contrast option above it will pick one at random and will be displayed in the message section at the start, furthermore the log2-fold changes represent the difference in expression between the "Control" group and the "Treatment" group, stated in the text on top of the deseq2 output.
# 
#            # res <- cbind(res, tax_table(ps_bracken_species_deseq))
# 
#             res$kefir_type <- i
#             res$category=j
# 
#             mm_deseq_cs[[i]][[j]] <- res
# 
#             deseq_output <- capture.output(res)
# 
# 
#             mm_deseq_cs[[i]][[j]]$log_comparision <-  deseq_output[1]
#             mm_deseq_cs[[i]][[j]]$wald_comparision <-  deseq_output[2]
# 
#             write.csv(  mm_deseq_cs[[i]][[j]],paste("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/DESeq2/DESeq2_cs_",i,"_",j,".csv",sep=""),quote = FALSE,row.names = TRUE)
#             #
# 
#         # mm_cs[[i]][[j]] <-run_lefse(ps_bracken_species, group = j)
#         # mm_cs[[i]][[j]]@marker_table$kefir_type=i
#         # mm_cs[[i]][[j]]@marker_table$category=j
#         # write.csv( mm_cs[[i]][[j]]@marker_table,paste("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/lefse/lefse_cs_",i,"_",j,".csv",sep=""),quote = FALSE,row.names = FALSE)
# 
#             #write.csv( mm_cs[[i]][[j]]@marker_table,paste("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/lefse/lefse_cs_",i,"_",j,".csv",sep=""),quote = FALSE,row.names = FALSE)
# 
# 
#       }#end of else
# 
#     }#end of else for seq_id
#   }#end of j
# } #end of i

  
  
  #############################################################################################################################
  # import the lefse data
  #############################################################################################################################
  
 #  temp <- c()
 # 
 #   setwd("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/lefse/")
 #  
 #  # temp[["mech"]] = list.files(pattern="*mech.tsv", recursive = TRUE, full.names = TRUE)
 #  # # #temp[["class"]] = list.files(pattern="*class.tsv", recursive = TRUE, full.names = TRUE)
 #  # # #temp[["group"]] = list.files(pattern="*group.tsv", recursive = TRUE, full.names = TRUE)
 # temp[["lefse_cs"]] = list.files(pattern="*_cs", recursive = TRUE, full.names = TRUE)
 # temp[["lefse_total"]] = list.files(pattern=".csv", recursive = TRUE, full.names = TRUE)
 # temp[["lefse_total"]] <-  temp[["lefse_total"]][-c(grep("_cs",temp[["lefse_total"]]))]
 # 
 #  
 #  lefse <- c()
 #  lefse_total <- c()
 #  # restistome <- c()
 #   for (i in names(temp)){
 #  # 
 #  # 
 #  # 
 #  lefse[[i]] = lapply( temp[[i]],read_csv)
 #  # 
 #  # 
 #  # 
 #  
 #  names(  lefse[[i]]) <-  gsub("./lefse_|.csv","",temp[[i]])
 #  
 #   
 #  
 #  
 #  
 #  
 #  
 #  
 #  
 # 
 #  
 #  
 #  for (j in names(lefse[[i]])){
 #  
 #    
 #    
 #    lefse_total <- rbind( lefse_total, lefse[[i]][[j]])
 #                          
 #  
 #    # lefse[[i]][[j]] %>% 
 #    # ggplot(aes(x=ef_lda,y=feature,fill=enrich_group))+
 #    # geom_col()+
 #    # facet_wrap(~metadata,scales = "free")
 #  
 #  
 #  
 #  }
 #   
 #  
 #   }
 #  
  # lefse_total$padj_vs <- p.adjust(lefse_total$pvalue, method = "bonferroni")
  # 
  # 
  # length(lefse_total$padj_vs)
  # 
  # 
  # length(lefse_total$padj_vs[which(lefse_total$padj_vs <= 0.05)])
  # 
  # 
  # View(lefse_total[which(lefse_total$padj_vs <= 0.05),])
  # 
  # 
  # 
  # 
  # lefse_total <-
  # lefse_total[which(lefse_total$padj_vs <= 0.05),]
  # 
  # 
  
  
  #############################################################################################################################
  # import the deseq data
  #############################################################################################################################
  
  

   temp <- c()

    setwd("Q:/H2020 Master/Citizen Science Project/Results/05_short_read_functional_profiling/05_resistome/DESeq2/")

   # temp[["mech"]] = list.files(pattern="*mech.tsv", recursive = TRUE, full.names = TRUE)
   # # #temp[["class"]] = list.files(pattern="*class.tsv", recursive = TRUE, full.names = TRUE)
   # # #temp[["group"]] = list.files(pattern="*group.tsv", recursive = TRUE, full.names = TRUE)
  temp[["deseq_cs"]] = list.files(pattern="*_cs", recursive = TRUE, full.names = TRUE)
  temp[["deseq_total"]] = list.files(pattern=".csv", recursive = TRUE, full.names = TRUE)
  temp[["deseq_total"]] <-  temp[["deseq_total"]][-c(grep("_cs",temp[["deseq_total"]]))]


  
   deseq <- c()
   deseq_total <- c()
   # restistome <- c()
    for (i in names(temp)){
   #
   #
   #
   deseq[[i]] = lapply( temp[[i]],read_csv)
   #
   #
   #

   names(  deseq[[i]]) <-  gsub("./DESeq2_|.csv","",temp[[i]])


   for (j in names(deseq[[i]])){


     deseq[[i]][[j]]$log_comparision <- gsub(	
       ".*group ",
       "",deseq[[i]][[j]]$log_comparision)
     
     
     deseq[[i]][[j]]$deseq_group <- sign(deseq[[i]][[j]]$log2FoldChange)
     
     
     for ( k in levels(as.factor(deseq[[i]][[j]]$log_comparision))){
    
     
    compare <-  str_split_fixed(deseq[[i]][[j]]$log_comparision[which(deseq[[i]][[j]]$log_comparision==k)]," vs ", 2)[1,1]
     contrast <- str_split_fixed(deseq[[i]][[j]]$log_comparision[which(deseq[[i]][[j]]$log_comparision==k)]," vs ", 2)[1,2]

     

     deseq[[i]][[j]]$deseq_group[which(deseq[[i]][[j]]$log_comparision==k)] <-   gsub("-1",contrast,deseq[[i]][[j]]$deseq_group[which(deseq[[i]][[j]]$log_comparision==k)])
     deseq[[i]][[j]]$deseq_group[which(deseq[[i]][[j]]$log_comparision==k)] <-   gsub("1",compare,deseq[[i]][[j]]$deseq_group[which(deseq[[i]][[j]]$log_comparision==k)])
     
      



     # lefse[[i]][[j]] %>%
     # ggplot(aes(x=ef_lda,y=feature,fill=enrich_group))+
     # geom_col()+
     # facet_wrap(~metadata,scales = "free")



     } # end of k
     
     deseq[[i]][[j]]$analysis_group <- i
     deseq_total <- rbind( deseq_total, deseq[[i]][[j]])
   
   }# end of j
    } # end of i
   
   
   deseq_total <- 
   deseq_total[which( deseq_total$padj<=0.05),]
  

   deseq_total$envfit_overlap <- NA
   deseq_total$log_comparision <- gsub(	
     ".*group ",
     "",deseq_total$log_comparision)
     
   deseq_total$mech <- gsub("[0-9]+\\.?[0-9]+", '',deseq_total$...1)
   
   deseq_total <- 
   deseq_total[-c(which(!( deseq_total$log_comparision %in% c("Lacticaseibacillus paracasei vs Lentilactobacillus hilgardii",
                                            
                                            "Lacticaseibacillus paracasei vs Zymomonas mobilis",
                                            
                                            "Lentilactobacillus hilgardii vs Zymomonas mobilis",
                                            
                                            
                                            "Lactobacillus kefiranofaciens vs Lactococcus lactis",
                                            
                                            "Lactobacillus kefiranofaciens vs Lactococcus cremoris",
                                            
                                            "Lactobacillus helveticus vs Lactobacillus kefiranofaciens",
                          
                                            "Lactobacillus helveticus vs Lactococcus lactis",
                                            
                                            "Lactococcus cremoris vs Lactococcus lactis",
                                            "Lactobacillus helveticus vs Lactococcus cremoris")
                                            
                                           
                                           ) & deseq_total$category=="species") )    ,]
   
  
   
   deseq_total$mech <- 
   gsub("[0-9]|\\.","",deseq_total$mech)
   
   #############################################################################################################################
   # What species features are differential abundant when looking at all species comparisons
   #############################################################################################################################
   
   
   
   data <- 
   deseq_total[which(deseq_total$category %in% c('species') &
                       deseq_total$analysis_group=="deseq_total"),]
   

   
   data$n_occurences <- paste(data$deseq_group,"_", data$mech,"_",data$kefir_type,sep="")
     

   look_up_matrix <- as.data.frame(     table( data$n_occurences ))
   
   look_up_matrix <-    as.data.frame(cbind(look_up_matrix, str_split_fixed(   look_up_matrix$Var1,"_", 3)))
   look_up_matrix <- 
   look_up_matrix[-c(
   which(look_up_matrix$`3`=="Water.kefir" &    look_up_matrix$Freq<2)),]
   
   look_up_matrix <- 
     look_up_matrix[-c(
       which(look_up_matrix$`3`=="Milk.kefir" &    look_up_matrix$Freq<3)),]
   
   data <- 
   data[
   which(data$n_occurences %in% look_up_matrix$Var1),]
   # note this needs to be plotted seperately work here tomorrow
   
   
   tester <- c()
   for (species in levels(as.factor(data$deseq_group))){
     
     
     tester <- data.frame(unit=gsub("\\.","",levels(as.factor(data$mech[which(data$deseq_group==species)]))))
     
     print(paste(species, "contains", nrow(tester), "differentially abundant features which includes",paste(tester$unit, collapse=", "),sep=" "))
     
   
   }
   
   
   
   
   
   
   tester <- c()

   
   for (kefir in c("Milk.kefir","Water.kefir")){
     
     
     env <- 
       as.list( envfit_output[[kefir]]$vectors)  
     top_five <-  names(which( rank(-env$r)<=5))
     print(paste("just a reminder the number of significant AMR classes for", kefir, "microbiome is", length( env$pvals_fdr[which(    env$pvals_fdr <=0.05)])))
     
     env_names <- 
       names( env$pvals_fdr[which(    env$pvals_fdr <=0.05)])
     
  
     
     
   for (condition in    
        
        levels(as.factor( deseq_total$deseq_group[which(deseq_total$category %in% c('conditions')  &
                                                        deseq_total$analysis_group=="deseq_total" &
                                                  deseq_total$kefir_type ==kefir)
                                                  ] ))
   ){
     
     
     tester <- data.frame(unit=
     deseq_total$mech[which(deseq_total$category %in% c('conditions')  &
                                     deseq_total$analysis_group=="deseq_total" &
                                     deseq_total$deseq_group==condition&
                        deseq_total$kefir_type ==kefir)])
     
     

     
    

     #print(paste("Within the", gsub("\\."," ", tolower(kefir)),"microbiome the fermentation paramater", gsub("\\."," ", tolower(condition)), "contains", nrow(tester), "differentially abundant features which includes",paste(tester$unit, collapse=", "),sep=" "))
     
     
     # work here 
     print(paste("Within the", gsub("\\."," ", tolower(kefir)),"microbiome the fermentation paramater", 
                 gsub("\\."," ", tolower(condition)), "contains", nrow(tester), "differentially abundant features which includes",paste(tester$unit, collapse=", "),
                 "of these differential abundant features",        
                 length(tester$unit [which(tester$unit %in% env_names)]), "are considered to be classes that contributed to differences across",
                 gsub("\\."," ", tolower(kefir)), "metagenomes, and includes", paste(tester$unit [which(tester$unit %in% top_five)],collapse=", "),
                 
                 "which represent ", length(tester$unit [which(tester$unit %in% top_five)]),
                 "AMR classes that mosr contributed to differences across",
                 gsub("\\."," ", tolower(kefir)), "metagenome." ,sep=" "))
     
     
     deseq_total$envfit_overlap[which(deseq_total$category %in% c('conditions')  &
                                        deseq_total$analysis_group=="deseq_total" &
                                        deseq_total$deseq_group==condition &
                                        deseq_total$kefir_type ==kefir &
                                        deseq_total$mech  %in% env_names)
                                
     ] <- "*"
     
     
   }
   
     
     

     
     
     
   }
      
  
   env_names <- c()
     env <- c()
   for (kefir in c("Milk.kefir","Water.kefir")){
     

     env <- 
       as.list( envfit_output[[kefir]]$vectors)  
   top_five <-  names(which( rank(-env$r)<=5))
     print(paste("just a reminder the number of significant AMR classes for", kefir, "microbiome is", length( env$pvals_fdr[which(    env$pvals_fdr <=0.05)])))
     
     env_names <- 
     names( env$pvals_fdr[which(    env$pvals_fdr <=0.05)])
     
     
     
     for (condition in    
          
          levels(as.factor( deseq_total$deseq_group[which(deseq_total$category %in% c('kefir_type_v2')  &
                                                          deseq_total$analysis_group=="deseq_cs" &
                                                          deseq_total$kefir_type ==kefir)
          ] ))
     ){
       
       
       tester <- data.frame(unit=
                              deseq_total$mech[which(deseq_total$category %in% c('kefir_type_v2')  &
                                                       deseq_total$analysis_group=="deseq_cs" &
                                                       deseq_total$deseq_group==condition&
                                                       deseq_total$kefir_type ==kefir)])
       
       
       

       # work here 
       print(paste("Within the", gsub("\\."," ", tolower(kefir)),"microbiome the fermentation paramater", 
                   gsub("\\."," ", tolower(condition)), "contains", nrow(tester), "differentially abundant features which includes",paste(tester$unit, collapse=", "),
                   "of these differential abundant features",        
                   length(tester$unit [which(tester$unit %in% env_names)]), "are considered to be classes that contributed to differences across",
                   gsub("\\."," ", tolower(kefir)), "metagenomes, and includes", paste(tester$unit [which(tester$unit %in% top_five)],collapse=", "),
                   
                "which represent ", length(tester$unit [which(tester$unit %in% top_five)]),
                  "AMR classes that mosr contributed to differences across",
                   gsub("\\."," ", tolower(kefir)), "metagenome." ,sep=" "))
       
       
       
       
       deseq_total$envfit_overlap[which(deseq_total$category %in% c('kefir_type_v2')  &
                                          deseq_total$analysis_group=="deseq_cs" &
                                          deseq_total$deseq_group==condition &
                                          deseq_total$kefir_type ==kefir &
                                          deseq_total$mech  %in% env_names)
                                  
       ] <- "*"
       
       
       
       
       
       deseq_total$envfit_overlap[which(deseq_total$category %in% c('kefir_type_v2')  &
                                          deseq_total$analysis_group=="deseq_cs" &
                                          deseq_total$deseq_group==condition &
                                          deseq_total$kefir_type ==kefir &
                                          deseq_total$mech  %in% top_five)
                                  
       ] <- "**"
       
       
       
       
       
     }
     
     # pvals<-as.data.frame(env$pvals)
     # pvals_fdr<-as.data.frame(env$pvals_fdr)
     # arrows<-as.data.frame(env$arrows*sqrt(env$r))
     # C<- cbind(arrows, pvals,pvals_fdr, as.data.frame(env$r))
     # 
     # Cred<-subset(C,pvals_fdr<0.05)
     # Cred <- cbind(Cred, Species = rownames(Cred))
     # 
     # colnames(Cred)[which(colnames(Cred)=="env$r")] <- "cor"
     # 
     # Cred<- Cred %>% dplyr::select(-c(Dim1,Dim2,))
     # 
     # Cred$deseq_category <- condition
     # Cred$deseq_category <- kefir
     # Cred$differntial_abundance <- NA
     # Cred$differntial_abundance[which(Cred$Species %in% env_names)] <- "*"
     
   }
     
     
     
     
     
#View(deseq_total[which(deseq_total$envfit_overlap=="*"),])
   
   
   
   #############################################################################################################################
   # plot deseq categories kefir_type_v2 and conditions 
   #############################################################################################################################
   
   
   pd <- 
   deseq_total[which(deseq_total$category %in% c('kefir_type_v2') &
                       deseq_total$analysis_group=="deseq_total"),]  %>% 
     mutate(kefir_type=gsub("\\."," ",kefir_type)) %>% 
     mutate(deseq_group=gsub("\\."," ",deseq_group)) %>% 
     mutate(category=str_to_title(category)) %>% 
     mutate(sample_unit=paste(kefir_type," - ",gsub("_type_v2"," type", category),sep="")) %>% 
     
     
     #  mutate(sample_unit=paste(...1, kefir_type)) %>% 
     #mutate(log2FoldChange = abs(as.numeric(reorder_within(log2FoldChange,kefir_type,category)))) %>%
     ggplot(aes(y=reorder_within(...1,log2FoldChange,sample_unit),x= log2FoldChange,fill=deseq_group))+ #fill=enrich_group
     # ggplot(aes(y= reorder(log2FoldChange,category),x= ...1,fill=deseq_group))+ #fill=enrich_group
     geom_col()+
     scale_y_reordered() +
     facet_wrap(~sample_unit,scales = "free")+
     #coord_flip() +
     theme_bw()+
     labs(fill="Covariate",y="Feature")+
     theme(
       axis.text.y =element_text(size=10.5),
       axis.text.x = element_blank(),
       strip.background = element_rect(
         color="black", fill="white"),
       strip.text = element_text(size=15.5),
       #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
       axis.title.y = element_text( size=35, face="bold",hjust = 0.5, vjust = 1.5),
       legend.key.size = unit(1, 'cm'), #change legend key size
       legend.key.height = unit(1, 'cm'), #change legend key height
       legend.key.width = unit(1, 'cm'), #change legend key width
       legend.direction = "vertical",
       # legend.position = "top",
       #legend.box = "horizontal",
       legend.box = "vertical",
       # legend.text = element_text(size=30),
       legend.title = element_text(size=17.5),
       legend.text=element_text(size=17.5),
       legend.position = c(0.12, 0.8))
   
   pe <- 
   
   deseq_total[which(deseq_total$category %in% c('conditions') &
                       deseq_total$analysis_group=="deseq_total"),]  %>% 
     mutate(kefir_type=gsub("\\."," ",kefir_type)) %>% 
     mutate(deseq_group=gsub("\\."," ",deseq_group)) %>% 
     mutate(category=str_to_title(category)) %>% 
     mutate(sample_unit=paste(kefir_type," - ",gsub("_type_v2"," type", category),sep="")) %>% 
     
     
     #  mutate(sample_unit=paste(...1, kefir_type)) %>% 
     #mutate(log2FoldChange = abs(as.numeric(reorder_within(log2FoldChange,kefir_type,category)))) %>%
     ggplot(aes(y=reorder_within(...1,log2FoldChange,sample_unit),x= log2FoldChange,fill=deseq_group))+ #fill=enrich_group
     # ggplot(aes(y= reorder(log2FoldChange,category),x= ...1,fill=deseq_group))+ #fill=enrich_group
     geom_col()+
     scale_y_reordered() +
     facet_wrap(~sample_unit,scales = "free")+
     #coord_flip() +
     theme_bw()+
     labs(fill="Fermentation parameter",y="AMR class")+
     theme(
       axis.text.y =element_text(size=10.5),
       axis.text.x = element_blank(),
       strip.background = element_rect(
         color="black", fill="white"),
       strip.text = element_text(size=15.5),
       #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
       axis.title.y = element_text( size=35, face="bold",hjust = 0.5, vjust = 1.5),
       legend.key.size = unit(1, 'cm'), #change legend key size
       legend.key.height = unit(1, 'cm'), #change legend key height
       legend.key.width = unit(1, 'cm'), #change legend key width
      # legend.direction = "horizontal",
      legend.direction = "vertical",
       # legend.position = "top",
       #legend.box = "horizontal",
      legend.box = "vertical",
       # legend.text = element_text(size=30),
       legend.title = element_text(size=17.5),
       legend.text=element_text(size=17.5),
       legend.position = c(0.12, 0.75))
   
  
   
   #############################################################################################################################
   # plot deseq category species 
   #############################################################################################################################
   
   
   
   data[which(data$category %in% c( 'species') &
                       data$analysis_group=="deseq_total"),]  %>% 
     mutate(kefir_type=gsub("\\."," ",kefir_type)) %>% 
     mutate(deseq_group=gsub("\\."," ",deseq_group)) %>% 
     mutate(category=str_to_title(category)) %>% 
     mutate(sample_unit=paste(kefir_type," - ",category,sep="")) %>% 
     
     
     #  mutate(sample_unit=paste(...1, kefir_type)) %>% 
     #mutate(log2FoldChange = abs(as.numeric(reorder_within(log2FoldChange,kefir_type,category)))) %>%
     ggplot(aes(y=reorder_within(...1,log2FoldChange,sample_unit),x= log2FoldChange,fill=deseq_group))+ #fill=enrich_group
     # ggplot(aes(y= reorder(log2FoldChange,category),x= ...1,fill=deseq_group))+ #fill=enrich_group
     geom_col()+
     scale_y_reordered() +
     facet_wrap(~sample_unit,scales = "free")+
     #coord_flip() +
     theme_bw()+
     labs(fill="Covariate",y="Feature")+
     theme(
       axis.text.y =element_text(size=10.5),
       axis.text.x = element_blank(),
       strip.background = element_rect(
         color="black", fill="white"),
       strip.text = element_text(size=15.5),
       #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
       axis.title.y = element_text( size=35, face="bold",hjust = 0.5, vjust = 1.5),
       legend.key.size = unit(1, 'cm'), #change legend key size
       legend.key.height = unit(1, 'cm'), #change legend key height
       legend.key.width = unit(1, 'cm'), #change legend key width
       legend.direction = "horizontal",
       legend.position = "top",
       legend.box = "horizontal",
       # legend.text = element_text(size=30),
       legend.title = element_text(size=17.5),
       legend.text=element_text(size=17.5))
   
   
   
   
   
  
   
   
   
  mech_total_resistome[which(mech_total_resistome$Mechanism=="MLS" &
                                mech_total_resistome$kefir_category=="Water.kefir"),] %>% 
    
    
    
    ggplot(aes(x=conditions,y=nor_hits,fill=conditions))+
  geom_boxplot()
  
  
  # lefse[[i]][[j]] %>% 
  # ggplot(aes(x=ef_lda,y=feature,fill=enrich_group))+
  # geom_col()+
  # facet_wrap(~metadata,scales = "free")
  
  
  
  
  ##############################################################################################################################
  # plot time baby
  #############################################################################################################################
  jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/CS_metagenomics/resistome_figure.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
  
  
  ggarrange(pa,
             ggarrange(   pcoa_data[[1]]$conditions,   pcoa_data[[2]]$conditions, nrow=1,ncol=2, labels=c("B.","C."),common.legend = FALSE,  font.label = list(size = 30)),
            pe,
             nrow=3,ncol=1, heights = c(9, 3,3),widths=c(10,1),labels=c("A.","","D."),  font.label = list(size = 30))#+theme(panel.spacing.x=unit(-10, "lines"))
 
  graphics.off()
  
  

  redefined_plot <- pd + theme(legend.position = "right")
  
  
  jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_metagenomics/Figures -CS_Metagenomics/supplementary_resistome_figure.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
  
  
  ggarrange(
            ggarrange(   pcoa_data[[1]]$kefir_type_v2,   pcoa_data[[2]]$kefir_type_v2, nrow=1,ncol=2, labels=c("A.", "B."),common.legend = FALSE,  font.label = list(size = 30)),
            redefined_plot,
            nrow=2,ncol=1, heights = c(5,5),labels=c("","C."),  font.label = list(size = 30))#+theme(panel.spacing.x=unit(-10, "lines"))
  
  graphics.off()
  
  
  

  #############################################################################################################################
  # junk code
  #############################################################################################################################
  
  


metadata <-  dplyr::select( t1, c(,category, kefir_type_v2,conditions,Stage)) # Stage needs to be included in cs only dataset

  
  s_tax_tab <-  metacache[,which( colnames(metacache) %in%      sub_metadata$merge_column[-c(which(is.na(
    getElement(sub_metadata,i))))])] %>% 
    
  
  
    
    
    ef.df<-as.data.frame(en_pcoa$vectors$arrows*sqrt(en_pcoa$vectors$r))
    ef.df$species<-rownames(ef.df)

    
  
    plot(en_pcoa)
    pcoa
    
    spp.scrs <- as.data.frame(scores(    en_pcoa, display = "vectors"))
    spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
    
    
    
    
    
    
    
    Ordination.model2 <- rda(t ~ conditions, data=t_names)
    summary(Ordination.model2)
    plot(Ordination.model2)
    
    anova.cca(    Ordination.model2, step = 1000)
    anova.cca(    Ordination.model2, step = 1000, by = "term")
    
    
    kruskal.test(as.matrix(t) ~ conditions, data=t_names)
    pairwise.wilcox.test(meta$chao, meta$AgeGroup, p.adjust.method="fdr")
    
    
   # The included environmental variables explain 73.41% of the variation in fish community composition across sites."

    
    ord <- cca(t)
    fit <- envfit(ord ~ conditions , metadata, perm = 1000)


    #devtools::install_github("gavinsimpson/ggvegan")
    #BiocManager::install("ggbio")
    

    en_pcoa$vectors
    chem.scores.envfit <- as.data.frame(scores(    en_pcoa, display = "vectors"))

    
    pc
    set.seed(123)
    en_pcoa <- vegan::envfit(pcoa,metadata, permutations = 1000, na.rm = TRUE, choices = c(1,2))

    en_pcoa$vectors
    envfit(pcoa~conditions,data=metadata, permutations = 1000, na.rm = TRUE)
    
    
    spp.scrs <- as.data.frame(scores(en_pcoa, display = "vectors")) #
    plot(en_pcoa)
    
    arrow<-data.frame( en_pcoa$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
    
    ###############################
    #attempt 2
    dune.mds <- metaMDS(x_matrix, distance = "bray", autotransform = F)
    dune.mds <- metaMDS(t, distance = "bray", autotransform = T)
    
    
    metadata[ which(rownames( metadata)%in% rownames(t)),]
    
    dune.envfit <- envfit(    dune.mds,     t1, permutations = 1000, na.rm = TRUE, display = "species")
    
    envfit(    dune.mds ~ conditions,     data=metadata, permutations = 1000, na.rm = TRUE), display = "species")
    
    
    
    dune.envfit$vectors
    
    
    site.scrs <- as.data.frame(scores(dune.mds, display = "sites"))
    
    plot(    dune.envfit)
   View(metadata_adonis)
    
    
   library(glmnet)
   
   tempcv <- cv.glmnet(x=as.matrix(t), y=metadata, family="multinomial", 
                       nfolds=20, alpha=0.5)
   coefsMin <- coef(tempcv, s="lambda.min")
   
   
    
    mech="Aminocoumarins" 
    
    
    
    
    data(QuickStartExample)
    x <- QuickStartExample$x
    y <- QuickStartExample$y
    
    
    
    
  for (mech in levels(as.factor(colnames(t)))){
    
    
    dplyr::select( t, -c(data_source, `kefir type`, Stage, Sample.y,category, kefir_category, kefir_type_v2))
    
    
    t2 <- t1[,which(colnames(t1) %in% c(mech,"conditions", "kefir type" ))]
    
    
    
    ggbetweenstats(data =t2,
                   x=conditions, 
                   y=Trimethoprim,
                   #title=type,
                   type = "nonparametric", # A
                   ggsignif.args    = list(textsize = 2, tip_length = 0.01)
                   
    )
    
    
    res<- cor.test(as.numeric(  getElement(t2[which(t2$conditions=="Laboratory controlled"),], mech)),
                                as.numeric(  getElement(t2[which(t2$conditions=="Household conditions"),], mech)),
                   method="kendall")
    
    
  }
   
    
    correlation_data[[i]] <- rbind(data.frame(
      group_type=as.character(i),
      mech_type=as.character(mech),
      p_value=as.numeric(res$p.value),
      r=as.numeric(res$estimate)),
      correlation_data[[i]])
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  group.labs <-
    iris.gg$data |>
    summarise(xvar = mean(xvar),
              yvar = mean(yvar), .by = groups)
  
  group.labs
  #>       groups   xvar   yvar
  #> 1     setosa -2.217 -0.288
  #> 2 versicolor  0.495  0.548
  #> 3  virginica  1.723 -0.260
  Now, just use geom_label to draw labels for the groups.
  
  iris.gg + geom_label(data = group.labs,
                       aes(x = xvar, y=yvar, label=groups),
                       size = 5) +
    theme(legend.position = "none")
  
  
  # results <- princomp(t)
  # 
  # biplot(results)
  

  
  # library(factoextra)
  # 
  # fviz_pca_biplot(pca),
  #                 label="var",
  #                 col.ind = "cos2",
  #                 col.var = "black",
  #                 gradient.cols = c("blue","green","red"))
  
  
  
  data(iris)
  
  
  library(ggplot2)
  library(ggbiplot)
  library(dplyr)
  library(corrplot)
  
  ggbiplot(pca,
         #  labels = crime$st ,
           circle = TRUE,
           varname.size = 4,
           varname.color = "red")

  
  

  ggbiplot(pca,
           groups = t$category,
           #labels = crime$st,
          # labels.size = 4,
           var.factor = 1.4,
           ellipse = TRUE, ellipse.level = 0.5, ellipse.alpha = 0.1,
           circle = TRUE,
           varname.size = 4,
           varname.color = "black") +
    labs(fill = "Region", color = "Region") +
    theme(legend.direction = 'horizontal', legend.position = 'top')
                  habillage = rownames(t))
  
