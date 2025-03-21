
########################################################################################################################
#Libraries used
########################################################################################################################

.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")
pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis)
pacman::p_load(readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)
library(vegan)
library(ggplot2)
library(grid)


########################################################################################################################
#Import metacache data
########################################################################################################################
metacache <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_species_profile.csv")
metacache_strain <-  read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_strain_profile.csv")

metacache$sample_id <- gsub("-","_",metacache$sample_id)
metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]


########################################################################################################################
#Import metadata
########################################################################################################################

global_mk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_milk_kefir_metadata_v1.csv")
global_wk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_water_kefir_metadata_v1.csv")
Citizen_Scientist_metadata_v8 <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/Citizen Scientist metadata_v8.csv")

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/kefir4all_sample_metadata_v2.csv")
kefir4all_metadata$merge_column <-  gsub("_host_removed_R..fastq.gz","",kefir4all_metadata$merge_column)
########################################################################################################################
# Merge metadata into one
########################################################################################################################
total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`))



########################################################################################################################
#Import prevalence data 
########################################################################################################################
metacache_prevalence<-c()
metacache_prevalence[["Milk.kefir"]] <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/prevalence/milk_metacache_prevalence.csv")


metacache_prevalence[["Water.kefir"]] <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/prevalence/water_metacache_prevalence.csv")


########################################################################################################################
# Look at microbes in the medium controls 
########################################################################################################################

mediums_controls_type <- c()
# total data
medium_controls <-   metacache[which(metacache$sample_id %in% 
                                       total_metadata$merge_column[which(
                                         total_metadata$`kefir type`=="Medium control")]),]

# get metadata 

medium_metadata <- total_metadata[which(
  total_metadata$`kefir type`=="Medium control"),]



medium_metadata$type <- "Water.kefir"

medium_metadata$type[which(medium_metadata$data_source=="Walsh et al 2023")] <- "Milk.kefir"

medium_metadata$type[grep("MLC",medium_metadata$merge_column)] <- "Milk.kefir"

mediums_controls_type[["Water.kefir"]] <- medium_controls$sample_id



# seperate in milk and water kefir medium

mediums_controls_type[["Water.kefir"]] <- metacache[which(metacache$sample_id %in% 
                                                            medium_metadata$merge_column[which(
                                                              medium_metadata$type=="Water.kefir")]),]



mediums_controls_type[["Milk.kefir"]] <- metacache[which(metacache$sample_id %in% 
                                                           medium_metadata$merge_column[which(
                                                             medium_metadata$type=="Milk.kefir")]),]



n1_species <- c()
maxab_species <- c()
i <-  "Milk.kefir"
medium_prevalence <- c()

species_medium_prevalence <- c()
environment_medium_prevalence <- c()
#filter medium control and determine both the prevalent microbes that may be medium derived and the environmentally microbes that may be medium derived.

for (i in names(mediums_controls_type)){
  
  mediums_controls_type[[i]] <- mediums_controls_type[[i]] %>% column_to_rownames("sample_id")
  mediums_controls_type[[i]][is.na(mediums_controls_type[[i]])] <-0
  maxab_species <- apply(mediums_controls_type[[i]],2, max, na.rm=TRUE)
  
  n1_species <-names(which(maxab_species > .1))
  
  
  mediums_controls_type[[i]] <- mediums_controls_type[[i]] [,which(colnames(mediums_controls_type[[i]] ) %in%  n1_species )]
  
  
  medium_prevalence[[i]] <- data.frame(specie=as.character(),
                                       sum_specie=as.numeric())
  species <- "Thermus thermophilus"  
  species <-"Lactobacillus kefiranofaciens"
  for (species in colnames(mediums_controls_type[[i]])){
    specie <- species
    
    sum_specie <-  length(which(mediums_controls_type[[i]][,species] >.1))
    
    data.prevalence <- cbind(specie,sum_specie)
    medium_prevalence[[i]] <- merge(data.prevalence,medium_prevalence[[i]], all=TRUE)
    
  }
  
  medium_prevalence[[i]]<- medium_prevalence[[i]][which(as.numeric(medium_prevalence[[i]]$sum_specie)> nrow(mediums_controls_type[[i]])*.1),]
  medium_prevalence[[i]]$kefir_type<- i
  
  mediums_controls_type[[i]]  <- mediums_controls_type[[i]][,which(colnames(mediums_controls_type[[i]]) %in%    medium_prevalence[[i]]$specie)]
  
  
  species_medium_prevalence[[i]] <- 
    
    metacache_prevalence[[i]]$specie[
      metacache_prevalence[[i]]$specie %in% 
        colnames(mediums_controls_type[[i]])]
  
  
  
  
}

#i <- "Water.kefir"


# environment_medium_prevalence[[i]] <- environmental_microbes$species[which(environmental_microbes$type==c("Occurs in both", i))][which(environmental_microbes$species[which(environmental_microbes$type==c("occurs in both", i))]%in% 
#   colnames(mediums_controls_type[[i]]))]
# work here tomorrow


#}









########################################################################################################################
#Get data from the global_MK study
########################################################################################################################

#metacache <- metacache[which(gsub("_S.*","", metacache$sample_id) %in% global_mk_metadata$Sample),]

########################################################################################################################
#PCOA used
########################################################################################################################


metacache [is.na(metacache )] <- 0

#metacache <- metacache[which(metacache$sample_id %in% 
#kefir4allmetadata$merge_column[grep("ID52", kefir4allmetadata$Sample)]), ]


#maxab_species <- apply(dplyr::select(metacache, -sample_id),2, max, na.rm=TRUE)

#n1_species <-which(maxab_species < 0.1)


#metacache <- 
#metacache %>% column_to_rownames("sample_id")


########################################################################################################################
#Alpha diversity total dataset
########################################################################################################################

my.files_summary_shannon <- diversity(metacache %>% column_to_rownames("sample_id") )
my.files_richness <- specnumber(metacache %>% column_to_rownames("sample_id"))
my.files_evenness<- my.files_summary_shannon/log(my.files_richness)
my.files_beta<- vegdist(metacache %>% column_to_rownames("sample_id"), method = "bray")

my.files_summary<- cbind(shannon = my.files_summary_shannon, richness = my.files_richness, pielou = my.files_evenness,site=metacache$sample_id)


my.files_summary<- as.data.frame(my.files_summary)
my.files_summary <- merge(my.files_summary,total_metadata,by.x="site",by.y="merge_column",all.x=TRUE)



my.files_summary$`kefir type`[which(is.na(my.files_summary$`kefir type`))] <- "ML"
my.files_summary$data_source[which(is.na(my.files_summary$`kefir type`))] <- "Walsh et al 2022" 



my.files_summary <- my.files_summary[-c(which(my.files_summary$site %in% total_metadata$merge_column[which(total_metadata$`kefir type` %in% c("Extraction control", "Medium control"))])),]

my.files_summary$data_source[which(is.na(my.files_summary$data_source))] <-  "Walsh et al 2023" 

my.files_summary$conditions <-  "Household conditions"

my.files_summary$conditions[which(my.files_summary$data_source %in% c("Walsh et al 2023","Mortensen et al 2023" ))] <- "Laboratory controlled"

setwd("Q:/H2020 Master/Citizen Science Project/Plots/Evolution")

#jpeg(filename='Alpha diversity_kefir_type.jpeg', width = 35*700, height=30*700,res=1700,pointsize = 15) #, width=2000, height=1950)
########################################################################################################################
#tstatistically test samples in sterile conditions (global studies) and kefir4all (non sterile conditions) uing wilcoxon 
########################################################################################################################

library(rstatix)
wilcoxon_mk <- wilcox_test(shannon ~ conditions, data = my.files_summary[which(my.files_summary$`kefir type`=="ML"),])#
wilcoxon_mk <- wilcoxon_mk %>% add_xy_position(x = "conditions", fun = "mean_se", scales = "free", step.increase = 0)
wilcoxon_mk$type <-"mk"


wilcoxon_wk <- wilcox_test(shannon ~ conditions, data = my.files_summary[which(my.files_summary$`kefir type`=="WL"),])#
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

my.files_summary$merge <- paste(gsub("ML","Milk Liquid",
                                     gsub("WL","Water Liquid",my.files_summary$`kefir type`)),"\n                   ",
                                          my.files_summary$conditions)


pb <-
  
  ggbetweenstats(data=my.files_summary %>%
         filter(`kefir type` %in% c("ML","WL")), x=merge, y=shannon, 
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


  library(ggpubr)
# pa <- 
# ggplot(my.files_summary %>%
#          filter(`kefir type` %in% c("ML","WL")), aes(x=`kefir type`, y=as.numeric(shannon),colour=conditions)) +
#   geom_boxplot() +
#   #labs(title= 'Alpha diversity of timepoints') +
#   geom_point(position=position_dodge(width=0.75),aes(group=conditions))+
#   theme_bw()+
#   xlab("Fermentation conditions")+
#   ylab("Alpha diversity values (Shannon)")+
#   labs(fill = "Timepoint")+
#   #guides(colour = guide_legend(override.aes = list(size=25)))+
#   #ylim(0,250)+
#   theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
#         legend.title = element_text( size=25, face="bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
#         axis.text.y = element_text(hjust = 1, size = 10),
#         axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
#         axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2))+
#   scale_x_discrete(labels=c('Milk Liquid', 'Water - Liquid'))+
#   stat_compare_means(label = "p.signif")
  #facet_wrap(~category.x+...2,scale="free")+
  
  #stat_pvalue_manual(   data = wilcoxon,
                        #step.increase = 0.1,
                        #label = "p_mark",
                        #bracket.nudge.y = 
                          #wilcoxon$y.position[which(wilcoxon$y.position ==max(wilcoxon$y.position))]*.5
                        #)

########################################################################################################################
#test with statitsics
########################################################################################################################
library(nortest)

hist(as.numeric(my.files_summary$shannon))
adtest <- ad.test(as.numeric(my.files_summary$shannon))
if(adtest$p.value<0.05){
  
  print(paste("ad test is less than 0.05, data is not normally distributed"))
}

#Check homogeneity of variances assumption
#library(car)
my.files_summary$shannon <- as.numeric(my.files_summary$shannon)
my.files_summary$`Sample ID` <- as.factor(my.files_summary$site )

library(car )
levene <- leveneTest(shannon ~ `Sample ID`, data = my.files_summary)

if(levene$`Pr(>F)`[1]  <0.05){
  
  print(paste("for ",file,"leveneTest is less than 0.05, data does not have equal variance."))
}




#As the p value obtained from the Shapiro-Wilk test and Levene's test is significant (p < 0.05), we conclude that the data is not normally distributed and does not have equal variance

krustal <- kruskal.test(shannon ~ `Sample ID`, data = my.files_summary)#
if(krustal$p.value  <0.05){
  
  print("for kruskal.test is less than 0.05, as the p value obtained from the Kruskal-Wallis test test is significant p < 0.05), we conclude that there are significant differences in alpha diversity values among the time point varieties.")
}else{
  print("for kruskal.test is greater than 0.05, as the p value obtained from the Kruskal-Wallis test test is not significant , we conclude that there are no significant differences in alpha diversity values among the time point varieties.")
}


########################################################################################################################
#Alpha diversity total dataset
########################################################################################################################

# ggplot(my.files_summary, aes(x=`kefir type`, y=as.numeric(shannon),fill=`kefir type`)) +
#   geom_boxplot() +
#   #labs(title= 'Alpha diversity of timepoints') +
#   geom_point()+
#   theme_bw()+
#   xlab("Time points")+
#   ylab("Alpha diversity values (Shannon)")+
#   labs(fill = "Timepoint")+
#   #guides(colour = guide_legend(override.aes = list(size=25)))+
#   #ylim(0,250)+
#   theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
#         legend.title = element_text( size=25, face="bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
#         axis.text.y = element_text(hjust = 1, size = 10),
#         axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
#         axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
#         legend.position = "none")+
#   scale_x_discrete(labels=  c('Milk - Grains', 'Milk Liquid', 'Water - Grains', 'Water - Liquid'))+
#   stat_compare_means(label = "p")+
#   stat_compare_means(label = "p.signif", method = "kruskal.test", comparisons =c('MG', 'ML', 'WG', 'WL'))
#                      #ref.group = ".all.") 
# 
# 
# 



library(rlang)
library(vctrs)
library(statsExpressions)
#install.packages("ggstatsplot", dependencies = TRUE, INSTALL_opts = '--no-lock')
library(ggstatsplot)
#graphics.off()
pa <- 

ggbetweenstats(
  data = my.files_summary,
  x = `kefir type`,
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
#ref.group = ".all.") 

  
  
  

my.files_summary_cs <- my.files_summary[which(my.files_summary$data_source=="This study"),]

ggplot(my.files_summary_cs, aes(x=`kefir type`, y=as.numeric(shannon),fill=`kefir type`)) +
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
                                                                                                                                                                                                         #Look into persitance across environmental micorbes 
#######################################################################################################################




###########################################################################################################################################################


# Sources of environmental microbes occurring at >.01 relative abundance

###########################################################################################################################################################
metacache <- metacache[which(metacache$sample_id %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$data_source=="This study")]),]



metacache <-  metacache %>% column_to_rownames("sample_id")

maxab_species <- apply(metacache,2, max, na.rm=TRUE)







# all species at >.1 % RA not jst prevalent micorbes
species_data  <- metacache [,-c(which(maxab_species < .1))] %>% rownames_to_column("sample") %>% pivot_longer(!sample,names_to = "clade_name",values_to = "relative_abundance")


species_data <-  
  merge(species_data,kefir4all_metadata,by.x="sample",by.y="merge_column",all.x=TRUE)

species_data <- species_data[-c(which(duplicated(species_data))),]

species_data <-  species_data[-c(grep("WLC",species_data$sample)),]

species_data <-  species_data[-c(grep("MLC",species_data$sample)),]


# species_data <- 
# species_data[-c(which(is.na(species_data$type))),]



species_data$type <- NA
#species_data$type <- gsub("milk","Milk kefir",species_data$type)
#species_data$type <- gsub("water","Water kefir",species_data$type)


#table(dplyr::select(species_data,clade_name,type))



species_data$type[which(species_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                         kefir4all_metadata$Stage!="T0")])] <- "Water.kefir"

species_data$type[which(species_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                         kefir4all_metadata$Stage!="T0")])] <- "Milk.kefir"

species_data$type[which(species_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                         kefir4all_metadata$Stage=="T0")])] <- "Water.kefir.T0"

species_data$type[which(species_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                         kefir4all_metadata$Stage=="T0")])] <- "Milk.kefir.T0"


# remove extraction controls
# species_data <- species_data[-c(which(is.na(
#    species_data$type))),]


compositional_differences_.1 <- data.frame(species=as.character( levels(as.factor(species_data$clade_name))),
                                           `Milk kefir`=as.character(NA),
                                           `Milk kefir T0`=as.character(NA),
                                           `Water kefir`=as.character(NA),
                                           `Water kefir T0`=as.character(NA)
)

species <- "Lactococcus lactis"
type_2 <- "Water.kefir.T0"
for (species in compositional_differences_.1$species){
  
  for (type_2 in levels(as.factor(species_data$type))){
    compositional_differences_.1[which(   compositional_differences_.1$species==species),which(colnames(   compositional_differences_.1 )==type_2)] <-
      nrow(subset(species_data, subset=clade_name==species &
                    type==type_2 &
                    relative_abundance >0.1))
  }
}

species_of_interest <- 
  subset( compositional_differences_.1, subset=
            compositional_differences_.1$Milk.kefir.T0==0 &
            compositional_differences_.1$Water.kefir.T0==0 &
            compositional_differences_.1$Milk.kefir>0 &
            compositional_differences_.1$Water.kefir>0
  ) %>% dplyr::select(species) %>% 
  mutate(type="Occurs in both")


# Identify potential environmental microbes in both milk and water kefir 


# Identify potential environmental microbes in milk kefir

t<- 
  compositional_differences_.1[which(compositional_differences_.1$Milk.kefir.T0==0 & compositional_differences_.1$Milk.kefir>0) ,] %>% 
  dplyr::select(species) %>% mutate(type="Milk.kefir")

t <- t[-c(which(t$species %in% species_of_interest$species)),]


# Identify potential environmental microbes in water kefir

t1 <-
  compositional_differences_.1[which(compositional_differences_.1$Water.kefir.T0==0 & compositional_differences_.1$Water.kefir>0) ,] %>% 
  dplyr::select(species)%>% mutate(type="Water.kefir")


t1 <- t1[-c(which(t1$species %in% species_of_interest$species)),]


species_of_interest <-  rbind(species_of_interest,t,t1)


print(paste("okay, we can see", nrow(species_of_interest), "microbes that may be environmentally derived. That's not the most interesting unless they persist. Lets look at that now"))

environmental_microbes  <- species_of_interest

#write.csv(environmental_microbes, "Q:/H2020 Master/Citizen Science Project/Results/00/environmental_microbes.txt",quote=FALSE,row.names = FALSE)
###########################################################################################################################################################


# Sources of environmental microbes occurring at >.01 relative abundance

###########################################################################################################################################################

species <- "Pseudomonas rhodesiae"  
type_2 <- "Milk.kefir" 
species <- "Leuconostoc mesenteroides"

#Pseudomonas putida look here in more details

species_of_interest_detected_in_grain <- data.frame("sample" =as.character(),
                                                    "clade_name" =as.character(),
                                                    "relative_abundance"=as.numeric(),
                                                    "Sample" =as.character(),
                                                    "kefir type"=as.character(),
                                                    "Stage"=as.character(),
                                                    "data_source"=as.character(),
                                                    "type"=as.character())

sample_of_interest <- c()
sample_of_interest_reduced <- c()
for (species in species_of_interest$species){
  
  type_2 <- species_of_interest$type[which(species_of_interest$species==species)]
  
  if(type_2== "Occurs in both"){
    type_2=c("Milk.kefir", "Water.kefir")
    grain_type <- c("MG","WG")
  }else if(type_2=="Water.kefir"){
    grain_type <- "WG"
  }else{
    grain_type <- "MG"
  }
  
  
  
  
  # if(length(levels(as.factor(subset(species_data, subset=clade_name==species &
  #           type %in% type_2 &
  #           relative_abundance >.1)$Stage)))>4){
  
  #identify what I would call persistence events where the micro appears to have been transferred to the grain, to do this I need to identify if its in the grain and then is it transfered
  
  
  for (grains in levels(as.factor(grain_type))){
    
    
    if(nrow(subset(species_data, subset=clade_name==species &
                   type == type_2 &
                   relative_abundance >.1&
                   `kefir type`==grains))==0){
      next
    }else{
      
      
      
      
      # identify the sample with the microbe of interest subset(species_data, subset=clade_name==species &
      sample_of_interest <-  subset(species_data, subset=clade_name==species &
                                      type == type_2 &
                                      relative_abundance >.1) %>% dplyr::select(sample,clade_name,relative_abundance,Sample, `kefir type`, Stage,data_source,type)
      
      
      # iDENTIFY DATAFRAMES WITH DUPLICATED SAMPLE NAMES AND COME FROM DIFFERENT STAGES
      
      if(TRUE %in% duplicated(sample_of_interest$Sample)){
        sample_of_interest <-    sample_of_interest[which(sample_of_interest$Sample %in% sample_of_interest$Sample[which(duplicated(sample_of_interest$Sample))]),]
        
        sample_of_interest_reduced <-   sample_of_interest[which( sample_of_interest$`kefir type`==grains),]
        
        
        if(
          
          TRUE %in% (as.numeric(gsub("T","",sample_of_interest_reduced $Stage[which(duplicated(sample_of_interest_reduced $Sample))][1]))+1==
                     as.numeric(gsub("T","", sample_of_interest_reduced $Stage[which(duplicated(sample_of_interest_reduced $Sample))][2])))
          
          
        ){
          
          
          
          species_of_interest_detected_in_grain <- rbind(     sample_of_interest, species_of_interest_detected_in_grain)
        }
        
        
        
        
      } 
    }
  }
}



persistance <- 
  species_of_interest_detected_in_grain %>% 
  mutate(`kefir type`=gsub("M","Milk ",
                           gsub(  "W","Water ",
                                  gsub( "L", "Liquid",
                                        gsub("G", "Grain",`kefir type`))))) %>% 
  mutate(Stage=gsub("T1","wk1",
                    gsub( "T2","wk5",
                          gsub( "T3","wk9",
                                gsub( "T4","wk13",
                                      gsub( "T5","wk17",
                                            gsub( "T6","wk21",Stage))))))) %>%
  mutate(Stage= factor(Stage, levels=c("wk0","wk1","wk5","wk9","wk13","wk17","wk21"))) %>% 
  
  ggplot(aes(x=Stage, y=clade_name, size=relative_abundance, color=relative_abundance)) + # change to number of occurrences 
  geom_point(alpha=0.5)  +
  scale_size(range=c(2, 20), name="Relative abundance (%)")+
  labs(x="",y="", colour="Species")+
  facet_wrap(~`kefir type`,dir="v",scales="free")+
  theme_bw() +
  theme(#plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5),size = 55), #element_text(color="red", size=14, face="bold",hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(size=17.5),
    axis.text.y = element_text(size=13.5,face="italic"),
    strip.background = element_rect(
      color="black", fill="white"),
    strip.text = element_text(size=20),
    #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
    axis.title.y = element_text( size=20, face="bold",hjust = 0.5, vjust = 1.5),
    legend.key.size = unit(2, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm'), #change legend key width
    legend.text = element_text(size=20,face = "italic"),
    legend.title = element_text(size=20, vjust = 0.5))+
  scale_colour_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
                         breaks= c(-5,0,10,30,70,90,100),
                         values= rescale(as.numeric(c(-5,0,10,30,70,90, 100))),
                         guide="colorbar",
                         name="Relative abundance (%)" ,
                         labels=c(-5,0,10,30,70,90,100))+
  guides(size = guide_legend(label.position = "bottom")) 
#title.position = "right", title.vjust = 0.1))

# guides(color= guide_legend(), size=guide_legend())






library(ggpubr)


jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/Evolution/manuscript figures/Figure4.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)


 ggarrange(
   ggarrange(pa,pb, nrow=1,ncol=2, labels=c("A.","B."),common.legend = FALSE,  font.label = list(size = 30)),
             persistance, nrow=2,ncol=1, heights = c(4, 6),widths=c(10,1),labels=c("","C."),  font.label = list(size = 30))
graphics.off()










































































########################################################################################################################
#Alpha diversity across fermentation categories
########################################################################################################################


#aes(, order = order)
ggplot(my.files_summary_cs, aes(x=fct_reorder2(category_confirmed, -as.numeric(shannon),`kefir type.x`), y=as.numeric(shannon),fill=category_confirmed)) +
  geom_boxplot() +
  #facet_wrap(~ medium)+
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
        legend.position = "none")
#


my.files_summary_cs %>%
  mutate(meanViz = mean(as.numeric(shannon), na.rm = TRUE)) %>%
  ggplot(aes(x=fct_reorder2(category_confirmed, as.numeric(meanViz),as.factor(`kefir type.x`)), y=as.numeric(shannon),fill=category_confirmed)) +
  geom_boxplot() +
  #facet_wrap(~ medium)+
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
        legend.position = "none")





my.files_summary_cs %>%
  mutate(meanViz = mean(as.numeric(shannon), na.rm = TRUE)) %>%
  arrange(desc(meanViz)) %>% 
  group_by(`kefir type.x`) %>% 
  
  
  
  ########################################################################################################################
#Alpha diversity across fermentation categories and timepoints
########################################################################################################################


#View(table(dplyr::select(my.files_summary_cs[-c(which(my.files_summary_cs$`kefir type.x`=="Media control")),], Stage, category_type)))



ggplot(my.files_summary_cs[-c(which(my.files_summary_cs$`kefir type.x`=="Media control")),], aes(x=Stage, y=as.numeric(shannon),fill=Stage)) +
  geom_boxplot() +
  facet_wrap(~ category_type)+
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
        legend.position = "none")