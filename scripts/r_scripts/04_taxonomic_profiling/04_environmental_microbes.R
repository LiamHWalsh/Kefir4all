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

﻿########################################################################################################################
#Research objectives in this script
# See if there is a differnece between sterile and non sterile conditions
# Identify pathogenic microbes
# see how household environment affects there prevalence 
################################################################################################################################################################################################################################################



                                                                                   # See if there is a differnece between sterile and non sterile conditions
#############################################################################################################################################################################################################################################

#Libraries used
########################################################################################################################

pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis)
pacman::p_load(readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)
library(vegan)
library(ggplot2)
library(grid)


########################################################################################################################
#Import metacache data
########################################################################################################################
metacache <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_species_profile_v2.csv")
metacache_strain <-  read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_strain_profile.csv")

metacache$sample_id <- gsub("-","_",metacache$sample_id)
metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]


########################################################################################################################
#Import metadata
########################################################################################################################

global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv")
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv")
Citizen_Scientist_metadata_v8 <- read_csv(CS_METADATA_PRIVATE)

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv")
kefir4all_metadata$merge_column <-  gsub("_host_removed_R..fastq.gz","",kefir4all_metadata$merge_column)
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

my.files_summary$conditions <-  "Non- sterile conditions"

my.files_summary$conditions[which(my.files_summary$data_source %in% c("Walsh et al 2023","Mortensen et al 2023" ))] <- "Sterile conditions"

setwd("Q:/H2020 Master/Citizen Science Project/Plots/Evolution")

#jpeg(filename='Alpha diversity_kefir_type.jpeg', width = 35*700, height=30*700,res=1700,pointsize = 15) #, width=2000, height=1950)
########################################################################################################################
#Plot alpha diversity between samples in sterile conditions (global studies) and kefir4all (non sterile conditions)
########################################################################################################################
pa <- 
  ggplot(my.files_summary %>% 
           filter(`kefir type` %in% c("ML","WL")), aes(x=`kefir type`, y=as.numeric(shannon),fill=conditions)) +
  geom_boxplot() +
  #labs(title= 'Alpha diversity of timepoints') +
  geom_point(position=position_dodge(width=0.75),aes(group=conditions))+
  theme_bw()+
  xlab("Fermentation conditions")+
  ylab("Alpha diversity values (Shannon)")+
  labs(fill = "Timepoint")+
  #guides(colour = guide_legend(override.aes = list(size=25)))+
  #ylim(0,250)+
  theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
        legend.title = element_text( size=25, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
        axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2))+
  scale_x_discrete(labels=c('Milk Liquid', 'Water - Liquid'))


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






################################################################################################################################################################################################################################################

                                                                                                            #Identify environmental microbes

#############################################################################################################################################################################################################################################




########################################################################################################################

########################################################################################################################





########################################################################################################################
#Libraries used
########################################################################################################################

pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis)
pacman::p_load(readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)
library(vegan)
library(ggplot2)
library(grid)


########################################################################################################################
#Import metacache data
########################################################################################################################
metacache <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_species_profile_v2.csv")


metacache_strain <-  read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_strain_profile_v2.csv")

metacache$sample_id <- gsub("-","_",metacache$sample_id)
metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]

########################################################################################################################
#Import environmetal microbes list generated in script 04_compositional_heatmap_metacache_05_02_23
########################################################################################################################
environmental_microbes <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/environmental_microbes.txt")

########################################################################################################################
#Import prevalence data 
########################################################################################################################
metacache_prevalence<-c()
metacache_prevalence[["Milk.kefir"]] <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/prevalence/milk_metacache_prevalence.csv")


metacache_prevalence[["Water.kefir"]] <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/prevalence/water_metacache_prevalence.csv")

########################################################################################################################
#Import and modify survey_responses 
########################################################################################################################
setwd("Q:/H2020 Master/Citizen Science Project/Results/00/")
temp = list.files(pattern=".xlsx", recursive = FALSE)
myfiles = lapply(temp,read_excel)
names(myfiles) <- gsub("file_names_survey_responses_|.xlsx","",temp)



#########################################################################################################################
#Import metadata
########################################################################################################################

global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv")
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv")
Citizen_Scientist_metadata_v8 <- read_csv(CS_METADATA_PRIVATE)

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv")

#kefir4all_metadata$`kefir type`[grep("MLC|WLC",kefir4all_metadata$merge_column)] <- "Liquid control"
kefir4all_metadata$merge_column <-  gsub("_host_removed_R..fastq.gz","",kefir4all_metadata$merge_column)
kefir4all_metadata <- kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$merge_column))),]


total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`))



print(paste("The citizen science project is composed of ",nrow(kefir4all_metadata))) 

print(paste("number of sample according to types in this study"))

print(table(kefir4all_metadata$`kefir type`))


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


mediums_controls_type[[i]] <-mediums_controls_type[[i]] [,which(colnames(mediums_controls_type[[i]] ) %in%  n1_species )]


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
#Import and modify fermented food survey_responses 
########################################################################################################################
setwd("Q:/H2020 Master/Citizen Science Project/Results/00/")
temp = list.files(pattern=".xlsx", recursive = FALSE)
myfiles = lapply(temp,read_excel)
names(myfiles) <- gsub("file_names_survey_responses_|.xlsx","",temp)



#########################################################################################################################
#Import getting started survey metadata
########################################################################################################################

#Teagasc_getting_started <- read_excel(file.path(DATA_DIR, "Teagasc - getting started.xlsx")
library(surveymonkey)
# Import the options command into your r_profile using the below command 
#usethis::edit_r_profile()
#options(sm_oauth_token = "vOcHPmfepVilX6cUC5JVr4axcPfCv4UxPKhWEd-rbegkTh58EF2I.SvzlZYBpNnQliGeObex30.GNzFQP8lUb.UVzWMjHVeli4G-pgueVnnuhu66IEzB5ovNjlYmP1YU")
getOption("sm_oauth_token")
Teagasc_getting_started <-  313643610 %>%
  fetch_survey_obj %>%
  parse_survey


Teagasc_getting_started$`Please enter your kefir ID` <- gsub(".*ID0|IDØ|Id0|1D0|Id","",Teagasc_getting_started$`Please enter your kefir ID`)

Teagasc_getting_started$`Please enter your kefir ID` <- gsub("ID|id|1D|ID |IDD|\\.|","",Teagasc_getting_started$`Please enter your kefir ID`)

Teagasc_getting_started$`Please enter your kefir ID` <- gsub(" ","",Teagasc_getting_started$`Please enter your kefir ID`)

#Citizen_Scientist_metadata



Teagasc_getting_started$`Please enter your kefir ID` <- paste("ID",Teagasc_getting_started$`Please enter your kefir ID`,sep = "")
Teagasc_getting_started$`Please enter your kefir ID` <- gsub("ID|ID0|ID00","",Teagasc_getting_started$`Please enter your kefir ID`)

Teagasc_getting_started$`Please enter your kefir ID`[which(as.numeric(Teagasc_getting_started$`Please enter your kefir ID`) %in% c(1:9))] <- paste( "ID00",Teagasc_getting_started$`Please enter your kefir ID`[which(as.numeric(Teagasc_getting_started$`Please enter your kefir ID`) %in% c(1:9))],sep="")

Teagasc_getting_started$`Please enter your kefir ID`[which(as.numeric(Teagasc_getting_started$`Please enter your kefir ID`) %in% c(10:99))] <- paste( "ID0",Teagasc_getting_started$`Please enter your kefir ID`[which(as.numeric(Teagasc_getting_started$`Please enter your kefir ID`) %in% c(10:99))],sep="")

Teagasc_getting_started$`Please enter your kefir ID`[which(as.numeric(Teagasc_getting_started$`Please enter your kefir ID`) >=100)] <- paste( "ID",Teagasc_getting_started$`Please enter your kefir ID`[which(as.numeric(Teagasc_getting_started$`Please enter your kefir ID`) >=100)],sep="")




Teagasc_getting_started <- Teagasc_getting_started %>% mutate(pets = coalesce(
 
    `Do you have any pets in your household? - no pets`,
    `Do you have any pets in your household? - yes, cat(s)`,
    `Do you have any pets in your household? - yes, dog(s)`,
    `Do you have any pets in your household? - yes, rabbit(s)/hamster(s)/...`,
    `Do you have any pets in your household? - yes, other (please specify)`))

Teagasc_getting_started <- Teagasc_getting_started %>% mutate(pets = coalesce(
  
  `Do you have any pets in your household? - no pets`,
  `Do you have any pets in your household? - yes, cat(s)`,
  `Do you have any pets in your household? - yes, dog(s)`,
  `Do you have any pets in your household? - yes, rabbit(s)/hamster(s)/...`,
  `Do you have any pets in your household? - yes, other (please specify)`))


Teagasc_getting_started <- Teagasc_getting_started %>% mutate(other_fermented_food=coalesce(
`Do you produce other fermented foods or beverages at home? - no`,
`Do you produce other fermented foods or beverages at home? - yes, a different kefir`,
`Do you produce other fermented foods or beverages at home? - yes, other (please specify)`,
`Do you produce other fermented foods or beverages at home? - yes, sauerkraut`,
`Do you produce other fermented foods or beverages at home? - yes, a different kefir`,
`Do you produce other fermented foods or beverages at home? - yes, kombucha`,
`Do you produce other fermented foods or beverages at home? - yes, kimchi`
))


#########################################################################################################################
#merge metadata 
########################################################################################################################



total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`))





#########################################################################################################################
#merge metadata 
########################################################################################################################

kefir4all_metadata <- merge(kefir4all_metadata, Teagasc_getting_started,by.x="Sample",by.y="Please enter your kefir ID",all.x=TRUE)
kefir4all_metadata$pets_type <- "NA" 
kefir4all_metadata$pets[
kefir4all_metadata$Stage=="T0"] <- "no pets"


kefir4all_metadata$pets_type[-c(which(
kefir4all_metadata$pets %in% c("no pets",NA)))] <- "yes"

kefir4all_metadata$pets_type[which(kefir4all_metadata$pets=="no pets")] <- "no"


kefir4all_metadata$other_fermented_food_type <- "NA"
kefir4all_metadata$other_fermented_food_type[which(kefir4all_metadata$other_fermented_food=="no")] <- "no"
# 
 kefir4all_metadata$other_fermented_food_type[
   kefir4all_metadata$Stage=="T0"] <- "no"


kefir4all_metadata$other_fermented_food_type[-c(which(kefir4all_metadata$other_fermented_food %in% c("no",NA)))] <- "yes"



kefir4all_metadata <- 
  merge(
    kefir4all_metadata, 
    rbind(dplyr::select(myfiles[["mk"]], merge_column,Sample, observations, category_confirmed),
          dplyr::select(myfiles[["wk"]], merge_column,Sample, observations, category_confirmed)),
    by="merge_column",
    all.x=TRUE)

kefir4all_metadata <- 
merge(kefir4all_metadata, 
      dplyr::select(myfiles[["wk"]],merge_column,water_type),
      by="merge_column",
      all.x=TRUE
      )





########################################################################################################################
# get grain samples from metacache 
########################################################################################################################

grain_microbes <- c()

rare_grain_microbes[["Milk.kefir"]] <- grain_microbes[["Milk.kefir"]] <- 
  metacache[
    which(metacache$sample_id %in%
            kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                    kefir4all_metadata$`kefir type` %in% c("MG"))]
    ),] %>% pivot_longer(!sample_id,names_to="species",values_to="relative_abundance") %>% 
  filter(relative_abundance<0.1)

rare_grain_microbes[["Milk.kefir"]] <- rare_grain_microbes[["Milk.kefir"]][which(rare_grain_microbes[["Milk.kefir"]]$species %in% 
                                            metacache_prevalence[["Milk.kefir"]]$specie)                                          
                                            ,]


grain_microbes[["Milk.kefir"]] <- 
  metacache[
    which(metacache$sample_id %in%
            kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                    kefir4all_metadata$`kefir type` %in% c("MG"))]
    ),] %>% pivot_longer(!sample_id,names_to="species",values_to="relative_abundance") %>% 
  filter(relative_abundance>0.1)





grain_microbes[["Water.kefir"]] <- 
  
  metacache[
    which(metacache$sample_id %in%
            kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                    kefir4all_metadata$`kefir type` %in% c("WG"))]
    ),]



rare_grain_microbes[["Water.kefir"]] <- 
  
  grain_microbes[["Water.kefir"]][which(grain_microbes[["Water.kefir"]]$sample_id=="TG_FL_T2_WG_0AN_S328"),]%>% pivot_longer(!sample_id,names_to="species",values_to="relative_abundance") %>% 
  filter(relative_abundance<0.1)
  
grain_microbes[["Water.kefir"]] <- 
  grain_microbes[["Water.kefir"]][which(grain_microbes[["Water.kefir"]]$sample_id=="TG_FL_T2_WG_0AN_S328"),]%>% pivot_longer(!sample_id,names_to="species",values_to="relative_abundance") %>% 
 filter(relative_abundance>0.1)


########################################################################################################################
#what medium micorbes were not found in the grain
########################################################################################################################


medium_prevalence[["Milk.kefir"]]$specie[which(medium_prevalence[["Milk.kefir"]]$specie%in%  metacache_prevalence[["Milk.kefir"]]$specie)]

medium_prevalence[["Water.kefir"]]$specie[which(medium_prevalence[["Water.kefir"]]$specie%in%  metacache_prevalence[["Water.kefir"]]$specie)]

                  


########################################################################################################################
#what micorbes can be identified in the grains that weren't their at baseline
########################################################################################################################

metacache_grain <- c()


base <- c()

for (type in names(grain_microbes)){
  
  if(type=="Milk.kefir"){
    
  
  base <- "MG"}else{
    base <- "WG"
    
  }
  metacache_grain <-   metacache[which(metacache$sample_id %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type`==base)]),] %>% 
                                   pivot_longer(!sample_id,names_to = "species",values_to = "relative_abundance") %>% 
                                  filter(relative_abundance>.01)
  grain_microbes[[type]]
}


########################################################################################################################
#What micorbs bloom in the kefir fermentations. 
########################################################################################################################





########################################################################################################################
#Modufy dataset into metacache citizen scientist and baseline used
########################################################################################################################

kefir4all_metadata <- kefir4all_metadata[-c(which(kefir4all_metadata$`kefir type`=="Medium control" )),]
kefir4all_metadata <- kefir4all_metadata[-c(which(kefir4all_metadata$`kefir type`=="Extraction control")),]

metacache [is.na(metacache )] <- 0

metacache <- metacache[which(metacache$sample_id %in% 
                               kefir4all_metadata$merge_column), ]

#n1_species <-which(maxab_species < 0.1)


metacache <- 
  metacache %>% column_to_rownames("sample_id")


my.files_summary_shannon <- diversity(metacache  )
my.files_richness <- specnumber(metacache)
my.files_evenness<- my.files_summary_shannon/log(my.files_richness)
my.files_beta<- vegdist(metacache, method = "bray")

my.files_summary<- cbind(shannon = my.files_summary_shannon, richness = my.files_richness, pielou = my.files_evenness,site=metacache$sample_id)


my.files_summary<- as.data.frame(my.files_summary)
my.files_summary <- 
merge(my.files_summary,kefir4all_metadata,by.x=0,by.y="merge_column",all.X=TRUE)

# Look at household pet infleunce in total data 


#my.files_summary <- my.files_summary[-c(which(my.files_summary$`kefir type`== "Liquid control")),]


# 
# 
# ########################################################################################################################
# #Maaslin for pets
# ########################################################################################################################
# 
# library(tidyverse)
# library(data.table)
# library(vegan)
# library(ape)
# library(cowplot)
# library(Maaslin2)
# library(dplyr)
# library(ggtext)
# 
# 
#   
#   
# species_profile <- c()
# sub_metadata <- c()
# type <- "mk"
# 
# 
# for (type in c("mk","wk")){
# 
#   if(type=="mk"){
#     sub_metadata <- kefir4all_metadata[which(kefir4all_metadata$`kefir type` %in% c("WG","WL")),]
#   }else{
#   
#   sub_metadata <- kefir4all_metadata[which(kefir4all_metadata$`kefir type` %in% c("WG","WL")),]
#   }
#   sub_metadata <- sub_metadata[-c(which(duplicated(sub_metadata$merge_column))),]
# 
#   #sub_metadata <- sub_metadata[-c(which(duplicated(sub_metadata$merge_column))),]
#   sub_metadata <- sub_metadata %>% remove_rownames()%>% column_to_rownames("merge_column")
#   
#   sub_metadata <- sub_metadata[-c(  
#     which(sub_metadata$pets_type=="NA")),]
#   
#   
#   species_profile[[type]] <- metacache[which(rownames(metacache) %in%  rownames(sub_metadata)),]#& 
#   
  
  #sub_metadata <- dplyr::select(myfiles[[type]], merge_column,Sample, Stage, Room_temperature, pH, grain_wash, Liquid_volume, grains_weight, observations, milk_type,total_fermentation_response_frequency)
  
  
  #sub_metadata <- sub_metadata[-c(which(duplicated(sub_metadata$merge_column))),]
  
  
  
  
 # unlink("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/Maaslin2",recursive = TRUE)
  # dir.create("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/Maaslin2/pets")
  # dir.create(paste("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/Maaslin2/pets/",type,sep=""))
  # out_dir <- paste("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/Maaslin2/pets/",type,sep="")
  # 
  # 
  # 
  # Maaslin2(
  #   species_profile[[type]],
  #   sub_metadata,
  #   out_dir,
  #   fixed_effects = c("pets_type"),
  #   reference = c("pets_type,no"),
  #   random_effects = c('Sample'),
  #   min_prevalence = 0.1
 # )
  # 
#}


# ########################################################################################################################
# #Maaslin for other feremented foods
# ########################################################################################################################
  
  # 
  # species_profile <- c()
  # sub_metadata <- c()
  # type <- "mk"
  # 
  # 
  # for (type in c("mk","wk")){
  #   
  #   if(type=="mk"){
  #     sub_metadata <- kefir4all_metadata[which(kefir4all_metadata$`kefir type` %in% c("WG","WL")),]
  #   }else{
  #     
  #     sub_metadata <- kefir4all_metadata[which(kefir4all_metadata$`kefir type` %in% c("WG","WL")),]
  #   }
  #   sub_metadata <- sub_metadata[-c(which(duplicated(sub_metadata$merge_column))),]
  #   
  #   #sub_metadata <- sub_metadata[-c(which(duplicated(sub_metadata$merge_column))),]
  #   sub_metadata <- sub_metadata %>% remove_rownames()%>% column_to_rownames("merge_column")
  # 
  #   
  #   sub_metadata <- sub_metadata[-c(  
  #     which(sub_metadata$other_fermented_food_type=="NA")),]
  #   
  #   
  #   species_profile[[type]] <- metacache[which(rownames(metacache) %in%  rownames(sub_metadata)),]#& 
  #   
    
    #sub_metadata <- dplyr::select(myfiles[[type]], merge_column,Sample, Stage, Room_temperature, pH, grain_wash, Liquid_volume, grains_weight, observations, milk_type,total_fermentation_response_frequency)
    
    
    #sub_metadata <- sub_metadata[-c(which(duplicated(sub_metadata$merge_column))),]
    
    
  
  # dir.create("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/Maaslin2/other_fermented_foods")
  # dir.create(paste("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/Maaslin2/other_fermented_foods/",type,sep=""))
  # out_dir <- paste("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/Maaslin2/other_fermented_foods/",type,sep="")
  # 
  
  # 
  # Maaslin2(
  #   species_profile[[type]],
  #   sub_metadata,
  #   out_dir,
  #   fixed_effects = c("other_fermented_food_type"),
  #   reference = c("other_fermented_food_type,no"),
  #   random_effects = c('Sample'),
  #   min_prevalence = 0.1
  # )
  # 
  
#}
  
  






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
  mutate(Stage= factor(Stage, levels=c("wk0","wk1","wk5","wk9","wk13","wk17","wk21"))
         
  ) %>% 
  
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



########################################################################################################################
#extract the environmental microbes 
########################################################################################################################




species_profile <- c()
sub_metadata <- c()



for (type in c("Milk.kefir","Water.kefir")){

  if(type=="Milk.kefir"){
    sub_metadata[[type]] <- kefir4all_metadata[which(kefir4all_metadata$`kefir type` %in% c("MG","ML")),]
  }else{

    sub_metadata[[type]] <- kefir4all_metadata[which(kefir4all_metadata$`kefir type` %in% c("WG","WL")),]
  }
  sub_metadata[[type]] <- sub_metadata[[type]][-c(which(duplicated(sub_metadata[[type]]$merge_column))),]

  #sub_metadata <- sub_metadata[-c(which(duplicated(sub_metadata$merge_column))),]
  sub_metadata[[type]] <- sub_metadata[[type]] %>% remove_rownames()%>% column_to_rownames("merge_column")
  
  
    species_profile[[type]] <- metacache[which(rownames(metacache) %in%  rownames(sub_metadata[[type]])),]
    
    
    
    species_profile[[type]] <- species_profile[[type]][, which(colnames( species_profile[[type]]) %in% environmental_microbes$species[which(environmental_microbes$type %in% c("Occurs in both", type))]) ] %>% as.data.frame() %>% 
      rownames_to_column("sample_id") %>% pivot_longer(!sample_id,names_to = "species",values_to = "relative_abundance")
    
    species_profile[[type]] <- merge(    species_profile[[type]], sub_metadata[[type]],by.x="sample_id",by.y=0,all.x=TRUE)
      
  
}




t <- 
species_profile[[type]][which(    species_profile[[type]]$relative_abundance>1),]

#dir.create("Q:/H2020 Master/Citizen Science Project/plots/Evolution/04_environmental_microbes_boxplots")
i<- "pets_type"
#other_fermented_food_type
#observations
p <- c()
for (type in c("Milk.kefir","Water.kefir")){
  
for (i in names(dplyr::select(species_profile[[type]],other_fermented_food_type ,pets_type, category_confirmed ))
){
  
  # 


  # p <- 
  # species_profile[[type]][c(which(getElement(species_profile[[type]],i)!="NA"
  #                                                                )),]%>% 
  #   filter(relative_abundance>.1) %>% 
  # 
  #      ggplot(aes_string(x="`kefir type`", y="relative_abundance",fill=i)) +
  #        geom_boxplot()+
  #        #labs(title= 'Alpha diversity of timepoints') +
  #        geom_point(position=position_dodge(width=0.75),aes_string(group=i))+
  #        theme_bw()+
  #        xlab("Time points")+
  #        ylab("Alpha diversity values (Shannon)")+
  #        labs(fill = "Timepoint")+
  #        # scale_x_discrete(labels=c('Milk - Grains', 'Milk Liquid', 'Water - Grains', 'Water - Liquid'))+
  #        facet_wrap(~species,scales = "free_y")+
  #        #guides(colour = guide_legend(override.aes = list(size=25)))+
  #        #ylim(0,250)+
  #        theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
  #              legend.title = element_text( size=25, face="bold"),
  #              axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
  #              axis.text.y = element_text(hjust = 1, size = 10),
  #              axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
  #              axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
  #              legend.position = "right")
  # jpeg(filename=paste('Q:/H2020 Master/Citizen Science Project/plots/Evolution/04_environmental_microbes_boxplots/',i,'_',type,'_changes_environmental_microbes_boxplot.jpeg',sep=""), width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
  # 
  # plot(p)
  # graphics.off()
  
    
  # [c(which(getElement(species_profile[[type]],i)!="NA"
  #                               )),] %>% 
  
  
}

}

# species_profile[[type]][c(which(getElement(species_profile[[type]],i)!="NA")),]

species <-  "Acinetobacter lwoffii"   
species <- "Leclercia adecarboxylata"
species <- "Lentilactobacillus hilgardii"
species <- 

type <- "Milk.kefir"

species <- 
"Raoultella terrigena"

t <- c()

for (type in c("Milk.kefir","Water.kefir")){
  environmental_microbes[paste(type,"_frequency",sep="")] <- NA
  environmental_microbes[paste(type,"_pets",sep="")] <- NA
  environmental_microbes[paste(type,"_other_fermented_foods",sep="")] <- NA

for (species in levels(as.factor(species_profile[[type]]$species))){
  
  t <- species_profile[[type]][which( species_profile[[type]]$species==species),] %>% 
    filter(relative_abundance>.1)

  
  environmental_microbes[which(environmental_microbes$species==species),which(colnames(environmental_microbes)==paste(type,"_frequency",sep=""))] <- nrow(t)
  environmental_microbes[which(environmental_microbes$species==species),which(colnames(environmental_microbes)==paste(type,"_pets",sep=""))] <- paste(levels(as.factor(t$pets)), sep=" ", collapse = " ")
  environmental_microbes[which(environmental_microbes$species==species),which(colnames(environmental_microbes)==paste(type,"_other_fermented_foods",sep=""))] <- paste(levels(as.factor(t$other_fermented_food)), sep=" ", collapse = " ")
  
  if(type=="water.kefir"){
    environmental_microbes[paste(type,"_water_type",sep="")] <- NA
    
    environmental_microbes[which(environmental_microbes$species==species),which(colnames(environmental_microbes)==paste(type,"_water_type",sep=""))] <- paste(levels(as.factor(t$water_type)), sep=" ",collaspe=",")
  }
}


}

environmental_microbes$species[which(environmental_microbes$Milk.kefir_frequency>=1)][

environmental_microbes$species[which(environmental_microbes$Milk.kefir_frequency>=1)] %in% colnames(mediums_controls_type[["Milk.kefir"]])]


i <- "water_type" 





p <- 
#[[type]][c(which(getElement(species_profile[[type]],i)!="NA"
#)),]%>% 
  species_profile[["Water.kefir"]][which( species_profile[[type]]$species=="Limosilactobacillus fermentum"),] %>% 
  #filter(relative_abundance>0.1) %>% 
  ggplot(aes_string(x="`kefir type`", y="relative_abundance",fill=i)) +
  geom_boxplot()+
  #labs(title= 'Alpha diversity of timepoints') +
  geom_point(position=position_dodge(width=0.75),aes_string(group=i))+
  theme_bw()+
  xlab("Time points")+
  ylab("Alpha diversity values (Shannon)")+
  labs(fill = "Timepoint")+
  # scale_x_discrete(labels=c('Milk - Grains', 'Milk Liquid', 'Water - Grains', 'Water - Liquid'))+
  facet_wrap(~species,scales = "free_y")+
  #guides(colour = guide_legend(override.aes = list(size=25)))+
  #ylim(0,250)+
  theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
        legend.title = element_text( size=25, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
        axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
        legend.position = "right")



#jpeg(filename=paste('Q:/H2020 Master/Citizen Science Project/plots/Evolution/04_environmental_microbes_boxplots/',i,'_',type,'changes_environmental_microbes_boxplot.jpeg',sep=""), width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)

#plot(p)
#graphics.off()






ggplot(env_microbes_long, aes(x=Stage, y=as.numeric(relative_abundance),fill=other_fermented_food_type)) +
  geom_boxplot()+
  #labs(title= 'Alpha diversity of timepoints') +
  geom_point(position=position_dodge(width=0.75),aes(group=other_fermented_food_type))+
  theme_bw()+
  xlab("Time points")+
  ylab("Alpha diversity values (Shannon)")+
  labs(fill = "Timepoint")+
  #  scale_x_discrete(labels=c('Milk - Grains', 'Milk Liquid', 'Water - Grains', 'Water - Liquid'))+
  facet_wrap(~species+`kefir type`,scales = "free_y")+
  #guides(colour = guide_legend(override.aes = list(size=25)))+
  #ylim(0,250)+
  theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
        legend.title = element_text( size=25, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
        axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
        legend.position = "right")

























# Now look at household influence in samples in wh

graphics.off()
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


my.files_summary_cs <- merge(my.files_summary_cs,kefir4all_total_metadata,by.x="site",by.y="Row.names",all.x=TRUE)



print(paste("Number of samples in metacache is",length(rownames(metacache))))

# samples  not accounted for by metacache- base_line typically did not have a compostional profile 






total_compositional_data <- metacache%>%
  rownames_to_column("sample") %>% pivot_longer(!sample  , names_to  ="clade_name" , values_to ="relative_abundance")

total_compositional_data <-  
  merge(total_compositional_data,kefir4all_metadata,by.x="sample",by.y="merge_column")



total_compositional_data$type <- NA


total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                                                 kefir4all_metadata$Stage!="T0")])] <- "Water kefir"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                                                 kefir4all_metadata$Stage!="T0")])] <- "Milk kefir"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                                                 kefir4all_metadata$Stage=="T0")])] <- "Water\n kefir\n T0"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                                                 kefir4all_metadata$Stage=="T0")])] <- "Milk\n kefir\n T0"


# remove extraction controls
total_compositional_data <- total_compositional_data[-c(which(is.na(
  total_compositional_data$type))),]

#total_compositional_data$type <- gsub("Milk kefir_baseline","Milk kefir - T0",total_compositional_data$type)
#total_compositional_data$type <- gsub("Water kefir_baseline","Water kefir - T0",total_compositional_data$type)


#print(paste(table(dplyr::select(total_compositional_data,`kefir type`))))

###########################################################################################################################################################

# Heatmap of all prevalent microbes

###########################################################################################################################################################


# jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/Evolution/Figure 1_v3.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)



total_compositional_data$kefirdataset <- NA
total_compositional_data$kefirdataset[which(total_compositional_data$`kefir type` %in% c("WL","WG"))] <- "Water kefir"

total_compositional_data$kefirdataset[which(total_compositional_data$`kefir type` %in% c("ML","MG"))] <- "Milk kefir"
#which(is.na(total_compositional_data$kefirdataset))



total_compositional_data %>% 
  ggplot(aes(x=type,y=relative_abundance,fill=type))+
  geom_boxplot()+
  facet_wrap(~clade_name, scales = "free")



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




# t <-
# subset( compositional_differences, subset=
#           compositional_differences$Milk.kefir.T0 >=0 &
#           compositional_differences$Water.kefir.T0 ==0 &
#           compositional_differences$Milk.kefir>=0&
#           compositional_differences$Water.kefir==0
# )$species

###########################################################################################################################################################


# Output prevalence into

###########################################################################################################################################################


pacman::p_load(officer,magrittr,officedown)

prevalence_doc <- read_docx() %>% 
  body_add_fpar(fpar(
    ftext(  " Look at the prevalence of both the milk and water kefir cs dataset, not including the global studies ",
            fp_text(bold = FALSE, font.size = 11,underlined=FALSE,font.family = "Calibri"))))
#dir.create("Q:/H2020 Master/Citizen Science Project/Manuscripts/officer_outputs/")
for (i in names(t)){
  
  
  prevalence_doc <- 
    
    prevalence_doc %>%
    body_add_fpar(  fpar( ftext( paste("We identified the prevalent microbial species in", i, " to be"),
                                 fp_text(bold = FALSE, font.size = 11,underlined=FALSE,font.family = "Calibri")) ,
                          
                          ftext(paste(t[[i]][1:length(t[[i]])-1],collapse=", "),
                                fp_text(bold = FALSE,italic = TRUE, font.size = 11,underlined=FALSE,font.family = "Calibri")),
                          ftext(  " and ",
                                  fp_text(bold = FALSE, font.size = 11,underlined=FALSE,font.family = "Calibri")) , 
                          ftext(paste(t[[i]][length(t[[i]])]),
                                fp_text(bold = FALSE,italic = TRUE, font.size = 11,underlined=FALSE,font.family = "Calibri"))))
}


#prevalence_doc %>%
#print(target = "Q:/H2020 Master/Citizen Science Project/Manuscripts/officer_outputs//cs_prevalence_output.docx")

###########################################################################################################################################################


# Of these prevlent microbes, what is the range of RA they occur at and whats their number of samples tht they dominate in 

###########################################################################################################################################################


total_compositional_data %>% 
  ggplot(aes(x=type,y=relative_abundance,fill=type))+
  geom_boxplot()+
  facet_wrap(~clade_name, scales="free")


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

for (sample_name in levels(as.factor(total_compositional_data$sample))){
  dominanting_species <- rbind(dominanting_species,
                               
                               total_compositional_data %>% 
                                 filter(sample== sample_name) %>% 
                                 filter(relative_abundance==max(relative_abundance))
                               
  )
}

# breakdown of dominate species per kefir type

dominate_breakdown <- as.data.frame(
  table(dplyr::select(  dominanting_species, clade_name,`kefir type`))) 

dominate_breakdown$Freq <- as.numeric(dominate_breakdown$Freq)

dominate_breakdown <- 
  dominate_breakdown%>% pivot_wider(names_from="kefir.type",values_from = "Freq")



# range of values per dominate type

dominate_breakdown$min <- NA
dominate_breakdown$max <- NA


for (sample_name in levels(as.factor(dominate_breakdown$clade_name))){
  
  
  
  
  dominate_breakdown$min[which(dominate_breakdown$clade_name==sample_name)] <-  getElement(   dominanting_species %>% 
                                                                                                filter(clade_name== sample_name) %>% 
                                                                                                filter(relative_abundance==min(relative_abundance)),"relative_abundance")
  
  
  dominate_breakdown$max [which(dominate_breakdown$clade_name==sample_name)]<-  getElement(   dominanting_species %>% 
                                                                                                filter(clade_name== sample_name) %>% 
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


dominate_breakdown$MG_percent <- (dominate_breakdown$MG / samples_grain_mk)*100

dominate_breakdown$ML_percent <- (dominate_breakdown$ML / samples_ml)*100
dominate_breakdown$WG_percent <- (dominate_breakdown$WG / samples_wg)*100
dominate_breakdown$WL_percent <- (dominate_breakdown$WL / samples_wl)*100

length(rownames(metacache)[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL"))])])







# = ggplot(pcoa, aes(x=PC1, y=PC2)) + geom_point(size=3) + theme_bw()#


#pcoa$kefir_dataset <- "Milk.kefir"
#pcoa$kefir_dataset[which(pcoa$`kefir type.x` %in% c("WL","WG"))] <-"Water kefir"

#pcoa <- pcoa[-c(which(pcoa$`kefir type.x`=="Liquid control")),]




table(dplyr::select(dominanting_species, clade_name, `kefir type`))
length(levels(as.factor(total_compositional_data$sample)))

#total_prevalence <- rbind(species_prevalence[["milk"]],species_prevalence[["water"]])



###########################################################################################################################################################


# Sources of environmental microbes occurring at >.01 relative abundance

###########################################################################################################################################################

maxab_species <- apply(metacache,2, max, na.rm=TRUE)

# all species at >.1 % RA not jst prevalent micorbes
species_data  <- metacache [,-c(which(maxab_species < .1))] %>% rownames_to_column("sample") %>% pivot_longer(!sample,names_to = "clade_name",values_to = "relative_abundance")


species_data <-  
  merge(species_data,kefir4all_metadata,by.x="sample",by.y="merge_column")


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
species_data <- species_data[-c(which(is.na(
  species_data$type))),]


compositional_differences_.1 <- data.frame(species=as.character( levels(as.factor(species_data$clade_name))),
                                           `Milk kefir`=as.character(NA),
                                           `Milk kefir T0`=as.character(NA),
                                           `Water kefir`=as.character(NA),
                                           `Water kefir T0`=as.character(NA)
)

species <- "Leuconostoc mesenteroides"
type_2 <- "Milk.kefir"
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
  compositional_differences[which(compositional_differences$Water.kefir.T0==0 & compositional_differences$Water.kefir>0) ,] %>% 
  dplyr::select(species)%>% mutate(type="Water.kefir")


t1 <- t1[-c(which(t1$species %in% species_of_interest$species)),]


species_of_interest <-  rbind(species_of_interest,t,t1)





print(paste("okay, we can see", nrow(species_of_interest), "microbes that may be environmentally derived. That's not the most interesting unless they persist. Lets look at that now"))

environmental_microbes  <- species_of_interest

write.csv(environmental_microbes, "Q:/H2020 Master/Citizen Science Project/Results/00/environmental_microbes.txt",quote=FALSE,row.names = FALSE)
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

# 
# jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/Evolution/Figure 1.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
# 
# 
# ggarrange(
#   ggarrange(p_heatmap[[1]],p_heatmap[[2]], nrow=1,ncol=2, labels=c("A.","B."),common.legend = TRUE,  font.label = list(size = 30)),
#             persistance, nrow=2,ncol=1, heights = c(6, 4),widths=c(10,1),labels=c("","C."),  font.label = list(size = 30))
# graphics.off()
# 
# 
# 
# total_compositional_data %>% 
#   ggplot(aes(x=type,y=relative_abundance,fill=type))+
#   geom_boxplot()+
#   facet_wrap(~clade_name)
# 
# graphics.off()




# aspect.ratio = 2/1)+

#scale_x_discrete(limits=species_data_long$Stage))


species_data %>% 
  ggplot(aes(x=type,y=relative_abundance,fill=type))+
  geom_boxplot()+
  facet_wrap(~clade_name)



##############
