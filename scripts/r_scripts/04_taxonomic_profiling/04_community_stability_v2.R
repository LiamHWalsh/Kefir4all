# Anchor to the project root before here() loads. RStudio Server resets the session cwd
# to $HOME after sourcing .Rprofile, and here() computes its root AT PACKAGE LOAD
# TIME, so force the project root before any here load. This block is inert on non-portable
for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {



  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}


# ==============================================================================
# PATH CONFIGURATION â€” portable relative paths via the here package
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
OUT_DIR <- here::here("output_reviewed", "04_taxonomic_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "Extended_Data_Figure_1")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
PRIMARY_SOURCE <- file.path(DATA_DIR, "04_metacache_total_species_profile_v2.csv")

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
#Libraries used 
########################################################################################################################


# .libPaths() removed â€” use default R library or renv
library(readxl)
library(tidyr)
library(readr)
library(devtools)
library(taxize)
library(ape)
library(ggtree)
library(ggdendro)
# column_to_rownames() below needs tibble; load it explicitly so this runs
# headless (a container has no interactive session to have already loaded it).
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
library(tibble)


########################################################################################################################
#Libraries used
########################################################################################################################

# .libPaths() removed â€” use default R library or renv
library(ggplot2)
library(plyr)
library(dplyr)
library(gplots)
library(vegan)
library(scales)
library(vegan)
library(ggplot2)
library(grid)
if (!requireNamespace("ggstatsplot", quietly=TRUE)) install.packages("ggstatsplot")
library(ggstatsplot)


########################################################################################################################
#Import metacache data
########################################################################################################################
metacache <- read_csv(file.path(DATA_DIR, "04_metacache_total_species_profile_v2.csv"))
metacache_strain <-  read_csv(file.path(DATA_DIR, "04_metacache_total_strain_profile_v2.csv"))

metacache$sample_id <- gsub("-","_",metacache$sample_id)
metacache$sample_id <- gsub(".*/","",metacache$sample_id)

#metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]

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
kefir4all_metadata <- 
  merge(
    kefir4all_metadata, 
    rbind(dplyr::select(myfiles[["mk"]], merge_column,Sample, observations, category_confirmed),
          dplyr::select(myfiles[["wk"]], merge_column,Sample, observations, category_confirmed)),
    by="merge_column",
    all.x=TRUE)



print(paste("The citizen science project is composed of ",nrow(kefir4all_metadata))) 

print(paste("number of sample according to types in this study"))

print(table(kefir4all_metadata$`kefir type`))

########################################################################################################################
#Get data from the global_MK study
########################################################################################################################

#metacache <- metacache[which(gsub("_S.*","", metacache$sample_id) %in% global_mk_metadata$Sample),]

########################################################################################################################
#Modufy dataset into metacache citizen scientist and baseline used
########################################################################################################################
my.files_beta<-c()

metacache [is.na(metacache )] <- 0
kefir4all_metadata <- 
  kefir4all_metadata[-c(which(
    kefir4all_metadata$`kefir type`=="Medium control")),]


metacache_total <- metacache

metacache <- metacache[which(metacache$sample_id %in% 
                               kefir4all_metadata$merge_column[
                                 which(kefir4all_metadata$`kefir type` %in% c("WG","MG"))
                               ]), ]



species_data <- c()

species_data[["Milk.kefir"]] <- metacache[which(metacache$sample_id %in% 
                               kefir4all_metadata$merge_column[
                                 which(kefir4all_metadata$`kefir type` %in% c("MG"))
                               ]), ] |> column_to_rownames("sample_id")

species_data[["Water.kefir"]] <- metacache[which(metacache$sample_id %in% 
                                                  kefir4all_metadata$merge_column[
                                                    which(kefir4all_metadata$`kefir type` %in% c("WG"))
                                                  ]), ] |> column_to_rownames("sample_id")



kefir4all_metadata$Stage <- gsub("T1","wk01",kefir4all_metadata$Stage )
kefir4all_metadata$Stage <- gsub("T2","wk05",kefir4all_metadata$Stage )
kefir4all_metadata$Stage <- gsub("T3","wk09",kefir4all_metadata$Stage )
kefir4all_metadata$Stage <- gsub("T4","wk13",kefir4all_metadata$Stage )
kefir4all_metadata$Stage <- gsub("T5","wk17",kefir4all_metadata$Stage )
kefir4all_metadata$Stage <- gsub("T6","wk21",kefir4all_metadata$Stage )

type <- "Water.kefir"
grain_name <- c()
bray_distance <- c()
for (type in names(species_data)){
  
  if(type=="Milk.kefir"){
    grain_name="MG"
  }else{
    grain_name="WG"
  }

# Get betadiversity distnace between To grain and grain of the other of the participant
my.files_beta[[type]]<- as.matrix(vegdist(species_data[[type]], method = "bray"))


my.files_beta[[type]] <- my.files_beta[[type]][which(rownames(my.files_beta[[type]])==kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                                                                                                                      kefir4all_metadata$`kefir type` %in% grain_name)][1]),]


# Get betadiversity between samples of the sample participant



# 
# 

my.files_beta[[type]] <- 
as.data.frame(my.files_beta[[type]] ) |> rownames_to_column("merge_column")
colnames(my.files_beta[[type]])[2] <- "distance"




my.files_beta[[type]] <- 
  merge(kefir4all_metadata,
        my.files_beta[[type]] , by="merge_column",all.y=TRUE)


my.files_beta[[type]]$type <- gsub("\\."," ",type)


bray_distance[[type]] <- 
  
  
  
  ggbetweenstats(data = my.files_beta[[type]] |> 
                   filter(Stage!="T0"),
                 x=Stage, 
                 y=distance,
                 title=type,
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
#


# my.files_beta[[type]] |> 
#   #facet_wrap(~ medium)+
#   #labs(title= 'Alpha diversity of timepoints') +
#   #guides(colour = guide_legend(override.aes = list(size=25)))+
#   #ylim(0,250)+
# #

for (timepoint in levels(as.factor(my.files_beta[[type]]$Stage))){
  
  
  print(paste(timepoint,"-",type, "min distace is ", min(my.files_beta[[type]]$distance[which(my.files_beta[[type]]$Stage==timepoint)])))
  
  
  print(paste(timepoint,"-",type, "mean distace is ", mean(my.files_beta[[type]]$distance[which(my.files_beta[[type]]$Stage==timepoint)])))
 
  print(paste(timepoint,"-",type, "max distace is ", max(my.files_beta[[type]]$distance[which(my.files_beta[[type]]$Stage==timepoint)])))
  
  
}


}




# 
# 
# 
#            bray_distance[[2]],
#            nrow=1,ncol=2, labels=c("A.","B."),common.legend = FALSE,  font.label = list(size = 30))
# 
# graphics.off()



#here
########################################################################################################################
#Look at beta diversity in liquid
########################################################################################################################

metacache_liquid<- metacache_total[which(metacache_total$sample_id %in% 
                               kefir4all_metadata$merge_column[
                                 which(kefir4all_metadata$`kefir type` %in% c("ML","WL"))
                               ]), ]



species_data_liquid <- c()

species_data_liquid[["Milk.kefir"]] <- metacache_liquid[which(metacache_liquid$sample_id %in% 
                                                  kefir4all_metadata$merge_column[
                                                    which(kefir4all_metadata$`kefir type` %in% c("ML"))
                                                  ]), ] |> column_to_rownames("sample_id")

species_data_liquid[["Water.kefir"]] <- metacache_liquid[which(metacache_liquid$sample_id %in% 
                                                   kefir4all_metadata$merge_column[
                                                     which(kefir4all_metadata$`kefir type` %in% c("WL"))
                                                   ]), ] |> column_to_rownames("sample_id")



type <- "Water.kefir"
liquid_name <- c()
bray_distance_liquid <- c()
my.files_beta_liquid <- c()
for (type in names(species_data_liquid)){
  
  if(type=="Milk.kefir"){
    liquid_name="ML"
  }else{
    liquid_name="WL"
  }
  
  # Get betadiversity distnace between To grain and grain of the other of the participant
  my.files_beta_liquid[[type]]<- as.matrix(vegdist(species_data_liquid[[type]], method = "bray"))
  
  
  my.files_beta_liquid[[type]] <- my.files_beta_liquid[[type]][which(rownames(my.files_beta_liquid[[type]])==kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                                                                                                kefir4all_metadata$`kefir type` %in% liquid_name)][1]),]
  
  
  # Get betadiversity between samples of the sample participant
  
  
  
  # 
  # 
  
  my.files_beta_liquid[[type]] <- 
    as.data.frame(my.files_beta_liquid[[type]] ) |> rownames_to_column("merge_column")
  colnames(my.files_beta_liquid[[type]])[2] <- "distance"
  
  
  
  
  my.files_beta_liquid[[type]] <- 
    merge(kefir4all_metadata,
          my.files_beta_liquid[[type]] , by="merge_column",all.y=TRUE)
  
  
  my.files_beta_liquid[[type]]$type <- gsub("\\."," ",type)
  
  
  bray_distance_liquid[[type]] <- 
    
    
    
    ggbetweenstats(data = my.files_beta_liquid[[type]] |> 
                     filter(Stage!="T0"),
                   x=Stage, 
                   y=distance,
                   title=type,
                   type = "nonparametric", # A
                   ggsignif.args    = list(textsize = 2.5, tip_length = 0.01)
                   
    )+
    theme_bw()+
    xlab("Timepoints")+
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
      #
  
  for (timepoint in levels(as.factor(my.files_beta_liquid[[type]]$Stage))){
    
    
    print(paste(timepoint,"-",type, "min distace is ", min(my.files_beta_liquid[[type]]$distance[which(my.files_beta_liquid[[type]]$Stage==timepoint)])))
    
    
    print(paste(timepoint,"-",type, "mean distace is ", mean(my.files_beta_liquid[[type]]$distance[which(my.files_beta_liquid[[type]]$Stage==timepoint)])))
    
    print(paste(timepoint,"-",type, "max distace is ", max(my.files_beta_liquid[[type]]$distance[which(my.files_beta_liquid[[type]]$Stage==timepoint)])))
    
    
  }
  
  
}


########################################################################################################################
#kefir grain stability over time 
########################################################################################################################



#get appropriate dataframe
########################################################################################################################


beta_diversity <- c()
total_beta_diversity <- c()

for (type in c("Water.kefir","Milk.kefir")){
  
  for (kefir_type in c("Liquid","Grain")){

    if(kefir_type=="Liquid"){
beta_diversity <-
as.matrix(vegdist(species_data_liquid[[type]], method = "bray"))
    }else{
      beta_diversity <-
        as.matrix(vegdist(species_data[[type]], method = "bray"))
}
  



beta_diversity <- as.data.frame(beta_diversity) |> rownames_to_column("sample_to") |> pivot_longer(!sample_to,names_to="samples_from",values_to = "beta_diversity_values")


beta_diversity <- 
merge(
#data.frame(sample_id=as.character(rownames(beta_diversity ))),
  beta_diversity,
           kefir4all_metadata,
           by.x="sample_to",
           by.y="merge_column",
           all.x=TRUE)


colnames(beta_diversity ) <- c("sample_to","samples_from", "beta_diversity_values", "id.to" ,"kefir type.to", "Stage.to", "data_source.to" ,"id_old.to",       "observations.to",        "category_confirmed.to")

beta_diversity <- 
  merge(
    #data.frame(sample_id=as.character(rownames(beta_diversity ))),
    beta_diversity,
    kefir4all_metadata,
    by.x="samples_from",
    by.y="merge_column",
    all.x=TRUE)


colnames(beta_diversity ) <- c("sample_to","samples_from", "beta_diversity_values", "id.to" ,"kefir type.to", "Stage.to", "data_source.to" ,"id_old.to",       "observations.to",        "category_confirmed.to","id.from" ,"kefir type.from", "Stage.from", "data_source.from" ,"id_old.from",       "observations.from",        "category_confirmed.from")




beta_diversity <-  beta_diversity[which(beta_diversity$id.to==beta_diversity$id.from),]
                                       

beta_diversity <- beta_diversity[-c(which(beta_diversity$Stage.to=="T0")),]

total_beta_diversity <- rbind(total_beta_diversity,beta_diversity)

 }# second for loop
  }# first for loop 



########################################################################################################################
# plot total comparisions with beta diversity

grain_bray_data=  total_beta_diversity |> 
  filter(`kefir type.to` %in% c("MG","WG")) |> 
  mutate(`kefir type.to`=gsub("MG","Milk kefir",
                              gsub("WG","Water kefir",`kefir type.to`))) 



grain_bray_data =grain_bray_data [-c(
                                                which(grain_bray_data$id.to %in% 
                                                                          names(which(
                                                                                    table(grain_bray_data$id.to)==1)))),]



pc <- 

  ggplot(grain_bray_data,aes(x = Stage.to , y =Stage.from)) +
  geom_tile(aes(fill = beta_diversity_values),#color = "white",
            #lwd = 1.5,
            linetype = 1) +
  labs(x="", y="", title="")+ #y="Feature"
  facet_wrap(~`kefir type.to`+id.to ,scales="free")+
  
    scale_fill_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
                       breaks= c(0,0.10,0.30,0.70,0.90,1),
                       values= rescale(as.numeric(c(0,0.10,0.30,0.70,0.90,1))),
                       guide="colorbar",
                       name="Bray Curtis dissimilarity" ,
                       labels=c(0,0.10,0.30,0.70,0.90,1))+
  theme_bw()+
  theme(#plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5),size = 55), #element_text(color="red", size=14, face="bold",hjust = 0.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(size=12.5),
    axis.text.y = element_text(size=12.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold"),
    #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
    axis.title.y = element_text( size=10, face="bold",hjust = 0.5, vjust = 1.5),
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm'), #change legend key width
    legend.text = element_text(size=20,face = "italic"),
    legend.title = element_text(size=20, vjust = 0.5))




########################################################################################################################
# now order the dataframe with the comparisons wk01 to wk05, wk05 to wk09,wk13,wk17,wk21

# remove single occurrences IDS 
total_beta_diversity <-
total_beta_diversity[-c(which(total_beta_diversity$id.to %in% c("ID013","ID079","ID096","ID069"))),]



total_beta_diversity_temporal <-
total_beta_diversity[c(which(total_beta_diversity$Stage.to=="wk01" &
                             total_beta_diversity$Stage.from=="wk05"),
                       which(total_beta_diversity$Stage.to=="wk05" &
                          total_beta_diversity$Stage.from=="wk09"),
                     which(total_beta_diversity$Stage.to=="wk09" &
                        total_beta_diversity$Stage.from=="wk13"),
                     which(total_beta_diversity$Stage.to=="wk13" &
                             total_beta_diversity$Stage.from=="wk17"),
                     which(total_beta_diversity$Stage.to=="wk17" &
                             total_beta_diversity$Stage.from=="wk21")
                       
                       ),]

  

########################################################################################################################
# plot the reordered the dataframe with the comparisons wk01 to wk05, wk05 to wk09,wk13,wk17,wk21


pd <- 

  total_beta_diversity_temporal |> 
  filter(`kefir type.to` %in% c("MG","WG")) |> 
  mutate(`kefir type.to`=gsub("MG","Milk kefir",
                              gsub("WG","Water kefir",`kefir type.to`))) |> 
  
  ggplot(aes(x = id.to , y =Stage.from)) +
  geom_tile(aes(fill = beta_diversity_values),#color = "white",
            #lwd = 1.5,
            linetype = 1) +
  labs(x="", y="", title="")+ #y="Feature"
  facet_wrap(~`kefir type.to` ,scales="free")+
  
    scale_fill_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
                       breaks= c(0,0.10,0.30,0.70,0.90,1),
                       values= rescale(as.numeric(c(0,0.10,0.30,0.70,0.90,1))),
                       guide="colorbar",
                       name="Bray Curtis dissimilarity" ,
                       labels=c(0,0.10,0.30,0.70,0.90,1))+
  theme_bw()+
  theme(#plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5),size = 55), #element_text(color="red", size=14, face="bold",hjust = 0.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(size=12.5),
    axis.text.y = element_text(size=12.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold"),
    #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
    axis.title.y = element_text( size=10, face="bold",hjust = 0.5, vjust = 1.5),
    legend.key.size = unit(2, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm'), #change legend key width
    legend.text = element_text(size=20,face = "italic"),
    legend.title = element_text(size=20, vjust = 0.5))




########################################################################################################################
# final plot 


library(ggpubr)



# Wrapped in tryCatch: ggbetweenstats panels embed a directlabels layer that
# can fail headless when a stage group is degenerate (e.g. single pair,
# min==max), which would otherwise halt the whole script before the
# Extended Data Figure 1 output below runs.
tryCatch(
  {
    tiff(filename=file.path(FIGURES_DIR, 'Figure_3.tiff'), width = 1800, height = 1200, units = "px", res = 300, compression = "lzw", type = "cairo", pointsize = 12) #, width=2000, height=1950)

    ggarrange(
      ggarrange( bray_distance[[1]],
                 bray_distance[[2]],
                 nrow=1,ncol=2, labels=c("A.","B."),common.legend = FALSE,  font.label = list(size = 30)),
      pc,
      nrow=2,ncol=1, labels=c("","C."),common.legend = FALSE,  font.label = list(size = 30),heights=c(.4,1))

    graphics.off()
  },
  error = function(e) {
    try(graphics.off(), silent = TRUE)
    message("V2: legacy Figure_3.tiff render skipped (", conditionMessage(e), ")")
  }
)

# The stability/prevalence/compositional analysis below is not required to
# build Extended Data Figure 1 (the mean_compositional_change and
# mean_gower_change panels it needs are already built above). Guarded in a
# tryCatch so a failure anywhere in this section cannot block that figure.
tryCatch(
  {

########################################################################################################################
#differences in stability over time
########################################################################################################################
title_data <- 
data.frame(type=as.character(c("MG","WG","ML","WL")),
                             title=as.character(c("Milk kefir grain", "Water kefir grain", "Milk kefir liquid","Water kefir liquid")))

t <- c()
title <- c()
mean_compositional_change <- c()



total_beta_diversity_temporal$same_time <- paste(total_beta_diversity_temporal$Stage.to," - ", total_beta_diversity_temporal$Stage.from,sep="")




for (type in c("MG","WG","ML","WL")){
  
 title <- title_data$title[title_data$type==type]
 
 mean_compositional_change[[type]] <- 
ggbetweenstats(data = total_beta_diversity_temporal[which(total_beta_diversity_temporal$`kefir type.to`==type),],
               x=same_time, 
               y=beta_diversity_values,
               title=paste(title, " - metagenomes",sep=""),
               type = "nonparametric", # A
               ggsignif.args    = list(textsize = 2, tip_length = 0.01))+
                 theme_bw()+
                 xlab("Timepoints")+
                 ylab("Bray Curtis distance")+
                 #labs(fill = "Timepoint")+
                 #guides(colour = guide_legend(override.aes = list(size=25)))+
                 #ylim(0,250)+
                 theme(#plot.title = element_text(hjust = 0.5,size=35,face="bold"),
                   legend.title = element_text( size=25, face="bold"),
                   axis.text.x = element_text(hjust = 1, size = 15),#angle = 45
                   axis.text.y = element_text(hjust = 1, size = 10),
                   axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
                   axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
                   legend.position = "none",
                   plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
                                                               box.color = "#748696",size=35, lineheight = 2))
               


}


#write.csv(total_beta_diversity_temporal, file.path(DATA_DIR, "bray_curtis_values_same_participant.csv"), quote = FALSE,row.names=FALSE )



pacman::p_load(rstatix)


t <- 

total_beta_diversity_temporal |>
  group_by(`kefir type.to`) |>
  dunn_test( beta_diversity_values ~  same_time)




                                                           ########################################################################################################################
                                                                                   #Gower metadata distance
                                                           ########################################################################################################################



########################################################################################################################
# Gower distance matrix

########################################################################################################################
# .libPaths() removed â€” use default R library or renv

########################################################################################################################
#Libraries used
########################################################################################################################
#install.packages("blastula")
library(readxl)
library(rlang)
library(ggplot2)
library(grid)
library(dplyr)
library(tidyr)

#devtools::install_github("tntp/surveymonkey")

# surveymonkey requires a SurveyMonkey OAuth token (not distributable).
# Install with: devtools::install_github("tntp/surveymonkey")
# Sections below that call surveymonkey functions will be skipped if not available.
if (requireNamespace("surveymonkey", quietly=TRUE)) library(surveymonkey)
# Import the options command into your r_profile using the below command
#usethis::edit_r_profile()
#options(sm_oauth_token = "vOcHPmfepVilX6cUC5JVr4axcPfCv4UxPKhWEd-rbegkTh58EF2I.SvzlZYBpNnQliGeObex30.GNzFQP8lUb.UVzWMjHVeli4G-pgueVnnuhu66IEzB5ovNjlYmP1YU")
getOption("sm_oauth_token")




library(readr)
SURVEY_MK <- file.path(PRIVATE_DIR, "fermentation_survey_mk_edited_dated.txt")
SURVEY_WK <- file.path(PRIVATE_DIR, "fermentation_survey_wk_edited_dated.txt")
if (!file.exists(SURVEY_MK) || !file.exists(SURVEY_WK)) {
  stop(
    "This script requires private fermentation survey files not distributed with the repo.\n",
    "Missing: ", SURVEY_MK, "\n",
    "Missing: ", SURVEY_WK, "\n",
    "See data/private/README.md for details."
  )
}
survey_mk <- read_delim(SURVEY_MK, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
survey_wk <- read_delim(SURVEY_WK, delim = "\t", escape_double = FALSE, trim_ws = TRUE)




survey_mk$previous_weightofgrains <- iconv(survey_mk$previous_weightofgrains, to = "UTF-8", sub = "")
survey_mk$previous_weightofgrains<- gsub("g.*|G.*", "",survey_mk$previous_weightofgrains)



survey_mk$previous_weightofgrains<- as.numeric(gsub("[^0-9.-]", "",survey_mk$previous_weightofgrains))

survey_wk$previous_weightofgrains <- iconv(survey_wk$previous_weightofgrains, to = "UTF-8", sub = "")
survey_wk$previous_weightofgrains<- as.numeric(gsub("[^0-9.-]", "",survey_wk$previous_weightofgrains))

survey_wk$previous_weightofgrains[
  which(survey_wk$previous_weightofgrains==
          "Don't have an exact measurment but about 1.2X")] <- NA




survey_mk$stage[which(is.na(survey_mk$stage))] <- "wk21"

survey_wk$stage[which(is.na(survey_wk$stage) &
                        survey_wk$Participant_ID=="ID85")] <- "wk13"


survey_wk$stage[which(is.na(survey_wk$stage))] <- "wk21"



survey_wk$sample_id <-  paste(survey_wk$Participant_ID,survey_wk$Fermentation_date,sep=" - ")
survey_mk$sample_id <-  paste(survey_mk$Participant_ID,survey_mk$Fermentation_date,sep=" - ")




survey_wk<- 
  survey_wk[-c(
    which(duplicated(survey_wk$sample_id))),]


survey_mk<- 
  survey_mk[-c(
    which(duplicated(survey_mk$sample_id))),]





pacman::p_load(StatMatch)
library(StatMatch)
gower_mk<-
  gower.dist(
    as.data.frame( unclass( dplyr::select(survey_mk,-c(final_comments,Fermentation_date,Participant_ID, School))),stringsAsFactors=TRUE)
  )


rownames(gower_mk) <- survey_mk$sample_id
colnames(gower_mk) <- survey_mk$sample_id



#gower_mk[which(rownames(gower_mk)=="ID9_25/04/2022"),which(colnames(gower_mk)=="ID9_25/04/2022")]


gower_wk<-
  gower.dist(
    as.data.frame( unclass( dplyr::select(survey_wk,-c(final_comments,Fermentation_date,Participant_ID, School))),stringsAsFactors=TRUE)
  )



rownames(gower_wk) <- survey_wk$sample_id
colnames(gower_wk) <- survey_wk$sample_id



library(dplyr)
library(tibble)
library(scales)

gower_wk <- 
  as.data.frame(gower_wk) |> rownames_to_column("sample_id") |> pivot_longer(!sample_id,names_to = "ID_to",values_to = "distance" )

gower_mk <- 
  as.data.frame(gower_mk) |> rownames_to_column("sample_id") |> pivot_longer(!sample_id,names_to = "ID_to",values_to = "distance" )




nrow(gower_wk)

gower_wk <- 
  merge(gower_wk, 
        dplyr::select(survey_wk,sample_id,Participant_ID,Fermentation_date,stage ),
        by="sample_id",
        all.x=TRUE
  )

nrow(gower_wk)




nrow(gower_wk)

gower_wk <- 
  merge(gower_wk, 
        dplyr::select(survey_wk,sample_id,Participant_ID,Fermentation_date,stage ),
        by.x="ID_to",
        by.y="sample_id",
        all.x=TRUE
  )

nrow(gower_wk)

nrow(gower_mk)

gower_mk <- 
  merge(gower_mk, 
        dplyr::select(survey_mk,sample_id,Participant_ID,Fermentation_date,stage ),
        by="sample_id",
        all.x=TRUE
  )

nrow(gower_mk)


nrow(gower_mk)

gower_mk <- 
  merge(gower_mk, 
        dplyr::select(survey_mk,sample_id,Participant_ID,Fermentation_date,stage ),
        by.x="ID_to",
        by.y="sample_id",
        all.x=TRUE
  )

nrow(gower_mk)



#note .x will refer to the sample_id
# .y will refer to the ID_to





gower_dist=c()

gower_dist[["Milk_kefir"]] <- gower_mk

gower_dist[["Water_kefir"]] <- gower_wk

mean_gower_change<-c() 

if (!requireNamespace("ggstatsplot", quietly=TRUE)) install.packages("ggstatsplot")
library(ggstatsplot)

for (type in names(gower_dist)){
  
  
  
  gower_dist[[type]]$stage_comparision <- paste(gower_dist[[type]]$stage.x, gower_dist[[type]]$stage.y,sep=" - ")
  
  
  gower_dist[[type]] <-   gower_dist[[type]][which(gower_dist[[type]]$stage_comparision %in% c("wk01 - wk05", "wk05 - wk09","wk09 - wk13", "wk13 - wk17", "wk17 - wk21")),]
  
  
  
  
  
  mean_gower_change[[type]] <- 
    
    
    ggbetweenstats(data = gower_dist[[type]][which(gower_dist[[type]]$Participant_ID.x==gower_dist[[type]]$Participant_ID.y),],
                   x=stage_comparision, 
                   y=distance,
                   title=paste(gsub("_"," ",type), "- fermenation surveys"),
                   type = "nonparametric", # A
                   ggsignif.args    = list(textsize = 2, tip_length = 0.01))+
    theme_bw()+
    xlab("Timepoints")+
    ylab("Gower distance")+
    #labs(fill = "Timepoint")+
    #guides(colour = guide_legend(override.aes = list(size=25)))+
    #ylim(0,250)+
    theme(#plot.title = element_text(hjust = 0.5,size=35,face="bold"),
      legend.title = element_text( size=25, face="bold"),
      axis.text.x = element_text( hjust = 1, size = 15), #angle = 45
      axis.text.y = element_text(hjust = 1, size = 10),
      axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
      axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
      legend.position = "none",
      plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
                                                  box.color = "#748696",size=35, lineheight = 2))
  
  
  
  
  
}


library(ggpubr)

# This downstream block produces Extended Data Figure 1 (the
# mean_compositional_change / mean_gower_change composite below).
#
# Rendered with Cairo::CairoTIFF() rather than base R's tiff(type = "cairo"):
# the base bitmap backend is not reliably available in this container, while
# the Cairo package's own device is. Wrapped in tryCatch with an explicit
# file-existence/size check so a silently empty or missing output file fails
# loudly instead of looking like a successful run.
#
# Canvas is 3600x2400: the panels' large title boxes
# (theme(plot.title = ggtext::element_textbox_simple(size = 35, ...))),
# combined with ggbetweenstats()'s own built-in subtitle, need this much room
# in a 3x2 grid -- a smaller canvas collapses a grid viewport to zero
# dimensions (reproduces identically with base png(), so it is not specific
# to the Cairo device).
supp_fig_2_path <- file.path(FIGURES_DIR, "Extended_Data_Figure_1.tiff")
supp_fig_2_plot <- ggarrange(
  mean_compositional_change[[1]],
  mean_compositional_change[[2]],
  mean_compositional_change[[3]],
  mean_compositional_change[[4]],
  mean_gower_change[[1]],
  mean_gower_change[[2]],
  nrow = 3, ncol = 2, labels = c("A.", "B.", "C.", "D.", "E.", "F."),
  common.legend = FALSE, font.label = list(size = 30))
tryCatch({
  Cairo::CairoTIFF(file = supp_fig_2_path, width = 3600, height = 2400,
                    dpi = 300, pointsize = 12, bg = "white")
  print(supp_fig_2_plot)
  grDevices::dev.off()
  if (!file.exists(supp_fig_2_path) || file.info(supp_fig_2_path)$size < 500) {
    stop("file was not written or is suspiciously small")
  }
  message("Wrote Extended_Data_Figure_1.tiff to ", FIGURES_DIR)
  # 24x16 in: save_figure()'s vector PDF leg (ggsave(..., device = cairo_pdf))
  # needs a larger canvas than its raster PNG leg to clear the same
  # zero-dimension-viewport issue described above for the main TIFF -- this
  # size comfortably covers both legs' thresholds.
  save_figure(supp_fig_2_plot, file.path(FIGURES_DIR, "Extended_Data_Figure_1"),
              width_mm = 24 * 25.4, height_mm = 16 * 25.4,
              source = PRIMARY_SOURCE, seed = NULL)
}, error = function(e) {
  try(grDevices::dev.off(), silent = TRUE)
  stop("Failed to write Extended_Data_Figure_1.tiff: ", conditionMessage(e))
})
  },
  error = function(e) {
    message("Non-essential downstream analysis failed (Extended Data Figure 1 is unaffected): ", conditionMessage(e))
  }
)
