
########################################################################################################################
#Libraries used 
########################################################################################################################


.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")
pacman::p_load(readxl,tidyr,readr,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,taxizedb )


########################################################################################################################
#Libraries used
########################################################################################################################

.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")
pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis)
pacman::p_load(readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)
library(vegan)
library(ggplot2)
library(grid)
library(ggstatsplot)


########################################################################################################################
#Import metacache data
########################################################################################################################
metacache <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_species_profile_v2.csv")
metacache_strain <-  read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_strain_profile_v2.csv")

metacache$sample_id <- gsub("-","_",metacache$sample_id)
metacache$sample_id <- gsub(".*/","",metacache$sample_id)

#metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]

Citizen_Scientist_metadata_v8 <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/Citizen Scientist metadata_v8.csv")

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )

########################################################################################################################
#Import and modify survey_responses 
########################################################################################################################
setwd("Q:/H2020 Master/Citizen Science Project/Results/00/")
temp = list.files(pattern=".xlsx", recursive = FALSE)
myfiles = lapply(temp,read_excel)
names(myfiles) <- gsub("file_names_survey_responses_|.xlsx","",temp)


kefir4all_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/kefir4all_sample_metadata_v2.csv")

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
                               ]), ] %>% column_to_rownames("sample_id")

species_data[["Water.kefir"]] <- metacache[which(metacache$sample_id %in% 
                                                  kefir4all_metadata$merge_column[
                                                    which(kefir4all_metadata$`kefir type` %in% c("WG"))
                                                  ]), ] %>% column_to_rownames("sample_id")



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


#)

# all_test <- 
#   adonis2(my.files_beta[[type]]~Stage, data=kefir4all_metadata ,permutations =  10000)
# 
# 

my.files_beta[[type]] <- 
as.data.frame(my.files_beta[[type]] ) %>% rownames_to_column("merge_column")
colnames(my.files_beta[[type]])[2] <- "distance"




my.files_beta[[type]] <- 
  merge(kefir4all_metadata,
        my.files_beta[[type]] , by="merge_column",all.y=TRUE)


my.files_beta[[type]]$type <- gsub("\\."," ",type)


bray_distance[[type]] <- 
  
  
  
  ggbetweenstats(data = my.files_beta[[type]] %>% 
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


# my.files_beta[[type]] %>% 
#   filter(Stage!="T0")%>% 
#   ggplot(aes(x=Stage, y=as.numeric(distance),fill=Stage)) +
#   geom_boxplot() +
#   #facet_wrap(~ medium)+
#   #labs(title= 'Alpha diversity of timepoints') +
#   geom_point(position=position_dodge(width=0.75),aes(group=Stage))+
#   theme_bw()+
#   xlab("Time points")+
#   ylab("Bray Curtis distance")+
#   labs(fill = "Timepoint")+
#   #guides(colour = guide_legend(override.aes = list(size=25)))+
#   #ylim(0,250)+
#   theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
#         legend.title = element_text( size=25, face="bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
#         axis.text.y = element_text(hjust = 1, size = 10),
#         axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
#         axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
#         legend.position = "none")
# #

for (timepoint in levels(as.factor(my.files_beta[[type]]$Stage))){
  
  
  print(paste(timepoint,"-",type, "min distace is ", min(my.files_beta[[type]]$distance[which(my.files_beta[[type]]$Stage==timepoint)])))
  
  
  print(paste(timepoint,"-",type, "mean distace is ", mean(my.files_beta[[type]]$distance[which(my.files_beta[[type]]$Stage==timepoint)])))
 
  print(paste(timepoint,"-",type, "max distace is ", max(my.files_beta[[type]]$distance[which(my.files_beta[[type]]$Stage==timepoint)])))
  
  
}


}




# 
# jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/Evolution/Figure compositional_change.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
# 
# 
# ggarrange( bray_distance[[1]],
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
                                                  ]), ] %>% column_to_rownames("sample_id")

species_data_liquid[["Water.kefir"]] <- metacache_liquid[which(metacache_liquid$sample_id %in% 
                                                   kefir4all_metadata$merge_column[
                                                     which(kefir4all_metadata$`kefir type` %in% c("WL"))
                                                   ]), ] %>% column_to_rownames("sample_id")



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
  
  
  #)
  
  # all_test <- 
  #   adonis2(my.files_beta[[type]]~Stage, data=kefir4all_metadata ,permutations =  10000)
  # 
  # 
  
  my.files_beta_liquid[[type]] <- 
    as.data.frame(my.files_beta_liquid[[type]] ) %>% rownames_to_column("merge_column")
  colnames(my.files_beta_liquid[[type]])[2] <- "distance"
  
  
  
  
  my.files_beta_liquid[[type]] <- 
    merge(kefir4all_metadata,
          my.files_beta_liquid[[type]] , by="merge_column",all.y=TRUE)
  
  
  my.files_beta_liquid[[type]]$type <- gsub("\\."," ",type)
  
  
  bray_distance_liquid[[type]] <- 
    
    
    
    ggbetweenstats(data = my.files_beta_liquid[[type]] %>% 
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
  



beta_diversity <- as.data.frame(beta_diversity) %>% rownames_to_column("sample_to") %>% pivot_longer(!sample_to,names_to="samples_from",values_to = "beta_diversity_values")


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

grain_bray_data=  total_beta_diversity %>% 
  filter(`kefir type.to` %in% c("MG","WG")) %>% 
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
                       breaks= c(-0.05,0,0.10,0.30,0.70,0.90,1),
                       values= rescale(as.numeric(c(-0.05,0,0.10,0.30,0.70,0.90,1))),
                       guide="colorbar",
                       name="Bray Curtis dissimilarity" ,
                       labels=c(-0.05,0,0.10,0.30,0.70,0.90,1))+
  theme_bw()+
  theme(#plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5),size = 55), #element_text(color="red", size=14, face="bold",hjust = 0.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(size=12.5),
    axis.text.y = element_text(size=12.5),
    strip.background = element_rect(
      color="black", fill="white"),
    strip.text = element_text(size=10),
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

  total_beta_diversity_temporal %>% 
  filter(`kefir type.to` %in% c("MG","WG")) %>% 
  mutate(`kefir type.to`=gsub("MG","Milk kefir",
                              gsub("WG","Water kefir",`kefir type.to`))) %>% 
  
  ggplot(aes(x = id.to , y =Stage.from)) +
  geom_tile(aes(fill = beta_diversity_values),#color = "white",
            #lwd = 1.5,
            linetype = 1) +
  labs(x="", y="", title="")+ #y="Feature"
  facet_wrap(~`kefir type.to` ,scales="free")+
  
  scale_fill_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
                       breaks= c(-0.05,0,0.10,0.30,0.70,0.90,1),
                       values= rescale(as.numeric(c(-0.05,0,0.10,0.30,0.70,0.90,1))),
                       guide="colorbar",
                       name="Relative abundance  (%)" ,
                       labels=c(-0.05,0,0.10,0.30,0.70,0.90,1))+
  theme_bw()+
  theme(#plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5),size = 55), #element_text(color="red", size=14, face="bold",hjust = 0.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(size=12.5),
    axis.text.y = element_text(size=12.5),
    strip.background = element_rect(
      color="black", fill="white"),
    strip.text = element_text(size=10),
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



jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 3_model.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)

 ggarrange(
 ggarrange( bray_distance[[1]],
            bray_distance[[2]],
            nrow=1,ncol=2, labels=c("A.","B."),common.legend = FALSE,  font.label = list(size = 30)),
 pc,
 nrow=2,ncol=1, labels=c("","C."),common.legend = FALSE,  font.label = list(size = 30),heights=c(.4,1))
#

graphics.off()



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


#write.csv(total_beta_diversity_temporal, "Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/bray_curtis_values_same_participant.csv", quote = FALSE,row.names=FALSE )



pacman::p_load(rstatix)


t <- 

total_beta_diversity_temporal %>%
  group_by(`kefir type.to`) %>%
  dunn_test( beta_diversity_values ~  same_time)




                                                           ########################################################################################################################
                                                                                   #Gower metadata distance
                                                           ########################################################################################################################



########################################################################################################################
# Gower distance matrix

########################################################################################################################
.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")

########################################################################################################################
#Libraries used
########################################################################################################################
#install.packages("blastula")
pacman::p_load(readxl,rlang,ggplot2,jpeg,patchwork,png,grid,ggimage,dplyr,tidyr,blastula)

#devtools::install_github("tntp/surveymonkey")

library(surveymonkey)
# Import the options command into your r_profile using the below command
#usethis::edit_r_profile()
#options(sm_oauth_token = "vOcHPmfepVilX6cUC5JVr4axcPfCv4UxPKhWEd-rbegkTh58EF2I.SvzlZYBpNnQliGeObex30.GNzFQP8lUb.UVzWMjHVeli4G-pgueVnnuhu66IEzB5ovNjlYmP1YU")
getOption("sm_oauth_token")




library(readr)
survey_mk <- read_delim("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/fermentation_survey_mk_edited_dated.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)


survey_wk <- read_delim("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/fermentation_survey_wk_edited_dated.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)




survey_mk$previous_weightofgrains<- gsub("g.*|G.*", "",survey_mk$previous_weightofgrains)



survey_mk$previous_weightofgrains<- as.numeric(gsub("[^0-9.-]", "",survey_mk$previous_weightofgrains))


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
  as.data.frame(gower_wk) %>% rownames_to_column("sample_id") %>% pivot_longer(!sample_id,names_to = "ID_to",values_to = "distance" )

gower_mk <- 
  as.data.frame(gower_mk) %>% rownames_to_column("sample_id") %>% pivot_longer(!sample_id,names_to = "ID_to",values_to = "distance" )




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

jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/CS_metagenomics/manuscript figures/supplementary figure 2_v2.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)

ggarrange(
  
           mean_compositional_change[[1]],
           mean_compositional_change[[2]],
           mean_compositional_change[[3]],
           mean_compositional_change[[4]],
           mean_gower_change[[1]],
           mean_gower_change[[2]],
                    nrow=3,ncol=2, labels=c("A.","B.","C.","D.","E.","F."),common.legend = FALSE,  font.label = list(size = 30))

graphics.off()


# identify outliers metagenomes with a higher dissimilarity value

total_beta_diversity_temporal[which(total_beta_diversity_temporal$`kefir type.to`=="MG" &
                                      total_beta_diversity_temporal$`Stage.from`=="wk21"),]


total_beta_diversity_temporal[which(total_beta_diversity_temporal$`kefir type.to`=="WG" &
                                      total_beta_diversity_temporal$`Stage.from`=="wk21"),]













































########################################################################################################################
#old junk script
########################################################################################################################




my.files_beta[["Milk.kefir"]]<- as.matrix(vegdist(species_data[["Milk.kefir"]], method = "bray"))

my.files_beta[["Milk.kefir"]] <- 
  as.data.frame(my.files_beta[["Milk.kefir"]] ) %>% rownames_to_column("merge_column")





my.files_beta[["Milk.kefir"]] <- 
merge(kefir4all_metadata,
      my.files_beta[["Milk.kefir"]] %>% pivot_longer(!merge_column,names_to = "sample",values_to = "distance"), by="merge_column",all.y=TRUE) 
  


  my.files_beta[["Milk.kefir"]] <- 
  merge(kefir4all_metadata,
        my.files_beta[["Milk.kefir"]],by.x="merge_column", by.y="sample",all.y=TRUE)
  
  



my.files_beta[["Milk.kefir"]] <- my.files_beta[["Milk.kefir"]][
  which(my.files_beta[["Milk.kefir"]]$Sample.x== my.files_beta[["Milk.kefir"]]$Sample.y),]


my.files_beta[["Milk.kefir"]]$x.axis <- paste(
  my.files_beta[["Milk.kefir"]]$Sample.x, 
  my.files_beta[["Milk.kefir"]]$Stage.x)

my.files_beta[["Milk.kefir"]]$y.axis <- paste(
  my.files_beta[["Milk.kefir"]]$Sample.y, 
  my.files_beta[["Milk.kefir"]]$Stage.y)


ggplot(my.files_beta[["Milk.kefir"]] ,aes(x=x.axis, y=y.axis,fill=as.numeric(distance))) +
  geom_tile() +
  scale_fill_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
                       breaks= c(1,.75,.5,.25,0,-5),
                       values= rescale(as.numeric(c(1,.75,.5,.25,0,-5))),
                       guide="colorbar",
                       name="")+
  
  facet_wrap(~Sample.x,scales = "free")
  #facet_wrap(~ medium)+
  #labs(title= 'Alpha diversity of timepoints') +
 # geom_point(position=position_dodge(width=0.75),aes(group=Stage))+
  theme_bw()+
  xlab("Kefir grain Samples")+
  ylab("Kefir grain Samples")+
  labs(fill = "Timepoint")+

    
  #guides(colour = guide_legend(override.aes = list(size=25)))+
  #ylim(0,250)+
  theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
        legend.title = element_text( size=25, face="bold"),
        axis.text.x = element_text( hjust = 1, size = 15),#angle = 45,
        axis.text.y = element_text(hjust = 1, size = 10),
        axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
        axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
        legend.position = "none")




dist_mi <- 1/dist_m # one over, as qgraph takes similarity matrices as input

rownames(as.data.frame(dist_mi))

library(qgraph)
qgraph(dist_mi, layout='spring', vsize=3)
     
#n1_species <-which(maxab_species < 0.1)

#https://www.statology.org/hierarchical-clustering-in-r/

library(cluster)
#define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
  agnes(species_data[["Milk.kefir"]] , method = x)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac)

library(ggtree)
clust <- agnes(species_data[["Milk.kefir"]], method = "ward")


murders_dist <- dist(species_data[["Milk.kefir"]], method="euclidean")  
fit <- hclust(murders_dist, method="ward.D") 

plot(fit, family="Arial")
rect.hclust(fit, k=4, border="cadetblue")

library(ggdendro)
ggdendrogram(fit, rotate = TRUE, theme_dendro = FALSE) +
  theme_minimal() + xlab("") + ylab("")



ggtree(as.hclust(clust))

pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram") 

gap_stat <- clusGap(species_data[["Milk.kefir"]], FUN = hcut, nstart = 25, K.max = 10, B = 50)

library(factoextra)
fviz_gap_stat(gap_stat)




###################
#Prevalence info
specie <- c()
sum_specie <- c()
species_prevalence <- c()


species_profile <- c()
species_profile[["water"]] <- metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL"))]),]#& 
#kefir4all_metadata$Stage!="T0")]),]

species_profile[["milk"]] <- metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML"))]),] 



#species_profile[["water"]] <- metacache[which(rownames(