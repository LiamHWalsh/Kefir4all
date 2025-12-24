
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
#Make a metadata file ONLY NEED TO RUN THIS ONCE
########################################################################################################################
# kefir4all_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/00/file_names.csv")
# 
# # RMEOVE LANE NAMES IN sample code AS IT WAS DESIGNED TO RENAME FILES ON THE CLUSTER PRIOR TO COMPUTATIONAL ANLAYSIS 
# 
# kefir4all_metadata$new_name <- gsub("_L00.","",kefir4all_metadata$new_name)
# 
# kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$new_name))),]
# 
# kefir4all_metadata <- cbind(kefir4all_metadata,
# str_split_fixed(kefir4all_metadata$new_name,"_",5))
# 
# kefir4all_metadata <-
# dplyr::select(kefir4all_metadata,new_name,`2`,`3`,`4`)%>%
#   dplyr::rename(merge_column=new_name,
#          Sample=2,
#          `kefir type`=3,
#          Stage=4)
# 
# 
# 
# kefir4all_metadata[-c(grep("ID",kefir4all_metadata$Sample)),]$Stage <- "T0"
# kefir4all_metadata[-c(grep("ID",kefir4all_metadata$Sample)),]$Stage
# 
# kefir4all_metadata[which(kefir4all_metadata$Stage=="T0"),]$`kefir type`[grep("ML",kefir4all_metadata[which(kefir4all_metadata$Stage=="T0"),"merge_column"])] <- "ML"
# kefir4all_metadata[which(kefir4all_metadata$Stage=="T0"),]$`kefir type`[grep("MG",kefir4all_metadata[which(kefir4all_metadata$Stage=="T0"),"merge_column"])] <- "MG"
# kefir4all_metadata[which(kefir4all_metadata$Stage=="T0"),]$`kefir type`[grep("WL",kefir4all_metadata[which(kefir4all_metadata$Stage=="T0"),"merge_column"])] <- "WL"
# kefir4all_metadata[which(kefir4all_metadata$Stage=="T0"),]$`kefir type`[grep("WG",kefir4all_metadata[which(kefir4all_metadata$Stage=="T0"),"merge_column"])] <- "WG"
# kefir4all_metadata$data_source <- "This study"
# 
# write.csv(kefir4all_metadata, "Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/kefir4all_sample_metadata_v2.csv", quote = FALSE,row.names = FALSE)
#note manual fix samll errors like extraction control
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
kefir4all_metadata <- kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$merge_column))),]
########################################################################################################################
# Merge metadata into one
########################################################################################################################
total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`))
total_metadata <- 
total_metadata[-c(which(total_metadata$`kefir type`=="Medium control")),]


########################################################################################################################
#Get data from the global_MK study
########################################################################################################################

metacache <- metacache[which( metacache$sample_id %in% total_metadata$merge_column),]

########################################################################################################################
#PCOA used
########################################################################################################################


metacache [is.na(metacache )] <- 0

#metacache <- metacache[which(metacache$sample_id %in% 
#kefir4allmetadata$merge_column[grep("ID52", kefir4allmetadata$Sample)]), ]


maxab_species <- apply(dplyr::select(metacache, -sample_id),2, max, na.rm=TRUE)

#n1_species <-which(maxab_species < 0.1)


#metacache <- 
#metacache %>% column_to_rownames("sample_id")





########################################################################################################################
#PCOA used for total data
########################################################################################################################

#kefir4all_metadata$merge_column[-c(which(kefir4all_metadata$merge_column %in% metacache$sample_id[which(metacache$sample_id %in% kefir4all_metadata$merge_column)]))]

#metacache <- metacache[which(metacache$sample_id %in% kefir4all_metadata$merge_column),]

data.bray=vegdist(metacache %>% remove_rownames() %>% column_to_rownames("sample_id"),na.rm = TRUE) #Bray-Curtis distance

#data.bray=vegdist(metacache, na.rm = FALSE) #Bray-Curtis distance


metacache$sample_id[which(metacache$`Leclercia adecarboxylata`==max(metacache$`Leclercia adecarboxylata`))]

data.b.pcoa=cmdscale(data.bray,k=2,eig=TRUE,add = TRUE) #ordination



pcoa = data.frame(PC1 = data.b.pcoa$points[,1], PC2 = data.b.pcoa$points[,2])
percent_explained <- 100* data.b.pcoa$eig/sum(data.b.pcoa$eig)


p = ggplot(pcoa, aes(x=PC1, y=PC2)) + geom_point(size=3) + theme_bw()

#############################################################################################################################
#Clustering by dataset
#############################################################################################################################
pcoa <- merge(pcoa,total_metadata,by.x=0,by.y="merge_column",all.x=TRUE)

#############################################################################################################################
#plot by dataset
#############################################################################################################################

pcoa$`kefir type`[which(is.na(pcoa$`kefir type`))] <- "ML"

pcoa$data_source[which(is.na(pcoa$data_source))] <- "Walsh et al 2023"

pcoa$data_source_specific <- pcoa$data_source
pcoa$data_source_specific[which(pcoa$data_source=="This study" &
                                  pcoa$`kefir type` %in% c("WL","WG"))] <- "Kefir4All - Water kefir"


pcoa$data_source_specific[which(pcoa$data_source=="This study" &
                                  pcoa$`kefir type` %in% c("ML","MG"))] <- "Kefir4All - Milk kefir"

pcoa$data_source_specific[which(pcoa$data_source=="Mortensen et al 2023")] <- "Breselge et al - Water kefir"

pcoa$data_source_specific[which(pcoa$data_source=="Walsh et al 2023")] <- "Walsh et al - Milk kefir"

centroid <- data.frame(data_source_specific=as.character(levels(as.factor(pcoa$data_source_specific))),
                       PC1=as.numeric(0),
                       PC2=as.numeric(0))

for (i in centroid$data_source){
  centroid$PC1[which(centroid$data_source==i)] <- mean(pcoa$PC1[which(pcoa$data_source_specific ==i)])
  centroid$PC2[which(centroid$data_source==i)] <- mean(pcoa$PC2[which(pcoa$data_source_specific==i)])
  
}




#merge pcoa by community types 
# 
# p_total <- 
# 
# ggplot(pcoa, aes(PC1, y=PC2,shape=data_source_specific ,colour=`kefir type`))+
#   geom_point(size=5) +
#   geom_point(data=centroid,size=10,shape=21, color="black",aes(fill=data_source_specific),  show.legend=FALSE)+ #aes(, colour=data_source)),size=10)+
#   #geom_convexhull(alpha=.1)+
#   stat_ellipse(geom = "polygon",
#                aes(fill=data_source_specific),
#                alpha = 0.25,
#                type = "norm",
#                show.legend=TRUE)+
#   #geom_convexhull(alpha=.1)+
#   #stat_ellipse(geom = "polygon",
#   #alpha = 0.25,
#   #type = "norm")+
#   #geom_text(colour="blue", check_overlap = TRUE, size=2.5, 
#   #hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.025)+
#   # directlabels::geom_dl(data=labels, aes(label = species), method = "smart.grid")+
#   # Filter data first
#   #geom_segment(data=species.long3, 
#   #aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
#   #colour="red", size=0.7, arrow=arrow()) +
# 
# ########
# labs(x=paste("PCoA1 - ", round(percent_explained[1]), "%", sep=""), y=paste("PCoA2 - ", round(percent_explained[2]), "%", sep=""), title="") +
#   coord_equal() +
#   theme_bw()+
#   theme(legend.position = "right",#c(.85,.2),#axis.text.x = element_blank(),  # remove x-axis text
#     #axis.text.y = element_blank(), # remove y-axis text
#     axis.ticks = element_blank(),  # remove axis ticks
#     axis.text = element_blank(),
#     axis.title.x = element_text(size=20), # remove x-axis labels
#     axis.title.y = element_text(size=20), # remove y-axis labels
#     panel.background = element_blank(), 
#     panel.grid.major = element_blank(),  #remove major-grid labels
#     panel.grid.minor = element_blank(),  #remove minor-grid labels
#     plot.background = element_blank(),
#     legend.text=element_text(size = 20),
#     legend.key.size = unit(2, 'cm'), #change legend key size
#     legend.key.height = unit(2, 'cm'), #change legend key height
#     legend.key.width = unit(2, 'cm'), #change legend key width
#     legend.title = element_text(size=20), #change legend title font size
#     #legend.key = element_rect(fill = "white", color = NA),
#     panel.border = element_rect(colour = "black", fill=NA))+
#    #legend.box.background = element_rect(colour = "red"))+
#   scale_colour_discrete(labels=c('Milk grain', 'Milk liquid', 'Water grain',"Water liquid"), name = "Kefir type")+
#   scale_shape_discrete(name="Data source")+
#    guides(colour = guide_legend(override.aes = list(size=10)),
#           shape = guide_legend(override.aes = list(size=10)))


p_total <- 
  
  ggplot(pcoa, aes(PC1, y=PC2,colour=data_source_specific))+
  geom_point(size=5, alpha=.25) +
  geom_point(data=centroid,size=10,shape=21,  alpha = 0.5, color="black",aes(fill=data_source_specific),  show.legend=FALSE)+ #aes(, colour=data_source)),size=10)+
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
  theme(legend.position = c(.76,.2),#legend.position = "",#,#axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        #legend.direction = "horizontal",
        #legend.box = "horizontal",
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text = element_blank(),
        axis.title.x = element_text(size=20), # remove x-axis labels
        axis.title.y = element_text(size=20), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.text=element_text(size = 20),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=20), #change legend title font size
        #legend.key = element_rect(fill = "white", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 2.5))+
  #legend.box.background = element_rect(colour = "red"))+
  #scale_colour_discrete()+
  scale_shape_discrete(name="Data source")+
  guides(colour = guide_legend(override.aes = list(size=10),title="Data source"))


library(cowplot)
library(dplyr)
library(ggplot2)


# Add density curves to y and x axis
xdens <- 
  axis_canvas(p_total, axis = "x") + 
  geom_density(data = pcoa, aes(x = PC1, fill = data_source_specific, colour = data_source_specific), alpha = 0.3)
ydens <-
  axis_canvas(p_total, axis = "y", coord_flip = TRUE) + 
  geom_density(data = pcoa, aes(x = PC2, fill = data_source_specific, colour = data_source_specific), alpha = 0.3) +
  coord_flip()

p_total_v2 <- 
p_total %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()





str(data.bray)
data.bray_metadata=data.frame(Row.names=rownames(as.matrix(data.bray)))

data.bray_metadata=
  merge(data.bray_metadata, pcoa, by="Row.names",all.x=TRUE, sort = F)

t <- 
  adonis2(data.bray~data.bray_metadata[,"data_source_specific"])


library(pairwiseAdonis)
#default is 999 permutations
res<-pairwiseAdonis::pairwise.adonis(data.bray,data.bray_metadata$data_source_specific)


res<-pairwiseAdonis::pairwise.adonis(data.bray~data.bray_metadata[,"data_source_specific"])





#subset_data1 <- as.matrix(data.bray)

# calculates the beta-dispersion for each group, when comparing 2 or more
dispersion <- betadisper(data.bray, group =data.bray_metadata$data_source_specific)


permutest(dispersion)

# tests if centroid distances are significantly different from each other
pathotype.disp.anova <- anova(dispersion) 



# test significance between each group
pathotype.disp.TukeyHSD <- TukeyHSD(dispersion)


pathotype.anosim <- anosim(data.bray, group =data.bray_metadata$data_source_specific)

#here 
units <- data.frame(units=as.character(c("Mortensen et al - Water kefir", "This study - Milk kefir",       "This study - Water kefir",     "Walsh et al - Milk kefir" )))

units <- as.data.frame(units)
subset_data1 <- as.matrix(data.bray)
aov_tab <- c()

for (i in c("Mortensen et al - Water kefir", "This study - Milk kefir",       "This study - Water kefir",     "Walsh et al - Milk kefir" )){
  
  
  for (j in units$units){
    
    
    if(i==j){
      next
    }else{
      subset_data1 <- as.matrix(data.bray)
      subset_data1 <- subset_data1[which(rownames(subset_data1) %in% data.bray_metadata$Row.names[which(data.bray_metadata$data_source_specific %in% c(i, j))]),
                                   
                                   which(colnames(subset_data1) %in% data.bray_metadata$Row.names[which(data.bray_metadata$data_source_specific %in% c(i,j))])
                                   
      ]
      
      
      subset_metadata <-
        data.bray_metadata[which(data.bray_metadata$data_source_specific%in% c(i, j)),]
      
      t <- 
        anosim(subset_data1, group =subset_metadata$data_source_specific)
      
      aov_tab <- rbind(aov_tab,data.frame(id=as.character(paste(i, "vs",j,sep=" ")),
                                          p_value=as.numeric(t$signif),
                                          R2=as.numeric(t$statistic)))
      
      
    }
    
    
  }
  
  units <- units[-c(which(units$units==j))]
}




###########################################################################################################################################################


# Get dominating species from each sample from both global and cs projects

###########################################################################################################################################################



dominanting_species=data.frame(sample_id=as.character(),
                               species=as.character(),
                               relative_abundance=as.numeric())
                              
metacache_long <- 
metacache %>% pivot_longer(!sample_id,names_to = "species",values_to = "relative_abundance")
  
  

for (sample_name in metacache$sample_id){
  
  
  dominanting_species <- 
    
    metacache_long %>% 
    filter(sample_id== sample_name) %>% 
    filter(relative_abundance==max(relative_abundance)) %>% 
    rbind(.,dominanting_species)
     
                              
}

# breakdown of the dominant species in each study
dominanting_species <- 
merge(dominanting_species,total_metadata,by.x="sample_id",by.y="merge_column",all.x=TRUE)
dominanting_species$group <- "milk"
  
dominanting_species$group[which(dominanting_species$`kefir type` %in% c("WL","WG"))] <- "water"

dominanting_species %>% 
dplyr::select(.,species,`kefir type`,data_source) %>% 
  table()



species_of_interest <- 
c(

  levels(as.factor(
dominanting_species$species[ c(which(dominanting_species$group=="water" &
                                      dominanting_species$species %in% names(which(table(dplyr::select(  dominanting_species[which(dominanting_species$group=="water"),], species)) >=  length(levels(as.factor( dominanting_species$sample_id[which(dominanting_species$group=="water")])))*.1))))]

)),

levels(as.factor(
dominanting_species$species[ c(which(dominanting_species$group=="milk" &
                                        dominanting_species$species %in% names(which(table(dplyr::select(  dominanting_species[which(dominanting_species$group=="milk"),], species))>=  length(levels(as.factor( dominanting_species$sample_id[which(dominanting_species$group=="milk")])))*.1))))]



)))





dominanting_species_v2<- dominanting_species[which(dominanting_species$species %in% species_of_interest),]


pcoa_v1 <- 
  merge(pcoa, dominanting_species,by.x="Row.names",by.y="sample_id",all.x=TRUE)




pcoa_v2 <- 
merge(pcoa, metacache_long[which(metacache_long$species %in% species_of_interest),] ,by.x="Row.names",by.y="sample_id",all.y=TRUE)





 
 pcoa_v1$species[ c(which(pcoa_v1$group=="water" &
                  pcoa_v1$species %in% names(which(table(dplyr::select(  dominanting_species, species))<=  length(levels(as.factor( pcoa$Row.names[which(pcoa_v1$group=="water")])))*.1))))] <- "Other"
# 
# 
# 
pcoa_v1$species[c(which(pcoa_v1$group=="milk" &
                         pcoa_v1$species %in% names(which(table(dplyr::select(  dominanting_species, species))<=  length(levels(as.factor( pcoa_v1$Row.names[which(pcoa_v1$group=="milk")])))*.1))))] <- "Other"
 




ggplot(pcoa_v1,aes(PC1, y=PC2,colour=species))+
  geom_point(size=5) +
  geom_point(data=centroid,size=10,shape=21, color="black",aes(fill=data_source_specific),  show.legend=FALSE)+ #aes(, colour=data_source)),size=10)+
  #geom_convexhull(alpha=.1)+
  stat_ellipse(geom = "polygon",
               aes(fill=species),
               alpha = 0.25,
               type = "norm",
               show.legend=FALSE)+
 # facet_wrap(~group)+

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
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=20), #change legend title font size
        #legend.key = element_rect(fill = "white", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 2.5))+
  #legend.box.background = element_rect(colour = "red"))+
  scale_colour_discrete(name="Data source")+
  scale_shape_discrete(name="Data source")+
  guides(colour = guide_legend(override.aes = list(size=10),ncol=2,byrow=TRUE))



# facet by species and RA
  
#   pcoa_v2 %>% 
#   
#   mutate(`kefir type`= gsub("MG|WG","Grain",
#                                    gsub("ML|WL","Liquid", `kefir type`)))%>% 
# ggplot(aes(PC1, y=PC2, label=relative_abundance,shape=`kefir type` ,colour=relative_abundance))+ #group=relative_abundance
#   geom_point(size=5)+
#   facet_wrap(~species)+
#   scale_color_distiller(palette = "RdYlBu")+
#   labs(colour="Relative abundance")+
#   
#   
#   ########
# #labs(x=paste("PCoA1 - ", round(percent_explained[1]), "%", sep=""), y=paste("PCoA2 - ", round(percent_explained[2]), "%", sep=""), title="") +
# #coord_equal() +
# theme_bw()+
#   theme(legend.position = c(.6,.2),#axis.text.x = element_blank(),  # remove x-axis text
#         legend.box = "horizontal",
#         legend.direction = "vertical",
#         axis.text = element_blank(), # remove y-axis text
#         axis.ticks = element_blank(),  # remove axis ticks
#         axis.title= element_blank(), # remove x-axis labels
#         panel.background = element_blank(), 
#         panel.grid.major = element_blank(),  #remove major-grid labels
#         panel.grid.minor = element_blank(),  #remove minor-grid labels
#         plot.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size = 2.5),
#         legend.text=element_text(face="plain", size = 15),
#         legend.key.size = unit(1.5, 'cm'), #change legend key size
#         legend.key.height = unit(1.5, 'cm'), #change legend key height
#         legend.key.width = unit(1.5, 'cm'), #change legend key width
#         legend.title = element_text(size=20), #change legend title font size
#         strip.background=element_rect(fill="white"),
#         strip.text = element_text(face="bold.italic", size=20))+
# guides(shape= guide_legend(override.aes = list(size=10)))
# 




library(viridis)
p_species <- 
pcoa_v2 %>% 
  mutate(`kefir type` = gsub("MG|WG", "Grain",
                             gsub("ML|WL", "Liquid", `kefir type`))) %>% 
  ggplot(aes(PC1, y = PC2, 
             label = relative_abundance, 
             shape = `kefir type`, 
             colour = relative_abundance)) +
  geom_point(size = 5) +
  facet_wrap(~species) +
  scale_color_viridis(option = "C", direction = -1) +  # Use a viridis palette with reversed colors
  labs(colour = "Relative Abundance") +
  theme_bw() +
  theme(legend.position = c(.6, .2),
        legend.box = "horizontal",
        legend.direction = "vertical",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2.5),
        legend.text = element_text(face = "plain", size = 15),
        legend.key.size = unit(1.5, 'cm'),
        legend.key.height = unit(1.5, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        legend.title = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold.italic", size = 20)) +
  guides(shape = guide_legend(override.aes = list(size = 10)))




# pcoa_v2 %>% 
#   mutate(`kefir type` = gsub("MG|WG", "Grain",
#                              gsub("ML|WL", "Liquid", `kefir type`))) %>% 
#   ggplot(aes(PC1, y = PC2, 
#              label = relative_abundance, 
#              shape = `kefir type`, 
#              colour = relative_abundance, 
#              alpha = relative_abundance)) + 
#   geom_point(size = 5) +
#   facet_wrap(~species) +
#   scale_color_distiller(palette = "RdYlBu") +
#   scale_alpha_continuous(range = c(0.3, 1), limits = c(20, 40)) +
#   labs(colour = "Relative Abundance", alpha = "Relative Abundance") +
#   theme_bw() +
#   theme(legend.position = c(.6, .2),
#         legend.box = "horizontal",
#         legend.direction = "vertical",
#         axis.text = element_blank(), 
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 2.5),
#         legend.text = element_text(face = "plain", size = 15),
#         legend.key.size = unit(1.5, 'cm'),
#         legend.key.height = unit(1.5, 'cm'),
#         legend.key.width = unit(1.5, 'cm'),
#         legend.title = element_text(size = 20),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(face = "bold.italic", size = 20)) +
#   guides(shape = "none")  # Hide shape legend
#  # guides(shape = guide_legend(override.aes = list(size = 8)))





library(ggpubr)
jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 5_model_v2.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)


ggarrange(p_total_v2,
          p_species,
          nrow=1,ncol=2,labels=c("A.","B."),  font.label = list(size = 20))
 graphics.off()

 
 
 
 library(Cairo)   # For high-quality PNG/PDF/SVG
 
 # Define figure size and resolution
 fig_width <- 10  # Inches (adjust based on journal requirements)
 fig_height <- 6  # Inches
 dpi_res <- 300   # High resolution (300 dpi for journal quality)
 
 
 
 combined_plot <- ggarrange(
   p_total_v2, p_species, 
   nrow = 1, ncol = 2, 
   labels = c("A.", "B."), 
   font.label = list(size = 20)
 )
 
 
  pdf("Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure5_model.pdf")#, width = fig_width, height = fig_height)
  print(combined_plot)
  dev.off()

 svg("Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure5_model.svg", width = fig_width, height = fig_height)
 print(combined_plot)
 dev.off()
 
 
 CairoPNG("Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure5_model.png", width = fig_width * dpi_res, height = fig_height * dpi_res, res = dpi_res)
 print(combined_plot)
 dev.off()
 
 
 
#write.csv(dominanting_species,"Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_dominating_species.csv",row.names = FALSE,quote = FALSE )

###########################################################################################################################################################

                                                                # Anaysis for figure 3 stops here 
# PCOA dominatating species on milk and water kefir seperately - old anlaysis used to make an old figure keeping as a reference 

###########################################################################################################################################################



metacache <- 
  metacache %>% remove_rownames() %>% column_to_rownames("sample_id")

beta_diversity<- c()




beta_diversity[["milk"]] <- vegdist(metacache[which(rownames(metacache) %in% total_metadata$merge_column[which(total_metadata$`kefir type` %in% c("MG","ML"))]),],
                                    na.rm = TRUE) #Bray-Curtis distance


beta_diversity[["water"]] <- vegdist( metacache[which(rownames(metacache) %in% total_metadata$merge_column[which(total_metadata$`kefir type` %in%c("WG","WL") )]),],
                                      na.rm = TRUE) #Bray-Curtis distance


data.b.pcoa=c()
pcoa <- c()

pcoa_plot <- c()
percent_explained <- 100* data.b.pcoa$eig/sum(data.b.pcoa$eig)

type <- "milk"
for (type in names(beta_diversity)){
  
  

  data.b.pcoa=cmdscale(beta_diversity[[type]],k=2,eig=TRUE,add = TRUE) #ordination
  
  pcoa[[type]] <- 
    data.frame(PC1 = data.b.pcoa$points[,1], PC2 = data.b.pcoa$points[,2])
  
  percent_explained <- 100* data.b.pcoa$eig/sum(data.b.pcoa$eig)

  
  pcoa[[type]] <- merge(pcoa[[type]],total_metadata,by.x=0,by.y="merge_column",all.x=TRUE)


  pcoa[[type]] <- merge( pcoa[[type]],dominanting_species,by.x="Row.names",by.y="sample_id",all.x=TRUE)
  

  pcoa[[type]]$species[-c(which(pcoa[[type]]$species %in% names(which(table(dplyr::select(  dominanting_species[which(dominanting_species$group==type),], species))>=  length(levels(as.factor( pcoa[[type]]$Row.names)))*.1))))] <- "Other"
  
  
  
  pcoa[[type]]$kefir_dataset <- type
  
  pcoa_plot[[type]] <- 
    
    pcoa[[type]] %>% 
    mutate(`kefir type.x`=gsub(".L", "Liquid",
                               gsub(".G","Grain",`kefir type.x`))) %>% 
    ggplot(aes(PC1, y=PC2, label=relative_abundance ,fill=species,colour=species,group=species,shape=`kefir type.x`))+
    geom_point(size=5) +
    #geom_convexhull(alpha=.1)+
    #stat_ellipse(geom = "polygon",
    #alpha = 0.25,
    #type = "norm")+

  ########
  labs(x=paste("PCoA1 - ", round(percent_explained[1]), "%", sep=""), y=paste("PCoA2 - ", round(percent_explained[2]), "%", sep=""), title="") +
    coord_equal() +
    theme_bw()+
    theme(legend.position = "right",#axis.text.x = element_blank(),  # remove x-axis text
      #axis.text.y = element_blank(), # remove y-axis text
      axis.ticks = element_blank(),  # remove axis ticks
      axis.text = element_blank(),
      axis.title = element_text(size = 20),
      #axis.title.x = element_text(size=18), # remove x-axis labels
      #axis.title.y = element_text(size=18), # remove y-axis labels
      panel.background = element_blank(), 
      panel.grid.major = element_blank(),  #remove major-grid labels
      panel.grid.minor = element_blank(),  #remove minor-grid labels
      plot.background = element_blank(),
      legend.text=element_text(face="italic",size = 20),
      legend.key.size = unit(2, 'cm'), #change legend key size
      legend.key.height = unit(2, 'cm'), #change legend key height
      legend.key.width = unit(2, 'cm'), #change legend key width
      legend.title = element_text(size=20) #change legend title font size
    )+
    #legend.box.background = element_rect(colour = "red"))+
    scale_colour_discrete( name = "Species detected")+
    scale_shape_discrete( name = "Kefir type")+
    guides(colour = guide_legend(override.aes = list(size=10)),
           fill = FALSE,
           shape=guide_legend(override.aes = list(size=10)))
  
  
  
}




library(ggpubr)

jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/Evolution/Figure 2_v3.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)


ggarrange(p_total,
  ggarrange(  pcoa_plot[[1]],  pcoa_plot[[2]], nrow=1,ncol=2, labels=c( "B.","C."),common.legend = FALSE,  font.label = list(size = 30)),
  nrow=2,ncol=1,labels=c("A.","C."),widths=c(10,1),heights = c(1,.5),  font.label = list(size = 30))
graphics.off()

