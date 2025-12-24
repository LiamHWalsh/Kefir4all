########################################################################################################################
#Research objectives in this script
#Describe the compositional profile of the CS dataset only
# Identify prevalent species
#identify dominating species
# track persistence of microbes in the grain that may be environmentally sourced
# Identify pathogenic microbes
########################################################################################################################


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


########################################################################################################################
#Import metacache data
########################################################################################################################
metacache <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_species_profile_v2.csv")
metacache_strain <-  read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_strain_profile_v2.csv")

metacache$sample_id <- gsub("-","_",metacache$sample_id)

metacache$sample_id <- gsub(".*/","",metacache$sample_id)


#metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]

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

print(paste("The citizen science project is composed of ",nrow(kefir4all_metadata))) 

print(paste("number of sample according to types in this study"))

print(table(kefir4all_metadata$`kefir type`))
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
#Modufy dataset into metacache citizen scientist and baseline used
########################################################################################################################


metacache [is.na(metacache )] <- 0
kefir4all_metadata <- 
kefir4all_metadata[-c(which(
kefir4all_metadata$`kefir type`=="Medium control")),]

metacache <- metacache[which(metacache$sample_id %in% 
kefir4all_metadata$merge_column), ]

#n1_species <-which(maxab_species < 0.1)


metacache <- 
metacache %>% column_to_rownames("sample_id")



print(paste("Number of samples in metacache is",length(rownames(metacache))))

# samples  not accounted for by metacache- base_line typically did not have a compostional profile 



kefir4all_metadata$merge_column[-c(which(kefir4all_metadata$merge_column %in% rownames(metacache)))]


###################
#Prevalence info
specie <- c()
sum_specie <- c()
species_prevalence <- c()


species_profile <- c()
species_profile[["water"]] <- metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL"))]),]#& 
                                                                                                               #kefir4all_metadata$Stage!="T0")]),]

species_profile[["milk"]] <- metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML"))]),] 
          
                                                                                                                                                                                                                

 #species_profile[["water"]] <- metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL")& 
                                                                                                              #kefir4all_metadata$Stage!="T0")]),]
                                                                                                                                                                                                                          
                                                                                                                                                                                                                           
#species_profile[["water_baseline"]] <- 
#metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                                                                                                                                                                       
        
#species_profile[["milk_baseline"]] <- 
  #metacache[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                   #kefir4all_metadata$Stage=="T0")]),]                             
  
  
set_entrez_key("34bbe7fc331d967e9da7f325f4a8b9e17509")

taxize_class <- c()

#taxize_name <- classification(tax_info$name, db = "ncbi")
taxize_tree <- c()

tax_phylo<- c()
pa <- c()
metadata <- c()
type <-  names(species_profile)[4]
maxab_species <- c()
total_compositional_data <- c()
type <- "water"

for (type in names(species_profile)){
  
  maxab_species <- apply(species_profile[[type]],2, max, na.rm=TRUE)
  
  n1_species <-names(which(maxab_species > .1))
  
  species_profile[[type]] <- species_profile[[type]][,which(colnames(species_profile[[type]]) %in%  n1_species )]
  
  species_prevalence[[type]] <- data.frame(specie=as.character(),
                                       sum_specie=as.numeric())
  
  #i <- "Lentilactobacillus hilgardii" 
for (i in colnames(species_profile[[type]])){
  specie <- i
  
  sum_specie <-  length(which(species_profile[[type]][,i] >.1))
  
  data.prevalence <- cbind(specie,sum_specie)
  species_prevalence[[type]] <- merge(data.prevalence,species_prevalence[[type]], all=TRUE)
  
}
  species_prevalence[[type]] <- species_prevalence[[type]][which(as.numeric(species_prevalence[[type]]$sum_specie)> nrow(species_profile[[type]])*.1),]
  species_prevalence[[type]]$kefir_type<- type
  species_profile[[type]] <- species_profile[[type]][,which(colnames(species_profile[[type]]) %in%    species_prevalence[[type]]$specie)]
  
  
  
   #taxize_class[[type]] <- taxize::classification(colnames(species_profile[[type]]), db = "ncbi")
  
  taxize_class[[type]] <- taxize::classification(species_prevalence[[type]]$specie, db = "ncbi")
  
  #taxize_name <- classification(tax_info$name, db = "ncbi")
  taxize_tree[[type]] <- taxize::class2tree(taxize_class[[type]], check = TRUE)
  
  tax_phylo[[type]] <- taxize_tree[[type]]$phylo
  
  pa[[type]] <- ggtree(tax_phylo[[type]], layout = "circular")#
  
metadata[[type]] <- data.frame(species=as.character(pa[[type]]$data$label),
                               genus=as.character(NA),
                               order=as.character(NA))


for (taxa in names(taxize_class[[type]])){
  metadata[[type]]$genus[which(metadata[[type]]$species==taxa)] <- 
    
    taxize_class[[type]][[taxa]]$name[which(taxize_class[[type]][[taxa]]$rank=="genus")]
  
  metadata[[type]]$order[which(metadata[[type]]$species==taxa)] <- 
    
    taxize_class[[type]][[taxa]]$name[which(taxize_class[[type]][[taxa]]$rank=="order")]
  
  
}

metadata[[type]] <- 
metadata[[type]][-c(which(is.na(metadata[[type]]$genus))),]

pa[[type]] <- pa[[type]] %<+% metadata[[type]]+
  # 
  geom_tippoint(mapping=aes(color=order), 
                size=5,
                show.legend=TRUE)+
  geom_tiplab(cex=.000010,offset=4,hjust=.05,linetype="dotted",size=4.5,face="Italic")+

  guides(colour=guide_legend(override.aes = list(size=10),title="Orders detected"))




species_profile[[type]] <- 
species_profile[[type]] %>% rownames_to_column("sample") %>% pivot_longer(!sample, names_to = "clade_name",values_to = "relative_abundance")




############################################################################################################################


############################################################################################################################
  pacman::p_load(ggnewscale,ggtreeExtra)
  pa[[type]]<- pa[[type]] + new_scale_fill()+
    geom_fruit(data=species_profile[[type]], geom=geom_tile,
               mapping=aes(x=sample,y=clade_name, fill=relative_abundance),
               offset = .7,size = 10,pwidth = 1,
               axis.params=list(axis="x",color = "white"))+
    #scale_fill_gradientn(colours = pal,breaks=col_breaks, name="Relative abundance" , labels=col_breaks_labels)+
    scale_fill_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
                         breaks= c(-5,0,10,30,70,90,100),
                         values= rescale(as.numeric(c(-5,0,10,30,70,90, 100))),
                         guide="colorbar",
                         name="Relative abundance" ,
                         labels=c(-5,0,10,30,70,90,100))+
    theme_void() +
    theme( legend.title = element_text(face="bold",size=25),
           legend.text = element_text(face="italic",size=18),
           legend.position="right",
           legend.key.size = unit(2, 'cm'))
  
  
  species_profile[[type]]$type <- type
total_compositional_data <- rbind(species_profile[[type]],total_compositional_data )

print(paste(type, "has", length(levels(as.factor( species_profile[[type]]$sample))),"number of sample"))



print(paste(type, "has", length(levels(as.factor( species_profile[[type]]$clade_name))),"prevalent microbes"))
  #.02
  
}

#write.csv(species_prevalence[["milk"]],"Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/milk_metacache_prevalence.csv" )
#write.csv(species_prevalence[["water"]],"Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/water_metacache_prevalence.csv" )

library(heatmaply)

#ggheatmap(
#metacache[,which(colnames(metacache) %in% levels(as.factor(total_compositional_data$clade_name)))]) 


total_compositional_data <- metacache[,which(colnames(metacache) %in% levels(as.factor(total_compositional_data$clade_name)))] %>%
  rownames_to_column("sample") %>% pivot_longer(!sample  , names_to  ="clade_name" , values_to ="relative_abundance")
                                     



kefir4all_metadata <- 
  merge(
    kefir4all_metadata, 
    rbind(dplyr::select(myfiles[["mk"]], merge_column,Sample, observations, category_confirmed),
          dplyr::select(myfiles[["wk"]], merge_column,Sample, observations, category_confirmed)),
    by="merge_column",
    all.x=TRUE)




total_compositional_data <-  
  merge(total_compositional_data,kefir4all_metadata,by.x="sample",by.y="merge_column")



total_compositional_data$type <- NA
#total_compositional_data$type <- gsub("milk","Milk kefir",total_compositional_data$type)
#total_compositional_data$type <- gsub("water","Water kefir",total_compositional_data$type)


#table(dplyr::select(total_compositional_data,clade_name,type))

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                                                 kefir4all_metadata$Stage!="T0")])] <- "Water kefir"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                                                 kefir4all_metadata$Stage!="T0")])] <- "Milk kefir"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
                                                                                                     kefir4all_metadata$Stage=="T0")])] <- "Water\n kefir\n T0"

total_compositional_data$type[which(total_compositional_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
                                                                                                                 kefir4all_metadata$Stage=="T0")])] <- "Milk\n kefir\n T0"


# remove extraction controls
#total_compositional_data <- total_compositional_data[-c(which(is.na(
    #total_compositional_data$type))),]

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




#total_compositional_data$`kefir type`[which(total_compositional_data$`kefir type` %in% c("MG","WG"))] <- "Grain"

#total_compositional_data$`kefir type`[which(total_compositional_data$`kefir type` %in% c("ML","WL"))] <- "Liquid"



#which(is.na(total_compositional_data$kefirdataset))
library(ggh4x)

p_heatmap <- c()
i <- "Milk kefir"



View(total_compositional_data[which(total_compositional_data$type=="Water\n kefir\n T0"),])

total_compositional_data <- total_compositional_data[-c(which(total_compositional_data$Sample.x %in% c("EC",   "Medium control") )),]



names(
species_prevalence)[1] <- "Water kefir"


names(
  species_prevalence)[2] <-"Milk kefir"




for (i in levels(as.factor((total_compositional_data$kefirdataset)))){
  
    p_heatmap[[i]] <-   total_compositional_data[which(total_compositional_data$kefirdataset==i &
                                                         total_compositional_data$clade_name  %in% species_prevalence[[i]]$specie
                                                         ),]%>%

      mutate(`kefir type`=gsub("ML|WL","Liquid",
                              gsub("MG|WG", "Grain",`kefir type` ))) %>% 
      mutate(type = factor(type, levels = rev(levels(factor(type))))) %>%  # Reverse the order of the levels
      mutate(sample_combined = paste(`kefir type`, sample, sep = "_")) %>%  # Combine kefir type and sample
      mutate(sample_combined = factor(sample_combined, levels = unique(sample_combined[order(`kefir type`, sample)]))) %>%  # Reorder the x-axis
       
      ggplot( aes(x = sample_combined  , y =clade_name)) +
       geom_tile(aes(fill = relative_abundance),#color = "white",
                 #lwd = 1.5,
                 linetype = 1) +
       #coord_fixed()+
       # scale_fill_gradientn(
       #colors =brewer.pal(3,"Set1")[2:1],breaks=c(0,1), labels=c("Absence","Presence"))+
       #guides(fill=guide_legend(title="Presence/\nAbsence"))+
       labs(x="", y="", title="")+ #y="Feature"
       #facet_wrap(~type,scales="free_x")+
       facet_grid(~type,margins=FALSE, space='free',scales="free")+
      # facet_grid2(.~CHR, scales = "free_x", space = "free_x", switch = "x",
      #             strip = strip_themed(
      #               background_x = elem_list_rect(
      #                 fill = rainbow(length(unique(df_temp$CHR)))))) +
       scale_fill_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
                            breaks= c(-5,0,10,30,70,90,100),
                            values= rescale(as.numeric(c(-5,0,10,30,70,90, 100))),
                            guide="colorbar",
                            name="Relative abundance" ,
                            labels=c(-5,0,10,30,70,90,100))+
       #facet_grid(Category~ cluster_info,margins=FALSE)+
       theme_bw()
    
    
    p_heatmap[[i]] <-     p_heatmap[[i]] +geom_point(data =   total_compositional_data[which(total_compositional_data$kefirdataset==i),]%>% 
                                                     mutate(`kefir type`=gsub("ML|WL","Liquid",
                                                                              gsub("MG|WG", "Grain",`kefir type` ))) %>% 
                                                       mutate(type = factor(type, levels = rev(levels(factor(type))))) %>%  # Reverse the order of the levels
                                                       mutate(sample_combined = paste(`kefir type`, sample, sep = "_")) %>%  # Combine kefir type and sample
                                                       mutate(sample_combined = factor(sample_combined, levels = unique(sample_combined[order(`kefir type`, sample)]))),  # Reorder the x-axis
                                                     aes(x=sample_combined,y=.35, col = `kefir type`), size = 4.6, shape = 15) +
                                         # geom_point(data =  total_compositional_data[which(total_compositional_data$kefirdataset==i),], 
      theme(legend.position = "top",
            strip.background = element_rect(
              color="black", fill="white"),
            strip.text = element_text(size=15.5),
            #strip.text.y = element_text(color = "white")
            legend.key.size = unit(2, 'cm'), #change legend key size
            legend.key.height = unit(2, 'cm'), #change legend key height
            legend.key.width = unit(2, 'cm'),
            legend.title = element_text( size=17.5, face="bold"),
            axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1, size = 10, face="italic"),
            legend.text = element_text(size = 17.5),
            axis.text.y = element_text(size=13.5, face="italic"),
            axis.title = element_text(face = "bold", size = 10.5)
      )+
      guides(color = guide_legend(override.aes = list(size=10),
                                  name="Relative abundance (%)"))
   
      
}


 #graphics.off()
 

 
 ###########################################################################################################################################################
 
 
 # Identify unique prevalent microbes in either milk or water kefir
 
 ###########################################################################################################################################################
 
  # 
   #species <-  "Leclercia adecarboxylata"  
   #species <- "Acetobacter aceti"
   
   #Look for microbes that are not present during the project but are present at T0
   
   
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
sample_name <- "TG_ID012_MG_T3_S189"




for (sample_name in levels(as.factor(total_compositional_data$sample))){
  if(sample_name %in% c("TG_ID012_MG_T3_S189", "TG_ID115_MG_T4_S96")){
    next
  }
  dominanting_species <- rbind(dominanting_species,
                               
  total_compositional_data %>% 
    filter(sample== sample_name) %>% 
    filter(relative_abundance==max(relative_abundance))
  
  )
}


dominanting_species[which(duplicated(dominanting_species$sample)),]

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


dominate_breakdown$MG_percent <- (dominate_breakdown$MG / samples_mk)*100

dominate_breakdown$ML_percent <- (dominate_breakdown$ML / samples_ml)*100
dominate_breakdown$WG_percent <- (dominate_breakdown$WG / samples_wg)*100
dominate_breakdown$WL_percent <- (dominate_breakdown$WL / samples_wl)*100

length(rownames(metacache)[which(rownames(metacache) %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL"))])])


View(dominate_breakdown)


sum(dominate_breakdown$MG_percent)

# = ggplot(pcoa, aes(x=PC1, y=PC2)) + geom_point(size=3) + theme_bw()#


#pcoa$kefir_dataset <- "Milk.kefir"
#pcoa$kefir_dataset[which(pcoa$`kefir type.x` %in% c("WL","WG"))] <-"Water kefir"

#pcoa <- pcoa[-c(which(pcoa$`kefir type.x`=="Liquid control")),]


#total_prevalence <- rbind(species_prevalence[["milk"]],species_prevalence[["water"]])





########################################################################################################################
# get grain samples from metacache 
########################################################################################################################




grain_microbes <- c()


grain_microbes[["Milk.kefir"]] <- 
  metacache[
    which(rownames(metacache) %in%
            kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                    kefir4all_metadata$`kefir type` %in% c("MG"))]
    ),]



grain_microbes[["Water.kefir"]]<- 
  
  metacache[
    which(rownames(metacache) %in%
            kefir4all_metadata$merge_column[which(kefir4all_metadata$Stage=="T0"&
                                                    kefir4all_metadata$`kefir type` %in% c("WG"))]
    ),]


grain_microbes[["Water.kefir"]] <- 
  grain_microbes[["Water.kefir"]][rownames(grain_microbes[["Water.kefir"]])=="TG_FL_T2_WG_0AN_S328",]


for (i in names(grain_microbes)){
  
 unclassified_rate=as.numeric(100-
                                                    rowSums(grain_microbes[[i]]))
  
  maxab_species <- apply( grain_microbes[[i]],2, max, na.rm=TRUE)
  
  n1_species <-names(which(maxab_species > .1))
  
  rare_biosphere_rate= rowSums(grain_microbes[[i]] [,-c(which(colnames( grain_microbes[[i]] ) %in%  n1_species ))])
    
  grain_microbes[[i]] <-  grain_microbes[[i]] [,which(colnames( grain_microbes[[i]] ) %in%  n1_species )]
  
  grain_microbes[[i]] <- grain_microbes[[i]] %>% rownames_to_column("sample_id") %>% pivot_longer(!sample_id,names_to="species",values_to = "relative_abundance")
  grain_microbes[[i]]$type <- i
  
  
  grain_microbes[[i]]=rbind(grain_microbes[[i]],
  data.frame(sample_id=as.character(unique(grain_microbes[[i]]$sample_id)),
             species="Unclassified",
             relative_abundance=unclassified_rate,
             type=i),
  
  
  data.frame(sample_id=as.character(unique(grain_microbes[[i]]$sample_id)),
             species="<0.1%" ,
             relative_abundance= rare_biosphere_rate,
             type=i)
  
  
  )
 
}



########################################################################################################################
# Plot grain species data 
########################################################################################################################




bubble_plot <- 
  ggplot(   
    rbind(  grain_microbes[[1]],
            grain_microbes[[2]]), aes(x=sample_id, y=species, size=relative_abundance, color=rev(relative_abundance))) +
  geom_point(alpha=0.5)  +
  labs(x="",y="", colour="Species")+
  facet_wrap(~type,scales="free")+
  theme_bw() +
  scale_size_continuous( breaks= c(30,20,10,1,0),
                         range= c(0,30),
                         name="Relative abundance (%)")+

  scale_color_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ),
                        breaks= c(30,20,10,1,0,-5),
                        values= rescale(as.numeric(c(-5,0,10,35))),
                         guide="colorbar",
                        name="")+

           theme(
             axis.text.x = element_blank(),
             axis.text.y = element_text(size=13.5,face="italic"),
             strip.background = element_rect(
               color="black", fill="white"),
             strip.text = element_text(size=15.5),
             #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
             axis.title.y = element_text( size=35, face="bold",hjust = 0.5, vjust = 1.5),
             legend.key.size = unit(2, 'cm'), #change legend key size
             legend.key.height = unit(2, 'cm'), #change legend key height
             legend.key.width = unit(2, 'cm'), #change legend key width
             legend.direction = "horizontal", 
             legend.position = "top",
             legend.box = "horizontal",
            # legend.text = element_text(size=30),
             legend.title = element_text(size=17.5),
             legend.text=element_text(size=17.5))
         


########################################################################################################################
# Identify microbes present in kefir grans >.1 RA
########################################################################################################################

dplyr::select(grain_microbes[[1]],species, relative_abundance)
dplyr::select(grain_microbes[[2]],species, relative_abundance)

########################################################################################################################
# Identify prevalent micorbes not in grains >.1 RA
########################################################################################################################


species_prevalence$`Milk kefir`[-c(which(species_prevalence$`Milk kefir`$specie %in% grain_microbes[["Milk.kefir"]]$species )),]

species_prevalence$`Water kefir`[-c(which(species_prevalence$`Water kefir`$specie %in% grain_microbes[["Water.kefir"]]$species )),]


print(paste(nrow(species_prevalence$`Milk kefir`[-c(which(species_prevalence$`Milk kefir`$specie %in% grain_microbes[["Milk.kefir"]]$species )),]),"prevelant species not detected at T0"))

print(paste(nrow(species_prevalence$`Water kefir`[-c(which(species_prevalence$`Water kefir`$specie %in% grain_microbes[["Water.kefir"]]$species )),]),"prevelant species not detected at T0"))



########################################################################################################################
#add this prevalence info to the p_heatmap data 
########################################################################################################################

label=
p_heatmap$`Milk kefir`$data$clade_name


label[-c(which(label %in% grain_microbes[["Milk.kefir"]]$species ))]=paste("*",label[-c(which(label %in% grain_microbes[["Milk.kefir"]]$species ))],sep="")


p_heatmap$`Milk kefir` <- p_heatmap$`Milk kefir` + scale_y_discrete(labels = label)

label=
  p_heatmap$`Water kefir`$data$clade_name


label[-c(which(label %in% grain_microbes[["Water.kefir"]]$species ))]=paste("*",label[-c(which(label %in% grain_microbes[["Water.kefir"]]$species ))],sep="")


p_heatmap$`Water kefir`<- p_heatmap$`Water kefir` + scale_y_discrete(labels = label)






########################################################################################################################
#plot TO compositional profile 
########################################################################################################################

total_compositional_data[which(
total_compositional_data$Stage=="T0"&
  total_compositional_data$`kefir type`=="ML"),] %>% 
  filter(relative_abundance>0.1)

ggplot(   
  grain_microbes[[1]], aes(x=sample_id, y=species, size=relative_abundance, color=species)) +
  geom_point(alpha=0.5)  +
 # scale_size(range=c(0, 55), name='Proportion (%)')+
  labs(x="",y="Species", colour="Species")+
  facet_wrap(~type,scales="free_x")+
  theme_bw() +
  theme(#plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5),size = 55), #element_text(color="red", size=14, face="bold",hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=30,face="italic"),
    strip.background=element_rect(fill="lightblue"),
    strip.text = element_text(size= 35 ),
    #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
    axis.title.y = element_text( size=35, face="bold",hjust = 0.5, vjust = 1.5),
    #legend.key.size = unit(3.5, 'cm'), #change legend key size
    #legend.key.height = unit(3.5, 'cm'), #change legend key height
    #legend.key.width = unit(3.5, 'cm'), #change legend key width
    legend.text = element_text(size=30,face = "italic"),
    legend.title = element_text(size=35))+
  scale_size_continuous( breaks= c(35,20,10,1,0),
                         range= c(0,95),
                         name="Proportion (%)")+

  # aspect.ratio = 2/1)+
  guides(color = guide_legend(override.aes = list(size =20)))
         





# 
# 
# ###########################################################################################################################################################
# 
# 
# # Sources of environmental microbes occurring at >.01 relative abundance
# 
# ###########################################################################################################################################################
# 
# maxab_species <- apply(metacache,2, max, na.rm=TRUE)
# 
# # all species at >.1 % RA not jst prevalent micorbes
# species_data  <- metacache [,-c(which(maxab_species < .1))] %>% rownames_to_column("sample") %>% pivot_longer(!sample,names_to = "clade_name",values_to = "relative_abundance")
# 
# 
# species_data <-  
#   merge(species_data,kefir4all_metadata,by.x="sample",by.y="merge_column")
# 
# 
# species_data$type <- NA
# #species_data$type <- gsub("milk","Milk kefir",species_data$type)
# #species_data$type <- gsub("water","Water kefir",species_data$type)
# 
# 
# #table(dplyr::select(species_data,clade_name,type))
# 
# species_data$type[which(species_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
#                                                                                                                  kefir4all_metadata$Stage!="T0")])] <- "Water.kefir"
# 
# species_data$type[which(species_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
#                                                                                                                  kefir4all_metadata$Stage!="T0")])] <- "Milk.kefir"
# 
# species_data$type[which(species_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("WG","WL") & 
#                                                                                                                  kefir4all_metadata$Stage=="T0")])] <- "Water.kefir.T0"
# 
# species_data$type[which(species_data$sample %in% kefir4all_metadata$merge_column[which(kefir4all_metadata$`kefir type` %in% c("MG","ML") & 
#                                                                                                                  kefir4all_metadata$Stage=="T0")])] <- "Milk.kefir.T0"
# 
# 
# # remove extraction controls
# #species_data <- species_data[-c(which(is.na(
#   #species_data$type))),]
# 
# 
# compositional_differences_.1 <- data.frame(species=as.character( levels(as.factor(species_data$clade_name))),
#                                         `Milk kefir`=as.character(NA),
#                                         `Milk kefir T0`=as.character(NA),
#                                         `Water kefir`=as.character(NA),
#                                         `Water kefir T0`=as.character(NA)
# )
# 
# species <- "Leuconostoc mesenteroides"
# type_2 <- "Milk.kefir"
# for (species in compositional_differences_.1$species){
#   
#   for (type_2 in levels(as.factor(species_data$type))){
#     compositional_differences_.1[which(   compositional_differences_.1$species==species),which(colnames(   compositional_differences_.1 )==type_2)] <-
#       nrow(subset(species_data, subset=clade_name==species &
#                     type==type_2 &
#                     relative_abundance >0.1))
#   }
# }
# 
# species_of_interest <- 
#   subset( compositional_differences_.1, subset=
#             compositional_differences_.1$Milk.kefir.T0==0 &
#             compositional_differences_.1$Water.kefir.T0==0 &
#             compositional_differences_.1$Milk.kefir>0 &
#             compositional_differences_.1$Water.kefir>0
#   ) %>% dplyr::select(species) %>% 
#   mutate(type="Occurs in both")
#   
# 
#    # Identify potential environmental microbes in both milk and water kefir 
# 
#    
#    # Identify potential environmental microbes in milk kefir
#    
# t<- 
#      compositional_differences_.1[which(compositional_differences_.1$Milk.kefir.T0==0 & compositional_differences_.1$Milk.kefir>0) ,] %>% 
#        dplyr::select(species) %>% mutate(type="Milk.kefir")
#      
#    t <- t[-c(which(t$species %in% species_of_interest$species)),]
#    
# 
#    # Identify potential environmental microbes in water kefir
#    
#   t1 <-
#      compositional_differences[which(compositional_differences$Water.kefir.T0==0 & compositional_differences$Water.kefir>0) ,] %>% 
#        dplyr::select(species)%>% mutate(type="Water.kefir")
#    
#    
#    t1 <- t1[-c(which(t1$species %in% species_of_interest$species)),]
#    
#    
#    species_of_interest <-  rbind(species_of_interest,t,t1)
#    
#    
# 
#    
# 
# print(paste("okay, we can see", nrow(species_of_interest), "microbes that may be environmentally derived. That's not the most interesting unless they persist. Lets look at that now"))
#    
# environmental_microbes  <- species_of_interest

#write.csv(environmental_microbes, "Q:/H2020 Master/Citizen Science Project/Results/00/environmental_microbes.txt",quote=FALSE,row.names = FALSE)

library(ggpubr)
# 


jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 2_v4.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)


ggarrange( bubble_plot,
  ggarrange(p_heatmap[[1]],p_heatmap[[2]], nrow=1,ncol=2, labels=c("B.","C."),common.legend = TRUE,  font.label = list(size = 30)),
            nrow=2,ncol=1, heights = c(5, 6),widths=c(10,1),labels=c("A.",""),  font.label = list(size = 30))
graphics.off()

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







