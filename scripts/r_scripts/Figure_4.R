#############################################################################################################################################################################################################################################

#Libraries used
#############################################################################################################################################################################################################################################


.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")

#update.packages(oldPkgs = "ggstatsplot")

pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis,tidytext,Hmisc)
pacman::p_load(readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales,ggrepel,ggside,ggnewscale,ggstatsplot,DESeq2,ltm,dplyr)
library(vegan)
library(ggplot2)
library(grid)

pacman::p_load(cooccur,visNetwork,igraph,indicspecies)


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



dominating_species <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/dominating species/04_metacache_dominating_species.csv")


dominating_species <- dominating_species[-c(which(duplicated(dominating_species$sample_id))),]

dominating_species <-
dominating_species[-c(which(dominating_species$sample_id== "TG_ID029_ML_T2_S16")),]

kefir4all_metadata <- kefir4all_metadata[-c(which(kefir4all_metadata$merge_column== "TG_ID029_ML_T2_S16")),]

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


# 
setwd("Q:/H2020 Master/Citizen Science Project/Results/07_metabolomics/")
# 
# 
# 
# temp <- c()
temp = list.files(pattern="*.xlsx", recursive = FALSE, full.names = TRUE)




myfiles = lapply(temp,read_excel)



names( myfiles) <- c("Milk.kefir","Water.kefir")


colnames(myfiles[["Water.kefir"]])[which(colnames(myfiles[["Water.kefir"]])=="Name")] <- "Compound"


colnames( myfiles[["Water.kefir"]])[
  c( grep("^[0-9][0-9]_T*", colnames( myfiles[["Water.kefir"]])))] <- paste("0",
                                                                            colnames( myfiles[["Water.kefir"]])[
                                                                              c( grep("^[0-9][0-9]_T*", colnames( myfiles[["Water.kefir"]])))],sep="")

nrow(kefir4all_metadata)
kefir4all_metadata <- merge(kefir4all_metadata,
                            dplyr::select(dominating_species,-c(data_source, `kefir type`)),
                            by.x="merge_column",
                            by.y="sample_id",
                            all.x=TRUE)
nrow(kefir4all_metadata)



nrow(kefir4all_metadata)
kefir4all_metadata <- merge(kefir4all_metadata,
                           Citizen_Scientist_metadata_v8,
                            by.x="Sample",
                            by.y="ID",
                            all.x=TRUE)
nrow(kefir4all_metadata)




kefir4all_metadata$category_confirmed[which(kefir4all_metadata$Stage=="T0")] <- 
           gsub("GO",	"Goat",
                gsub( "FU",	"Cow milk (whole)",
                      gsub("LO", 	"Cow milk (low fat)",
                           gsub("SK",	"Cow milk (skim)",
                                gsub("FG",	"White sugar-dried fig",
                                     gsub("FL" , "White sugar-dried fig and fresh lemon",
                                          gsub( "AP",	"White sugar-dried apricot",
                                                gsub( "BR",	"Brown sugar-None",
                                                      kefir4all_metadata$Sample[which(kefir4all_metadata$Stage=="T0")]     
                                                      
                                                      
                                                ))))))))



sample_type <- c()   
i="Water.kefir"
myfiles_wide <- c()
prevalent_class <- c()
metadata_metabolomics  <- c()

prevalent_compounds <- c()


core_compounds <- c()

t <- c()
t_names <- c()
means <- c()


compound_frequency <- c()


for (i in names(myfiles)){
  
  
  metadata <-kefir4all_metadata
  
  if(i=="Milk.kefir"){
    sample_type <- "ML"
  }else{
    sample_type <- "WL"
  }
  
  colnames( myfiles[[i]])[
    grep("_T*", colnames( myfiles[[i]]))] <- 
    colnames( myfiles[[i]])[
      grep("_T*", colnames( myfiles[[i]]))] %>% gsub("t","T",.) %>% gsub("__","_",.) %>% gsub("_0AET","_0AE",.) #%>% paste("ID",.,sep="")
  
  
  metadata <-  metadata[which( metadata$`kefir type` == sample_type ),]
  metadata$merge_column_vs <- paste(gsub("ID","",   metadata$Sample), metadata$Stage,sep="_")
  
  metadata$merge_column_vs[which( metadata$Stage=="T0" )] <-   gsub("_S[0-9].*|_.L|TG_","",metadata$merge_column[which( metadata$Stage=="T0" )]) #gsub("TG_|_S[0-9]+\\.?[0-9]+.*|_.L","",metadata$merge_column[which( metadata$Stage=="T0" )])
  metadata_metabolomics[[i]] <- metadata
  
  
  myfiles_wide[[i]] <-   myfiles[[i]]
  
  
  compound_frequency[[i]] <- 
  
  as.data.frame(table(myfiles_wide[[i]]$Class))
  
  
  
  #print(paste(i, " has ", length(myfiles_wide[[i]]$Compound), " volatiles", "This consisted of",  paste(  compound_frequency[[i]]$Var1,"(",   compound_frequency[[i]]$Freq, ")",collaspe="",sep=""), sep=""))
  
  print(paste(i, " has ", length(myfiles_wide[[i]]$Compound), " volatiles", " This consisted of ",   paste(compound_frequency[[i]]$Var1,"(",   compound_frequency[[i]]$Freq, ")",collapse=", "), sep=""))
  
 # what the average number of volatiles
 

  myfiles[[i]] <- 
    myfiles[[i]] %>% pivot_longer(!c(  Class,   Compound,       CAS,      RI, `Ref RI`, `Reference odour descriptor` ),names_to = "sample_id",values_to = "concentrations")
  
  t1 <- 
    nrow(myfiles[[i]])
  
  
  myfiles[[i]] <- 
    merge( myfiles[[i]], 
           metadata, 
           by.x="sample_id",
           by.y="merge_column_vs",
           all.x=TRUE)
  
  t2 <- 
    nrow(myfiles[[i]])
  
  if(t1!=t2){
    print(paste("merge command for",i," is not correct revise now",sep=""))
    #break
  }
  
  

  prevalent_class[[i]] <- 
  names(which(table(myfiles_wide[[i]]$Class) > round(length(levels(as.factor(myfiles[[i]]$sample_id)))*.10)))
  
  prevalent_compounds[[i]] <- 
  names(which(table(myfiles[[i]]$Compound[-c(
    which(myfiles[[i]]$concentrations==0.000))]) > round(length(levels(as.factor(myfiles[[i]]$sample_id)))*.10)))
  
t_names <- 
  
    names(which(table(myfiles[[i]]$Compound[-c(
      which(myfiles[[i]]$concentrations==0.000))]) >= length(levels(as.factor(myfiles[[i]]$sample_id)))))
  

  t <- 
data.frame(core_compounds=as.character(t_names),
           min=as.numeric(NA),
           max=as.numeric(NA),
           mean=as.numeric(NA))


for (j in t_names){
  
  
  t$min[which(t$core_compounds==j)] <- min(myfiles[[i]]$concentrations[which(myfiles[[i]]$Compound==j)])
  t$max[which(t$core_compounds==j)] <- max(myfiles[[i]]$concentrations[which(myfiles[[i]]$Compound==j)])
  t$mean[which(t$core_compounds==j)] <- mean(myfiles[[i]]$concentrations[which(myfiles[[i]]$Compound==j)])
}

  t <- 
  merge(t, 
  
  dplyr::select(myfiles_wide[[i]],Class,  Compound),
  by.x="core_compounds",
  by.y="Compound",
  all.x=TRUE)
  

  
  
t$rank <- 
  rank(-t$mean)

  core_compounds[[i]] <- t


  
  
  means[[i]] <- aggregate(concentrations ~  Compound , data = myfiles[[i]], mean)
  
  
  
  means[[i]] <- 
    
    merge(means[[i]], 
          
          dplyr::select(myfiles_wide[[i]],Class,  Compound),
          by="Compound",
          all.x=TRUE)
  

  
  means[[i]]$rank <- rank(-means[[i]]$concentrations)
  
 
  
  means[[i]]$core <- NA
  
  means[[i]]$core[which( means[[i]]$Compound %in% core_compounds[[i]]$core_compounds)] <- "*"
  
  
 print(paste(  core_compounds[[i]]$core_compounds," (",  core_compounds[[i]]$Class,")" ,collapse = ", ",sep=""))
  
  

}

########################################################################################################################
#Alpha diversity total dataset
########################################################################################################################



myfiles[["Milk.kefir"]]$group <- "Milk kefir"
myfiles[["Water.kefir"]]$group <- "Water kefir"



total_tile_plot_data <- 
rbind(
myfiles[[1]],

myfiles[[2]]) 

total_tile_plot_data$Class[-c(which(total_tile_plot_data$Class %in% levels(as.factor(c(prevalent_class[[1]],prevalent_class[[2]])))))] <- "Other"


levels(as.factor(total_tile_plot_data$sample_id[which(is.na(total_tile_plot_data$Stage))]))

total_tile_plot_data$Stage[which(is.na(total_tile_plot_data$Stage))] <- "T2"
total_tile_plot_data %>%  
 mutate(reference_odour=gsub(",.*","",`Reference odour descriptor`))

  table(total_tile_plot_data$reference_odour)
  
 
  pa <-  total_tile_plot_data %>% 
    mutate(Stage=
  
  gsub("T1","wk01",
    gsub("T2","wk05",
      gsub("T3","wk09",
        gsub("T4","wk13", 
          gsub("T5","wk17", 
           gsub("T6","wk21",Stage))))))) %>% 
  
  
    mutate(sample_stage_unit=paste(Stage,merge_column,sep="_")) %>% 
  #  filter(`kefir type` %in% c("MG","ML","WG","WL")) %>% 
  # arrange(Sample,Mechanism) %>% 
  ggplot(aes(y = sample_stage_unit,x =Compound)) +
  geom_tile(aes(fill = concentrations),#color = "white",
            #lwd = 1.5,
            linetype = 1) +
  labs(x="Samples", y="Volatiles", title="", fill="Abundance")+ #y="Feature"
  #geom_ysidetile(aes(x = "Fermentation conditions", yfill = conditions))+
    scale_y_reordered() +
  facet_wrap(~group ,scales="free",
             )+

scale_fill_gradientn(colours = c("Dark blue","yellow","yellow2","orange", "red","darkred" ))+
 # scale_fill_gradient(low = "green", high = "red") + 
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
    axis.text.x = element_text(size=7.5,angle = 90,hjust = .7),
    strip.background = element_rect(
      color="black", fill="white"),
    strip.text = element_text(size=25),
    #axis.title.x = element_text( size=35, face="bold",hjust = 0.5,vjust = -2),
    axis.title.y = element_text( size=10, face="bold",hjust = 0.5, vjust = 1.5),
    legend.key.size = unit(1.25, 'cm'), #change legend key size
    legend.key.height = unit(1.25, 'cm'), #change legend key height
    legend.key.width = unit(1.25, 'cm'), #change legend key width
    legend.text = element_text(size=15,angle=45, hjust=1),
    legend.title = element_text( size=25, face="bold",vjust = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    new_scale_fill()+
  geom_ysidetile(aes(x = "Fermentation conditions", yfill = Stage))+
      guides(yfill = guide_legend(nrow = 1,label.theme = element_text(angle = 0)))
    
  


  


########################################################################################################################
#Alpha diversity total dataset
########################################################################################################################

  my.files_summary_shannon <- c()
  my.files_richness <- c()
  my.files_evenness <- c()
  my.files_beta <- c()
  my.files_summary<-c()
    
  for (i in names(myfiles)){

my.files_summary_shannon[[i]] <- vegan::diversity(  t(dplyr::select(myfiles_wide[[i]], -c("Class", "CAS","RI","Ref RI","Reference odour descriptor"  )) %>% column_to_rownames("Compound")) )

my.files_richness[[i]] <- specnumber(  t(dplyr::select(myfiles_wide[[i]], -c("Class", "CAS","RI","Ref RI","Reference odour descriptor"  )) %>% column_to_rownames("Compound")))
my.files_evenness[[i]] <- my.files_summary_shannon[[i]]/log(my.files_richness[[i]])
my.files_beta[[i]] <- vegdist(  t(dplyr::select(myfiles_wide[[i]], -c("Class", "CAS","RI","Ref RI","Reference odour descriptor"  )) %>% column_to_rownames("Compound")), method = "bray")

my.files_summary[[i]]<- cbind(shannon = my.files_summary_shannon[[i]], richness = my.files_richness[[i]], pielou = my.files_evenness[[i]],site= rownames(t(dplyr::select(myfiles_wide[[i]], -c("Class", "CAS","RI","Ref RI","Reference odour descriptor"  )) %>% column_to_rownames("Compound"))))


my.files_summary[[i]]<- as.data.frame(my.files_summary[[i]])
my.files_summary[[i]] <- merge(my.files_summary[[i]], metadata_metabolomics[[i]],by.x="site",by.y="merge_column_vs",all.x=TRUE)



my.files_summary[[i]]$shannon <- as.numeric( my.files_summary[[i]]$shannon)
#setwd("Q:/H2020 Master/Citizen Science Project/Plots/Evolution")

#jpeg(filename='Alpha diversity_kefir_type.jpeg', width = 35*700, height=30*700,res=1700,pointsize = 15) #, width=2000, height=1950)


}




metabolites_plot_data <- rbind( my.files_summary[[1]],
                                my.files_summary[[2]] )

metabolites_plot_data$timeframe <- "Early stage"


metabolites_plot_data$timeframe[which(metabolites_plot_data$Stage %in% c("T4", "T5","T6"))] <- "Late stage"


  my.files_summary[["Water.kefir"]]$group <- "Water kefir"
  my.files_summary[["Milk.kefir"]]$group <- "Milk kefir"
  
  
  

  
  grouped_ggbetweenstats(
    data =   
      metabolites_plot_data[which(metabolites_plot_data$species %in% c("Lactobacillus helveticus","Lactococcus cremoris","Lactobacillus kefiranofaciens","Lactococcus lactis","Lacticaseibacillus paracasei","Zymomonas mobilis","Lentilactobacillus hilgardii" )),],
    x = species,
    y = shannon,
    fill=species,
    color=species,
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
  


  
  grouped_ggbetweenstats(
    data =   
      metabolites_plot_data,
    x = timeframe,
    y = shannon,
    fill=timeframe,
    color=timeframe,
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




  
  #############################################################################################################################
  
  
  #############################################################################################################################
  
  t <- c()
  my.files_beta <-c()
  my.files_beta_cmdscale <- c()
  pcoa <- c()
  en_pcoa <- c()
  envfit_output <-c()
  pcoa_data <- c()
  
  env_correlations <- c()
  base <- c()
  for (i in names(myfiles)){
    
    
    t <- t(myfiles_wide[[i]][,which(colnames(myfiles_wide[[i]] ) %in% c("Compound", metadata_metabolomics[[i]]$merge_column_vs ))] %>% column_to_rownames("Compound"))
  
    my.files_beta[[i]] <- vegdist(t)
  
  
   
    my.files_beta_cmdscale[[i]] <-  cmdscale(    my.files_beta[[i]],k=2,eig=TRUE,add = TRUE) #ordination
    
    
    
    pcoa[[i]] = data.frame(PC1 =  my.files_beta_cmdscale[[i]]$points[,1], PC2 =  my.files_beta_cmdscale[[i]]$points[,2])
    percent_explained <- 100*  my.files_beta_cmdscale[[i]]$eig/sum( my.files_beta_cmdscale[[i]]$eig)
    
    
    pcoa[[i]] <- merge(    pcoa[[i]],metadata_metabolomics[[i]],by.x=0,by.y="merge_column_vs",all.x=TRUE)
    
    t1= nrow(t)
    
    t <- merge(t,
               dplyr::select(metadata_metabolomics[[i]], Stage, species, category_confirmed,merge_column_vs ),
    by.x=0,
    by.y="merge_column_vs",
    all.x=TRUE)
    
    
    
   t2= nrow(t)
    if(t1!=t2){
      print(paste("merge command for",i," is not correct revise now",sep=""))
      #break
    }
    
    
    
    t <- t %>% column_to_rownames("Row.names")
    en_pcoa <- vegan::envfit(   my.files_beta_cmdscale[[i]],   t,  permutations = 1000,na.rm = TRUE)
    
    en_pcoa$vectors$pvals_fdr <-  p.adjust(en_pcoa$vectors$pvals,method = "bonferroni")
    
    
    envfit_output[[i]] <- en_pcoa
    
    
    
    
    # ggplot(   pcoa[[i]], aes(PC1, y=PC2,colour= Stage,fill= Stage))+
    #   geom_point(size=5) +
    # 
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
    #   theme(#legend.position = "none",#axis.text.x = element_blank(),  # remove x-axis text
    #     #axis.text.y = element_blank(), # remove y-axis text
    #     axis.ticks = element_blank(),  # remove axis ticks
    #     axis.text = element_blank(),
    #     axis.title.x = element_text(size=18), # remove x-axis labels
    #     axis.title.y = element_text(size=18), # remove y-axis labels
    #     panel.background = element_blank(), 
    #     panel.grid.major = element_blank(),  #remove major-grid labels
    #     panel.grid.minor = element_blank(),  #remove minor-grid labels
    #     plot.background = element_blank(),
    #     legend.text=element_text(face="italic",size = 15),
    #     legend.key.size = unit(2, 'cm'), #change legend key size
    #     legend.key.height = unit(2, 'cm'), #change legend key height
    #     legend.key.width = unit(2, 'cm'), #change legend key width
    #     legend.title = element_text(size=20) #change legend title font size
    #   )
    
    
    
    
    
    A <- as.list(en_pcoa$vectors)
    pvals<-as.data.frame(A$pvals)
    pvals_fdr<-as.data.frame(A$pvals_fdr)
    arrows<-as.data.frame(A$arrows*sqrt(A$r))
    C<- cbind(arrows, pvals,pvals_fdr, as.data.frame(A$r))
    
    Cred<-subset(C,pvals<0.05)
    Cred <- cbind(Cred, Species = rownames(Cred))
    
    colnames(Cred)[which(colnames(Cred)=="A$r")] <- "cor"
    
    
    
    

    
    for (target in c("Stage", "species", "category_confirmed")){
      
      data.scores_plot <-     
        pcoa[[i]] 
      
      
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
        centroid$PC1[which(centroid$target==e)] <- mean(data.scores_plot$PC1 [which(data.scores_plot$target ==e)])
        centroid$PC2[which(centroid$target==e)] <- mean(data.scores_plot$PC2[which(data.scores_plot$target==e)])
        
      }#end of e
      
      
      
      #test <- 
      #dplyr::select( data.scores, c( Row.names, V1, V2, levels(as.factor( mech_total_resistome$Mechanism)))) %>% pivot_longer(!c( Row.names, V1, V2),values_to = "RA", names_to = "mech") 
      
      
      library(ggnewscale)
      
      colnames(data.scores_plot)[which(colnames(data.scores_plot)=="target")] <- "Fermentation parameter"
      
      pcoa_data[[i]][[target]] <- 
        
        ggplot(data = data.scores_plot, aes(x = PC1, y = PC2)) + 
        geom_point(data = data.scores_plot, aes(colour = `Fermentation parameter`), size = 3, alpha = 0.5) + 
        geom_point(data=centroid,size=10,shape=21, color="black",aes(fill=target),  show.legend=FALSE)+ 
        stat_ellipse(geom = "polygon",
                     aes(fill= `Fermentation parameter`),
                     alpha = 0.25,
                     type = "norm",
                     show.legend=FALSE)+
        #scale_colour_manual(values = c("orange", "steelblue"))  + 
        #scale_fill_manual(values = c("orange", "steelblue"))  + 
        new_scale_color()+
        # geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
        #              data =  Cred, size =1, alpha = 0.5, colour = "grey30") +
        geom_point(data = Cred, aes(x = Dim1, y = Dim2, colour = cor), 
                   shape = "diamond", size = 4, alpha = 0.6) +
        geom_text_repel(data = Cred, aes(x = Dim1, y = Dim2, colour = cor), 
                        label = row.names(Cred), fontface = "bold",min.segment.length=0.01) + 
        #geom_text(data = Cred, aes(x = Dim1, y = Dim2), colour = "grey30", 
        #fontface = "bold", label = row.names(Cred)) + 
        theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
              plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
                                                          box.color = "#748696",size=35, lineheight = 2),
              panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
              axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
              legend.title = element_text(size = 15), 
              legend.text = element_text(size = 15),
              legend.position = c(.075,.8),
              legend.direction="vertical")+
        #legend.box="Horizontal") + 
        labs(colour = "Fermentation conditions",
             x="PC-1",
             y="PC-2", 
             title = gsub("\\."," ",i))+
        scale_color_distiller(palette = "RdYlBu")
      
      
      
      
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
    
    
    #paste(env_correlations[[i]]$Species[which(env_correlations[[i]]$rank<=5)],sep="",collapse = ", ")
    
    print(paste("in", tolower(gsub("\\."," ",i)), length(env_correlations[[i]]$Species), "statistically significant volatiles contributed to differences across", tolower(gsub("\\."," ",i)) ,
                "metagenomes, the most distinguishing of which included ",
                
                paste(tolower(env_correlations[[i]]$Species[which(env_correlations[[i]]$rank<=5)])," (", round(env_correlations[[i]]$cor[which(env_correlations[[i]]$rank<=5)],2), ")",collapse = ", ",sep="")
                
    ))
    
    
    base[[i]] <- 
      
      dplyr::select(myfiles[[i]],sample_id, category_confirmed)
    
    base[[i]] <- 
      base[[i]][-c(which(duplicated(base[[i]]$sample_id))),]
    
    

    
    

    
  }#end of i

  table(base[["Milk.kefir"]]$category_confirmed)
  table(base[["Water.kefir"]]$category_confirmed)


  
  
  ########################################################################################################################
  #Import metacache data
  ########################################################################################################################
  
  metacache <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_species_profile_v2.csv")
  metacache_strain <-  read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_strain_profile.csv")
  
  metacache$sample_id <- gsub("-","_",metacache$sample_id)
  #metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]
  
  
  metacache [is.na(metacache )] <- 0
  
  
  metacache$sample_id <- 
    gsub(".*/","",metacache$sample_id)
  
  

  #metacache$sample_id <- gsub("-","_",metacache$sample_id)
  metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]

  
  metacache <- metacache %>% remove_rownames() %>% column_to_rownames("sample_id")
  species_correlations <- c()
  species_profile <- c()
  
  metabolite_profile <-c() 
  
  merged_data <- c()
  cor_matrix <- c()
  correlations <- c()
  combo <- c()
  
  ## function to tidy the outputs

  
  tidy_corr <- function(x, value) {
    
    x %>%
      #cor_matrix[[i]]$r %>% 
     as.data.frame() %>% 
      rownames_to_column("species") %>%
      filter(species %in% colnames(species_profile[[i]])) %>%
      dplyr::select(species, all_of(colnames(metabolite_profile[[i]]))) %>%
      pivot_longer(!species, names_to = "functional_feature", values_to = value)
    
  }
  
  
  
  count_gt_one <- function(x) {
    sum(x > .1)
  }
  
  
  
  
  for (i in names(myfiles)){
    
  species_profile[[i]] <- metacache[which(rownames(metacache) %in%   levels(as.factor(myfiles[[i]]$merge_column))),]#& 
  
 
  

  
  species_profile[[i]] <-   na.omit(    species_profile[[i]] )
  
  
  

  
  
  maxab_species <- apply(species_profile[[i]],2, max, na.rm=TRUE)
  
  n1_species <-names(which(maxab_species >= .1))
  
  species_profile[[i]] <- species_profile[[i]][,which(colnames(species_profile[[i]]) %in%  n1_species )]
  
  
  maxab <- apply(species_profile[[i]],2,   count_gt_one)
  
  #n2_species <-names(which(maxab >= 6)) #round(length(levels(as.factor(rownames(species_profile[[i]]))))*.1)))
  
  
  #species_profile[[i]] <- species_profile[[i]][,which(colnames(species_profile[[i]]) %in%  n2_species )]
  
  
  if(i=="Water.kefir"){
  metabolite_profile[[i]] <- as.data.frame(dplyr::select(myfiles[[i]][-c(which(myfiles[[i]]$sample_id=="113_T2")),], merge_column, Compound, concentrations) %>% pivot_wider(names_from = Compound, values_from = concentrations)) %>% column_to_rownames("merge_column" )
  
  }else{
    
    metabolite_profile[[i]] <- as.data.frame(dplyr::select(myfiles[[i]], merge_column, Compound, concentrations) %>% pivot_wider(names_from = Compound, values_from = concentrations)) %>% column_to_rownames("merge_column" )
    
    
}
  metabolite_profile[[i]] <- metabolite_profile[[i]][match(rownames(  species_profile[[i]]),rownames(metabolite_profile[[i]])), ]
  
  

  if(
 identical(rownames(metabolite_profile[[i]]),  rownames(species_profile[[i]] ))!=TRUE){
    
    print(paste("species and metabolic dataframes do not align for",i,"crashing now"))
    break
  }
  
  
  #cor_matrix[[i]] <- cor(  species_profile[[i]],   metabolite_profile[[i]])
  
  
  cor_matrix[[i]] <- rcorr(as.matrix(species_profile[[i]]),  as.matrix( metabolite_profile[[i]]), "spearman")
  
  
  
  ## get the p-values and do p-value correction
  p_values <- cor_matrix[[i]]$P %>%
    tidy_corr(., "p") %>%
    mutate(fdr = p.adjust(p, method = "fdr"))
  
  ## get the r-values
  r_values <-
  cor_matrix[[i]]$r%>%
    tidy_corr(., "r")
  
  
  combo <- inner_join(r_values, p_values, by = c("species", "functional_feature")) %>%
    arrange(desc(r * r))
  
  combo$name <- gsub("\\."," ",i)


  
  correlations[[i]] <- combo[which(combo$p <= 0.05), ]
  #write.csv(  correlations_total[[i]], paste("Q:/H2020 Master/Citizen Science Project/Results/07_metabolomics/",i,"species_correlations_cs_data_set_0.1.csv"), row.names = FALSE , quote=FALSE)

 #  combo[which(combo$fdr <=.1), ]

  
  }
  
  #############################################################################################################################
  
  # Community change shifts and metabolites
  #############################################################################################################################
  # just a correlation statistic like the metadata one used. 

  
  table(correlations[["Milk.kefir"]]$species)
  
  table(correlations[["Water.kefir"]]$species)
  
  
  
  levels(as.factor(correlations[["Milk.kefir"]]$functional_feature))
  
  
  levels(as.factor(correlations[["Water.kefir"]]$functional_feature))
  
  
  
  #merged_data <- merge(  species_profile[[i]], metabolite_profile[[i]] , by.x=0,by.y = "merge_column",all.x=TRUE)
  
  

  
  
  #############################################################################################################################
  
  # Correlate metabolics by species RA
  # Correlate mutation rate by volatile change 
  #############################################################################################################################
  


  
  
  
  
  
  
  
  #############################################################################################################################
  
  # import all metabolics datasets and modify the exisitig kefir4all metabolote datasets
  #############################################################################################################################
  
  # 
  setwd("Q:/H2020 Master/Citizen Science Project/Results/07_metabolomics/other_datasets")
  # 
  # 
  # 
  # temp <- c()
  temp = list.files(pattern="*.xlsx", recursive = TRUE, full.names = TRUE)
  
  
  
  
  metabolic_data = lapply(temp,read_excel)
  
  
  #"Kefir4all.Milk.kefir","Kefir4all.Water.kefir",
  names(   metabolic_data) <- c("Gethins.Milk.kefir","Walsh.Milk.kefir","Breselge.Water.kefir")
  
  
metabolic_data[["Breselge.Water.kefir"]]$Sample <- gsub("-","_",metabolic_data[["Breselge.Water.kefir"]]$Sample )
  


 kefir_name <- c()
 volatiles <- c()
 volatiles_long <- c()
 

 metadata_metabolomics_datasets <- c()
 
 metadata_metabolomics_datasets[["Kefir4All_Milk.kefir"]] <-  metadata_metabolomics[["Milk.kefir"]]
 
 metadata_metabolomics_datasets[["Kefir4All_Milk.kefir"]]$conditions="Household conditions"
 
 
 metadata_metabolomics_datasets[["Kefir4All_Milk.kefir"]]$conditions[which(metadata_metabolomics_datasets[["Kefir4All_Milk.kefir"]]$Stage=="T0")] <- "Laboratory controlled"
 
 
 metadata_metabolomics_datasets[["Kefir4All_Water.kefir"]] <-  metadata_metabolomics[["Water.kefir"]]
 
 metadata_metabolomics_datasets[["Kefir4All_Water.kefir"]]$conditions="Household conditions"
 
 
 
 metadata_metabolomics_datasets[["Kefir4All_Water.kefir"]]$conditions[which(metadata_metabolomics_datasets[["Kefir4All_Water.kefir"]]$Stage=="T0")] <- "Laboratory controlled"

 
 
 
 metadata_metabolomics[[ "Walsh.Milk.kefir"]] <- 
read_delim("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/mk_walsh_et_al.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
 
colnames(metadata_metabolomics[[ "Walsh.Milk.kefir"]])[which(colnames( metadata_metabolomics[[ "Walsh.Milk.kefir"]])=="sample_id")] <- "Sample"
 metadata_metabolomics[[  "Breselge.Water.kefir"]] <- 
 global_wk_metadata 
 
 
 colnames(metadata_metabolomics[[  "Breselge.Water.kefir"]])[which( colnames(metadata_metabolomics[[  "Breselge.Water.kefir"]])=="Sample")] <- "Sample_id"
 
 metadata_metabolomics[[  "Breselge.Water.kefir"]]$Sample<- 
   
   paste(str_split_fixed( metadata_metabolomics[[ "Breselge.Water.kefir"]]$merge_column,  "_",5)[,1], str_split_fixed( metadata_metabolomics[[  "Breselge.Water.kefir"]]$merge_column,  "_",5)[,2], str_split_fixed( metadata_metabolomics[[  "Breselge.Water.kefir"]]$merge_column,  "_",5)[,3], str_split_fixed( metadata_metabolomics[[  "Breselge.Water.kefir"]]$merge_column,  "_",5)[,4],sep="_") #, str_split_fixed( metadata_metabolomics[[  "Breselge.Water.kefir"]]$Sample,  "_",5)[,5],sep="_")
 
 
 
 #gsub("_S\\d.*","",metadata_metabolomics[[  "Breselge.Water.kefir"]]$merge_column)
 
 #colnames(metadata_metabolomics[[  "Breselge.Water.kefir"]])[which( colnames(metadata_metabolomics[[  "Breselge.Water.kefir"]])=="merge_column")] <- "Sample"
 
# SampleSheet_KefirDanone2 <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/SampleSheet__KefirDanone2.csv")

 
#metadata_metabolomics[[  "Breselge.Water.kefir"]]$Sample <- 
 #gsub("-","_",metadata_metabolomics[[  "Breselge.Water.kefir"]]$Sample )
 

 library(readr)
 global_wk_raw_metabolomics <- read_csv("global_wk_raw_metabolomics.csv")
 
 for (kefir_name in c("Milk.kefir", "Water.kefir")){
   
   volatiles[[kefir_name]] <- 
     dplyr::select(myfiles_wide[[kefir_name]],-c(CAS, RI, `Ref RI`, `Reference odour descriptor`))
   
   
   
  for (i in names(metabolic_data)[grep(kefir_name, names(metabolic_data))]){
    
   # kefir_name <- paste(str_split_fixed(i,"\\.",3)[2],  str_split_fixed(i,"\\.",3)[3],sep=".")

 
    
    #any(is.na(volatiles[[kefir_name]]))
      
    tester <- 
t(metabolic_data[[i]] %>%  column_to_rownames("Sample"))
 
# 
# 
#     missing_1 <- 
# rownames(tester)[-c(
#         which(rownames(tester) %in% volatiles[[kefir_name]]$Compound))]
#    
#     
#       missing <-    
#  global_wk_raw_metabolomics[ 
# which(   
# global_wk_raw_metabolomics$`Name (common name brackets)` %in% rownames(tester)[-c(
#   which(rownames(tester) %in% volatiles[[kefir_name]]$Compound))]),c("Name (common name brackets)","CAS")]
#   
#       missing_1[-c(
#      which( missing_1 %in% missing$`Name (common name brackets)`))]
#   
#     nrow(missing)
# 
#     missing[
#     which(missing$CAS %in% myfiles[["Water.kefir"]]$CAS),]
    
 volatiles[[kefir_name]] <- 
   
   merge(
 tester,
  volatiles[[kefir_name]],
   by.x=0,by.y="Compound",
   all=TRUE) %>% relocate(Class,.before="Row.names")
 

 


 colnames( volatiles[[kefir_name]])[which(colnames( volatiles[[kefir_name]])=="Row.names")] <- "Compound"
 
 # if(i=="Breselge.Water.kefir"){
 #   
 #   # volatiles_long[[kefir_name]]$Sample <- 
 #   # gsub("-","_",volatiles_long[[kefir_name]]$Sample )
 #   
 #   metabolic_data[[i]]$Sample <- gsub("-","_", metabolic_data[[i]]$Sample)
 # }
 # 
 
 metadata_metabolomics_datasets[[i]] <- dplyr::select(metabolic_data[[i]],Sample) %>% mutate(conditions=as.character("Laboratory controlled"))
 
 if(i=="Gethins.Milk.kefir"){
   next
 }else{
   metadata_metabolomics_datasets[[i]] <- 
   merge(metadata_metabolomics_datasets[[i]],  metadata_metabolomics[[ i]] ,by="Sample",all.x=TRUE)
   
 }
 
 
 
 

  }
 
   
   volatiles_long[[kefir_name]] <- 
     volatiles[[kefir_name]] %>% pivot_longer(!c(Class,Compound),names_to = "Sample",values_to = "Concentrations")
   
   
   
 }
 



 metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$`kefir type` <-"ML" 
 metadata_metabolomics_datasets[["Gethins.Milk.kefir"]]$`kefir type` <-"ML" 
 
 metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$`kefir type`[ which(metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$grain=="Milk")] <- "Medium control"
 
 metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$`kefir type`[ which(metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$source=="kefir grain")] <- "MG"
 
 
 metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$`kefir type`[ which(is.na( metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$`kefir type`))] <- "Medium control"
 
 metadata_metabolomics_datasets[["Breselge.Water.kefir"]]
 

 metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$merge_column <-  metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$Sample

 metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$merge_column <-  paste( metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$merge_column,"_L001",sep="")

 metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$merge_column_vs= metadata_metabolomics_datasets[["Walsh.Milk.kefir"]]$Sample
 
 metadata_metabolomics_datasets[["Gethins.Milk.kefir"]]$merge_column_vs= metadata_metabolomics_datasets[["Gethins.Milk.kefir"]]$Sample
 
 
 metadata_metabolomics_datasets[["Gethins.Milk.kefir"]]$merge_column= NA
 metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$merge_column_vs <-  metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$Sample

 

 
 
 volatile_metadata <- c()
 
 
 
 volatile_metadata[["Milk.kefir"]] <-  rbind(
   
   dplyr::select(metadata_metabolomics_datasets[["Walsh.Milk.kefir"]],Sample, `kefir type`, merge_column_vs, conditions,merge_column),
 dplyr::select(metadata_metabolomics_datasets[["Gethins.Milk.kefir"]],Sample, `kefir type`,merge_column_vs,conditions,merge_column),
               dplyr::select(metadata_metabolomics_datasets[["Kefir4All_Milk.kefir"]],Sample, `kefir type`, merge_column_vs,conditions,merge_column)
               
 )
 
 metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$merge_column[
which(is.na(metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$merge_column))] <- metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$Sample[which(is.na(metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$merge_column))]
 
 
 metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$data_source[
   which(is.na(metadata_metabolomics_datasets[["Breselge.Water.kefir"]]$merge_column))] <- "Mortensen et al 2023" 
 
 
 volatile_metadata[["Water.kefir"]] <-  rbind(
   
   dplyr::select(metadata_metabolomics_datasets[["Breselge.Water.kefir"]],Sample, `kefir type`,merge_column_vs,conditions,merge_column),
   dplyr::select(metadata_metabolomics_datasets[["Kefir4All_Water.kefir"]],Sample, `kefir type`, merge_column_vs,conditions,merge_column)
   
 )
 
 volatile_metadata[["Water.kefir"]] <-  rbind(volatile_metadata[["Water.kefir"]],
                                              c("ID113","WL","113_T2","Household conditions", ""))

 
 #############################################################################################################################
 
 # alpha diversity time baby
 #############################################################################################################################
 my.files_summary_shannon <- c()
 my.files_richness <- c()
 my.files_evenness <- c()
 my.files_beta <- c()
 my.files_summary<-c()
 
 

 mets <- c()
 
 
 
 for (i in  c("Milk.kefir", "Water.kefir")){
   
   
   volatiles[[i]] <- 
   volatiles[[i]][, which(colnames( volatiles[[i]]) %in%  c("Class",     "Compound",  volatile_metadata[[i]]$merge_column_vs [which( volatile_metadata[[i]]$`kefir type` %in% c("ML","WL")  )] )) ]
   
   
   mets <- dplyr::select( volatiles[[i]], Class,Compound)
   
   volatiles[[i]] <- 
   dplyr::select( volatiles[[i]],-c(Class)) %>% column_to_rownames("Compound")

  volatiles[[i]][is.na(volatiles[[i]] )] <- 0
   my.files_summary_shannon[[i]] <- vegan::diversity(  t(   volatiles[[i]]))
   
   my.files_richness[[i]] <- specnumber(   t(   volatiles[[i]]) )
   my.files_evenness[[i]] <- my.files_summary_shannon[[i]]/log(my.files_richness[[i]])
   my.files_beta[[i]] <- vegdist(   t(   volatiles[[i]]) , method = "bray")
   
   my.files_summary[[i]]<- cbind(shannon = my.files_summary_shannon[[i]], richness = my.files_richness[[i]], pielou = my.files_evenness[[i]],site= rownames( t(   volatiles[[i]])))
   
   
   my.files_summary[[i]]<- as.data.frame(my.files_summary[[i]])
   my.files_summary[[i]] <- merge(my.files_summary[[i]], volatile_metadata[[i]],by.x="site",by.y="merge_column_vs",all.x=TRUE)
   
   
   
   my.files_summary[[i]]$shannon <- as.numeric( my.files_summary[[i]]$shannon)
   #setwd("Q:/H2020 Master/Citizen Science Project/Plots/Evolution")
   
   #jpeg(filename='Alpha diversity_kefir_type.jpeg', width = 35*700, height=30*700,res=1700,pointsize = 15) #, width=2000, height=1950)
   
   
 }
 
 
 my.files_summary[["Water.kefir"]]$group <- "Water kefir"
 my.files_summary[["Milk.kefir"]]$group <- "Milk kefir"
 
 metabolites_plot_data <- rbind( my.files_summary[[1]],
                                 my.files_summary[[2]] )
 

 
 
 p_alpha <- 
  grouped_ggbetweenstats(
   data = metabolites_plot_data,
   x = conditions,
   y = shannon,
   fill=conditions,
   color=conditions,
   grouping.var =group,
   type = "nonparametric", # ANOVA or Kruskal-Wallis
   plot.type = "box",
   pairwise.comparisons = TRUE,
   pairwise.display = "significant",
   centrality.plotting = FALSE,
   ggsignif.args    = list(textsize = 4, tip_length = 0.01),
   bf.message = FALSE, 
   xlab="Time points",
   ylab="Alpha diversity values (Shannon)",
   fill = "Timeframe",
   ggtheme = ggplot2::theme_bw(),
   ggplot.component = list(theme(
     plot.title = 
       ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
                                      box.color = "#748696",size=20, lineheight = 2),
     legend.title = element_text( size=25, face="bold"),
     axis.text.x = element_text(hjust = 1, size = 15),
     axis.text.y = element_text(hjust = 1, size = 10),
     axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
     axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
     legend.position = "none"),
     ggplot2::scale_color_manual(values = c("#F8766D", "#00BFC4")) ))      
 
   
   
   
   
   
   
   
   
   
   
 
 
 
 
 
 

 #############################################################################################################################
 
 # beta diversity time baby
 #############################################################################################################################

 t <- c()
 my.files_beta_total <-c()
 my.files_beta_cmdscale_total <- c()
 pcoa_total <- c()
 en_pcoa_total <- c()
 envfit_output_total <-c()
 pcoa_data_total <- c()
 
 env_correlations_total <- c()
 base_total <- c()
 adonis_total <- c()
 
 
 wilcoxon_data <- c()
 
 
 wilcoxon <- c()
 
 wilcoxon_plot <- c()
 
 plot <- c()
 
 #volatiles[["Water.kefir"]] <- 
 #volatiles[["Water.kefir"]] %>% dplyr::select(-c(Class)) %>% column_to_rownames("Compound")
 
 
 
 lowest <- c()

 
 highest <- c()
   
 volatile_details <- 
  rbind( 
 dplyr::select(myfiles[[1]],Class, Compound),
 
 dplyr::select(myfiles[[2]],Class, Compound), 
 
 dplyr::select(global_wk_raw_metabolomics,Class, Compound))
 
 volatile_details <- 
 volatile_details[-c(which(duplicated(volatile_details$Compound))),]
 
 
 for (i in names(myfiles)){
   
   

   my.files_beta_total[[i]] <- vegdist( t(volatiles[[i]]))
   
   
   
   my.files_beta_cmdscale_total[[i]] <-  cmdscale(    my.files_beta_total[[i]],k=2,eig=TRUE,add = TRUE) #ordination
   
   
   
   pcoa_total[[i]] = data.frame(PC1 =  my.files_beta_cmdscale_total[[i]]$points[,1], PC2 =  my.files_beta_cmdscale_total[[i]]$points[,2])
   percent_explained <- 100*  my.files_beta_cmdscale_total[[i]]$eig/sum( my.files_beta_cmdscale_total[[i]]$eig)
   
   percent_explained[1:2]
   pcoa_total[[i]] <- merge(    pcoa_total[[i]], volatile_metadata[[i]],by.x=0,by.y="merge_column_vs",all.x=TRUE)
   
   t1= nrow(t(volatiles[[i]]))
   
   t <- merge(t(volatiles[[i]]),
              dplyr::select(volatile_metadata[[i]],-c(`kefir type`, Sample)),
              by.x=0,
              by.y="merge_column_vs",
              all.x=TRUE)
   
   
   
   t2= nrow(t)
   if(t1!=t2){
     print(paste("merge command for",i," is not correct revise now",sep=""))
     #break
   }
   
   
   
   t <- t %>% column_to_rownames("Row.names")
   
   
   
   
   en_pcoa_total <- vegan::envfit(   my.files_beta_cmdscale_total[[i]],   t,  permutations = 1000,na.rm = TRUE)
   
   en_pcoa_total$vectors$pvals_fdr <-  p.adjust(en_pcoa_total$vectors$pvals,method = "bonferroni")
   
   
   envfit_output_total[[i]] <- en_pcoa_total
   
   
   
   
   # ggplot(   pcoa[[i]], aes(PC1, y=PC2,colour= Stage,fill= Stage))+
   #   geom_point(size=5) +
   # 
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
   #   theme(#legend.position = "none",#axis.text.x = element_blank(),  # remove x-axis text
   #     #axis.text.y = element_blank(), # remove y-axis text
   #     axis.ticks = element_blank(),  # remove axis ticks
   #     axis.text = element_blank(),
   #     axis.title.x = element_text(size=18), # remove x-axis labels
   #     axis.title.y = element_text(size=18), # remove y-axis labels
   #     panel.background = element_blank(), 
   #     panel.grid.major = element_blank(),  #remove major-grid labels
   #     panel.grid.minor = element_blank(),  #remove minor-grid labels
   #     plot.background = element_blank(),
   #     legend.text=element_text(face="italic",size = 15),
   #     legend.key.size = unit(2, 'cm'), #change legend key size
   #     legend.key.height = unit(2, 'cm'), #change legend key height
   #     legend.key.width = unit(2, 'cm'), #change legend key width
   #     legend.title = element_text(size=20) #change legend title font size
   #   )
   
   
   A <- as.list(en_pcoa_total$vectors)
   pvals<-as.data.frame(A$pvals)
  # pvals_fdr<-as.data.frame(A$pvals_fdr)
   arrows<-as.data.frame(A$arrows*sqrt(A$r))
   C<- cbind(arrows, pvals, as.data.frame(A$r))
   
   Cred<-subset(C,pvals<0.05)
   Cred <- cbind(Cred, Species = rownames(Cred))
   
   colnames(Cred)[which(colnames(Cred)=="A$r")] <- "cor"

   for (target in c("conditions")){
     
     
     t_n <-   as.data.frame(as.matrix(my.files_beta_total[[i]]))
     t_metadata <- volatile_metadata[[i]]
     
     
     t_n <-  
     t_n[c(which(rownames(t_n )%in% t_metadata$merge_column_vs )),]
     
     t_metadata <-  t_metadata[which(t_metadata$merge_column_vs %in% rownames(t_n )), ]
     
     
     t_metadata <-  t_metadata[match(rownames(t_n ),t_metadata$merge_column_vs), ]
     

     colnames(     t_metadata)[which(colnames( t_metadata)==target)] <- "interest"
    # t_n <- t_n[which(t_n)]
    
     
     adonis_total[[i]][[target]] <- 
       adonis2( t_n~interest, data=   t_metadata  ,permutations =  10000)
     
     
     
     #library(ltm)
     
     mock_data <- data.frame(volatiles=as.character(colnames(t)[-c(which(colnames(t) %in% c("conditions", "merge_column") ))]),
     Biserial=as.numeric(NA),
     Wilcoxon_p=as.numeric(NA),
     Wilcoxon_W=as.numeric(NA))
     
   wilcoxon <- c()
     for (v in colnames(t)[-c(which(colnames(t) %in% c("conditions", "merge_column") ))]){
       
       

       mock_data$Biserial[which(mock_data$volatiles==v)] <-     biserial.cor(getElement(t,v), getElement(t,"conditions"))
       
       wilcoxon <-
       wilcox.test( t[which(t$conditions=="Household conditions"),which(colnames(t)==v)], 
                    t[which(t$conditions=="Laboratory controlled"),which(colnames(t)==v)], alternative = "two.sided")
       
       
       mock_data$Wilcoxon_p[which(mock_data$volatiles==v)] <- wilcoxon$p.value
         
       mock_data$Wilcoxon_W[which(mock_data$volatiles==v)] <- wilcoxon$statistic
         
       # 
       plot <- t
       
       colnames(plot)[which(colnames(plot)==v)] <-"metabolite"
       wilcoxon_plot[[i]][[target]][[v]] <- ggbetweenstats( # independent samples
         data = plot,
         x = conditions,
         y = metabolite,
         title=v,
         plot.type = "box", # for boxplot
         type = "nonparametric", # for wilcoxon
         centrality.plotting = FALSE # remove median
       )

       
       wilcox.test( t[which(t$conditions=="Household conditions"),which(colnames(t)==v)], 
                       t[which(t$conditions=="Laboratory controlled"),which(colnames(t)==v)], alternative = "two.sided")
       
       
       
   
       
       
     }
     
     mock_data$fdr <- p.adjust(mock_data$Wilcoxon_p, method = "bonferroni")
   wilcoxon_data[[i]] <- mock_data
   
     
     metabolite_profile[[i]] <- metabolite_profile[[i]][match(rownames(  species_profile[[i]]),rownames(metabolite_profile[[i]])), ]
     
  
     data.scores_plot <-     
       pcoa_total[[i]] 
     
     
     colnames(data.scores_plot)[which(      colnames(data.scores_plot)==target)] <- "target"
     
     #something to consider in a minute 
     
     # if(target=="species"){
     #   if (i=="Milk.kefir"){
     #     
     #     data.scores_plot$target[-c(which(data.scores_plot$target %in% c("Lactobacillus helveticus","Lactobacillus kefiranofaciens","Lactococcus cremoris", "Lactococcus lactis")))] <- "Other"
     #     #work here 
     #     
     #     
     #     
     #   }else{
     #     data.scores_plot$target[-c(which(data.scores_plot$target %in%  c("Lacticaseibacillus paracasei", "Lentilactobacillus hilgardii", "Zymomonas mobilis")))] <- "Other"
     #     
     #   }
     # }
     
     centroid <- data.frame(target=as.character(levels(as.factor(data.scores_plot$target))),
                            V1=as.numeric(0),
                            V2=as.numeric(0))
     
     for (e in centroid$target){
       centroid$PC1[which(centroid$target==e)] <- mean(data.scores_plot$PC1 [which(data.scores_plot$target ==e)])
       centroid$PC2[which(centroid$target==e)] <- mean(data.scores_plot$PC2[which(data.scores_plot$target==e)])
       
     }#end of e
     
     
     
     #test <- 
     #dplyr::select( data.scores, c( Row.names, V1, V2, levels(as.factor( mech_total_resistome$Mechanism)))) %>% pivot_longer(!c( Row.names, V1, V2),values_to = "RA", names_to = "mech") 
     
     
     library(ggnewscale)
     
     colnames(data.scores_plot)[which(colnames(data.scores_plot)=="target")] <- "Fermentation parameter"
     
     pcoa_data_total[[i]][[target]] <- 
       
       # ggplot(data = data.scores_plot, aes(x = PC1, y = PC2)) + 
       # geom_point(data = data.scores_plot, aes(colour = `Fermentation parameter`), size = 3, alpha = 0.5) + 
       # geom_point(data=centroid,size=10,shape=21, color="black",aes(fill=target),  show.legend=FALSE)+ 
       # stat_ellipse(geom = "polygon",
       #              aes(fill= `Fermentation parameter`),
       #              alpha = 0.25,
       #              type = "norm",
       #              show.legend=FALSE)+
       # #scale_colour_manual(values = c("orange", "steelblue"))  + 
       # #scale_fill_manual(values = c("orange", "steelblue"))  + 
       # new_scale_color()+
       # # geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
       # #              data =  Cred, size =1, alpha = 0.5, colour = "grey30") +
       # geom_point(data = Cred, aes(x = Dim1, y = Dim2, colour = cor), 
       #            shape = "diamond", size = 4, alpha = 0.6) +
       # geom_text_repel(data = Cred[which(Cred$`A$pvals`<0.05),], aes(x = Dim1, y = Dim2, colour = cor), 
       #                 label = row.names(Cred), fontface = "bold",min.segment.length=0.01) + 
       # #geom_text(data = Cred, aes(x = Dim1, y = Dim2), colour = "grey30", 
       # #fontface = "bold", label = row.names(Cred)) + 
       # theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
       #       panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
       #       axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
       #       legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
       #       legend.text = element_text(size = 9, colour = "grey30"),
       #       legend.position = c(.075,.8),
       #       legend.direction ="horizontal")+
       #       #legend.direction="vertical")+
       # #legend.box="Horizontal") + 
       # labs(colour = "Fermentation conditions",
       #      x="PC-1",
       #      y="PC-2")+
       # scale_color_distiller(palette = "RdYlBu")
       # 
     
     
     
     ggplot(data = data.scores_plot, aes(x = PC1, y = PC2)) + 
       geom_point(aes(colour = `Fermentation parameter`), size = 3, alpha = 0.5) + 
       geom_point(data=centroid, size=10, shape=21, color="black", aes(fill=target), show.legend=FALSE)+ 
       stat_ellipse(geom = "polygon",
                    aes(fill= `Fermentation parameter`),
                    alpha = 0.25,
                    type = "norm",
                    show.legend=FALSE) +
       new_scale_color() +
       geom_point(data = Cred, aes(x = Dim1, y = Dim2, colour = cor), 
                  shape = "diamond", size = 4, alpha = 0.6) +
       geom_text_repel(data = Cred[Cred$`A$pvals` < 0.05,], aes(x = Dim1, y = Dim2, colour = cor), 
                       label = row.names(Cred), fontface = "bold", min.segment.length = 0.01)+
       
       
       theme( axis.title.x = element_text(size = 15, face = "bold", hjust = 0.5, vjust = -1),
              axis.title.y = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 2),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(hjust = 1, size = 10),
             panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
             axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
             legend.title = element_text(size = 15),
             legend.text = element_text(size = 15),
             legend.position = "top",
             legend.box = "horizontal",
             legend.direction = "horizontal",
             
             legend.key.size = unit(1, 'cm'), # change legend key size
             legend.key.height = unit(1, 'cm'), # change legend key height
             legend.key.width = unit(1, 'cm'), # change legend key width
             plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
                                                         box.color = "#748696",size=20, lineheight = 2)) +
       labs(colour = "Fermentation conditions",
            x=paste("PCoA1 - ", round(percent_explained[1]), "%", sep=""), y=paste("PcoA2 - ", round(percent_explained[2]), "%", sep=""),
            title=gsub("\\."," ",i)) +
       scale_color_distiller(palette = "RdYlBu") 
     
     
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
   
   
   env_correlations_total[[i]] <- Cred
   
   
   #paste(env_correlations[[i]]$Species[which(env_correlations[[i]]$rank<=5)],sep="",collapse = ", ")
   
   # print(paste("in", tolower(gsub("\\."," ",i)), length(env_correlations[[i]]$Species), "statistically significant volatiles contributed to differences across", tolower(gsub("\\."," ",i)) ,
   #             "metagenomes, the most distinguishing of which included ",
   #             
   #             paste(tolower(env_correlations[[i]]$Species[which(env_correlations[[i]]$rank<=5)])," (", round(env_correlations[[i]]$cor[which(env_correlations[[i]]$rank<=5)],2), ")",collapse = ", ",sep="")
   #             
   # ))
   
   
   # base[[i]] <- 
   #   
   #   dplyr::select(myfiles[[i]],sample_id, category_confirmed)
   # 
   # base[[i]] <- 
   #   base[[i]][-c(which(duplicated(base[[i]]$sample_id))),]
   # 
   
   
   wilcoxon_data[[i]] <- 
   wilcoxon_data[[i]][which(wilcoxon_data[[i]]$fdr<=0.05),]
   
   
  
   
   wilcoxon_data[[i]] <- 
     merge(
       wilcoxon_data[[i]] , 
       volatile_details,
       by.x="volatiles",
       by.y="Compound",
       all.x=TRUE)
   
     lowest <- 
 wilcoxon_data[[i]] %>% 
     filter(Biserial < 0) %>%
     slice_min(order_by = Biserial, n = 5) 
  
     highest <- 
   wilcoxon_data[[i]] %>% 
     filter(Biserial > 0) %>%
     slice_max(order_by = Biserial, n = 5)
     
     
     n <- as.data.frame(table(wilcoxon_data[[i]]$Class[which(wilcoxon_data[[i]]$Biserial< 0)] ))
     
     n <- n[order(n$Freq, decreasing = TRUE),]

     
     print(paste(i, "- laboratory has", nrow(   wilcoxon_data[[i]] %>% 
                                     filter(Biserial < 0)), "volatile"," belonging to the following classes", paste(n$Var1,"(",
                                                                                                 n$Freq,")",collapse=", "),
   
   paste( lowest$class, "-", lowest$volatiles, "(Rank biserial =", -(round(lowest$Biserial,2)),", p <0.01 )" ,collapse = ", "), "were significantly higher in",i, "samples produced under household conditions compared to laboratory conditions")
     )
     
     
     
     
     n <- as.data.frame(table(wilcoxon_data[[i]]$Class[which(wilcoxon_data[[i]]$Biserial> 0)] ))
     
     n <- n[order(n$Freq, decreasing = TRUE),]
     
     # problem is here
     print(paste(i, " - household has", nrow(   wilcoxon_data[[i]][which(wilcoxon_data[[i]]$Biserial>0),]), "volatile"," belonging to the following classes", paste(n$Var1,"(",
                                                                                                 n$Freq,")",collapse=", "),
   paste( highest$volatiles, "(Rank biserial =", (round(highest$Biserial,2)),", p <0.01 )" ,collapse = ", "),  "between",i ,"samples produced under laboratory conditions compared to household conditions.")
     )
 
     
     


    
     
     }# end of for loop for i
 
 
 

 # plus is for household and minus is for lab 
 
 wilcoxon
 wilcoxon_data$Milk.kefir[
 which(wilcoxon_data$Milk.kefir$fdr<=0.05),]
 
 
 

 
 
 
 wilcoxon_data$Water.kefir[
   which(wilcoxon_data$Water.kefir$fdr<=0.05),]
 
 
 
 
 

 
 
 

 volatiles_long[["Milk.kefir"]]
 
which(is.na(total_tile_plot_data$concentrations[which(
 total_tile_plot_data$`kefir type`=="ML")]))
 


which(is.na(total_tile_plot_data$concentrations[which(
  total_tile_plot_data$`kefir type`=="WL")]))

levels(as.factor(
total_tile_plot_data$Compound[which(
  total_tile_plot_data$`kefir type`=="ML" & 
    total_tile_plot_data$concentrations>0
    )]))



levels(as.factor(
  total_tile_plot_data$Compound[which(
    total_tile_plot_data$`kefir type`=="WL" & 
      total_tile_plot_data$concentrations>0
  )]))


combined_data <-  rbind(
  wilcoxon_data[["Milk.kefir"]] %>% 
    dplyr::mutate(kefir_type = "Milk kefir"),
  wilcoxon_data[["Water.kefir"]] %>% 
    dplyr::mutate(kefir_type = "Water kefir"))


combined_data$conditions <- NA
combined_data$conditions[which(combined_data$Biserial>0)] = 'Household conditions'

combined_data$conditions[which(combined_data$Biserial<0)] =   'Laboratory controlled'



# 
#   
# combined_data$of_interest <- NA
# 
# combined_data$of_interest[which(combined_data$kefir_type=="Milk kefir" & !combined_data$volatiles %in% 
#                                   levels(as.factor(
#                                     total_tile_plot_data$Compound[which(
#                                       total_tile_plot_data$`kefir type`=="ML" & 
#                                         total_tile_plot_data$concentrations>0
#                                     )]))
# )] <- NA

   
   library(ggplot2)
 library(dplyr)
 library(ggforce)

p_biserial=
 combined_data %>% 
   dplyr::filter((Biserial > 0.3 & fdr <= 0.05) | 
                   (Biserial < -0.3 & fdr <= 0.05)) %>% 
   
   ggplot(aes(y = reorder_within(volatiles, Biserial, kefir_type), x = as.numeric(Biserial),fill = conditions)) +
   geom_bar(stat = "identity", color = "black", width = 1, size = 1) +
   ggforce::facet_col(kefir_type ~ ., scales = "free", space = "free") +
   scale_y_reordered() +
   theme_bw() +
   labs(y = "Volatiles", x = "Point-Biserial Correlation Coefficient", fill = "Correlation") +
   theme(
     legend.title = element_text(size = 15),
     legend.text = element_text(size = 15),
     axis.text.x = element_text(size = 9),
     axis.text.y = element_text(hjust = 1, size = 9),
     axis.title.x = element_text(size = 15, face = "bold", hjust = 0.5, vjust = -1),
     axis.title.y = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 2),
     legend.position = "top",
     strip.background = element_rect(color = "black", fill = "white"),
     strip.text = element_text(size = 20), 
     
     legend.key.size = unit(1, 'cm'), # change legend key size
     legend.key.height = unit(1, 'cm'), # change legend key height
     legend.key.width = unit(1, 'cm') # change legend key width
     
   ) 






   # Adding annotations for Laboratory controlled and Household
  #Unicode characters \u2192 (???) and \u2190 (???) are used to add right and left arrows respectively. 
  #To make the arrows bigger in your ggplot2 annotations, you can use larger Unicode characters for the arrows, such as:
  
 # \u27A1 (??????) for a right arrow.
#\u2B05 (??????) for a left arrow.
  
  
  # annotate("text", x = 0.5, y = Inf, label = "\u2192 Household", color = "blue", size = 6, hjust = 0, vjust = 9) +
  # annotate("text", x = -0.5, y = Inf, label = "Laboratory controlled \u2190", color = "red", size = 6, hjust = 1, vjust = 6)
  # 
  # 
  # 
 
 

 ########################################################################################################################
 #Import metacache data for correlation analysis total
 ########################################################################################################################
 metacache <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_species_profile_v2.csv")

 metacache$sample_id <- gsub("04_Metacache_individual/|04_Metacache_indiviudal_other_datasets/","",metacache$sample_id) 
 metacache$sample_id <- gsub("-","_",metacache$sample_id)
 metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]
 
 metacache [is.na(metacache )] <- 0
 
 
 metacache <- metacache %>% remove_rownames() %>% column_to_rownames("sample_id")
 species_correlations_total <- c()
 species_profile_total <- c()
 
 metabolite_profile_total <-c() 
 
 merged_data <- c()
 cor_matrix_total <- c()
 correlations_total <- c()
 combo_total <- c()
 
 ## function to tidy the outputs
 clus.boot <- c()
 
 tidy_corr <- function(x, value) {
   
   x %>%
     as.data.frame() %>%
     rownames_to_column("species") %>%
     filter(species %in% colnames(species_profile_total[[i]])) %>%
     dplyr::select(species, all_of(colnames(metabolite_profile_total[[i]]))) %>%
     pivot_longer(!species, names_to = "functional_feature", values_to = value)
   
 }
 
 
 
 count_gt_one <- function(x) {
   sum(x > .1)
 }
 
 
 cluster_quality_data <- c()

 indval <- c()
 chiSquare <- c()
 cooccur<-c()
 
 
 colnames(volatile_metadata[["Milk.kefir"]])[which(colnames(volatile_metadata[["Milk.kefir"]]) == "Row.names")] <- "merge_column"
 
  for (i in names(myfiles)){
    
    ########################################################################################################################
#data adjustment
    metabolite_profile_total[[i]] <-     t(   volatiles[[i]] ) #%>% column_to_rownames("Compound"))
    
    if(i=="Milk.kefir"){
    mets <- 
      volatile_metadata[[i]][-c(which(is.na(    volatile_metadata[[i]]$merge_column))),]
    }else{
    
    
    mets <-      volatile_metadata[[i]]
    
    }
    mets <- mets[which(mets$`kefir type` != "Medium control"),]
    
    
    mets <-
    mets[c(which(mets$merge_column_vs %in%  rownames( metabolite_profile_total[[i]] ))),]
    
    # "10_S8_L001"       

   species_profile_total[[i]] <- metacache[which(rownames(metacache) %in%       levels(as.factor(mets$merge_column))),]#& 
   
   species_profile_total[[i]] <-   na.omit(    species_profile_total[[i]] )
   
   
   maxab_species <- apply(species_profile_total[[i]],2, max, na.rm=TRUE)
   
   n1_species <-names(which(maxab_species >= .1))
   
   species_profile_total[[i]] <- species_profile_total[[i]][,which(colnames(species_profile_total[[i]]) %in%  n1_species )]
   
   
   maxab <- apply(species_profile_total[[i]],2,   count_gt_one)
   
  # n2_species <-names(which(maxab >= 6)) #round(length(levels(as.factor(rownames(species_profile[[i]]))))*.1)))
   
   
   #species_profile_total[[i]] <- species_profile_total[[i]][,which(colnames(species_profile_total[[i]]) %in%  n2_species )]

   
   metabolite_profile_total[[i]] <-  metabolite_profile_total[[i]][which( rownames(metabolite_profile_total[[i]]) %in% mets$merge_column_vs ), ]

   
   
   
   t2=   nrow(metabolite_profile_total[[i]])
   
   metabolite_profile_total[[i]] <- merge(   metabolite_profile_total[[i]],
                                            dplyr::select(mets,merge_column_vs,merge_column),
                                             by.x=0,
                                             by.y="merge_column_vs",
                                             all.x=TRUE
                                             
                                             )
   
   
   t3=   nrow(metabolite_profile_total[[i]])
   
   
   
   if(t2!=t3){
     print(paste("merge command for",i," is not correct revise now",sep=""))
     #break
   }
   
   
   
   metabolite_profile_total[[i]] <- dplyr::select(   metabolite_profile_total[[i]], -c(Row.names)) %>% column_to_rownames("merge_column")
   
   

   volatile_metadata[[i]] <- mets
   #rownames( metabolite_profile_total[[i]] )[-c( which(rownames( metabolite_profile_total[[i]] ) %in% mets$merge_column_vs)) ]
   
   
   

  # mets$merge_column_vs[-c(which(mets$merge_column_vs %in%  rownames( metabolite_profile_total[[i]] )))]
   

   
   metabolite_profile_total[[i]]  <- metabolite_profile_total[[i]][match(rownames(   
                                                                        species_profile_total[[i]]   ),rownames(metabolite_profile_total[[i]])), ]
   #mistake here 
   ########################################################################################################################
   #data correlation   
   if(
     identical(rownames(metabolite_profile_total[[i]]),  rownames(species_profile_total[[i]] ))!=TRUE){
     
     print(paste("species and metabolic dataframes do not align for",i,"crashing now"))
     break
   }
   
   
   #cor_matrix[[i]] <- cor(  species_profile[[i]],   metabolite_profile[[i]])
   
   
   cor_matrix_total[[i]] <- rcorr(as.matrix(species_profile_total[[i]]),  as.matrix( metabolite_profile_total[[i]]), "spearman")
   
   
   ## get the p-values and do p-value correction
   p_values <- cor_matrix_total[[i]]$P %>%
     tidy_corr(., "p") %>%
     mutate(fdr = p.adjust(p, method = "fdr"))
   
   ## get the r-values
   r_values <-
     cor_matrix_total[[i]]$r%>%
     tidy_corr(., "r")
   
   
   combo <- inner_join(r_values, p_values, by = c("species", "functional_feature")) %>%
     arrange(desc(r * r))
   
   combo$name <- gsub("\\."," ",i)
   
   
   
   correlations_total[[i]] <- combo[which(combo$p <= 0.05), ]
   
   # build a co_occurece network 
   
   

 
   

 
 # create absence/presence matrix, 1 denotes present
 pres_volatile <- ifelse( t(metabolite_profile_total[[i]]) ==0,0,1)

 pres_species <- ifelse( t(species_profile_total[[i]]) <0.01,0,1)
 
 
 threshold=0
 
 
 co_occurrence <- 
 pres_species %*%  t(pres_volatile)
 

 # 
 # 
 # 
 # # Function to create co-occurrence network
 # create_co_occurrence_network <- function(data1, data2, threshold = 0) {
 #   # Calculate co-occurrence matrix
 #   co_occurrence <- data1 %*% t(data2)
 #   
 #   # Apply threshold if necessary
 #   co_occurrence[co_occurrence <= threshold] <- 0
 #   
 #   # Convert co-occurrence matrix to binary adjacency matrix
 #   adjacency_matrix <- ifelse(co_occurrence > 0, 1, 0)
 #   
 #   # Create graph object
 #   graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
 #   
 #   return(graph)
 # }
 # 
 # 
 # co_occurrence_network <- create_co_occurrence_network(pres_species, pres_volatile)
 # 

 
 ########################################################################################################################
 #cooccur package
 
 if(identical(colnames(pres_volatile), colnames(pres_species))){
   
   
  # Calculate the observed co-occurrence matrix (correlation matrix)
   #observed <- cor(t(rbind(pres_volatile,pres_species)), method = "spearman")
  #  observed <- cor(cbind(metabolite_profile_total[[i]],species_profile_total[[i]]), method = "spearman")
  #  
  #  
  #  #observed_vec <- as.vector(observed)
  #  
  #  
  #  # Generate null distributions by permuting the columns (metabolite data) of the dataset
  # # rbind(pres_volatile,pres_species)
  #  null_distributions <- replicate(1000, {
  #    shuffled_data <- apply(cbind(metabolite_profile_total[[i]],species_profile_total[[i]]), 2, sample)
  #    cor(shuffled_data, method = "spearman")
  #  })
  #  
  #  observed_vec <- as.vector(observed)
  #  
  #  
  #  
  #  p_values <- sapply(seq_along(observed), function(e) {
  #    sum(null_distributions[e, ] >= observed[e]) / length(null_distributions[e, ])
  #  })
  #  
  #  # Create a data frame with sample IDs and corresponding p-values
  #  p_values_df <- data.frame(sample_id = sample_ids[common_samples], p_value = p_values)
  #  
  #  # Print the data frame
  #  print(p_values_df)
  #  
  #  
  # 
  #  
  #  p_values <- colMeans(null_distributions >= as.vector(observed) )
  #  
  #  colnames(p_values) <- rownames(observed)
  #  
  #  
  #  p_values_df <- data.frame(sample_id = sample_ids, p_value = p_values)
    
 co <- cooccur(rbind(pres_volatile,pres_species), spp_names = TRUE)

 
 co$results$sp1_name_group="volatile"
 
 co$results$sp1_name_group[ which(co$results$sp1_name %in% rownames( pres_species))] <- "species"

 
 co$results$sp2_name_group="volatile"
 
 co$results$sp2_name_group[ which(co$results$sp2_name %in% rownames( pres_species))] <- "species"
 
 
 
 co$results$compare_group <- "volatile_vs_volatile"
 
 co$results$compare_group[
 which(co$results$sp1_name_group=="species" &
   co$results$sp2_name_group=="volatile" |  

     co$results$sp2_name_group=="species" &
     co$results$sp1_name_group=="volatile" )] <- "species_vs_volatile"
 
 
 co$results$compare_group[
   which(co$results$sp1_name_group=="species" &
           co$results$sp2_name_group=="species")] <- "species_vs_species"
 
 

 co$results$p_lt_fdr <- 
 p.adjust(co$results$p_lt ,method = "bonferroni")
 

 co$results$p_gt_fdr <- 
   p.adjust(co$results$p_gt ,method = "bonferroni")
 
 
 co$results[
 which(co$results$p_gt_fdr<=0.05 |
         co$results$p_lt_fdr<=0.05),]
 
 
 
 co_sv <- co$results[which(co$results$compare_group == "species_vs_volatile"),]
 
 
 co_sv$p_lt_fdr <- 
   p.adjust( co_sv$p_lt ,method = "bonferroni")
 
 
 
 co_sv$p_gt_fdr <- 
   p.adjust( co_sv$p_gt ,method = "bonferroni")
 
 
 
 co$results$p_gt_fdr <- 
   p.adjust(co$results$p_gt ,method = "bonferroni")
 
 
 
 co_sv[
   which(co_sv$p_gt<=0.05 |
           co_sv$p_lt<=0.05),]
 
 
 co_sv[
   which(co_sv$p_gt_fdr<=0.05 |
           co_sv$p_lt_fdr<=0.05),]
 
 
 
 cooccur[[i]] <- co
 ########################################################################################################################
 #Hierarchial clustering 
 pacman::p_load(fpc,factoextra)
 
 test <- my.files_beta_total[[i]] 
 my.files_beta_total[[i]][is.na(my.files_beta_total[[i]])] <- 0
 
 

 beta <- 
 vegdist(  metabolite_profile_total[[i]])
 beta[is.na(beta)] <- 0
 
 

 hc <- hclust( beta,method="ward.D2")
                   
 

 
 # n <-
 # fviz_nbclust(as.matrix(beta[[i]]), hcut, method = "silhouette",print.summary = TRUE,verbose = TRUE)
 # 
 # 
 # n <- which(n$data$y==max(n$data$y)) # E
 # 
 
 
 cluster_quality_data[[i]] <-  data.frame(cluster=as.numeric(1:10),
                                     Milk.kefir=as.numeric(0),
                                     Water.kefir=as.numeric(0),
                                     mean=as.numeric(0),
                                     dissolution=as.numeric(0),
                                     recovery =as.numeric(0) )
   

 
 
 
 for (j in 2:10){
k=j
 
   
    

   clus.boot[[i]][[j]] <- clusterboot(hc,
                            distances=TRUE,
                            B=1000, # Number of bootstrap resamples
                            clustermethod=hclustCBI, # for hierarchical clustering 
                            method="ward.D2", # use what we used in "hclust"
                            k=j, 
                            count=FALSE) # 
   
   
   
   
   
   # print(  clus.boot[[i]][[j]],statistics=c("mean","dissolution","recovery")
   #         )
   
  
  #mean <-   
   cluster_quality_data[[i]][which( cluster_quality_data[[i]]$cluster==j), "mean"] <-  paste(clus.boot[[i]][[j]]$bootmean,collapse=",")
   cluster_quality_data[[i]][which( cluster_quality_data[[i]]$cluster==j), "recovery"] <-  paste(clus.boot[[i]][[j]]$bootrecover,collapse=",")
   cluster_quality_data[[i]][which( cluster_quality_data[[i]]$cluster==j), "dissolution" ] <-  paste(clus.boot[[i]][[j]]$bootbrd,collapse=",")
   
if(k==2){
  next
}else{
  k <- k-1
   
   if (   
     mean(clus.boot[[i]][[k]]$bootmean) >=
     mean(clus.boot[[i]][[j]]$bootmean)){
     print(paste("optimal solution reach for", i,"which is,",k=j, "cluster"))
     #n <- j
     break
    
   }else{
    
     
     
     
   }
   
   # mean(clus.boot[[i]][[3]]$bootmean)
   # mean(clus.boot[[i]][[4]]$bootmean)
   # 
   # 
   # 
   # 
   # mean(clus.boot[[i]][[3]]$bootrecover)
   # mean(clus.boot[[i]][[4]]$bootrecover)
   # 
   # mean(clus.boot[[i]][[3]]$bootbrd)
   # mean(clus.boot[[i]][[4]]$bootrecover)
   
} 
 
 }#end of j
   
 
 
 ########################################################################################################################
 #differnetial abundnace 
 
 
 
 
 ########################################################################################################################
 #cooccurence 
 
 
 indval[[i]] <- multipatt( species_profile_total[[i]],  clus.boot[[i]][[j]]$result$partition, 
                     control = how(nperm=999)) 
 

 
   volatile_metadata[[i]] <- 
 merge(
data.frame( clus.boot[[i]][[j]]$result$partition) %>% dplyr::rename("clusters"=1),
volatile_metadata[[i]],
by.x=0,
by.y="merge_column",
all=TRUE)
 
   
   
   
   

   #plot <- as.data.frame( xtabs(~ conditions+clusters, volatile_metadata[[i]] )) %>% pivot_wider(names_from = "clusters",values_from = "Freq")

   
   
   chiSquare[[i]] <- 
   chisq.test(as.matrix(xtabs(~ conditions+clusters, volatile_metadata[[i]] )))
   
  
                     
 }else{
   print(paste("For the cooccurnece in ",i," the volatile and species data doesn't align"))
 }
# co_occurrence_matrix <- table(species_profile_total[[i]]),as.matrix(  metabolite_profile_total[[i]])))

   
   #  combo[which(combo$fdr <=.1), ]
   
   
  # write.csv(  correlations_total[[i]], paste("Q:/H2020 Master/Citizen Science Project/Results/07_metabolomics/species_correlations",i,"species_correlations_total_data_set_01.csv"), row.names = FALSE , quote=FALSE)
  }
 
 
 
 hc <- hclust( my.files_beta_total[["Milk.kefir"]] ,method="ward.D2")

   
   clus.boot_total_mk<- clusterboot(hc,
                                      distances=TRUE,
                                      B=1000, # Number of bootstrap resamples
                                      clustermethod=hclustCBI, # for hierarchical clustering 
                                      method="ward.D2", # use what we used in "hclust"
                                      k=3, 
                                      count=FALSE) # 
   
   
   clus.boot_total_mk<- 
   
   merge(
     data.frame(clus.boot_total_mk$result$partition) %>% dplyr::rename("clusters"=1), 
   
 
   rbind(
     
     dplyr::select(metadata_metabolomics_datasets[["Walsh.Milk.kefir"]],Sample, `kefir type`, merge_column_vs, conditions,merge_column),
     dplyr::select(metadata_metabolomics_datasets[["Gethins.Milk.kefir"]],Sample, `kefir type`,merge_column_vs,conditions,merge_column),
     dplyr::select(metadata_metabolomics_datasets[["Kefir4All_Milk.kefir"]],Sample, `kefir type`, merge_column_vs,conditions,merge_column)
     
   ),
     by.x=0,
   by.y="merge_column_vs",
   all.x=TRUE)

 
 
 
   chisq.test(as.matrix(xtabs(~ conditions+clusters, clus.boot_total_mk )))
   
   
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 correlations_total[["Milk.kefir"]][which(correlations_total[["Milk.kefir"]]$fdr<= 0.05),]
 
 correlations_total[["Milk.kefir"]][which(correlations_total[["Milk.kefir"]]$fdr<= 0.05 &
                                            correlations_total[["Milk.kefir"]]$functional_feature=="Hexanal"),]
 
 
 correlations_total[["Water.kefir"]][which(correlations_total[["Water.kefir"]]$fdr<= 0.05),]
 
 
 # correlations_total[["Milk.kefir"]] %>% 
 #   ggplot(aes(x=species,col=species))+
 #   geom_bar()
 

 table(   correlations_total[["Milk.kefir"]]$species)
 
 
 
 
 
 table(   correlations[["Milk.kefir"]]$species)
 
 
 
 
 
 
 
 correlations_total[[1]][ which(correlations_total[[1]]$fdr<=0.05),]
 
 

 
 
 ########################################################################################################################
 #Compostional changes result in metabolic changes 
 ########################################################################################################################
 
 
 
 species_bray_curtis <- c()
 volatile_bray_curtis <- c()
  abund_temp=c()
 
  

  for (i in names(myfiles)){
    
    species_bray_curtis[["cs"]][[i]] <-  as.data.frame(as.matrix(vegdist( species_profile[[i]])))
    
    
    species_bray_curtis[["total"]][[i]] <-  as.data.frame(as.matrix(vegdist( species_profile_total[[i]])))
    
  
    
    volatile_bray_curtis[["cs"]][[i]] <- as.data.frame(as.matrix(vegdist( metabolite_profile[[i]])))
    
    volatile_bray_curtis[["total"]][[i]] <- as.data.frame(as.matrix(vegdist( metabolite_profile_total[[i]])))
    
  for (j in c("cs", "total")){
    
   



 if(identical(rownames( volatile_bray_curtis[[j]][[i]] ), rownames( species_bray_curtis[[j]][[i]]))==FALSE){
   
   print(paste("rownames for", i ,"do not align crashing now" )) 
   break
   
 }else{
   
  
   abund_temp[[j]][[i]] = mantel(as.matrix( species_bray_curtis[[j]][[i]]), as.matrix(volatile_bray_curtis[[j]][[i]]), method = "spearman", permutations = 9999, na.rm = TRUE)
   
   
   
 }
 
 
  }
    
  }
 
 
 ########################################################################################################################
 #Volatile changes over time 
 ########################################################################################################################

  
  mets <-   rbind(metadata_metabolomics[["Milk.kefir"]],
  metadata_metabolomics[["Water.kefir"]])
  
  
  mets$Stage <- gsub("T1","wk01",  mets$Stage )
  mets$Stage <- gsub("T2","wk05",  mets$Stage )
  mets$Stage <- gsub("T3","wk09",  mets$Stage )
  mets$Stage <- gsub("T4","wk13",  mets$Stage )
   mets$Stage <- gsub("T5","wk17",  mets$Stage )
   mets$Stage <- gsub("T6","wk21",  mets$Stage )
  
   brayplot_data <- c()
   brayplot_data_to <- c()
 for (i in names(myfiles)){

   t <- 
   
   volatile_bray_curtis[["cs"]][[i]] %>% rownames_to_column("sample_to") %>% pivot_longer(!sample_to,names_to = "sample_from",values_to = "bray")
   
 
   t1=nrow(t) 
   t <- 
   merge(t,
        mets,
         by.x="sample_to",
         by.y="merge_column",
         all.x=TRUE)

  t2=nrow(t)
  
  
  if(t1!=t2){
    print(paste("merge command for",i," is not correct revise now",sep=""))
    break
  }
  
  
   
  
  t1=nrow(t)
  
  
   t <- 
     merge(t,
            mets,
           by.x="sample_from",
           by.y="merge_column",
           all.x=TRUE)

   t2=nrow(t)
   
   if(t1!=t2){
     print(paste("merge command for",i," is not correct revise now",sep=""))
     #break
   }
  
   
   
   t$timecategory <- paste(t$Stage.y, t$Stage.x,sep="_")
   colnames(t)[which(colnames(t)=="Stage.x")] = "Stage.to"
   
   
   colnames(t)[which(colnames(t)=="Stage.y")] = "Stage.from"
   
   
   
   brayplot_data_to[[i]] <-  t[which(t$timecategory %in% c("T0_wk01","T0_wk05", "T0_wk09","T0_wk013", "T0_wk17","T0_wk21")), ]
   
   
   brayplot_data[[i]] <- t[which(t$timecategory %in% c("T0_wk01","wk01_wk05", "wk05_wk09","wk09_wk013", "wk13_wk17","wk17_wk21")), ]
   
   
   
   # ggbetweenstats(data =  brayplot_data[[i]] ,
   #                x=timecategory , 
   #                y=bray,
   #                title=i,
   #                type = "nonparametric", # A
   #                ggsignif.args    = list(textsize = 2, tip_length = 0.01)
   #                
   # )+
   #   theme_bw()+
   #   xlab("Time points")+
   #   ylab("Bray Curtis distance")+
   #   labs(fill = "Timepoint")+
   #   #guides(colour = guide_legend(override.aes = list(size=25)))+
   #   #ylim(0,250)+
   #   theme(#plot.title = element_text(hjust = 0.5,size=35,face="bold"),
   #     legend.title = element_text( size=25, face="bold"),
   #     axis.text.x = element_text( hjust = 1, size = 15),#angle = 45,
   #     axis.text.y = element_text(hjust = 1, size = 10),
   #     axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
   #     axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
   #     legend.position = "none",
   #     plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
   #                                                 box.color = "#748696",size=35, lineheight = 2))
   # 
   
   


   
   
   
   # y is sample to
  # x is sample from
   
 }
   
   plot_t0_change <- 
     
  
   grouped_ggbetweenstats(
     data =    
       rbind(
         brayplot_data_to[[1]],
         brayplot_data_to[[2]]) %>% mutate(Kefir.x=paste(Kefir.x," Kefir",sep="")),
     x=timecategory , 
     y=bray,
     fill=timecategory,
     color=timecategory,
     grouping.var =Kefir.x,
     #type = "nonparametric", # ANOVA or Kruskal-Wallis
     plot.type = "box",
     pairwise.comparisons = TRUE,
     pairwise.display = "significant",
     #centrality.plotting = FALSE,
     ggsignif.args    = list(textsize = 4, tip_length = 0.01),
     type = "nonparametric", # ANOVA or Kruskal-Wallis
  
   centrality.plotting = FALSE,
   bf.message = FALSE,
   xlab="Time points",
   ylab="Alpha diversity values (Shannon)",
   fill = "Timepoint",
   ggtheme = ggplot2::theme_bw(),
   ggplot.component = list(theme(
     plot.title = 
       ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
                                      box.color = "#748696",size=20, lineheight = 2),
     legend.title = element_text( size=25, face="bold"),
     axis.text.x = element_text(hjust = 1, size = 15),
     axis.text.y = element_text(hjust = 1, size = 10),
     axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
     axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
     legend.position = "none")))


     
     
     
 
   
   grouped_ggbetweenstats(
     data =    
       rbind(
         brayplot_data[[1]],
         brayplot_data[[2]])%>% mutate(Kefir.x=paste(Kefir.x," kefir",sep="")),
     x=timecategory , 
    y=bray,
     fill=timecategory,
     color=timecategory,
     grouping.var =Kefir.x,
     type = "nonparametric", # ANOVA or Kruskal-Wallis
     plot.type = "box",
     pairwise.comparisons = TRUE,
     pairwise.display = "significant",
     centrality.plotting = FALSE,
     ggsignif.args    = list(textsize = 4, tip_length = 0.01),
     bf.message = FALSE,
     xlab="Time points",
     ylab="Alpha diversity values (Shannon)",
     fill = "Timepoint",
     ggtheme = ggplot2::theme_bw(),
   ggplot.component = list(theme(
     plot.title = 
     ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
                                    box.color = "#748696",size=20, lineheight = 2),
   legend.title = element_text( size=25, face="bold"),
   axis.text.x = element_text(hjust = 1, size = 15),
   axis.text.y = element_text(hjust = 1, size = 10),
   axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
   axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
   legend.position = "none")))
   

   
   
   
   plot_data <-  rbind(
     brayplot_data_to[[1]],
     brayplot_data_to[[2]]) %>% mutate(Kefir.x=paste(Kefir.x," kefir",sep=""))
   
   
   for (type in c("Water kefir","Milk kefir")){
   for (timepoint in levels(as.factor(plot_data$Stage.to))){
     
     

     print(paste(timepoint,"-",type, "min distace is ", min(plot_data$bray[which(plot_data$Stage.to==timepoint)])))
     
     
     print(paste(timepoint,"-",type, "mean distace is ", mean(plot_data$bray[which(plot_data$Stage.to==timepoint)])))
     
     print(paste(timepoint,"-",type, "max distace is ", max(plot_data$bray[which(plot_data$Stage.to==timepoint)])))
     
     
   }
   
   }
   
   
   
   ggbetweenstats(
     data =  preprocessing_summary[[i]][which(preprocessing_summary[[i]]$n_of_reads >=r1 &
                                                
                                                preprocessing_summary[[i]]$n_of_reads<=r2),],
     
     x =study_accession,
     y = n_of_reads,
     fill=merge_column,
     color=merge_column,
     #grouping.var =diversity_metric,
     type = "nonparametric", # ANOVA or Kruskal-Wallis
     plot.type = "box",
     pairwise.comparisons = TRUE,
     pairwise.display = "significant",
     centrality.plotting = FALSE,
     ggsignif.args    = list(textsize = 2.5, tip_length = 0.01),
     bf.message = FALSE,
     
     xlab = "Study accession",
     ylab="Distribution of reads",
     ggtheme = ggplot2::theme_bw(),
     ggplot.component = list(theme(text = element_text(size = 20))))#+
   # ggplot2::scale_y_continuous(
   #   limits = c(r1, r2) #,
   # breaks = seq(from = 35, to = 85, by = 5)
   
   
   
   
 
 
   
   plot_data_same_participant  =    rbind(brayplot_data[[1]],
   brayplot_data[[2]])
   
   plot_data_same_participant <- plot_data_same_participant[which(plot_data_same_participant$Sample.x==plot_data_same_participant$Sample.y),]
  
   
   grouped_ggbetweenstats(
     data =    
       plot_data_same_participant,
     x=timecategory , 
     y=bray,
     fill=timecategory,
     color=timecategory,
     grouping.var =Kefir.x,
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
     ylab("Bray Curtis distance")+
     labs(fill = "Timepoint")+
     #guides(colour = guide_legend(override.aes = list(size=25)))+
     #ylim(0,250)+
     theme(#plot.title = element_text(hjust = 0.5,size=35,face="bold"),
       legend.title = element_text( size=25, face="bold"),
       axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
       axis.text.y = element_text(hjust = 1, size = 10),
       axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
       axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
       legend.position = "none",
       plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
                                                   box.color = "#748696",size=55, lineheight = 2))
   
     #guides(colour = guide_legend(override.aes = list(size=25)))+
     #ylim(0,250)+
     theme(plot.title = element_text(hjust = 0.5,size=35,face="bold"),
           legend.title = element_text( size=25, face="bold"),
           axis.text.x = element_text(hjust = 1, size = 15),
           axis.text.y = element_text(hjust = 1, size = 10),
           axis.title.x = element_text( size=15, face="bold",hjust = 0.5,vjust = -1),
           axis.title.y = element_text( size=15, face="bold",hjust = 0.5, vjust = 2),
           legend.position = "none")
   
   
   
   
   # library(ggpubr)
   # 
   # jpeg("Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 8.jpeg", width = 7864, height=5200,res =300,pointsize = 15) 
   # # 
   # # 
   # ggarrange(pa,
   #           ggarrange(    pcoa_data_total[["Milk.kefir"]][["conditions"]],   pcoa_data_total[["Water.kefir"]][["conditions"]], nrow=1,ncol=2, labels=c("B.","C."),common.legend = FALSE,  font.label = list(size = 30)),
   #      
   #           
   #           
   #           plot_t0_change,
     
   #           
   #           
   #           nrow=3,ncol=1,labels=c("A.","","D."),  font.label = list(size = 30), heights=c(1,.5,1))
   #          #ncol=1,nrow=2,labels = c("A.","","C."), font.label = list(size = 20),heights=c(.7,1))
   # 
   # # 
   # graphics.off()
   
   
   
   
   ########################################################################################################################
   #where does metabolic chnage result in chnages in volatiles
   ########################################################################################################################
   
   
   ########################################################################################################################
   # correlate bray curtis distance vs metabolic concentrations
   ########################################################################################################################
   t <- c()
   
   correlation_data <- c()
   res <- c()
   base_id <- c()
   for (type in c("Water.kefir","Milk.kefir")){
     
    t0 <-  
      myfiles[[type]][which(myfiles[[type]]$Stage=="T0"),]
    
    t0 <- 
    t0[-c(
   which( t0$concentrations==0)),]
    
     t <- 
       merge(
         brayplot_data_to[[type]],
         
        # my.files_beta[[type]], 
         #myfiles_superfocus[["p3"]][which(myfiles_superfocus[["p3"]]$sample_id %in%  my.files_beta[[type]]$merge_column),],
        myfiles[[type]],
         by.x="merge_column_vs.x",
         by.y="sample_id",
         all.x=TRUE)
     
     
     
     
    # t <- dplyr::select(t,merge_column,bray_distance, `Subsystem Level 3`,`reads_%`, merge_column.y )
     
     
     
     # if( 
     #   
     #   nrow(  myfiles[[type]][which(  myfiles[[type]]$sample_id %in%   brayplot_data_to[[type]]$merge_column_vs.x),]) !=
     #   nrow(t)
     # ){
     #   print("Bug in the merge command gonna crash till you fix it")
     #   stop()
     # }
     # 
     
     for (i in levels(as.factor(t$Compound))){
       
       
       if( i %in% t0$Compound){
         
         base_id <- "Detected in baseline"
       }else{
         paste(i," of ",type," is not detected in baseline",sep="")
         base_id <- "Not detected in baseline"
       }
       
       
     
     
     
     
     
    
       
      
       if(length(as.numeric( t$bray[which( t$Compound==i)]))<4){
         
         correlation_data <- rbind(data.frame(compound=as.character(i),
                                              kefir_type=as.character(type),
                                              p_value=as.character("n to small"),
                                              r=as.character("n to small"),
                                              base_detection=as.character(base_id)),
                                   correlation_data)
         
         
       }else{
         
         res <- cor.test(as.numeric(t$bray[which( t$Compound==i)]),
                         
                         as.numeric(t$concentrations[which( t$Compound==i)]),
                         method="spearman")
         
         
         
         
         correlation_data <- rbind(data.frame(compound=as.character(i),
                                              kefir_type=as.character(type),
                                              p_value=as.numeric(res$p.value),
                                              r=as.numeric(res$estimate),
                                              base_detection=as.character(base_id)),
                                   correlation_data)
         
         
       }
       
     }
     
     
   }
   
   
   ########################################################################################################################
   # breakdown correlations based on strength
   ########################################################################################################################
   
   # correlation_data <- 
   #   correlation_data[-c(which(correlation_data$p_value=="n to small")),] 
   # 
   correlation_data$fdr <- p.adjust(as.numeric(correlation_data$p_value),method = "bonferroni")
   
   
   
   
   
   
   strong_correlations_graph <- correlation_data[which(
                                                                                        correlation_data$r>=.8 & correlation_data$fdr<=0.05 |
                                                                                        correlation_data$r<0 &
                                                                                        correlation_data$r<= -.8 &
                                                                                        correlation_data$fdr<=0.05 
                                                                                        ),]
   
   
   correlation_data <- 
     correlation_data[-c(
       which(is.na( correlation_data$fdr))),]
   
   correlations_graph <- correlation_data[which(
     correlation_data$r>.3 & correlation_data$fdr<=0.05 |
       correlation_data$r<0 &
       correlation_data$r< -.3 &
       correlation_data$fdr<=0.05 
   ),]
   
   

   
   
   table( correlations_graph$kefir_type)
   
   table( correlations_graph$kefir_type[which(correlations_graph$r>=.3 & correlations_graph$fdr<=0.05  )])
   
   
correlations_graph[which(correlations_graph$r>.3 & correlations_graph$fdr<=0.05  ),]

 range(correlations_graph$r[which(correlations_graph$r>.3 & correlations_graph$fdr<=0.05  )])

   table( correlations_graph$kefir_type[which( correlations_graph$r<0 & correlations_graph$r< -.3 & correlations_graph$fdr<=0.05  )])
   
   correlations_graph[which( correlations_graph$r<0 & correlations_graph$r< -.3 & correlations_graph$fdr<=0.05  ),]
   
   
   range( correlations_graph$r[which( correlations_graph$r<0 & correlations_graph$r< -.3 & correlations_graph$fdr<=0.05  )])
   

   
   ########################################################################################################################
   # find the top 10 highest positive and negative correlations
   ########################################################################################################################
   
   
   correlations_graph$rank_mk <- NA
  
   
    correlations_graph $rank_mk[which( correlations_graph $kefir_type=="Milk.kefir")] <- 
   rank(-  correlations_graph $r[

  which( correlations_graph $kefir_type=="Milk.kefir")], )
   
   
   
    correlations_graph $rank_wk <- NA
   
   
    correlations_graph $rank_wk[which( correlations_graph $kefir_type=="Water.kefir")] <- 
     rank(-  correlations_graph $r[
       
       which( correlations_graph $kefir_type=="Water.kefir")], )
   
   

   
wk_up_range <- 
   
   range(
    correlations_graph $rank_wk[-c(which(is.na(  correlations_graph $rank_wk)))])[1]
   
   
wk_down_range <- 
   range(
      correlations_graph $rank_wk[-c(which(is.na(  correlations_graph $rank_wk)))])[2]


wk_up_range:c(wk_up_range+9)

wk_down_range:(wk_down_range -9)


                    

mk_up_range <- 
  
  range(
     correlations_graph $rank_mk[-c(which(is.na(  correlations_graph $rank_mk)))])[1]


mk_down_range <- 
  range(
     correlations_graph $rank_mk[-c(which(is.na(  correlations_graph $rank_mk)))])[2]


mk_up_range:c(mk_up_range+9)

mk_down_range:(mk_down_range -9)



cor_plot_data <-  correlations_graph [which( correlations_graph $rank_mk %in% c(mk_up_range:c(mk_up_range+9), mk_down_range:(mk_down_range -9)) | 
                                          
                                           correlations_graph $rank_wk %in% c(wk_up_range:c(wk_up_range+9), wk_down_range:(wk_down_range -9)) 
),                      ]
   
     # merge(correlation_data[which(correlation_data$fdr<=0.05 &
     #                                correlation_data$r>=.8 | 
     #                                correlation_data$r>= -.8 & 
     #                                correlation_data$fdr<=0.05 &
     #                                correlation_data$r<0),],t,
     #       by.x="compound",
     #       by.y="Compound",
     #       all.x="TRUE"
     # )
   
   
cor_plot_data$r <- as.numeric(cor_plot_data$r)
   range(   strong_correlations_graph$r)
   
   library(ggforce)
   

   
   
   p=
     cor_plot_data%>% 
     mutate(kefir_type=gsub("\\."," ",kefir_type)) %>% 
     #mutate(pathway=gsub("EC 5.1.2.- ","",pathway)) %>% 
     ggplot(
       # aes(r,reorder(metabolite,r))) +
       aes(y=reorder_within(compound,r,kefir_type), x=as.numeric(r)))+
     geom_bar(stat="identity",fill="coral1",color="black",width = 1,size=1)+
     
     ggforce::facet_col(kefir_type~ ., scales = "free", space = "free")+
     scale_y_reordered() +
     #facet_wrap(~ kefir_type, scales = "free",ncol=1) +
     theme_bw() +
     labs(y = "Subsystem", x = "Spearman's rank correlation coefficient", fill = "Correlation") +
     theme(
       legend.title = element_text(size = 25, face = "bold"),
       axis.text.x = element_text(size = 10),
       axis.text.y = element_text(hjust = 1, size = 10),
       axis.title.x = element_text(size = 15, face = "bold", hjust = 0.5, vjust = -1),
       axis.title.y = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 2),
       legend.position = "top",
       strip.background = element_rect(color = "black", fill = "white"),
       strip.text = element_text(size = 12, face = "bold")
     )
   
   
   
  # or 

   
   
   
   ########################################################################################################################
   # plot all correlations of moderate strength note there is no strong correlations
   ########################################################################################################################
   
   
 
   
   
   p=
   correlations_graph %>% 
     mutate(kefir_type = gsub("\\.", " ", kefir_type)) %>% 
     ggplot(aes(y = reorder_within(compound, r, kefir_type), x = as.numeric(r))) +
     geom_bar(stat = "identity", fill = "coral1", color = "black", width = 1, size = 1) +
     ggforce::facet_col(kefir_type ~ ., scales = "free", space = "free") +
     scale_y_reordered() +
     geom_text(data = correlations_graph %>%
                 mutate(kefir_type = gsub("\\.", " ", kefir_type)) %>% filter(base_detection == "Not detected in baseline"),
               aes(label = "*", y = reorder_within(compound, r, kefir_type), x = as.numeric(r)),
               hjust = -0.5, vjust = .8 ,size = 10, color = "blue") +
     theme_bw() +
     labs(y = "Volatiles", x = "Spearman's rank correlation coefficient", fill = "Correlation") +
     theme(
       legend.title = element_text(size = 25, face = "bold"),
       axis.text.x = element_text(size = 10),
       axis.text.y = element_text(hjust = 1, size = 10),
       axis.title.x = element_text(size = 15, face = "bold", hjust = 0.5, vjust = -1),
       axis.title.y = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 2),
       legend.position = "top",
       strip.background = element_rect(color = "black", fill = "white"),
       strip.text = element_text(size = 2)
     )
   
   
   
   
   # 
   # library(ggpubr)
   # 
   # jpeg("Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 8_v2.jpeg", width = 7864, height=5200,res =300,pointsize = 15)
   # #
   # #
   # ggarrange(pa,
   #       
   # 
   # 
   #           #plot_t0_change,
   #          p,
   # 
   # 
   # 
   #           nrow=2,ncol=1,labels=c("A.","B."),  font.label = list(size = 30))#, heights=c(1,.5,1))
   # 
   # 
   # #
   # graphics.off()
   # 
   # 
   
   

   
   
   library(ggpubr)
   
   jpeg("Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 9_v2.jpeg", width = 7864, height=5200,res =300,pointsize = 15)
   #
   #
   ggarrange(p_alpha,ggarrange( pcoa_data_total[["Milk.kefir"]][["conditions"]],
                        pcoa_data_total[["Water.kefir"]][["conditions"]],
                        ncol=2,nrow=1,labels=c("B.", ""), font.label = list(size = 30)),
             
            
             #plot_t0_change,
            
             p_biserial, 
             
             
             nrow=3,ncol=1,labels=c("A.","D."),  font.label = list(size = 30), heights=c(.5,.5,1))
   
   
   #
   graphics.off()
   
   

   
   
   