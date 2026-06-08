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

#remove.packages("ggtree")
pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis)
pacman::p_load(ggpubr,readxl,devtools,taxize,rotl,ape,treeio,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales,ggnewscale)

library(BiocManager)
#BiocManager::install("ggtree")
#remotes::install_github("YuLab-SMU/tidytree")
#remove.packages("ggtree")
#remotes::install_github("YuLab-SMU/ggtree")
library(ggtree)
#might inlcude number of Cds and tRNAs in outer circle
#######################################################################################################
#importing tree file after the gtdbtk infer step. The tree was imported to itol were another tree file was downloaded 
  #snake.tree <- read.tree("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/gtdbtk.bac120.decorated.tree")
# 
# tree <- ggtree(snake.tree, layout = "circular", open.angle = 0)
# id <- snake.tree$tip.label
# data_labels <- as.data.frame(id)
# 



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



SampleSheet_KefirDanone2 <- read_csv(file.path(DATA_DIR, "SampleSheet__KefirDanone2.csv")

colnames(SampleSheet_KefirDanone2) <- SampleSheet_KefirDanone2[c(13),]
SampleSheet_KefirDanone2 <- SampleSheet_KefirDanone2[-c(1:13),]
library(readr)
mk_walsh_et_al <- read_delim(file.path(DATA_DIR, "mk_walsh_et_al.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


########################################################################################################################
# Merge metadata into one
########################################################################################################################
total_metadata <- rbind(dplyr::select(kefir4all_metadata, data_source, merge_column, `kefir type`,Stage,Sample),
                        dplyr::select(global_mk_metadata, data_source, merge_column, `kefir type`,Stage,Sample),
                        dplyr::select(global_wk_metadata, data_source, merge_column, `kefir type`,Stage,Sample))


total_metadata $category <- NA
total_metadata $category[which(total_metadata $`kefir type` %in% c("WL","WG"))] <- "Water.kefir"
total_metadata $category[which(total_metadata $`kefir type` %in% c("ML","MG"))] <- "Milk.kefir"


global_wk_metadata$country <- NA



colnames(global_wk_metadata)
colnames(global_mk_metadata)
global_wk_metadata <- dplyr::select(global_wk_metadata,"Sample" ,      "kefir type",   "Stage",        "Timepoint",    "country",      "data_source",  "merge_column")


global_metadata_total <- 
rbind(global_wk_metadata,
      global_mk_metadata)

global_metadata_total $category <- NA
global_metadata_total $category[which(global_metadata_total $`kefir type` %in% c("WL","WG"))] <- "Water.kefir"
global_metadata_total $category[which(global_metadata_total $`kefir type` %in% c("ML","MG"))] <- "Milk.kefir"



#######################################################################################################
# Importing checkm data
#######################################################################################################
setwd("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification")

temp = list.files(pattern="*quality_report.tsv", recursive = TRUE, full.names = TRUE)
myfiles_checkm = lapply(temp,read_table2, comment="--")

names_checkm <- temp
#names_checkm <- as.data.frame(names_checkm)

names_checkm <- gsub( "./|quality_report.tsv|02", "",names_checkm)


names(myfiles_checkm) <- names_checkm


checkm_total <- c()

#######################################################################################################
# merge checkm data togther 
#######################################################################################################

for (i in names(myfiles_checkm)){


  myfiles_checkm[[i]]$type <- i
  
  checkm_total <- rbind(checkm_total , myfiles_checkm[[i]] )
  
}

library(data.table)
#fwrite(  checkm_total, "Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/check_total_MAGs_v2.csv")

#######################################################################################################
#remove less than 50% complete 
#checkm_total <- 
#checkm_total[which(checkm_total$Completeness>=50),]

# checkm_1 <- checkm_total[which(checkm_total$Completeness>=50 &
#                                  checkm_total$Contamination<=10 ),]

checkm_1 <-   checkm_total[-c(which(checkm_total$type %in% c( "initial_cs_bin_metawrap_merge", "initial_cs_bin_metawrap_unmerge" ,    "initial_global_bin_metawrap_merge" ,  "initial_global_bin_metawrap_unmerge"))),]

t <- 
dplyr::select(checkm_1, Name, Completeness, Contamination)

colnames(t) <- c("genome","completeness","contamination")

t$genome <- paste(t$genome,".fa",sep="")
#write.csv(t,"Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/drep/drep_checkm_file.csv",row.names = FALSE,quote = FALSE)
#######################################################################################################
# Importing classification data 
#######################################################################################################

setwd("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification")
temp_classification = list.files(pattern="*bin_taxonomy.tab", recursive = TRUE, full.names = TRUE)
myfiles_classification = lapply(temp_classification,read_delim, 
                                delim = "\t", escape_double = FALSE, 
                                col_names = FALSE, trim_ws = TRUE)



names_classification  <- gsub( "./|bin_taxonomy.tab|02", "",temp_classification )


names(myfiles_classification ) <- names_classification 


classification_total <-c()
# note getting just eukaroytes here 

for (i in names(myfiles_classification)){
  
  
  myfiles_classification[[i]]$type <- i
  #myfiles_classification[[i]] <-  myfiles_classification[[i]][grep("Eukaryota;",myfiles_classification[[i]]$X2),]  
  
  classification_total <- rbind(classification_total , myfiles_classification[[i]] )
  
}



classification_total$X1 <- gsub(".fa","",classification_total$X1)

# 
# classification_total <-  merge(checkm_total,classification_total,by.x="Name",by.y="X1",all.y=TRUE)
# 
# 
# 
# 
# unit <- as.data.frame(str_split_fixed(classification_total$X2,";",7))
# 
# levels(as.factor(classification_total$phylum))
# 
# 
# colnames(unit) <- c("Domain", "phylum","class","order", "family","genus","species" )
# classification_total <- cbind(classification_total,unit)
# 
# 
# 
# classification_total <- classification_total[which(classification_total$phylum %in% c("Basidiomycota",	
#                                                            "Ascomycota" ,
#                                                            "Euglenozoa",
#                                                            "Microsporidia",
#                                                            "Zoopagomycota",
#                                                            "Mucoromycota",
#                                                            "Oomycota",
#                                                            "Tubulinea")),]
#                                                            
# 
# View(
# classification_total[
# which(
# 
# classification_total$Completeness>40),])

#######################################################################################################
# pick a represnetative and see how weel busco aligns to it 
#######################################################################################################
species="Kluyveromyces marxianus" 

dataset <- "initial_global_bin_metawrap_unmerge"
t <- c()
representative_dataset <- c()
for (dataset in levels(as.factor(classification_total$type.x))){
  
  t <- 
  classification_total[which(classification_total$type.x==dataset),]
  
  for (species in levels(as.factor(t$species))){
   
    t1 <-t[which(t$species==species &
                   t$Completeness==max(t$Completeness[which(t$species==species)])[1]),]
   
    representative_dataset <- rbind(t1,representative_dataset)
  }
  
  
}






  #write.csv(classification_total,"Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/busco/representative_mags.csv",quote = FALSE,row.names = FALSE )                                                         
                                                           
                                             

names_classification  <- gsub( "./|bin_taxonomy.tab|02", "",temp_classification )


names(myfiles_classification ) <- names_classification



#######################################################################################################
# Importing busco data 
#######################################################################################################



setwd("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/busco/")
temp_busco = list.files(pattern="*txt", recursive = TRUE, full.names = TRUE)

temp_busco <- as.data.frame(temp_busco)
temp_busco$merge_id <- gsub(".*\\/|.txt","",temp_busco$temp_busco)
  


temp_busco_rep <-
temp_busco[which(
temp_busco$merge_id %in%representative_dataset$Name ),]
  
  

myfiles_busco = lapply(temp_busco_rep$temp_busco,read_delim,
                                delim = "\t", escape_double = FALSE,
                                col_names = FALSE, trim_ws = TRUE)




#Interpretation of Results
# Complete (C): This indicates the presence of full-length single-copy orthologs. Higher percentages of complete BUSCOs suggest better genome completeness.
# Single-Copy (S): Single-copy orthologs that are not duplicated.
# Duplicated (D): Orthologs found more than once, indicating potential gene duplication or assembly issues.
# Fragmented (F): Partial orthologs, which can suggest incomplete assembly or fragmented genes.
# Missing (M): Orthologs that are not found, indicating gaps in the assembly.
# 5. Thresholds and Benchmarks
# Completeness Threshold: A commonly used benchmark for good completeness in eukaryotic MAGs is having more than 90% of BUSCO genes complete. However, this threshold can vary depending on the specific goals of your study.
# Contamination Threshold: Ideally, you should aim for less than 5% of BUSCO genes duplicated to indicate low contamination and high-quality assembly.

names(myfiles_busco ) <-temp_busco_rep$temp_busco


busco_total_data <- c()

for (i in names(myfiles_busco)){

  busco_total_data <-   
  rbind(busco_total_data,
  
  data.frame(complete=

str_split(
    gsub("\\[|\\]|,", "", myfiles_busco[[i]][7,]),"%",7)[[1]][1], 

duplicated = str_split(
  gsub("\\[|\\]|,", "", myfiles_busco[[i]][7,]),"%",7)[[1]][4] ,

sample=i

)
)
}









busco_total_data$complete <- gsub("C:","",busco_total_data$complete)


busco_total_data$complete <- as.numeric(
  busco_total_data$complete )


busco_total_data$duplicated <- as.numeric(gsub("F:","",busco_total_data$duplicated))

busco_total_data <- 
  
  cbind(busco_total_data,
str_split_fixed(
busco_total_data$sample, "\\/",4)[,c(2,3,4)])



busco_total_data$`3` <- gsub(".txt","",busco_total_data$`3`)

busco_total_data <- merge(
                        busco_total_data,
                                classification_total,
                                                      by.x="3",
                                                                by.y="X1",all.x=TRUE)





busco_total_data_top <- 

busco_total_data[
which(busco_total_data$complete >=90 & busco_total_data$duplicated <=5),]


nrow(busco_total_data_top )

View(busco_total_data_top
     )
  




# representative_dataset_v2 <- representative_dataset
# 
# representative_dataset_v2$busco_total <- NA
# 
# for (i in names(myfiles_busco)){
#   
#   
#   representative_dataset_v2$busco_total[
#   which(representative_dataset_v2$Name==i)] <- 
#   myfiles_busco[[i]][[7,"X1"]]
# }






str_split(representative_dataset_v2$busco_total,":",7)

representative_dataset_v2$busco_completeness <- 
gsub("%\\[.*|C:","",representative_dataset_v2$busco_total)

View(representative_dataset_v2[
which(representative_dataset_v2$busco_completeness>50),])

i="WK003-2-G-48h-01C6_S3_metawrap_bins_maxbin2_bins_bin.2"

representative_dataset_v2$busco_total
names(myfiles_busco) <- temp_busco_rep$merge_id



#######################################################################################################
# Importing GTDB-TK data 
#######################################################################################################

# need to look at busco results and see if they are complementary. 
             
             # &
             # t$Contamination<=10),]
setwd("Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification")

temp_gtdbtk = list.files(pattern="gtdbtk.bac120.summary",recursive=TRUE)

myfiles_gtdbtk = lapply(temp_gtdbtk,read_delim, 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

names(myfiles_gtdbtk) <- gsub( "./|gtdbtk.bac120.summary.tsv|02", "",temp_gtdbtk )

gtdbtk_total <-c()



#######################################################################################################
# Mash the metadata together 
#######################################################################################################

metadata_linker <- data.frame(datasource=as.character(names(myfiles_gtdbtk)),
                              metadata_file=as.character(c("kefir4all_metadata","kefir4all_metadata","global_metadata_total","global_metadata_total")))

metadata <- c()

metadata[["kefir4all_metadata"]] <- kefir4all_metadata
metadata[["global_metadata_total"]] <- global_metadata_total




#metadata[["global_metadata_total"]]$Sample <- paste(metadata[["global_metadata_total"]]$Sample,"_",gsub(".*_","",metadata[["global_metadata_total"]]$merge_column),sep="")



# ensure the metadata is compatable
metadata_for_loop <- c()

#maybe I just need to add metadata to the unmerged
#_S[[:digit:]]+.*
i <- "cs_bin_metawrap_unmerge" 
i <- "global_bin_metawrap_merge"



#i <-  "cs_bin_metawrap_unmerge" 
for (i in names(myfiles_gtdbtk)){

  # make user_genome compatible with merge column in metadata
  
  #WK069_2_L_48h_05D2_S81
  
  #B11_S23
  #myfiles_gtdbtk[[i]]$user_genome <- gsub("-","_",myfiles_gtdbtk[[i]]$user_genome)# unhash this is mistake occurs

  if(length(grep("unmerge", i))==1){
    myfiles_gtdbtk[[i]]$original_name <-  myfiles_gtdbtk[[i]]$user_genome
    myfiles_gtdbtk[[i]]$user_genome <- gsub("-","_",myfiles_gtdbtk[[i]]$user_genome)# mistake will be here
    myfiles_gtdbtk[[i]]$Sample <- 
      gsub("_metawrap_bin_reassembly.*|MK_","",myfiles_gtdbtk[[i]]$user_genome)
    

  myfiles_gtdbtk[[i]] <- 
  merge(myfiles_gtdbtk[[i]], metadata[[metadata_linker$metadata_file[which(metadata_linker$datasource==i)]]],by.x="Sample",by.y= "merge_column",all.x=TRUE)
 
  myfiles_gtdbtk[[i]]$Sample <- 
    gsub("_metawrap_bin_reassembly.*|MK_","",myfiles_gtdbtk[[i]]$user_genome)
  
   myfiles_gtdbtk[[i]]$type <- i
  
  
  myfiles_gtdbtk[[i]]$base_name <-  myfiles_gtdbtk[[i]]$Sample
  }else{
    myfiles_gtdbtk[[i]]$original_name <-  myfiles_gtdbtk[[i]]$user_genome
    myfiles_gtdbtk[[i]]$Sample <-     gsub("_T_","_T.",myfiles_gtdbtk[[i]]$user_genome) 
    myfiles_gtdbtk[[i]]$Sample <-     gsub("IRL_","IRL-", myfiles_gtdbtk[[i]]$Sample) 
    myfiles_gtdbtk[[i]]$Sample  <-    gsub("_metawrap_bin_reassembly.*|MK_","",myfiles_gtdbtk[[i]]$Sample)
    


    #myfiles_gtdbtk[[i]]$Sample <- 
    #str_split_fixed(myfiles_gtdbtk[[i]]$Sample,"_",4) [,3]
    
    #for (type in levels(as.factor(metadata[[metadata_linker$metadata_file[which(metadata_linker$datasource==i)]]]$category))){
      
      #myfiles_gtdbtk[[i]]$Sample <-  
      #str_split_fixed(gsub("TG_","",myfiles_gtdbtk[[i]]$user_genome),"_",4)[,2]
      
      
      
      
      myfiles_gtdbtk[[i]]$`kefir type` <- "merged"
      
      if(i=="global_bin_metawrap_merge"){
        myfiles_gtdbtk[[i]]$Stage<-    "T1"
        myfiles_gtdbtk[[i]]$Sample <-     gsub(".*_","",myfiles_gtdbtk[[i]]$Sample) 
        myfiles_gtdbtk[[i]]$data_source<- gsub("MK","Walsh et al",
                                               gsub("WK", "Samuel et al", 
                                                    gsub("[[:digit:]]+.*","",str_split_fixed(myfiles_gtdbtk[[i]]$user_genome,"_",4)[,1])))
      }else{

        myfiles_gtdbtk[[i]]$Stage<-    str_split_fixed(gsub("TG_|_merged","",myfiles_gtdbtk[[i]]$user_genome),"_",4)[,3]
        myfiles_gtdbtk[[i]]$data_source<- "This study"
      }
      
      myfiles_gtdbtk[[i]]$type <- i
      
      myfiles_gtdbtk[[i]]$base_name <-  myfiles_gtdbtk[[i]]$Sample
      myfiles_gtdbtk[[i]]$Timepoint<- "merged"
      myfiles_gtdbtk[[i]]$country <-      myfiles_gtdbtk[[i]]$Sample
      

      
  
     
        
       
    
      myfiles_gtdbtk[[i]]$Sample.y <-  myfiles_gtdbtk[[i]]$Sample
      
      myfiles_gtdbtk[[i]] <-  myfiles_gtdbtk[[i]] %>% 
        relocate(Sample.y,.after=warnings)
      
      myfiles_gtdbtk[[i]] <- myfiles_gtdbtk[[i]] %>% 
        relocate(Sample,.before=user_genome)
      
      myfiles_gtdbtk[[i]] <- myfiles_gtdbtk[[i]] %>% 
        relocate(Sample,.before=user_genome)
      
      myfiles_gtdbtk[[i]]$category <- gsub("MK","Milk.kefir",
                                           gsub("WK", "Water.kefir", 
                                                gsub("[[:digit:]]+.*","",str_split_fixed(myfiles_gtdbtk[[i]]$user_genome,"_",4)[,1])))
      
      

      
  }
  
  
  
}






# make sure the metadata columns align between the global and cs datasets 
myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]]$Timepoint <- "24x1"
myfiles_gtdbtk[["cs_bin_metawrap_merge"]]$Timepoint <- "24x1"

myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]]$country <- "IRL-cs"#


myfiles_gtdbtk[["cs_bin_metawrap_merge"]]$country <- "IRL-cs"


myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]]$category <- NA
myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]]$category[which(myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]]$`kefir type` %in% c("WL","WG"))] <- "Water.kefir"
myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]]$category[which(myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]]$`kefir type` %in% c("ML","MG"))] <- "Milk.kefir"



myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]] <- myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]] %>% 
  relocate(data_source,.before = category)



myfiles_gtdbtk[["global_bin_metawrap_merge"]]$country[grep("metawrap",myfiles_gtdbtk[["global_bin_metawrap_merge"]]$country)] <- "NA"

#View(myfiles_gtdbtk[["global_bin_metawrap_merge"]])

#View(myfiles_gtdbtk[["cs_bin_metawrap_unmerge"]])

# work from here
  
  gtdbtk_total <-  
    
    rbind(myfiles_gtdbtk[[1]],
          myfiles_gtdbtk[[2]],
          myfiles_gtdbtk[[3]],
          myfiles_gtdbtk[[4]]
    )
  

#######################################################################################################
#merge with checkm data 
#######################################################################################################
  
 # View(as.data.frame(levels(as.factor(checkm_total$Name[which(checkm_total$type== "global_bin_metawrap_unmerge")])))
checkm_total$Name[which(checkm_total$type== "global_bin_metawrap_unmerge")] <- gsub("-1-","_1_",  checkm_total$Name[which(checkm_total$type== "global_bin_metawrap_unmerge")])

checkm_total$Name[which(checkm_total$type== "global_bin_metawrap_unmerge")] <- gsub("-2-","_2_",  checkm_total$Name[which(checkm_total$type== "global_bin_metawrap_unmerge")])
  
gtdbtk_total <- 
merge(gtdbtk_total, checkm_total,by.x="user_genome",by.y="Name",all.x=TRUE)
  
gtdbtk_total$data_type <- gsub("_bin_metawrap.*","",gtdbtk_total$type.x)


# now its time to add one final piece of metadata that will link merged and unmerged
gtdbtk_total$linker <- gtdbtk_total$country



gtdbtk_total$linker[which( gtdbtk_total$linker=="IRL-cs" )] <- NA


gtdbtk_total$linker[which(is.na(gtdbtk_total$linker )
)] <- gtdbtk_total$Sample[which(is.na(gtdbtk_total$linker ))]

gtdbtk_total$linker <- gsub("_S[[:digit:]]+.*|_merged|","",gtdbtk_total$linker)
                            
gtdbtk_total$linker <- gsub("_WL_|_ML_|_WG_|_MG_","_",gtdbtk_total$linker)

gtdbtk_total$linker <- gsub("_1_L.*|_2_L.*|_1_G.*|_2_G.*","",gtdbtk_total$linker)
gtdbtk_total$linker <- gsub("WK_","",gtdbtk_total$linker)





#write.csv(gtdbtk_total,file.path(DATA_DIR, "mag_metadata_v3.csv.gz"),row.names = FALSE)
#write.table(gtdbtk_total,"Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/mag_metadata_v3.txt",quote = FALSE,row.names = FALSE)


#######################################################################################################
# compare genome numbers and quality across datasets and identify represntative genomes of both types
#######################################################################################################
  
type <- "global"
category <- "global_bin_metawrap_merge"  

t <- c()
link <- "TG_ID113_T6"


number_of_species_recovered <- data.frame(linker=as.character(),
                                          type=as.character(),
                                          n_species=as.numeric())

# this section gets te number of species detected and the representative genomes across both merged and unmerged datasets

t1 <- c()
s <- c()
t2 <- c()
t3 <- c()
gtdbtk_representatives <-c()
for (link in levels(as.factor(gtdbtk_total$linker))){
  
  if(link %in% c( "A4","A5","H3","T10", "T9" ,"U_2_10","U_2_11","U_2_12", "U_2_7","U_2_8", "U_2_9", "U10", "U11" , "U12","U9","Y3"    )){
    next
  }
  
  t <- gtdbtk_total[which(gtdbtk_total$linker==link #&
                          #gtdbtk_total$type.x==category
  ),]
  
  t1 <-   as.data.frame(table(dplyr::select(t,type.x,classification)))
  
  if(0 %in% t1$Freq){
  t1 <- t1[-c(which(t1$Freq==0)),]
  }
  
  for (type in levels(as.factor(gtdbtk_total$type.x[which(gtdbtk_total$linker==link )]))){
                                

    
 s <-    data.frame(linker=as.character(link),
               type=as.character(type),
               n_species=as.numeric(length(t1$type.x[which(t1$type.x==type)])))
    
    
    
    number_of_species_recovered <- rbind(s,
                                         number_of_species_recovered)
    

  
  #quality of genomes recovered
 
   for (species in levels(as.factor(t1$classification))){
    
 t2 <-    t[which(t$classification==species &
                     t$type.x==type),]
    
 

 
 t2$Completeness_rank <- rank(-t2$Completeness)
 t2$Contamination_rank <- rank(t2$Contamination) # lower contamination is advantageous hence no - sign
 
 t2$rank_score <-  t2$Completeness_rank+   t2$Contamination_rank 
 
 t2 <- t2[which(t2$rank_score==   min(t2$rank_score)),]
 
   if(nrow(t2)>1){ # in cases where the ranking is equal priortise the one with the greater change
     
     net_change_completeness <- t2$Completeness[1]-t2$Completeness[2]
     net_change_contamination <- t2$Contamination[1]-t2$Contamination[2]
     
     if(net_change_completeness>net_change_contamination){
       
       t2 <-  t2[which(t2$Completeness_rank==min(t2$Completeness_rank)),] 
     }else{
       t2 <-  t2[which(t2$Contamination_rank ==min(t2$Contamination_rank )),] 
     }
     
   }
   #t2 <-  t2[which(t2$Completeness_rank==min(t2$Completeness_rank)),]
 
 
 if(nrow(t2)>1){ # in cases where the they are basically the same MAG
   
   t2 <-  t2[1,]
   
 }
 # save you represnetatives genomes of both merged and unmerged in one dataframe
 gtdbtk_representatives <- rbind(gtdbtk_representatives,t2)
  }
 
    }

      }




gtdbtk_representatives$classification <- gsub(".*s__","",gtdbtk_representatives$classification)



# this section tell which of the reference genomes the merged or unmerged is superior 

link <- "FR-A"


mag_quaility_comparision <- data.frame(linker=as.character(),
                                       species=as.character(),
                                       merged=as.numeric(),
                                       unmerged=as.numeric()
                                       )


gtdbtk_representatives$merge_type <- gsub(".*_","",gtdbtk_representatives$type.x)



gtdbtk_representatives$path <- NA
  
gtdbtk_representatives$path[which(gtdbtk_representatives$type.x=="cs_bin_metawrap_merge" )] <- "/data/Food/analysis/R6564_NGS/kefir4all/02_assembly/02_metawrap/02_metawrap_merged/finished_cs_bins/"

gtdbtk_representatives$path[which(gtdbtk_representatives$type.x=="cs_bin_metawrap_unmerge" )] <- "/data/Food/analysis/R6564_NGS/kefir4all/02_assembly/02_metawrap/02_metawrap_unmerged/finished_cs_bins/"
  
gtdbtk_representatives$path[which(gtdbtk_representatives$type.x=="global_bin_metawrap_merge" )] <- "/data/Food/analysis/R6564_NGS/kefir4all/02_assembly/02_metawrap/02_metawrap_merged/finished_global_bins/"

gtdbtk_representatives$path[which(gtdbtk_representatives$type.x=="global_bin_metawrap_unmerge" )] <- "/data/Food/analysis/R6564_NGS/kefir4all/02_assembly/02_metawrap/02_metawrap_unmerged/finished_global_bins/"
  
gtdbtk_representatives$path <- paste(gtdbtk_representatives$path,gtdbtk_representatives$original_name,".fa",sep="")
levels(as.factor(gtdbtk_representatives$type.x))



gtdbtk_representatives_high_quaility <- 
  gtdbtk_representatives[which(gtdbtk_representatives$Completeness >80 &
                                 gtdbtk_representatives$Contamination<5),]


gtdbtk_representatives$classification <- gsub("Lactococcus lactis_E","Lactococcus cremoris",gtdbtk_representatives$classification)

gtdbtk_representatives$classification <-  gsub("_.","",gtdbtk_representatives$classification)


n_mags <- 
as.data.frame(table(gtdbtk_representatives$classification))


gtdbtk_representatives_high_quaility$classification <- gsub("Lactococcus lactis_E","Lactococcus cremoris",gtdbtk_representatives_high_quaility$classification)

gtdbtk_representatives_high_quaility$classification <-  gsub("_.","",gtdbtk_representatives_high_quaility$classification)


n_mags_high_quaility <- 
  as.data.frame(table(gtdbtk_representatives_high_quaility$classification))



#write.csv(n_mags[-c(
  #which(n_mags$Freq<10)),], "Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/high_quaility_mags_detected.csv",row.names=FALSE)





                                              

#write.csv(gtdbtk_representatives[which(gtdbtk_representatives$classification %in% 
                                         #n_mags$Var1[-c(which(n_mags$Freq<10))]),],
          #"Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/gtdbtk_representatives_for_roary.csv",row.names=FALSE)



#species <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Liquorilactobacillus;s__Liquorilactobacillus ghanensis" 



link <- "TG_ID092_T3"

net_change_completeness <- c()
net_change_contamination <- c()

mag_quaility_comparision <- data.frame(linker=as.character(),
                                       species=as.character(),
                                       merged=as.numeric(),
                                       unmerged=as.numeric()
)


best_gtdbtk_representatives <- c()

for (link in levels(as.factor( gtdbtk_representatives$linker))){
  
  t <- 
  gtdbtk_representatives[which(gtdbtk_representatives$linker==link),]
  
  
  for (species in levels(as.factor(t$classification))){
    t2 <- t[which(t$classification==species),]
    
    
    
    t2$Completeness_rank <- rank(-t2$Completeness)
    t2$Contamination_rank <- rank(t2$Contamination) # lower contamination is advantageous hence no - sign
    
    t2$rank_score <-  t2$Completeness_rank+   t2$Contamination_rank 
    
    t2 <- t2[which(t2$rank_score==   min(t2$rank_score)),]
    

    if(nrow(t2)>1){ # in cases where the ranking is equal priortise the one with the greater change
      
      net_change_completeness <- t2$Completeness[1]-t2$Completeness[2]
      net_change_contamination <- t2$Contamination[1]-t2$Contamination[2]
      
      if(net_change_completeness>net_change_contamination){
       
        t2 <-  t2[which(t2$Completeness_rank==min(t2$Completeness_rank)),] 
      }else{
        t2 <-  t2[which(t2$Contamination_rank ==min(t2$Contamination_rank )),] 
      }
      
    }
    
    if(nrow(t2)>1){ # in cases where the they are basically the same MAG
      
      t2 <-  t2[1,]
      
    }
    
    
    
    # now that the represnetative is sorted figure out where it can from and record this in a seperate dataframe
    if(t2$merge_type=="unmerge"){
      mag_quaility_comparision <- 
        rbind(mag_quaility_comparision,
              data.frame(linker=as.character(link),
                         species=as.character(species),
                         type=as.character("unmerge"),
                         value=as.numeric(1)),
              
              data.frame(linker=as.character(link),
                         species=as.character(species),
                         type=as.character("merge"),
                         value=as.numeric(0)
    )
      )
    }else if (t2$merge_type=="merge"){
      mag_quaility_comparision <- 
        
        rbind(mag_quaility_comparision,
              data.frame(linker=as.character(link),
                         species=as.character(species),
                         type=as.character("unmerge"),
                         value=as.numeric(0)),
              
              data.frame(linker=as.character(link),
                         species=as.character(species),
                         type=as.character("merge"),
                         value=as.numeric(1)
      )
      )
    }else{
      print(paste("something went wrong with", species,link,sep=" "))
    }
    
    best_gtdbtk_representatives <- rbind(t2,best_gtdbtk_representatives)
}

  
  
  
  
}

          
  # write.table(    best_gtdbtk_representatives, "Q:/H2020 Master/Citizen Science Project/Results/03_mag_classification/per_grain_gtdbtk_representatives.txt",row.names=FALSE,quote=FALSE,sep="\t")       
    
          

#######################################################################################################
# plot time
#######################################################################################################


#
# 
#View()
mag_quaility_comparision$species_v2 <-   gsub(".*;s__","",mag_quaility_comparision$species)

mag_quaility_comparision$species_v2[
which(mag_quaility_comparision$species_v2 =="")] <- paste( "Putatively novel -", gsub(".*;g__|;s__","",
                                                                                mag_quaility_comparision$species[
                                                                                  which(mag_quaility_comparision$species_v2 =="")]))



t1 <- 
dplyr::select(gtdbtk_representatives,linker, data_source, category)

t1$data_source <- gsub(" 2023","",t1$data_source)


t1 <- 
t1[-c(which(duplicated(t1))),]

mag_quaility_comparision <- merge(
mag_quaility_comparision,
t1,by="linker",all.x=TRUE)

mag_quaility_comparision <- 

mag_quaility_comparision[-c(which(is.na(mag_quaility_comparision$category))),]

mag_quaility_comparision$type <- str_to_title(mag_quaility_comparision$type)


library(forcats)

 pacman::p_load(tidytext)
 
p1 <- 
ggplot(data=mag_quaility_comparision[which(
 mag_quaility_comparision$value!=0 &
   mag_quaility_comparision$category=="Milk.kefir"
   ),], aes(x= fct_infreq(species_v2), fill=type))+
  #mag_quaility_comparision$value!=0),], aes(x= species_v2, fill=type))+
  geom_bar(state="identity")+
  #labs(title = "8hours")+
  #+
  facet_wrap(.~as.factor(category),scales = "free_x")+
  scale_x_reordered()+
  theme_bw()+
  theme(#legend.position = "none",
    axis.text.x = element_text(size = 11.5, colour = "black",angle = 90,hjust=1,face="italic"),
    
    axis.title.y = element_text(size = 20, face = "bold"),
  #  legend.title = element_text(size = 16, face = "bold"), 
   # legend.text = element_text(size = 8, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black" ,vjust = 0.5),
    axis.ticks.y =element_blank(),
    legend.text=element_text(size = 20),
    legend.key.size = unit(2, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm'), #change legend key width
    legend.title = element_text(size=20),
    strip.background = element_rect(
      color="black", fill="white"),
    strip.text = element_text(size=15.5))+
  #   plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
  #                                               box.color = "#748696",size=55, lineheight = 2))+ 
  # scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Sample", y = "Number of species sample specific representative MAGs", fill="Assembly type")


p2 <- 
  ggplot(data=mag_quaility_comparision[which(
    mag_quaility_comparision$value!=0 &
      mag_quaility_comparision$category=="Water.kefir"
  ),], aes(x= fct_infreq(species_v2), fill=type))+
  #mag_quaility_comparision$value!=0),], aes(x= species_v2, fill=type))+
  geom_bar(state="identity")+
  #labs(title = "8hours")+
  #+
  facet_wrap(.~as.factor(category),scales = "free_x")+

  theme_bw()+
  theme(#legend.position = "none",
    axis.text.x = element_text(size = 10.5, colour = "black",angle = 90,hjust=1,face="italic"),
    
    axis.title.y = element_text(size = 20, face = "bold"),
    #  legend.title = element_text(size = 16, face = "bold"), 
    # legend.text = element_text(size = 8, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black" ,vjust = 0.5),
    axis.ticks.y =element_blank(),
    legend.text=element_text(size = 20),
    legend.key.size = unit(2, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm'), #change legend key width
    legend.title = element_text(size=20),
    strip.background = element_rect(
      color="black", fill="white"),
    strip.text = element_text(size=15.5))+
  #   plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
  #                                               box.color = "#748696",size=55, lineheight = 2))+ 
  # scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Sample", y = "Number of MAGs", fill="Assembly type")


number_of_species_recovered <- 
merge(number_of_species_recovered, t1,by="linker",all.x=TRUE)

number_of_species_recovered <- 
number_of_species_recovered[-c(which(is.na(number_of_species_recovered$category))),]

number_of_species_recovered$type_v2 <-    str_to_title(gsub(".*_","",number_of_species_recovered$type))


table(
mag_quaility_comparision[which(
  mag_quaility_comparision$value!=0 &
    mag_quaility_comparision$category=="Milk.kefir"
),"species"])



table(
  mag_quaility_comparision[which(
    mag_quaility_comparision$value!=0 &
      mag_quaility_comparision$category=="Water.kefir"
  ),"species"])


# ggplot(data=number_of_species_recovered, aes(x= fct_reorder(linker,n_species),y=n_species, fill=type_v2))+
#   #mag_quaility_comparision$value!=0),], aes(x= species_v2, fill=type))+
#   geom_col()+
#   #labs(title = "8hours")+
#   #+
#   facet_wrap(.~as.factor(category),scales = "free_x")+
#   scale_x_reordered()+
#   theme_bw()+
#   theme(#legend.position = "none",
#     axis.text.x = element_blank(),
#     
#     axis.title.y = element_text(size = 20, face = "bold"),
#     #  legend.title = element_text(size = 16, face = "bold"), 
#     # legend.text = element_text(size = 8, face = "bold", colour = "black"),
#     axis.text.y = element_text(size = 10, colour = "black" ,vjust = 0.5),
#     axis.ticks.y =element_blank(),
#     legend.text=element_text(size = 20),
#     legend.key.size = unit(2, 'cm'), #change legend key size
#     legend.key.height = unit(2, 'cm'), #change legend key height
#     legend.key.width = unit(2, 'cm'), #change legend key width
#     legend.title = element_text(size=20),
#     strip.background = element_rect(
#       color="black", fill="white"),
#     strip.text = element_text(size=15.5))+
#   #   plot.title = ggtext::element_textbox_simple(face="bold",halign  = 0.5,linetype = 1, # turn on border
#   #                                               box.color = "#748696",size=55, lineheight = 2))+ 
#   # scale_y_continuous(expand = c(0,0)) + 
#   labs(x = "Sample", y = "Number of MAGs", fill="Assembly type")





ggplot(data=number_of_species_recovered, aes(x= type,y=n_species, fill=type_v2))+
  #mag_quaility_comparision$value!=0),], aes(x= species_v2, fill=type))+
  geom_boxplot()+
  #labs(title = "8hours")+
  #+
  facet_wrap(.~as.factor(category),scales = "free_x")




# jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/Evolution/Figure 5_MAG_comparisons.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
# 
# 
#ggarrange( bubble_plot,
#            ggarrange(p1,p2, nrow=1,ncol=2, labels=c("B.","C."),common.legend = TRUE,  font.label = list(size = 30)),
#            nrow=2,ncol=1, heights = c(5, 6),widths=c(10,1),labels=c("A.",""),  font.label = list(size = 30))
# graphics.off()

#######################################################################################################
# Select the reference genomes for downstream anlaysis 
#######################################################################################################


mag_breakdown <- as.data.frame(
table(dplyr::select(
mag_quaility_comparision, species_v2, data_source, category,type)))


mag_breakdown <- 
mag_breakdown[-c(which(mag_breakdown$Freq ==0)),]

t <- mag_quaility_comparision[which(
  mag_quaility_comparision$value!=0),]

length(t$species_v2[which(
t$species_v2=="Acetobacter fabarum")])






length(mag_quaility_comparision$species_v2[which(
  mag_quaility_comparision$value!=0)])

length()

mag_quaility_comparision 
number_of_species_recovered 

#######################################################################################################
# modify gtdbtk data 
#######################################################################################################
# make id column more informative and comparable 

gtdbtk_total$merge_column <- 
gsub("_metawrap_bin_reassembly|reassembled_|original_","",gtdbtk_total$user_genome)
     

gtdbtk_total$merge_column <- 
gsub(".orig|.strict|.permissive","",gtdbtk_total$merge_column)
     
     
     




# total_metadata 
# 
# 
# levels(as.factor(gtdbtk_total$type.x))



# 
# s1 <- as.data.frame(str_split_fixed(gtdbtk_total$user_genome[which(gtdbtk_total$type.x=="cs_bin_metawrap_merge")] ,"_",7))
# 
# s1 <- 
# merge(s1, total_metadata,by.x="V3",by.y="Sample",all.x=TRUE)
# 
# s1$merge_column <-NA
#   
# s1$merge_column[which( s1$Stage=="T0")] <-   paste(s1$V1,s1$V3,s1$V5,s1$V4,S1$,sep="_")
# 
# s1$merge_column[which( s1$Stage!="T0")] <-   paste(s1$V1,s1$V7,sep="_")
# 
# 
# 
# s <- as.data.frame(
# str_split_fixed(gsub("TG_","",gtdbtk_total$user_genome) ,"_",5))
# 
# 
# 
# str_split_fixed("ID087_WL_T3_S171_reassembled_bin.1.orig","_",4)
# 
# 
# view(as.data.frame(gtdbtk_total$user_genome))
# 
# s$V5 <- gsub("reassembled_|original_|.orig|.strict|.permissive","",s$V4)
# 
# s$V6 <- 
# paste(s$V2,s$V3,s$V5,sep="_")
# 
# 
# s <-
# dplyr::select(s, V1,V6,V4)
# 
# 
# gtdbtk_total <- 
#   cbind(s,gtdbtk_total) 

#to ensure all mags are captured
#gtdbtk_total$user_genome[-c(which(t$user_genome %in% gtdbtk_total$user_genome))]



gtdbtk_total <- gtdbtk_total[which(gtdbtk_total$Completeness>=50 &
                          gtdbtk_total$Contamination<=10),]



# Select one MAG from the same ID with the best stats
i <- "A11_S11_bin.1"

t <- c()

#gtdbtk_representatives <- c()
for (i in levels(as.factor(gtdbtk_total$merge_column))){
 
   t <- 
  
  gtdbtk_total[which(gtdbtk_total$merge_column==i),]
   
   t$Completeness_rank <- rank(-t$Completeness)
   t$Contamination_rank <- rank(t$Contamination) # lower contamination is advantageous hence no - sign
   
   t$rank_score <-  t$Completeness_rank+   t$Contamination_rank 
   
   t <- t[which(t$rank_score==   min(t$rank_score)),]
   
   if(nrow(t)>1){ # in cases where the ranking is equal priortise genome completeness 
     
    t <-  t[which(t$Completeness_rank==min(t$Completeness_rank)),]
   }
   
   if(nrow(t)>1){ # in cases where the they are basically the same MAG
     
     t <-  t[1,]
    
   }

   gtdbtk_representatives <- rbind(gtdbtk_representatives,t)
}


gtdbtk_representatives <- 
gtdbtk_representatives %>% 
  relocate(merge_column,.before = user_genome)




#######################################################################################################
# Improve the current metadata

#######################################################################################################



# as we don't have a metadata file for the merged data we need to outline which are milk kefir and which are water kefir, not necessary with the unmerged group was the metadata exists
gtdbtk_representatives$kefir_category <- NA




gtdbtk_representatives$kefir_category[which(gtdbtk_representatives$type.x %in% c("global_bin_metawrap_merge", "cs_bin_metawrap_merge"))][
  grep("MK",
       gtdbtk_representatives$merge_column[which(gtdbtk_representatives$type.x %in% c("global_bin_metawrap_merge", "cs_bin_metawrap_merge"))])] <- "Milk.Kefir"




gtdbtk_representatives$kefir_category[which(gtdbtk_representatives$type.x %in% c("global_bin_metawrap_merge", "cs_bin_metawrap_merge"))][
grep("WK",
gtdbtk_representatives$merge_column[which(gtdbtk_representatives$type.x %in% c("global_bin_metawrap_merge", "cs_bin_metawrap_merge"))])] <- "Water.Kefir"


View(gtdbtk_representatives)


gtdbtk_representatives$merge_column[which(gtdbtk_representatives$kefir_category=="Milk.Kefir" &
                                      gtdbtk_representatives$type.x=="global_bin_metawrap_merge")]




gtdbtk_representatives$kefir_category[grep("WK",gtdbtk_representatives$merge_column)] <- 


gtdbtk_representatives$merge_column[which(gtdbtk_representatives$type.x=="global_bin_metawrap_merge")]









gtdbtk_representatives$merge_column[which(gtdbtk_representatives$type.x=="cs_bin_metawrap_merge")]
  
  
gtdbtk_representatives$merge_column[which(gtdbtk_representatives$type.x=="global_bin_metawrap_merge")]

  
  


gtdbtk_representatives$assembly_type <- "Unmerged"
gtdbtk_representatives$assembly_type[grep("_merged_",gtdbtk_representatives$merge_column)] <- "Merged"

gtdbtk_representatives$base_name <- gsub("MK_|WK_|_bins.*|_merged_.*","",gtdbtk_representatives$merge_column)


t <- 
merge(total_metadata,gtdbtk_representatives,by.x="merge_column",by.y="base_name",all.y=TRUE)
#######################################################################################################
# plot mags recovered across categories

#######################################################################################################
#which(duplicated()

levels(as.factor(gtdbtk_representatives$classification))







#View(   gtdbtk_representatives )
gtdbtk_representatives$species <- 
gsub(".*;s__","",gtdbtk_representatives$classification)





gtdbtk_representatives$kefir_type <- "Milk.kefir"
gtdbtk_representatives$kefir_type [grep("WK",gtdbtk_representatives$merge_column)]


gtdbtk_representatives$kefir_type [grep("WK|WL|WG",gtdbtk_representatives$merge_column)] <- "Water.kefir"


ggplot(gtdbtk_representatives,aes(x=species,fill=type.x))+
  geom_bar()+
  facet_wrap(~kefir_type,scales="free" )+
  # geom_boxplot()+
  #geom_point()+
  labs(x="Species", y="strains comparisons", title="",fill="Evolutionary change") +
  #coord_equal() +
  theme_bw()+
  theme(legend.position = "right",#axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.y = element_text(size=15,face="italic"),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=7.5,face = "italic",angle=90,hjust = .5,vjust = .6),
        #axis.text.x = element_text(size=10,angle = 45,hjust = 1), # remove x-axis labels
        #axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.text=element_text(face="italic",size = 20),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=20),
        strip.background = element_rect(
          color="black", fill="white"),
        strip.text = element_text(size=15.5,face="italic"))
  
  




n_mags_data <- as.data.frame(
  #species_without_change_participant
  table(dplyr::select(       gtdbtk_representatives  ,species,type.x,)))




p_1 <- 
  n_participants_with_evolution %>% 
  mutate(species=gsub("_"," ",species)) %>% 
  mutate(change=gsub("FALSE", "Evolutionary change",
                     gsub("TRUE","No evolutionary change",change ))) %>% 
  ggplot(.,aes(x=reorder(species,Freq),y=Freq,fill=change))+
  geom_col()+
  # geom_boxplot()+
  #geom_point()+
  labs(x="Species", y="strains comparisons", title="",fill="Evolutionary change") +
  #coord_equal() +
  theme_bw()+
  theme(legend.position = "right",#axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.y = element_text(size=15,face="italic"),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=15,face = "italic",angle=90,hjust = .5,vjust = .6),
        #axis.text.x = element_text(size=10,angle = 45,hjust = 1), # remove x-axis labels
        #axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.text=element_text(face="italic",size = 20),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=20),
        strip.background = element_rect(
          color="black", fill="white"),
        strip.text = element_text(size=15.5,face="italic"))









View(as.data.frame(
gtdbtk_total$user_genome))

#names_gtdbtk <- as.data.frame(names_gtdbtk)

#######################################################################################################

sample_name.1_gtdbtk <- c()
sample_name_gtdbtk <- c()
bin_name_gtdbtk <- c()
bin_name.1_gtdbtk <- c()
classification.1 <- c()
classification <- c()
novelty <- c()
fastani <- c()

gtdbtk_data <-data.frame(sample_name_gtdbtk=character(),
                         bin_name_gtdbtk=character(),
                         classification=character(),
                         novelty=character(),
                         fastani=character()
)

for (i in 1:length(myfiles_gtdbtk)){
  length_gtdbtk <- length(getElement(myfiles_gtdbtk[[i]], "user_genome"))
  sample_name_gtdbtk <- replicate(length_gtdbtk,names_gtdbtk[[i]])
  
  bin_name_gtdbtk <- getElement(myfiles_gtdbtk[[i]], "user_genome")
  
  
  classification <- getElement(myfiles_gtdbtk[[i]], "classification")
  
  fastani <- getElement(myfiles_gtdbtk[[i]], "fastani_ani")
  
  novelty <- getElement(myfiles_gtdbtk[[i]], "note")
  
  data.1 <- cbind(sample_name_gtdbtk,bin_name_gtdbtk,classification,fastani,novelty )
  
  gtdbtk_data  <- merge(gtdbtk_data ,data.1,all=TRUE )
  
}


gtdbtk_data <-as.data.frame(gtdbtk_data)




gtdbtk_data$species <- gtdbtk_data$classification
gtdbtk_data$genus <- gtdbtk_data$classification
gtdbtk_data$species <- gsub(".*s__","",gtdbtk_data$species)
gtdbtk_data$genus <- gsub(".*g__","",gtdbtk_data$genus)

gtdbtk_data$genus <- gsub(";.*","",gtdbtk_data$genus)

#names_checkm <- gsub( "checkm_quality_tab_assessment.*", "",names_checkm$names_checkm)



#########################################################################################################

#remove controls 
gtdbtk_data$sample_name_gtdbtk <- gsub("_.*","",gtdbtk_data$sample_name_gtdbtk)
controls_to_remove <- c("A1","A2","A3","A4","A5","T6","T7","T8","T9","T10","T11", "T12", "U1","U2","U3","U4","U5","U6","U9","U10","U11","U12", "W11","W12","Y1","Y2","Y3","Y4")
c.1 <-c() 
c.2 <- c()
for (i in controls_to_remove){
  c.1 <-  which(gtdbtk_data$sample_name_gtdbtk==i)
  c.2 <- c(c.1,c.2)}

gtdbtk_data <- gtdbtk_data[-c(c.2),]




checkm_data$base<-checkm_data$id
checkm_data$base <- gsub("_S.*","",checkm_data$base)

#View(gtdbtk_data[which(gtdbtk_data$novelty=="N/A"),])
c.3 <-c() 
c.4 <- c()
for (i in controls_to_remove){
  c.3 <-  which(checkm_data$base==i)
  c.4 <- c(c.3,c.4)}

checkm_data <- checkm_data[-c(c.4),]

checkm_data <- checkm_data[-c(grep("H3_S87_merge_host_removed_bin",checkm_data$id)),]

#########################################################################################################
#Getting prevalence info of the number of genomes recovered

frequency_species <- data.frame(x=as.character(),
                                freq=as.numeric())
data.1 <- c()


for (i in levels(as.factor(gtdbtk_data$species))){
  data.1 <-  count(gtdbtk_data$species[which(gtdbtk_data$species==i)])
  frequency_species <- merge(data.1,frequency_species,all=TRUE)
}

#write.csv(frequency_species,"E:/Excel/Metagenomic data/metawrap/species_frequency_data.csv", row.names = FALSE)



#######################################################################################################

gtdbtk_data$base <-  sub('^([^.]+.[^.]+).*', '\\1', gtdbtk_data$bin_name_gtdbtk)

gtdbtk_data$base <- paste(gtdbtk_data$sample_name_gtdbtk,"_",gtdbtk_data$base,sep="")
checkm_data$base <- gsub("_S\\w+\\_merge_host_removed","", checkm_data$id)

##Duplicated H3 work here tomorrow

total_data_full<- merge(gtdbtk_data,checkm_data,by="base")

total_data <- dplyr::select(total_data_full,id, species, genus, Completeness.1, Contamination.1,fastani,novelty,Strain_heterogeneity.1 )


total_data <-  total_data[!duplicated(total_data$id), ]

total_data$genus <- gsub("_E","",total_data$genus)

#levels(as.factor(total_data$genus))
#subset(temp_test, !(temp_test$base %in% temp_gtdbtk_test$base))

#######################################################################################################
# changing the tip labels to show genus level info
#Note this section isn't necessary






#######################################################################################################
# Removing control data 

id <- gsub("_S.*","",id)
controls_to_remove <- c("A1","A2","A3","A4","A5","T6","T7","T8","T9","T10","T11", "T12", "U1","U2","U3","U4","U5","U6","U9","U10","U11","U12", "W11","W12","Y1","Y2","Y3","Y4")
c.1 <-c() 
c.2 <- c()
for (i in controls_to_remove){
  c.1 <-  which(id==i)
  c.2 <- c(c.1,c.2)}
#Drop control tips
snake.tree <- drop.tip(snake.tree,snake.tree$tip.label[c.2])
#Drop contaminated tips
snake.tree <- drop.tip(snake.tree, snake.tree$tip.label[grep("H3_S87_merge_host_removed_bin",snake.tree$tip.label)])

#c.3 <- c()
#c.4 <- c()


#total_data$id <- sub('^([^.]+.[^.]+).*', '\\1', total_data$id)


#for (i in c.2){
#c.3 <- which(total_data$id==i)
#c.4 <- c(c.3,c.4)}

#total_data <-   total_data[-c(c.4),]
#######################################################################################################
pacman::p_load(phylotools,ape)

#merge(data_labels,total_data, by="id")

tip_label <- dplyr::select(total_data,"id", "genus")


tip_label$id <-sub('^([^.]+.[^.]+).*', '\\1', tip_label$id)


ntree <- sub.taxa.label(snake.tree,tip_label)
#Make sure this is correct tomorrow

#View(merge(test, tip_label, by.x="base",by.y="id",all.x=TRUE))



######################################################################################################
#Summary info on MAGS:

length(total_data$id) # Numbers of MAGs recovered




#total_data[which(total_data$Completeness.1>=90, total_data$Contamination.1 <=5),]

High_mags <-subset(x=total_data,
                   subset=as.numeric(Completeness.1) >= 90 &
                     as.numeric(Contamination.1) <= 5)


##############
#Note this metadata file will be used for clustering
# High_mags$id <- paste0(High_mags$id,".fa",sep="")
# 
# High_mags<- dplyr::select(High_mags,id,Completeness.1,Contamination.1,Strain_heterogeneity.1)
# colnames(High_mags) <- c("genome","completeness","contamination","strain_heterogeneity")
# High_mags$species <- gsub("_.*","",High_mags$species)
# write.csv(High_mags,"E:/Excel/Metagenomic data/metawrap/High_mags_reference_genomes.csv",row.names = FALSE,quote = FALSE,sep=",")
############


#Lactococcus_lactis_reference <- High_mags[grep("Lactococcus lactis",High_mags$species),1]
#Lactococcus_lactis_reference <- as.data.frame(Lactococcus_lactis_reference)

#Lactococcus_lactis_reference$bin <- Lactococcus_lactis_reference$Lactococcus_lactis_reference

#Lactococcus_lactis_reference$bin <- gsub(".*_","",Lactococcus_lactis_reference$bin)

#Lactococcus_lactis_reference$Lactococcus_lactis_reference <- gsub("_bin.*","",Lactococcus_lactis_reference$Lactococcus_lactis_reference)


#write.table(Lactococcus_lactis_reference,"E:/Excel/Metagenomic data/gtdbtk/Lactococcus_lactis_reference_genomes.txt",row.names = FALSE,quote = FALSE)



med_mags <- subset(x=total_data,  
                   subset=as.numeric(Completeness.1) >= 80 &
                     as.numeric(Completeness.1) < 90 &
                     as.numeric(Contamination.1) <= 5)





novel_MAGS <- subset(x=High_mags,
                     subset=novelty=="N/A")

novel_MAGS_mid <- subset(x=med_mags,
                         subset=novelty=="N/A")





#get the bins with novel  species 


#write.csv(
#rbind(novel_MAGS,novel_MAGS_mid),"E:/Excel/Metagenomic data/metawrap/putative_novel_MAGS.csv")

combined_med_high_Mags <- rbind(High_mags,med_mags)
'%!in%' <- function(x,y)!('%in%'(x,y)) 
snake.tree <- drop.tip(snake.tree, snake.tree$tip.label[which(snake.tree$tip.label %!in% combined_med_high_Mags$id)])




##############
#Note this section is metadata about the number of species per assembly used in the manuscript
#combined_med_high_Mags$species <- 
# gsub("_E","",combined_med_high_Mags$species)
# combined_med_high_Mags$species <- 
#   gsub("_A","",combined_med_high_Mags$species)
# 
# 
# combined_med_high_Mags$species[combined_med_high_Mags$species==""] <- "Putative novel species"
# 
# write.csv(combined_med_high_Mags,"E:/STORE N GO/Writings/kefir-metagenomics/Tables/Supplementary table 3S.csv",row.names = FALSE )

##############

#####################################################################################################
# genus and species info 
ngenus <- levels(as.factor(combined_med_high_Mags$genus))

# Break down of the number of genuses
genus_breakdown <- data.frame(genus_name=as.character(),
                              length=as.numeric())
length <- c()
genus_name <- c()

for (i in (levels(as.factor(combined_med_high_Mags$genus)))){
  length <-length(which(total_data$genus==i))
  genus_name <- as.character(i)
  g_data <- cbind(length,genus_name)
  genus_breakdown <- merge(g_data,genus_breakdown,all=TRUE)
}


#length(which(total_data$genus=="Acetobacter"))

#####################################################################################################


#scale_color_manual(values=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"))


#setwd("E:/STORE N GO/R/Plots")
#png(filename='ggtree_MAGS_gtdbtk.png', width = 35*700, height=30*700,res=1700,pointsize = 15) #, width=2000, height=1950)

#plot(p)

#graphics.off()



###################################################################################
#Changing the geom_tippoint labels
# You can remove this step if you don't want the genus name fllowed by n and a number

# y <- 1
# for (i in genus_breakdown$genus_name){
#   total_data$genus[which(total_data$genus==i)] <- paste0(genus_breakdown$genus_name[y]," -n(",genus_breakdown$length[y],")")
#   y <- y+1
# }

#######################################################################################################
combined_med_high_Mags$putative <- combined_med_high_Mags$novelty

na_location <- grep("N/A", combined_med_high_Mags$putative )

combined_med_high_Mags$putative[-c(na_location)] <- "Existing species"
combined_med_high_Mags$putative[na_location] <- "Putative novel species"


#######################################################################################################
#total_data$id <- sub('^([^.]+.[^.]+).*', '\\1',total_data$id )

#Plotting section


View(combined_med_high_Mags)

p <- c()
genus <- c("Lactobacillus" )
for (genus in levels(as.factor(combined_med_high_Mags$genus ))){
  
  reduced_tree <- drop.tip(snake.tree, snake.tree$tip.label[which(snake.tree$tip.label %!in% combined_med_high_Mags$id[which(combined_med_high_Mags$genus==genus)])])
  
  tree <- ggtree(  reduced_tree)#, layout = "circular"), open.angle = 0)
  
  p[[genus]] <- tree %<+% combined_med_high_Mags + geom_tippoint(aes(color=genus,shape=putative), size=5)+ #, shape=putative
    scale_shape_manual(values=c(16,18))+
    scale_size_manual(values=c(2,6))
  
  
  pal_complete <- colorRampPalette(c("blue" ))(100)
  pal_contamination <- colorRampPalette(c("green" ))(100)
  
  pal_strain_heterogeneity <- colorRampPalette(c("red" ))(100)
  
  
  
  
  heatmap_data_completeness <- dplyr::select(combined_med_high_Mags[which(combined_med_high_Mags$genus==genus),],"id", "Completeness.1") %>% remove_rownames() %>%  column_to_rownames("id")
  heatmap_data_contamination <- dplyr::select(combined_med_high_Mags[which(combined_med_high_Mags$genus==genus),], "id","Contamination.1")%>% remove_rownames() %>%  column_to_rownames("id")
  
  heatmap_data_pal_strain_heterogeneity <- dplyr::select(combined_med_high_Mags[which(combined_med_high_Mags$genus==genus),],"id","Strain_heterogeneity.1")%>% remove_rownames() %>%  column_to_rownames("id")
  
  
  heatmap_data_completeness$Completeness.1 <- as.numeric(heatmap_data_completeness$Completeness.1)
  
  heatmap_data_contamination$Contamination.1 <- as.numeric(heatmap_data_contamination$Contamination.1)
  heatmap_data_pal_strain_heterogeneity$Strain_heterogeneity.1 <- as.numeric(heatmap_data_pal_strain_heterogeneity$Strain_heterogeneity.1)
  
  
  # rownames(heatmap_data_completeness) <- combined_med_high_Mags$id[which(combined_med_high_Mags$genus==genus)]
  # rownames(heatmap_data_contamination) <- combined_med_high_Mags$id[which(combined_med_high_Mags$genus==genus)]
  # 
  # rownames(heatmap_data_pal_strain_heterogeneity) <- combined_med_high_Mags$id[which(combined_med_high_Mags$genus==genus)]
  # 
  
  
  
  
  
  
  
  
  
  
  
  #######################################################################################################
  library(ggnewscale)
  
  p1 <-gheatmap(p[[genus]],heatmap_data_completeness, offset =.03, width = .2,low = "lightblue", high = "blue", 
                colnames = FALSE,
                font.size=2,
                colnames_angle=90, 
                colnames_offset_y = .25,
                hjust = 0.5,
                legend_title = "Completeness")
  
  #merge(t, heatmap_data_completeness, by.x="id",by.y=0,all=TRUE)
  # <- data.frame(id=as.character(reduced_tree$tip.label))
  
  
  
  p2 <- p1 + new_scale_fill()
  
  p3 <- gheatmap(p2,heatmap_data_contamination, offset =.01, width = .2,low = "white", high = "green", 
                 colnames = FALSE,
                 font.size=2,
                 colnames_angle=90, 
                 colnames_offset_y = .25,
                 hjust = 0.5,
                 legend_title = "Contamination")
  
  
  
  
  p4 <- p3 + new_scale_fill()
  
  p[[genus]] <- gheatmap(p4,heatmap_data_pal_strain_heterogeneity, offset =.05, width = .2,low = "white", high = "red", 
                         colnames = FALSE,
                         font.size=2,
                         colnames_angle=90, 
                         colnames_offset_y = .25,
                         hjust = 0.5,
                         legend_title = "Strain heterogeneity")+
    ggtitle(genus)+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1.5, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(face="plain", size=20),
          legend.position="right",
          plot.title = ggtext::element_textbox_simple(face="bold.italic",halign  = 0.5,linetype = 1, # turn on border
                                                      box.color = "#748696",size=20, lineheight = 2))+# border color
    guides(color = guide_legend(title="Genus",
                                override.aes = list(size = 10),
                                label.theme = element_text(angle = 0, face = "italic",size = 20)),
           shape=guide_legend(title="Classification",
                              override.aes = list(size = 10),
                              label.theme = element_text(angle = 0,size = 20)))
  
  #    
}







# color=NULL, must be included in the gheatmap
#+scale_fill_viridis_c(option="A")


setwd("E:/STORE N GO/R/Plots/Phylogenetic trees")
png(filename='ggtree_metawrap_putative_gtdbtk_TESTER.png', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)

ggarrange(p[["Acetobacter"]], p[["Lactobacillus"]], nrow=1,ncol=2,labels=c("A.","B."),common.legend=TRUE)


graphics.off()



###########################################################################################################################################
#Incorporating the dram annotations

###########################################################################################################################################
