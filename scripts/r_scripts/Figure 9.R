





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
total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 1" )] ="Lactococcus lactis" 
total_prevalence$...5[which(total_prevalence$...5=="Lactococcus_lactis subcluster 2")] ="Lactococcus cremoris"

total_prevalence$...5[which(total_prevalence$...5=="Zymomonas_mobilis_subcluster 1")] ="Zymomonas mobilis"




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


#we applied the authors' recommended thresholds of 99.999% identity and 50% coverage.

colnames( global_mk_metadata)[which(colnames( global_mk_metadata)=="country")]="Sample_grain"


global_wk_metadata$Sample_grain = global_wk_metadata$Sample




global_metadata <- c()
global_metadata[["Milk.kefir"]] <- global_mk_metadata

global_metadata[["Water.kefir"]] <- global_wk_metadata






species="Lactococcus cremoris"



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


milk_taxonomic_profile_prevalence

# type="Milk.kefir"
# 
# t1 <- instrain[which(instrain$classification %in% 
#                       gsub("_"," ", total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]) &
#                        instrain$category.y==type&
#                        instrain$popANI_reference>.98),]
# 
# t1$cluster=
# gsub(".*_","",t1$cluster)
# 
# 
# levels(as.factor(t1$classification))
# 
# 
# #country
# data_long <- 
#   dplyr::select(t1,sample_id, classification, cluster) %>% 
#   mutate(present = 1) %>%
#   pivot_wider(names_from = cluster, values_from = present, values_fill = 0) %>%
#   pivot_longer(cols = c(`4`, `5`, `10`, `6`, `9`, `7`, `11`, `12`, `3`, `1`, `2`), 
#                names_to = "Strain_cluster", values_to = "Presence") %>%
#   mutate(Strain = as.numeric(Strain_cluster)) # Convert Strain to numeric for correct ordering
# 
# 
# ggplot(data_long, aes(x = Strain_cluster, y = sample_id, fill = Presence)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient(low = "white", high = "steelblue") +
#   facet_wrap(~ classification, scales = "free_y") + 
#   theme_minimal() +
#   labs(x = "Strain", y = "Sample ID", fill = "Presence", title = "Heatmap of Strain Presence by Sample") +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# 
# 
# 
#  # pivot_longer(cols = starts_with("`"), names_to = "Strain", values_to = "Presence")
#   pivot_longer(!classification, names_to = "cluster", values_to = "Freq")
#   #select(-classification) %>%
#   #column_to_rownames("sample_id")
# 
# 
# 
# 
# 
# 
# 




library(igraph)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(ggalluvial)
library(widyr )

multi_strain_detection=c()
unique_combinations_df=c()
occurence_data=c()


species="Lactococcus lactis"
type="Milk.kefir"

cluster_size_total=c()

cooccurrence_plot=c()
for (species in  levels(as.factor(instrain$classification))){
  
  for (type in c("Milk.kefir", "Water.kefir")){
    
    t1 <- instrain[which(instrain$classification==species &
      instrain$category.y==type&
        instrain$popANI_reference>.98),]
    
    if(nrow(t1)==0){
      
      next
    }
    t1=
      dplyr::select(t1,"sample_id" ,                 "genome"  ,                   "coverage"  ,                 "breadth", "Stage.y" ,   "Sample.y", "category.y","data_source.y","cluster", "classification"    )
    
    
    t1=
      
      merge(global_metadata[[type]] ,t1, by.y="sample_id",by.x="merge_column",all.y=TRUE )
    
    
    t1$Sample_grain[grep("ID",t1$Sample.y)] <-    t1$Sample.y[grep("ID",t1$Sample.y)] 
    
 
  
    
   t2=as.data.frame(table(t1$merge_column))

    
    multi_strain_detection=rbind(multi_strain_detection,
                                 data.frame(species=species,
                                            kefir=type,
                                            n_total_detections=    nrow(t2), 
                                            n_multi_detections= length(which(t2$Freq>1))
                                 
                                            
                                            ))
    
    
    if(length(which(t2$Freq>1))==0){
      
      print(paste( " We did not detect more than one strain cluster in more than one sample for the",species, "within the", type, "microbiome"))
      next     }
    
    
    
    
    
    
    
    
    cooccurrence_matrix <- 
      dplyr::select(t1, merge_column, classification, cluster) %>% 
      mutate(cluster=paste(classification,cluster)) %>% 
      pairwise_count(cluster, merge_column) %>% 
      spread(item2, n, fill = 0) %>%
      column_to_rownames(var = "item1")
    
    cooccurrence_long <- as.data.frame(cooccurrence_matrix) %>% rownames_to_column("Var2") %>% pivot_longer(!Var2,names_to = "Var1",values_to = "Freq")
    
    cooccurrence_plot[[species]] <- 
      # Assuming long format data
      ggplot(cooccurrence_long, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
      geom_alluvium(aes(fill = Var1)) +
      geom_stratum() +
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
      theme_void()
    
    
    pacman::p_load(widyr)
    cluster_size <-
    dplyr::select(t1,merge_column, classification, cluster) %>%
      dplyr::count(cluster)
    # 
    cooccurrence_matrix <-
    dplyr::select(t1,merge_column, classification, cluster) %>%
      mutate(cluster=paste(classification,cluster)) %>% 
      # group_by(merge_column) %>%
      # summarise(clusters = list(cluster)) %>%
      # unnest(clusters) %>%
      pairwise_count(cluster, merge_column) %>%  spread(item2, n, fill = 0) %>%
      column_to_rownames(var = "item1")

    # cooccurrence_matrix <- as.matrix(cooccurrence_matrix)
    # cooccurrence_matrix <- cooccurrence_matrix + t(cooccurrence_matrix)
    # graph <- graph_from_adjacency_matrix(cooccurrence_matrix, mode = "undirected", weighted = TRUE)
    # 
    # # Add cluster sizes as node attributes
    # V(graph)$size <- cluster_size$n[match(V(graph)$name, cluster_size$cluster)]
    # 
    # # Plot the network with ggraph
    # ggraph(graph, layout = "fr") +
    #   geom_edge_link(aes(edge_width = weight), edge_color = "blue", alpha = 0.6) +
    #   geom_node_point(aes(size = size), color = "red") +  # Size of the nodes based on the cluster size
    #   geom_node_text(aes(label = name), repel = TRUE) +  # Add labels for clusters
    #   theme_void() +
    #   ggtitle("Cluster Co-Occurrence Network") +
    #   scale_edge_width(range = c(0.5, 2)) +  # Adjust the edge width range
    #   theme(legend.position = "right")
    # 
    
    # #country
    # cooccurrence_matrix <- 
    #   dplyr::select(t1,merge_column, classification, cluster) %>% 
    #   mutate(present = 1) %>%
    #   pivot_wider(names_from = merge_column, values_from = present, values_fill = 0) %>%
    #   select(-classification) %>%
    #   column_to_rownames("cluster")
    # 
    
    
    
    
#     
#     
#     cooccurrence <- crossprod(as.matrix(cooccurrence_matrix))
#     graph <- graph_from_adjacency_matrix(cooccurrence, mode = "undirected", weighted = TRUE, diag = FALSE)
#     V(graph)$label <- rownames(cooccurrence_matrix)
#     
#     # Plot the network with vertex labels as clusters
#     plot(graph, 
#          vertex.label = V(graph)$label, # Use cluster names for labels
#          vertex.size = 10,              # Size of the nodes
#          edge.width = E(graph)$weight,  # Width of edges based on co-occurrence frequency
#          main = "Co-Occurrence Network Based on Clusters")
# 
#     # Plot the network
#     plot(graph, vertex.label = V(graph)$name, vertex.size = 10, edge.width = E(graph)$weight)
#     
#     
#     
#     melt(    cooccurrence_matrix %>% rownames_to_column("clust")) %>% 
#       ggplot(aes(variable, clust, fill = value)) +
#       geom_tile(color = "white") +
#       scale_fill_gradient(low = "white", high = "red") +
#       theme_minimal() +
#       labs(title = "Co-occurrence Heatmap of Clusters",
#            x = "Cluster",
#            y = "Cluster",
#            fill = "Co-occurrence") +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))
#     
#     
#     
#     
#     cooccurrence <- crossprod(as.matrix(cooccurrence_matrix))
#     
#     # Set the diagonal to 0 (no self-co-occurrence)
#     diag(cooccurrence) <- 0
#     
#     # Create a graph object from the co-occurrence matrix
#     graph <- graph_from_adjacency_matrix(cooccurrence, mode = "undirected", weighted = TRUE)
#     
#     # Convert igraph object to tidygraph
#     tidy_graph <- as_tbl_graph(graph)
#     
#     # Step 2: Plot the Co-occurrence Network
#     set.seed(42)  # For reproducibility of the layout
#     p1=  ggraph(tidy_graph, layout = 'fr') + 
#       geom_edge_link(aes(edge_width = weight), edge_colour = "lightblue", alpha = 0.7) +
#       geom_node_point(size = 6, colour = "steelblue") +
#       geom_node_text(aes(label = name), repel = TRUE, size = 4, colour = "black") +
#       scale_edge_width(range = c(0.2, 2)) +  # Adjust edge thickness
#       theme_void() +
#       labs(title = "Co-occurrence Network of Lactococcus cremoris Strains") +
#       theme(
#         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#         legend.position = "none"
#       )
#     
#     
#     jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/CS_metagenomics/cooccurnce_network_milk.kefir.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
#     
#     
# plot(p1)
#     graphics.off()
    
    cooccurrence_long$species=species
    cooccurrence_long$kefir=type
    
    
    cluster_size$species=species
    cluster_size$kefir=type
    
    cluster_size$var= paste(species,cluster_size$cluster)
    
    
    
occurence_data=rbind(               cooccurrence_long,
                                    occurence_data)
    
    cluster_size_total <- rbind( cluster_size,
                                 cluster_size_total)



    
    
  }
  
  
}



occurence_data_v2=
  occurence_data[which(occurence_data$Freq>=5),]


# occurence_data[
# which(occurence_data$Var1 %in%   cluster_size_total$var[which(cluster_size_total$n>5)] &
#         
#         occurence_data$Var2 %in%   cluster_size_total$var[which(cluster_size_total$n>5)]
#         
#         
#         
#         ),]
  


occurence_data_v2=
rbind(
occurence_data_v2[which(occurence_data_v2$species %in% gsub("_"," ",total_prevalence$...5[which(total_prevalence$kefir_type=="water")])&
                       
                       occurence_data_v2$kefir=="Water.kefir"
                       ),],




occurence_data_v2[which(occurence_data_v2$species %in% gsub("_"," ",total_prevalence$...5[which(total_prevalence$kefir_type=="milk")])&
                       
                       occurence_data_v2$kefir=="Milk.kefir"
),] ) 
  
  
 plot_data=  
   occurence_data_v2
 
 
  plot_data$Var1_c=
  
  str_split_fixed(plot_data$Var1,"_",2)[,2]


plot_data$Var2_c=
  
  str_split_fixed(plot_data$Var2,"_",2)[,2]

plot_data$Var1_primary=
  
  gsub(".* ","",plot_data$Var1)
 
plot_data$Var2_primary=
  
  gsub(".* ","",plot_data$Var2)

plot_data$kefir <- gsub("\\."," ",plot_data$kefir )

  
  
b=

ggplot(plot_data, aes(axis1 = Var1_primary, axis2 = Var2_primary, y = Freq)) +
  geom_alluvium(aes(fill = factor(Var1_c, levels = mixedsort(unique(Var1_c))))) +
  geom_stratum() +
  labs(title = "",
       x = "",
       y = "Number of strain Detections",
       fill = "Secondary cluster") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  facet_wrap(~species + kefir, scales = "free") +
  theme_bw() +
  theme(
    legend.position = c(0.78, 0.1),  # Adjust position to bottom left corner
    legend.justification = "left",  # Align legend to the left
    legend.direction = "horizontal", # Make the legend horizontal
    legend.box = "horizontal",  # Arrange in a horizontal box
    legend.box.just = "left",  # Align legend box to the left
    legend.spacing.x = unit(0.5, 'cm'), # Adjust spacing between items
    legend.spacing.y = unit(0.5, 'cm'), # Adjust spacing between rows
    legend.text = element_text(face = "italic", size = 20),
    legend.key.size = unit(1.5, 'cm'),  # Adjust legend key size
    legend.key.height = unit(1.5, 'cm'),  # Adjust legend key height
    legend.key.width = unit(1.5, 'cm'),  # Adjust legend key width
    legend.title = element_text(size = 20,hjust = 0.5),
    #legend.title.align = 0.5,  # Center-align title
    axis.text.y = element_text(size = 15, face = "italic"),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid labels
    panel.grid.minor = element_blank(),  # Remove minor grid labels
    plot.background = element_blank(),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(size = 15.5, face = "italic")
  )+
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5))


# jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/test2_v2.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
# 
# plot(strain_coccurence_plot)
# graphics.off()
# #



  #theme_void()




multi_strain_detection_prevelant=
  rbind(
    multi_strain_detection[which(multi_strain_detection$species %in% gsub("_"," ",total_prevalence$...5[which(total_prevalence$kefir_type=="water")])&
                              
                                   multi_strain_detection$kefir=="Water.kefir"
    ),],
    
    
    multi_strain_detection[which(multi_strain_detection$species %in% gsub("_"," ",total_prevalence$...5[which(total_prevalence$kefir_type=="milk")])&
                                   
                                   multi_strain_detection$kefir=="Milk.kefir"
    ),]
     ) 




a=
multi_strain_detection_prevelant %>% 
  mutate(kefir=gsub("\\."," ",kefir)) %>% 
  filter(n_multi_detections !=0) %>% 
 # filter(species %in% unique(plot_data$species)) %>% 
tidyr::gather(key = "type", value = "count", n_total_detections, n_multi_detections) %>% 
  mutate(type=
                gsub( "n_multi_detections", "Multi-Strain Detections",
                      gsub("n_total_detections" ,"Total Detections",type))) %>% 



  ggplot(aes(x = reorder(species, -count), y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Species",
       y = "Number of Detections",
       fill = "Detection Type") +
  coord_flip() +
  facet_wrap(~ kefir, scales = "free_y") +
  theme_bw()+
  theme( legend.position = c(0.35, .7), 
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.y = element_text(size=15,face="italic"),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
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
        strip.text = element_text(size=15.5))




# jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/test2_v2.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
# 
# plot(strain_coccurence_plot)
# graphics.off()
# #


library(ggpubr)

jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/v2/Figure 11.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)


ggarrange(a, b,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(size = 30),heights = c(.4,1), common.legend = FALSE)
graphics.off()





    
    total_prevalence
    
    length(names( table(t1$sample_id)))
    
    length(names( table(t1$sample_id)))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # #country
    # cooccurrence_matrix <- 
    #   dplyr::select(t1,sample_id, classification, cluster) %>% 
    #   mutate(present = 1) %>%
    #   pivot_wider(names_from = sample_id, values_from = present, values_fill = 0) %>%
    #   select(-classification) %>%
    #   column_to_rownames("cluster")
    # 
    # 
    
    country
    cooccurrence_matrix <-
      dplyr::select(t1,sample_id, classification, cluster) %>%
      mutate(present = 1) %>%
      pivot_wider(names_from = cluster, values_from = present, values_fill = 0) %>%
      select(-classification) %>%
      column_to_rownames("sample_id")

    
    
    
    
    
    cooccurrence <- crossprod(as.matrix(cooccurrence_matrix))
    
    # Set the diagonal to 0 (no self-co-occurrence)
    diag(cooccurrence) <- 0
    
    # Create a graph object from the co-occurrence matrix
    graph <- graph_from_adjacency_matrix(cooccurrence, mode = "undirected", weighted = TRUE)
    
    # Convert igraph object to tidygraph
    tidy_graph <- as_tbl_graph(graph)
    
    # Step 2: Plot the Co-occurrence Network
    set.seed(42)  # For reproducibility of the layout
    p1=  ggraph(tidy_graph, layout = 'fr') + 
      geom_edge_link(aes(edge_width = weight), edge_colour = "lightblue", alpha = 0.7) +
      geom_node_point(size = 6, colour = "steelblue") +
      geom_node_text(aes(label = name), repel = TRUE, size = 4, colour = "black") +
      scale_edge_width(range = c(0.2, 2)) +  # Adjust edge thickness
      theme_void() +
      labs(title = "Co-occurrence Network of Lactococcus cremoris Strains") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "none"
      )
    
    
    
    jpeg(filename='Q:/H2020 Master/Citizen Science Project/plots/CS_metagenomics/cooccurnce_network_milk.kefir.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
    
    
    plot(p1)
    graphics.off()
    
    
    library(igraph)
    
    # Calculate the co-occurrence matrix
    cooccurrence <- crossprod(as.matrix(cooccurrence_matrix))
    
    
    
    # Create a graph object
    graph <- graph.adjacency(cooccurrence, mode = "undirected", diag = FALSE)
    
    # Plot the network
    plot(graph, vertex.label = V(graph)$name, vertex.size = 10, edge.width = E(graph)$weight)
    
    
    
 
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



clus_breakdown_prevalent %>% 
  mutate(clust=gsub(".*_","",clust)) %>% 
  
  
  #ggplot(aes(x=  reorder(...2, count),y=count,fill=detection_category))+
  ggplot( aes(x = fct_reorder(species, Freq, .fun = sum),y=Freq,fill=clust))+
  geom_col()+
  facet_wrap(~type,scales="free")+
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












instrain_prevelant_cs_clust <- 
  instrain[which(
    instrain$popANI_reference>.98 &
      
      instrain$sample_id %in% kefir4all_metadata$merge_column),]




instrain_prevelant_cs_clust <- 
  rbind(
    instrain_prevelant_cs_clust[which(instrain_prevelant_cs_clust$category.y=="Milk.kefir" &
                                        instrain_prevelant_cs_clust$classification %in% total_prevalence$...5[which(total_prevalence$kefir_type=="milk")]),],
    
    instrain_prevelant_cs_clust[which(instrain_prevelant_cs_clust$category.y=="Water.kefir" &
                                        instrain_prevelant_cs_clust$classification %in% total_prevalence$...5[which(total_prevalence$kefir_type=="water")]),]
  )













library(ggalluvial)



instrain_prevelant_cs_clust <- 
  instrain_prevelant_cs_clust %>% 
  mutate(clusters=gsub(".*_","",cluster)) %>% 
  
  mutate(  Stage =  gsub("T1", "wk01",
                         gsub( "T2", "wk05",
                               gsub( "T3", "wk09",
                                     gsub("T4", "wk13",
                                          gsub("T5", "wk17",
                                               gsub( "T6", "wk21", Stage.y)))))),
           Stage = factor(Stage, levels = c("T0", "wk01", "wk05", "wk09", "wk13", "wk17","wk21")))



plot <- 
  as.data.frame(
    xtabs(~Stage+clusters+ classification+category.y,instrain_prevelant_cs_clust )) %>% 
  filter(Freq !=0) %>% 
  
  mutate(
    Stage = factor(Stage, levels = c("T0", "wk01", "wk05", "wk09", "wk13", "wk17","wk21")),
    clusters=factor(clusters,levels=c(1:20))
    
  )%>% 
  #mutate(c
  
  ggplot( aes(x = Stage, stratum = clusters, alluvium = clusters, y = Freq, fill = clusters)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum() +
  # scale_fill_brewer(type = "qual")+#, palette = "Set1") +
  theme_bw() +
  labs(title = "",
       x = "Stage",
       y = "Number of strains detected")+
  
  facet_wrap(~category.y+classification,scales = "free_y")+
  theme(legend.position = "right",#axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=15,,hjust = .5,vjust = .6),
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
        strip.text = element_text(size=15.5,face="italic"))+
  scale_colour_discrete(l = 40)




#########################################


#########################################







jpeg(filename='Q:/H2020 Master/Citizen Science Project/Manuscripts/CS_Metagenomics/Figures -CS_Metagenomics/instrain_strain_clusters_cs_v3.jpeg', width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)


plot(plot)
graphics.off()
