# Anchor to the project root before here() loads. RStudio Server resets the session cwd
# to $HOME after sourcing .Rprofile, and here() computes its root AT PACKAGE LOAD
# TIME, so force the project root before any here load. This block is inert on non-portable
for (.proj in c("/home/rstudio/kefir4all", "/workspaces/kefir4all_review_history")) {



  if (dir.exists(file.path(.proj, "scripts")) && getwd() != .proj) { setwd(.proj); break }
}


# ==============================================================================
# PATH CONFIGURATION — portable relative paths via the here package
# Run install.packages("here") if not installed.
# For files listed in data/private/README.md, provide your own copies.
# ==============================================================================
if (!requireNamespace("here",   quietly = TRUE)) install.packages("here")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(here)
source(here::here("scripts", "r_scripts", "figure_style.R"))
kefir_use_style()
library(tibble)  # needed for column_to_rownames() used below
DATA_DIR    <- here::here("data")                         # public data files
PRIVATE_DIR <- here::here("data", "private")              # private/sensitive files
OUT_DIR <- here::here("output", "04_taxonomic_profiling")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIGURES_DIR <- here::here("figures", "04_taxonomic_profiling")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# Private metadata — not distributed with repo (participant identifiers)
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

# .libPaths() removed — use default R library or renv
library(ComplexHeatmap)
library(ape)
library(devtools)
library(ggtree)
library(gplots)
library(plyr)
library(readr)
library(scales)
library(taxize)
library(tidyr)

library(vegan)
library(ggplot2)
library(grid)


########################################################################################################################
#Import metacache data
########################################################################################################################
metacache <- read_csv(file.path(DATA_DIR, "04_metacache_total_species_profile.csv"))
metacache_strain <-  read_csv(file.path(DATA_DIR, "04_metacache_total_strain_profile.csv"))

metacache$sample_id <- gsub("-","_",metacache$sample_id)
metacache <- metacache[-c(which(metacache$sample_id=="TG_EC_S72")),]

########################################################################################################################
#Make a metadata file ONLY NEED TO RUN THIS ONCE
########################################################################################################################
# 
# # RMEOVE LANE NAMES IN sample code AS IT WAS DESIGNED TO RENAME FILES ON THE CLUSTER PRIOR TO COMPUTATIONAL ANLAYSIS 
# 
# kefir4all_metadata$new_name <- gsub("_L00.","",kefir4all_metadata$new_name)
# 
# kefir4all_metadata[-c(which(duplicated(kefir4all_metadata$new_name))),]
# 
# 
# dplyr::select(kefir4all_metadata,new_name,`2`,`3`,`4`)%>%
# 
# 
# 
# kefir4all_metadata[-c(grep("ID",kefir4all_metadata$Sample)),]$Stage <- "T0"
# kefir4all_metadata[-c(grep("ID",kefir4all_metadata$Sample)),]$Stage
# 
# kefir4all_metadata$data_source <- "This study"
# 
# write.csv(kefir4all_metadata, file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"), quote = FALSE,row.names = FALSE)
#note manual fix samll errors like extraction control
########################################################################################################################
#Import metadata
########################################################################################################################

global_mk_metadata <- read_csv(file.path(DATA_DIR, "global_milk_kefir_metadata_v1.csv"))
global_wk_metadata <- read_csv(file.path(DATA_DIR, "global_water_kefir_metadata_v1.csv"))
Citizen_Scientist_metadata_v8 <- read_csv(CS_METADATA_PRIVATE)

Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] <- gsub("ID","ID00",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==3)] )
Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] <- gsub("ID","ID0",Citizen_Scientist_metadata_v8$ID[which(nchar(Citizen_Scientist_metadata_v8$ID)==4)] )



kefir4all_metadata <- read_csv(file.path(DATA_DIR, "kefir4all_sample_metadata_v2.csv"))
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
# 
#   # Filter data first
# 
# ########


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
        panel.border = element_blank())+
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





########################################################################################################################
# PCoA split by kefir type (milk vs water) -- addressing the reviewer concern that pooling both types in
# one ordination mostly just recovers the (expected) milk-vs-water separation. Each panel below is a FRESH
# Bray-Curtis + PCoA computed only on that type's samples, not just a subset of the combined ordination's
# points, so the axes reflect within-type structure rather than between-type structure.
########################################################################################################################

# Fixed palette so each dataset keeps the same colour in both panels (ggplot's default hue
# assignment is per-panel/alphabetical, so "Kefir4All" was coming out red in the milk panel but
# teal in the water panel). Four genuinely distinct, colourblind-safe hues (Okabe-Ito-derived) --
# NOT a warm/cool pairing, since reusing "red vs blue" as the two-way contrast in both the milk
# and the water panel reads as if the same two things are being compared both times.
dataset_palette <- c(
  "Kefir4All - Milk kefir"       = "#D55E00",  # vermillion
  "Walsh et al - Milk kefir"     = "#56B4E9",  # sky blue
  "Kefir4All - Water kefir"      = "#009E73",  # bluish green
  "Breselge et al - Water kefir" = "#CC79A7"   # reddish purple
)

# Captures each panel's data aspect ratio (PC1 range / PC2 range) so the combined figure can give
# each coord_equal() panel a column width that actually matches what it needs -- otherwise cowplot
# hands it a column sized from a guess (rel_widths), coord_equal shrinks the panel to preserve the
# 1:1 data scaling, and the mismatch shows up as blank gutters on the sides (the "white space
# between A and B").
panel_aspects <- list()

split_pcoa_panel <- function(sample_ids, label) {
  mat <- metacache %>% remove_rownames() %>% column_to_rownames("sample_id")
  mat <- mat[rownames(mat) %in% sample_ids, , drop = FALSE]
  cat("\n=== ", label, " PCoA: ", nrow(mat), " samples ===\n", sep = "")

  dist <- vegdist(mat, na.rm = TRUE)
  ord <- cmdscale(dist, k = 2, eig = TRUE, add = TRUE)
  pct <- 100 * ord$eig / sum(ord$eig)
  cat(label, "percent variance explained: PC1 =", round(pct[1], 1), "%  PC2 =", round(pct[2], 1), "%\n")

  df <- data.frame(Row.names = rownames(mat), PC1 = ord$points[, 1], PC2 = ord$points[, 2])
  df <- merge(df, pcoa[, c("Row.names", "data_source_specific")], by = "Row.names")
  panel_aspects[[label]] <<- diff(range(df$PC1)) / diff(range(df$PC2))
  cat(label, "PC1/PC2 data aspect ratio:", round(panel_aspects[[label]], 2), "\n")
  # Fixed level order (Kefir4All before the published comparison dataset) so the legend
  # ordering is the same in both panels, matching the fixed colour palette.
  df$data_source_specific <- factor(df$data_source_specific, levels = names(dataset_palette))
  df$panel_label <- label
  cat(label, "group sizes:\n")
  print(table(df$data_source_specific))

  centroid <- aggregate(cbind(PC1, PC2) ~ data_source_specific, data = df, FUN = mean)

  ggplot(df, aes(PC1, PC2, colour = data_source_specific)) +
    geom_point(size = 5, alpha = .25) +
    geom_point(data = centroid, size = 10, shape = 21, alpha = 0.5, color = "black",
               aes(fill = data_source_specific), show.legend = FALSE) +
    stat_ellipse(geom = "polygon", aes(fill = data_source_specific), alpha = 0.25,
                 type = "norm", show.legend = FALSE) +
    scale_colour_manual(values = dataset_palette) +
    scale_fill_manual(values = dataset_palette) +
    # facet_wrap on a single-value column, rather than labs(title=...), gives a title styled
    # EXACTLY like Figure 2's actual live panels (grain_heatmap_mk/wk and the C/D heatmaps, all
    # of which use facet strips, not plot.title) -- bold, size 22, centred, no background box.
    # The grain_bar_mk/wk code elsewhere in the repo with plot.title size 15.5 is dead code, not
    # what Figure 2 actually renders; copying its style was the mistake made further up this chat.
    facet_wrap(~ panel_label) +
    labs(x = paste0("PCoA1 - ", round(pct[1]), "%"), y = paste0("PCoA2 - ", round(pct[2]), "%")) +
    coord_equal() + theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 22, face = "bold", colour = "black"),
          strip.background = element_blank(),
          plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
          axis.ticks = element_blank(), axis.text = element_blank(),
          axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          panel.background = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.background = element_blank(),
          # margin(l=4) on legend.text puts visible air between each circular key and its label
          # (the default gap reads as near-zero once the key is enlarged to size 10).
          legend.text = element_text(size = 14, margin = margin(l = 4)),
          legend.title = element_text(size = 16),
          legend.key.height = unit(1.1, "cm"),
          legend.spacing.x = unit(0.3, "cm"),
          panel.border = element_blank()) +
    # ncol = 1 stacks the two legend entries vertically instead of side-by-side -- the long
    # published-study labels (e.g. "Breselge et al - Water kefir") were overflowing the panel's
    # column width and getting clipped when laid out in a single horizontal row.
    # alpha = 1 here matters: the scatter layer is drawn at alpha = .25, and without this override
    # the legend swatches inherit that and render as pale/washed-out versions of the palette colours
    # -- a mismatch against panel C's legend, which shows the same hex values at full saturation.
    guides(colour = guide_legend(override.aes = list(size = 10, alpha = 1), title = "Data source", ncol = 1))
}

milk_ids  <- pcoa$Row.names[pcoa$data_source_specific %in% c("Kefir4All - Milk kefir", "Walsh et al - Milk kefir")]
water_ids <- pcoa$Row.names[pcoa$data_source_specific %in% c("Kefir4All - Water kefir", "Breselge et al - Water kefir")]
cat("\n=== Split-PCoA sample counts ===\n")
cat("Milk kefir ids:", length(milk_ids), "  Water kefir ids:", length(water_ids), "\n")

p_milk  <- split_pcoa_panel(milk_ids,  "Milk kefir")
p_water <- split_pcoa_panel(water_ids, "Water kefir")

# Panel labels ("A.", "B.", size 20) match the house style used for the other combined figures
# in this script (e.g. Figure_5's ggarrange(..., labels = c("A.", "B."), font.label = list(size = 20))).
# rel_widths uses each panel's own captured data aspect ratio (panel_aspects), so the column width
# cowplot hands it already matches what coord_equal() needs -- no guessed ratio, no blank gutters.
p_split <- cowplot::plot_grid(p_milk, p_water, nrow = 1, labels = c("A.", "B."), label_size = 20,
                               rel_widths = c(panel_aspects[["Milk kefir"]], panel_aspects[["Water kefir"]]),
                               align = "hv", axis = "tb")
ggsave(file.path(FIGURES_DIR, "Figure_4_split_pcoa.tiff"), p_split, width = 16, height = 8, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Figure_4_split_pcoa_preview.png"), p_split, width = 16, height = 8, dpi = 150, bg = "white")
cat("\nWrote:", file.path(FIGURES_DIR, "Figure_4_split_pcoa.tiff"), "\n")


str(data.bray)
data.bray_metadata=data.frame(Row.names=rownames(as.matrix(data.bray)))

data.bray_metadata=
  merge(data.bray_metadata, pcoa, by="Row.names",all.x=TRUE, sort = F)

t <- 
  adonis2(data.bray~data.bray_metadata[,"data_source_specific"])


if (!requireNamespace("pairwiseAdonis", quietly=TRUE)) devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", upgrade="never")
library(pairwiseAdonis)
#default is 999 permutations
res<-pairwiseAdonis::pairwise.adonis(data.bray,data.bray_metadata$data_source_specific)


res<-pairwiseAdonis::pairwise.adonis(data.bray, data.bray_metadata[,"data_source_specific"])





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
      
      t <- tryCatch(
        anosim(subset_data1, group = subset_metadata$data_source_specific),
        error = function(e) { message("anosim skipped (", i, " vs ", j, "): ", e$message); NULL }
      )
      if (!is.null(t))
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

########################################################################################################################
# Figure 4 panel C: the community types actually observed (by dominant species) in milk and water
# kefir, and how many samples from each dataset -- including Kefir4All's own ("this study") samples --
# fall into each one. Reuses the dominant-species-per-sample table above (already thresholded at
# >=10% prevalence within milk/water) rather than any new per-sample clustering.
########################################################################################################################

dominanting_species_v2$kefir_group <- ifelse(dominanting_species_v2$group == "milk", "Milk kefir", "Water kefir")
dominanting_species_v2$dataset_label <- dplyr::case_when(
  dominanting_species_v2$group == "milk"  & dominanting_species_v2$data_source == "This study" ~ "Kefir4All - Milk kefir",
  dominanting_species_v2$group == "milk"  & dominanting_species_v2$data_source != "This study" ~ "Walsh et al - Milk kefir",
  dominanting_species_v2$group == "water" & dominanting_species_v2$data_source == "This study" ~ "Kefir4All - Water kefir",
  dominanting_species_v2$group == "water" & dominanting_species_v2$data_source != "This study" ~ "Breselge et al - Water kefir"
)

# Re-apply the >=10% prevalence threshold PER GROUP. species_of_interest above is a single list
# shared across milk and water, so a species that clears 10% in milk (e.g. Lactococcus lactis/
# cremoris) slips into the water panel too even though it's only a handful of water samples'
# dominant species there -- not a real water-kefir community type, just noise from the shared list.
group_totals <- dominanting_species %>%
  dplyr::distinct(sample_id, group) %>%
  dplyr::count(group, name = "n_group_total")

species_prevalence <- dominanting_species_v2 %>%
  dplyr::count(group, species, name = "n_species") %>%
  dplyr::left_join(group_totals, by = "group") %>%
  dplyr::mutate(prevalence = n_species / n_group_total)

cat("\n=== Species prevalence within group (dominant-species basis) -- kept if >= 10% ===\n")
print(species_prevalence %>% dplyr::arrange(group, dplyr::desc(prevalence)))

keep_species <- species_prevalence %>% dplyr::filter(prevalence >= 0.1)
dominanting_species_v3 <- dominanting_species_v2 %>%
  dplyr::semi_join(keep_species, by = c("group", "species"))

cat("\n=== Community types observed (dominant species), sample counts by dataset ===\n")
print(table(dominanting_species_v3$species, dominanting_species_v3$dataset_label))

community_type_counts <- dominanting_species_v3 %>%
  dplyr::count(kefir_group, dataset_label, species, name = "n_samples")

species_order <- community_type_counts %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(total = sum(n_samples), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::pull(species)
community_type_counts$species <- factor(community_type_counts$species, levels = species_order)
community_type_counts$dataset_label <- factor(community_type_counts$dataset_label, levels = names(dataset_palette))

# Styled to match panels A/B: white background, black box border, same type sizes.
# The legend is driven by a dummy, invisible geom_point layer (not the bar fill) so the keys
# render as solid circles like A/B's point legend -- a fill legend from geom_col always draws
# as rectangles no matter what override.aes is given, since the key glyph is tied to the geom.
p_community_types <-
  ggplot(community_type_counts, aes(x = species, y = n_samples)) +
  geom_col(aes(fill = dataset_label), colour = "black", linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(x = species, y = 0, colour = dataset_label), size = 0, alpha = 0) +
  facet_wrap(~ kefir_group, scales = "free", nrow = 1) +
  scale_fill_manual(values = dataset_palette) +
  scale_colour_manual(values = dataset_palette) +
  labs(x = NULL, y = "Number of samples", colour = "Data source") +
  # Wrapped, horizontal (not rotated) labels instead of fighting rotated-text overflow: at 40
  # degrees, "Lactobacillus kefiranofaciens" projects roughly 2 inches horizontally past its own
  # tick, which clip = "off" / wider margins / wider panel.spacing all failed to fully contain
  # (it kept getting clipped somewhere -- by the panel, by the facet, or by cowplot's own
  # viewport). Wrapping onto 2 short horizontal lines avoids the problem rather than chasing it.
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 14)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, face = "italic", colour = "black"),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 20),
        axis.ticks.x = element_blank(),
        # Facet strip styled to match the EXACT house convention used for the equivalent faceted
        # panels (C/D heatmaps) in Figure 2: bold, size 22, no background box
        # (heatmap_theme in 04_taxonomic_profiling.R: strip.text size=22 face="bold", strip.background
        # = element_blank()) -- same semantic "Milk kefir"/"Water kefir" label as there.
        strip.text = element_text(size = 22, face = "bold", colour = "black"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Extra top margin so the "C." cowplot panel label (drawn above the plot, same corner as
        # the first facet's strip text) doesn't collide with "Milk kefir". No longer need the
        # wide left/right margin or panel.spacing that rotated-label overflow required -- wrapped
        # horizontal labels stay inside their own panel.
        plot.margin = margin(t = 30, r = 5, b = 5, l = 5),
        plot.background = element_blank(),
        legend.position = "bottom",
        # margin(l=4) matches A/B: visible air between each circular key and its label.
        legend.text = element_text(size = 14, margin = margin(l = 4)),
        legend.title = element_text(size = 16),
        legend.spacing.x = unit(0.3, "cm"),
        legend.key.size = unit(1.1, "cm"),
        panel.border = element_blank(),
        panel.spacing = unit(1.5, "lines")) +
  # size = 10, alpha = 1 in the override brings the (invisible, size-0/alpha-0) dummy points back
  # up to the same visual weight as A/B's circular legend keys. nrow = 1: single row across all
  # four entries, not stacked 2x2.
  guides(colour = guide_legend(nrow = 1,
                                override.aes = list(size = 10, alpha = 1)))

# Plain "A."/"B."/"C." labels only -- matching the house style used elsewhere in this script
# (e.g. Figure_5's labels = c("A.", "B.")). Each panel keeps its own in-plot "Milk kefir"/"Water
# kefir" title (bold facet strips), matching Figure 2's convention.
# All three panels share the SAME "Data source" colour scheme, so one combined legend (pulled from
# panel C, the only one with all 4 categories present) replaces the three separate, partially
# redundant per-panel legends -- cuts the repeated vertical whitespace those were adding up.
legend_combined <- cowplot::get_legend(p_community_types + theme(legend.position = "bottom"))

# A single 22-inch-wide row made all three panels shallow and left large empty areas around the
# coord_equal() PCoAs. Put A/B together on the top row and give C the full lower row instead.
# A/B retain aspect-derived relative widths, while the deeper canvas lets their plotting regions
# grow substantially without sacrificing Panel C's long species labels.
# TOP_ROW_HEIGHT must be big enough to hold BOTH the non-panel chrome (top margin + the size-22
# facet-strip title + x-axis title + bottom margin, measured via ggplotGrob() introspection at
# ~1.05in) AND the full coord_equal() panel height implied by FIG_WIDTH / sum(panel_aspects)
# (~7.49in here) -- otherwise the panel is height-constrained below its width allocation, coord_equal
# shrinks it to fit, and the leftover space shows up as a visible gap between the strip title and the
# panel (and as gutters beside it) instead of the title sitting snug against the panel like Figure 3's.
FIG_WIDTH         <- 15.5
TOP_ROW_HEIGHT    <- 8.6
TYPES_ROW_HEIGHT  <- 4.3
LEGEND_HEIGHT     <- 0.65

top_grid <- cowplot::plot_grid(
  p_milk + theme(legend.position = "none"),
  p_water + theme(legend.position = "none"),
  nrow = 1,
  labels = c("A.", "B."),
  label_size = 20,
  rel_widths = c(panel_aspects[["Milk kefir"]], panel_aspects[["Water kefir"]]),
  align = "hv",
  axis = "tb"
)

types_grid <- cowplot::plot_grid(
  p_community_types + theme(legend.position = "none"),
  ncol = 1,
  labels = "C.",
  label_size = 20
)

panels_grid <- cowplot::plot_grid(
  top_grid,
  types_grid,
  ncol = 1,
  rel_heights = c(TOP_ROW_HEIGHT, TYPES_ROW_HEIGHT)
)

p_split_types <- cowplot::plot_grid(
  panels_grid,
  legend_combined,
  ncol = 1,
  rel_heights = c(TOP_ROW_HEIGHT + TYPES_ROW_HEIGHT, LEGEND_HEIGHT)
)
ggsave(file.path(FIGURES_DIR, "Figure_4_split_pcoa_community_types.tiff"), p_split_types,
       width = FIG_WIDTH, height = TOP_ROW_HEIGHT + TYPES_ROW_HEIGHT + LEGEND_HEIGHT,
       dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Figure_4_split_pcoa_community_types_preview.png"), p_split_types,
       width = FIG_WIDTH, height = TOP_ROW_HEIGHT + TYPES_ROW_HEIGHT + LEGEND_HEIGHT,
       dpi = 150, bg = "white")
cat("\nWrote:", file.path(FIGURES_DIR, "Figure_4_split_pcoa_community_types.tiff"), "\n")

########################################################################################################################
# Kefir4All community-type drift over time (ggalluvial): does a citizen-science kefir line's
# starting dominant-species community type stay put, or drift into a different type, across its
# own repeat timepoints? Kefir4All ONLY -- Walsh et al / Breselge et al are cross-sectional
# (single timepoint per sample), so there's no "drift" to show for them.
# Stage -> week recode and participant ID ("Sample" column) follow the same convention already
# used in 06_instrain_temporal_alluvial.R and 04_community_stability_v2.R.
########################################################################################################################

library(ggalluvial)

stage_levels <- c("T0", "wk01", "wk05", "wk09", "wk13", "wk17", "wk21")
recode_stage <- function(x) {
  dplyr::recode(x, T1 = "wk01", T2 = "wk05", T3 = "wk09", T4 = "wk13", T5 = "wk17", T6 = "wk21", T0 = "T0")
}

alluvial_data <- dominanting_species_v3 %>%
  dplyr::filter(data_source == "This study") %>%
  dplyr::inner_join(
    kefir4all_metadata %>% dplyr::select(merge_column, Sample, Stage),
    by = c("sample_id" = "merge_column")
  ) %>%
  dplyr::filter(!is.na(Stage), Stage %in% c("T0", "T1", "T2", "T3", "T4", "T5", "T6")) %>%
  dplyr::mutate(Stage = factor(recode_stage(Stage), levels = stage_levels),
                kefir_group = ifelse(group == "milk", "Milk kefir", "Water kefir"),
                # Grain (MG/WG) and liquid (ML/WL) are different biological compartments for the
                # SAME participant at the SAME stage -- collapsing both into one alluvium per
                # participant creates duplicate (alluvium, x) pairs, which geom_flow() rejects
                # ("Data is not in a recognized alluvial form"). Track them as separate lines.
                compartment = ifelse(`kefir type` %in% c("MG", "WG"), "Grain", "Liquid"),
                participant_line = paste(Sample, compartment, sep = "_"))

cat("\n=== Kefir4All alluvial: sample counts per Stage x kefir type ===\n")
print(table(alluvial_data$Stage, alluvial_data$kefir_group))

participant_counts <- alluvial_data %>% dplyr::distinct(participant_line, kefir_group, compartment) %>%
  dplyr::count(kefir_group, compartment)
cat("\nDistinct participant lines per kefir type x compartment:\n")
print(participant_counts)

# How many participant lines actually have >1 timepoint (otherwise there's no "drift" to draw)?
line_stage_counts <- alluvial_data %>% dplyr::distinct(participant_line, Stage, kefir_group) %>%
  dplyr::count(participant_line, kefir_group, name = "n_timepoints")
cat("\nParticipant lines with >1 timepoint (the ones an alluvial ribbon actually tracks):\n")
print(line_stage_counts %>% dplyr::filter(n_timepoints > 1) %>% dplyr::count(kefir_group))

# Safety net: if any (participant_line, Stage) pair still has >1 row (e.g. a genuine duplicate
# sample), keep just one so geom_flow() doesn't choke -- report it rather than silently assuming.
dupes <- alluvial_data %>% dplyr::count(participant_line, Stage) %>% dplyr::filter(n > 1)
if (nrow(dupes) > 0) {
  cat("\nWARNING:", nrow(dupes), "(participant_line, Stage) pairs still have >1 row -- deduplicating (keeping first):\n")
  print(dupes)
  alluvial_data <- alluvial_data %>% dplyr::group_by(participant_line, Stage) %>% dplyr::slice(1) %>% dplyr::ungroup()
}

species_palette <- setNames(
  scales::hue_pal()(length(unique(alluvial_data$species))),
  sort(unique(alluvial_data$species))
)

p_alluvial <-
  ggplot(alluvial_data, aes(x = Stage, stratum = species, alluvium = participant_line, fill = species)) +
  ggalluvial::geom_flow(stat = "alluvium", alpha = 0.7, colour = "darkgray") +
  ggalluvial::geom_stratum(width = 0.5) +
  # Grain and liquid are separate panels, not just separate lines within one panel -- they're
  # different compartments (the kefir grains themselves vs. the fermented milk/water around them),
  # so mixing their flows in the same axis would read as one undifferentiated community.
  facet_grid(compartment ~ kefir_group, scales = "free_x") +
  scale_fill_manual(values = species_palette, name = "Community type\n(dominant species)") +
  labs(x = NULL, y = "Number of samples") +
  theme_bw(base_size = 11) +
  # Same facet-strip convention as the rest of the figure (and Figure 2): bold, size 22, centred,
  # no background box. Previously knocked down to 16 to squeeze into a 2x2 grid -- restored properly
  # and given the figure more height below instead of shrinking the text.
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        strip.text = element_text(size = 22, face = "bold", colour = "black"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(face = "italic"))

ggsave(file.path(FIGURES_DIR, "Figure_4_community_type_alluvial.tiff"), p_alluvial,
       width = 14, height = 11, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Figure_4_community_type_alluvial_preview.png"), p_alluvial,
       width = 14, height = 11, dpi = 150, bg = "white")
cat("\nWrote:", file.path(FIGURES_DIR, "Figure_4_community_type_alluvial.tiff"), "\n")

########################################################################################################################
# Full Figure 4: ordination (A/B) + Kefir4All community-type drift (C). Panel C (the community
# type counts bar chart, still saved separately above as Figure_4_split_pcoa_community_types.tiff)
# is dropped from this combined figure -- the alluvial plot already shows the same community types
# per dataset, plus their dynamics over time, so the static count chart was redundant here.
# The alluvial keeps its own legend (species, not "Data source") since it's a different
# categorical variable -- two legends is the right call, not one merged legend.
########################################################################################################################

AB_grid <- cowplot::plot_grid(p_milk + theme(legend.position = "none"),
                               p_water + theme(legend.position = "none"),
                               nrow = 1, labels = c("A.", "B."), label_size = 20,
                               rel_widths = c(panel_aspects[["Milk kefir"]], panel_aspects[["Water kefir"]]),
                               align = "hv", axis = "tb")

AB_with_legend <- cowplot::plot_grid(AB_grid, legend_combined, ncol = 1,
                                      rel_heights = c(TOP_ROW_HEIGHT, LEGEND_HEIGHT))

c_grid <- cowplot::plot_grid(p_alluvial, labels = "C.", label_size = 20)

AB_HEIGHT <- TOP_ROW_HEIGHT + LEGEND_HEIGHT
C_HEIGHT  <- 11

p_figure4_full <- cowplot::plot_grid(AB_with_legend, c_grid, ncol = 1,
                                      rel_heights = c(AB_HEIGHT, C_HEIGHT))

ggsave(file.path(FIGURES_DIR, "Figure_4_full.tiff"), p_figure4_full,
       width = FIG_WIDTH, height = AB_HEIGHT + C_HEIGHT, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(FIGURES_DIR, "Figure_4_full_preview.png"), p_figure4_full,
       width = FIG_WIDTH, height = AB_HEIGHT + C_HEIGHT, dpi = 150, bg = "white")
cat("\nWrote:", file.path(FIGURES_DIR, "Figure_4_full.tiff"), "\n")


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

  ########
labs(x=paste("PCoA1 - ", round(percent_explained[1]), "%", sep=""), y=paste("PCoA2 - ", round(percent_explained[2]), "%", sep=""), title="") +
  coord_equal() +
  theme_bw()+
  theme(legend.position = "top",#c(.85,.2),#axis.text.x = element_blank(),  # remove x-axis text
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
        panel.border = element_blank())+
  #legend.box.background = element_rect(colour = "red"))+
  scale_colour_discrete(name="Data source")+
  scale_shape_discrete(name="Data source")+
  guides(colour = guide_legend(override.aes = list(size=10),ncol=2,byrow=TRUE))



# facet by species and RA
  
#   pcoa_v2 %>% 
#   
#   
#   
#   ########
# #labs(x=paste("PCoA1 - ", round(percent_explained[1]), "%", sep=""), y=paste("PCoA2 - ", round(percent_explained[2]), "%", sep=""), title="") +
# #coord_equal() +
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
        panel.border = element_blank(),
        legend.text = element_text(face = "plain", size = 15),
        legend.key.size = unit(1.5, 'cm'),
        legend.key.height = unit(1.5, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        legend.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 20)) +
  guides(shape = guide_legend(override.aes = list(size = 10)))




# pcoa_v2 %>% 





library(ggpubr)
tiff(filename=file.path(FIGURES_DIR, 'Figure_5.tiff'), width = 1800, height = 1200, units = "px", res = 300, compression = "lzw", type = "cairo", pointsize = 12) #, width=2000, height=1950)


ggarrange(p_total_v2,
          p_species,
          nrow=1,ncol=2,labels=c("A.","B."),  font.label = list(size = 20))
 graphics.off()

 
 
 
 library(Cairo)   # For high-quality PNG/PDF/SVG
 
 # Define figure size and resolution
 fig_width <- 6  # Inches (adjust based on journal requirements)
 fig_height <- 4  # Inches
 dpi_res <- 300   # High resolution (300 dpi for journal quality)
 
 
 
 combined_plot <- ggarrange(
   p_total_v2, p_species, 
   nrow = 1, ncol = 2, 
   labels = c("A.", "B."), 
   font.label = list(size = 20)
 )
 
 
  pdf(file.path(FIGURES_DIR, "Figure5_model.pdf"), width = fig_width, height = fig_height)#, width = fig_width, height = fig_height)
  print(combined_plot)
  dev.off()

 svg(file.path(FIGURES_DIR, "Figure5_model.svg"), width = fig_width, height = fig_height)
 print(combined_plot)
 dev.off()
 
 
 tiff(file.path(FIGURES_DIR, "Figure5_model.tiff"), width = fig_width, height = fig_height, units = "in", res = dpi_res, compression = "lzw", type = "cairo")
 print(combined_plot)
 dev.off()
 
 
 

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
    theme(legend.position = "top",#axis.text.x = element_blank(),  # remove x-axis text
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

tiff(filename=file.path(FIGURES_DIR, 'Figure_5_v2.tiff'), width = 1800, height = 1200, units = "px", res = 300, compression = "lzw", type = "cairo", pointsize = 12) #, width=2000, height=1950)


ggarrange(p_total,
  ggarrange(  pcoa_plot[[1]],  pcoa_plot[[2]], nrow=1,ncol=2, labels=c( "B.","C."),common.legend = FALSE,  font.label = list(size = 30)),
  nrow=2,ncol=1,labels=c("A.","C."),widths=c(10,1),heights = c(1,.5),  font.label = list(size = 30))
graphics.off()
