# 
# # Update all CRAN packages
# update.packages(ask = FALSE, checkBuilt = TRUE)
# 
# # Update all Bioconductor packages
# BiocManager::install(ask = FALSE)
# 
# # Restart R session
# .rs.restartR()
# Clear global environment ####
rm(list=ls())

# 1. Load required libraries ####
install.packages("BiocManager",version=3.18)
# install.packages("ggplot2")
library(ggplot2)
library(ggpubr)
# BiocManager::install("ggtree")
library(ggtree)
packageVersion("ggtree")
# ‘3.16.3’
library(RColorBrewer)
# install.packages("ggfun", dependencies = TRUE)
library(ggfun)
library(stringr)
library(tidyverse)
library(pheatmap)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtreeExtra")
library(ggtreeExtra)

# Set directory 
dir.home <- here::here() 

#Load in data
Object <- readRDS("/Users/hhampson/Library/CloudStorage/GoogleDrive-hhampson@usc.edu/Shared drives/ATHLETE_CHS_Environment_Microbiome/0_data/CHS/Formatted Data/FormattedObject.RDS")
table <- data.frame(Object$Table)
#Remove mitochondria
table <- table[,-149]

taxa_tab2 <- as.data.frame(t(table))
taxa_tab2 <- taxa_tab2[c(1:206),]

# Change rownames (taxa names) to column
taxa_tab2 <- tibble::rownames_to_column(taxa_tab2, "V1")

# 2. Clean Taxa Table Names (JG Comment: There was an issue with this code) ----
species_names_og <- taxa_tab2[,1]
phylogenic_tree_names_df <- data.frame(species_names_og) |> 
  tidylog::mutate(
    species_names_full = str_replace(
      species_names_og, "_p_", "--p_") |>
      str_replace("_c_", "--c_") |>
      str_replace("_o_", "--o_") |>
      str_replace("_f_", "--f_") |>
      str_replace("_g_", "--g_") |>
      str_replace("_s_", "--s_"), 
    genus_name_full  = sapply(strsplit(species_names_full,"--s"),"[[",1),
    order_name_full  = sapply(strsplit(species_names_full,"--f"),"[[",1),
    family_name_full = sapply(strsplit(species_names_full,"--g"),"[[",1),
    class_name_full  = sapply(strsplit(species_names_full,"--o"),"[[",1),
    phylum_name_full = sapply(strsplit(species_names_full,"--c"),"[[",1),
    species_names = sapply(strsplit(species_names_full,"--"),"[[",7), 
    genus_name    = sapply(strsplit(species_names_full,"--"),"[[",6),
    family_name   = sapply(strsplit(species_names_full,"--"),"[[",5),
    order_name    = sapply(strsplit(species_names_full,"--"),"[[",4),
    class_name    = sapply(strsplit(species_names_full,"--"),"[[",3),
    phylum_name   = sapply(strsplit(species_names_full,"--"),"[[",2), 
    domain_name   = sapply(strsplit(species_names_full,"--"),"[[",1))

# get data frame of individual names only (not used later, only for reference)
individual_names_only <- phylogenic_tree_names_df %>%
  dplyr::select(species_names,genus_name,
                order_name,family_name,class_name,phylum_name)

# in new data, -- separates phylum, __ separates the level identifier from the name 
species_names <- phylogenic_tree_names_df$species_names_full
genus_names   <- unique(phylogenic_tree_names_df$genus_name_full)
family_names  <- unique(phylogenic_tree_names_df$family_name_full)
order_names   <- unique(phylogenic_tree_names_df$order_name_full)
class_names   <- unique(phylogenic_tree_names_df$class_name_full)
phylum_names  <- unique(phylogenic_tree_names_df$phylum_name_full)


## a. Defining taxa level names (JG- the issue above was causing issues here) ----
names_all <- c(species_names,phylum_names,class_names,
               order_names,family_names,genus_names) 
phylo_level <- sapply(strsplit(names_all,"--"),length) 

# Sanity check- are there the correct number of species:phylum
table(phylo_level)

## b. Defining branch tree text ----
# Create an empty list 
l1=list()  


# 3. Read in data to generate node size, direction and location ------
#Location of New Results - Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/2_Results
#BaH-ZING Means Results: Formatted_BaHZING_Mean_Results_03_08_2024.csv
#ZING Means Results: Formatted_ZING_Means_Results_03_09_2024.csv
# m1 <- read.csv("/Users/hehampson/USC Lab Dropbox/Chatzi Active Projects/MOVED Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/2_Results/Fiber_Sensitivity_BaH_ZING_Results_03_23_2024.csv") %>% 
m1 <- read.csv("/Users/hhampson/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/Chatzi Projects and Analyses/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/2_Results/Fiber_Sensitivity_BaH_ZING_Results_03_23_2024.csv") %>% 
# m1 <- read.csv(file = fs::path(dir.home %>% dirname(),
#                                "2_Results",
#                                "Fiber_Sensitivity_BaH_ZING_Results_03_23_2024.csv")) %>% 
  tidylog::filter(Exposure!="Dispersion") %>% 
  tidylog::filter(Exposure!="Omega") %>% 
  tidylog::mutate(
    Direction = dplyr::case_when(
      Component=="Probability" & Mean>1 ~ "Positive",
      Component=="Probability" ~ "Negative",
      Component=="Means" & Mean > 0 ~ "Positive",
      Component=="Means" ~ "Negative"),
    Mean_sig_only=ifelse(grepl("\\*",sig),Mean,NA),
    Direction_sig_only=ifelse(grepl("\\*",sig),Direction,NA),
    se = Mean/SD) %>% 
  mutate(Taxa = str_replace(Taxa, "_p_", "--p_")) %>% 
  mutate(Taxa = str_replace(Taxa,"_c_", "--c_")) %>%
  mutate(Taxa = str_replace(Taxa,"_o_", "--o_")) %>% 
  mutate(Taxa = str_replace(Taxa,"_f_", "--f_")) %>% 
  mutate(Taxa = str_replace(Taxa,"_g_", "--g_")) %>% 
  mutate(Taxa = str_replace(Taxa,"_s_", "--s_"))

# Filter to mixture effect only
psi <- m1 %>% dplyr::filter(Exposure=="Mixture")
betas <- m1 %>% dplyr::filter(Exposure!="Mixture")
# Reset row names after filtering
rownames(psi) <- NULL
rownames(betas) <- NULL

## Summarize betas ----
betas_summary <- betas |> 
  dplyr::group_by(Taxa, Component) |> 
  dplyr::summarise(n_sig = sum(sig=="*")) |> 
  dplyr::ungroup()

# Defining directions from data - individual effects
beta_dd <- betas_summary %>% dplyr::filter(Component=="Means")
beta_dd <- beta_dd %>% tidylog::filter(n_sig > 0)

# Defining directions from data - PSI
psi_dd <- psi %>% dplyr::filter(Component=="Means")
# psi_dd$Mean=abs(psi_dd$Mean)
psi_dd$group <- (psi_dd$Direction)
psi_dd <- psi_dd %>% 
  tidylog::filter(Taxa %in% beta_dd$Taxa | sig == "*") |>
  tidylog::select(Taxa,Mean,group, se, Direction_sig_only) 
# 
sig.only <- psi_dd %>%
  filter(grepl("s__",Taxa) | grepl("g__Dorea--s",Taxa)) %>%
  filter(!is.na(Direction_sig_only)| grepl("g__Dorea--s",Taxa))
  # filter(!is.na(Direction_sig_only))
sig.only <- sig.only$Taxa

#Filter taxa_tab2 to only significant species
colnames(table) <- str_replace(colnames(table), "_p_", "--p_")
colnames(table) <- str_replace(colnames(table),"_c_", "--c_") 
colnames(table) <- str_replace(colnames(table),"_o_", "--o_")
colnames(table) <- str_replace(colnames(table),"_f_", "--f_") 
colnames(table) <- str_replace(colnames(table),"_g_", "--g_")
colnames(table) <- str_replace(colnames(table),"_s_", "--s_")

# table <- table[sig.only]

taxa_tab2 <- as.data.frame(t(table))
taxa_tab2 <- taxa_tab2[c(1:206),]

# Change rownames (taxa names) to column
taxa_tab2 <- tibble::rownames_to_column(taxa_tab2, "V1")

# 2. Clean Taxa Table Names (JG Comment: There was an issue with this code) ----
species_names_og <- taxa_tab2[,1]
phylogenic_tree_names_df <- data.frame(species_names_og) |> 
  tidylog::mutate(
    species_names_full = str_replace(
      species_names_og, "_p_", "--p_") |>
      str_replace("_c_", "--c_") |>
      str_replace("_o_", "--o_") |>
      str_replace("_f_", "--f_") |>
      str_replace("_g_", "--g_") |>
      str_replace("_s_", "--s_"), 
    genus_name_full  = sapply(strsplit(species_names_full,"--s"),"[[",1),
    order_name_full  = sapply(strsplit(species_names_full,"--f"),"[[",1),
    family_name_full = sapply(strsplit(species_names_full,"--g"),"[[",1),
    class_name_full  = sapply(strsplit(species_names_full,"--o"),"[[",1),
    phylum_name_full = sapply(strsplit(species_names_full,"--c"),"[[",1),
    species_names = sapply(strsplit(species_names_full,"--"),"[[",7), 
    genus_name    = sapply(strsplit(species_names_full,"--"),"[[",6),
    family_name   = sapply(strsplit(species_names_full,"--"),"[[",5),
    order_name    = sapply(strsplit(species_names_full,"--"),"[[",4),
    class_name    = sapply(strsplit(species_names_full,"--"),"[[",3),
    phylum_name   = sapply(strsplit(species_names_full,"--"),"[[",2), 
    domain_name   = sapply(strsplit(species_names_full,"--"),"[[",1))

# get data frame of individual names only (not used later, only for reference)
individual_names_only <- phylogenic_tree_names_df %>%
  dplyr::select(species_names,genus_name,
                order_name,family_name,class_name,phylum_name)

# in new data, -- separates phylum, __ separates the level identifier from the name 
species_names <- phylogenic_tree_names_df$species_names_full
genus_names   <- unique(phylogenic_tree_names_df$genus_name_full)
family_names  <- unique(phylogenic_tree_names_df$family_name_full)
order_names   <- unique(phylogenic_tree_names_df$order_name_full)
class_names   <- unique(phylogenic_tree_names_df$class_name_full)
phylum_names  <- unique(phylogenic_tree_names_df$phylum_name_full)


## a. Defining taxa level names (JG- the issue above was causing issues here) ----
names_all <- c(species_names,phylum_names,class_names,
               order_names,family_names,genus_names) 
phylo_level <- sapply(strsplit(names_all,"--"),length) 

# Sanity check- are there the correct number of species:phylum
table(phylo_level)

## b. Defining branch tree text ----
# Create an empty list 
l1=list()  

# For testing:
# n=7
# m=names_all[which(phylo_level==n-1)][39]
# table(phylo_level)

## c. For loop to get Newick format tree structure -----
for (n in 7:3){
  l2=list()
  # For genus level vars:
  if(n==7){
    for (m in names_all[which(phylo_level==n-1)]){
      names2=names_all[phylo_level==n]
      l1[m]=paste0("(",paste(names2[grep(paste0(m, "-"),names2)],collapse=","),")",m)
    }
  } else{
    for (m in names_all[which(phylo_level==n-1)]){
      names2=names_all[phylo_level==n]
      l2[m]=paste0("(",paste(l1[names2[grep(m,names2)]],collapse=","),")",m)
    }
    l1=l2
  }
  print(paste0(n, ": ", length(l1), "; ", length(l2)))
}

## d. Defining parameters for branch tree text ----
# This code requires Newick format 
l3 <- paste0("(",paste(l1,collapse=","),")Bacteria;")

## e. Plot Base  ----
tree <- ape::read.tree(text = l3)
# library(igraph)
# igraph_tree <- ape::as.igraph.phylo(tree, directed = TRUE, ignore.root = TRUE, use.labels = TRUE)
# plot(igraph_tree, layout = layout_with_kk, 
#      vertex.size = 5,
#      vertex.label = NA,
#      # edge.curved = 0.1,
#      # vertex.label.cex = 0.5, 
#      edge.arrow.size = 0.1) 

(ggtree::ggtree(tree, layout = "fan", open.angle= 20))
# ggtree(tree)
# ggtree(tree, layout="roundrect")
# ggtree(tree, layout="slanted")
# ggtree(tree, layout="ellipse")
# ggtree(tree, layout="circular")
# ggtree(tree, layout="fan", open.angle=0)
# ggtree(tree, layout="equal_angle")
# ggtree(tree, layout="daylight")
# ggtree(tree, branch.length='none')
# ggtree(tree, layout="ellipse", branch.length="none")
# ggtree(tree, branch.length='none', layout='circular')
# ggtree(tree, layout="daylight", branch.length = 'none')





# Sanity check:
paste0("Number of nodes is ", length(tree$node.label),". It should be ", sum(table(phylo_level))-table(phylo_level)[6]+1, 
       " (the number of unique genus_names, family_names, order_names, class_names, phylum_names, and domain_names.")
paste0("Number of tips is ", length(tree$tip.label) ,". It should be ", length(species_names))


## Trouble shooting more than expected nodes
# str_count(l3, "s__unclassified113,")
# x <- purrr::map(l1, ~str_count(., "s__unclassified113,") ) |> unlist()
# str_count(l3, "s__unclassified1,")
# temp_df <- data.frame(tip.label = tree$tip.label) 
# temp_df_dup <- janitor::get_dupes(temp_df, tip.label)

# Creating empty vectors for nodes and tip labels
node_num <- vector()
tip_num <- vector()

# Assigning labels to the correct locations
for (i in 1:length(phylum_names)) {
  print(i)
  tips_matching_taxon <- tree$tip.label[grep(as.character(phylum_names[i]), tree$tip.label)]
  
  if (length(tips_matching_taxon) == 0) {
    # No tips matching the taxon, set to NA
    node_num[i] <- NA
    tip_num[i] <- NA
  } else if (length(tips_matching_taxon) < 5) {
    # Less than 5 tips matching the taxon, set to NA
    node_num[i] <- NA
    tip_num[i] <- length(tips_matching_taxon)
  } else {
    # Find MRCA of the tips matching the taxon
    nodem <- try(ggtree::MRCA(tree, tips_matching_taxon))
    
    if (inherits(nodem, "try-error")) {
      node_num[i] <- NA
      tip_num[i] <- NA
    } else {
      node_num[i] <- nodem
      tip_num[i] <- length(tips_matching_taxon)
    }
  }
}

# Creating a dataframe to reference for the final base ggtree plot
dd3 <- data.frame(phylum_names, node_num, tip_num)
dd3$node_num <- as.numeric(dd3$node_num)
dd3$Phylum <- sapply(strsplit(dd3$phylum_names, "--"), "[[", 2)
dd3 <- dd3[!is.na(dd3$node_num), ]
dd3 <- dd3[order(dd3$node_num), ]

## f. Create Base ggtree plot ----
p_no_highlight <- ggtree::ggtree(tree, layout = "fan", open.angle= 100 )
p_highlight <- p_no_highlight +
  geom_hilight(data=dd3,
               mapping=aes(fill=Phylum, node=node_num),
               alpha=0.2,
               extend=0.5) +
  scale_fill_manual(values=rainbow(nrow(dd3))[c(2,1,4,3)])

p_highlight

#Significant & Non-significant on plot
p1 <- (ggtree::ggtree(tree, layout = "fan", open.angle= 20)) %<+% psi_dd + 
    geom_point(aes(color=Direction_sig_only,size=abs(Mean)), alpha=1) + 
    scale_colour_manual(values=c("darkgreen","red"),na.translate = F)

p2 <- (ggtree::ggtree(tree, layout = "fan", open.angle= 20)) %<+% beta_dd + 
    geom_point(aes(color=n_sig, size = n_sig), alpha=1) + #, size = Mean , size = 4
    scale_color_distiller(palette = "Reds", na.value = NA, direction = 1) + #, name = "Significant Individual PFAS"
    scale_size_continuous(range = c(1,5)) 

#Load all results for each individual pfas
betas.ind <- betas %>% 
  filter(Component=="Means") %>% 
  select(Taxa,Exposure,Mean)
betas.ind <- betas.ind %>% 
  filter(grepl("s__",Taxa))
betas.ind.wide <- betas.ind %>% 
  pivot_wider(names_from = Exposure,values_from = Mean)
betas.ind.wide <- data.frame(betas.ind.wide)
rownames(betas.ind.wide) <- betas.ind.wide$Taxa
betas.ind.wide <- betas.ind.wide %>% 
  select(-Taxa)
betas.ind.wide <- data.matrix(betas.ind.wide)
# mycol <- colorpanel(n=100, low="blue",mid="white",high="red")
# heatmap_res <- heatmap.2(x=betas.ind.wide,
#                          # col=mycol,
#                          key.title="Color Key",
#                          key.xlab="",
#                          key.ylab="",
#                          Rowv=TRUE,
#                          Colv=TRUE,
#                          # cellnote=pmatrix,
#                          cexRow=.75,
#                          cexCol=1,
#                          notecex=0.5,
#                          notecol="black",
#                          trace="none",
#                          breaks = seq(-0.25,0.25,length.out = 201),
#                          margins = c(11,11))


pfna <- data.frame(betas.ind.wide) %>% 
  select(pfna)
pfpes <- data.frame(betas.ind.wide) %>% 
  select(pfpes)
pfos <- data.frame(betas.ind.wide) %>% 
  select(pfos)
pfhps <- data.frame(betas.ind.wide) %>% 
  select(pfhps)
pfoa <- data.frame(betas.ind.wide) %>% 
  select(pfoa)
pfhxs <- data.frame(betas.ind.wide) %>% 
  select(pfhxs)
pfda <- data.frame(betas.ind.wide) %>% 
  select(pfda)


# min_value <- min(betas.ind$Mean)
min_value <- -3
# max_value <- max(betas.ind$Mean)
max_value <- 3
#Mixtures
plot1 <- gheatmap(p1, pfna, offset = 0, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot2 <- gheatmap(plot1, pfpes, offset = 0.6, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot3 <- gheatmap(plot2, pfos, offset = 1.2, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot4 <- gheatmap(plot3, pfhps, offset = 1.8, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot5 <- gheatmap(plot4, pfoa, offset = 2.4, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot6 <- gheatmap(plot5, pfhxs, offset = 3.0, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot7 <- gheatmap(plot6, pfda, offset = 3.6, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
mixture.plot <- plot7+ theme(legend.position = "none") 
#Individuals

plot8 <- gheatmap(p2, pfna, offset = 0, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot9 <- gheatmap(plot8, pfpes, offset = 0.6, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot10 <- gheatmap(plot9 , pfos, offset = 1.2, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot11 <- gheatmap(plot10, pfhps, offset = 1.8, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot12 <- gheatmap(plot11, pfoa, offset = 2.4, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot13 <- gheatmap(plot12, pfhxs, offset = 3.0, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot14 <- gheatmap(plot13, pfda, offset = 3.6, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = c("#006400","darkgreen"), mid = c("darkgreen","white","red"), high = c("red","#8B0000"), midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
  # scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
individual.plot <- plot14+theme(legend.position = "none")

library(cowplot)
combined.plot <- plot_grid(mixture.plot,individual.plot,ncol=2)
ggsave(plot = combined.plot, "/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/2_Results/Figure 1.jpeg",width=30,height=30)
# pdf("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/2_Results/Results_Dendrograms_03_26.pdf",
#     width=30,height=20)
# ggarrange(mixture.plot,individual.plot,nrow=1,ncol = 2,labels=c("A.","B."))
# dev.off()

# ggsave(plot=plot1,"Legends2.jpeg")
names <- ggtree(tree, layout = "fan", open.angle= 20) %<+%  psi_dd + 
  geom_point(aes(color=Direction_sig_only,size=abs(Mean)), alpha=1) + 
  scale_colour_manual(values=c("darkgreen","red"),na.translate = F)+
  geom_tiplab()

ggsave(plot=names,"Names of sig points.jpeg",width=40,height=40)

# +theme_tree()+
#   # geom_label(aes(label = label))
#   geom_label(aes(label = label, fill = ifelse(label %in% positive, "Significant Positive",
#                                               ifelse(label %in% negative, "Significant Negative","Non-Significant"))),size=6,family="Times")+
#   scale_fill_manual(values = c("Significant Positive" = "#CC0066","Significant Negative" = "springgreen4", "Non-Significant" = "lightgray"),name = NULL)+
#   theme(legend.position = "none")

 
tree$tip.label <- str_replace(tree$tip.label,".*g__","")
# 
# 
# 
# p1+geom_tiplab()
# which(tree$tip.label %in% sig.only)
# # Create the bar chart
# bar_chart <- ggplot(data = psi_dd, aes(x = Taxa, y = abs(Mean))) +
#   geom_bar(stat = "identity",aes(fill=Mean)) +  # Customize fill color as needed
#   theme_void()+  # Remove background and gridlines
#   scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50", limits = c(-3.0, 3.0))
# 
# 
# p3 <- (ggtree::ggtree(tree, layout = "fan", open.angle= 20)) %<+% psi_dd +
#   geom_bar(data = psi_dd, stat = "identity",aes(x = Taxa, y = abs(Mean), fill=Mean))+
#   theme_void()+  # Remove background and gridlines
#   scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50", limits = c(-3.0, 3.0))
# 
# (ggtree::ggtree(tree, layout = "fan", open.angle= 20))%<+% psi_dd+
#   ggplot(data = psi_dd, aes(x = Taxa, y = abs(Mean))) +
#   geom_bar(stat = "identity",aes(fill=Mean)) +  # Customize fill color as needed
#   theme_void()+  # Remove background and gridlines
#   scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50", limits = c(-3.0, 3.0))
# 
# 
# ggtree(tree) + bar_chart 
# 
# bar_chart <- ggplot(data = psi_dd, aes(x = Taxa, y = abs(Mean))) +
#   geom_bar(stat = "identity",aes(fill=Mean)) +  # Customize fill color as needed
#   theme_void()+  # Remove background and gridlines
#   scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50", limits = c(-3.0, 3.0))+
#   coord_flip()
# 
# # library(ggtree)
# # library(phyloseq)
# # library(dplyr)
# 
# combined_plot <- ggtree(tree) + 
#   bar_chart 
# 
# # Wrap the combined plot in a circle
# combined_plot +
#   coord_polar(theta = "y")
# 
# 
# bar_chart +coord_polar(theta = "y")
# ggtree(tree)
# 
# tree$tip.label
# 

# # Cretrees# # Create the heatmap plot
# heatmap_plot <- gheatmap(p1, pfna, offset = 0, width = 0.1, colnames_offset_y = -5) +
#   scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0, na.value = "grey50", limits = c(min_value, max_value))
# 
# # Combine the plots
# combined_plot <- heatmap_plot + annotation_custom(ggplotGrob(bar_chart), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
# 
# # Display the combined plot
# print(combined_plot)
















#Significant
colnames(table) <-  str_replace(colnames(table),"_p_", "--p_")
colnames(table) <-  str_replace(colnames(table),"_c_", "--c_")
colnames(table) <-  str_replace(colnames(table),"_o_", "--o_")
colnames(table) <-  str_replace(colnames(table),"_f_", "--f_")
colnames(table) <-  str_replace(colnames(table),"_g_", "--g_")
colnames(table) <-  str_replace(colnames(table),"_s_", "--s_")

table <- table[sig.only]

taxa_tab2 <- as.data.frame(t(table))

# Change rownames (taxa names) to column
taxa_tab2 <- tibble::rownames_to_column(taxa_tab2, "V1")

# 2. Clean Taxa Table Names (JG Comment: There was an issue with this code) ----
species_names_og <- taxa_tab2[,1]
phylogenic_tree_names_df <- data.frame(species_names_og) |> 
  tidylog::mutate(
    species_names_full = str_replace(
      species_names_og, "_p_", "--p_") |>
      str_replace("_c_", "--c_") |>
      str_replace("_o_", "--o_") |>
      str_replace("_f_", "--f_") |>
      str_replace("_g_", "--g_") |>
      str_replace("_s_", "--s_"), 
    genus_name_full  = sapply(strsplit(species_names_full,"--s"),"[[",1),
    order_name_full  = sapply(strsplit(species_names_full,"--f"),"[[",1),
    family_name_full = sapply(strsplit(species_names_full,"--g"),"[[",1),
    class_name_full  = sapply(strsplit(species_names_full,"--o"),"[[",1),
    phylum_name_full = sapply(strsplit(species_names_full,"--c"),"[[",1),
    species_names = sapply(strsplit(species_names_full,"--"),"[[",7), 
    genus_name    = sapply(strsplit(species_names_full,"--"),"[[",6),
    family_name   = sapply(strsplit(species_names_full,"--"),"[[",5),
    order_name    = sapply(strsplit(species_names_full,"--"),"[[",4),
    class_name    = sapply(strsplit(species_names_full,"--"),"[[",3),
    phylum_name   = sapply(strsplit(species_names_full,"--"),"[[",2), 
    domain_name   = sapply(strsplit(species_names_full,"--"),"[[",1))

# get data frame of individual names only (not used later, only for reference)
individual_names_only <- phylogenic_tree_names_df %>%
  dplyr::select(species_names,genus_name,
                order_name,family_name,class_name,phylum_name)

# in new data, -- separates phylum, __ separates the level identifier from the name 
species_names <- phylogenic_tree_names_df$species_names_full
genus_names   <- unique(phylogenic_tree_names_df$genus_name_full)
family_names  <- unique(phylogenic_tree_names_df$family_name_full)
order_names   <- unique(phylogenic_tree_names_df$order_name_full)
class_names   <- unique(phylogenic_tree_names_df$class_name_full)
phylum_names  <- unique(phylogenic_tree_names_df$phylum_name_full)


## a. Defining taxa level names (JG- the issue above was causing issues here) ----
names_all <- c(species_names,phylum_names,class_names,
               order_names,family_names,genus_names) 
phylo_level <- sapply(strsplit(names_all,"--"),length) 

# Sanity check- are there the correct number of species:phylum
table(phylo_level)

## b. Defining branch tree text ----
# Create an empty list 
l1=list()  

# For testing:
# n=7
# m=names_all[which(phylo_level==n-1)][39]
# table(phylo_level)

## c. For loop to get Newick format tree structure -----
for (n in 7:3){
  l2=list()
  # For genus level vars:
  if(n==7){
    for (m in names_all[which(phylo_level==n-1)]){
      names2=names_all[phylo_level==n]
      l1[m]=paste0("(",paste(names2[grep(paste0(m, "-"),names2)],collapse=","),")",m)
    }
  } else{
    for (m in names_all[which(phylo_level==n-1)]){
      names2=names_all[phylo_level==n]
      l2[m]=paste0("(",paste(l1[names2[grep(m,names2)]],collapse=","),")",m)
    }
    l1=l2
  }
  print(paste0(n, ": ", length(l1), "; ", length(l2)))
}

## d. Defining parameters for branch tree text ----
# This code requires Newick format 
l3 <- paste0("(",paste(l1,collapse=","),")Bacteria;")

## e. Plot Base  ----
tree2 <- ape::read.tree(text = l3)
# library(igraph)
# igraph_tree <- ape::as.igraph.phylo(tree, directed = TRUE, ignore.root = TRUE, use.labels = TRUE)
# plot(igraph_tree, layout = layout_with_kk, 
#      vertex.size = 5,
#      vertex.label = NA,
#      # edge.curved = 0.1,
#      # vertex.label.cex = 0.5, 
#      edge.arrow.size = 0.1) 

tree2$tip.label
tree2$node.label

tree2$node.label <- str_replace(tree2$node.label, ".*--", "")
tree2$node.label <- str_replace(tree2$node.label,".*?__","")
tree2$node.label <- str_replace_all(tree2$node.label,"_"," ")
tree2$tip.label <- str_replace(tree2$tip.label, "^.*?--g__", "")
tree2$tip.label <- str_replace_all(tree2$tip.label,"--s__unclassified"," sp.")
tree2$tip.label <- str_replace_all(tree2$tip.label,"_"," ")
tree2$tip.label <- str_replace_all(tree2$tip.label,"unclassified246","Genus 246")
tree2$node.label <- str_replace_all(tree2$node.label,"unclassified246","Genus 246")

sig.psi <- psi_dd %>% 
  filter(!is.na(Direction_sig_only))
sig.psi$Taxa <- str_replace(sig.psi$Taxa, "^.*?--g__", "")
sig.psi$Taxa <- str_replace(sig.psi$Taxa,"d__Bacteria--p__Firmicutes--c__Clostridia--o__Lachnospirales--f__Lachnospiraceae","Lachnospiraceae")
sig.psi$Taxa <- str_replace(sig.psi$Taxa,"--s__unclassified"," sp.")
sig.psi$Taxa <- str_replace(sig.psi$Taxa,"unclassified246","Genus 246")
sig.psi$Taxa <- str_replace_all(sig.psi$Taxa,"_"," ")
# sig.psi$Label <- sig.psi$Taxa
positive <- sig.psi %>% 
  filter(Direction_sig_only=="Positive")
negative <- sig.psi %>% 
  filter(Direction_sig_only=="Negative")
negative <- negative$Taxa
positive <- positive$Taxa

p <- ggtree(tree2) +theme_tree()+
  # geom_label(aes(label = label))
  geom_label(aes(label = label, fill = ifelse(label %in% positive, "Significant Positive",
                                              ifelse(label %in% negative, "Significant Negative","Non-Significant"))),size=6,family="Times")+
  scale_fill_manual(values = c("Significant Positive" = "#CC0066","Significant Negative" = "springgreen4", "Non-Significant" = "lightgray"),name = NULL)+
  theme(legend.position = "none")
  # theme(element_text(family="Times",face = "italic"))
p
p <- p + theme(plot.margin = margin(2, 2, 2, 2, "cm"))
pdf("Significant Dendrogram.pdf",width=30,height=20)
p 
dev.off()


# theme(element_text(family="Times",face = "italic"))

# pdf("Significant_Tree.pdf",width=100,height=100)
install.packages("htmlwidgets")
library(htmlwidgets)
# svg("Significant_Dendrogram.html")
saveWidget(p, file = "SignificantDendrogram.html")
class(p)
# dev.off()

# Read the SVG file
svg_content <- readLines("Significant_Dendrogram.svg")
writeLines("<mxfile host='app.diagrams.net' modified='2024-03-24T00:00:00.000Z' agent='5.0 (Macintosh; Intel Mac OS X 10_15_7)'><diagram name='Untitled Diagram' id='Id'><mxGraphModel><root><mxCell/></root></mxGraphModel></diagram></mxfile>", "draw.io.xml")


# ggtree(tree)
# ggtree(tree, layout="roundrect")
# ggtree(tree, layout="slanted")
# ggtree(tree, layout="ellipse")
# ggtree(tree, layout="circular")
# ggtree(tree, layout="fan", open.angle=0)
# ggtree(tree, layout="equal_angle")
# ggtree(tree, layout="daylight")
# ggtree(tree, branch.length='none')
# ggtree(tree, layout="ellipse", branch.length="none")
# ggtree(tree, branch.length='none', layout='circular')
# ggtree(tree, layout="daylight", branch.length = 'none')





# Sanity check:
paste0("Number of nodes is ", length(tree$node.label),". It should be ", sum(table(phylo_level))-table(phylo_level)[6]+1, 
       " (the number of unique genus_names, family_names, order_names, class_names, phylum_names, and domain_names.")
paste0("Number of tips is ", length(tree$tip.label) ,". It should be ", length(species_names))


## Trouble shooting more than expected nodes
# str_count(l3, "s__unclassified113,")
# x <- purrr::map(l1, ~str_count(., "s__unclassified113,") ) |> unlist()
# str_count(l3, "s__unclassified1,")
# temp_df <- data.frame(tip.label = tree$tip.label) 
# temp_df_dup <- janitor::get_dupes(temp_df, tip.label)

# Creating empty vectors for nodes and tip labels
node_num <- vector()
tip_num <- vector()

# Assigning labels to the correct locations
for (i in 1:length(phylum_names)) {
  print(i)
  tips_matching_taxon <- tree$tip.label[grep(as.character(phylum_names[i]), tree$tip.label)]
  
  if (length(tips_matching_taxon) == 0) {
    # No tips matching the taxon, set to NA
    node_num[i] <- NA
    tip_num[i] <- NA
  } else if (length(tips_matching_taxon) < 5) {
    # Less than 5 tips matching the taxon, set to NA
    node_num[i] <- NA
    tip_num[i] <- length(tips_matching_taxon)
  } else {
    # Find MRCA of the tips matching the taxon
    nodem <- try(ggtree::MRCA(tree, tips_matching_taxon))
    
    if (inherits(nodem, "try-error")) {
      node_num[i] <- NA
      tip_num[i] <- NA
    } else {
      node_num[i] <- nodem
      tip_num[i] <- length(tips_matching_taxon)
    }
  }
}

# Creating a dataframe to reference for the final base ggtree plot
dd3 <- data.frame(phylum_names, node_num, tip_num)
dd3$node_num <- as.numeric(dd3$node_num)
dd3$Phylum <- sapply(strsplit(dd3$phylum_names, "--"), "[[", 2)
dd3 <- dd3[!is.na(dd3$node_num), ]
dd3 <- dd3[order(dd3$node_num), ]

## f. Create Base ggtree plot ----
p_no_highlight <- ggtree::ggtree(tree, layout = "fan", open.angle= 100 )
p_highlight <- p_no_highlight +
  geom_hilight(data=dd3,
               mapping=aes(fill=Phylum, node=node_num),
               alpha=0.2,
               extend=0.5) +
  scale_fill_manual(values=rainbow(nrow(dd3))[c(2,1,4,3)])

p_highlight

#Significant & Non-significant on plot
p1 <- (ggtree::ggtree(tree, layout = "fan", open.angle= 20)) %<+% psi_dd + 
  geom_point(aes(color=Direction_sig_only,size=abs(Mean)), alpha=1) + 
  scale_colour_manual(values=c("darkgreen","red"),na.translate = F)

p2 <- (ggtree::ggtree(tree, layout = "fan", open.angle= 20)) %<+% beta_dd + 
  geom_point(aes(color=n_sig, size = n_sig), alpha=1) + #, size = Mean , size = 4
  scale_color_distiller(palette = "Reds", na.value = NA, direction = 1) + #, name = "Significant Individual PFAS"
  scale_size_continuous(range = c(1,5)) 

#Load all results for each individual pfas
betas.ind <- betas %>% 
  filter(Component=="Means") %>% 
  select(Taxa,Exposure,Mean)
betas.ind <- betas.ind %>% 
  filter(grepl("s__",Taxa))
betas.ind.wide <- betas.ind %>% 
  pivot_wider(names_from = Exposure,values_from = Mean)
betas.ind.wide <- data.frame(betas.ind.wide)
rownames(betas.ind.wide) <- betas.ind.wide$Taxa
betas.ind.wide <- betas.ind.wide %>% 
  select(-Taxa)
betas.ind.wide <- data.matrix(betas.ind.wide)
# mycol <- colorpanel(n=100, low="blue",mid="white",high="red")
# heatmap_res <- heatmap.2(x=betas.ind.wide,
#                          # col=mycol,
#                          key.title="Color Key",
#                          key.xlab="",
#                          key.ylab="",
#                          Rowv=TRUE,
#                          Colv=TRUE,
#                          # cellnote=pmatrix,
#                          cexRow=.75,
#                          cexCol=1,
#                          notecex=0.5,
#                          notecol="black",
#                          trace="none",
#                          breaks = seq(-0.25,0.25,length.out = 201),
#                          margins = c(11,11))


pfna <- data.frame(betas.ind.wide) %>% 
  select(pfna)
pfpes <- data.frame(betas.ind.wide) %>% 
  select(pfpes)
pfos <- data.frame(betas.ind.wide) %>% 
  select(pfos)
pfhps <- data.frame(betas.ind.wide) %>% 
  select(pfhps)
pfoa <- data.frame(betas.ind.wide) %>% 
  select(pfoa)
pfhxs <- data.frame(betas.ind.wide) %>% 
  select(pfhxs)
pfda <- data.frame(betas.ind.wide) %>% 
  select(pfda)


# min_value <- min(betas.ind$Mean)
min_value <- -1.5
# max_value <- max(betas.ind$Mean)
max_value <- 1.5
#Mixtures
plot1 <- gheatmap(p1, pfna, offset = 0, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot2 <- gheatmap(plot1, pfpes, offset = 0.6, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot3 <- gheatmap(plot2, pfos, offset = 1.2, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot4 <- gheatmap(plot3, pfhps, offset = 1.8, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot5 <- gheatmap(plot4, pfoa, offset = 2.4, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot6 <- gheatmap(plot5, pfhxs, offset = 3.0, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot7 <- gheatmap(plot6, pfda, offset = 3.6, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
mixture.plot <- plot7+ theme(legend.position = "none") 
#Individuals

plot8 <- gheatmap(p2, pfna, offset = 0, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot9 <- gheatmap(plot8, pfpes, offset = 0.6, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot10 <- gheatmap(plot9 , pfos, offset = 1.2, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot11 <- gheatmap(plot10, pfhps, offset = 1.8, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot12 <- gheatmap(plot11, pfoa, offset = 2.4, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot13 <- gheatmap(plot12, pfhxs, offset = 3.0, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
plot14 <- gheatmap(plot13, pfda, offset = 3.6, width = 0.1, colnames_offset_y = -5) +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "black", midpoint = 0, na.value = "grey50",limits = c(min_value, max_value))
individual.plot <- plot14+theme(legend.position = "none") 

# pdf("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/2_Results/Results_Dendrograms.pdf",
#     width=30,height=20)
# ggarrange(mixture.plot,individual.plot,nrow=1,ncol = 2,labels=c("A.","B."))
# dev.off()



p <- ggtree(tree) +theme_tree()+
  # geom_label(aes(label = label))+
  geom_label(aes(label = label, fill = ifelse(label %in% positive, "Significant Positive", 
                                              ifelse(label %in% negative, "Significant Negative","Non-Significant"))),size=2)+
  scale_fill_manual(values = c("Significant Positive" = "#CC0066","Significant Negative" = "#3399FF", "Non-Significant" = "lightgray"),name = NULL)


p <- ggtree(tree) +theme_tree()+
  # geom_label(aes(label = label))+
  geom_label(aes(label = label, fill = ifelse(label %in% positive, "Significant Positive", 
                                              ifelse(label %in% negative, "Significant Negative","Non-Significant"))),size=2)+
  scale_fill_manual(values = c("Significant Positive" = "#CC0066","Significant Negative" = "#3399FF", "Non-Significant" = "lightgray"),name = NULL)


# 4. Creating final plots ----
library(treeio)
library(ggtree)
# tree$node.label <- paste("Node", 1:length(tree$node.label))
# ggtree(tree, aes(label = node.label))
# ggtree(tree) + geom_nodelab(aes(label = node.label), geom = "label")
# tree$node.label <- paste("Node", 1:length(tree$node.label))

# Print the tree with labels
# p <- ggtree(tree)

# Add node labels
# p <- p + geom_text2(aes(subset = isTip, label = ifelse(isTip, tree$tip.label, node.label)), vjust = 0)

label <- c(tree$node.label,tree$tip.label)
label <- str_replace(label,"Oscillospiraceae--g__UCG_005","Oscillospiraceae--g__Oscillospiraceae_UCG_005")
label <- str_replace(label,"--s__unclassified","_sp.")
label <- str_remove(label,".+__")
label <- str_replace_all(label,"_"," ")
tree$node.label
tree$tip.label
tree$node.label <- c(tree$node.label)
tree$node.label <- str_replace(tree$node.label,"Oscillospiraceae--g__UCG_005","Oscillospiraceae--g__Oscillospiraceae_UCG_005")
# tree$node.label <- str_replace(tree$node.label,"--s__unclassified","_sp.")
tree$node.label <- str_remove(tree$node.label,".+__")
tree$node.label <- str_replace_all(tree$node.label,"_"," ")
tree$tip.label <- c(tree$tip.label)
tree$tip.label <- str_replace(tree$tip.label,"Oscillospiraceae--g__UCG_005","Oscillospiraceae--g__Oscillospiraceae_UCG_005")
tree$tip.label <- str_replace(tree$tip.label,"--s__unclassified","_sp.")
tree$tip.label <- str_remove(tree$tip.label,".+__")
tree$tip.label <- str_replace_all(tree$tip.label,"_"," ")

ggtree(tree) + geom_tiplab() +
  # geom_label(aes(x=branch, label=), fill='lightgreen') +
  geom_label(aes(label=label), fill='orange') 
# geom_text(aes(label=tree$node.label), hjust=-.5)
# pdf("TreeSigResults.pdf",width=30,height=30)
# ggtree(tree) +
#     geom_label(aes(label=label), fill='#CC0066') 
# dev.off()

#Get significant only names
sig <- psi_dd %>% 
  filter(Direction_sig_only=="Negative")
sig.pos <- psi_dd %>% 
  filter(Direction_sig_only=="Positive")
sig$Taxa
sig$Taxa <- str_replace(sig$Taxa,"Oscillospiraceae--g__UCG_005","Oscillospiraceae--g__Oscillospiraceae_UCG_005")
sig$Taxa <- str_replace(sig$Taxa,"--s__unclassified","_sp.")
sig$Taxa <- str_remove(sig$Taxa,".+__")
sig$Taxa <- str_replace_all(sig$Taxa,"_"," ")
sig.pos$Taxa
sig.pos$Taxa <- str_replace(sig.pos$Taxa,"Oscillospiraceae--g__UCG_005","Oscillospiraceae--g__Oscillospiraceae_UCG_005")
sig.pos$Taxa <- str_replace(sig.pos$Taxa,"--s__unclassified","_sp.")
sig.pos$Taxa <- str_remove(sig.pos$Taxa,".+__")
sig.pos$Taxa <- str_replace_all(sig.pos$Taxa,"_"," ")

# colored_labels <- data.frame(label = sig$Taxa,
#                            color = c("#CC0066"))
positive <- sig.pos$Taxa
pos.cor <- sig.pos$Mean
negative <- sig$Taxa
neg.cor <- sig$Mean
effect_estimates <- data.frame(label = c(positive,negative),
                               effect = c(pos.cor,neg.cor))  # Adjust with your actual effect estimates


# tree$color <- label_colors 
# Merge the label_colors data frame with the tree's node labels
# tree_data <- merge(tree$data, label_colors, by = "label", all.x = TRUE)

# Plot the tree with labels and custom colors
p <- ggtree(tree) +theme_tree()+
  # geom_label(aes(label = label))+
  geom_label(aes(label = label, fill = ifelse(label %in% positive, "Significant Positive", 
                                              ifelse(label %in% negative, "Significant Negative","Non-Significant"))),size=2)+
  scale_fill_manual(values = c("Significant Positive" = "#CC0066","Significant Negative" = "#3399FF", "Non-Significant" = "lightgray"),name = NULL)


# pdf("TreeSigResults.pdf",width=20,height=10)   
# p
# dev.off()

# ggsave("plot.pdf", p, width = 10, height = 8, units = "in", device = cairo_pdf, 
#        bg = "white", limitsize = FALSE, dpi = 300,  # Additional options for PDF output
#        expand = c(0.5, 0.5)) # Adjust the margin values as needed

# 
# gplots <- list()
# # Mixture effects, significant only
# (gplots[[1]]=  p %<+% psi_dd + 
#         geom_point(aes(color=Direction_sig_only, size = abs(correlation)), alpha=1) + 
#         scale_colour_manual(values=c("blue","red"),na.translate = F))
# gplots[[1]]


