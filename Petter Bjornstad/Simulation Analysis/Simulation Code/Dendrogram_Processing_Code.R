#Dendrogram Code
#1. Rename Simulation Variables ----
# #Load formatted object
# temp_file <- tempfile(fileext = ".RDS") 
# s3$download_file(bucket,"FormattedObject.RDS", temp_file)
# FormattedObject <- readRDS(temp_file)
# Table <- data.frame(FormattedObject$Table)
# 
# Y <- Table[1:207]
# Y <- Y[-149]
# N <- nrow(Y)
# P.s <- ncol(Y)

# Use your Z-matrices to get names at each level
species_names <- rownames(Z.s.g)
genus_names   <- colnames(Z.s.g)
family_names  <- colnames(Z.g.f)
order_names   <- colnames(Z.f.o)
class_names   <- colnames(Z.o.c)
phylum_names  <- colnames(Z.c.p)

# Build Newick tree from Z-matrices
# Species -> Genus
genus_list <- list()
for(g in genus_names){
  children <- species_names[Z.s.g[, g] == 1]
  genus_list[[g]] <- paste0("(", paste(children, collapse=","), ")", g)
}

# Genus -> Family
family_list <- list()
for(f in family_names){
  children <- genus_names[Z.g.f[, f] == 1]
  family_list[[f]] <- paste0("(", paste(unlist(genus_list[children]), collapse=","), ")", f)
}

# Family -> Order
order_list <- list()
for(o in order_names){
  children <- family_names[Z.f.o[, o] == 1]
  order_list[[o]] <- paste0("(", paste(unlist(family_list[children]), collapse=","), ")", o)
}

# Order -> Class
class_list <- list()
for(cl in class_names){
  children <- order_names[Z.o.c[, cl] == 1]
  class_list[[cl]] <- paste0("(", paste(unlist(order_list[children]), collapse=","), ")", cl)
}

# Class -> Phylum
phylum_list <- list()
for(p in phylum_names){
  children <- class_names[Z.c.p[, p] == 1]
  phylum_list[[p]] <- paste0("(", paste(unlist(class_list[children]), collapse=","), ")", p)
}

# Full Newick string
newick <- paste0("(", paste(unlist(phylum_list), collapse=","), ")root;")

# Build and plot tree
tree <- ape::read.tree(text = newick)
ggtree::ggtree(tree, layout = "fan", open.angle = 40)

#####################

taxa_tab2 <- as.data.frame(t(Y))
taxa_tab2 <- tibble::rownames_to_column(taxa_tab2, "Taxa")
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

all.taxa <- phylogenic_tree_names_df
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

# Plot Base Trees
tree <- ape::read.tree(text = l3)
(ggtree::ggtree(tree, layout = "fan", open.angle= 40))
labels <- data.frame(tree$tip.label)
labels <- labels %>% 
  filter(grepl("p__Bacteroidota",tree.tip.label))
phy.labels <- labels$tree.tip.label
# pdf("GetTaxaNames.pdf",width=75,height=75)
# p <- (ggtree::ggtree(tree, layout = "fan", open.angle= 40)) +
#   geom_tiplab(aes(label=label, color=ifelse(label %in% phy.labels, "highlight", "normal")))
# dev.off()


#Causal scenarios 
#Figure out the causal taxa in each scenario
#Load simulation scenarios
temp_file <- tempfile(fileext = ".rds") 
s3$download_file(bucket,"Simulation_Scenarios_03_12_25.rds", temp_file)
simulation.parameters <- readRDS(temp_file)
# simulation.parameters <- readRDS(fs::path(dir.code,"Simulation_Scenarios_03_12_25.rds"))

#Load all simulatin parameters
temp_file <- tempfile(fileext = ".csv") 
s3$download_file(bucket,"Simulation_Parameters.csv", temp_file)
sim.par.all <- read.csv(temp_file)
# sim.par.all <- read.csv(fs::path(dir.code,"Simulation_Parameters.csv")) %>% 
sim.par.all <- sim.par.all %>% 
  dplyr::select(-X) %>% 
  dplyr::select(c("P.e.causal","P.s.causal.all","N","OR.exposure")) %>% 
  dplyr::rename(P.s.scenario=P.s.causal.all)

#Set causal scenarios
scenario1 <- simulation.parameters[[sim.par.all$P.s.scenario[1]]]
scenario2 <- simulation.parameters[[sim.par.all$P.s.scenario[2]]]
scenario3 <- paste0("species",simulation.parameters[[sim.par.all$P.s.scenario[3]]])
scenario4 <- paste0("species",simulation.parameters[[sim.par.all$P.s.scenario[4]]])
scenario5 <- paste0("species",simulation.parameters[[sim.par.all$P.s.scenario[5]]])
scenario6 <- paste0("species",simulation.parameters[[sim.par.all$P.s.scenario[6]]])
scenario7 <- paste0("species",simulation.parameters[[sim.par.all$P.s.scenario[7]]])
scenario8 <- paste0("species",simulation.parameters[[sim.par.all$P.s.scenario[8]]])
scenario9 <- paste0("species",simulation.parameters[[sim.par.all$P.s.scenario[9]]])


#Make Second base tree - only causal 
taxa_tab2 <- as.data.frame(t(Y))
# taxa_tab2 <- as.data.frame(t(Table[1:206]))
taxa_tab2 <- tibble::rownames_to_column(taxa_tab2, "Taxa")
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

#Filter to only significant causal species
phylogenic_tree_names_df <- phylogenic_tree_names_df %>% 
  filter(species_names_full %in% scenario8)

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

# Plot Base Trees
tree2 <- ape::read.tree(text = l3)
(ggtree::ggtree(tree2, layout = "fan", open.angle= 240))

#Scenarios
# scenario.all <- scenario8
# p <- ggtree(tree2, layout = "fan", open.angle= 240)
# ggtree(tree2, layout = "fan", open.angle= 240)  +theme_tree()+
#   geom_label(aes(label = label, fill = ifelse(label %in% scenario, "Significant","Non-Significant")),size=2)+
#   scale_fill_manual(values = c("Significant" = "#CC0066","Non-Significant" = "#3399FF"),name = NULL)

#Scenario 2
# scenario <- scenario2
scenario.fxn <- function(scenario){
  causal.or <- 1.5
  noncausal.or <- 1.0
  #Species
  dat <- data.frame(species=phylogenic_tree_names_df$species_names_full) 
  dat <- dat %>% 
    mutate(ExpVal= ifelse(species %in% scenario,log(causal.or),log(noncausal.or)))
  #Genus
  genus.names <- data.frame(genus = phylogenic_tree_names_df$genus_name_full,species=phylogenic_tree_names_df$species_names_full)
  genus.names <- tidylog::left_join(dat,genus.names,by="species")
  genus.names <- genus.names %>% 
    select(-species) %>% 
    group_by(genus) %>% 
    summarize(ExpVal=mean(ExpVal)) %>% 
    ungroup() %>% 
    rename(Taxa=genus)
  #Family
  fam.names <- genus.names
  fam.names <- fam.names %>% 
    rename(family=Taxa)
  fam.names$family <- sub("--g__.*$","",fam.names$family)
  fam.names <- fam.names %>% 
    group_by(family) %>% 
    summarize(ExpVal=mean(ExpVal)) %>% 
    ungroup()%>% 
    rename(Taxa=family)
  #Order
  ord.names <- fam.names
  ord.names <- ord.names %>% 
    rename(order=Taxa)
  ord.names$order <- sub("--f__.*$","",ord.names$order)
  ord.names <- ord.names %>% 
    group_by(order) %>% 
    summarize(ExpVal=mean(ExpVal)) %>% 
    ungroup()%>% 
    rename(Taxa=order)
  #Class
  class.names <- ord.names
  class.names <- class.names %>% 
    rename(class=Taxa)
  class.names$class <- sub("--o__.*$","",class.names$class)
  class.names <- class.names %>% 
    group_by(class) %>% 
    summarize(ExpVal=mean(ExpVal)) %>% 
    ungroup()%>% 
    rename(Taxa=class)
  #Phylum
  phy.names <- class.names
  phy.names <- phy.names %>% 
    rename(phylum=Taxa)
  phy.names$phylum <- sub("--c__.*$","",phy.names$phylum)
  phy.names <- phy.names %>% 
    group_by(phylum) %>% 
    summarize(ExpVal=mean(ExpVal)) %>% 
    ungroup()%>% 
    rename(Taxa=phylum)
  
  dat <- dat %>% 
    rename(Taxa=species)
  expected.data <- rbind(dat,genus.names)
  expected.data <- rbind(expected.data,fam.names)
  expected.data <- rbind(expected.data,ord.names)
  expected.data <- rbind(expected.data,class.names)
  expected.data <- rbind(expected.data,phy.names)
  return(expected.data)
}

#Scenario 1
scenario1.data <- scenario.fxn(scenario1)
#Scenario 2
scenario2.data <- scenario.fxn(scenario2)
#Scenario 3
scenario3.data <- scenario.fxn(scenario3)
#Scenario 4
scenario4.data <- scenario.fxn(scenario4)
#Scenario 5
scenario5.data <- scenario.fxn(scenario5)
#Scenario 6
scenario6.data <- scenario.fxn(scenario6)
#Scenario 7
scenario7.data <- scenario.fxn(scenario7)
#Scenario 8
scenario8.data <- scenario.fxn(scenario8)
# 
# #Plot scenarios
min <- 0
max <- log(1.5)


# # Define custom breaks for size scale based on your condition
custom_sizes <- ifelse(max== 0, c(0.05, 0.075, 0.1), c(1, 2, 3))
# 
# #Scenario 1
# #FUll tree

scen1a <- ggtree(tree, layout = "fan", open.angle= 40)  %<+% scenario1.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient(low = "#CCCCCC", high = "#CC0066", limits = c(min_val, max_val)) +
  # scale_color_gradient2(low = "gray", mid = "orange", high = "#CC0066",
  #                       midpoint = max/2, limits = c(min,max), na.value = NA)+
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")

scen1b <- ggtree(tree2,layout="fan",open.angle=225) %<+% scenario1.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient2(low = "gray", mid = "orange", high = "#CC0066",
  #                       midpoint = max/2, limits = c(min,max), na.value = NA)+
  # scale_color_gradient(low = "#3399FF", high = "#CC0066",limits=c(0,0.5))
  # scale_color_gradient(low = "#3399FF", high = "#CC0066", limits = c(min_val, max_val)) +
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")

#Scenario2 - 3 of 6 causal exposures
# scenario2.data <- scenario2.data %>% 
#   mutate(ExpValExp = ifelse(ExpVal!=0,mean(ExpVal)))
# # mean(scenario2.data$ExpVal*1  +scenario2.data$ExpVal*1 +scenario2.data$ExpVal*1 +scenario2.data$ExpVal*0+scenario2.data$ExpVal*0+scenario2.data$ExpVal*0)
# 
scen2 <- ggtree(tree2,layout="fan",open.angle=225) %<+% scenario2.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient2(low = "gray", mid = "orange", high = "#CC0066",
  #                       midpoint = max/2, limits = c(min,max), na.value = NA)+
  # scale_color_gradient(low = "#3399FF", high = "#CC0066",limits=c(0,0.5))+
  # scale_color_gradient(low = "black", high = "lightskyblue", limits = c(min_val, max_val)) +
  # scale_color_gradient(low = "#3399FF", high = "#CC0066", limits = c(min_val, max_val)) +
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")
# scale_size_continuous(range = c(0,0.5), limits = c(0,0.5))
# scale_colour_gradient2(values=c(low="#3399FF",mid = "white",high="#CC0066"),na.value = F)


#Scenario3
scen3 <- ggtree(tree2,layout="fan",open.angle=225) %<+% scenario3.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient2(low = "gray", mid = "orange", high = "#CC0066",
  #                       midpoint = max/2, limits = c(min,max), na.value = NA)+
  # scale_color_gradient(low = "#3399FF", high = "#CC0066",limits=c(0,0.5))
  # scale_color_gradient(low = "#3399FF", high = "#CC0066", limits = c(min_val, max_val)) +
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")


#Scenario4
scen4 <- ggtree(tree2,layout="fan",open.angle=225) %<+% scenario4.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient2(low = "gray", mid = "orange", high = "#CC0066",
  #                       midpoint = max/2, limits = c(min,max), na.value = NA)+
  # scale_color_gradient(low = "#3399FF", high = "#CC0066",limits=c(0,0.5))
  # scale_color_gradient(low = "#3399FF", high = "#CC0066", limits = c(min_val, max_val)) +
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")

#Scenario5
scen5 <- ggtree(tree2,layout="fan",open.angle=225) %<+% scenario5.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient(low = "#3399FF", high = "#CC0066",limits=c(0,0.5))
  # scale_color_gradient(low = "#3399FF", high = "#CC0066", limits = c(min_val, max_val)) +
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")

#Scenario6
scen6 <- ggtree(tree2,layout="fan",open.angle=225) %<+% scenario6.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient(low = "#3399FF", high = "#CC0066",limits=c(0,0.5))
  # scale_color_gradient(low = "#3399FF", high = "#CC0066", limits = c(min_val, max_val)) +
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")

#Scenario7
scen7 <- ggtree(tree2,layout="fan",open.angle=225) %<+% scenario7.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient(low = "#3399FF", high = "#CC0066",limits=c(0,0.5))
  # scale_color_gradient(low = "#3399FF", high = "#CC0066", limits = c(min_val, max_val)) +
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")

#Scenario8
scen8 <- ggtree(tree2,layout="fan",open.angle=225) %<+% scenario8.data +
  geom_point(aes(color=ExpVal, size = ExpVal), alpha=1)+
  scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066", midpoint = max/2, limits = c(min, max), na.value = NA) +
  # scale_color_gradient(low = "#3399FF", high = "#CC0066",limits=c(0,0.5))
  # scale_color_gradient(low = "#3399FF", high = "#CC0066", limits = c(min_val, max_val)) +
  scale_size_continuous(breaks = custom_sizes, limits = c(0, max(custom_sizes)))+
  theme(legend.position = "none")

# pdf(fs::path(here::here() %>% dirname(),"2_Results","Scenario_Plots_03_31_Figure 1.pdf"),height=10,width=20)
# # ggarrange(scen1a,scen2,scen3,scen4,scen5,scen6,scen7,scen8,ncol=4,nrow=2,labels = c("A.", "B.","C.","D.","E.","F.","G.","H."))+
# #   theme(family="Times")
# plot_grid(scen1b,scen2,scen3,scen4,scen5,scen6,scen7,scen8,ncol=4,labels = c("A. Scenario 1", "B. Scenarios 2 & 3","C. Scenarios 4 & 5","D. Scenarios 6 & 7","E. Scenarios 8 & 9","F. Scenarios 10 & 11","G. Scenarios 12 & 13","H. Scenarios 14 & 15"))+
# theme(element_text(family="Times"))
# dev.off()




# # P.e <- 7
# P.e.causal <- c(0,3,6)
# N <- 100
# OR.exposure <- c(1,1.5)
# P.s.scenario <- c(1:8)
# Iterations <- 100
# Iter <- 1:100
# sim.par <- expand.grid(P.e.causal = P.e.causal, P.s.scenario = P.s.scenario, N = N, OR.exposure = OR.exposure) #36 combinations
# sim.par <- sim.par %>%
#   filter((OR.exposure==1 & P.s.scenario==1 & P.e.causal==0) | (OR.exposure==1.5 & P.e.causal!=0) & OR.exposure==1.5 & P.s.scenario!=1) #20
# # saveRDS(sim.par,"/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Simulation_Parameters_03_21.rds")
# # saveRDS(sim.par,fs::path(here::here(),"Formatted Data","Simulation_Parameters_all_04_14.rds"))

# for (i in 1:15) {
#   sim.par.new <- sim.par[i,]
#   sim.par.new.a <- replicate(Iterations, sim.par.new, simplify = FALSE)
#   sim.par.new <- do.call(rbind, sim.par.new.a)
#   sim.par.new <- sim.par.new %>%
#     mutate(Iterations=Iter)
#   saveRDS(sim.par.new,paste0("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Simulation_Parameter_V",i,"_04_14.rds"))
# }

# rm(sim.par.new,sim.par.new.a,tree,tree2)

#Create taxa key
species.names <- rownames(Z.s.g)
species.key <- data.frame(species=species.names,sim.num=paste0("species",c(1:length(species.names))))
genus.names <- colnames(Z.s.g)
genus.key <- data.frame(genus=genus.names,sim.num=paste0("genus",c(1:length(genus.names))))
family.names <- colnames(Z.g.f)
family.key <- data.frame(family=family.names,sim.num=paste0("family",c(1:length(family.names))))
order.names <- colnames(Z.f.o)
order.key <- data.frame(order=order.names,sim.num=paste0("order",c(1:length(order.names))))
class.names <- colnames(Z.o.c)
class.key <- data.frame(class=class.names,sim.num=paste0("class",c(1:length(class.names))))
phylum.names <- colnames(Z.c.p)
phylum.key <- data.frame(phylum=phylum.names,sim.num=paste0("phylum",c(1:length(phylum.names))))

#Rename matrices
#Z.s.g
rownames(Z.s.g) <- species.key$sim.num
colnames(Z.s.g) <- genus.key$sim.num
#Z.g.f
rownames(Z.g.f) <- genus.key$sim.num
colnames(Z.g.f) <- family.key$sim.num
#Z.f.o
rownames(Z.f.o) <- family.key$sim.num
colnames(Z.f.o) <- order.key$sim.num
#Z.o.c
rownames(Z.o.c) <- order.key$sim.num
colnames(Z.o.c) <- class.key$sim.num
#Z.c.p
rownames(Z.c.p) <- class.key$sim.num
colnames(Z.c.p) <- phylum.key$sim.num

#Y - matrix
colnames(Y) <- species.key$sim.num
# saveRDS(Y,fs::path(here::here(),"Y_matrix.rds"))
# 
# #Save Matrices for Simulation
# # saveRDS(Z.s.g,fs::path(here::here(),"Formatted Data","Z.s.g.RDS"))
# # saveRDS(Z.g.f,fs::path(here::here(),"Formatted Data","Z.g.f.RDS"))
# # saveRDS(Z.f.o,fs::path(here::here(),"Formatted Data","Z.f.o.RDS"))
# # saveRDS(Z.o.c,fs::path(here::here(),"Formatted Data","Z.o.c.RDS"))
# # saveRDS(Z.c.p,fs::path(here::here(),"Formatted Data","Z.c.p.RDS"))
# 
#Simulation Scenarios Format
#Species Key
species.key <- species.key %>%
  mutate(species.new=str_replace_all(species,"_p__","--p__"),
         species.new=str_replace_all(species.new,"_c__","--c__"),
         species.new=str_replace_all(species.new,"_o__","--o__"),
         species.new=str_replace_all(species.new,"_f__","--f__"),
         species.new=str_replace_all(species.new,"_g__","--g__"),
         species.new=str_replace_all(species.new,"_s__","--s__"))

scenario1.index <- scenario1
scenario2.index <- which(species.key$species.new %in% scenario2)
scenario3.index <- which(species.key$species.new %in% scenario3)
scenario4.index <- which(species.key$species.new %in% scenario4)
scenario5.index <- which(species.key$species.new %in% scenario5)
scenario6.index <- which(species.key$species.new %in% scenario6)
scenario7.index <- which(species.key$species.new %in% scenario7)
scenario8.index <- which(species.key$species.new %in% scenario8)

simulation.scenarios <- list(scenario1.index,scenario2.index,scenario3.index,scenario4.index,scenario5.index,scenario6.index,
                             scenario7.index,scenario8.index)
# # saveRDS(simulation.scenarios,"/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Simulation_Scenarios_03_23.rds")
# saveRDS(simulation.scenarios,"/Users/hhampson/USC Lab Dropbox/Chatzi Active Projects/MOVED Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Simulation_Scenarios_05_11.rds")

