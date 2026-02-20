#Clean Data & Calculate Variables for Analysis ----
#Libraries 
library(purrr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(missForest)
library(phyloseq)
library(ggpubr)
library(ggtree)
packageVersion("ggtree")
# ‘3.6.2’
library(RColorBrewer)
# install.packages("ggfun", dependencies = TRUE)
library(ggfun)
library(stringr)
library(tidyverse)
library(egg)

#Load data
FormattedObject <- readRDS(fs::path(here::here(),"Formatted Data","FormattedObject.RDS"))
Table <- data.frame(FormattedObject$Table)
# ori.taxa <- colnames(Table)[1:207]
# #Remove mitochondria
# ori.taxa <- ori.taxa[-149] #206 taxa remain
# Tab
#Check if characteristics differ by those with and without microbiome 
# demo <- readRDS("/Users/hehampson/USC Lab Dropbox/Chatzi Active Projects/MOVED Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/0_data_processing/SOL CHS Diet Data Processing/CHS/2_cleaned_data/temp/Nutrients_Data_Formatted.rds") %>% 
#   select(id,m_total_dietary_fiber_g,m_energy_kcal) %>% 
#   rename(m_fiber=m_total_dietary_fiber_g,
#          m_kcal=m_energy_kcal)
# cov <-  readRDS("/Users/hehampson/USC Lab Dropbox/Chatzi Active Projects/MOVED Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/1_analysis_ready_data/CHS MetaAir MetaChem Cleaned Redcap and Exposures Outcome Data long V3.rds") %>%
#   filter(study=="MetaAir") %>%
#   select(id,ma_age,ma_sex,ma_race_eth,ma_parental_edu,ma_bmi)
# ids <- data.frame(id=Table$id) %>%
#   mutate(id=as.numeric(id))
# demo <- demo %>% 
#   mutate(CompleteData=ifelse(id %in% ids$id,"Yes","No"))
# cov <- cov %>% 
#   mutate(CompleteData=ifelse(id %in% ids$id,"Yes","No"))
# length(which(cov$CompleteData=="Yes")) #71 with CompleteData
# length(which(cov$CompleteData=="Yes")) #103 of 156 without CompleteData
# table(cov$CompleteData)  
# t.test(ma_age~CompleteData,data=cov)
# t.test(ma_bmi~CompleteData,data=cov)
# 
# table(cov$ma_sex, cov$CompleteData)
# chisq.test(table(cov$ma_sex,cov$CompleteData))
# chisq.test(table(cov$ma_race_eth,cov$CompleteData))
# chisq.test(table(cov$ma_parental_edu,cov$CompleteData))
# 
# length(which(demo$CompleteData=="Yes")) #70 with CompleteData
# length(which(demo$CompleteData=="Yes")) #86 of 156 without CompleteData
# table(demo$CompleteData)  
# t.test(m_fiber~CompleteData,data=demo)
# t.test(m_kcal~CompleteData,data=demo)

# # Combine the datasets
# cov <- cov %>%
#   mutate(id = as.numeric(id))  # Convert id to numeric
# 
# # Join and rename columns
# data <- cov %>%
#   left_join(demo, by = c("CompleteData","id")) %>%
#   dplyr::select(CompleteData,
#     Age = ma_age,
#     Sex = ma_sex,
#     `Race/Ethnicity` = ma_race_eth,
#     `Parental Education` = ma_parental_edu,
#     BMI = ma_bmi,
#     `Fiber Intake` = m_fiber,
#     `Energy Intake` = m_kcal
#   )
# 
# # Create a formatted summary table
# summary_table <- data %>%
#   tbl_summary(
#     by = CompleteData,  # Compare by CompleteData (Yes/No)
#     statistic = list(
#       all_continuous() ~ "{mean} ({sd})",   # Mean and SD for continuous variables
#       all_categorical() ~ "{n} ({p}%)"      # Count and percentages for categorical variables
#     ),
#     missing = "no"
#   ) %>%
#   add_p() %>%  # Add p-values for comparisons
#   modify_header(label = "**Participant Characteristics**") %>%
#   modify_spanning_header(everything() ~ "**Complete Data**") %>% 
#   bold_labels()
# 
# # Print the table
# summary_table


#Load fiber data 
# "/Users/hhampson/USC Lab Dropbox/Chatzi Active Projects/MOVED Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis"
fib <- readRDS("/Users/hehampson/USC Lab Dropbox/Chatzi Active Projects/MOVED Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/0_data_processing/SOL CHS Diet Data Processing/CHS/2_cleaned_data/temp/Nutrients_Data_Formatted.rds") %>% 
  select(id,m_total_dietary_fiber_g,m_energy_kcal) %>% 
  rename(m_fiber=m_total_dietary_fiber_g,
         m_kcal=m_energy_kcal)
ids <- data.frame(id=Table$id) %>% 
  mutate(id=as.numeric(id))
fib <- tidylog::left_join(ids,fib,by="id") #70 out of 71 have energy and fiber intake
#Load other covariates for imputation - age, sex, race/ethnicity, parental education
# cov <- readRDS("/Users/hhampson/USC Lab Dropbox/Chatzi Active Projects/MOVED Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/1_analysis_ready_data/CHS MetaAir MetaChem Cleaned Redcap and Exposures Outcome Data long V3.rds") %>% 
cov <-  readRDS("/Users/hehampson/USC Lab Dropbox/Chatzi Active Projects/MOVED Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/1_analysis_ready_data/CHS MetaAir MetaChem Cleaned Redcap and Exposures Outcome Data long V3.rds") %>%
  filter(study=="MetaAir") %>% 
  select(id,ma_age,ma_sex,ma_race_eth,ma_parental_edu)
cov$id <- as.numeric(cov$id)
cov <- tidylog::left_join(ids,cov,by="id")
fib <- tidylog::full_join(fib,cov,by="id")
#Impute energy and fiber for 1 id 
fib$ma_sex <- as.factor(fib$ma_sex)
fib$ma_race_eth <- as.factor(fib$ma_race_eth)
fib$ma_parental_edu <- as.factor(fib$ma_parental_edu)
fib$id <- as.numeric(fib$id)

###i. Energy Intake ----
missForest_imputed <- data.frame(
  original = fib$m_kcal,
  imputed_missForest = missForest(fib)$ximp$m_kcal
)
fib$m_kcal_imp <- missForest_imputed$imputed_missForest

###ii. Fiber ----
missForest_imputed <- data.frame(
  original = fib$m_fiber,
  imputed_missForest = missForest(fib)$ximp$m_fiber
)
fib$m_fiber_imp <- missForest_imputed$imputed_missForest

#Filter to newly imputed variables
fib <- fib %>% 
  select(-m_fiber,-m_kcal) %>% 
  rename(m_fiber=m_fiber_imp,
         m_kcal=m_kcal_imp)

#1. Calculate Energy-Adjusted Fiber Residual ----
#a. Regress fiber (y) on energy intake (x)
mod <- lm(m_fiber ~ m_kcal,data=fib)
mean(fib$m_kcal) #Mean total energy = 1974.538
mean.energy <- data.frame(m_kcal=mean(fib$m_kcal) )
#Predicted fiber intake at mean energy
predict(mod,mean.energy) #17.48132 = predicted fiber intake at mean energy intake 
#Add constant (predicted fiber intake at mean energy intake) to make fiber residuals "intake values"
fiber.resid = resid(mod)+predict(mod,mean.energy)
fib$fib.adj <- fiber.resid 
fib$id <- as.character(fib$id)
fib <- fib %>% 
  rename(fiber=fib.adj,
         total.energy=m_kcal) %>% 
  select(id,fiber,total.energy)
#Merge into data
Table <- tidylog::left_join(Table,fib,by="id")
#save fiber
# saveRDS(fib,fs::path(here::here(),"Formatted Data","Fiber_Adjusted.RDS"))
rm(cov,fib,ids,mean.energy,missForest_imputed,mod)

#2. Calculate Z-matrices ----
#Set up for analysis ----
## Outcomes ----
Y <- Table[1:207]
Y <- Y[-149]
N <- nrow(Y)
P.s <- ncol(Y)

## Exposure ----
X <- Table %>% 
  select(pfda:pfpes)
X.q <- apply(X, 2, function(v) { cut(v, breaks=quantile(v, probs=c(0, 0.25, 0.5, 0.75, 1)), include.lowest=TRUE, labels=FALSE) })
P.e <- ncol(X.q)
profiles <- rbind(rep(-0.5,P.e), rep(0.5, P.e)) 

##Covariates ----
Age <- Table$age
Sex <- Table$sex
Education <- Table[c("education1","education2")]
Ethnicity <- Table[c("race1","race2")]
Energy <- Table$total.energy
Fiber <- Table$fiber
W <- data.frame(Age, Sex, Education, Ethnicity,Energy,Fiber)
Q <- ncol(W)

#Library Size ----
L <- Y
L <- L %>%
  mutate(LibrarySize=rowSums(across(everything())))
L <- L %>%
  select(LibrarySize)

#Z Matrices ----
## Genus ----
GenusData <- data.frame(t(FormattedObject$Species.Genus.Matrix))
speciesNames <- names(Table)[grep("d__", names(Table))]
GenusMatrix <- as.matrix(GenusData[match(speciesNames, row.names(GenusData)),])
GenusIndicator <- GenusMatrix%*%(1:ncol(GenusMatrix))
#Fix replication errors due to similar names
GenusMatrix[c(6, 21, 25, 54, 62, 82, 92, 97, 106, 107, 130, 144, 153, 154, 155, 157, 176, 181, 187, 199, 205), c(39, 20)] <- 0
GenusMatrix <- GenusMatrix[-149,] #Remove mitochondria from rownames
GenusMatrix <- GenusMatrix[,-79] #Remove mitochondria from colnames
P.g <- ncol(GenusMatrix)
#Check only 1 genus per species (species cant be in multiple genera)
numGenusPerSpecies <- as.numeric(apply(GenusMatrix, 1, sum)) # note that one species has two genus and two have no genus
#Check mitochondria removed
grep("Rickettsiales_f__Mitochondria",colnames(GenusMatrix))
grep("Rickettsiales_f__Mitochondria",rownames(GenusMatrix))
GenusData <- data.frame(GenusMatrix)


#Family
FamilyData <- data.frame(t(FormattedObject$Genus.Family.Matrix))
FamilyData <- FamilyData[-79,] #Remove mitochondria from rownames
FamilyData <- FamilyData[,-28] #Remove mitochondria from colnames
P.f <- ncol(FamilyData)
numFamilyPerGenus <- as.numeric(apply(FamilyData, 1, sum))
grep("Rickettsiales_f__Mitochondria",colnames(FamilyData))
grep("Rickettsiales_f__Mitochondria",rownames(FamilyData))

#Order
OrderData <- data.frame(t(FormattedObject$Family.Order.Matrix))
grep("Rickettsiales",colnames(OrderData))
grep("Rickettsiales",rownames(OrderData))
OrderData <- OrderData[-28,]
OrderData <- OrderData[,-21]
P.o <- ncol(OrderData)
numOrderPerorder <- as.numeric(apply(OrderData, 1, sum))

# Class
ClassData <- data.frame(t(FormattedObject$Order.Class.Matrix))
ClassData <- ClassData[-21,]
# ClassData[,4] <- rep(0, nrow(ClassData))
# ClassData[,10] <- rep(0, nrow(ClassData))
grep("Rickettsiales",colnames(ClassData))
grep("Rickettsiales",rownames(ClassData))
P.c <- ncol(ClassData)
numClassPerOrder <- as.numeric(apply(ClassData, 1, sum))

# Phylum
PhylumData <- data.frame(t(FormattedObject$Class.Phylum.Matrix))
P.p <- ncol(PhylumData)
numPhylumPerClass <- as.numeric(apply(PhylumData, 1, sum))

rm(GenusIndicator,GenusMatrix,Education,Ethnicity)

Z.s.g <- GenusData
Genus.R <- ncol(Z.s.g)
#Family
Z.g.f <- FamilyData
Family.R <- ncol(Z.g.f)
#Order
Z.f.o <- OrderData
Order.R <- ncol(Z.f.o)
#Class
Z.o.c <- ClassData
Class.R <- ncol(Z.o.c)
#Phylum
Z.c.p <- PhylumData
Phylum.R <- ncol(Z.c.p)


#3. Rename Simulation Variables ----
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
p <- (ggtree::ggtree(tree, layout = "fan", open.angle= 40)) +
  geom_tiplab(aes(label=label, color=ifelse(label %in% phy.labels, "highlight", "normal")))
# dev.off()


#Causal scenarios
scenario1 <- c()
scenario2 <- c(phy.labels[which(grepl("g__Alistipes",phy.labels))],phy.labels[which(grepl("g__Parabacteroides",phy.labels))],
               phy.labels[which(grepl("g__Odoribacter",phy.labels))])
scenario3 <- c(scenario2,phy.labels[which(grepl("g__Prevotellaceae_UCG_001",phy.labels))],phy.labels[which(grepl("g__Prevotellaceae_NK3B31_group",phy.labels))],
               phy.labels[which(grepl("g__Prevotella--s__unclassified946",phy.labels))],phy.labels[which(grepl("g__Paraprevotella",phy.labels))])
scenario4 <- c(scenario3,phy.labels[which(grepl("g__Prevotella_9",phy.labels))])
scenario5 <- c(scenario4,"d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified1",   
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified446", 
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified557", 
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified668", 
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified746")
scenario6 <- c(scenario5,"d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified779", 
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified1113",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified1335",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified2112",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified2223")
scenario7 <- c(scenario6,"d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified2555",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified2828",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified2884",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified3282",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified4101")
scenario8 <- c(scenario7,"d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified5541",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified5542",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified6541",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified7319",
               "d__Bacteria--p__Bacteroidota--c__Bacteroidia--o__Bacteroidales--f__Bacteroidaceae--g__Bacteroides--s__unclassified8874")

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




# P.e <- 7
P.e.causal <- c(0,3,6)
N <- 100
OR.exposure <- c(1,1.5)
P.s.scenario <- c(1:8)
Iterations <- 100
Iter <- 1:100
sim.par <- expand.grid(P.e.causal = P.e.causal, P.s.scenario = P.s.scenario, N = N, OR.exposure = OR.exposure) #36 combinations
sim.par <- sim.par %>%
  filter((OR.exposure==1 & P.s.scenario==1 & P.e.causal==0) | (OR.exposure==1.5 & P.e.causal!=0) & OR.exposure==1.5 & P.s.scenario!=1) #20
# saveRDS(sim.par,"/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Simulation_Parameters_03_21.rds")
# saveRDS(sim.par,fs::path(here::here(),"Formatted Data","Simulation_Parameters_all_04_14.rds"))

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

