#Test code for Laura 
#Load Packages
library(lme4)
library(lmerTest)
library(emmeans)
library(Seurat)
library(tidyr)
library(dplyr)

#Load in subsetted data
so <- readRDS('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/Test_Seurat_Object_for_Laura.rds')

#Run Mixed Effect Model - 10 test genes in aPT cells only
# Extract the gene expression data for all genes
gene_expression <- as.data.frame(GetAssayData(so, layer = "data"))

#Transpose gene expression dataset to merge with metadata
gene_expression <- t(gene_expression)
gene_expression <- data.frame(gene_expression) #Make a dataframe again after transposing

#Set gene list
gene_list <- colnames(gene_expression)
gene_expression$cellname <- rownames(gene_expression) 
rownames(gene_expression) <- NULL

# Extract the metadata
metadata <- so@meta.data
# metadata <- metadata %>%
#   mutate(across(everything(),~ifelse(.==".",NA,.)))
metadata$cellname <- rownames(metadata)
rownames(metadata) <- NULL

# Combine the gene expression data and metadata
data <- tidylog::full_join(metadata,gene_expression,by="cellname")
rm(metadata,gene_expression)

#Subset to T2D only to compare medication groups in T2D
data_subset <- data %>%
  filter(group=="Type_2_Diabetes") 

#Make sure med group is faster and that no_med is the reference group
data_subset$medication <- factor(data_subset$medication)

#Prepare dataframe for results
full_results <- data.frame()

# Total number of genes
total_genes <- length(gene_list)

# Calculate the batch size
batch_size <- round(total_genes / 5)

# Simulate a vector of genes (replace this with your actual gene data)
genes <- 1:total_genes

# Loop through the genes in batches
for (i in seq(1, total_genes, by = batch_size)) {
  
  # Get the current batch (subsetting the gene vector)
  batch <- genes[i:min(i + batch_size - 1, total_genes)]
  for (gene in gene_list[batch]) { #tester genes
    # m0 <- as.formula(paste0(gene," ~ group + (1 | kit_id)"))
    # model <- lmer(m0,data=data_subset)
    #Adjust for key covariates
    m1 <- as.formula(paste0(gene," ~ medication + age + sex + bmi + (1 | kit_id)"))
    model1 <- lmer(m1,data=data_subset)
    
    emm_options(pbkrtest.limit = 4241)
    # Compute estimated marginal means
    emm <- emmeans(model1, ~ medication)
    
    # Define and test the contrasts
    contrasts <- contrast(emm, list(
      "GLP1 vs No Med" = c(0, 1, 0, -1),  # GLP1 vs No Med
      "SGT2 vs No Med" = c(0, 0, 1, -1),  # SGT2 vs No Med
      "GLP1+SGT2 vs No Med" = c(0, 0, 0, 1),  # GLP1+SGT2 vs No Med
      "GLP1+SGT2 vs GLP1" = c(0, -1, 1, 0),  # GLP1+SGT2 vs GLP1
      "GLP1+SGT2 vs SGT2" = c(0, 0, -1, 1)  # GLP1+SGT2 vs SGT2
    ))
    print(contrasts)

    contrasts_new <- contrast(emm, list(
      "GLP1 vs No Med" = c(-1, 0, 1, 0),  # GLP1 vs No Med
      "SGT2 vs No Med" = c(-1, 1, 0, 0),  # SGT2 vs No Med
      "GLP1+SGT2 vs No Med" = c(-1, 0, 0, 1),  # GLP1+SGT2 vs No Med
      "GLP1+SGT2 vs GLP1" = c(0, 0, -1, 1),  # GLP1+SGT2 vs GLP1
      "GLP1+SGT2 vs SGT2" = c(0, -1, 0, 1)  # GLP1+SGT2 vs SGT2
    ))
    print(contrasts_new)
    
  }
}


