library(arsenal)
library(Biobase)
library(BiocGenerics)
library(BiocParallel)
library(broom.mixed)
library(colorspace)
library(cowplot)
library(data.table)
library(DirichletReg)
library(dplyr)
library(edgeR)
library(emmeans)
library(enrichR)
library(foreach)
library(future)
library(future.apply)
library(GSEABase)
library(ggdendro)
library(ggpubr)
library(glmmTMB)
library(harmony)
library(jsonlite)
library(kableExtra)
library(limma)
library(MAST)
library(Matrix)
library(msigdbr)
library(muscat)
library(NMF)
library(nebula)
library(patchwork)
library(pheatmap)
library(readxl)
library(REDCapR)
library(reshape2)
library(rstatix)
library(SAVER)
# library(scDC)
library(scater)
library(scran)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(UpSetR)
library(WriteXLS)
library(parallel)
library(doParallel)
library(quantreg)
library(aws.s3)

# Set up environment for Kopah
keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# Pull in PB90 Seurat object
pb90 <- s3readRDS(object = "Kidney transcriptomics/Single cell RNA seq/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds", 
                  bucket = "scrna",  region = "")

# Pull in clinical dataset
rh_rh2_croc_improve_unique <- 
  s3readRDS(object = "data_clean/rh_rh2_croc_improve_unique.RDS", 
            bucket = "harmonized.dataset",
            region = "")

# subset to T2D and OB from RH/RH2/IMPROVE Baseline, HC from CROC
pb90_subset <- subset(pb90, group != "Type_1_Diabetes" & T2D_HC_Phil != "HC_igA" & visit == "baseline")
length(unique(pb90_subset$record_id)) # 49 people with biopsies

# add metadata to PB90 subset
meta_subset <- left_join(subset(pb90_subset@meta.data, select = -c(age, eGFR_CKD_epi, hba1c)), 
                         rh_rh2_croc_improve_unique,
                         by = "kit_id") 
length(unique(meta_subset$record_id.x)) # matches 49 people with biopsies
rownames(meta_subset) <- meta_subset$barcode
pb90_subset@meta.data <- meta_subset

# add general cell type labels
celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),   # grouped thin limbs
  DCT_CNT = c("DCT", "dDCT", "CNT"),   # grouped distal tubule/connecting tubule
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A"), 
  Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",                         # mesangial cells
  PEC = "PEC",                       # parietal epithelial cells
  Schwann = "SchwannCells",
  Other = c("non-specific")          # catchall
)

map_celltype_to_general <- function(celltype, celltype_groups) {
  for (group_name in names(celltype_groups)) {
    if (celltype %in% celltype_groups[[group_name]]) {
      return(group_name)
    }
  }
  return("Other")  # For any celltype not found in groups
}


pb90_subset$KPMP_celltype_general <- sapply(pb90_subset$KPMP_celltype, 
                                            map_celltype_to_general, 
                                            celltype_groups = celltype_groups)


# filter out rare gene
ncol(pb90_subset) # Number of cells before filtering (117859)
nrow(pb90_subset) # Number of genes before filtering (31332)

# Filter out rare genes expressed in less than "gene_pct" of cells
expr_matrix <- as.matrix(GetAssayData(pb90_subset, assay = "RNA", layer = "counts"))

# Get cell type information (adjust column name if needed)
cell_types <- pb90_subset@meta.data$KPMP_celltype_general

# Calculate proportion of cells expressing each gene within each cell type
gene_pct_threshold <- 0.05  # 5% threshold

# Initialize a matrix to store proportions for each gene in each cell type
unique_cell_types <- unique(cell_types)
gene_proportions_by_type <- matrix(0, 
                                   nrow = nrow(expr_matrix), 
                                   ncol = length(unique_cell_types))
rownames(gene_proportions_by_type) <- rownames(expr_matrix)
colnames(gene_proportions_by_type) <- unique_cell_types

# Calculate proportion for each cell type
for(ct in unique_cell_types) {
  # Get cells of this type
  cells_of_type <- which(cell_types == ct)
  
  # Calculate proportion of cells expressing each gene in this cell type
  expr_subset <- expr_matrix[, cells_of_type]
  num_cells_expressing <- rowSums(expr_subset > 0)
  total_cells_in_type <- length(cells_of_type)
  
  gene_proportions_by_type[, ct] <- num_cells_expressing / total_cells_in_type
}

# Keep genes expressed in at least 5% of cells in AT LEAST ONE cell type
max_proportion_per_gene <- apply(gene_proportions_by_type, 1, max)
genes_to_keep <- names(max_proportion_per_gene[max_proportion_per_gene >= gene_pct_threshold])

# Check results
cat("Total genes:", nrow(expr_matrix), "\n") # Total genes: 31332 
cat("Genes to keep:", length(genes_to_keep), "\n") # Genes to keep: 16286 
cat("Proportion kept:", round(length(genes_to_keep)/nrow(expr_matrix)*100, 1), "%\n") # Proportion kept: 52 %

# Optional: Check distribution of max proportions
hist(max_proportion_per_gene)

# Filter the Seurat object
pb90_subset <- subset(pb90_subset, features = genes_to_keep)


#Check the number of Mitochondrial genes to start
sum(grepl("^MT-", rownames(pb90_subset))) 

# Identify mitochondrial genes (human: start with "MT-")
mito_genes <- grep("^MT-", rownames(pb90_subset), value = TRUE)
pb90_subset <- subset(pb90_subset, features = setdiff(rownames(pb90_subset), mito_genes))

#Check the number of Mitochondrial genes after filtering to ensure filtering step was successful
sum(grepl("^MT-", rownames(pb90_subset))) #Should be 0

# Identify ribosomal genes
ribo_genes <- c(
  "RPL22", "RPL11", "RPS8", "RPL5", "RPS27", "RPS7", "RPS27A", "RPL31", "RPL37A", "RPL32", "RPL15", "RPL14", "RPL29",
  "RPL24", "RPL22L1", "RPL35A", "RPL9", "RPL34", "RPS3A", "RPL37", "RPS23", "RPS14", "RPS18", "RPS10", "RPL10A", 
  "RPS20", "RPL7", "RPL30", "RPL8", "RPS6", "RPL35", "RPL12", "RPL7A", "RPS24", "RPLP2", "RPL27A", "RPS13", "RPS3",
  "RPS25", "RPS26", "RPL41", "RPL6", "RPLP0", "RPL21", "RPS29", "RPL4", "RPLP1", "RPS17", "RPS2", "RPS15A", "RPL13",
  "RPL26", "RPL23A", "RPL23", "RPL19", "RPL27", "RPL38", "RPL17", "RPS15", "RPL36", "RPS28", "RPL18A", "RPS16", 
  "RPS19", "RPL18", "RPL13A", "RPS11", "RPS9", "RPL28", "RPS5", "RPS21", "RPL3", "RPS4X", "RPL36A", "RPL39", 
  "RPL10", "RPS4Y1"
) # grep("^RPL|^RPS", rownames(attempt_so_raw), value = TRUE) captures some none ribosomal genes

pb90_subset <- subset(pb90_subset, features = setdiff(rownames(pb90_subset), ribo_genes))


# Check if data are normalized
Assays(pb90_subset)  # Should show e.g., "RNA" with "counts" and "data"
pb90_subset@assays # Another method to see RNA assays
head(GetAssayData(pb90_subset, layer = "counts")[, 1:5])  # Raw counts
head(GetAssayData(pb90_subset, layer = "data")[, 1:5])    # Normalized data


# Exploration of data
# Randomly select 100 genes from the Seurat object
genes <- sample(rownames(pb90_subset), 100)

# Set up a 2x3 plotting layout so you can plot multiple histograms in one figure
par(mfrow = c(2, 3))

# Loop over each randomly selected gene
for (g in genes) {
  
  # Plot a histogram of the expression values for gene 'g'
  # using normalized expression values from the "data" slot
  hist(GetAssayData(pb90_subset, layer = "data")[g, ],
       main = g,                # Title of the plot = gene name
       xlab = "Normalized Expression")     # Label for x-axis
  
  # using raw counts expression values from the "counts" slot
  hist(GetAssayData(pb90_subset, layer = "counts")[g, ],
       main = g,                # Title of the plot = gene name
       xlab = "Raw Counts Expression")     # Label for x-axis
}

# Save Seurat object for downstream analysis
s3saveRDS(pb90_subset, object = "data_clean/subset/pb90_ckd_analysis_subset.rds", bucket = "scrna", region = "")
