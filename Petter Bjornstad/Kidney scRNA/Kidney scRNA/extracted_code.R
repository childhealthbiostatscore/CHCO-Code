## ----include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------
library(reprex)
library(tidyverse)
library(BiocManager)        
library(arsenal)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(future)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(venn)
library(rstatix)
library(table1)
library(Biobase)
library(ReactomeGSA)
library(GSEABase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(SingleCellExperiment)
library(fgsea)
library(EnhancedVolcano)
library(openxlsx)
library(BiocManager)
# library(MAST)
library(ggrepel)
# library(qpcR)
library(ggpubr)
library(openxlsx)
library(ggplot2)
library(GGally)
library(GSEABase)
library(limma)
library(reshape2)
library(data.table)
library(knitr)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
#library(NMF)
library(rsvd)
library(RColorBrewer)
# library(MAST)
library(devtools)
# install_github("Sun-lab/ideas",force=T)
#library(ideas)
library(foreach)
library(parallel)
library(doRNG)
library(doParallel)
library(fs)
# registerDoParallel(cores = 6)
library(VennDiagram)
library(janitor)
# devtools::install_github('immunogenomics/presto')
# library(presto)
library(knitr)
library(lme4)
library(lmerTest)
#install.packages("glmmTMB")
# Reinstall glmmTMB from source
# install.packages("glmmTMB", type = "source")
library(glmmTMB)
# Install DoubletFinder (if not already installed)
# devtools::install_github("chris-mcginnis-ucsf/DoubletFinder",force=T)
# Load the package
# Install DoubletFinder from GitHub (use devtools to install)
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("chris-mcginnis-ucsf/DoubletFinder",force=T)
# library(DoubletFinder)
# install.packages("emmeans")
library(emmeans)
library(pheatmap)
library(enrichplot)
library(enrichR)
dbs <- c("GO_Biological_Process_2023", 
         "KEGG_2021_Human",
         # "Reactome_2022", 
         "Reactome_Pathways_2024",
         # "MSigDB_Oncogenic_Signatures",
         # "MSigDB_Computational",
         "MSigDB_Hallmark_2020")
dbs_celltype <- c(
  "GO_Biological_Process_2023",
  "KEGG_2021_Human",
  "Reactome_Pathways_2024",
  "MSigDB_Hallmark_2020",
  # Cell type specific databases
  "CellMarker_2024",  # Cell type marker genes
  "Azimuth_Cell_Types_2021",  # Cell type signatures
  "PanglaoDB_Augmented_2021",  # Single-cell markers
  "Descartes_Cell_Types_and_Tissue_2021",  # Developmental cell types
  "Human_Gene_Atlas",  # Tissue/cell expression
  "ARCHS4_Tissues"  # Tissue expression
)
# BiocManager::install("edgeR",force=T)
library(edgeR)
library(devtools)
# install_github("lhe17/nebula")
# library(nebula)

# remove.packages("boot")  # Remove broken version
# install.packages("boot", type = "source")  # Reinstall from source
library(boot)
library(furrr)
library(future)
# BiocManager::install("scran")
library(scran)
library(BiocParallel)
# library(DESeq2)
if (!require("GSA", quietly = TRUE)) {
  install.packages("GSA")
}
library(GSA)
# if (!require("reticulate", quietly = TRUE)) {
#   install.packages("reticulate")
# }

# #Mac laptop
# dir.dat <- c("//Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive")
# dir.dat2 <- c("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/scRNA")
# dir.code <- c("/Users/hhampson/Documents/CHCO-Code/Petter Bjornstad/Liver analysis/Liver scRNAseq")
# dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Kidney scRNAseq Project/Organoid Results/NEBULA")

# #Mac Studio File Path
# dir.dat <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive")
# dir.dat2 <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/scRNA")
# dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Kidney scRNAseq Project/Organoid Results/NEBULA")

##c. Load functions ----
source("Kidney_functions_sc.R")
source("Rename_Seurat_Function.R")
source("Celltype_Annotation_Function.R")


## ----Cyberduck setup,include=F,echo=F------------------------------------------------------------------------------------------------------------------------------------
# install.packages("reticulate",force=T)
library(reticulate)
reticulate::use_python("/home/hhampson/miniconda3/bin/python") # replace with your username

## Load boto3 and pandas
boto3 <- reticulate::import("boto3")
pd <- reticulate::import("pandas")

## Create an S3 client
# install.packages("jsonlite")  # Install if not already installed
library(jsonlite)  # Load the package

keys <- fromJSON("/home/hhampson/keys.json") # replace with your Lambda username
session <- boto3$session$Session(
  aws_access_key_id = keys$MY_ACCESS_KEY,
  aws_secret_access_key = keys$MY_SECRET_KEY
)

## Create an S3 client with the session
s3 <- session$client("s3", endpoint_url = "https://s3.kopah.uw.edu")

gc()
bucket <- "scrna" # bucket name in Kopah
temp_file <- tempfile(fileext = ".rds") # need to create a temporary file
s3$download_file(bucket, "Kidney transcriptomics/Single cell RNA seq/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds", temp_file)
so_kpmp_sc <- readRDS(temp_file)
# # read file
# bucket <- "scrna" # bucket name in Kopah
# temp_file <- tempfile(fileext = ".rds") # need to create a temporary file
# s3$download_file(bucket,"Liver transcriptomics/Single nucleus RNA seq/NoRef_PetterLiver_ClinData_Labels_Char_041924.rds", temp_file)
# so_liver_sn <- readRDS(temp_file)

# DefaultAssay(so_liver_sn) <- "RNA"
# dim(so_liver_sn)#36601 130124
# nrow(so_liver_sn) #36,601 genes
# ncol(so_liver_sn) # 130,124 nuceli
# invisible(gc())

# #Load metadata
# gc()
# bucket <- "scrna" # bucket name in Kopah
# temp_file <- tempfile(fileext = ".csv") # need to create a temporary file
# s3$download_file(bucket, "Liver transcriptomics/liver_biopsy_metadata_PN.csv", temp_file)
# meta_liver_raw <- read.csv(temp_file)
# gc()
gc()
bucket <- "harmonized.dataset" # bucket name in Kopah
temp_file <- tempfile(fileext = ".rds") # need to create a temporary file
s3$download_file(bucket, "harmonized_dataset.csv", temp_file)
harm_meta_data <- read.csv(temp_file)


# Filter out the gmt files for KEGG, Reactome and GOBP
# list.files(bg_path)
# gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
# gmt_files
# kegg_legacy <- prepare_gmt(gmt_files[1], unique(full_results$Gene), savefile = FALSE)
# reactome <- prepare_gmt(gmt_files[3], unique(full_results$Gene), savefile = FALSE)
# go <- prepare_gmt(gmt_files[4], unique(full_results$Gene), savefile = FALSE)
# 
# 
# temp_file <- tempfile(fileext = ".gmt")
# s3$download_file("gsea", "gmt_files/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt", temp_file)
# kegg <- GSA.read.gmt(temp_file)

temp_file <- tempfile(fileext = ".gmt")
s3$download_file("gsea", "gmt_files/c8.all.v2025.1.Hs.symbols.gmt", temp_file)
c8 <- GSA.read.gmt(temp_file)

#Load formatted organoid data
gc()
bucket <- "scrna" # bucket name in Kopah
temp_file <- tempfile(fileext = ".rds") # need to create a temporary file
s3$download_file(bucket, "Kidney organoids/organoids_processed_matrix.rds", temp_file)
sc_org <- readRDS(temp_file)


#Load kidney cell type annotations
gc()
bucket <- "scrna" # bucket name in Kopah
temp_file <- tempfile(fileext = ".csv") # need to create a temporary file
s3$download_file(bucket, "Kidney organoids/top_20_markers_15.csv", temp_file)
annot_15 <- read.csv(temp_file)
annot_15 <- annot_15 %>% 
  dplyr::select(c("cluster","Celltype")) %>% 
  distinct()
  
gc()
bucket <- "scrna" # bucket name in Kopah
temp_file <- tempfile(fileext = ".csv") # need to create a temporary file
s3$download_file(bucket, "Kidney organoids/top_20_markers_30.csv", temp_file)
annot_30 <- read.csv(temp_file)
annot_30 <- annot_30 %>% 
  dplyr::select(c("cluster","Celltype")) %>% 
  distinct()



## ----include=T,echo=T----------------------------------------------------------------------------------------------------------------------------------------------------
#Human Biopsy samples in PB90 origional number of cells
cells_all <-  ncol(so_kpmp_sc) #211218 cells before processing
genes_all <- nrow(so_kpmp_sc) # 31332 genes before processing
print(paste0(cells_all," cells and ",genes_all," genes before processing in PB90"))

#Set ids for organoid samples
ids <- c("CRC-10","CRC-11","CRC-03","RH-50-T","RH-72-T","RH-62-T","IT_19")

#Filter to organoid samples only
so_kpmp_sc <- subset(so_kpmp_sc,record_id %in% ids)
cells_org <-  ncol(so_kpmp_sc) #14348 cells before processing
genes_org <- nrow(so_kpmp_sc) #31332 genes before processing
print(paste0(cells_org," cells and ",genes_org," genes after filtering to 7 participants with organoid data"))



