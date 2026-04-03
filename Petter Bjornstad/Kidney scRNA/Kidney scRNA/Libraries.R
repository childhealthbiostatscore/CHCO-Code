library(reticulate)
use_python("/home/hhampson/miniconda3/bin/python", required = TRUE)
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
# BiocManager::install("speckle")
library(speckle)
library(limma)
library(edgeR)  # often used together for normalization

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
library(nebula)