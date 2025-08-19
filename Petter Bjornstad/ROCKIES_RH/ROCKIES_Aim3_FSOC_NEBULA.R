library(scran)
library(future)
library(future.apply)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(REDCapR)
library(data.table)
library(emmeans)
library(NMF)
library(pheatmap)
library(UpSetR)
library(enrichR)
library(WriteXLS)
library(SAVER)
library(readxl)
library(limma)
library(edgeR)
library(BiocGenerics)
library(GSEABase)
library(slingshot)
library(SingleCellExperiment)
library(MAST)
library(muscat)
library(scater)
library(Seurat)
library(jsonlite)
library(dplyr)
library(glmmTMB)
library(reshape2)
library(broom.mixed)
library(nebula)



#load in correct data, ensure FSOC accuracy. Only Type 2's 


load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')

dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/'

load("C:/Users/netio/Downloads/TCA_genes.txt")
load('C:/Users/netio/Downloads/OxPhos_genes.txt')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

test <- so_subset@meta.data

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat <- dat %>% dplyr::select(record_id, epic_sglti2_1, starts_with('fsoc')) %>% 
  filter(!is.na(epic_sglti2_1))

#Fix up FSOC data
test$epic_sglti2_1 <- NULL
test$fsoc_l_cortex <- NULL
test$fsoc_l_kidney <- NULL
test$fsoc_l_medulla <- NULL
test$fsoc_r_cortex <- NULL
test$fsoc_r_kidney <- NULL
test$fsoc_r_medulla <- NULL

find_fsoc_averages <- function(data){
  tmp_data <- data %>% dplyr::select(starts_with('fsoc'))
  fsoc_full_combined <- rowMeans(tmp_data, na.rm=T)
  tmp_data <- data %>% dplyr::select(starts_with('fsoc_l_'))
  fsoc_l_combined <- tmp_data %>% rowMeans(na.rm=T)
  tmp_data <- data %>% dplyr::select(starts_with('fsoc_r_'))
  fsoc_r_combined <- tmp_data %>% rowMeans(na.rm=T)
  
  tmp_df <- cbind(fsoc_l_combined, fsoc_r_combined, fsoc_full_combined)
  return(tmp_df)
  
}

tmp_results <- find_fsoc_averages(dat)




test <- test %>% left_join(dat, by='record_id')

so_subset@meta.data$epic_sglti2_1 <- test$epic_sglti2_1
so_subset@meta.data$fsoc_l_cortex <- test$fsoc_l_cortex
so_subset@meta.data$fsoc_l_kidney <- test$fsoc_l_kidney
so_subset@meta.data$fsoc_l_medulla <- test$fsoc_l_medulla
so_subset@meta.data$fsoc_r_cortex <- test$fsoc_r_cortex
so_subset@meta.data$fsoc_r_kidney <- test$fsoc_r_kidney
so_subset@meta.data$fsoc_r_medulla <- test$fsoc_r_medulla
so_subset@meta.data$fsoc_l_combined <- test$fsoc_l_combined
so_subset@meta.data$fsoc_r_combined <- test$fsoc_r_combined
so_subset@meta.data$fsoc_full_combined <- test$fsoc_full_combined




rm(test)



so_subset$celltype1 <- case_when(grepl("PT-",so_subset$celltype_rpca)~"PT",
                                 grepl("TAL-",so_subset$celltype_rpca)~"TAL",
                                 grepl("EC-",so_subset$celltype_rpca)~"EC",
                                 grepl("POD",so_subset$celltype_rpca)~"POD",
                                 grepl("MAC",so_subset$celltype_rpca)~"MAC",
                                 grepl("MON",so_subset$celltype_rpca)~"MON",
                                 grepl("PC-",so_subset$celltype_rpca)~"PC",
                                 grepl("FIB",so_subset$celltype_rpca)~"FIB_MC_VSMC",
                                 grepl("DTL",so_subset$celltype_rpca)~"DTL",
                                 so_subset$celltype_rpca=="DCT"~"DCT",
                                 so_subset$celltype_rpca=="ATL"~"ATL",
                                 so_subset$celltype_rpca=="B"~"B",
                                 so_subset$celltype_rpca=="T"~"T")
so_subset$celltype1 <- as.character(so_subset$celltype1)

so_subset$KPMP_celltype2 <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2 <- ifelse(so_subset$KPMP_celltype=="aPT" | 
                                so_subset$KPMP_celltype=="PT-S1/S2" | 
                                so_subset$KPMP_celltype == "PT-S3","PT",
                              ifelse(grepl("TAL",so_subset$KPMP_celltype),"TAL",
                                     ifelse(grepl("EC-",so_subset$KPMP_celltype),"EC",so_subset$KPMP_celltype2)))


so_subset$DCT_celltype <- ifelse((so_subset$KPMP_celltype=="DCT" | 
                                    so_subset$KPMP_celltype=="dDCT"), "DCT","Non-DCT")
















#filter data to only TAL cells 


so_subset <- subset(so_subset, )






#gene lists
tal_genes 
tal_genes <- c("SLC12A1", "KCNJ1", "CLCNKB", "BSND", "CLDN16", 
               "ATP1A1", "ATP1B1", "ATP1B3", "FXYD2", 
              "PPARGC1A", "NRF1", "TFAM", "COX4I1", "COX6A1", "COX7A1", "ATP5F1A", "ATP5F1B",
    "UQCRB", "UQCRC1","NDUFA4", "NDUFB5", "NDUFS1", 
    "UMOD", "CLDN10", "WNK1", "WNK4","STK39","OXSR1")








counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round



design_fsoc <- model.matrix(
  ~ FSOC_medulla 
  design_fsoc 
  FSOC_medulla + FSOC_cortex 
  count_mat))
FSOC_cortex + group + age + sex + eGFR + BMI,
data = metadata
)


































