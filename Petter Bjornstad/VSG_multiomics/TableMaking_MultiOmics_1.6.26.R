##### Making tables for multiomics VSG paper 


library(readr)
library(dplyr)
library(gridExtra)
library(Seurat)
library(lme4)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(stats) 
#
library(lmerTest)
library(broom.mixed)
library(emmeans)
#
library(UpSetR)
library(tidyr)
#
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggridges)
library(tidytext)
library(msigdbr)
library(fgsea)
library(stringr)
library(purrr)
library(rlang)
#
library(KEGGREST)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)
library(edgeR)
library(limma)
library(ggrepel)
#

base_dir <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/scrna/ "

scrna <- readRDS(paste0(base_dir, "data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds"))
scrna <- subset(scrna, subset = percent.mt < 50 & nFeature_RNA < 5000 & nFeature_RNA > 500) 

#
clin <- read_csv(paste0(base_dir, "data_clean/pb90_rhrh2improve_clinical_subset.csv", show_col_types = FALSE)) %>%
  dplyr::select(record_id, visit,acr_u, gfr_bsa_plasma, gfr_raw_plasma, bmi, eGFR_CKiD_U25_avg, eGFR_fas_cr_cysc) %>%
  dplyr::mutate(visit = dplyr::recode(visit,"baseline" = "pre","12_months_post_surgery" = "post"))


#############
## IMPROVE ##
#############
improve <- subset(
  scrna, 
  subset = cohort == "IMPROVE" | 
    record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T"))
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$pre_post <- 'pre'
improve@meta.data <- improve@meta.data %>%
  tibble::rownames_to_column("cell") %>%   
  left_join(clin, by = c("record_id" = "record_id", "pre_post" = "visit")) %>%
  tibble::column_to_rownames("cell") 
improve@meta.data[improve@meta.data$record_id == "RH-59-T", ]$record_id <- 'IT_07'
improve@meta.data[improve@meta.data$record_id == "RH-60-T", ]$record_id <- 'IT_08'
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$record_id <- 'IT_09'
improve@meta.data[improve@meta.data$record_id == "RH-66-T", ]$record_id <- 'IT_10'

########
## RH ##
########
rh <- subset(
  scrna,
  subset = cohort %in% c('RENAL HEIR', 'RENAL HEIRITAGE') &
    !scrna$record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T", "RH2-07-O")
)
rh$pre_post <- ifelse(rh$cohort == "RENAL HEIR", "pre", "post")
rh@meta.data <- rh@meta.data %>%
  tibble::rownames_to_column("cell") %>%   
  left_join(clin, by = c("record_id" = "record_id")) %>%
  tibble::column_to_rownames("cell") 
#
rh@meta.data[rh@meta.data$record_id == "RH2-14-T", ]$record_id <- 'RH-23-T'
rh@meta.data[rh@meta.data$record_id == "RH2-19-T", ]$record_id <- 'RH-67-T'

#
improve@meta.data$treatment <- "VSG"
rh@meta.data$treatment <- "Standard"

healthy <- subset(scrna, 
                  subset = record_id %in% c("CRC-03","CRC-14","CRC-02","CRC-13",
                                            "CRC-10","CRC-39","CRC-40","CRC-11",
                                            "CRC-46","CRC-58","CRC-56","CRC-54"))
healthy@meta.data$treatment <- "Healthy"
healthy@meta.data$pre_post <- "healthy"
#
combined <- merge(improve, rh)






































