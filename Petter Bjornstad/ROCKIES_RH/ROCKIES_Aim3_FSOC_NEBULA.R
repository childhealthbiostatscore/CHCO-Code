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
library(doParallel)



#load in correct data, ensure FSOC accuracy. Only Type 2's 


load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')


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
  
  fsoc_medulla <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_medulla)
  fsoc_cortex <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_cortex)
  fsoc_kidney <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_kidney)
  
  tmp_df <- cbind(fsoc_l_combined, fsoc_r_combined, fsoc_medulla, fsoc_cortex, fsoc_kidney, fsoc_full_combined)
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
so_subset@meta.data$fsoc_medulla <- test$fsoc_medulla
so_subset@meta.data$fsoc_cortex <- test$fsoc_cortex
so_subset@meta.data$fsoc_kidney <- test$fsoc_kidney





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


so_subset <- subset(so_subset, celltype2 == 'TAL')


##Get all FSOC data 

harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')

#harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat2 <- dat %>% filter(visit == 'baseline') %>% 
  filter(study %in% c('RENAL-HEIR', 'RENAL-HEIRitage', 'CROCODILE') | record_id == 'IT_19') %>%
  #  filter(group != 'Obese Control') %>% 
  dplyr::select(mrn, record_id, study, visit, group, starts_with('fsoc'), bmi, 
                epic_sglti2_1, age, sex, epic_mfm_1, epic_insulin_1, epic_glp1ra_1)

dat2 <- dat2[-which(dat2$study == 'CROCODILE' & dat2$group == 'Type 1 Diabetes'),]


tests <- c('fsoc_l_cortex', 'fsoc_r_cortex', 
           'fsoc_l_kidney', 'fsoc_r_kidney', 
           'fsoc_l_medulla', 'fsoc_r_medulla', 
           'fsoc_l_combined', 'fsoc_r_combined',
           'fsoc_medulla', 'fsoc_cortex', 'fsoc_kidney',
           'fsoc_full_combined')




tmp_df <- dat2 %>% dplyr::select(starts_with('fsoc'))



find_fsoc_averages <- function(data){
  tmp_data <- data %>% dplyr::select(starts_with('fsoc'))
  fsoc_full_combined <- rowMeans(tmp_data, na.rm=T)
  
  tmp_data <- data %>% dplyr::select(starts_with('fsoc_l_'))
  fsoc_l_combined <- tmp_data %>% rowMeans(na.rm=T)
  
  tmp_data <- data %>% dplyr::select(starts_with('fsoc_r_'))
  fsoc_r_combined <- tmp_data %>% rowMeans(na.rm=T)
  
  fsoc_medulla <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_medulla) %>%
    rowMeans(na.rm=T)
  fsoc_cortex <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_cortex) %>%
    rowMeans(na.rm=T)
  fsoc_kidney <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_kidney) %>%
    rowMeans(na.rm=T)
  
  tmp_df <- cbind(fsoc_l_combined, fsoc_r_combined, fsoc_medulla, fsoc_cortex, fsoc_kidney, fsoc_full_combined)
  return(tmp_df)
  
}

tmp_results <- find_fsoc_averages(dat2)


dat_results <- dat2 %>% bind_cols(tmp_results)


meta.data <- so_subset@meta.data
dat_results <- dat_results %>% filter(record_id %in% meta.data$record_id)


table1::table1(~age + sex + epic_mfm_1 + epic_insulin_1 + epic_glp1ra_1 +
                 fsoc_l_cortex + fsoc_l_kidney + fsoc_l_medulla + 
                 fsoc_r_cortex + fsoc_r_kidney + fsoc_r_medulla + 
                 fsoc_l_combined + fsoc_r_combined + fsoc_medulla + 
                 fsoc_cortex + fsoc_kidney + fsoc_full_combined | group, data = dat_results)



names(meta.data)[which(names(meta.data) %in% c('fsoc_l_kidney', 'fsoc_r_kidney', 
                              'fsoc_l_medulla', 'fsoc_r_medulla', 'fsoc_l_cortex', 'fsoc_r_cortex',
                              'fsoc_l_combined', 'fsoc_r_combined',
                              'fsoc_medulla', 'fsoc_cortex', 'fsoc_kidney',
                              'fsoc_full_combined'))]

meta.data <- meta.data %>% dplyr::select(-fsoc_l_kidney, -fsoc_l_medulla, 
                                         -fsoc_r_kidney, -fsoc_r_medulla, 
                                         -fsoc_l_cortex, -fsoc_r_cortex) %>% 
  left_join(dat_results, by='record_id')



so_subset@meta.data$fsoc_l_kidney <- meta.data$fsoc_l_kidney
so_subset@meta.data$fsoc_l_medulla <- meta.data$fsoc_l_medulla
so_subset@meta.data$fsoc_l_cortex <- meta.data$fsoc_l_cortex

so_subset@meta.data$fsoc_r_kidney <- meta.data$fsoc_r_kidney
so_subset@meta.data$fsoc_r_medulla <- meta.data$fsoc_r_medulla
so_subset@meta.data$fsoc_r_cortex <- meta.data$fsoc_r_cortex

so_subset@meta.data$fsoc_l_combined <- meta.data$fsoc_l_combined
so_subset@meta.data$fsoc_r_combined <- meta.data$fsoc_r_combined
so_subset@meta.data$fsoc_medulla <- meta.data$fsoc_medulla
so_subset@meta.data$fsoc_cortex <- meta.data$fsoc_cortex
so_subset@meta.data$fsoc_kidney <- meta.data$fsoc_kidney
so_subset@meta.data$fsoc_full_combined <- meta.data$fsoc_full_combined








#filter data to only TAL cells 





#gene lists

tal_genes <- c("SLC12A1", "KCNJ1", "CLCNKB", "BSND", "CLDN16", 
               "ATP1A1", "ATP1B1", "ATP1B3", "FXYD2", 
              "PPARGC1A", "NRF1", "TFAM", "COX4I1", "COX6A1", "COX7A1", "ATP5F1A", "ATP5F1B",
    "UQCRB", "UQCRC1","NDUFA4", "NDUFB5", "NDUFS1", 
    "UMOD", "CLDN10", "WNK1", "WNK4","STK39","OXSR1")




dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/FSOC_NEBULA/'




counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round

sce <- SingleCellExperiment(assays = list(counts = counts_path))
# sce <- computeSumFactors(sce)
sce <- computeSumFactors(sce)
# View size factors
# sizeFactors(sce)
# STEP 3: Calculate offset → log(size factors)
pooled_offset <- sizeFactors(sce)
so_subset@meta.data$pooled_offset <- pooled_offset



genes_list <- tal_genes

so_subset <- subset(so_subset, features = tal_genes)

counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round
count_gene <- counts_path
meta_gene <- subset(so_subset)@meta.data


complete_idx <- complete.cases(meta_gene$fsoc_medulla)
cat("Cells with complete data:", sum(complete_idx), "\n")

# Step 3: Filter all your data to only include cells with complete predictor data
meta_gene <- meta_gene[complete_idx, ]
count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells


# Step 4: Create prediction matrix from the complete data
pred_gene <- model.matrix(~fsoc_medulla, data = meta_gene)


    # library <- meta_gene$library_size
library <- meta_gene$pooled_offset
    data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
    
    if (is.null(data_g_gene)) {
      data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
    }
    
    #With offset
    result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                     pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
    


#Make dataframe of final results
full_results <- as.data.frame(result)
#Calculate number of genes filtered out for low expression 

# print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
full_results <- full_results %>%
  mutate(fdr=p.adjust(`summary.p_fsoc_medulla`,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$`summary.p_fsoc_medulla`, 1e-10))  # Avoid log(0)

write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TALGeneList_TAL_cells_medulla_T2D_pooledoffset.csv")))



#Left Medulla



counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round
count_gene <- counts_path
meta_gene <- subset(so_subset)@meta.data


complete_idx <- complete.cases(meta_gene$fsoc_l_medulla)
cat("Cells with complete data:", sum(complete_idx), "\n")

# Step 3: Filter all your data to only include cells with complete predictor data
meta_gene <- meta_gene[complete_idx, ]
count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells


# Step 4: Create prediction matrix from the complete data
pred_gene <- model.matrix(~fsoc_l_medulla, data = meta_gene)


# library <- meta_gene$library_size
library <- meta_gene$pooled_offset
data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)

if (is.null(data_g_gene)) {
  data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
}

#With offset
result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                 pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)



#Make dataframe of final results
full_results <- as.data.frame(result)
#Calculate number of genes filtered out for low expression 

# print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
full_results <- full_results %>%
  mutate(fdr=p.adjust(`summary.p_fsoc_l_medulla`,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$`summary.p_fsoc_l_medulla`, 1e-10))  # Avoid log(0)

write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TALGeneList_TAL_cells_l_medulla_T2D_pooledoffset.csv")))

















#Right Medulla



counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round
count_gene <- counts_path
meta_gene <- subset(so_subset)@meta.data


complete_idx <- complete.cases(meta_gene$fsoc_r_medulla)
cat("Cells with complete data:", sum(complete_idx), "\n")

# Step 3: Filter all your data to only include cells with complete predictor data
meta_gene <- meta_gene[complete_idx, ]
count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells


# Step 4: Create prediction matrix from the complete data
pred_gene <- model.matrix(~fsoc_r_medulla, data = meta_gene)


# library <- meta_gene$library_size
library <- meta_gene$pooled_offset
data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)

if (is.null(data_g_gene)) {
  data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
}

#With offset
result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                 pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)



#Make dataframe of final results
full_results <- as.data.frame(result)
#Calculate number of genes filtered out for low expression 

# print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
full_results <- full_results %>%
  mutate(fdr=p.adjust(`summary.p_fsoc_r_medulla`,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$`summary.p_fsoc_r_medulla`, 1e-10))  # Avoid log(0)

write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TALGeneList_TAL_cells_r_medulla_T2D_pooledoffset.csv")))














##Lean Control Analysis 

remove(list=ls())






load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)



tal_genes <- c("SLC12A1", "KCNJ1", "CLCNKB", "BSND", "CLDN16", 
               "ATP1A1", "ATP1B1", "ATP1B3", "FXYD2", 
               "PPARGC1A", "NRF1", "TFAM", "COX4I1", "COX6A1", "COX7A1", "ATP5F1A", "ATP5F1B",
               "UQCRB", "UQCRC1","NDUFA4", "NDUFB5", "NDUFS1", 
               "UMOD", "CLDN10", "WNK1", "WNK4","STK39","OXSR1")




dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/FSOC_NEBULA/'









harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')

#harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat2 <- dat %>% filter(visit == 'baseline') %>% 
  filter(study %in% c('RENAL-HEIR', 'RENAL-HEIRitage', 'CROCODILE') | record_id == 'IT_19') %>%
  #  filter(group != 'Obese Control') %>% 
  dplyr::select(mrn, record_id, study, visit, group, starts_with('fsoc'), bmi, 
                epic_sglti2_1, age, sex, epic_mfm_1, epic_insulin_1, epic_glp1ra_1)

dat2 <- dat2[-which(dat2$study == 'CROCODILE' & dat2$group == 'Type 1 Diabetes'),]


tests <- c('fsoc_l_cortex', 'fsoc_r_cortex', 
           'fsoc_l_kidney', 'fsoc_r_kidney', 
           'fsoc_l_medulla', 'fsoc_r_medulla', 
           'fsoc_l_combined', 'fsoc_r_combined',
           'fsoc_medulla', 'fsoc_cortex', 'fsoc_kidney',
           'fsoc_full_combined')




tmp_df <- dat2 %>% dplyr::select(starts_with('fsoc'))



find_fsoc_averages <- function(data){
  tmp_data <- data %>% dplyr::select(starts_with('fsoc'))
  fsoc_full_combined <- rowMeans(tmp_data, na.rm=T)
  
  tmp_data <- data %>% dplyr::select(starts_with('fsoc_l_'))
  fsoc_l_combined <- tmp_data %>% rowMeans(na.rm=T)
  
  tmp_data <- data %>% dplyr::select(starts_with('fsoc_r_'))
  fsoc_r_combined <- tmp_data %>% rowMeans(na.rm=T)
  
  fsoc_medulla <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_medulla) %>%
    rowMeans(na.rm=T)
  fsoc_cortex <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_cortex) %>%
    rowMeans(na.rm=T)
  fsoc_kidney <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_kidney) %>%
    rowMeans(na.rm=T)
  
  tmp_df <- cbind(fsoc_l_combined, fsoc_r_combined, fsoc_medulla, fsoc_cortex, fsoc_kidney, fsoc_full_combined)
  return(tmp_df)
  
}

tmp_results <- find_fsoc_averages(dat2)


dat_results <- dat2 %>% bind_cols(tmp_results)


meta.data <- so_subset@meta.data
dat_results <- dat_results %>% filter(record_id %in% meta.data$record_id)



names(meta.data)[which(names(meta.data) %in% c('fsoc_l_kidney', 'fsoc_r_kidney', 
                                               'fsoc_l_medulla', 'fsoc_r_medulla', 'fsoc_l_cortex', 'fsoc_r_cortex',
                                               'fsoc_l_combined', 'fsoc_r_combined',
                                               'fsoc_medulla', 'fsoc_cortex', 'fsoc_kidney',
                                               'fsoc_full_combined'))]

meta.data <- meta.data %>% dplyr::select(-fsoc_l_kidney, -fsoc_l_medulla, 
                                         -fsoc_r_kidney, -fsoc_r_medulla, 
                                         -fsoc_l_cortex, -fsoc_r_cortex) %>% 
  left_join(dat_results, by='record_id')



so_subset@meta.data$fsoc_l_kidney <- meta.data$fsoc_l_kidney
so_subset@meta.data$fsoc_l_medulla <- meta.data$fsoc_l_medulla
so_subset@meta.data$fsoc_l_cortex <- meta.data$fsoc_l_cortex

so_subset@meta.data$fsoc_r_kidney <- meta.data$fsoc_r_kidney
so_subset@meta.data$fsoc_r_medulla <- meta.data$fsoc_r_medulla
so_subset@meta.data$fsoc_r_cortex <- meta.data$fsoc_r_cortex

so_subset@meta.data$fsoc_l_combined <- meta.data$fsoc_l_combined
so_subset@meta.data$fsoc_r_combined <- meta.data$fsoc_r_combined
so_subset@meta.data$fsoc_medulla <- meta.data$fsoc_medulla
so_subset@meta.data$fsoc_cortex <- meta.data$fsoc_cortex
so_subset@meta.data$fsoc_kidney <- meta.data$fsoc_kidney
so_subset@meta.data$fsoc_full_combined <- meta.data$fsoc_full_combined



#filter data to only TAL cells 





#gene lists

tal_genes <- c("SLC12A1", "KCNJ1", "CLCNKB", "BSND", "CLDN16", 
               "ATP1A1", "ATP1B1", "ATP1B3", "FXYD2", 
               "PPARGC1A", "NRF1", "TFAM", "COX4I1", "COX6A1", "COX7A1", "ATP5F1A", "ATP5F1B",
               "UQCRB", "UQCRC1","NDUFA4", "NDUFB5", "NDUFS1", 
               "UMOD", "CLDN10", "WNK1", "WNK4","STK39","OXSR1")




dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/FSOC_NEBULA/'



so_subset <- subset(so_subset, group == 'Lean_Control')
so_subset <- subset(so_subset, features = tal_genes)

counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round
count_gene <- counts_path
meta_gene <- subset(so_subset)@meta.data


complete_idx <- complete.cases(meta_gene$fsoc_medulla)
cat("Cells with complete data:", sum(complete_idx), "\n")

# Step 3: Filter all your data to only include cells with complete predictor data
meta_gene <- meta_gene[complete_idx, ]
count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells


# Step 4: Create prediction matrix from the complete data
pred_gene <- model.matrix(~fsoc_medulla, data = meta_gene)


# library <- meta_gene$library_size
library <- meta_gene$pooled_offset
data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)

if (is.null(data_g_gene)) {
  data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
}

#With offset
result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                 pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)



#Make dataframe of final results
full_results <- as.data.frame(result)
#Calculate number of genes filtered out for low expression 

# print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
full_results <- full_results %>%
  mutate(fdr=p.adjust(`summary.p_fsoc_medulla`,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$`summary.p_fsoc_medulla`, 1e-10))  # Avoid log(0)

write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TALGeneList_TAL_cells_medulla_LC_pooledoffset.csv")))



#Left Medulla



counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round
count_gene <- counts_path
meta_gene <- subset(so_subset)@meta.data


complete_idx <- complete.cases(meta_gene$fsoc_l_medulla)
cat("Cells with complete data:", sum(complete_idx), "\n")

# Step 3: Filter all your data to only include cells with complete predictor data
meta_gene <- meta_gene[complete_idx, ]
count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells


# Step 4: Create prediction matrix from the complete data
pred_gene <- model.matrix(~fsoc_l_medulla, data = meta_gene)


# library <- meta_gene$library_size
library <- meta_gene$pooled_offset
data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)

if (is.null(data_g_gene)) {
  data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
}

#With offset
result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                 pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)



#Make dataframe of final results
full_results <- as.data.frame(result)
#Calculate number of genes filtered out for low expression 

# print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
full_results <- full_results %>%
  mutate(fdr=p.adjust(`summary.p_fsoc_l_medulla`,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$`summary.p_fsoc_l_medulla`, 1e-10))  # Avoid log(0)

write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TALGeneList_TAL_cells_l_medulla_LC_pooledoffset.csv")))

















#Right Medulla



counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round
count_gene <- counts_path
meta_gene <- subset(so_subset)@meta.data


complete_idx <- complete.cases(meta_gene$fsoc_r_medulla)
cat("Cells with complete data:", sum(complete_idx), "\n")

# Step 3: Filter all your data to only include cells with complete predictor data
meta_gene <- meta_gene[complete_idx, ]
count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells


# Step 4: Create prediction matrix from the complete data
pred_gene <- model.matrix(~fsoc_r_medulla, data = meta_gene)


# library <- meta_gene$library_size
library <- meta_gene$pooled_offset
data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)

if (is.null(data_g_gene)) {
  data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
}

#With offset
result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                 pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)



#Make dataframe of final results
full_results <- as.data.frame(result)
#Calculate number of genes filtered out for low expression 

# print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
full_results <- full_results %>%
  mutate(fdr=p.adjust(`summary.p_fsoc_r_medulla`,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$`summary.p_fsoc_r_medulla`, 1e-10))  # Avoid log(0)

write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TALGeneList_TAL_cells_r_medulla_LC_pooledoffset.csv")))









#Graphing

remove(list=ls())

tests <- c('medulla', 'l_medulla', 'r_medulla')
dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/FSOC_NEBULA/'


for(i in c(1:length(tests))){
  
  t2d <- data.table::fread(paste0(dir.results, 'NEBULA_TALGeneList_TAL_cells_', tests[i], '_T2D_pooledoffset.csv'))
  names(t2d)[str_which(names(t2d), 'summary.logFC_[:alpha:]+')] <- 'logFC'
  names(t2d)[str_which(names(t2d), 'summary.p_[:alpha:]+')] <- 'pvalue'
  t2d <- t2d %>%
    dplyr::select(gene = summary.gene, logFC, pvalue) %>% 
    mutate(group = 'T2D')
  
  lc <- data.table::fread(paste0(dir.results, 'NEBULA_TALGeneList_TAL_cells_', tests[i], '_LC_pooledoffset.csv')) 
  names(lc)[str_which(names(lc), 'summary.logFC_[:alpha:]+')] <- 'logFC'
  names(lc)[str_which(names(lc), 'summary.p_[:alpha:]+')] <- 'pvalue'
  lc <- lc %>% 
    dplyr::select(gene = summary.gene, logFC, pvalue) %>%
    mutate(group = 'LC')
  
  full_results <- bind_rows(t2d, lc)
  
  full_results$color1 <- ifelse(full_results$pvalue < 0.05, "lightcoral", "gray")
  
  # Identify significant points (fdr < 0.05)
 
  max <- max(full_results$logFC)
  # max <- 3.1
  min <- min(full_results$logFC)
  
  if(tests[i] == 'medulla'){
    text <- 'Medulla'
    
  }else if(tests[i] == 'l_medulla'){
    text <- 'Left Medulla'
  }else if(tests[i] == 'r_medulla'){
    text <- 'Right Medulla'
  }
  
  dot_plot <- ggplot(full_results, aes(
    y = reorder(gene, logFC),
    x = logFC,
    color = color1,
    shape = group,
    size = abs(logFC)
  )) +
    geom_point(alpha = 0.7) +
    scale_color_identity(name = 'Significance (pvalue)', labels = c('> 0.05', '<0.05'), guide = 'legend')+
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed Genes in ",  text, ': TAL Cells'),
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0(
        "Participant Number: T2D: 9,", 'LC: 11',
        "; Genes = 28", 
        ", Cells = 17,902 & 12,682"
      )
    ) +
    theme(plot.caption = element_text(size = 8), 
          plot.title = element_text(hjust = 0),
          axis.text.y = element_text(size = 8),
          # axis.text.x = element_text(angle = 0, hjust = 1),
          axis.line = element_line(color = "black", size = 0.5),
          axis.ticks.x = element_line(color = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()
    )
  
  png(paste0(dir.results, 'NEBULA_FSOC_', tests[i], '_LCvsT2D_Dotplots.png'), 
      width =600, height = 800)
  print(dot_plot)
  dev.off()
  
  
}










