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





#insulin clamp integration


dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Insulin_Clamp/'


insulin_genes <- c("INSR", "IGF1R", "INSRR",
                     "IRS1", "IRS2", "IRS4", "GRB2", "SHC1",
                     "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2",
  "AKT1", "AKT2", "AKT3",
  "PDPK1", "PTEN", "INPP4B",
  "MTOR", "RPTOR", "RICTOR", "MLST8",
  "AKT1S1", "TSC1", "TSC2", "RHEB",
  "EIF4EBP1", "RPS6KB1", "RPS6", "EIF4E",
  "SLC2A1", "SLC2A2", "SLC2A4",
  "SLC5A1", "SLC5A2",
  "HK1", "HK2", "GCK",
  "PFKFB2", "PFKFB3", "PKM", "LDHA",
  "PCK1", "PCK2", "G6PC", "FBP1", "FBP2",
  "FASN", "ACACA", "ACACB", "SCD",
  "SREBF1", "SREBP2", "INSIG1", "INSIG2",
  "CPT1A", "CPT2", "ACOX1",
  "PPARGC1A", "PPARGC1B", "NRF1", "TFAM",
  "ATP5F1A", "COX4I1", "NDUFS1", "SDHA",
  "ATP5F1A", "COX4I1", "NDUFS1", "SDHA",
  "PTPN1", "PTPRF", "PTPN11",
  "SOCS1", "SOCS3", "SOCS7",
  "PRKAA1", "PRKAA2",
  "FOXO1", "FOXO3", "FOXO4")




#load in correct data, ensure FSOC accuracy. Only Type 2's 


load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')



so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

test <- so_subset@meta.data

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat2 <- dat2 %>% dplyr::select(record_id, epic_sglti2_1, airg, acprg) %>% 
  filter(!is.na(epic_sglti2_1))

#Fix up FSOC data
test$airg <- NULL
test$acprg<- NULL



test <- test %>% left_join(dat2, by='record_id')

so_subset@meta.data$airg <- test$airg
so_subset@meta.data$acprg <- test$acprg




demo_table <- dat %>% filter(visit == 'baseline') %>% 
  filter(record_id %in% so_subset@meta.data$record_id) %>% 
  filter(!is.na(airg))

demo_table$group2 <- ifelse(demo_table$epic_sglti2_1 == 'Yes', 'T2D-SGLT2', 'T2D-No SGLT2')


table1::table1(~age + sex + bmi + airg + acprg + epic_sglti2_1 | group2, data = demo_table)




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





#TAL 
so_subset <- subset(so_subset, celltype2 == 'TAL')



tests <- c('airg', 'acprg')




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

celltypes <- c('PT', 'TAL', 'EC')




#gene lists


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






























#PT 








#TAL 






#EC 













