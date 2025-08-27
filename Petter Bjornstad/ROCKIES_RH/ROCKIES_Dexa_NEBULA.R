#DEXA Analysis
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




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat2 <- dat %>% dplyr::select(record_id, mrn, visit, epic_sglti2_1, starts_with('dexa_'))
dat2 <- dat2 %>% filter(!duplicated(record_id))

#dat2 <- dat2 %>% filter(!is.na(dexa_body_fat))






dex_var <- c("dexa_ag_ratio",  "dexa_body_fat",            "dexa_bone_mineral_density", "dexa_est_vat",             
             "dexa_fat_kg",               "dexa_lean_kg",              "dexa_lean_mass",            "dexa_trunk_kg",            
            "dexa_trunk_mass"  )




gene_list <- c(
      "LEPR", "ADIPOR1", "ADIPOR2",
      "PPARG", "PPARA", "PPARD",
      "SREBF1", "SREBF2",
      "STAT3", "JAK2", "SOCS3",
      "TNF", "IL6", "IL1B", "CCL2",
      "NFKB1", "RELA", "IKBKB",
      "TLR2", "TLR4", "MYD88",
      "NLRP3", "CASP1", "IL18",
      "FASN", "ACACA", "SCD",
      "CPT1A", "CPT2", "ACOX1",
      "LDLR", "SCARB1", "ABCA1",
      "APOE", "APOB", "LPL",
      "DGAT1", "DGAT2", "PLIN1",
      "TGFB1", "TGFB2", "SMAD2", "SMAD3",
      "COL1A1", "COL3A1", "COL4A1",
      "FN1", "CTGF", "ACTA2",
      "MMP2", "MMP9", "TIMP1",
      "NOX1", "NOX2", "NOX4",
      "SOD1", "SOD2", "CAT", "GPX1",
      "NRF2", "KEAP1", "NQO1",
      "HMOX1", "GCLC", "GCLM",
      "HMOX1", "GCLC", "GCLM",
      "IRS1", "IRS2", "INSR",
      "GLUT4", "GLUT2",
      "SOCS1", "SOCS3", "PTPN1",
      "PRKAA1", "PRKAA2",
      "HSPA5", "EIF2AK3", "ERN1", "ATF6", # UPR sensors
      "XBP1", "ATF4", "DDIT3",
      "EDEM1", "DNAJB9", "HYOU1"
    )



#cell-type specific analyses 


 
"PT" = c(
                        # Fatty acid uptake and metabolism
                        "CD36", "FABP1", "FABP3", "SLC27A2", "SLC27A4",
                        # Glucose/lipid interaction
                        "PPARGC1A", "FOXO1", "SIRT1", "AMPK",
                        # PT-specific transporters affected by lipids
                        "SLC5A2", "SLC9A3", "SLC34A1"
                      )
              "POD" = c(
                        # Podocyte lipotoxicity
                        "NPHS1", "NPHS2", "PODXL", "SYNPO", "CD2AP",
                        # Cholesterol efflux
                        "ABCA1", "ABCG1", "APOE",
                      "CHOP", "BIP", "XBP1s"

)
"Endothelial" = c(
  # Endothelial dysfunction
  "VCAM1", "ICAM1", "SELE", "SELP",
  # NO synthesis
  "NOS3", "GUCY1A3", "GUCY1B3",
  # Angiogenesis
  "VEGFA", "VEGFR2", "ANGPT2"
)
"Immune" = c(
  # Macrophage polarization
  "CD68", "CD86", "CD163", "CD206",
  # Cytokines
  "IL10", "IL12", "TGFB1",
  # Chemokines
  "CCL2", "CCL3", "CCL5", "CXCL10"
)
"TAL" = c(
  # Metabolic adaptation
  "PPARGC1A", "NRF1", "TFAM",
  # Transport regulation
  "SLC12A1", "KCNJ1", "ATP1A1")

"DCT" = c(
  # Calcium handling
  "SLC12A3", "TRPM6", "CALB1",
  # Metabolic genes
  "HIF1A", "EPAS1"
)







#LC vs. T2D (no SLGT2)



load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

#dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
#dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))

#so_subset <- subset(so_subset, record_id == dat_groups$record_id)
test <- so_subset@meta.data

test <- test %>% dplyr::select(-starts_with('dexa_'))


test <- test %>% left_join(dat2, by='record_id') %>% 
  dplyr::select(record_id, starts_with('dexa_'))


so_subset@meta.data$dexa_ag_ratio <- test$dexa_ag_ratio
so_subset@meta.data$dexa_body_fat <- test$dexa_body_fat
so_subset@meta.data$dexa_bone_mineral_density <- test$dexa_bone_mineral_density
so_subset@meta.data$dexa_est_vat <- test$dexa_est_vat
so_subset@meta.data$dexa_fat_kg <- test$dexa_fat_kg
so_subset@meta.data$dexa_lean_kg <- test$dexa_lean_kg
so_subset@meta.data$dexa_lean_mass <- test$dexa_lean_mass
so_subset@meta.data$dexa_trunk_kg <- test$dexa_trunk_kg
so_subset@meta.data$dexa_trunk_mass <- test$dexa_trunk_mass


#Make sure exposure/independent/x variable or group variable is a factor variable
so_subset$group <- factor(so_subset$group)
#Make sure to set reference level
so_subset$group  <- relevel(so_subset$group ,ref="Lean_Control")

counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round


dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/dexa/'

#function
T2D_LC_Dexa_Analysis <- function(so_subset, dir.results, celltype, genes){
  if(celltype == 'All'){
    so_celltype <- so_subset 
  }else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype != 'DCTall'){
    so_celltype <- subset(so_subset,KPMP_celltype==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype=='DCT')
  }
  
  
  
  
  
  if(celltype %in% c('PT', 'PT-S1/S2', 'PT-S3', 'aPT')){
    tmp_gene_list <- c(genes, PT)
  }else if(celltype == 'POD'){
    tmp_gene_list <- c(genes, POD)
  }else if(celltype %in% c('EC', 'EC-AEA', 'EC-AVR', 'EC-GC', 'EC-PTC')){
    tmp_gene_list <- c(genes, Endothelial)
  }else if(celltype %in% c('B')){
    tmp_gene_list <- c(genes, Immune)
  }else if(celltype %in% c('TAL', 'C-TAL-1', 'C-TAL-2', 'dTAL')){
    tmp_gene_list <- c(genes, TAL)
  }else if(celltype == 'DCT'){
    tmp_gene_list <- c(genes, DCT)
  }else{
    tmp_gene_list <- genes
  }
  
  
  so_celltype <- subset(so_celltype, features = tmp_gene_list)
  
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  

  so_celltype$group <- factor(so_celltype$group)
  so_celltype$group  <- relevel(so_celltype$group ,ref="Lean_Control")
  
  
  print(paste0(celltype2, ' is running.'))
  #FOR LOOP OVER VARIABLES 
  for(iter in c(1:length(dex_var))){
    print(paste0('Working on ', dex_var[iter]))
  
    so_celltype@meta.data$Variable <- so_celltype@meta.data[,which(names(so_celltype@meta.data) == dex_var[iter])]
    counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
    count_gene <- counts_path
    meta_gene <- subset(so_celltype)@meta.data
  
  
  
  complete_idx <- complete.cases(meta_gene$Variable)
  cat("Cells with complete data:", sum(complete_idx), "\n")
  
  # Step 3: Filter all your data to only include cells with complete predictor data
  meta_gene <- meta_gene[complete_idx, ]
  count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells
  
  if(length(unique(meta_gene$group)) != 2){
    print(paste0('Skipped ', dex_var[iter], ' in ', celltype2, ' Cells.'))
    next
  }
  
  
  num_cells <- nrow(meta_gene)
  num_part <- unique(meta_gene$record_id) %>% length()
  
  tmp_df <- meta_gene %>% dplyr::select(record_id, group, epic_sglti2_1) %>% filter(!duplicated(record_id))
  
  num_t2d <- tmp_df %>% filter(group == 'Type_2_Diabetes') %>% nrow()
  num_lc <- tmp_df %>% filter(group == 'Lean_Control') %>% nrow()
  
  # Step 4: Create prediction matrix from the complete data
  pred_gene <- model.matrix(~Variable*group, data = meta_gene)
  
  
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
  
  
  full_results$num_cells <- num_cells
  full_results$num_t2d <- num_t2d
  full_results$num_lc <- num_lc
  
  write.table(full_results,paste0(dir.results,"NEBULA_", 
                                  celltype2, "_cells_", dex_var[iter], "_LC_vs_T2D_pooledoffset.csv"),
              row.names=F, quote=F, sep=',')
  }
  
  
  
  
  print(paste0(celltype2, ' is done.'))
  
  
}

celltypes_vec <- c('All', 'PT', 'PT-S1/S2', 'PT-S3', 'aPT', 'POD', 
               'TAL', 'C-TAL-1','C-TAL-2', 'dTAL', 'DCT', 'dDCT',
               'EC', 'EC-AEA', 'EC-AVR', 'EC-GC', 'EC-PTC', 
               "cDC",
               "cycT",
               #              "CD4+ T",
               #             "CD8+ T",
               "NK",
               "B",
               "MON",
               "MAC",
               "MC")

for(i  in 1:length(celltypes_vec)){
  T2D_LC_Dexa_Analysis(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}












#T2D SGLT2 Analysis 







load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

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




test <- so_subset@meta.data

test <- test %>% dplyr::select(-starts_with('dexa_'))


test <- test %>% left_join(dat2, by='record_id') %>% 
  dplyr::select(record_id, starts_with('dexa_'), group, epic_sglti2_1.y)
test$group2 <- 'T2D-No SGLT2'
test$group2[which(test$epic_sglti2_1.y == 'Yes')] <- 'T2D-SGLT2'


so_subset@meta.data$dexa_ag_ratio <- test$dexa_ag_ratio
so_subset@meta.data$dexa_body_fat <- test$dexa_body_fat
so_subset@meta.data$dexa_bone_mineral_density <- test$dexa_bone_mineral_density
so_subset@meta.data$dexa_est_vat <- test$dexa_est_vat
so_subset@meta.data$dexa_fat_kg <- test$dexa_fat_kg
so_subset@meta.data$dexa_lean_kg <- test$dexa_lean_kg
so_subset@meta.data$dexa_lean_mass <- test$dexa_lean_mass
so_subset@meta.data$dexa_trunk_kg <- test$dexa_trunk_kg
so_subset@meta.data$dexa_trunk_mass <- test$dexa_trunk_mass
so_subset@meta.data$group2 <- test$group2


#Make sure exposure/independent/x variable or group variable is a factor variable
so_subset$group2 <- factor(so_subset$group2)
#Make sure to set reference level
so_subset$group2  <- relevel(so_subset$group2 ,ref="T2D-No SGLT2")

counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round


dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/dexa/'

#function
T2D_SGLT2_Dexa_Analysis <- function(so_subset, dir.results, celltype, genes){
  if(celltype == 'All'){
    so_celltype <- so_subset 
  }else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype != 'DCTall'){
    so_celltype <- subset(so_subset,KPMP_celltype==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype=='DCT')
  }
  
  
  
  
  
  if(celltype %in% c('PT', 'PT-S1/S2', 'PT-S3', 'aPT')){
    tmp_gene_list <- c(genes, PT)
  }else if(celltype == 'POD'){
    tmp_gene_list <- c(genes, POD)
  }else if(celltype %in% c('EC', 'EC-AEA', 'EC-AVR', 'EC-GC', 'EC-PTC')){
    tmp_gene_list <- c(genes, Endothelial)
  }else if(celltype %in% c('B')){
    tmp_gene_list <- c(genes, Immune)
  }else if(celltype %in% c('TAL', 'C-TAL-1', 'C-TAL-2', 'dTAL')){
    tmp_gene_list <- c(genes, TAL)
  }else if(celltype == 'DCT'){
    tmp_gene_list <- c(genes, DCT)
  }else{
    tmp_gene_list <- genes
  }
  
  
  so_celltype <- subset(so_celltype, features = tmp_gene_list)
  
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  
  print(paste0(celltype2, ' is running.'))
  #FOR LOOP OVER VARIABLES 
  for(iter in c(1:length(dex_var))){
    print(paste0('Working on ', dex_var[iter]))
    
    so_celltype@meta.data$Variable <- so_celltype@meta.data[,which(names(so_celltype@meta.data) == dex_var[iter])]
    counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
    count_gene <- counts_path
    meta_gene <- subset(so_celltype)@meta.data
    
    
    
    complete_idx <- complete.cases(meta_gene$Variable)
    cat("Cells with complete data:", sum(complete_idx), "\n")
    
    # Step 3: Filter all your data to only include cells with complete predictor data
    meta_gene <- meta_gene[complete_idx, ]
    count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells
    
    if(length(unique(meta_gene$group2)) < 2){
      print(paste0('Skipped ', dex_var[iter], ' in ', celltype2, ' Cells.'))
      next
    }
    
    
    num_cells <- nrow(meta_gene)
    num_part <- unique(meta_gene$record_id) %>% length()
    
    tmp_df <- meta_gene %>% dplyr::select(record_id, group2, epic_sglti2_1) %>% filter(!duplicated(record_id))
    
    num_sglt2 <- tmp_df %>% filter(group2 == 'T2D-SGLT2') %>% nrow()
    num_nosglt2 <- tmp_df %>% filter(group2 == 'T2D-No SGLT2') %>% nrow()
    
    # Step 4: Create prediction matrix from the complete data
    pred_gene <- model.matrix(~Variable*group2, data = meta_gene)
    
    
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
    
    
    full_results$num_cells <- num_cells
    full_results$num_sglt2<- num_sglt2
    full_results$num_nosglt2 <- num_nosglt2
    
    write.table(full_results,paste0(dir.results,"NEBULA_", 
                                    celltype2, "_cells_", dex_var[iter], "_T2D_SGLT2_pooledoffset.csv"),
                row.names=F, quote=F, sep=',')
  }
  
  
  
  
  print(paste0(celltype2, ' is done.'))
  
  
}

celltypes_vec <- c('All', 'PT', 'PT-S1/S2', 'PT-S3', 'aPT', 'POD', 
                   'TAL', 'C-TAL-1','C-TAL-2', 'dTAL', 'DCT', 'dDCT',
                   'EC', 'EC-AEA', 'EC-AVR', 'EC-GC', 'EC-PTC', 
                   "cDC",
                   "cycT",
                   #              "CD4+ T",
                   #             "CD8+ T",
                   "NK",
                   "B",
                   "MON",
                   "MAC",
                   "MC")

for(i  in 1:length(celltypes_vec)){
  T2D_SGLT2_Dexa_Analysis(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}



