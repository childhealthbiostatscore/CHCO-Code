#### Sex-Specific Analyses in Lean Controls




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
#library(table1)
library(clusterProfiler)
library('org.Hs.eg.db')







load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

#dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
#dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))

#so_subset <- subset(so_subset, record_id == dat_groups$record_id)
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))








dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/'


so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Lean_Control')
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))






harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )





## demographics tables 

dat_small <- dat %>% filter(record_id %in% test$record_id)



library(gt)
library(gtsummary)

desc_table1_fixed <- dat_small %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study) %>%
  tbl_summary(
    by = sex,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      race_ethnicity ~ "categorical",
      study ~ "categorical"
      
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave(paste0(dir.results, "LC_demographics.png"), 
         vwidth = 1200, vheight = 800)












#### Celltype Analyses




celltypes_vec <- c('All', 'PT', 'TAL', 'EC', 'DCTall', 'IC')



#function
LC_NEBULA_Analysis <- function(so_subset, dir.results, celltype, genes){
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
  }else if(celltype == 'IC'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c(
      "cDC",
      "cycT",
      "CD4+ T",
      "CD8+ T",
      "NK",
      "B",
      "MON",
      "MAC",
      "MC"))
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

celltypes_vec <- c('All', 'PT', 'TAL', 'EC', 'IC')

for(i  in 1:length(celltypes_vec)){
  LC_NEBULA_Analysis(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}





















##### GSEA





















































