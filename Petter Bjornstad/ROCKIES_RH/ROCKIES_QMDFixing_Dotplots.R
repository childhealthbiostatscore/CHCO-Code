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



##Fixing the order references for LC vs. T2D 
load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

remove(so_kpmp_sc)

#dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
#dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))

#so_subset <- subset(so_subset, record_id == dat_groups$record_id)
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

#load('C:/Users/netio/Downloads/TCA_genes.txt')
#load('C:/Users/netio/Downloads/OxPhos_genes.txt')


dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/'





#Make sure exposure/independent/x variable or group variable is a factor variable
so_subset$group <- factor(so_subset$group)
#Make sure to set reference level
so_subset$group  <- relevel(so_subset$group ,ref="Lean_Control")

counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round







missing_celltypes <- c('TAL', 'C-TAL-1', 'C-TAL-2', 'dTAL', 'aTAL',
                       'DCT', 'dDCT', 'EC',
                       "EC/VSMC","EC-AVR","EC-PTC","EC-AEA","EC-LYM","EC-GC",
                       "cDC",
                         "cycT",
                         "CD4+ T",
                         "CD8+ T",
                         "NK",
                         "B",
                         "MON",
                         "MAC",
                         "MC", 
                       'POD')





for (celltype in missing_celltypes) {
  print(paste0('Starting ', celltype))
  
  if(celltype %in% c('TAL', 'EC', 'POD')){
    so_celltype <- subset(so_subset, celltype2 == celltype)
  }else if(celltype == 'DCT'){
    so_celltype <- subset(so_subset, DCT_celltype == celltype)
  }else{
    so_celltype <- subset(so_subset, KPMP_celltype == celltype)
  }
  
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or group variable is a factor variable
  so_celltype$group <- factor(so_celltype$group)
  #Make sure to set reference level
  so_celltype$group  <- relevel(so_celltype$group ,ref="Lean_Control")
  
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  

  
  #Ox_Phos Cycle
  # List of genes
  genes_list <- ox_phos_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  test2 <- test %>% filter(record_id %in% unique(so_celltype@meta.data$record_id))
  t2d_count <- test2 %>% filter(group == 'Type_2_Diabetes') %>% nrow()
  lc_count <- test2 %>% filter(group == 'Lean_Control') %>% nrow()
  
  
  start_time <- Sys.time()
  
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- subset(so_celltype,features=g)@meta.data
      pred_gene <- model.matrix(~group, data = meta_gene)
      # library <- meta_gene$library_size
      library <- meta_gene$pooled_offset
      data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
      
      if (is.null(data_g_gene)) {
        data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
      }
      
      #With offset
      result <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
      
      list(gene = g, result = result)  # return both gene name and result
      
    }, error = function(e) {
      NULL
    })
  }
  
  stopCluster(cl)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # set the names of results based on gene names
  nebula_results_list <- Filter(Negate(is.null), nebula_results_list)  # remove NULLs first
  names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)  # set names
  nebula_results_list <- lapply(nebula_results_list, function(x) x$result)  # clean list back to just results
  
  PT_nebula_converged <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      converged <- nebula_results_list[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  nebula_summaries <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      df <- nebula_results_list[[gene_name]]$summary
      df <- df %>% mutate(Gene = gene_name)
      return(df)
    }
  )
  nonconverge_genes <- unique(PT_nebula_converged$Gene[which(PT_nebula_converged$Convergence_Code==-40)]) 
  
  #Make dataframe of final results
  full_results <- as.data.frame(nebula_summaries)
  #Calculate number of genes filtered out for low expression 
  low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results <- full_results %>%
    mutate(fdr=p.adjust(`p_groupType_2_Diabetes`,method="fdr"))  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results$`p_groupType_2_Diabetes`, 1e-10))  # Avoid log(0)
  
  write.csv(full_results,fs::path(dir.results,paste0("NEBULA_OX_PHOS_cycle_",celltype2,"_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv")))
  
  print(paste0(celltype, ' is done.'))
}



























####Fixing the T2D (SGLT2i) Analyses 



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


dat <- dat %>% dplyr::select(record_id, epic_sglti2_1) %>% 
  filter(!is.na(epic_sglti2_1))

test$epic_sglti2_1 <- NULL

test <- test %>% left_join(dat, by='record_id')

so_subset@meta.data$epic_sglti2_1 <- test$epic_sglti2_1



counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round


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













missing_celltypes <- c('TAL', 'C-TAL-1', 'C-TAL-2', 'dTAL', 'aTAL',
                       'DCT', 'dDCT', 'EC',
                       "EC/VSMC","EC-AVR","EC-PTC","EC-AEA","EC-LYM","EC-GC",
                       "cDC",
                       "cycT",
                       "CD4+ T",
                       "CD8+ T",
                       "NK",
                       "B",
                       "MON",
                       "MAC",
                       "MC", 'POD')

for (celltype in missing_celltypes) {
  #Filter to PT Cells

  print(paste0('Starting ', celltype))
  
  if(celltype %in% c('TAL', 'EC', 'POD')){
    so_celltype <- subset(so_subset, celltype2 == celltype)
  }else if(celltype == 'DCT'){
    so_celltype <- subset(so_subset, DCT_celltype == celltype)
  }else{
    so_celltype <- subset(so_subset, KPMP_celltype == celltype)
  }
  
  
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or group variable is a factor variable
  so_celltype$epic_sglti2_1 <- factor(so_celltype$epic_sglti2_1)
  #Make sure to set reference level
  so_celltype$epic_sglti2_1  <- relevel(so_celltype$epic_sglti2_1 ,ref="No")
  
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  test2 <- test %>% filter(record_id %in% unique(so_celltype@meta.data$record_id))
  sglt2_count <- test2 %>% filter(epic_sglti2_1 == 'Yes') %>% nrow()
  nosglt2_count <- test2 %>% filter(epic_sglti2_1 == 'No') %>% nrow()
  
  
  start_time <- Sys.time()
  
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- subset(so_celltype,features=g)@meta.data
      pred_gene <- model.matrix(~epic_sglti2_1, data = meta_gene)
      # library <- meta_gene$library_size
      library <- meta_gene$pooled_offset
      data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
      
      if (is.null(data_g_gene)) {
        data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
      }
      
      #With offset
      result <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
      
      list(gene = g, result = result)  # return both gene name and result
      
    }, error = function(e) {
      NULL
    })
  }
  
  stopCluster(cl)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # set the names of results based on gene names
  nebula_results_list <- Filter(Negate(is.null), nebula_results_list)  # remove NULLs first
  names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)  # set names
  nebula_results_list <- lapply(nebula_results_list, function(x) x$result)  # clean list back to just results
  
  PT_nebula_converged <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      converged <- nebula_results_list[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  nebula_summaries <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      df <- nebula_results_list[[gene_name]]$summary
      df <- df %>% mutate(Gene = gene_name)
      return(df)
    }
  )
  nonconverge_genes <- unique(PT_nebula_converged$Gene[which(PT_nebula_converged$Convergence_Code==-40)]) 
  
  #Make dataframe of final results
  full_results <- as.data.frame(nebula_summaries)
  #Calculate number of genes filtered out for low expression 
  low_exp <- length(tca_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  nebula_nonconverged_percent <- paste0(round((1-(length(tca_genes)-length(nonconverge_genes))/length(tca_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results <- full_results %>%
    mutate(fdr=p.adjust(`p_epic_sglti2_1Yes`,method="fdr"))  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results$`p_epic_sglti2_1Yes`, 1e-10))  # Avoid log(0)
  
  write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TCA_cycle_",celltype2,"_cells_T2D_SGLT2_unadjusted_pooled_offset.csv")))
  
 
  
  #Ox_Phos Cycle
  # List of genes
  genes_list <- ox_phos_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  test2 <- test %>% filter(record_id %in% unique(so_celltype@meta.data$record_id))
  sglt2_count <- test2 %>% filter(epic_sglti2_1 == 'Yes') %>% nrow()
  nosglt2_count <- test2 %>% filter(epic_sglti2_1 == 'No') %>% nrow()
  
  
  start_time <- Sys.time()
  
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- subset(so_celltype,features=g)@meta.data
      pred_gene <- model.matrix(~epic_sglti2_1, data = meta_gene)
      # library <- meta_gene$library_size
      library <- meta_gene$pooled_offset
      data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
      
      if (is.null(data_g_gene)) {
        data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
      }
      
      #With offset
      result <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
      
      list(gene = g, result = result)  # return both gene name and result
      
    }, error = function(e) {
      NULL
    })
  }
  
  stopCluster(cl)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # set the names of results based on gene names
  nebula_results_list <- Filter(Negate(is.null), nebula_results_list)  # remove NULLs first
  names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)  # set names
  nebula_results_list <- lapply(nebula_results_list, function(x) x$result)  # clean list back to just results
  
  PT_nebula_converged <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      converged <- nebula_results_list[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  nebula_summaries <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      df <- nebula_results_list[[gene_name]]$summary
      df <- df %>% mutate(Gene = gene_name)
      return(df)
    }
  )
  nonconverge_genes <- unique(PT_nebula_converged$Gene[which(PT_nebula_converged$Convergence_Code==-40)]) 
  
  #Make dataframe of final results
  full_results <- as.data.frame(nebula_summaries)
  #Calculate number of genes filtered out for low expression 
  low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results <- full_results %>%
    mutate(fdr=p.adjust(`p_epic_sglti2_1Yes`,method="fdr"))  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results$`p_epic_sglti2_1Yes`, 1e-10))  # Avoid log(0)
  
  write.csv(full_results,fs::path(dir.results,paste0("NEBULA_OX_PHOS_cycle_",celltype2,"_cells_T2D_SGLT2_unadjusted_pooled_offset.csv")))
  
  names(full_results)[which(names(full_results) == 'p_epic_sglti2_1Yes')] <- 'pvalue' 
  full_results$color1 <- ifelse(full_results$fdr < 0.05, "lightcoral", "gray")
  full_results$color2 <- ifelse(full_results$pvalue < 0.05, "lightcoral", "gray")
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_epic_sglti2_1Yes`)
  # max <- 3.1
  min <- min(full_results$`logFC_epic_sglti2_1Yes`)
  
  
  order <- data.table::fread(paste0("C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/NEBULA_OX_PHOS_cycle_", 
                                    celltype2, '_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv')) %>%
    arrange(`logFC_groupType_2_Diabetes`)
  
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y = factor(gene, levels = order$gene),
    x = `logFC_epic_sglti2_1Yes`,
    color = color1,
    size = abs(`logFC_epic_sglti2_1Yes`)
  )) +
    geom_point(alpha = 0.7) +
    scale_color_identity(name = 'Significance (FDR)', labels = c('> 0.05', '<0.05'), guide = 'legend')+
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "SGLT2i vs. No SGLT2i (T2D Only), Unadjusted (Pooled Offset)",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0(       "Participants on SGLT2: ", sglt2_count, ', No SGLT2: ', nosglt2_count,
                              "; Genes = ", Genes,
                              ", Cells = ", Cell,
                              ", Non-Convergence Rate: ", Nonconvergence_Rate,
                              ", Genes Filtered out for Low Expression: ", low_exp
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
  
  
  names(full_results)[which(names(full_results) == 'p_epic_sglti2_1Yes')] <- 'pvalue' 
  dot_plot_ox_phos_2 <- ggplot(full_results, aes(
    y = factor(gene, levels = order$gene),
    x = `logFC_epic_sglti2_1Yes`,
    color = color2,
    size = abs(`logFC_epic_sglti2_1Yes`)
  )) +
    geom_point(alpha = 0.7) +
    scale_color_identity(name = 'Significance (P-value)', labels = c('> 0.05', '<0.05'), guide = 'legend') + 
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OxPhos Cycle Genes in ",celltype," Cells"),
      subtitle = "SGLT2i vs. No SGLT2i (T2D Only), Unadjusted (Pooled Offset)",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0("Participants on SGLT2: ", sglt2_count, ', No SGLT2: ', nosglt2_count, 
                       "; Genes = ", Genes,
                       ", Cells = ", Cell,
                       ", Non-Convergence Rate: ", Nonconvergence_Rate,
                       ", Genes Filtered out for Low Expression: ", low_exp
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
  

  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
  print(paste0(celltype, ' is done.'))
}
















##If data is there, just need to plot 








missing_celltypes <- c('TAL', 'C-TAL-1', 'C-TAL-2', 'dTAL', 'aTAL',
                       'DCT', 'dDCT', 'EC',
                       "EC/VSMC","EC-AVR","EC-PTC","EC-AEA","EC-LYM","EC-GC",
                       "cDC",
                       "cycT",
                       "CD4+ T",
                       "CD8+ T",
                       "NK",
                       "B",
                       "MON",
                       "MAC",
                       "MC", 'POD')

cell_counts <- data.table::fread("C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/Cell_Counts.csv")

for (celltype in missing_celltypes) {
  #Filter to PT Cells
  
  print(paste0('Starting ', celltype))
  
  if(celltype %in% c('TAL', 'EC', 'POD')){
    so_celltype <- subset(so_subset, celltype2 == celltype)
  }else if(celltype == 'DCT'){
    so_celltype <- subset(so_subset, DCT_celltype == celltype)
  }else{
    so_celltype <- subset(so_subset, KPMP_celltype == celltype)
  }
  
  
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
full_results <- data.table::fread(paste0('C:/Users/'))
  
  print(paste0(celltype, ' is done.'))
}











