#All Dotplots Function
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






















NEBULA_LC_T2D_dotplots <- function(so_subset, celltype, dir.results){
  
  if(celltype == 'All'){
   so_celltype <- so_subset 
  }else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype != 'DCTall'){
    so_celltype <- subset(so_subset,KPMP_celltype==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype==celltype)
  }
  
  
  
  
  
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or group variable is a factor variable
  so_celltype$group <- factor(so_celltype$group)
  #Make sure to set reference level
  so_celltype$group  <- relevel(so_celltype$group ,ref="Lean_Control")
  
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
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
  low_exp <- length(tca_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  nebula_nonconverged_percent <- paste0(round((1-(length(tca_genes)-length(nonconverge_genes))/length(tca_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results <- full_results %>%
    mutate(fdr=p.adjust(`p_groupType_2_Diabetes`,method="fdr"))  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results$`p_groupType_2_Diabetes`, 1e-10))  # Avoid log(0)
  
  write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TCA_cycle_",celltype2,"_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv")))
  
  names(full_results)[6] <- 'pvalue' 
  full_results$color1 <- ifelse(full_results$fdr < 0.05, "lightcoral", "gray")
  full_results$color2 <- ifelse(full_results$pvalue < 0.05, "lightcoral", "gray")
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  # full_results$color3 <- ifelse(full_results$fdr3 < 0.2 & full_results$`logFC_SGLT2SGLT2i`3 > 0, "lightcoral",
  #                               ifelse(full_results$fdr3 < 0.2 & full_results$`logFC_SGLT2SGLT2i`3 < 0, "lightblue", "gray"))
  # 
  # # Identify significant points (fdr < 0.05)
  # significant_df3 <- full_results[full_results$fdr3 < 0.2, ]
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_tca <- ggplot(full_results, aes(
    y = reorder(gene, `logFC_groupType_2_Diabetes`),
    x = `logFC_groupType_2_Diabetes`,
    color = color1,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
    scale_color_identity(name = 'Significance (FDR)', labels = c('> 0.05', '<0.05'), guide = 'legend')+
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "LC vs. T2D (No SGLT2), Unadjusted (Pooled Offset)",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0("Participant Number: T2D: ", t2d_count, ', LC: ', lc_count, 
                       "; Genes = ", Genes,
                       ", Cells = ", Cell
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
  
  names(full_results)[6] <- 'pvalue' 
  dot_plot_tca_2 <- ggplot(full_results, aes(
    y = reorder(gene, `logFC_groupType_2_Diabetes`),
    x = `logFC_groupType_2_Diabetes`,
    color = color2,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
    scale_color_identity(name = 'Significance (P-value)', labels = c('> 0.05', '<0.05'), guide = 'legend')+
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "LC vs. T2D (No SGLT2), Unadjusted (Pooled Offset)",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0("Participant Number: T2D: ", t2d_count, ', LC: ', lc_count, 
                       "; Genes = ", Genes,
                       ", Cells = ", Cell
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
  
  names(full_results)[6] <- 'pvalue' 
  full_results$color1 <- ifelse(full_results$fdr < 0.05, "lightcoral", "gray")
  full_results$color2 <- ifelse(full_results$pvalue < 0.05, "lightcoral", "gray")
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y = reorder(gene, `logFC_groupType_2_Diabetes`),
    x = `logFC_groupType_2_Diabetes`,
    color = color1,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
    scale_color_identity(name = 'Significance (FDR)', labels = c('> 0.05', '<0.05'), guide = 'legend')+
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "LC vs. T2D (No SGLT2), Unadjusted (Pooled Offset)",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0(       "Participant Number: T2D: ", t2d_count, ', LC: ', lc_count,
                              "; Genes = ", Genes,
                              ", Cells = ", Cell
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

  
  names(full_results)[6] <- 'pvalue' 
  dot_plot_ox_phos_2 <- ggplot(full_results, aes(
    y = reorder(gene, `logFC_groupType_2_Diabetes`),
    x = `logFC_groupType_2_Diabetes`,
    color = color2,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
    scale_color_identity(name = 'Significance (P-value)', labels = c('> 0.05', '<0.05'), guide = 'legend') + 
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OxPhos Cycle Genes in ",celltype," Cells"),
      subtitle = "LC vs. T2D (No SGLT2), Unadjusted (Pooled Offset)",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0("Participant Number: T2D: ", t2d_count, ', LC: ', lc_count, 
                       "; Genes = ", Genes,
                       ", Cells = ", Cell
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
  

  # dot_plot
  # 
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
  
  
}









NEBULA_T2D_SGLT2_dotplots <- function(data, celltype, dir.results){
  
  if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype != 'DCT'){
    so_celltype <- subset(so_subset,KPMP_celltype==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'DCT'){
    so_celltype <- subset(so_subset, DCT_celltype==celltype)
  }
  
  
  
  
  
}










####Lean control analysis
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




cell_vector <- c('All', 
                 'PT', "PT-S1/S2","PT-S3","aPT", 
                 'TAL', "C-TAL-1","C-TAL-2","aTAL","dTAL", 
                 'DCTall', "DCT","dDCT", 
                 'EC', "EC/VSMC","EC-AVR","EC-PTC","EC-AEA","EC-LYM","EC-GC",
                 'POD',"cDC","cycT","CD4+ T", "CD8+ T","NK","B","MON", "MAC","MC")


for(i in c(2:length(cell_vector))){
  NEBULA_LC_T2D_dotplots(so_subset, celltype = cell_vector[i], dir.results = dir.results)
  print(cell_vector[i])
}














#### T2D analysis 
















