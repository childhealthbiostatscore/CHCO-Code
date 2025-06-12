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










load('C:/Users/netio/Documents/UofW/Rockies/Line4875_Rockies.RData')



kidneyimaging_analysis <- function(celltype, genes, gene_list_name = 'TCA', median = F, adjustment = NULL){
  if(median == F){
  k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
  }else{
    k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
  }
  
  if(celltype == 'PT'){
  so_celltype <- subset(so_subset,celltype2 == celltype)
  cat('PT Cells')
  }else if(celltype == 'TAL'){
    so_celltype <- subset(so_subset, TAL_celltype == celltype)
    cat('TAL Cells')
  }else if(celltype == 'DCT'){
    so_celltype <- subset(so_subset, DCT_celltype == celltype)
    cat('DCT Cells')
  }else{
    so_celltype <- subset(so_subset, KPMP_celltype == celltype)
    cat('Other celltypes')
  }
  # so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
  k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
  so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #4926 PT cells
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- genes

  total_results <- data.frame()
  for (exposure in k2_vars) {
    cl <- makeCluster(1)
    registerDoParallel(cl)
    
    start_time <- Sys.time()
    nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
      tryCatch({
        count_gene <- counts_path[g, , drop = FALSE]
        meta_gene <- subset(so_celltype,features=g)@meta.data
        
        if(!is.null(adjustment)){
          tmp.formula <- as.formula(paste0('~', exposure, '+', adjustment))
        }else{
          tmp.formula <- as.formula(paste0('~', exposure))
        }
        
        pred.formula <- as.formula(tmp.formula)
        pred_gene <- model.matrix(pred.formula, data = meta_gene)
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
        list(gene = g, summary = NA, overdispersion = NA, convergence = NA, algorithm = NA, covariance = NA, random_effect = NA)
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
    
    
    if(is.null(adjustment)){
      colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept",
                                  "p_value","gene_id","gene","Gene")
    }else{
      colnames(full_results) <- c("logFC_Intercept","logFC",'logFC_SGLTi2', "se_Intercept",
                                  "se",'se_SGLTi2', "p_Intercept","p_value",'p_value_SGLTi2', "gene_id","gene",
                                  "Gene")
    }
    
    
    
    full_results$Variable <- exposure
    #Calculate number of genes filtered out for low expression 
    full_results$low_exp <- length(tca_genes)-length(full_results$gene)
    #Filter out non-converging genes
    full_results <- full_results %>% 
      filter(!gene %in%  nonconverge_genes)
    #Calculate nonconvergence rate
    full_results$nebula_nonconverged_percent <- paste0(round((1-(length(tca_genes)-length(nonconverge_genes))/length(tca_genes))*100,3),"%")
    # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
    # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
    full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
    # mutate(fdr3=p.adjust(PValue3,method="fdr"))
    full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
    total_results <- rbind(total_results,full_results)
  }
  
  if(is.null(adjustment)){
    file.name <- paste0('C:/Users/netio/Documents/UofW/Rockies/NEBULA_', gene_list_name, 
                        '_median', median, '_', celltype, '_cells_PET_unadjusted_pooled_offset_T2D.txt')
  }else{
    file.name <- paste0('C:/Users/netio/Documents/UofW/Rockies/NEBULA_', gene_list_name, 
                        '_median', median, '_', celltype, '_cells_PET_adjusted_pooled_offset_T2D.txt')
  }
  
  
  write.table(total_results,file.name, row.names=F, quote=F, sep='\t')
  # total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
  
  # Define significance stars
  total_results <- total_results %>%
    mutate(signif = case_when(
      fdr < 0.01 ~ "**",
      fdr < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  # Select only the needed columns and rename LogFC for clarity
  if(median == T){
  subtitle1 <- paste0(celltype, ' (Median), ')
  }else{
    subtitle1 <- paste0(celltype, ', ')
  }
  
  if(!is.null(adjustment)){
    subtitle2 <- paste0('with adjustment for ', adjustment)
  }else{
    subtitle2 <- paste0('with no adjustment')
  }
  subtitle <- paste0(subtitle1, subtitle2)
  
  heatmap_data <- total_results %>%
    dplyr::select(Gene, Variable, logFC, signif)
  
  custom_order <- c("avg_c_k2_med", "avg_m_k2_med", "avg_c_f_med", "avg_m_f_med", "avg_c_k2_f_med", "avg_m_k2_f_med")
  heatmap_data$Variable <- factor(heatmap_data$Variable, levels = custom_order)
  custom_labels <- c("Average Cortical K2","Average Medulla K2","Average Cortical F","Average Medulla F",
                     "Average Cortical K2/F","Average Medulla K2/F")
  
  heat_map_p <- ggplot(heatmap_data, aes(x = Variable, y = Gene, fill = logFC)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = signif), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "#264653", mid = "white", high = "darkred",
      midpoint = 0,
      name = "LogFC"
    ) +
    theme_minimal() +
    labs(title = paste0(gene_list_name, " Genes vs. PET Variables (T2D)"),
         subtitle = subtitle,
         x = "Exposure",
         y = "Gene") +
    scale_x_discrete(labels = setNames(custom_labels, custom_order))+
    theme(
      text = element_text(face="bold"),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank()
    )
  
  # custom_colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "darkred")

  if(is.null(adjustment)){
  file.name <- paste0('C:/Users/netio/Documents/UofW/Rockies/NEBULA_Heatmap_', gene_list_name, 
                      'NEBULA_','median', median, '_', celltype, '_PET_unadjusted_pooled_offset_T2D.png')
  }else{
    file.name <- paste0('C:/Users/netio/Documents/UofW/Rockies/NEBULA_Heatmap_', gene_list_name, 
                        'NEBULA_', 'median', median, '_', celltype, '_PET_adjusted_pooled_offset_T2D.png')  
  }
  png(file.name, 
      width = 1500, height = 2000, res = 300)
  print(heat_map_p)
  dev.off()
  
}

#Finishing PT Analysis
kidneyimaging_analysis('PT', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

#PT Subtypes
kidneyimaging_analysis('PT-S1/S2', median = F, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('PT-S1/S2', median = F, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

kidneyimaging_analysis('PT-S1/S2', median = T, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('PT-S1/S2', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

#S3
kidneyimaging_analysis('PT-S3', median = F, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('PT-S3', median = F, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

kidneyimaging_analysis('PT-S3', median = T, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('PT-S3', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

#aPT
kidneyimaging_analysis('aPT', median = F, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('aPT', median = F, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

kidneyimaging_analysis('aPT', median = T, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('aPT', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

#TAL
kidneyimaging_analysis('TAL', median = F, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('TAL', median = F, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

kidneyimaging_analysis('TAL', median = T, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('TAL', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

#C-TAL-1
kidneyimaging_analysis('C-TAL-1', median = F, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('C-TAL-1', median = F, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

kidneyimaging_analysis('C-TAL-1', median = T, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('C-TAL-1', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

#C-TAL-2
kidneyimaging_analysis('C-TAL-2', median = F, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('C-TAL-2', median = F, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

kidneyimaging_analysis('C-TAL-2', median = T, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('C-TAL-2', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

#dTAL
kidneyimaging_analysis('dTAL', median = F, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('dTAL', median = F, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

kidneyimaging_analysis('dTAL', median = T, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('dTAL', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

#DCT
kidneyimaging_analysis('DCT', median = F, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('DCT', median = F, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')

kidneyimaging_analysis('DCT', median = T, genes = tca_genes, 
                       gene_list_name = 'TCA', adjustment = 'epic_sglti2_1')
kidneyimaging_analysis('DCT', median = T, genes = ox_phos_genes, 
                       gene_list_name = 'Ox-Phos', adjustment = 'epic_sglti2_1')












