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







#For T2D Only 

load('C:/Users/netio/Documents/UofW/Rockies/Line4875_Rockies.RData')

#For Both 



#Function
kidneyimaging_analysis <- function(celltype, genes, gene_list_name = 'TCA', T2D_Only = TRUE, 
                                   dir.results, cl_number = 1, cpc = 0.005, 
                                   set_cutoff = F, logFC_thresh = 10, pvalue_thresh = 0.1){
  
  if(celltype == 'All'){
    so_celltype <- sO_subset
  }else if(celltype == 'PT'){
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
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #4926 PT cells
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  
  celltype <- str_replace_all(celltype, pattern='/', replacement = '_')
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- genes
  
  total_results <- data.frame()

  nebula_results_list <- list()
  for(i in c(1:length(genes_list))){
    tryCatch({
      g <- genes_list[i]
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- subset(so_celltype,features=g)@meta.data
      
      if(T2D_Only == TRUE){
        tmp.formula <- as.formula(paste0('~epic_sglti2_1'))
      }else{
        tmp.formula <- as.formula(paste0('~group'))
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
      result <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, 
                       ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,
                       offset=data_g_gene$library, cpc= cpc)
      
      nebula_results_list[[i]] <- result
#      names(nebula_results_list[[i]]) <- g
    })
      # return both gene name and result
      
#    }, error = function(e) {
#      list(gene = g, summary = NA, overdispersion = NA, convergence = NA, 
#           algorithm = NA, covariance = NA, random_effect = NA)
#    })
#  }
#    nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
#      tryCatch({
#        count_gene <- counts_path[g, , drop = FALSE]
#        meta_gene <- subset(so_celltype,features=g)@meta.data
#        
#        if(T2D_Only == TRUE){
#          tmp.formula <- as.formula(paste0('~epic_sglti2_1'))
#        }else{
#          tmp.formula <- as.formula(paste0('~group'))
#        }
#        
#        pred.formula <- as.formula(tmp.formula)
#        pred_gene <- model.matrix(pred.formula, data = meta_gene)
#        # library <- meta_gene$library_size
#        library <- meta_gene$pooled_offset
#        data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
#        
#        if (is.null(data_g_gene)) {
#          data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
#        }
#        
#        #With offset
#        result <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, 
#                         ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,
#                         offset=data_g_gene$library, cpc= cpc)
#        
#        list(gene = g, result = result)  # return both gene name and result
#        
#      }, error = function(e) {
#        list(gene = g, summary = NA, overdispersion = NA, convergence = NA, 
#             algorithm = NA, covariance = NA, random_effect = NA)
#      })
#    }
    
#    stopCluster(cl)
#    end_time <- Sys.time()
#    print(end_time - start_time)
    print(paste0('Done with ', g, ' Gene'))
}
    
    # set the names of results based on gene names
    nebula_results_list <- Filter(Negate(is.null), nebula_results_list)  # remove NULLs first
    names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)  # set names
    nebula_results_list <- lapply(nebula_results_list, function(x) x$result)  # clean list back to just results
    
    PT_nebula_converged <- map_dfr(
      names(nebula_results_list),
      function(gene_name) {
        # Safely extract convergence code
        converged <- tryCatch({
          conv <- nebula_results_list[[gene_name]]$convergence
          if (is.null(conv) || length(conv) == 0) NA else conv
        }, error = function(e) NA)
        
        data.frame(Gene = gene_name,
                   Convergence_Code = converged)
      }
    )
    
    nebula_summaries <- map_dfr(
      names(nebula_results_list),
      function(gene_name) {
        tryCatch({
          # Check if the result exists and has summary info
          result <- nebula_results_list[[gene_name]]
          
          if (is.null(result) || is.null(result$summary)) {
            # If no summary info, return NULL (will be filtered out by map_dfr)
            return(NULL)
          } else {
            df <- result$summary
            
            # Check if summary is empty or not a data.frame
            if (is.null(df) || nrow(df) == 0 || !is.data.frame(df)) {
              return(NULL)
            } else {
              df <- df %>% mutate(Gene = gene_name)
              return(df)
            }
          }
        }, error = function(e) {
          # If any error occurs, return NULL (will be filtered out)
          cat("Error processing gene", gene_name, ":", e$message, "\n")
          return(NULL)
        })
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
    
    write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TCA_cycle_",celltype2,"_cells_SGLT2_T2D_Only_unadjusted_pooled_offset.csv")))

    total_results <- rbind(total_results,full_results)

    file.name <- paste0(dir.results, 'NEBULA_', gene_list_name, 
                        '_', celltype, '_cells_PET_adjusted_pooled_offset_T2D_SGLT2i_dotplotdata.txt')

  
  
  write.table(total_results,file.name, row.names=F, quote=F, sep='\t')
  
  # total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
  
  # Define significance stars
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_epic_sglti2_1Yes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_epic_sglti2_1Yes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  # full_results$color3 <- ifelse(full_results$fdr3 < 0.2 & full_results$`logFC_epic_sglti2_1Yes`3 > 0, "lightcoral",
  #                               ifelse(full_results$fdr3 < 0.2 & full_results$`logFC_epic_sglti2_1Yes`3 < 0, "lightblue", "gray"))
  # 
  # # Identify significant points (fdr < 0.05)
  # significant_df3 <- full_results[full_results$fdr3 < 0.2, ]
  
  max <- max(full_results$`logFC_epic_sglti2_1Yes`)
  # max <- 3.1
  min <- min(full_results$`logFC_epic_sglti2_1Yes`)
  
  dot_plot <- ggplot(full_results, aes(
    y = gene,
    x = `logFC_epic_sglti2_1Yes`,
    color = fdr,
    size = abs(`logFC_epic_sglti2_1Yes`)
  )) +
    geom_point(alpha = 0.7) +
    binned_scale(aesthetics = 'color', scale_name = 'stepsn', 
                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),
                 breaks = c(0.05, 0.1, 0.15, 0.2, 1),
                 limits = c(0, 1), 
                 show.limits = TRUE, guide = 'colorsteps') +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed ", gene_list_name, " Cycle Genes in ",celltype," Cells"),
      subtitle = "SGLT2i vs. No SGLT2i (T2D Only), Pooled Offset",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0(
        "FDR < 0.05, Genes = ", Genes,
        ", Cells = ", Cell,
        ", Non-Convergence Rate: ", Nonconvergence_Rate,
        ", Genes Filtered out for Low Expression: ", low_exp
      )
    ) +
    theme(
      plot.title = element_text(hjust = 0),
      axis.text.y = element_text(size = 8),
      # axis.text.x = element_text(angle = 0, hjust = 1),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks.x = element_line(color = "black"),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  
  
  
  
  # Select only the needed columns and rename LogFC for clarity
    subtitle1 <- paste0(celltype, ', ')
    subtitle2 <- paste0('with no adjustment')
  subtitle <- paste0(subtitle1, subtitle2)
  
  

  # custom_colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "darkred")
    file.name <- paste0(dir.results, 'NEBULA_', gene_list_name, 
                        '_', celltype, '_pooled_offset_T2D_SGLT2i_dotplot.png')  

  
  
  # Ensure clean graphics device state
  while(dev.cur() > 1) dev.off()
  
  # Remove existing file if it exists
  if (file.exists(file.name)) {
    file.remove(file.name)
  }
  
  # Create PNG with error handling
  tryCatch({
    png(file.name, 
        width = 1500, height = 2000, res = 300)
    print()
  }, finally = {
    dev.off()
  })
  
  # Verify file was created
  if(file.exists(file.name)) {
    cat("PNG file successfully created:", file.name, "\n")
  } else {
    cat("Failed to create PNG file:", file.name, "\n")
  }
  
} 
    
    
    
    kidneyimaging_analysis('PT', genes = tca_genes, T2D_Only = TRUE,
                           gene_list_name = 'TCA', dir.results = 'C:/Users/netio/Documents/UofW/Rockies/Increased_N/')




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #REFERENCE DATA 
    



celltypes <- c(
  "cDC",
  "cycT",
  "CD4+ T",
  "CD8+ T",
  "NK",
  "B",
  "MON",
  "MAC",
  "MC")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  celltype2 <- str_replace_all(celltype2," ","")
  celltype2 <- str_replace_all(celltype2,"\\+","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
  so_celltype$epic_sglti2_1 <- factor(so_celltype$epic_sglti2_1)
  #Make sure to set reference level
  so_celltype$epic_sglti2_1  <- relevel(so_celltype$epic_sglti2_1 ,ref="No")
  
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
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
  
  write.csv(full_results,fs::path(dir.results,paste0("NEBULA_TCA_cycle_",celltype2,"_cells_SGLT2_T2D_Only_unadjusted_pooled_offset.csv")))
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_epic_sglti2_1Yes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_epic_sglti2_1Yes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  # full_results$color3 <- ifelse(full_results$fdr3 < 0.2 & full_results$`logFC_epic_sglti2_1Yes`3 > 0, "lightcoral",
  #                               ifelse(full_results$fdr3 < 0.2 & full_results$`logFC_epic_sglti2_1Yes`3 < 0, "lightblue", "gray"))
  # 
  # # Identify significant points (fdr < 0.05)
  # significant_df3 <- full_results[full_results$fdr3 < 0.2, ]
  
  max <- max(full_results$`logFC_epic_sglti2_1Yes`)
  # max <- 3.1
  min <- min(full_results$`logFC_epic_sglti2_1Yes`)
  
  dot_plot_tca <- ggplot(full_results, aes(
    y = reorder(gene, `logFC_epic_sglti2_1Yes`),
    x = `logFC_epic_sglti2_1Yes`,
    color = fdr,
    size = abs(`logFC_epic_sglti2_1Yes`)
  )) +
    geom_point(alpha = 0.7) +
    binned_scale(aesthetics = 'color', scale_name = 'stepsn', 
                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),
                 breaks = c(0.05, 0.1, 0.15, 0.2, 1),
                 limits = c(0, 1), 
                 show.limits = TRUE, guide = 'colorsteps') +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "SGLT2i vs. No SGLT2i (T2D Only), Unadjusted (Pooled Offset)",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0(
        "FDR < 0.05, Genes = ", Genes,
        ", Cells = ", Cell,
        ", Non-Convergence Rate: ", Nonconvergence_Rate,
        ", Genes Filtered out for Low Expression: ", low_exp
      )
    ) +
    theme(
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
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
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
  
  write.csv(full_results,fs::path(dir.results,paste0("NEBULA_OX_PHOS_cycle_",celltype2,"_cells_SGLT2_T2D_Only_unadjusted_pooled_offset.csv")))
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_epic_sglti2_1Yes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_epic_sglti2_1Yes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_epic_sglti2_1Yes`)
  # max <- 3.1
  min <- min(full_results$`logFC_epic_sglti2_1Yes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y = reorder(gene, `logFC_epic_sglti2_1Yes`),
    x = `logFC_epic_sglti2_1Yes`,
    color = fdr,
    size = abs(`logFC_epic_sglti2_1Yes`)
  )) +
    geom_point(alpha = 0.7) +
    binned_scale(aesthetics = 'color', scale_name = 'stepsn', 
                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),
                 breaks = c(0.05, 0.1, 0.15, 0.2, 1),
                 limits = c(0, 1), 
                 show.limits = TRUE, guide = 'colorsteps') +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "SGLT2i vs. No SGLT2i (T2D Only), Unadjusted (Pooled Offset)",
      x = "Log Fold Change",
      y = "Gene",
      caption = paste0(
        "FDR < 0.05, Genes = ", Genes,
        ", Cells = ", Cell,
        ", Non-Convergence Rate: ", Nonconvergence_Rate,
        ", Genes Filtered out for Low Expression: ", low_exp
      )
    ) +
    theme(
      plot.title = element_text(hjust = 0),
      axis.text.y = element_text(size = 8),
      # axis.text.x = element_text(angle = 0, hjust = 1),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks.x = element_line(color = "black"),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  
  comb_plot <- dot_plot_tca + dot_plot_ox_phos
  # dot_plot
  # 
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_SGLT2_T2D_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
}

