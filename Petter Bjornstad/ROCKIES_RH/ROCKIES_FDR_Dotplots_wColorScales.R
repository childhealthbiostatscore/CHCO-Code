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


test <- harmonized_data %>% filter(kit_id %in% small_meta.data_unique$kit_id)
test2 <- harmonized_data %>% filter(mrn %in% test$mrn)

test3 <- harmonized_data %>% filter(mrn %in% test2$mrn) %>% 
  dplyr::select(record_id, kit_id, mrn, lc_k2, rc_k2, lc_f, rc_f, 
                lm_k2, rm_k2, lm_f, rm_f) %>% group_by(mrn) %>% 
  summarize(lc_k2 = mean(lc_k2, na.rm=T), rc_k2 = mean(rc_k2, na.rm=T), 
            lm_k2 = mean(lm_k2, na.rm=T), rm_k2 = mean(rm_k2, na.rm=T),
            lc_f = mean(lc_f, na.rm=T), rc_f = mean(rc_f, na.rm=T),
            lm_f = mean(lm_f, na.rm=T), rm_f = mean(rm_f, na.rm=T))

write.table(test3, 'C:/Users/netio/Documents/UofW/Rockies/T2D_rawimagingdata.txt', row.names=F, quote=F, sep='\t')

fixed_data <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/T2D_rawimagingdata.txt')

# #Calculate K2 and F variables
fixed_data <- fixed_data %>%
  mutate(avg_c_k2 = (lc_k2+rc_k2)/2) %>%
  mutate(avg_m_k2 = (lm_k2+rm_k2)/2) %>%
  mutate(avg_c_f = (lc_f+rc_f)/2) %>%
  mutate(avg_m_f = (lm_f+rm_f)/2)
fixed_data <- fixed_data %>%
  rowwise() %>%
  mutate(avg_c_k2_f = (avg_c_k2/avg_c_f)) %>%
  mutate(avg_m_k2_f = (avg_m_k2/avg_m_f))

med_c_k2 <- median(fixed_data$avg_c_k2, na.rm=T)
med_m_k2 <- median(fixed_data$avg_m_k2, na.rm=T)
med_c_f <- median(fixed_data$avg_c_f, na.rm=T)
med_m_f <- median(fixed_data$avg_m_f, na.rm=T)
med_c_k2_f <- median(fixed_data$avg_c_k2_f, na.rm=T)
med_m_k2_f <- median(fixed_data$avg_m_k2_f, na.rm=T)


fixed_data <- fixed_data %>% 
  rowwise() %>% 
  mutate(avg_c_k2_med = ifelse(avg_c_k2 >= med_c_k2, 'Above Median', 'Below Median'),
         avg_m_k2_med = ifelse(avg_m_k2 >= med_m_k2, 'Above Median', 'Below Median'),
         avg_c_f_med = ifelse(avg_c_f >= med_c_f, 'Above Median', 'Below Median'),
         avg_m_f_med = ifelse(avg_m_f >= med_m_f, 'Above Median', 'Below Median'),
         avg_c_k2_f_med = ifelse(avg_c_k2_f >= med_c_k2_f, 'Above Median', 'Below Median'),
         avg_m_k2_f_med = ifelse(avg_m_k2_f >= med_m_k2_f, 'Above Median', 'Below Median'))



#Try to add this data from dat 
so_subset@meta.data <- so_subset@meta.data[, !colnames(so_subset@meta.data) %in% c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 
                                                                                   'avg_m_f', 
                                                                                   'avg_c_k2_f', 'avg_m_k2_f',
                                                                                   'lc_k2', 'rc_k2', 'lm_k2', 
                                                                                   'rm_k2', 'lm_f', 'rm_f',
                                                                                   'avg_c_k2_med', 'avg_m_k2_med',
                                                                                   'avg_c_f_med', 'avg_m_f_med', 
                                                                                   'avg_c_k2_f_med', 'avg_m_k2_f_med')]

so_subset@meta.data <- so_subset@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(fixed_data, by = "mrn") %>%
  tibble::column_to_rownames("cell_id")


so_subset$avg_c_k2_med <- factor(so_subset$avg_c_k2_med)
so_subset$avg_c_k2_med <- relevel(so_subset$avg_c_k2_med,"Below Median")

so_subset$avg_m_k2_med <- factor(so_subset$avg_m_k2_med)
so_subset$avg_m_k2_med <- relevel(so_subset$avg_m_k2_med,"Below Median")

so_subset$avg_c_f_med <- factor(so_subset$avg_c_f_med)
so_subset$avg_c_f_med <- relevel(so_subset$avg_c_f_med,"Below Median")

so_subset$avg_m_f_med <- factor(so_subset$avg_m_f_med)
so_subset$avg_m_f_med <- relevel(so_subset$avg_m_f_med,"Below Median")

so_subset$avg_c_k2_f_med <- factor(so_subset$avg_c_k2_f_med)
so_subset$avg_c_k2_f_med <- relevel(so_subset$avg_c_k2_f_med,"Below Median")

so_subset$avg_m_k2_f_med <- factor(so_subset$avg_m_k2_f_med)
so_subset$avg_m_k2_f_med <- relevel(so_subset$avg_m_k2_f_med,"Below Median")














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

