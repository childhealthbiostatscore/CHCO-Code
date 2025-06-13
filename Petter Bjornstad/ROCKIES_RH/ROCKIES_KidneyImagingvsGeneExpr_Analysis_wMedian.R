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


##a. PT Cells 
###i. Continuous 
#### TCA Cycle 
#### Unadjusted

#Filter to PT Cells
rm(so_celltype)
so_celltype <- subset(so_subset,celltype2=="PT")
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
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- subset(so_celltype,features=g)@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "Proximal Tubule Cells",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

######Adjusted

#Filter to PT Cells
rm(so_celltype)
so_celltype <- subset(so_subset,celltype2=="PT")
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
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- subset(so_celltype,features=g)@meta.data
      pred.formula <- as.formula(paste0("~",exposure,"+epic_sglti2_1"))
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
  colnames(full_results) <- c("logFC_Intercept","logFC",'logFC_SGLTi2', "se_Intercept",
                              "se",'se_SGLTi2', "p_Intercept","p_value",'p_value_SGLTi2', "gene_id","gene",
                              "Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_adjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "Proximal Tubule Cells, Adj. SGLT2i",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_PET_adjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




##### Ox Phos
######Unadjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,celltype2=="PT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

# With parallelization
#Ox Phos
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "Proximal Tubule Cells",
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

png(fs::path(dir.results, "Heatmap_ox_phos_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

######Adjusted

#Filter to PT Cells
rm(so_celltype)
so_celltype <- subset(so_subset,celltype2=="PT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

# With parallelization
#OX PHOS Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- subset(so_celltype,features=g)@meta.data
      pred.formula <- as.formula(paste0("~",exposure,"+epic_sglti2_1"))
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
  colnames(full_results) <- c("logFC_Intercept","logFC",'logFC_SGLTi2', "se_Intercept",
                              "se",'se_SGLTi2', "p_Intercept","p_value",'p_value_SGLTi2', "gene_id","gene",
                              "Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_OX_PHOS_cycle_PT_cells_PET_Variables_adjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_OX_PHOS_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "OX PHOS Cycle Genes vs. PET Variables (T2D)",
       subtitle = "Proximal Tubule Cells, Adj. SGLT2i",
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

png(fs::path(dir.results, "Heatmap_OX_PHOS_cycle_NEBULA_PT_PET_adjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()





###ii. Median Split
#### TCA Cycle 
#####Unadjusted

#Filter to PT Cells
rm(so_celltype)
so_celltype <- subset(so_subset,celltype2=="PT")
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
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()

for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept",
                              "se","p_Intercept","p_value","gene_id","gene",
                              "Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset_Median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "Proximal Tubule Cells (T2)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#####Adjusted

#Filter to PT Cells
rm(so_celltype)
so_celltype <- subset(so_subset,celltype2=="PT")
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
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()

for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure,"+epic_sglti2_1"))
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
  colnames(full_results) <- c("logFC_Intercept","logFC",'logFC_SGLTi2', "se_Intercept",
                              "se",'se_SGLTi2', "p_Intercept","p_value",'p_value_SGLTi2', "gene_id","gene",
                              "Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_adjusted_pooled_offset_Median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "Proximal Tubule Cells (T2D), Adj. SGLT2",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_PET_adjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#### Ox Phos
#####Unadjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,celltype2=="PT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

# With parallelization
#Ox Phos
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept",
                              "se", "p_Intercept","p_value","gene_id","gene",
                              "Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "Proximal Tubule Cells (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_NEBULA_PT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#####Adjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,celltype2=="PT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

# With parallelization
#Ox Phos
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure,"+epic_sglti2_1"))
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
  colnames(full_results) <- c("logFC_Intercept","logFC",'logFC_SGLTi2', "se_Intercept",
                              "se",'se_SGLTi2', "p_Intercept","p_value",'p_value_SGLTi2', "gene_id","gene",
                              "Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "Proximal Tubule Cells (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_NEBULA_PT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


### ii. PT-S1/2
### a. Continuous
##### TCA Cycle 
######Unadjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S1/S2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_S1_S2_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "PT-S1/S2",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes (Type 2 Diabetes)",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


##### Ox Phos
######Unadjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S1/S2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_S1_S2_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "PT-S1/S2",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

######Adjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S1/S2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure,"+epic_sglti2_1"))
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
  colnames(full_results) <- c("logFC_Intercept","logFC",'logFC_SGLTi2', "se_Intercept",
                              "se",'se_SGLTi2', "p_Intercept","p_value",'p_value_SGLTi2', "gene_id","gene",
                              "Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_S1_S2_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "PT-S1/S2",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




### b. Median Split
#### TCA Cycle 
#####Unadjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S1/S2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_S1_S2_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "PT-S1/S2 (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()
# 
# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes (Type 2 Diabetes)",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#####Adjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S1/S2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_S1_S2_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "PT-S1/S2 (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()
# 
# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes (Type 2 Diabetes)",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#### Ox Phos
#####Unadjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S1/S2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_S1_S2_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "PT-S1/S2 (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()
# 
# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#####Adjusted

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S1/S2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_S1_S2_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "PT-S1/S2 (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()
# 
# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_PT_S1_S2_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


### iii. PT-S3
####a. Continuous
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S3")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_S3_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "PT-S3",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_S3_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S3")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_S3_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "PT-S3",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_PT_S3_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


### b. Median Split
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S3")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_S3_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "PT-S3 (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_S3_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()
# 
# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes (Type 2 Diabetes)",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_PT_S3_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="PT-S3")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_S3_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "PT-S3 (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_PT_S3_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()
# 
# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_PT_S3_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


### iv. aPT
###a. Continuos 
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="aPT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_aPT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "aPT",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_aPT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="aPT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_aPT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "aPT",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_aPT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




### b. Median Split
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="aPT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_aPT_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "aPT (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_aPT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()
# 
# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes (Type 2 Diabetes)",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_aPT_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="aPT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_aPT_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "aPT (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_aPT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()
# 
# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_aPT_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


###i. TAL
####a. Continuous
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,TAL_celltype=="TAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_TAL_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "TAL",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes (Type 2 Diabetes)",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,TAL_celltype=="TAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_TAL_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "TAL",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




#### b. Median Split
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,TAL_celltype=="TAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_TAL_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "TAL (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed TCA Cycle Genes (Type 2 Diabetes)",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,TAL_celltype=="TAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_TAL_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "TAL (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


###ii. C-TAL-1
####a. Continuous
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="C-TAL-1")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_C_TAL_1_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "C-TAL-1",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_C_TAL_1_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="C-TAL-1")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_C_TAL_1_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "C-TAL-1",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_C_TAL_1_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



#### b. Median Split
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="C-TAL-1")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_C_TAL_1_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "C-TAL-1 (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_C_TAL_1_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="C-TAL-1")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_C_TAL_1_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "C-TAL-1 (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_C_TAL_1_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


###iii. C-TAL-2
####a. Continuous
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="C-TAL-2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_C_TAL_2_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "C-TAL-2",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_C_TAL_2_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="C-TAL-2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_C_TAL_2_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "C-TAL-2",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_C_TAL_2_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


####b. Median
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="C-TAL-2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_C_TAL_2_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "C-TAL-2 (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_C_TAL_2_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="C-TAL-2")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_C_TAL_2_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "C-TAL-2 (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_C_TAL_2_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


### iv. dTAL
####a. Continuous
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="dTAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_dTAL_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "dTAL",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_dTAL_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="dTAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_dTAL_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "dTAL",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_dTAL_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


####b. Median
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="dTAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_dTAL_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "dTAL (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_dTAL_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="dTAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_dTAL_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "dTAL (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_dTAL_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


### v. aTAL
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="aTAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
      if(is.null(converged)){
        converged <- NA
      }
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  tryCatch({
  nebula_summaries <- map_dfr(
    names(nebula_results_list),
    function(gene_name) {
      df <- nebula_results_list[[gene_name]]$summary
      df <- df %>% mutate(Gene = gene_name)
      return(df)
    }
  )
  })
  
  nonconverge_genes <- unique(PT_nebula_converged$Gene[which(PT_nebula_converged$Convergence_Code==-40)]) 
  
  #Make dataframe of final results
  full_results <- as.data.frame(nebula_summaries)
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_aTAL_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "aTAL",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_aTAL_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="aTAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_aTAL_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "aTAL",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_aTAL_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



####b. Median
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="aTAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_aTAL_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "aTAL (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_aTAL_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="aTAL")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_aTAL_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "aTAL (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_aTAL_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


###i. DCT
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,DCT_celltype=="DCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_DCT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "DCT",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_DCT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,DCT_celltype=="DCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#ox phos
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_DCT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "DCT",
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

png(fs::path(dir.results, "Heatmap_ox_phos_NEBULA_DCT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))

# Identify significant points (fdr < 0.05)
significant_df <- total_results[total_results$fdr < 0.05, ]

Genes <- length(unique(total_results$gene))
Nuclei <- ncol(so_celltype)
# Nonconvergence_Rate <- nebula_nonconverged_percent
# total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
#                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df3 <- total_results[total_results$fdr3 < 0.2, ]

# Your custom facet labels
# custom_labels <- c(
#   "Average Cortical K2", "Average Medulla K2",
#   "Average Cortical F", "Average Medulla F",
#   "Average Cortical K2/F", "Average Medulla K2/F"
# )
# 
# # Ensure 'Variable' is a factor in the correct order
# custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
total_results <- total_results %>%
  group_by(Variable) %>%
  mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
  ungroup()

# Plot
dot_plot <- ggplot(total_results, aes(
  # y = reorder(gene, logFC),
  y = gene_ordered,
  x = logFC,
  color = color,
  size = abs(logFC)
)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
  scale_color_identity() +
  scale_size(range = c(2, 6), name = "|LogFC|") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Differentially Expressed ox phos Genes (Type 2 Diabetes)",
    subtitle = "PET Variables in DCT Cells, Unadjusted (Pooled Offset)",
    x = "Log Fold Change",
    y = "Gene",
    caption = paste0(
      "FDR < 0.05, Genes = ", Genes,
      ", Nuclei = ", Nuclei
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.y = element_text(size = 8),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

# Print the plot
dot_plot


png(fs::path(dir.results, "DotPlot_ox_phos_NEBULA_DCT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 2500, height = 3000, res = 300)
print(dot_plot)
dev.off()


####b. Median
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,DCT_celltype=="DCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_DCT_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "DCT (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_DCT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="DCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_DCT_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "DCT (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_DCT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()


###ii. DCT1
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="DCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_DCT_1_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "DCT-1",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_DCT1_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()




#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="DCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#ox phos
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_DCT1_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "DCT-1",
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

png(fs::path(dir.results, "Heatmap_ox_phos_NEBULA_DCT1_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



####b. Median
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="DCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_DCT_1_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "DCT-1 (T2D)",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_DCT_1_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="DCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_DCT_1_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "DCT-1 (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_DCT_1_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()

###iii. dDCT
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="dDCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_dDCT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "dDCT",
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

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_dDCT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))

# Identify significant points (fdr < 0.05)
significant_df <- total_results[total_results$fdr < 0.05, ]

Genes <- length(unique(total_results$gene))
Nuclei <- ncol(so_celltype)
# Nonconvergence_Rate <- nebula_nonconverged_percent
# total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
#                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df3 <- total_results[total_results$fdr3 < 0.2, ]

# Your custom facet labels
# custom_labels <- c(
#   "Average Cortical K2", "Average Medulla K2",
#   "Average Cortical F", "Average Medulla F",
#   "Average Cortical K2/F", "Average Medulla K2/F"
# )
# 
# # Ensure 'Variable' is a factor in the correct order
# custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
total_results <- total_results %>%
  group_by(Variable) %>%
  mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
  ungroup()

# Plot
dot_plot <- ggplot(total_results, aes(
  # y = reorder(gene, logFC),
  y = gene_ordered,
  x = logFC,
  color = color,
  size = abs(logFC)
)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
  scale_color_identity() +
  scale_size(range = c(2, 6), name = "|LogFC|") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Differentially Expressed TCA Cycle Genes (Type 2 Diabetes)",
    subtitle = "PET Variables in dDCT Cells, Unadjusted (Pooled Offset)",
    x = "Log Fold Change",
    y = "Gene",
    caption = paste0(
      "FDR < 0.05, Genes = ", Genes,
      ", Nuclei = ", Nuclei
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.y = element_text(size = 8),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

# Print the plot
dot_plot


png(fs::path(dir.results, "DotPlot_TCA_cycle_NEBULA_dDCT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 2500, height = 3000, res = 300)
print(dot_plot)
dev.off()

#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="dDCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#ox phos
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2","avg_m_k2","avg_c_f","avg_m_f","avg_c_k2_f","avg_m_k2_f")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_dDCT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
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
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "dDCT",
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

png(fs::path(dir.results, "Heatmap_ox_phos_NEBULA_dDCT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))

# Identify significant points (fdr < 0.05)
significant_df <- total_results[total_results$fdr < 0.05, ]

Genes <- length(unique(total_results$gene))
Nuclei <- ncol(so_celltype)
# Nonconvergence_Rate <- nebula_nonconverged_percent
# total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
#                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df3 <- total_results[total_results$fdr3 < 0.2, ]

# Your custom facet labels
# custom_labels <- c(
#   "Average Cortical K2", "Average Medulla K2",
#   "Average Cortical F", "Average Medulla F",
#   "Average Cortical K2/F", "Average Medulla K2/F"
# )
# 
# # Ensure 'Variable' is a factor in the correct order
# custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
total_results <- total_results %>%
  group_by(Variable) %>%
  mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
  ungroup()

# Plot
dot_plot <- ggplot(total_results, aes(
  # y = reorder(gene, logFC),
  y = gene_ordered,
  x = logFC,
  color = color,
  size = abs(logFC)
)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
  scale_color_identity() +
  scale_size(range = c(2, 6), name = "|LogFC|") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Differentially Expressed ox phos Genes (Type 2 Diabetes)",
    subtitle = "PET Variables in dDCT Cells, Unadjusted (Pooled Offset)",
    x = "Log Fold Change",
    y = "Gene",
    caption = paste0(
      "FDR < 0.05, Genes = ", Genes,
      ", Nuclei = ", Nuclei
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.y = element_text(size = 8),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

# Print the plot
dot_plot


png(fs::path(dir.results, "DotPlot_ox_phos_NEBULA_dDCT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 2500, height = 3000, res = 300)
print(dot_plot)
dev.off()

####b. Median
##### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="dDCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- tca_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
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

write.csv(total_results,fs::path(dir.results,"NEBULA_TCA_cycle_dDCT_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "dDCT (T2D)",
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

#Try diff plot
heat_map_p <- ggplot(heatmap_data, aes(x = Variable, y = Gene, fill = logFC)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = signif), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#264653", mid = "white", high = "darkred",
    midpoint = 0,
    limits = c(-1, 1),  # Adjust these values to suit your data
    oob = scales::squish,  # Squish out-of-bounds values
    name = "LogFC"
  ) +
  theme_minimal() +
  labs(title = "TCA Cycle Genes vs. PET Variables (Median)",
       subtitle = "dDCT (T2D)",
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
heat_map_p
# custom_colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "darkred")

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_dDCT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()



##### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,KPMP_celltype=="dDCT")
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
k2_ids <- unique(so_celltype$kit_id[which(!is.na(so_celltype$avg_c_k2))])
so_celltype <- subset(so_celltype, kit_id %in% k2_ids)
DefaultAssay(so_celltype) <- "RNA" 

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

nrow(so_celltype) #34 genes
ncol(so_celltype) #4926 PT cells
# With parallelization
#TCA Cycle
# List of genes
genes_list <- ox_phos_genes

k2_vars <- c("avg_c_k2_med","avg_m_k2_med","avg_c_f_med","avg_m_f_med","avg_c_k2_f_med","avg_m_k2_f_med")
total_results <- data.frame()
for (exposure in k2_vars) {
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  start_time <- Sys.time()
  nebula_results_list <- foreach(g = genes_list, .packages = c("nebula", "Matrix")) %dopar% {
    tryCatch({
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene <- so_celltype@meta.data
      pred.formula <- as.formula(paste0("~",exposure))
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
  colnames(full_results) <- c("logFC_Intercept","logFC","se_Intercept","se","p_Intercept","p_value","gene_id","gene","Gene")
  full_results$Variable <- exposure
  #Calculate number of genes filtered out for low expression 
  full_results$low_exp <- length(ox_phos_genes)-length(full_results$gene)
  #Filter out non-converging genes
  full_results <- full_results %>% 
    filter(!gene %in%  nonconverge_genes)
  #Calculate nonconvergence rate
  full_results$nebula_nonconverged_percent <- paste0(round((1-(length(ox_phos_genes)-length(nonconverge_genes))/length(ox_phos_genes))*100,3),"%")
  # nebula_nonconverged_percent <- (length(rownames(counts_path))-length(unique(full_results$gene)))/length(rownames(counts_path))
  # print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
  full_results$fdr <- p.adjust(full_results[,6],method="fdr")  
  # mutate(fdr3=p.adjust(PValue3,method="fdr"))
  full_results$PValue10 <- -log10(pmax(full_results[,6], 1e-10))
  total_results <- rbind(total_results,full_results)
}

write.csv(total_results,fs::path(dir.results,"NEBULA_ox_phos_cycle_dDCT_cells_PET_Variables_unadjusted_pooled_offset_median.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
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
  labs(title = "Ox Phos Genes vs. PET Variables (Median)",
       subtitle = "dDCT (T2D)",
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

png(fs::path(dir.results, "Heatmap_ox_phos_cycle_NEBULA_dDCT_PET_unadjusted_pooled_offset_T2D_median.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()

# 
# total_results$color <- ifelse(total_results$fdr < 0.05 & total_results$`logFC` > 0, "darkred",
#                              ifelse(total_results$fdr < 0.05 & total_results$`logFC` < 0, "#264653", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df <- total_results[total_results$fdr < 0.05, ]
# 
# Genes <- length(unique(total_results$gene))
# Nuclei <- ncol(so_celltype)
# # Nonconvergence_Rate <- nebula_nonconverged_percent
# # total_results$color3 <- ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 > 0, "lightcoral",
# #                               ifelse(total_results$fdr3 < 0.2 & total_results$`logFC`3 < 0, "lightblue", "gray"))
# # 
# # # Identify significant points (fdr < 0.05)
# # significant_df3 <- total_results[total_results$fdr3 < 0.2, ]
# 
# # Your custom facet labels
# # custom_labels <- c(
# #   "Average Cortical K2", "Average Medulla K2",
# #   "Average Cortical F", "Average Medulla F",
# #   "Average Cortical K2/F", "Average Medulla K2/F"
# # )
# # 
# # # Ensure 'Variable' is a factor in the correct order
# # custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
# # total_results$Variable <- factor(total_results$Variable, levels = custom_order, labels = custom_labels)
# total_results <- total_results %>%
#   group_by(Variable) %>%
#   mutate(gene_ordered = fct_reorder2(gene, Variable, logFC)) %>%
#   ungroup()
# 
# # Plot
# dot_plot <- ggplot(total_results, aes(
#   # y = reorder(gene, logFC),
#   y = gene_ordered,
#   x = logFC,
#   color = color,
#   size = abs(logFC)
# )) +
#   geom_point(alpha = 0.7) +
#   facet_wrap(~Variable, scales = "free_x",ncol=2) +  # allow each facet its own x-axis range
#   scale_color_identity() +
#   scale_size(range = c(2, 6), name = "|LogFC|") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   labs(
#     title = "Differentially Expressed Ox Phos Genes among Youth with Type 2 Diabetes",
#     subtitle = "PET Variables in PT Cells, Unadjusted (Pooled Offset)",
#     x = "Log Fold Change",
#     y = "Gene",
#     caption = paste0(
#       "FDR < 0.05, Genes = ", Genes,
#       ", Nuclei = ", Nuclei
#     )
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0),
#     axis.text.y = element_text(size = 8),
#     axis.line = element_line(color = "black", size = 0.5),
#     axis.ticks.x = element_line(color = "black"),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # Print the plot
# dot_plot
# 
# 
# png(fs::path(dir.results, "DotPlot_ox_phos_cycle_NEBULA_TAL_PET_unadjusted_pooled_offset_T2D.png"), 
#     width = 2500, height = 3000, res = 300)
# print(dot_plot)
# dev.off()