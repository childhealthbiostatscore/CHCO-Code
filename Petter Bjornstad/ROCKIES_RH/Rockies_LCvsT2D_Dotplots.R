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


load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_LCvsT2D_data.RData')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

load('C:/Users/netio/Downloads/TCA_genes.txt')
load('C:/Users/netio/Downloads/OxPhos_genes.txt')


dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Increased_N/DotPlots/LC_T2D/'

#All cells 
so_celltype <- so_subset
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #13534 PT cells

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

write.csv(full_results,fs::path(dir.results,"NEBULA_TCA_cycle_ALL_cells_LC_T2D_NoMed_unadjusted_pooled_offset_no_IT_08.csv"))

full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                             ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))

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

dot_plot <- ggplot(full_results, aes(
  y =   gene,
  x = `logFC_groupType_2_Diabetes`,
  color = fdr,
  size = abs(`logFC_groupType_2_Diabetes`)
)) +
  geom_point(alpha = 0.7) +
   binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 
                palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                
                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                
                limits = c(0, 1),                 
                show.limits = TRUE, guide = 'colorsteps')  +
  scale_size(range = c(2, 6), name = "|LogFC|") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +  # Retains grid lines
  labs(
    title = "Differentially Expressed TCA Cycle Genes in All Cell Types",
    subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
dot_plot

png(fs::path(dir.results, "Plot_TCA_cycle_NEBULA_All_Cells_T2D_LC_NoMed_unadjusted_pooled_offset_no_IT_08.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()
 
### Ox Phos
   
#Filter to PT Cells
so_celltype <- so_subset
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #13534 PT cells

#Make sure exposure/independent/x variable or group variable is a factor variable
so_celltype$group <- factor(so_celltype$group)
#Make sure to set reference level
so_celltype$group  <- relevel(so_celltype$group ,ref="Lean_Control")

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

# With parallelization
#Ox Phos Cycle
# List of genes
genes_list <- ox_phos_genes

cl <- makeCluster(1)
registerDoParallel(cl)

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

write.csv(full_results,fs::path(dir.results,"NEBULA_Ox_Phos_cycle_ALL_cells_LC_T2D_NoMed_unadjusted_pooled_offset_no_IT_08.csv"))

full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                             ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))

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

dot_plot <- ggplot(full_results, aes(
  y =   gene,
  x = `logFC_groupType_2_Diabetes`,
  color = fdr,
  size = abs(`logFC_groupType_2_Diabetes`)
)) +
  geom_point(alpha = 0.7) +
   binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
  scale_size(range = c(2, 6), name = "|LogFC|") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +  # Retains grid lines
  labs(
    title = "Differentially Expressed Ox Phos Cycle Genes in All Cell Types",
    subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
dot_plot

png(fs::path(dir.results, "Plot_Ox_Phos_cycle_NEBULA_All_Cells_T2D_LC_NoMed_unadjusted_pooled_offset_no_IT_08.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()
 

### i. All PT Cells 
#### TCA Cycle 
   
#Filter to PT Cells
PT_cells <- c('PT-1', 'PT-2', 'PT-3', 'PT-4', 'PT-5', 'aPT')
so_celltype <- subset(so_subset,celltype_harmony ==PT_cells)
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #13534 PT cells

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

write.csv(full_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_LC_T2D_NoMed_unadjusted_pooled_offset_no_IT_08.csv"))

full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                             ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))

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

dot_plot <- ggplot(full_results, aes(
  y =   gene,
  x = `logFC_groupType_2_Diabetes`,
  color = fdr,
  size = abs(`logFC_groupType_2_Diabetes`)
)) +
  geom_point(alpha = 0.7) +
   binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
  scale_size(range = c(2, 6), name = "|LogFC|") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +  # Retains grid lines
  labs(
    title = "Differentially Expressed TCA Cycle Genes in PT Cells",
    subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
dot_plot

png(fs::path(dir.results, "Plot_TCA_cycle_NEBULA_PT_Cells_T2D_LC_NoMed_unadjusted_pooled_offset_no_IT_08.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()
 

#### Ox Phos
   
#Filter to PT Cells
PT_cells <- c('PT-1', 'PT-2', 'PT-3', 'PT-4', 'PT-5', 'aPT')
so_celltype <- subset(so_subset,celltype_harmony ==PT_cells)
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #13534 PT cells

#Make sure exposure/independent/x variable or group variable is a factor variable
so_celltype$group <- factor(so_celltype$group)
#Make sure to set reference level
so_celltype$group  <- relevel(so_celltype$group ,ref="Lean_Control")


counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

# With parallelization
#Ox_Phos Cycle
# List of genes
genes_list <- ox_phos_genes

cl <- makeCluster(1)
registerDoParallel(cl)

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

write.csv(full_results,fs::path(dir.results,"NEBULA_OX_PHOS_cycle_PT_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv"))

full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                             ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))

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

dot_plot <- ggplot(full_results, aes(
  y =   gene,
  x = `logFC_groupType_2_Diabetes`,
  color = fdr,
  size = abs(`logFC_groupType_2_Diabetes`)
)) +
  geom_point(alpha = 0.7) +
   binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
  scale_size(range = c(2, 6), name = "|LogFC|") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +  # Retains grid lines
  labs(
    title = "Differentially Expressed OX PHOS Cycle Genes in PT Cells",
    subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
dot_plot

png(fs::path(dir.results, "Plot_OX_PHOS_cycle_NEBULA_PT_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()
 


### ii. PT Subtypes
   
celltypes <- c("PT-S1/S2","PT-S3","aPT")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
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
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
}
 

###i. TAL
   
celltypes <- c("TAL")
TAL_cells <- c('C-TAL-1', 'C-TAL-2', 'aTAL', 'dTAL')
so_celltype <- subset(so_subset,KPMP_celltype == TAL_cells)
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
 

### ii. TAL Subtypes
   
celltypes <- c("C-TAL-1","C-TAL-2","aTAL","dTAL")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
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
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
}
 


###i. DCT
   
celltypes <- c("DCT", 'dDCT')
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype == celltype)
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
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
 
### ii. DCT Subtypes
   
celltypes <- c("DCT","dDCT")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
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
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
}
 

### i. EC Cells
   
celltypes <- c("EC/VSMC","EC-AVR","EC-PTC","EC-AEA","EC-LYM","EC-GC")
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype ==celltype)
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
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
 
### ii. EC Subtypes
   
celltypes <- c("EC/VSMC","EC-AVR","EC-PTC","EC-AEA","EC-LYM","EC-GC")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
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
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
}
 

### i. Podocytes
   
celltypes <- c("POD")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,celltype2==celltype)
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
  
  # With parallelization
  #TCA Cycle
  # List of genes
  genes_list <- tca_genes
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
}
 

#Immune cells
   
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
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
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed TCA Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
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
  
  full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` > 0, "lightcoral",
                               ifelse(full_results$fdr < 0.05 & full_results$`logFC_groupType_2_Diabetes` < 0, "lightblue", "gray"))
  
  # Identify significant points (fdr < 0.05)
  significant_df <- full_results[full_results$fdr < 0.05, ]
  
  Genes <- length(unique(full_results$gene))
  Cell <- ncol(so_celltype)
  Nonconvergence_Rate <- nebula_nonconverged_percent
  
  
  max <- max(full_results$`logFC_groupType_2_Diabetes`)
  # max <- 3.1
  min <- min(full_results$`logFC_groupType_2_Diabetes`)
  
  dot_plot_ox_phos <- ggplot(full_results, aes(
    y =   gene,
    x = `logFC_groupType_2_Diabetes`,
    color = fdr,
    size = abs(`logFC_groupType_2_Diabetes`)
  )) +
    geom_point(alpha = 0.7) +
     binned_scale(aesthetics = 'color', scale_name = 'stepsn',                 palette = function(x) c('red', 'orange', 'purple', 'blue', 'grey'),                breaks = c(0.05, 0.1, 0.15, 0.2, 1),                limits = c(0, 1),                 show.limits = TRUE, guide = 'colorsteps')  +
    scale_size(range = c(2, 6), name = "|LogFC|") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +  # Retains grid lines
    labs(
      title = paste0("Differentially Expressed OX PHOS Cycle Genes in ",celltype," Cells"),
      subtitle = "T2D vs. LC (No Medication), Unadjusted (Pooled Offset)",
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
  png(fs::path(dir.results, paste("Plot_NEBULA_",celltype2,"_Cells_T2D_LC_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(comb_plot)
  dev.off()
  
}