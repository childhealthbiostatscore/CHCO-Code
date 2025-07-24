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









load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')

dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/'

load("C:/Users/netio/Downloads/TCA_genes.txt")
load('C:/Users/netio/Downloads/OxPhos_genes.txt')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

test <- so_subset@meta.data

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(screen_date) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat <- dat %>% dplyr::select(record_id, epic_sglti2_1) %>% 
  filter(!is.na(epic_sglti2_1))

test$epic_sglti2_1 <- NULL

test <- test %>% left_join(dat, by='record_id')

so_subset@meta.data$epic_sglti2_1 <- test$epic_sglti2_1



#Filter to PT Cells
so_celltype <- so_subset
# so_celltype <- subset(so_celltype,group=="Type_2_Diabetes")
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #68909 PT cells

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

test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
test <- test[-which(duplicated(test)),]


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

write.csv(full_results,fs::path(dir.results,"NEBULA_TCA_cycle_ALL_cells_SGLT2_T2D_Only_unadjusted_pooled_offset.csv"))

full_results$color1 <- ifelse(full_results$fdr < 0.05, "lightcoral", "gray")
full_results$color2 <- ifelse(full_results$pvalue < 0.05, "lightcoral", "gray")

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

order_tca <- data.table::fread("C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/NEBULA_TCA_cycle_ALL_cells_LC_T2D_NoMed_unadjusted_pooled_offset_no_IT_08.csv") %>%
  arrange(`logFC_groupType_2_Diabetes`)

dot_plot <- ggplot(full_results, aes(
  y = factor(gene, levels = order_tca$x),
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
    title = "Differentially Expressed TCA Cycle Genes in All Cell Types",
    subtitle = "SGLT2i vs. No SGLT2i (T2D Only), Unadjusted (Pooled Offset)",
    x = "Log Fold Change",
    y = "Gene",
    caption = paste0("Participants on SGLT2: ", sglt2_count, ', No SGLT2: ', nosglt2_count, 
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

png(fs::path(dir.results, "Plot_TCA_cycle_NEBULA_All_Cells_SGLT2_T2D_unadjusted_pooled_offset.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()

#### Ox Phos

#Filter to PT Cells
so_celltype <- so_subset
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #68909 PT cells

#Make sure exposure/independent/x variable or group variable is a factor variable
so_celltype$epic_sglti2_1 <- factor(so_celltype$epic_sglti2_1)
#Make sure to set reference level
so_celltype$epic_sglti2_1  <- relevel(so_celltype$epic_sglti2_1 ,ref="No")

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

# With parallelization
#OX_PHOS Cycle
# List of genes
genes_list <- ox_phos_genes

cl <- makeCluster(1)
registerDoParallel(cl)

test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
test <- test[-which(duplicated(test)),]


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

write.csv(full_results,fs::path(dir.results,"NEBULA_OX_PHOS_cycle_ALL_cells_SGLT2_T2D_Only_unadjusted_pooled_offset.csv"))

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
    title = "Differentially Expressed OX PHOS Cycle Genes in All Cell Types",
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
dot_plot

png(fs::path(dir.results, "Plot_OX_PHOS_cycle_NEBULA_All_Cells_SGLT2_T2D_unadjusted_pooled_offset.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()

##b. PT Cells 
#### TCA Cycle 

#Filter to PT Cells
so_celltype <- subset(so_subset,celltype2=="PT")
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #13534 PT cells

#Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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

test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
test <- test[-which(duplicated(test)),]


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

write.csv(full_results,fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_SGLT2_T2D_Only_unadjusted_pooled_offset.csv"))

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
    title = "Differentially Expressed TCA Cycle Genes in PT Cells",
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
dot_plot

png(fs::path(dir.results, "Plot_TCA_cycle_NEBULA_PT_Cells_SGLT2_T2D_unadjusted_pooled_offset.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()

#### Ox Phos

#Filter to PT Cells
so_celltype <- subset(so_subset,celltype2=="PT")
DefaultAssay(so_celltype) <- "RNA" 

nrow(so_celltype) #34 genes
ncol(so_celltype) #13534 PT cells

#Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
so_celltype$epic_sglti2_1 <- factor(so_celltype$epic_sglti2_1)
#Make sure to set reference level
so_celltype$epic_sglti2_1  <- relevel(so_celltype$epic_sglti2_1 ,ref="No")


counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round

# With parallelization
#Ox_Phos Cycle
# List of genes
genes_list <- ox_phos_genes

cl <- makeCluster(1)
registerDoParallel(cl)

test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
test <- test[-which(duplicated(test)),]


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

write.csv(full_results,fs::path(dir.results,"NEBULA_OX_PHOS_cycle_PT_cells_SGLT2_T2D_Only_unadjusted_pooled_offset.csv"))

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
    title = "Differentially Expressed OX PHOS Cycle Genes in PT Cells",
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
dot_plot

png(fs::path(dir.results, "Plot_OX_PHOS_cycle_NEBULA_PT_Cells_SGLT2_T2D_unadjusted_pooled_offset.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()

###i. PT Subtypes

celltypes <- c("PT-S1/S2","PT-S3","aPT")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
    y = gene,
    x = `logFC_epic_sglti2_1Yes`,
    color = fdr,
    size = abs(`logFC_epic_sglti2_1Yes`)
  )) +
    geom_point(alpha = 0.7) +
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
  
}


##c. TAL

celltypes <- c("TAL")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,celltype2==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
  
}

###i. TAL Subtypes

celltypes <- c("C-TAL-1","C-TAL-2","aTAL","dTAL")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
}


##d. DCT

celltypes <- c("DCT")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,DCT_celltype==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
  
}

### i. DCT Subtypes

celltypes <- c("DCT","dDCT")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
  
}


##e. EC Cells

celltypes <- c("EC")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,celltype2==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
}

###i. EC Subtypes

celltypes <- c("EC/VSMC","EC-AVR","EC-PTC","EC-AEA","EC-LYM","EC-GC")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,KPMP_celltype==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
  
}


##f. Podocytes

celltypes <- c("POD")
for (celltype in celltypes) {
  #Filter to PT Cells
  so_celltype <- subset(so_subset,celltype2==celltype)
  DefaultAssay(so_celltype) <- "RNA" 
  
  nrow(so_celltype) #34 genes
  ncol(so_celltype) #13534 PT cells
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  #Make sure exposure/independent/x variable or SGLT2 variable is a factor variable
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
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
}


##g. Immune cells

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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  
  cl <- makeCluster(1)
  registerDoParallel(cl)
  
  test <- so_subset@meta.data %>% dplyr::select(record_id, group, epic_sglti2_1)
  test <- test[-which(duplicated(test)),]
  
  
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
  png(fs::path(dir.results, paste("fdr/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_TCA_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_tca_2)
  dev.off()
  
  png(fs::path(dir.results, paste("fdr/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos)
  dev.off()
  
  png(fs::path(dir.results, paste("pvalue/Plot_OxPhos_NEBULA_",celltype2,"_Cells_T2D_SGLT2_NoMed_unadjusted_pooled_offset.png")),
      width = 5000, height = 2000, res = 300)
  print(dot_plot_ox_phos_2)
  dev.off()
  
  
  
  
}