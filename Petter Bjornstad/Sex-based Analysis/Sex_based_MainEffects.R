#Load in packages for functions we need. 

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
library(GSEABase)
library(clusterProfiler)
library('org.Hs.eg.db')


load('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line265.RData')



results.dir <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/MainEffects_Sex/'
#LC Sex Analyses

so <- subset(so, group == 'Type_2_Diabetes')

#All Analysis
so_subset <- so

cellcounts <- ncol(so_subset)
partcounts <- length(unique(so_subset@meta.data$kit_id))

celltype2 <- 'All'
full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(full_analysis)
full_analysis <- subset(so_subset, features = hvgs)
full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 


meta_gene <- full_analysis@meta.data
pred_gene <- model.matrix(~sex, data = meta_gene)
data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                          offset = full_analysis$pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary

write.table(result_allcells, paste(results.dir, 'NEBULA_LC_SexDifferences_AllCells.txt'), row.names=F,
            quote=F, sep='\t')


#Plotting 
sig_results <- as.data.frame(result_allcells) %>% 
  dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)

tmp_df <- sig_results
tmp_df$diffexp <- 'No'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'

tmp_df <- tmp_df %>% arrange(Pvalue)
tmp_df$label <- NA
tmp_df$label[1:20] <- tmp_df$Gene[1:20]

tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)

if(length(unique(tmp_df$diffexp)) > 1){
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('orange', 'grey', 'purple'),
                       labels = c('Downregulated', 'Not significant', 'Upregulated'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                         title = paste0('Lean Controls Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
}else{
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('grey'),
                       labels = c('Not significant'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                         title = paste0('Lean Controls Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  
}
pdf(paste0(results.dir, 'NEBULA_LC_SexDifferences_Volcano_AllCells.pdf'))
print(tmp_graph)
dev.off()





#PT Cells
  so_subset <- subset(so,celltype2 == 'PT')
  
  cellcounts <- ncol(so_subset)
  partcounts <- length(unique(so_subset@meta.data$kit_id))
  
  celltype2 <- 'PT'
  full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(full_analysis)
  full_analysis <- subset(so_subset, features = hvgs)
  full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 
  
  
  meta_gene <- full_analysis@meta.data
  pred_gene <- model.matrix(~sex, data = meta_gene)
  data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
  result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                            offset = full_analysis$pooled_offset,
                            ncore = 1, output_re = T, covariance = T,
                            reml = T, model = "NBLMM")
  
  result_allcells$summary
  
  write.table(result_allcells, paste(results.dir, 'NEBULA_LC_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
              quote=F, sep='\t')
  
  
  #Plotting 
  sig_results <- as.data.frame(result_allcells) %>% 
    dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                  Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)
  
  tmp_df <- sig_results
  tmp_df$diffexp <- 'No'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
  
  tmp_df <- tmp_df %>% arrange(Pvalue)
  tmp_df$label <- NA
  tmp_df$label[1:20] <- tmp_df$Gene[1:20]
  
  tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)
  
  if(length(unique(tmp_df$diffexp)) > 1){
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                           title = paste0('Lean Controls Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  }else{
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                           title = paste0('Lean Controls Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
    
  }
  pdf(paste0(results.dir, 'NEBULA_LC_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
  print(tmp_graph)
  dev.off()
  
  
  
  
  
  
  #EC Cells
  so_subset <- subset(so,celltype2 == 'EC')
  
  cellcounts <- ncol(so_subset)
  partcounts <- length(unique(so_subset@meta.data$kit_id))
  
  celltype2 <- 'EC'
  full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(full_analysis)
  full_analysis <- subset(so_subset, features = hvgs)
  full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 
  
  
  meta_gene <- full_analysis@meta.data
  pred_gene <- model.matrix(~sex, data = meta_gene)
  data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
  result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                            offset = full_analysis$pooled_offset,
                            ncore = 1, output_re = T, covariance = T,
                            reml = T, model = "NBLMM")
  
  result_allcells$summary
  
  write.table(result_allcells, paste(results.dir, 'NEBULA_T2D_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
              quote=F, sep='\t')
  
  
  #Plotting 
  sig_results <- as.data.frame(result_allcells) %>% 
    dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                  Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)
  
  tmp_df <- sig_results
  tmp_df$diffexp <- 'No'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
  
  tmp_df <- tmp_df %>% arrange(Pvalue)
  tmp_df$label <- NA
  tmp_df$label[1:20] <- tmp_df$Gene[1:20]
  
  tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)
  
  if(length(unique(tmp_df$diffexp)) > 1){
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                           title = paste0('Type 2 Diabetes Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  }else{
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                           title = paste0('Type 2 Diabetes Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
    
  }
  pdf(paste0(results.dir, 'NEBULA_T2D_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
  print(tmp_graph)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  

#TAL Cells
  so_subset <- subset(so, TAL_celltype =='TAL')
  
  cellcounts <- ncol(so_subset)
  partcounts <- length(unique(so_subset@meta.data$kit_id))
  
  celltype2 <- 'TAL'
  full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(full_analysis)
  full_analysis <- subset(so_subset, features = hvgs)
  full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 
  
  
  meta_gene <- full_analysis@meta.data
  pred_gene <- model.matrix(~sex, data = meta_gene)
  data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
  result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                            offset = full_analysis$pooled_offset,
                            ncore = 1, output_re = T, covariance = T,
                            reml = T, model = "NBLMM")
  
  result_allcells$summary
  
  write.table(result_allcells, paste(results.dir, 'NEBULA_LC_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
              quote=F, sep='\t')
  
  
  #Plotting 
  sig_results <- as.data.frame(result_allcells) %>% 
    dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                  Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)
  
  tmp_df <- sig_results
  tmp_df$diffexp <- 'No'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
  
  tmp_df <- tmp_df %>% arrange(Pvalue)
  tmp_df$label <- NA
  tmp_df$label[1:20] <- tmp_df$Gene[1:20]
  
  tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)
  
  if(length(unique(tmp_df$diffexp)) > 1){
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                           title = paste0('Lean Controls Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  }else{
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                           title = paste0('Lean Controls Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
    
  }
  pdf(paste0(results.dir, 'NEBULA_LC_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
  print(tmp_graph)
  dev.off()
  
  
  
  
  
  
  

#DCT Cells
  so_subset <- subset(so, celltype2 == 'DCT')
  
  cellcounts <- ncol(so_subset)
  partcounts <- length(unique(so_subset@meta.data$kit_id))
  
  celltype2 <- 'DCT'
  full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(full_analysis)
  full_analysis <- subset(so_subset, features = hvgs)
  full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 
  
  
  meta_gene <- full_analysis@meta.data
  pred_gene <- model.matrix(~sex, data = meta_gene)
  data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
  result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                            offset = full_analysis$pooled_offset,
                            ncore = 1, output_re = T, covariance = T,
                            reml = T, model = "NBLMM")
  
  result_allcells$summary
  
  write.table(result_allcells, paste(results.dir, 'NEBULA_LC_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
              quote=F, sep='\t')
  
  
  #Plotting 
  sig_results <- as.data.frame(result_allcells) %>% 
    dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                  Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)
  
  tmp_df <- sig_results
  tmp_df$diffexp <- 'No'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
  
  tmp_df <- tmp_df %>% arrange(Pvalue)
  tmp_df$label <- NA
  tmp_df$label[1:20] <- tmp_df$Gene[1:20]
  
  tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)
  
  if(length(unique(tmp_df$diffexp)) > 1){
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                           title = paste0('Lean Controls Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  }else{
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                           title = paste0('Lean Controls Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
    
  }
  pdf(paste0(results.dir, 'NEBULA_LC_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
  print(tmp_graph)
  dev.off()
  
  
  

  
  #Cell Subtypes
  celltypes <- unique(so@meta.data$KPMP_celltype)
  
  for(celltype in celltypes){
  so_subset <- subset(so, KPMP_celltype == celltype)
  
  cellcounts <- ncol(so_subset)
  partcounts <- length(unique(so_subset@meta.data$kit_id))
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(full_analysis)
  full_analysis <- subset(so_subset, features = hvgs)
  full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 
  
  
  meta_gene <- full_analysis@meta.data
  pred_gene <- model.matrix(~sex, data = meta_gene)
  data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
  result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                            offset = full_analysis$pooled_offset,
                            ncore = 1, output_re = T, covariance = T,
                            reml = T, model = "NBLMM")
  
  result_allcells$summary
  
  write.table(result_allcells, paste(results.dir, 'NEBULA_LC_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
              quote=F, sep='\t')
  
  
  #Plotting 
  sig_results <- as.data.frame(result_allcells) %>% 
    dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                  Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)
  
  tmp_df <- sig_results
  tmp_df$diffexp <- 'No'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
  
  tmp_df <- tmp_df %>% arrange(Pvalue)
  tmp_df$label <- NA
  tmp_df$label[1:20] <- tmp_df$Gene[1:20]
  
  tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)
  
  if(length(unique(tmp_df$diffexp)) > 1){
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                           title = paste0('Lean Controls Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  }else{
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                           title = paste0('Lean Controls Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
    
  }
  pdf(paste0(results.dir, 'NEBULA_LC_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
  print(tmp_graph)
  dev.off()
  
  print(celltype2)
  
  
}












#T2D Analyses
remove(list=ls())


load('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line265.RData')



results.dir <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/MainEffects_Sex/'
  
  
  
so <- subset(so, group == 'Type_2_Diabetes')



#All Analysis
so_subset <- so

cellcounts <- ncol(so_subset)
partcounts <- length(unique(so_subset@meta.data$kit_id))

celltype2 <- 'All'
full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(full_analysis)
full_analysis <- subset(so_subset, features = hvgs)
full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 


meta_gene <- full_analysis@meta.data
pred_gene <- model.matrix(~sex, data = meta_gene)
data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                          offset = full_analysis$pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary

write.table(result_allcells, paste(results.dir, 'NEBULA_T2D_SexDifferences_AllCells.txt'), row.names=F,
            quote=F, sep='\t')


#Plotting 
sig_results <- as.data.frame(result_allcells) %>% 
  dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)

tmp_df <- sig_results
tmp_df$diffexp <- 'No'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'

tmp_df <- tmp_df %>% arrange(Pvalue)
tmp_df$label <- NA
tmp_df$label[1:20] <- tmp_df$Gene[1:20]

tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)

if(length(unique(tmp_df$diffexp)) > 1){
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('orange', 'grey', 'purple'),
                       labels = c('Downregulated', 'Not significant', 'Upregulated'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                         title = paste0('Type 2 Diabetes Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
}else{
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('grey'),
                       labels = c('Not significant'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                         title = paste0('Type 2 Diabetes Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  
}
pdf(paste0(results.dir, 'NEBULA_T2D_SexDifferences_Volcano_AllCells.pdf'))
print(tmp_graph)
dev.off()





#PT Cells
so_subset <- subset(so,celltype2 == 'PT')

cellcounts <- ncol(so_subset)
partcounts <- length(unique(so_subset@meta.data$kit_id))

celltype2 <- 'PT'
full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(full_analysis)
full_analysis <- subset(so_subset, features = hvgs)
full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 


meta_gene <- full_analysis@meta.data
pred_gene <- model.matrix(~sex, data = meta_gene)
data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                          offset = full_analysis$pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary

write.table(result_allcells, paste(results.dir, 'NEBULA_T2D_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
            quote=F, sep='\t')


#Plotting 
sig_results <- as.data.frame(result_allcells) %>% 
  dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)

tmp_df <- sig_results
tmp_df$diffexp <- 'No'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'

tmp_df <- tmp_df %>% arrange(Pvalue)
tmp_df$label <- NA
tmp_df$label[1:20] <- tmp_df$Gene[1:20]

tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)

if(length(unique(tmp_df$diffexp)) > 1){
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('orange', 'grey', 'purple'),
                       labels = c('Downregulated', 'Not significant', 'Upregulated'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                         title = paste0('Type 2 Diabetes Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
}else{
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('grey'),
                       labels = c('Not significant'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                         title = paste0('Type 2 Diabetes Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  
}
pdf(paste0(results.dir, 'NEBULA_T2D_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
print(tmp_graph)
dev.off()








#TAL Cells
so_subset <- subset(so, TAL_celltype =='TAL')

cellcounts <- ncol(so_subset)
partcounts <- length(unique(so_subset@meta.data$kit_id))

celltype2 <- 'TAL'
full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(full_analysis)
full_analysis <- subset(so_subset, features = hvgs)
full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 


meta_gene <- full_analysis@meta.data
pred_gene <- model.matrix(~sex, data = meta_gene)
data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                          offset = full_analysis$pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary

write.table(result_allcells, paste(results.dir, 'NEBULA_T2D_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
            quote=F, sep='\t')


#Plotting 
sig_results <- as.data.frame(result_allcells) %>% 
  dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)

tmp_df <- sig_results
tmp_df$diffexp <- 'No'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'

tmp_df <- tmp_df %>% arrange(Pvalue)
tmp_df$label <- NA
tmp_df$label[1:20] <- tmp_df$Gene[1:20]

tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)

if(length(unique(tmp_df$diffexp)) > 1){
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('orange', 'grey', 'purple'),
                       labels = c('Downregulated', 'Not significant', 'Upregulated'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                         title = paste0('Type 2 Diabetes Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
}else{
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('grey'),
                       labels = c('Not significant'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                         title = paste0('Type 2 Diabetes Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  
}
pdf(paste0(results.dir, 'NEBULA_T2D_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
print(tmp_graph)
dev.off()








#DCT Cells
so_subset <- subset(so, celltype2 == 'DCT')

cellcounts <- ncol(so_subset)
partcounts <- length(unique(so_subset@meta.data$kit_id))

celltype2 <- 'DCT'
full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(full_analysis)
full_analysis <- subset(so_subset, features = hvgs)
full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 


meta_gene <- full_analysis@meta.data
pred_gene <- model.matrix(~sex, data = meta_gene)
data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                          offset = full_analysis$pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary

write.table(result_allcells, paste(results.dir, 'NEBULA_T2D_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
            quote=F, sep='\t')


#Plotting 
sig_results <- as.data.frame(result_allcells) %>% 
  dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)

tmp_df <- sig_results
tmp_df$diffexp <- 'No'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'

tmp_df <- tmp_df %>% arrange(Pvalue)
tmp_df$label <- NA
tmp_df$label[1:20] <- tmp_df$Gene[1:20]

tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)

if(length(unique(tmp_df$diffexp)) > 1){
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('orange', 'grey', 'purple'),
                       labels = c('Downregulated', 'Not significant', 'Upregulated'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                         title = paste0('Type 2 Diabetes Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
}else{
  tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('grey'),
                       labels = c('Not significant'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                         title = paste0('Type 2 Diabetes Sex Differences'), 
                         subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  
}
pdf(paste0(results.dir, 'NEBULA_T2D_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
print(tmp_graph)
dev.off()





#Cell Subtypes
celltypes <- unique(so@meta.data$KPMP_celltype)

for(celltype in celltypes){
  so_subset <- subset(so, KPMP_celltype == celltype)
  
  cellcounts <- ncol(so_subset)
  partcounts <- length(unique(so_subset@meta.data$kit_id))
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  full_analysis <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(full_analysis)
  full_analysis <- subset(so_subset, features = hvgs)
  full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 
  
  
  meta_gene <- full_analysis@meta.data
  pred_gene <- model.matrix(~sex, data = meta_gene)
  data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
  result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                            offset = full_analysis$pooled_offset,
                            ncore = 1, output_re = T, covariance = T,
                            reml = T, model = "NBLMM")
  
  result_allcells$summary
  
  write.table(result_allcells, paste(results.dir, 'NEBULA_T2D_SexDifferences_', celltype2, 'Cells.txt'), row.names=F,
              quote=F, sep='\t')
  
  
  #Plotting 
  sig_results <- as.data.frame(result_allcells) %>% 
    dplyr::select(Gene = summary.gene, LogFC = summary.logFC_sexMale, 
                  Pvalue = summary.p_sexMale) #%>% filter(Pvalue < 0.05)
  
  tmp_df <- sig_results
  tmp_df$diffexp <- 'No'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
  
  tmp_df <- tmp_df %>% arrange(Pvalue)
  tmp_df$label <- NA
  tmp_df$label[1:20] <- tmp_df$Gene[1:20]
  
  tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)
  
  if(length(unique(tmp_df$diffexp)) > 1){
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                           title = paste0('Type 2 Diabetes Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
  }else{
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                           title = paste0('Type 2 Diabetes Sex Differences'), 
                           subtitle = paste0(cellcounts, ' ', celltype2, ' Cells; ', partcounts, ' Participants'))
    
  }
  pdf(paste0(results.dir, 'NEBULA_T2D_SexDifferences_Volcano_', celltype2, 'Cells.pdf'))
  print(tmp_graph)
  dev.off()
  
  print(celltype2)
  
  
}






#Clean results 


summary_files <- list.files('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/MainEffects_Sex/', pattern='txt')


results.dir <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/MainEffects_Sex/'
for(i in c(1:length(summary_files))){
  tmp_df <- data.table::fread(paste0(results.dir, summary_files[i]))
new_file_name <- str_replace(summary_files[i], '.txt', '')

  tmp_df <- tmp_df %>% 
    dplyr::select(gene = summary.gene, 
                  logFC = summary.logFC_sexMale, 
                  pvalue = summary.p_sexMale)
  
  write.table(tmp_df, paste0(results.dir, new_file_name, '_cleaned.csv'),
              row.names=F, quote = F, sep=',')
  
  sig_genes_small <- tmp_df %>% filter(pvalue < 0.05)
  if(nrow(sig_genes_small) == 0){
    next
  }
  enrich_GO_BP <- enrichGO(gene = sig_genes_small$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP") %>% 
    #  as.data.frame() %>% #dplyr::select(Description, GeneRatio, p.adjust, Count) %>% 
    filter(p.adjust < 0.05)
  enrich_GO_MF <- enrichGO(gene = sig_genes_small$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF") %>% 
    filter(p.adjust < 0.05)
  enrich_GO_CC <- enrichGO(gene = sig_genes_small$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC") %>% 
    filter(p.adjust < 0.05)
  
  #dotplot
  jpeg(paste0(results.dir, 'GSEA/', new_file_name, '_significant_GSEA.jpeg'), height = 800, width = 600)
  if(nrow(as.data.frame(enrich_GO_BP) > 0)){
  dotplot(enrich_GO_BP, showCategory = 20)
  }
  if(nrow(as.data.frame(enrich_GO_MF) > 0)){
  dotplot(enrich_GO_MF, showCategory = 20)
  }
  if(nrow(as.data.frame(enrich_GO_CC) > 0)){
  dotplot(enrich_GO_CC, showCategory = 20)
  }
  
  dev.off()
  
  
  
  
  
  
  
  
  print(i)
  
}








