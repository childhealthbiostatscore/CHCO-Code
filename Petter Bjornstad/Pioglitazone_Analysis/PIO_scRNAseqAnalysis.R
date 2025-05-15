#scRNAseq analysis
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
library(foreach)
library(doParallel)
library(nebula)
library(Matrix)
library(DirichletReg)



setwd('C:/Users/netio/Documents/Harmonized_data/')


#Demographics and identifying medication use 

harmonized_data <- data.table::fread("harmonized_dataset.csv")
harmonized_data <- harmonized_data %>% filter(group == 'Type 2 Diabetes')

medications <- readxl::read_xlsx("Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>% dplyr::select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), pattern = '_1', replacement = '')








#Load Data and take a look 

scrna <- readRDS('../../Downloads/PB_90samples_Harmony_rpca_Fadhl_PhilApproved_091024.RDS')

meta.data <- scrna@meta.data

#Get recordID filters 

harmonized_data <- harmonized_data %>% 
  filter(rh_id %in% meta.data$record_id | croc_id %in% meta.data$record_id | improve_id %in% meta.data$record_id | penguin_id %in% meta.data$record_id | rh2_id %in% meta.data$record_id) %>%
  dplyr::select(record_id, mrn, group)

medications <- medications %>% dplyr::select(mrn, tzd) %>% filter(mrn %in% harmonized_data$mrn)
harmonized_data$mrn <- as.character(harmonized_data$mrn)

final_df <- medications %>% left_join(harmonized_data, by='mrn')
final_df$combined_id <- paste0(final_df$mrn, '_', final_df$record_id)
final_df <- final_df %>% filter(!duplicated(combined_id))

#Quality Control 


scrna_small <- subset(x = scrna,
                      record_id %in% final_df$record_id)

scrna_small <- subset(scrna_small, 
                      subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)


#Assigning the groups 
meta.data <- scrna_small@meta.data
new_df <- data.frame(ID = meta.data$record_id)

group_data <- final_df %>% dplyr::select(ID=record_id, group=tzd)

new_df <- new_df %>% left_join(group_data)

scrna_small$group_labels <- new_df$group
meta.data <- scrna_small@meta.data



#Traditional Analysis 
scrna_small <- NormalizeData(scrna_small)
scrna_small <- FindVariableFeatures(scrna_small, selection.method = 'vst', nfeatures=2000)




#
genes_list <- rownames(scrna_small)
scrna_small <- scrna_small %>% ScaleData(features=genes_list)


scrna_small <- RunUMAP(scrna_small, features = VariableFeatures(scrna_small))

#Plotting the UMAP 
pdf('../UofW/pio_scrna/AllCellTypes/Traditional_UMAP_plot.pdf', width=20, height=15)
DimPlot(scrna_small, reduction = 'umap')
dev.off()



#Need to code in putting the cell types in each box so it's more readable
pdf("../UofW/pio_scrna/AllCellTypes/CellType_Distributions.pdf", width=20, height=15)
ggplot(meta.data, aes(x=group_labels, fill=celltype_rpca, label = celltype_rpca))+
  geom_bar(position='fill')+theme_classic()+labs(x='Pio Treatment', fill='Cell Type')

dev.off()

#Add code to only show ones I think are significantly different

#Testing significant difference 
chisq.test(meta.data$group_labels, meta.data$celltype_rpca)

#Full Comparison irrespective of cell type 
sig_markers <- FindMarkers(scrna_small, ident.1 = 'Yes', 
                           ident.2 = 'No', group.by = 'group_labels')

sig_markers$diffexp <- 'No'
sig_markers$diffexp[sig_markers$avg_log2FC > 0.6 & sig_markers$p_val_adj < 0.05] <- 'Up'
sig_markers$diffexp[sig_markers$avg_log2FC < -0.6 & sig_markers$p_val_adj < 0.05] <- 'Down'

sig_markers$label <- NA
sig_markers$label[1:10] <- rownames(sig_markers)[1:10]

sig_markers$gene <- rownames(sig_markers)

#Making graph
pdf('../UofW/pio_scrna/AllCellTypes/SignificantExpressionDifferences_AllCells.pdf')
ggplot(sig_markers, aes(x= avg_log2FC, y=-log10(p_val_adj), col = diffexp, label=label))+
  geom_point()+
  geom_text(size=2, vjust = 2, color='black')+
  scale_color_manual(values = c('orange', 'grey', 'purple'),
                     labels = c('Downregulated', 'Not significant', 'Upregulated'))+
  geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
  geom_vline(xintercept = c(-0.6, 0.6), col='blue', linetype ='dashed')+
  theme_classic()+labs(x='Log2FC', y='-log10 pvalue', col ='Differential Expression', 
                       title = 'Use of Pioglitazone: Changes in Expression in All Cells')

dev.off()

write.table(sig_markers %>% filter(diffexp != 'No'), 
            '../UofW/pio_scrna/AllCellTypes/SignificantMarkers_AllCellTypes.txt', 
            row.names=F, quote=F, sep='\t')




#Each cell type specifically
scrna_small$celltype_group <- paste0(scrna_small$celltype_rpca, 
                                     sep='-', scrna_small$group_labels)
Idents(scrna_small) <- "celltype_group"

#Graph all the PT subtypes 
graph_list <- list()

celltype_list <- unique(scrna_small$celltype_rpca)







#NEBULA

counts_layer <- round(GetAssayData(scrna_small, layer = 'counts'))
library_size <- Matrix::colSums(round(GetAssayData(scrna_small, layer = 'counts')))
scrna_small$library_size <- library_size
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_layer))
 sce <- computeSumFactors(sce)
 # View size factors
sizeFactors(sce)
## Calculate offset → (size factors)
pooled_offset <- (sizeFactors(sce))






scrna_small <- FindVariableFeatures(scrna_small, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(scrna_small)

#Perform remaining steps on top 2000 hvgs
# Subset Seurat object to only HVGs
scrna_small_hvg <- subset(scrna_small, features = hvgs)
counts_layer <- round(GetAssayData(scrna_small_hvg, layer = 'counts'))

pdf('../../Documents/UofW/pio_scrna/AllCellTypes/UMAP_GroupLabels_PIO.pdf')
DimPlot(scrna_small, reduction = 'umap.harmony', group.by = 'celltype_harmony')
dev.off()



meta_gene <- scrna_small_hvg@meta.data
pred_gene <- model.matrix(~group_labels, data = meta_gene)
data_g_gene <- group_cell(count = counts_layer, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, 
                          offset = pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary

result_allcells <- as.data.frame(result_allcells)

write.table(result_allcells, '../../Documents/UofW/pio_scrna/AllCellTypes/AllCells_NEBULA_PIO_analysis.txt',
            row.names=F, quote=F, sep='\t')

#Plot the Overall Scores 

sig_markers <- result_allcells %>% dplyr::select(avg_log2FC = summary.logFC_group_labelsYes, 
                                          p_val_raw = summary.p_group_labelsYes,
                                          gene = summary.gene)

sig_markers <- sig_markers %>% mutate(p_val_adj = p.adjust(p_val_raw, method='BH')) %>%
  arrange(p_val_raw)

sig_markers$diffexp <- 'No'
sig_markers$diffexp[sig_markers$avg_log2FC > 0.6 & sig_markers$p_val_adj < 0.05] <- 'Up'
sig_markers$diffexp[sig_markers$avg_log2FC < -0.6 & sig_markers$p_val_adj < 0.05] <- 'Down'

sig_markers$label <- NA
sig_markers$label[1:15] <- sig_markers$gene[1:15]

#Raw P-values 
sig_markers_raw <- sig_markers
sig_markers_raw$diffexp <- 'No'
sig_markers_raw$diffexp[sig_markers_raw$avg_log2FC > 0.6 & sig_markers_raw$p_val_raw < 0.05] <- 'Up'
sig_markers_raw$diffexp[sig_markers_raw$avg_log2FC < -0.6 & sig_markers_raw$p_val_raw < 0.05] <- 'Down'

sig_markers_raw$label <- NA
sig_markers_raw$label[1:15] <- sig_markers_raw$gene[1:15]


#Making graph
pdf('../../Documents/UofW/pio_scrna/AllCellTypes/NEBULA_SignificantGenes_PIO_noCorrections.pdf')
ggplot(sig_markers_raw, aes(x= avg_log2FC, y=-log10(p_val_raw), col = diffexp, label=label))+
  geom_point()+
  geom_text(size=2, vjust = 2, color='black')+
  scale_color_manual(values = c('orange', 'grey', 'purple'),
                     labels = c('Downregulated', 'Not significant', 'Upregulated'))+
  geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
  geom_vline(xintercept = c(-0.6, 0.6), col='blue', linetype ='dashed')+
  theme_classic()+labs(x='Log2FC', y='-log10 pvalue', col ='Differential Expression', 
                       title = 'Expression Changes In Pioglitazone Use In All Cell Types')

dev.off()

pdf('../../Documents/UofW/pio_scrna/AllCellTypes/NEBULA_SignificantGenes_PIO_BHCorrected.pdf')
ggplot(sig_markers, aes(x= avg_log2FC, y=-log10(p_val_adj), col = diffexp, label=label))+
  geom_point()+
  geom_text(size=2, vjust = 2, color='black')+
  scale_color_manual(values = c('orange', 'grey', 'purple'),
                     labels = c('Downregulated', 'Not significant', 'Upregulated'))+
  geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
  geom_vline(xintercept = c(-0.6, 0.6), col='blue', linetype ='dashed')+
  theme_classic()+labs(x='Log2FC', y='-log10 pvalue', col ='Differential Expression', 
                       title = 'Expression Changes In Pioglitazone Use In All Cell Types')

dev.off()





write.table(sig_markers %>% filter(diffexp != 'No'), 
            '../../Documents/UofW/pio_scrna/AllCellTypes/NEBULA_PIO_SignificantMarkers_AllCellTypes.txt', 
            row.names=F, quote=F, sep='\t')

remove(sig_markers)


#Cell-Type Specific 

celltypes <- as.character(unique(scrna_small$celltype_harmony))

counts_layer <- round(GetAssayData(scrna_small, layer = 'counts'))
library_size <- Matrix::colSums(round(GetAssayData(scrna_small, layer = 'counts')))
scrna_small$library_size <- library_size
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_layer))
sce <- computeSumFactors(sce)
# View size factors
sizeFactors(sce)
## Calculate offset → (size factors)
scrna_small$pooled_offset <- (sizeFactors(sce))


for(i in c(1:length(celltypes))){
  tmp_celltype <- celltypes[i]
  if(str_detect(tmp_celltype, 'lowQuality')){
    next
  }
  celltype_name <- tmp_celltype
  if(str_detect(celltype_name, pattern = '/')){
    celltype_name <- str_replace_all(celltype_name, pattern='/', replacement = '_')
  }
  scrna_temp <- subset(scrna_small, celltype_harmony == tmp_celltype)
  scrna_temp <- FindVariableFeatures(scrna_temp, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(scrna_temp)
  counts_layer <- round(GetAssayData(scrna_temp, layer = 'counts'))
  
  #Perform remaining steps on top 2000 hvgs
  # Subset Seurat object to only HVGs
  scrna_temp_hvg <- subset(scrna_temp, features = hvgs)
  
  pdf(paste0('../../Documents/UofW/pio_scrna/CelllTypeSpecific/UMAP_PIO_', 
             celltype_name, '.pdf'))
  DimPlot(scrna_temp, reduction = 'umap.harmony', group.by = 'celltype_harmony')
  dev.off()
  
  
  #NEBULA Analysis 
  meta_gene <- scrna_temp_hvg@meta.data
  pred_gene <- model.matrix(~group_labels, data = meta_gene)
  data_g_gene <- group_cell(count = counts_layer, id = meta_gene$record_id, pred = pred_gene)
  result_allcells <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, 
                            offset = scrna_temp$pooled_offset,
                            ncore = 1, output_re = T, covariance = T,
                            reml = T, model = "NBLMM")
  
  result_allcells <- as.data.frame(result_allcells)
  
  write.table(result_allcells, 
              paste0('../../Documents/UofW/pio_scrna/CelllTypeSpecific/NEBULA_PIO_', 
                     celltype_name, '_analysis.txt'),
              row.names=F, quote=F, sep='\t')
  
  
  #Plotting for significant markers 
  sig_markers <- result_allcells %>% dplyr::select(avg_log2FC = summary.logFC_group_labelsYes, 
                                                   p_val_raw = summary.p_group_labelsYes,
                                                   gene = summary.gene)
  
  sig_markers <- sig_markers %>% mutate(p_val_adj = p.adjust(p_val_raw, method='BH')) %>%
    arrange(p_val_raw)
  
  sig_markers$diffexp <- 'No'
  sig_markers$diffexp[sig_markers$avg_log2FC > 0.6 & sig_markers$p_val_adj < 0.05] <- 'Up'
  sig_markers$diffexp[sig_markers$avg_log2FC < -0.6 & sig_markers$p_val_adj < 0.05] <- 'Down'
  
  sig_markers$label <- NA
  sig_markers$label[1:15] <- sig_markers$gene[1:15]
  
  #Raw P-values 
  sig_markers_raw <- sig_markers
  sig_markers_raw$diffexp <- 'No'
  sig_markers_raw$diffexp[sig_markers_raw$avg_log2FC > 0.6 & sig_markers_raw$p_val_raw < 0.05] <- 'Up'
  sig_markers_raw$diffexp[sig_markers_raw$avg_log2FC < -0.6 & sig_markers_raw$p_val_raw < 0.05] <- 'Down'
  
  sig_markers_raw$label <- NA
  sig_markers_raw$label[1:15] <- sig_markers_raw$gene[1:15]
  
  
  #Making graph
  pdf(paste0('../../Documents/UofW/pio_scrna/CelllTypeSpecific/NEBULA_SignificantGenes_PIO_noCorrections_', 
  celltype_name, '.pdf'))
  ggplot(sig_markers_raw %>% filter(!is.na(p_val_raw)), aes(x= avg_log2FC, y=-log10(p_val_raw), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('orange', 'grey', 'purple'),
                       labels = c('Downregulated', 'Not significant', 'Upregulated'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(-0.6, 0.6), col='blue', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                         title = paste0('Expression Changes In Pioglitazone Use In ', tmp_celltype))
  
  dev.off()
  
  pdf(paste0('../../Documents/UofW/pio_scrna/CelllTypeSpecific/NEBULA_SignificantGenes_PIO_BHCorrected_', 
  celltype_name, '.pdf'))
  if(length(unique(sig_markers$diffexp )) > 1){
  tmp_graph <- ggplot(sig_markers %>% filter(!is.na(p_val_adj)), aes(x= avg_log2FC, y=-log10(p_val_adj), col = diffexp, label=label))+
    geom_point()+
    geom_text(size=2, vjust = 2, color='black')+
    scale_color_manual(values = c('orange', 'grey', 'purple'),
                       labels = c('Downregulated', 'Not significant', 'Upregulated'))+
    geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
    geom_vline(xintercept = c(-0.6, 0.6), col='blue', linetype ='dashed')+
    theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                         title = paste0('Expression Changes In Pioglitazone Use In ', tmp_celltype))
  }else{
    tmp_graph <- ggplot(sig_markers %>% filter(!is.na(p_val_adj)), aes(x= avg_log2FC, y=-log10(p_val_adj), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not Significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(-0.6, 0.6), col='blue', linetype ='dashed')+
      theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                           title = paste0('Expression Changes In Pioglitazone Use In ', tmp_celltype))
  }
  tmp_graph
  dev.off()
  
  
  
  
  
  write.table(sig_markers %>% filter(diffexp != 'No'), 
              paste0('../../Documents/UofW/pio_scrna/CelllTypeSpecific/NEBULA_PIO_SignificantMarkers_', 
                     celltype_name, '.txt'), 
              row.names=F, quote=F, sep='\t')
  
  remove(sig_markers)
  
  print(i)
  
  
  
}

#Combined Cell types 

#EC

#PT

#TAL 





