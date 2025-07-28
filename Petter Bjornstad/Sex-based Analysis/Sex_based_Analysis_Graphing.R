#UMAP Plots

genes_list <- rownames(so)
so <- so %>% ScaleData(features=genes_list)


so <- RunUMAP(so, features = VariableFeatures(so))





# Extract UMAP coordinates and metadata
umap_coords <- Embeddings(so, reduction = "umap")
metadata <- so@meta.data

# Combine into a data frame
plot_data <- data.frame(
  UMAP_1 = umap_coords[,1],
  UMAP_2 = umap_coords[,2],
  sex = metadata$sex,
  group = metadata$group,
  # Add any other variables you want to color by
  clusters = Idents(so)
)

# Create the faceted plot

pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Sex_Group_Based_UMAP.pdf',
    width = 15, height = 20)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
  geom_point(size = 0.5, alpha = 0.7) +
  facet_grid(sex ~ group) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")
dev.off()


DimPlot(so, reduction = 'umap')









DimPlot(so, reduction = 'umap')+
  facet_grid(sex ~ group, 
             labeller = labeller(sex = label_both, 
                                 group = label_both))+
  theme(strip.text = element_text(size=10))














#Cell-specific distributions

meta.data <- so@meta.data
meta.data$group_sex <- paste0(meta.data$group, '_', meta.data$sex)

ggplot(meta.data %>% filter(celltype_rpca != 'ATL'), 
       aes(x=group_sex, fill=KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='Cell Type')


ggplot(meta.data %>% filter(celltype_rpca != 'ATL') %>% filter(celltype2 == 'PT'), 
       aes(x=group_sex, fill=KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='PT Cell Type')

ggplot(meta.data %>% filter(celltype_rpca != 'ATL') %>% filter(celltype2 == 'EC'), 
       aes(x=group_sex, fill=KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='EC Type')

ggplot(meta.data %>% filter(celltype_rpca != 'ATL') %>% filter(TAL_celltype == 'TAL'), 
       aes(x=group_sex, fill=KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='TAL Cell Type')





##PT Heatmaps 

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


interaction <- data.table::fread("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/ NEBULA_fullanalysis_offset.txt")
interaction <- interaction %>% dplyr::select(gene = summary.gene,
                                             interaction_logFC = summary.logFC_groupType_2_Diabetes.sexMale,
                                             interaction_pvalue = summary.p_groupType_2_Diabetes.sexMale) %>% 
  filter(interaction_pvalue < 0.05)







#LC Sex Analyses

so <- subset(so, group == 'Lean_Control')


#PT Cells
#so_subset <- subset(so,celltype2 == 'PT')
so_subset <- so 
remove(so)

cellcounts <- ncol(so_subset)
partcounts <- length(unique(so_subset@meta.data$kit_id))

#celltype2 <- 'PT'

counts_t1d_hc <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round
counts_t1d_hc_mtap <- counts_t1d_hc[interaction$gene,]

genes_to_analyze <- interaction$gene
count_data <- GetAssayData(so_subset, slot = "counts")[genes_to_analyze, ]
meta_gene <- so_subset@meta.data


pred_gene <- model.matrix(~sex, data = meta_gene)
data_g_gene <- list(count = counts_t1d_hc_mtap, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, offset = so_subset$pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary



lc <- result_allcells %>% as.data.frame() %>% 
  dplyr::select(gene = summary.gene, 
                lc_logFC = summary.logFC_sexMale,
                lc_pvalue = summary.p_sexMale)



remove(so_subset)

##T2D 
load('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line265.RData')
so <- subset(so, group == 'Type_2_Diabetes')

#PT Cells
so_subset <- subset(so,celltype2 == 'PT')

cellcounts <- ncol(so_subset)
partcounts <- length(unique(so_subset@meta.data$kit_id))

celltype2 <- 'PT'

counts_t1d_hc <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round
counts_t1d_hc_mtap <- counts_t1d_hc[interaction$gene,]

genes_to_analyze <- interaction$gene
count_data <- GetAssayData(so_subset, slot = "counts")[genes_to_analyze, ]
meta_gene <- so_subset@meta.data


pred_gene <- model.matrix(~sex, data = meta_gene)
data_g_gene <- list(count = counts_t1d_hc_mtap, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = data_g_gene$count, id = data_g_gene$id, pred = data_g_gene$pred, offset = so_subset$pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary



t2d <- result_allcells %>% as.data.frame() %>% 
  dplyr::select(gene = summary.gene, 
                t2d_logFC = summary.logFC_sexMale,
                t2d_pvalue = summary.p_sexMale)





combined <- interaction %>% left_join(lc, by='gene') %>% 
  left_join(t2d, by='gene') %>% filter(interaction_pvalue < 0.05)




library(pheatmap)
heatmap_data <- combined %>% dplyr::select(t2d_logFC, lc_logFC, interaction_logFC)   # exclude gene column
rownames(heatmap_data) <- combined$gene

# Create heatmap

pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/LogFCHeatmap_MainEffectsvsInteraction.pdf', 
    height=12, width=10)
pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = FALSE,  # since only 3 conditions
         scale = "none",        # or "row" to scale by gene
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101),  # adjust range as needed
         fontsize_row = 8,      # adjust gene label size
         fontsize_col = 10,     # adjust condition label size
         angle_col = 45,        # rotate column labels
         main = "LogFC Heatmap")
dev.off()


enrich_GO_BP <- enrichGO(gene = interaction$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP") %>% 
  #  as.data.frame() %>% #dplyr::select(Description, GeneRatio, p.adjust, Count) %>% 
  filter(p.adjust < 0.05)

#dotplot
pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/GSEA_SignificantInteraction_dotplot.pdf',
    width=10, height=15)
dotplot(enrich_GO_BP, showCategory = 20)
dev.off()







