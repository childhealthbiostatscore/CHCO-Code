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
#library(table1)
library(clusterProfiler)
library('org.Hs.eg.db')



#Load in single-cell RNAseq data and do QC


so <- readRDS("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/scRNA/data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds")
mito_genes <- grep("^MT-", rownames(so), value = TRUE)
so <- subset(so, features = setdiff(rownames(so), mito_genes))

ribo_genes <- c(
  "RPL22", "RPL11", "RPS8", "RPL5", "RPS27", "RPS7", "RPS27A", "RPL31", "RPL37A", "RPL32", "RPL15", "RPL14", "RPL29",
  "RPL24", "RPL22L1", "RPL35A", "RPL9", "RPL34", "RPS3A", "RPL37", "RPS23", "RPS14", "RPS18", "RPS10", "RPL10A", 
  "RPS20", "RPL7", "RPL30", "RPL8", "RPS6", "RPL35", "RPL12", "RPL7A", "RPS24", "RPLP2", "RPL27A", "RPS13", "RPS3",
  "RPS25", "RPS26", "RPL41", "RPL6", "RPLP0", "RPL21", "RPS29", "RPL4", "RPLP1", "RPS17", "RPS2", "RPS15A", "RPL13",
  "RPL26", "RPL23A", "RPL23", "RPL19", "RPL27", "RPL38", "RPL17", "RPS15", "RPL36", "RPS28", "RPL18A", "RPS16", 
  "RPS19", "RPL18", "RPL13A", "RPS11", "RPS9", "RPL28", "RPS5", "RPS21", "RPL3", "RPS4X", "RPL36A", "RPL39", 
  "RPL10", "RPS4Y1")

so<- subset(so, features = setdiff(rownames(so), ribo_genes))



so <- NormalizeData(so)
so <- ScaleData(so, features = VariableFeatures(so))





#Organize data groups

meta.data <- so@meta.data

so$kit_id[which(so$kit_id=="KI-0014643")] <- "KL-0014643"
so$kit_id[which(so$kit_id=="kl-0023998")] <- "KL-0023998"


harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')
dat <- harmonized_data %>%
  arrange(screen_date) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat <- dat %>% filter(visit == 'baseline')

dat <- dat %>% 
  filter(group=="Type 2 Diabetes" | group=="Lean Control" | group == 'Type 1 Diabetes' | group == 'Obese Control')
length(unique(dat$mrn))#22
length(unique(dat$record_id))#22

dat <- dat %>% 
  filter(kit_id %in% so$kit_id)
length(unique(dat$mrn))#4
length(unique(dat$record_id))#4

ids <- c(dat$kit_id) #32 total ids, 13 LC, 19 T2D no SGLT2
length(ids) #32

#Filter pb90 to these 32 kit ids for aim 1
so <- subset(so, kit_id %in% ids)





save.image('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line195.RData')

#Sex-specific analyses 


counts_layer <- round(GetAssayData(so, layer = 'counts'))
library_size <- Matrix::colSums(round(GetAssayData(so, layer = 'counts')))
so$library_size <- library_size
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_layer))
sce <- computeSumFactors(sce)

## Calculate offset → (size factors)
so$pooled_offset <- (sizeFactors(sce))



save.image('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line210.RData')

so$celltype1 <- case_when(grepl("PT-",so$celltype_rpca)~"PT",
                                  grepl("TAL-",so$celltype_rpca)~"TAL",
                                  grepl("EC-",so$celltype_rpca)~"EC",
                                  grepl("POD",so$celltype_rpca)~"POD",
                                  grepl("MAC",so$celltype_rpca)~"MAC",
                                  grepl("MON",so$celltype_rpca)~"MON",
                                  grepl("PC-",so$celltype_rpca)~"PC",
                                  grepl("FIB",so$celltype_rpca)~"FIB_MC_VSMC",
                                  grepl("DTL",so$celltype_rpca)~"DTL",
                                  so$celltype_rpca=="DCT"~"DCT",
                                  so$celltype_rpca=="ATL"~"ATL",
                                  so$celltype_rpca=="B"~"B",
                                  so$celltype_rpca=="T"~"T")
so$celltype1 <- as.character(so$celltype1)

so$KPMP_celltype2 <- as.character(so$KPMP_celltype)
so$celltype2 <- ifelse(so$KPMP_celltype=="aPT" | 
                                 so$KPMP_celltype=="PT-S1/S2" | 
                                 so$KPMP_celltype == "PT-S3","PT",
                               ifelse(grepl("TAL",so$KPMP_celltype),"TAL",
                                      ifelse(grepl("EC-",so$KPMP_celltype),"EC",so$KPMP_celltype2)))

so$KPMP_celltype2 <- as.character(so$KPMP_celltype)
so$celltype2 <- ifelse(so$KPMP_celltype=="aPT" | 
                                 so$KPMP_celltype=="PT-S1/S2" | 
                                 so$KPMP_celltype == "PT-S3","PT",
                               ifelse(grepl("TAL",so$KPMP_celltype),"TAL",
                                      ifelse(grepl("EC-",so$KPMP_celltype),"EC",so$KPMP_celltype2)))
#Make sure TAL includes all types
unique(so$KPMP_celltype)
so$TAL_celltype <- ifelse((so$KPMP_celltype=="C-TAL-1" | 
                                     so$KPMP_celltype=="C-TAL-2"|
                                     so$KPMP_celltype=="dTAL" |
                                     so$KPMP_celltype=="aTAL"), "TAL","Non-TAL")






save.image('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line255.RData')



#load('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line255.RData')


#Expression filtering

#expr_matrix <- as.matrix(GetAssayData(so, assay = "RNA", layer = "counts"))
# expr_matrix <- so_kpmp_sc@assays$RNA@counts
# Calculate the proportion of cells expressing each gene
#num_cells_per_gene <- rowSums(expr_matrix > 0)  # Count nonzero cells per gene
#total_cells <- ncol(expr_matrix)  # Total number of cells
#gene_proportion <- num_cells_per_gene / total_cells  # Fraction of cells expressing each gene
#remove(expr_matrix)
# Keep genes expressed in at least "gene_pct" of cells
#genes_to_keep <- names(gene_proportion[gene_proportion >= 0.05])
#so <- subset(so, features = genes_to_keep)




library(Matrix)

# Get sparse count matrix (don't convert to dense!)
counts_sparse <- GetAssayData(so, layer = "counts")

# Calculate detection rate per gene (% cells expressing)
detection_rate <- rowMeans(counts_sparse > 0)

# Calculate mean expression per gene
mean_expr <- rowMeans(counts_sparse)

# Filter genes based on detection rate and mean expression
keep_genes <- detection_rate >= 0.05 & mean_expr >= 0.1  # Adjust thresholds

# Subset your Seurat object
so <- so[keep_genes, ]




save.image('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line265.RData')



##Start here
load('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line265.RData')

pdf("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/CelltypeDistribution_bySex.pdf")
ggplot(meta.data %>% filter(celltype_rpca != 'ATL'), aes(x=group_labels, fill=celltype_rpca))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='Cell Type')

dev.off()


so <- subset(so, group %in% c('Type_2_Diabetes', 'Lean_Control'))

#Function
scrna_analysis <- function(status = 'Both', Controls = 'Both', genelist = NULL,
                          results.dir, celltype= 'All'){
  
  
  
  if(status == 'Both'){
    group_label <- c('Type_1_Diabetes', 'Type_2_Diabetes')
  }else if(status == 'T2D'){
    group_label <- c('Type_2_Diabetes')
  }else if(status == 'T1D'){
    group_label <- c('Type_1_Diabetes')
  }else{
    print('Error: No disease status selected')
    stop()
  }
  
  
  if(Controls == 'Both'){
    group_label <- c(group_label, 'Lean_Control', 'Obese_Control')
  }else if(Controls == 'LC'){
    group_label <- c(group_label, 'Lean_Control')
  }else if(Controls == 'OC'){
    group_label <- c(group_label, 'Obese_Control')
  }
  
  so_subset <- subset(so, group %in% group_label)
  remove(so)
  
  if(celltype == 'All'){
    so_celltype <- so_subset
    }else if(celltype == 'PT'){
    so_celltype <- subset(so_subset,celltype2 == celltype)
    cat('PT Cells')
    }else if(celltype == 'EC'){
   so_celltype <- subset(so_subset, celltype2 == celltype)
   cat("Endothelial Cells")
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
  
  
  
  
full_analysis <- FindVariableFeatures(so_celltype, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(full_analysis)
full_analysis <- subset(so_celltype, features = hvgs)
full_counts <- round(GetAssayData(full_analysis, layer = "counts")) 


meta_gene <- full_analysis@meta.data
pred_gene <- model.matrix(~group*sex, data = meta_gene)
data_g_gene <- list(count = full_analysis, id = meta_gene$record_id, pred = pred_gene)
result_allcells <- nebula(count = full_counts, id = full_analysis$record_id, pred = data_g_gene$pred, 
                          offset = full_analysis$pooled_offset,
                          ncore = 1, output_re = T, covariance = T,
                          reml = T, model = "NBLMM")

result_allcells$summary

write.table(result_allcells, paste(results.dir, 'NEBULA_fullanalysis_offset.txt'), row.names=F,
            quote=F, sep='\t')


#plotting_function(as.data.frame(result_allcells), 'Top2000HVGs')
if(status == 'T2D'){
sig_results <- as.data.frame(result_allcells) %>% 
  dplyr::select(Gene = summary.gene, LogFC = summary.logFC_groupType_2_Diabetes, 
                pvalue = summary.p_groupType_2_Diabetes) %>% filter(pvalue < 0.05)
}else if(status == 'T1D'){
  sig_results <- as.data.frame(result_allcells) %>% 
    dplyr::select(Gene = summary.gene, LogFC = summary.logFC_groupType_1_Diabetes, 
                  pvalue = summary.p_groupType_1_Diabetes) %>% filter(pvalue < 0.05)
}

}


scrna_analysis(status = 'T2D', Controls = 'LC', 
               results.dir = 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/', 
               celltype='PT')

scrna_analysis(status = 'T2D', Controls = 'LC', 
               results.dir = 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_All/', 
               celltype='All')

scrna_analysis(status = 'T2D', Controls = 'LC', 
               results.dir = 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_TAL/', 
               celltype='TAL')

scrna_analysis(status = 'T2D', Controls = 'LC', 
               results.dir = 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_EC/', 
               celltype='EC')






pdf("C:/Documents and Settings/netio/Documents/UofW/Projects/Sex_based_Analysis/CellDistributions_byGroup.pdf",
    width=15, height=12)
ggplot(meta.data, 
       aes(x=group, fill=KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+
  labs(x='Condition Group', fill='Cell Type')+
  theme(text = element_text(size = 15))

dev.off()




meta.data <- meta.data %>% 
  mutate(group_sex=paste0(group,'_',sex)) %>% 
  filter(group == 'Lean_Control' | group == 'Type_2_Diabetes')


pdf("C:/Documents and Settings/netio/Documents/UofW/Projects/Sex_based_Analysis/CellDistributions_bySex.pdf",
    width=15, height=12)
ggplot(meta.data, 
       aes(x=group_sex, fill=KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+
  labs(x='Sex and Diabetes Status', fill='Cell Type')+
  theme(text = element_text(size = 15))

dev.off()










#Plotting functions

plotting_function <- function(file, output_file_prefix, results.dir = 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/'){
  
  sig_markers <- data.table::fread(file)
  sig_markers <- sig_markers %>% dplyr::select(gene = summary.gene,
                                               logfc_sex = summary.logFC_sexMale, 
                                               logfc_group = summary.logFC_groupType_2_Diabetes, 
                                               logfc_interaction = summary.logFC_groupType_2_Diabetes.sexMale, 
                                               p_sex = summary.p_sexMale,
                                               p_group = summary.p_groupType_2_Diabetes, 
                                               p_interaction = summary.p_groupType_2_Diabetes.sexMale)
  
  variable_names <- c('Sex Effects', 'Diagnosis Group Effects', 'Interaction Effects')
  
  
  for(i in c(1:3)){
    tmp_logfc <- i + 1
    tmp_p <- i + 4
    tmp_df <- sig_markers %>% dplyr::select(gene, tmp_logfc, tmp_p)
    tmp_label <- variable_names[i]
    names(tmp_df) <- c('Gene', 'LogFC', 'Pvalue')
    tmp_df <- tmp_df %>% filter(!is.na(Pvalue))
    tmp_df$diffexp <- 'No'
    tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
    tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
    tmp_df <- tmp_df %>% arrange(Pvalue)
    tmp_df$label <- NA
    tmp_df$label[1:10] <- tmp_df$Gene[1:10]
    
    list_to_remove <- tmp_df %>% filter(Pvalue > 0.05 & abs(LogFC) > 8)
    
    if(nrow(list_to_remove) > 0){
      tmp_df <- tmp_df %>% filter(!Gene %in% list_to_remove$Gene)
    }
    
    #Making graph
    if(length(unique(tmp_df$diffexp)) > 1){
      tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
        geom_point()+
        geom_text(size=2, vjust = 2, color='black')+
        scale_color_manual(values = c('orange', 'grey', 'purple'),
                           labels = c('Downregulated', 'Not significant', 'Upregulated'))+
        geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
        geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
        theme_classic()+labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
                             title = paste0(tmp_label, ' in NEBULA Analysis'), 
                             subtitle = output_file_prefix)
    }else{
      tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
        geom_point()+
        geom_text(size=2, vjust = 2, color='black')+
        scale_color_manual(values = c('grey'),
                           labels = c('Not significant'))+
        geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
        geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
        theme_classic()+labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
                             title = paste0(tmp_label, ' in NEBULA Analysis'),
                             subtitle = output_file_prefix)
      
    }
    pdf(paste0(results.dir, 'NEBULA_VolcanoPlot_', output_file_prefix, '_', 
               tmp_label, '.pdf'))
    print(tmp_graph)
    dev.off()
  }
}






sig_genes <- data.table::fread("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/NEBULA_fullanalysis_offset.csv")

plotting_function("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/ NEBULA_fullanalysis_offset.txt",
                  'T2D_vs_LC_PT', results.dir = 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/')


so <- subset(so, features = sig_genes$gene)

meta <- so@meta.data
meta <- meta %>% filter(group %in% c("Type_2_Diabetes", 'Lean_Control'))

so <- subset(so, kit_id %in% meta$kit_id)

so@meta.data <- so@meta.data %>% 
  mutate(group_sex=paste0(group,'_',sex))


library(Seurat)
library(ggplot2)
library(patchwork)


so <- subset(so, celltype2 == 'PT')


for(gene in sig_genes_small$summary.gene) {
  p <- VlnPlot(so, 
               features = gene, 
               group.by = 'group_sex',
               pt.size = 0.1) + 
    ggtitle(paste("Expression of", gene))
  
  # Save each plot
  ggsave(filename = paste0("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/violin_", gene, ".png"), 
         plot = p, 
         width = 8, 
         height = 6)
}




DotPlot(so, features = sig_genes_small$summary.gene)



enrich_GO_BP <- enrichGO(gene = sig_genes_small$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP") %>% 
#  as.data.frame() %>% #dplyr::select(Description, GeneRatio, p.adjust, Count) %>% 
  filter(p.adjust < 0.05)

#dotplot
pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/GSEA_SignificantInteraction_dotplot.pdf',
    width=10, height=15)
dotplot(enrich_GO_BP, showCategory = 20)
dev.off()

# Bar plot
pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/GSEA_SignificantInteraction_barplot.pdf',
    width=10, height=15)
barplot(enrich_GO_BP, showCategory = 20)
dev.off()

# Network plot
pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/GSEA_SignificantInteraction_network.pdf',
    width=10, height=15)
cnetplot(enrich_GO_BP, categorySize = "pvalue", foldChange = NULL)
dev.off()
# Enrichment map
pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/GSEA_SignificantInteraction_map.pdf',
    width=10, height=15)
emapplot(enrich_GO_BP, showCategory = 30)
dev.off()
# Tree plot
pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_LC_PT/GSEA_SignificantInteraction_tree.pdf',
    width=10, height=15)
treeplot(enrich_GO_BP)
dev.off()


library(PANTHER.db)

pthOrganisms(PANTHER.db) <- "HUMAN"





GO_pathways <- function(data, results.dir, label){
  data <- data.table::fread(data)
  data <- data %>% as.data.frame() %>% 
    dplyr::select()
  
  
  
}



