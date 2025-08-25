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





#insulin clamp integration


dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Insulin_Clamp/'


insulin_genes <- c("INSR", "IGF1R", "INSRR",
                     "IRS1", "IRS2", "IRS4", "GRB2", "SHC1",
                     "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2",
  "AKT1", "AKT2", "AKT3",
  "PDPK1", "PTEN", "INPP4B",
  "MTOR", "RPTOR", "RICTOR", "MLST8",
  "AKT1S1", "TSC1", "TSC2", "RHEB",
  "EIF4EBP1", "RPS6KB1", "RPS6", "EIF4E",
  "SLC2A1", "SLC2A2", "SLC2A4",
  "SLC5A1", "SLC5A2",
  "HK1", "HK2", "GCK",
  "PFKFB2", "PFKFB3", "PKM", "LDHA",
  "PCK1", "PCK2", "G6PC", "FBP1", "FBP2",
  "FASN", "ACACA", "ACACB", "SCD",
  "SREBF1", "SREBP2", "INSIG1", "INSIG2",
  "CPT1A", "CPT2", "ACOX1",
  "PPARGC1A", "PPARGC1B", "NRF1", "TFAM",
  "ATP5F1A", "COX4I1", "NDUFS1", "SDHA",
  "ATP5F1A", "COX4I1", "NDUFS1", "SDHA",
  "PTPN1", "PTPRF", "PTPN11",
  "SOCS1", "SOCS3", "SOCS7",
  "PRKAA1", "PRKAA2",
  "FOXO1", "FOXO3", "FOXO4")




#load in correct data, ensure FSOC accuracy. Only Type 2's 


load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')



so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

test <- so_subset@meta.data

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat2 <- dat %>% dplyr::select(record_id, epic_sglti2_1, airg, acprg) %>% 
  filter(!is.na(epic_sglti2_1))

#Fix up FSOC data
test$airg <- NULL
test$acprg<- NULL



test <- test %>% left_join(dat2, by='record_id')

so_subset@meta.data$airg <- test$airg
so_subset@meta.data$acprg <- test$acprg




demo_table <- dat %>% filter(visit == 'baseline') %>% 
  filter(record_id %in% so_subset@meta.data$record_id) %>% 
  filter(!is.na(airg))

demo_table$group2 <- ifelse(demo_table$epic_sglti2_1 == 'Yes', 'T2D-SGLT2', 'T2D-No SGLT2')


table1::table1(~age + sex + bmi + airg + acprg + epic_sglti2_1 | group2, data = demo_table)




so_subset$celltype1 <- case_when(grepl("PT-",so_subset$celltype_rpca)~"PT",
                                 grepl("TAL-",so_subset$celltype_rpca)~"TAL",
                                 grepl("EC-",so_subset$celltype_rpca)~"EC",
                                 grepl("POD",so_subset$celltype_rpca)~"POD",
                                 grepl("MAC",so_subset$celltype_rpca)~"MAC",
                                 grepl("MON",so_subset$celltype_rpca)~"MON",
                                 grepl("PC-",so_subset$celltype_rpca)~"PC",
                                 grepl("FIB",so_subset$celltype_rpca)~"FIB_MC_VSMC",
                                 grepl("DTL",so_subset$celltype_rpca)~"DTL",
                                 so_subset$celltype_rpca=="DCT"~"DCT",
                                 so_subset$celltype_rpca=="ATL"~"ATL",
                                 so_subset$celltype_rpca=="B"~"B",
                                 so_subset$celltype_rpca=="T"~"T")
so_subset$celltype1 <- as.character(so_subset$celltype1)

so_subset$KPMP_celltype2 <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2 <- ifelse(so_subset$KPMP_celltype=="aPT" | 
                                so_subset$KPMP_celltype=="PT-S1/S2" | 
                                so_subset$KPMP_celltype == "PT-S3","PT",
                              ifelse(grepl("TAL",so_subset$KPMP_celltype),"TAL",
                                     ifelse(grepl("EC-",so_subset$KPMP_celltype),"EC",so_subset$KPMP_celltype2)))


so_subset$DCT_celltype <- ifelse((so_subset$KPMP_celltype=="DCT" | 
                                    so_subset$KPMP_celltype=="dDCT"), "DCT","Non-DCT")




counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round

sce <- SingleCellExperiment(assays = list(counts = counts_path))
# sce <- computeSumFactors(sce)
sce <- computeSumFactors(sce)
# View size factors
# sizeFactors(sce)
# STEP 3: Calculate offset → log(size factors)
pooled_offset <- sizeFactors(sce)
so_subset@meta.data$pooled_offset <- pooled_offset




dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Insulin_Clamp/'
celltypes <- c('PT', 'TAL', 'EC')

for(i in c(1:length(celltypes))){
#TAL 
  celltype <- celltypes[i]
  
so_celltype <- subset(so_subset, celltype2 == celltype)

num_cells <- nrow(so_celltype)
test <- so_celltype@meta.data %>% filter(!is.na(airg))
num_part <- length(unique(test$record_id))


#gene lists


genes_list <- insulin_genes

so_celltype <- subset(so_celltype, features = genes_list)

counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
count_gene <- counts_path
meta_gene <- subset(so_celltype)@meta.data


complete_idx <- complete.cases(meta_gene$airg)
cat("Cells with complete data:", sum(complete_idx), "\n")

# Step 3: Filter all your data to only include cells with complete predictor data
meta_gene <- meta_gene[complete_idx, ]
count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells


# Step 4: Create prediction matrix from the complete data
pred_gene <- model.matrix(~airg, data = meta_gene)


# library <- meta_gene$library_size
library <- meta_gene$pooled_offset
data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)

if (is.null(data_g_gene)) {
  data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
}

#With offset
result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                 pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)



#Make dataframe of final results
full_results <- as.data.frame(result)
#Calculate number of genes filtered out for low expression 

# print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
full_results <- full_results %>%
  mutate(fdr=p.adjust(`summary.p_airg`,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$`summary.p_airg`, 1e-10))  # Avoid log(0)

full_results$num_cells <- num_cells
full_results$num_part <- num_part

write.table(full_results,paste0(dir.results,"NEBULA_", 
                                celltype, "_cells_AIRg_T2D_pooledoffset.csv"),
            row.names=F, quote=F, sep=',')





#ACPRg

num_cells <- nrow(so_celltype)
test <- so_celltype@meta.data %>% filter(!is.na(airg))
num_part <- length(unique(test$record_id))


counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
count_gene <- counts_path
meta_gene <- subset(so_celltype)@meta.data


complete_idx <- complete.cases(meta_gene$acprg)
cat("Cells with complete data:", sum(complete_idx), "\n")

# Step 3: Filter all your data to only include cells with complete predictor data
meta_gene <- meta_gene[complete_idx, ]
count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells


# Step 4: Create prediction matrix from the complete data
pred_gene <- model.matrix(~acprg, data = meta_gene)


# library <- meta_gene$library_size
library <- meta_gene$pooled_offset
data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)

if (is.null(data_g_gene)) {
  data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
}

#With offset
result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                 pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)



#Make dataframe of final results
full_results <- as.data.frame(result)
#Calculate number of genes filtered out for low expression 

# print(paste0(nebula_nonconverged_percent*100, "% failed to converge"))
full_results <- full_results %>%
  mutate(fdr=p.adjust(`summary.p_acprg`,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$`summary.p_acprg`, 1e-10))  # Avoid log(0)

full_results$num_cells <- num_cells
full_results$num_part <- num_part

write.table(full_results,paste0(dir.results,"NEBULA_", 
                                                   celltype, "_cells_ACPRg_T2D_pooledoffset.csv"),
            row.names=F, quote=F, sep=',')

}







#Plotting the results 

remove(list=ls())

celltypes <- c('PT', 'TAL', 'EC')

for(i in c(1:length(celltypes))){
  celltype <- celltypes[i]
  
}

names(full_results)[6] <- 'pvalue' 
full_results$color1 <- ifelse(full_results$fdr < 0.05, "lightcoral", "gray")
full_results$color2 <- ifelse(full_results$pvalue < 0.05, "lightcoral", "gray")

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
  y = reorder(gene, `logFC_groupType_2_Diabetes`),
  x = `logFC_groupType_2_Diabetes`,
  color = color1,
  size = abs(`logFC_groupType_2_Diabetes`)
)) +
  geom_point(alpha = 0.7) +
  scale_color_identity(name = 'Significance (FDR)', labels = c('> 0.05', '<0.05'), guide = 'legend')+
  scale_size(range = c(2, 6), name = "|LogFC|") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +  # Retains grid lines
  labs(
    title = "Differentially Expressed TCA Cycle Genes in All Cell Types",
    subtitle = "LC vs. T2D (No SGLT2), Unadjusted (Pooled Offset)",
    x = "Log Fold Change",
    y = "Gene",
    caption = paste0(
      "Participant Number: T2D: ", t2d_count, ', LC: ', lc_count, 
      "; Genes = ", Genes,
      ", Cells = ", Cell
    )
  ) +
  theme(plot.caption = element_text(size = 8), 
        plot.title = element_text(hjust = 0),
        axis.text.y = element_text(size = 8),
        # axis.text.x = element_text(angle = 0, hjust = 1),
        axis.line = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()
  )
dot_plot

png(fs::path(dir.results, "fdr/Plot__TCA_cycle_NEBULA_All_Cells_T2D_LC_NoMed_unadjusted_pooled_offset_no_IT_08.png"), 
    width = 2500, height = 2000, res = 300)
print(dot_plot)
dev.off()




















