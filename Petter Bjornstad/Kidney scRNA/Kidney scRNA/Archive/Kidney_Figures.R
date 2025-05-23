#Figures
# #By sex
# jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap_by_Sex.jpeg"))
# so_kidney_sc <- RunUMAP(so_kidney_sc, dims = 1:15)
# DimPlot(so_kidney_sc, reduction = "umap", group.by = "sex") 
# dev.off()
# #By sglti2
# jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap_by_SGLT2i.jpeg"))
# DimPlot(so_kidney_sc, reduction = "umap", group.by = "sglt2i_ever") 
# dev.off()
# #by disease status
# jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap_by_disease_status.jpeg"))
# DimPlot(so_kidney_sc, reduction = "umap", group.by = "group") 
# dev.off()
# #General 
# jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap.jpeg"))
# DimPlot(so_kidney_sc, reduction = "umap")
# dev.off()

#Diff exp by sex
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySex_allcells_all_disease.pdf"),width=15,height=10)
de.markers(so_kidney_sc, NULL, "sex", id1 = "Female", id2 = "Male", NULL, "_top")
m_top <- m_top %>% head(2000)
significant_genes <- m_top %>% filter(p_val_adj < 0.05)

# Select the top 10 positive and top 10 negative log2FC genes that are significant
top_genes <- rbind(
  significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
  significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
)

labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
p <- EnhancedVolcano(m_top,
                     lab = labels,
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     title = paste0("Differentially Expressed Genes by Sex (Pseudobulk)"),
                     subtitle = paste0("Positive Log2 FC = Greater Expression in Female vs. Male\n",
                                       "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                     pCutoff = 0.05,
                     FCcutoff = 1.4,
                     labFace = 'bold',
                     pointSize = 4,
                     labSize = 5,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black',
                     legendPosition=NULL,
                     boxedLabels = TRUE,
                     max.overlaps=20)
plot(p)
dev.off()

#SGLT2is
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySGLT2i_allcells_all_disease.pdf"),width=15,height=10)
de.markers(so_kidney_sc, NULL, "sglt2i_ever", id1 = "Yes", id2 = "No", NULL, "_top")
m_top <- m_top %>% head(2000)
significant_genes <- m_top %>% filter(p_val_adj < 0.05)

# Select the top 10 positive and top 10 negative log2FC genes that are significant
top_genes <- rbind(
  significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
  significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
)

labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
p <- EnhancedVolcano(m_top,
                lab = labels,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Differentially Expressed Genes by SGLT2i (Pseudobulk)',
                subtitle = paste0("Positive Log2 FC = Greater Expression in SGLT2i Yes vs. No\n",
                "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                # pCutoff = 0.05,
                # FCcutoff = 0,
                pointSize = 1.5,
                pCutoff = 0.05,
                FCcutoff = 1.4,
                labFace = 'bold',
                labSize = 5,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                legendPosition=NULL,
                boxedLabels = TRUE,
                max.overlaps=20)
plot(p)
dev.off()

#By Cell Types
#Sex
cell_types <- c("MON","MAC","EC-GC","EC-AEA","EC-PTC")
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySex_byCelltypes.pdf"),width=15,height=8)
for (cell in cell_types){
  de.markers(so_kidney_sc, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  m_top <- m_top %>% head(2000)
  top_genes <- rbind(
    m_top %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    m_top %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  # Filter m_top for significant genes with adj p-value < 0.05
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)

  # Select the top 10 positive and top 10 negative log2FC genes that are significant
  top_genes <- rbind(
    significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = paste0("Differentially Expressed Genes by Sex in ",cell," cells"),
                       subtitle = paste0("Positive Log2 FC = Greater Expression in Female vs. Male\n",
                       "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                       pCutoff = 0.05,
                       FCcutoff = 1.4,
                       labFace = 'bold',
                       pointSize = 4,
                       labSize = 5,
                       drawConnectors = TRUE,
                       widthConnectors = 1.0,
                       colConnectors = 'black',
                       legendPosition=NULL,
                       boxedLabels = TRUE,
                       max.overlaps=20)
  
  plot(p)
}
dev.off()


#SGLT2is
cell_types <- c("MON","MAC","EC-GC","EC-AEA","EC-PTC")
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySGLT2i_byCelltypes.pdf"),width=15,height=8)
for (cell in cell_types){
  de.markers(so_kidney_sc, genes, "sglt2i_ever", id1 = "Yes", id2 = "No", cell, "_top")
  m_top <- m_top %>% head(2000)
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and top 10 negative log2FC genes that are significant
  top_genes <- rbind(
    significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = paste0("Differentially Expressed Genes by SGLT2i in ",cell," cells"),
                       subtitle = paste0("Positive Log2 FC = Greater Expression in SGLT2i Yes vs. No\n",
                                         "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                       pCutoff = 0.05,
                       FCcutoff = 1.4,
                       labFace = 'bold',
                       pointSize = 4,
                       labSize = 5,
                       drawConnectors = TRUE,
                       widthConnectors = 1.0,
                       colConnectors = 'black',
                       legendPosition=NULL,
                       boxedLabels = TRUE,
                       max.overlaps=20)

  plot(p)
}
dev.off()


#By Disease Status
#Sex
#HC
cell_types <- c("MON","MAC","EC-GC","EC-AEA","EC-PTC")
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySex_byCelltypes_HC.pdf"),width=15,height=8)
for (cell in cell_types){
  de.markers(so_HC, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  m_top <- m_top %>% head(2000)
  top_genes <- rbind(
    m_top %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    m_top %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  # Filter m_top for significant genes with adj p-value < 0.05
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and top 10 negative log2FC genes that are significant
  top_genes <- rbind(
    significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = paste0("Differentially Expressed Genes in Healthy Controls by Sex in ",cell," cells"),
                       subtitle = paste0("Positive Log2 FC = Greater Expression in Female vs. Male\n",
                                         "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                       pCutoff = 0.05,
                       FCcutoff = 1.4,
                       labFace = 'bold',
                       pointSize = 4,
                       labSize = 5,
                       drawConnectors = TRUE,
                       widthConnectors = 1.0,
                       colConnectors = 'black',
                       legendPosition=NULL,
                       boxedLabels = TRUE,
                       max.overlaps=20)
  
  plot(p)
}
dev.off()

#OB
cell_types <- c("MON","MAC","EC-GC","EC-AEA","EC-PTC")
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySex_byCelltypes_OB.pdf"),width=15,height=8)
for (cell in cell_types){
  de.markers(so_OB, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  m_top <- m_top %>% head(2000)
  top_genes <- rbind(
    m_top %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    m_top %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  # Filter m_top for significant genes with adj p-value < 0.05
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and top 10 negative log2FC genes that are significant
  top_genes <- rbind(
    significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = paste0("Differentially Expressed Genes in Obese Controls by Sex in ",cell," cells"),
                       subtitle = paste0("Positive Log2 FC = Greater Expression in Female vs. Male\n",
                                         "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                       pCutoff = 0.05,
                       FCcutoff = 1.4,
                       labFace = 'bold',
                       pointSize = 4,
                       labSize = 5,
                       drawConnectors = TRUE,
                       widthConnectors = 1.0,
                       colConnectors = 'black',
                       legendPosition=NULL,
                       boxedLabels = TRUE,
                       max.overlaps=20)
  
  plot(p)
}
dev.off()

#T1
cell_types <- c("MON","MAC","EC-GC","EC-AEA","EC-PTC")
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySex_byCelltypes_T1.pdf"),width=15,height=8)
for (cell in cell_types){
  de.markers(so_T1, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  m_top <- m_top %>% head(2000)
  top_genes <- rbind(
    m_top %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    m_top %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  # Filter m_top for significant genes with adj p-value < 0.05
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and top 10 negative log2FC genes that are significant
  top_genes <- rbind(
    significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = paste0("Differentially Expressed Genes in Type 1 by Sex in ",cell," cells"),
                       subtitle = paste0("Positive Log2 FC = Greater Expression in Female vs. Male\n",
                                         "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                       pCutoff = 0.05,
                       FCcutoff = 1.4,
                       labFace = 'bold',
                       pointSize = 4,
                       labSize = 5,
                       drawConnectors = TRUE,
                       widthConnectors = 1.0,
                       colConnectors = 'black',
                       legendPosition=NULL,
                       boxedLabels = TRUE,
                       max.overlaps=20)
  
  plot(p)
}
dev.off()

#T2
#Sex
cell_types <- c("MON","MAC","EC-GC","EC-AEA","EC-PTC")
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySex_byCelltypes_T2.pdf"),width=15,height=8)
for (cell in cell_types){
  de.markers(so_T2, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  m_top <- m_top %>% head(2000)
  top_genes <- rbind(
    m_top %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    m_top %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  # Filter m_top for significant genes with adj p-value < 0.05
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and top 10 negative log2FC genes that are significant
  top_genes <- rbind(
    significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = paste0("Differentially Expressed Genes in Type 2 by Sex in ",cell," cells"),
                       subtitle = paste0("Positive Log2 FC = Greater Expression in Female vs. Male\n",
                                         "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                       pCutoff = 0.05,
                       FCcutoff = 1.4,
                       labFace = 'bold',
                       pointSize = 4,
                       labSize = 5,
                       drawConnectors = TRUE,
                       widthConnectors = 1.0,
                       colConnectors = 'black',
                       legendPosition=NULL,
                       boxedLabels = TRUE,
                       max.overlaps=20)
  
  plot(p)
}
dev.off()

#SGLT2is
cell_types <- c("MON","MAC","EC-GC","EC-AEA","EC-PTC")
pdf(file = fs::path(dir.home,"Results_and_Figures","Volcano_bySGLT2i_byCelltypes_T2.pdf"),width=15,height=8)
for (cell in cell_types){
  de.markers(so_T2, genes, "sglt2i_ever", id1 = "Yes", id2 = "No", cell, "_top")
  m_top <- m_top %>% head(2000)
  top_genes <- rbind(
    m_top %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    m_top %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  # Filter m_top for significant genes with adj p-value < 0.05
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and top 10 negative log2FC genes that are significant
  top_genes <- rbind(
    significant_genes %>% arrange(desc(avg_log2FC)) %>% head(10),  # Top 10 positive log2FC
    significant_genes %>% arrange(avg_log2FC) %>% head(10)         # Top 10 negative log2FC
  )
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = paste0("Differentially Expressed Genes in Type 2 by SGLT2i in ",cell," cells"),
                       subtitle = paste0("Positive Log2 FC = Greater Expression in SGLT2i Yes vs. No\n",
                                         "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),                       pCutoff = 0.05,
                       FCcutoff = 1.4,
                       labFace = 'bold',
                       pointSize = 4,
                       labSize = 5,
                       drawConnectors = TRUE,
                       widthConnectors = 1.0,
                       colConnectors = 'black',
                       legendPosition=NULL,
                       boxedLabels = TRUE,
                       max.overlaps=20)
  
  plot(p)
}
dev.off()

#Save all results to excel 
de.markers(so_kidney_sc, NULL, "sex", id1 = "Female", id2 = "Male", NULL, "_top")
r1 <- m_top 
de.markers(so_kidney_sc, NULL, "sglt2i_ever", id1 = "Yes", id2 = "No", NULL, "_top")
r2 <- m_top

wb <- createWorkbook()
addWorksheet(wb,"By_Sex")
writeData(wb,"By_Sex",r1,rowNames = T)
addWorksheet(wb,"By_SGLT2i")
writeData(wb,"By_SGLT2i",r2,rowNames = T)
cell_types <- c("MON","MAC","EC-GC","EC-AEA","EC-PTC")
#By cell for sex
for (cell in cell_types){
  de.markers(so_kidney_sc, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  r <- m_top 
  name <- paste0("By_Sex_in_",cell,"_cells")
  addWorksheet(wb,name)
  writeData(wb,name,r,rowNames = T)
}
#By cell for SGLT2is
for (cell in cell_types){
  de.markers(so_kidney_sc, genes, "sglt2i_ever", id1 = "Yes", id2 = "No", cell, "_top")
  r <- m_top 
  name <- paste0("By_SGLT2i_in_",cell,"_cells")
  addWorksheet(wb,name)
  writeData(wb,name,r,rowNames = T)
}
#HC by cell for sex
for (cell in cell_types){
  de.markers(so_HC, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  r <- m_top 
  name <- paste0("HC_By_Sex_in_",cell,"_cells")
  addWorksheet(wb,name)
  writeData(wb,name,r,rowNames = T)
}
#OB by cell for sex
for (cell in cell_types){
  de.markers(so_OB, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  r <- m_top 
  name <- paste0("OB_By_Sex_in_",cell,"_cells")
  addWorksheet(wb,name)
  writeData(wb,name,r,rowNames = T)
}
#T1 by cell for sex
for (cell in cell_types){
  de.markers(so_T1, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  r <- m_top 
  name <- paste0("T1_By_Sex_in_",cell,"_cells")
  addWorksheet(wb,name)
  writeData(wb,name,r,rowNames = T)
}
#T2 by cell for sex
for (cell in cell_types){
  de.markers(so_T2, genes, "sex", id1 = "Female", id2 = "Male", cell, "_top")
  r <- m_top 
  name <- paste0("T2_By_Sex_in_",cell,"_cells")
  addWorksheet(wb,name)
  writeData(wb,name,r,rowNames = T)
}
#T2 by cell for SGLT2i
for (cell in cell_types){
  de.markers(so_T2, genes, "sglt2i_ever", id1 = "Yes", id2 = "No", cell, "_top")
  r <- m_top 
  name <- paste0("T2_By_SGLT2i_in_",cell,"_cells")
  addWorksheet(wb,name)
  writeData(wb,name,r,rowNames = T)
}

saveWorkbook(wb,fs::path(dir.home,"Results_and_Figures","Kidney_scRNA_Diff_Exp.xlsx"),overwrite = TRUE)



# Assuming your code prior to this point remains the same
library(stringr)

#HC
kidney_meta_raw_sc <- so_HC@meta.data %>%
  filter(cryostor_id != "")
form <- paste("sex", paste(colnames(kidney_meta_raw_sc)[c(18:27)], collapse = " + "), sep = " ~ ")

# Perform the tableby analysis
table_summary <- arsenal::tableby(formula = as.formula(form), data = kidney_meta_raw_sc, test = TRUE)

# Convert the summary to a data frame
summary_df <- as.data.frame(summary(table_summary))
colnames(summary_df)[1] <- "Variable"
summary_df$Variable <- str_remove_all(summary_df$Variable,"\\**")
summary_df$Variable <- str_remove_all(summary_df$Variable,"\\&nbsp;&nbsp;&nbsp;")
summary_df$Variable <- str_to_title(summary_df$Variable)

# Save the summary to a CSV file
write.csv(summary_df, "HC_Descriptives_by_Sex.csv", row.names = FALSE)

#T2
kidney_meta_raw_sc <- so_T2@meta.data %>%
  filter(cryostor_id != "")
form <- paste("sex", paste(colnames(kidney_meta_raw_sc)[c(18:27)], collapse = " + "), sep = " ~ ")

# Perform the tableby analysis
table_summary <- arsenal::tableby(formula = as.formula(form), data = kidney_meta_raw_sc, test = TRUE)

# Convert the summary to a data frame
summary_df <- as.data.frame(summary(table_summary))
colnames(summary_df)[1] <- "Variable"
summary_df$Variable <- str_remove_all(summary_df$Variable,"\\**")
summary_df$Variable <- str_remove_all(summary_df$Variable,"\\&nbsp;&nbsp;&nbsp;")
summary_df$Variable <- str_to_title(summary_df$Variable)

# Save the summary to a CSV file
write.csv(summary_df, "T2_Descriptives_by_Sex.csv", row.names = FALSE)
