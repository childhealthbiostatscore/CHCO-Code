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

<<<<<<< Updated upstream
#Compare effect esitmates
pdf(file = fs::path(dir.results,"Senesence Results","ComparisonPlots_ALT.pdf"),width=10,height=15)
plot(pvalue)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","ComparisonPlots_ALT2.pdf"),width=10,height=15)
plot(plot)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","ComparisonPlots_ALT3.pdf"),width=20,height=20)
plot(plot)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","ComparisonPlots_ALT4.pdf"),width=20,height=20)
plot(plot)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","Volcano_Diabetes.pdf"),width=15,height=10)
plot(p)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","Barchart_Diabetes.pdf"),width=8,height=5)
plot(b)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","Volcano_MASLD.pdf"),width=15,height=10)
plot(p)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","Barchart_MASLD.pdf"),width=8,height=5)
plot(b)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","Volcano_GLP_1.pdf"),width=15,height=10)
plot(p)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","Barchart_GLP_1.pdf"),width=8,height=5)
plot(b)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","ALT.pdf"),width=8,height=5)
plot(b)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","AST.pdf"),width=8,height=5)
plot(b)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","Fibrosis_Stage.pdf"),width=8,height=5)
plot(b)
dev.off()

pdf(file = fs::path(dir.results,"Senesence Results","Steatosis_Grade.pdf"),width=8,height=5)
=======
pdf(file = fs::path(dir.dat,"AST.pdf"),width=8,height=5)
>>>>>>> Stashed changes
plot(b)
dev.off()

#Diff exp by diabetes status
pdf(file = fs::path(dir.results,"Results_and_Figures","Volcano_byDiabetes_heptatocytes.pdf"),width=15,height=10)
de.markers(so_liver_sn_hep, genes, "diagnosis_of_diabetes", id2 = "No", id1 = "Yes", "Hepatocyte", "")
write.csv(m,fs::path(dir.code,"Differential_Expression_Diabetes.csv"))

m_top <- m
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
                     title = paste0("Differentially Expressed Genes by Diabetes Status (Pseudobulk All Cell Types)"),
                     subtitle = paste0("Positive Log2 FC = Greater Expression in Diabetes vs. No Diabetes\n",
                                       "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                     pCutoff = 0.05,
                     FCcutoff = 0.5,
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

#Print senesence figures
#Type 2
# Define custom colors based on significance and log2 fold change direction
colCustom <- ifelse(
  m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "red",           # Significant & positive FC
  ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "blue",    # Significant & negative FC
         "lightgray")                                               # Non-significant
)
colCustom[is.na(colCustom)] <- "lightgray"
names(colCustom)[colCustom =="red"] <- "Upregulated in Diabetes"
names(colCustom)[colCustom =="blue"] <- "Downregulated in Diabetes"
names(colCustom)[colCustom=="lightgray"] <- "Unchanged"

# Create the volcano plot
p <- EnhancedVolcano(m_top,
                     lab = labels,
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     title = paste0("Differentially Expressed Senescence Genes by Diabetes (Pseudobulk)"),
                     subtitle = paste0("Significant at FDR-P<0.05"),
                     pCutoff = 0.05,         # Set p-value cutoff for significance
                     FCcutoff = 0,           # No effect size threshold
                     labFace = 'bold',
                     pointSize = 4,
                     labSize = 5,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black',
                     legendPosition = NULL,
                     boxedLabels = TRUE,
                     max.overlaps = 30,
                     colCustom = colCustom)  # Apply custom colors
pdf(file = "Volcano_Senesence_Diabetes.pdf",width=15,height=8)
plot(p)
dev.off()

#MASLD
# Define custom colors based on significance and log2 fold change direction
colCustom <- ifelse(
  m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "red",           # Significant & positive FC
  ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "blue",    # Significant & negative FC
         "lightgray")                                               # Non-significant
)
colCustom[is.na(colCustom)] <- "lightgray"
names(colCustom)[colCustom =="red"] <- "Upregulated in MASLD"
names(colCustom)[colCustom =="blue"] <- "Downregulated in MASLD"
names(colCustom)[colCustom=="lightgray"] <- "Unchanged"

# Create the volcano plot
p <- EnhancedVolcano(m_top,
                     lab = labels,
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     title = paste0("Differentially Expressed Senescence Genes by MASLD (Pseudobulk)"),
                     subtitle = paste0("Significant at FDR-P<0.05"),
                     pCutoff = 0.05,         # Set p-value cutoff for significance
                     FCcutoff = 0,           # No effect size threshold
                     labFace = 'bold',
                     pointSize = 4,
                     labSize = 5,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black',
                     legendPosition = NULL,
                     boxedLabels = TRUE,
                     max.overlaps = 20,
                     colCustom = colCustom)  # Apply custom colors
pdf(file = "Volcano_Senesence_MASLD.pdf",width=15,height=8)
plot(p)
dev.off()

# Define custom colors based on significance and log2 fold change direction
colCustom <- ifelse(
  m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "red",           # Significant & positive FC
  ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "blue",    # Significant & negative FC
         "lightgray")                                               # Non-significant
)
colCustom[is.na(colCustom)] <- "lightgray"
names(colCustom)[colCustom =="red"] <- "Upregulated in GLP-1"
names(colCustom)[colCustom =="blue"] <- "Downregulated in GLP-1"
names(colCustom)[colCustom=="lightgray"] <- "Unchanged"
# Create the volcano plot
p <- EnhancedVolcano(m_top,
                     lab = labels,
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     title = paste0("Differentially Expressed Senescence Genes by GLP-1 Agonists (Pseudobulk)"),
                     subtitle = paste0("Significant at FDR-P<0.05"),
                     pCutoff = 0.05,         # Set p-value cutoff for significance
                     FCcutoff = 0,           # No effect size threshold
                     labFace = 'bold',
                     pointSize = 4,
                     labSize = 5,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black',
                     legendPosition = NULL,
                     boxedLabels = TRUE,
                     max.overlaps = 20,
                     colCustom = colCustom)  # Apply custom colors
pdf(file = "Volcano_Senesence_GLP1.pdf",width=15,height=8)
plot(p)
dev.off()

#MAST Results ----
m_top <- fcHurdle[,c("primerid","coef","fdr")]
m_top <- as.data.frame(m_top)
rownames(m_top) <- m_top$primerid
m_top <- m_top %>% 
  dplyr::select(-primerid)
significant_genes <- m_top %>% filter(fdr < 0.05)

# Select the top 10 positive and top 10 negative log2FC genes that are significant
top_genes <- rbind(
  significant_genes %>% arrange(desc(coef)) %>% head(40),  # Top 10 positive log2FC
  significant_genes %>% arrange(coef) %>% head(40)         # Top 10 negative log2FC
)

# Define custom colors based on significance and log2 fold change direction
colCustom <- ifelse(
  m_top$fdr < 0.05 & m_top$coef > 0, "red",           # Significant & positive FC
  ifelse(m_top$fdr < 0.05 & m_top$coef < 0, "blue",    # Significant & negative FC
         "lightgray")                                               # Non-significant
)
colCustom[is.na(colCustom)] <- "lightgray"
names(colCustom)[colCustom =="red"] <- "Upregulated"
names(colCustom)[colCustom =="blue"] <- "Downregulated"
names(colCustom)[colCustom=="lightgray"] <- "Unchanged"

labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
p <- EnhancedVolcano(m_top,
                     lab = labels,
                     x = 'coef',
                     y = 'fdr',
                     title = paste0("Fibrosis Stage"),
                     # subtitle = paste0("Positive Log2 FC = Greater Expression in Female vs. Male\n",
                     #                   "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                     # xlab = "Beta Coefficient",
                     ylab = "-Log10 Adjusted P-Value",
                     # Adjust x-axis scale
                     xlim = c(-0.25, 0.25),
                     pCutoff = 0.05,
                     FCcutoff = 0.005,
                     labFace = 'bold',
                     pointSize = 4,
                     labSize = 5,
                     drawConnectors = TRUE,
                     widthConnectors = 0.7,
                     colConnectors = 'black',
                     legendPosition=NULL,
                     boxedLabels = TRUE,
                     max.overlaps=50,
                     colCustom = colCustom)

plot(p)
  
pdf(file="Adj_Fibrosis_Grade_Plot.pdf",width=15,height=10)
plot(p)
dev.off()

<<<<<<< Updated upstream


table <- table1(~ age + nih_sex + nih_race + diagnosis_of_MASLD + diagnosis_of_diabetes + bmi+ tg + ast + alt + ggt + steatosis_percent + lobular_inflammation_percent + fibrosis_stage | glp1agonist, data = dat)
write.csv(table, fs::path(dir.results,"Table1_ByGLP1.csv"))

=======
pdf(file="Steatosis_Grade_Plot_Discrete.pdf",width=10,height=5)
plot(bubble_plot_D)
dev.off()

pdf(file=fs::path(dir.results,"UMAP_Cell_Types.pdf"))
DimPlot(so_liver_sn, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 4)+
  ggtitle("UMAP Plot with Cell Type Labels")
dev.off()


pdf(file = fs::path(dir.results,"Volcano_DEG_Liver_MASLT_All_Cell_Types.pdf"),width=20,height=20)
plot(p)
dev.off()

pdf(file=fs::path(dir.results,"GSEA_Top_20_Liver_Diabetes_All_Cell_Types.pdf"),width=10,height=5)
plot(enrich.p)
dev.off()
>>>>>>> Stashed changes







grid.draw(
  venn.diagram(
    x = list(
      "ACR" = acr_neg,
      "GLP1" = glp1_neg,
      "SGLT2" = sglt2_neg,
      "Type2" = type2_neg
    ),
    category.names = c("ACR", "GLP1", "SGLT2", "Type2"),
    filename = NULL
  )
)
