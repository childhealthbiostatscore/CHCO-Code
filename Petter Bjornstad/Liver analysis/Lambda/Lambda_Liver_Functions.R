#Define gene sets -----
sens_genes <- c("ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2","BMP6","C3","CCL1","CCL13",
                "CCL16","CCL2","CCL20","CCL24","CCL26","CCL3","CCL4","CCL5","CCL7","CCL8","CD55",
                "CD9","CSF1","CSF2","CSF2RB","CST4","CTNNB1","CTSB","CXCL1","CXCL10","CXCL12",
                "CXCL16","CXCL2","CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR","EREG","ESM1",
                "ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG","HGF","HMGB1","ICAM1","ICAM3",
                "IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13",
                "IL15","IL18","IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7","INHA","IQGAP2","ITGA2",
                "ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13","MMP14","MMP2",
                "MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT","PLAU","PLAUR",
                "PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG","SEMA3F","SERPINB4","SERPINE1",
                "SERPINE2","SPP1","SPX","TIMP2","TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B",
                "TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2")
#sens_genes2 <- c("GLB1","TP53", "CDKN2A", "RB1", "CDK4", "CDK6",
#                  "IL6", "IL8", "CXCL1", "CXCL10", "CXCL12", "ICAM1", "MMP1", "MMP2", "MMP3", "TIMP1",
#                  "H2AFX", "ATM", "ATR", "CHEK1", "CHEK2", "TP53BP1", "MDC1",
#                  "H3K9me3", "CBX1", "CBX3", "H2AFY","NOX4", "NOX5", "SOD1", "SOD2", "GPX1", 
#                  "GPX4", "NFE2L2", "PRDX1", "PRDX3", "PRDX6","BMI1", "EZH2","TGFBR2", "FN1", 
#                  "LGALS1", "SERPINE1", "CDKN2B", "IGFBP3", "IGFBP7", "LMNB1")
# Find the intersection
#full_sens_genes <- union(sens_genes, sens_genes2)
#sens_genes <- full_sens_genes

# # #Metabolism Pathway Gene Sets
# tca<- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","tca_cycle_pathway.csv")) %>%
#   pull(genesymbol)
# # gluco <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","gluconeogenesis_pathway.csv"))%>%
# #   pull(genesymbol)
# glyc <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","glycolysis_pathway.csv"))%>%
#   pull(genesymbol)
# beta_ox <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","beta_ox_pathway.csv"))%>%
#   pull(genesymbol)
# ox_phos <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","ox_phos_pathway.csv"))%>%
#   pull(genesymbol)
# fa_met <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","fa_met_pathway.csv"))%>%
#   pull(genesymbol)
# energy_met <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","energy_metabolism_pathway.csv"))%>%
#   pull(genesymbol)
# met <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","metabolism_pathway.csv"))%>%
#   pull(genesymbol)
# 
# #Inflammation Pathway Gene Sets
# cytokines <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","cytokines_inflammatory_response_pathway.csv"))%>%
#   pull(genesymbol)
# inflamatory <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","inflamatory_response_pathway.csv"))%>%
#   pull(genesymbol)
# chemokines <- read.csv(fs::path(dir.dat,"Liver project","Pathway_Genes","infl_chemokine_cytokine_pathway.csv"))%>%
#   pull(genesymbol)

#DEGS ----
# de.markers(so_liver_sn, genes, "diagnosis_of_diabetes", id2 = "No", id1 = "Yes", NULL, "_top")
# seurat_object <- so_kidney_sc
# group.by <- "epic_glp1ra_1"
# id1 <- "No"
# id2 <- "Yes"
# celltype="PT"
# genes <- sens_genes
# extension="_top"
# # function for de.markers
de.markers <- function(seurat_object, genes, group.by, id1, id2, celltype, extension,logfc.threshold,min.pct){
  m = FindMarkers(seurat_object, features = genes,group.by = group.by,ident.1 = id1, 
                  ident.2 = id2, subset.ident = celltype, verbose = F, logfc.threshold=0.001,
                  min.pct = 0.001)
  m$p_val_adj = p.adjust(m$p_val,method = "fdr")
  m <- m %>% 
    rownames_to_column('gene') %>%
    arrange(p_val) %>%
    column_to_rownames('gene') %>%
    dplyr::select(avg_log2FC,pct.1,pct.2,p_val,p_val_adj) %>%
    filter(!is.na(p_val))
  
  genes_subset <- rownames(m)[m$p_val <= 0.05]
  
  if (length(genes_subset) > 0){
    assign(paste0("genes_subset", extension), genes_subset, envir = .GlobalEnv)
  }
  assign(paste0("m", extension), m, envir = .GlobalEnv)
  return(knitr::kable(m, digits = 3
  ))
}

#DEG & GSEA Function ----
##a. All Genes ----
degs_fxn <- function(so,cell,exposure,gene_set,exp_group,ref_group,enrichment,top_gsea) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  condition <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group," vs. ",ref_group,")")
  
  #Differential Expression by Group
  if (!is.null(cell)) {
    Idents(so) <- so$celltype2
    cell_name <- str_replace_all(cell,"/","_")
  }
  # sens_genes <- c(sens_genes,"CDKN1A")
  # de.markers(so, gene_set, "group", id2 = "neither", id1 = "both", "PT", "")
  de.markers(so, gene_set, exposure, id2 = ref_group, id1 = exp_group, cell, "")
  colnames(m)[2] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group,")")
  colnames(m)[3] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",ref_group,")")
  
  # Filter for significant genes
  m_top <- m
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and 10 negative log2FC genes based on the largest magnitude of fold change
  top_positive_by_fc <- significant_genes %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
    head(10)  # Top 10 positive fold changes
  
  top_negative_by_fc <- significant_genes %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
    head(10)  # Top 10 negative fold changes
  
  # Select the top 10 positive and 10 negative log2FC genes based on significance (p_val_adj)
  top_positive_by_significance <- significant_genes %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
    head(10)  # Top 10 most significant positive fold changes
  
  top_negative_by_significance <- significant_genes %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
    head(10)  # Top 10 most significant negative fold changes
  
  # Combine top fold-change based and significance-based genes into a final list
  top_genes <- rbind(
    top_positive_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
    top_negative_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
    top_positive_by_significance %>% mutate(Selection = "Top 10 by Significance"),
    top_negative_by_significance %>% mutate(Selection = "Top 10 by Significance")
  )
  
  
  if (!is.null(cell)){
    title <- paste0("DEGs in ",cell_name," cells for ",condition)
  } else {
    title <- paste0("Bulk DEGs for ",condition)
  }
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       subtitle = paste0("Positive Log2 FC = Greater Expression in ", condition,"\n",
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
                       max.overlaps=60)
  if (!is.null(cell)){
    filename <- paste0("DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
  } else {
    filename <- paste0("Bulk_DEGs_for_",condition,".pdf") 
  }
  pdf(fs::path(dir.results,filename),width=20,height=10)
  plot(p)
  dev.off()
  
  #GSEA
  if (enrichment=="Yes") {
    #Gene set enrichment analysis
    gc()
    sce_sn_hep <- as.SingleCellExperiment(so)
    ## C2 category is according to canonical pathways: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
    geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
    ### filter background to only include genes that we assessed
    geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
    geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_sn_hep),]
    m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
    stats <- m$p_val_adj
    names(stats) <- rownames(m)
    eaRes <- fgsea(pathways = m_list, stats = na.omit(stats))
    #ooEA <- order(eaRes$pval, decreasing = FALSE)
    #kable(head(eaRes[ooEA, 1:7], n = 20))
    # Convert the leadingEdge column to comma-separated strings
    eaRes$leadingEdge <- sapply(eaRes$leadingEdge, function(x) paste(x, collapse = ", "))
    gc()
    #Significant pathways
    sig <- eaRes %>% 
      filter(padj<0.05)
    
    #Plot top pathways
    # Subset top pathways for visualization
    top_pathways <- eaRes[order(-abs(eaRes$NES)), ][1:top_gsea, ]  # Top 10 pathways based on adjusted p-value
    
    # Define a significance threshold
    significance_threshold <- 0.05
    
    # Add a significance column for coloring
    top_pathways$significance <- ifelse(
      top_pathways$padj > significance_threshold, "Non-significant",
      ifelse(top_pathways$NES > 0, "Positive Significant", "Negative Significant")
    )
    
    if (!is.null(cell)) {
      title <- paste0("Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition)
    } else {
      title <- paste0("Top Enriched Pathways by NES (GSEA) for ",condition," (Bulk)")
    }
    # Create a bar plot
    gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
      geom_bar(stat = "identity") +
      coord_flip() +  # Flip coordinates for better readability
      scale_fill_manual(
        values = c(
          "Positive Significant" = "red",
          "Negative Significant" = "blue",
          "Non-significant" = "gray"
        ),
        name = "Significance"
      ) +
      labs(
        title = title,
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")
      )
    
    
    if (!is.null(cell)){
      filename <- paste0("GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
    } else {
      filename <- paste0("Bulk_GSEA_for_",condition,".pdf") 
    }
    pdf(fs::path(dir.results,filename),width=20,height=20)
    plot(gsea_plot)
    dev.off()
    
  }
  #Make sure gene nemaes are in printed file
  deg_results <- m 
  deg_results$Gene <- rownames(deg_results)
  
  #save results to excel file
  write_multiple_sheets <- function(output_file, sheet_data) {
    # Create a new workbook
    wb <- createWorkbook()
    
    # Loop over the list of data frames
    for (sheet_name in names(sheet_data)) {
      # Add a new sheet to the workbook
      addWorksheet(wb, sheet_name)
      
      # Write the data frame to the sheet
      writeData(wb, sheet_name, sheet_data[[sheet_name]])
    }
    # Save the workbook to the specified file
    saveWorkbook(wb, file = output_file, overwrite = TRUE)
    
    # Print a message
    message("Excel file with multiple sheets created: ", output_file)
  } 
  # Example usage
  df1 <- deg_results
  
  if (enrichment=="Yes") {
    df2 <- eaRes
  }
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".xlsx"))  
  } else {
    output_file <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".xlsx"))
  }
  
  if (enrichment=="Yes") {
    sheet_data <- list(
      "DEG" = df1,
      "Pathway_Results" = df2
    )
  } else {
    sheet_data <- list(
      "DEG" = df1
    ) 
  }
  
  # Call the function
  write_multiple_sheets(output_file, sheet_data)
  
  return(p)
}

##b. Specific Pathways----
degs_fxn_pathway <- function(so,cell,exposure,gene_set,exp_group,ref_group,pathway) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  condition <- paste0(pathway,str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group," vs. ",ref_group,")")
  
  #Differential Expression by Group
  if (!is.null(cell)) {
    Idents(so) <- so$celltype2
    cell_name <- str_replace_all(cell,"/","_")
  }
  # sens_genes <- c(sens_genes,"CDKN1A")
  # de.markers(so, gene_set, "group", id2 = "neither", id1 = "both", "PT", "")
  de.markers(so, gene_set, exposure, id2 = ref_group, id1 = exp_group, cell, "")
  colnames(m)[2] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group,")")
  colnames(m)[3] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",ref_group,")")
  
  # Filter for significant genes
  m_top <- m
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and 10 negative log2FC genes based on the largest magnitude of fold change
  top_positive_by_fc <- significant_genes %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
    head(10)  # Top 10 positive fold changes
  
  top_negative_by_fc <- significant_genes %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
    head(10)  # Top 10 negative fold changes
  
  # Select the top 10 positive and 10 negative log2FC genes based on significance (p_val_adj)
  top_positive_by_significance <- significant_genes %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
    head(10)  # Top 10 most significant positive fold changes
  
  top_negative_by_significance <- significant_genes %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
    head(10)  # Top 10 most significant negative fold changes
  
  # Combine top fold-change based and significance-based genes into a final list
  top_genes <- rbind(
    top_positive_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
    top_negative_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
    top_positive_by_significance %>% mutate(Selection = "Top 10 by Significance"),
    top_negative_by_significance %>% mutate(Selection = "Top 10 by Significance")
  )
  
  
  if (!is.null(cell)){
    title <- paste0(str_to_title(pathway)," DEGs in ",cell_name," cells for ",condition)
  } else {
    title <- paste0(str_to_title(pathway)," Bulk DEGs for ",condition)
  }
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       subtitle = paste0("Positive Log2 FC = Greater Expression in ", condition,"\n",
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
                       max.overlaps=60)
  if (!is.null(cell)){
    filename <- paste0("DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
  } else {
    filename <- paste0("Bulk_DEGs_for_",condition,".pdf") 
  }
  pdf(fs::path(dir.results,filename),width=20,height=10)
  plot(p)
  dev.off()
  
  # Filter for significant genes
  m_top <- m
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and 10 negative log2FC genes based on the largest magnitude of fold change
  top_positive_by_fc <- significant_genes %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
    head(10)  # Top 10 positive fold changes
  
  top_negative_by_fc <- significant_genes %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
    head(10)  # Top 10 negative fold changes
  
  # Select the top 10 positive and 10 negative log2FC genes based on significance (p_val_adj)
  top_positive_by_significance <- significant_genes %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
    head(10)  # Top 10 most significant positive fold changes
  
  top_negative_by_significance <- significant_genes %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
    head(10)  # Top 10 most significant negative fold changes
  
  # Combine top fold-change based and significance-based genes into a final list
  top_genes <- rbind(
    top_positive_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
    top_negative_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
    top_positive_by_significance %>% mutate(Selection = "Top 10 by Significance"),
    top_negative_by_significance %>% mutate(Selection = "Top 10 by Significance")
  )
  if (!is.null(cell)){
    title <- paste0(str_to_title(pathway)," DEGs in ",cell_name," cells for ",condition)
  } else {
    title <- paste0(str_to_title(pathway)," Bulk DEGs for ",condition)
  }
  
  # Add a column to classify points based on significance and direction
  m_top$Significance <- ifelse(
    m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "Significant Positive",
    ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "Significant Negative", "Non-Significant")
  )
  
  # Add labels only for significant points
  m_top$label <- ifelse(m_top$Significance != "Non-Significant", rownames(m_top), NA)
  
  # Define custom colors
  colors <- c("Significant Positive" = "red", 
              "Significant Negative" = "blue", 
              "Non-Significant" = "gray")
  
  # Create the ggplot volcano plot
  p2 <- ggplot(m_top, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
    geom_point(size = 4, alpha = 0.8) +  # Points with size and transparency
    scale_color_manual(values = colors) +  # Set custom colors
    labs(
      title = title,
      subtitle = paste0("Positive Log2 FC = Greater Expression in ", condition, "\n",
                        "(Significant at FDR-P < 0.05)"),
      x = "Log2 Fold Change (avg_log2FC)",
      y = "-Log10 Adjusted P-Value (p_val_adj)"
    ) +
    theme_minimal() +  # Clean plot theme
    theme(
      text = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "top"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance threshold line
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Center line for FC
    geom_text_repel(
      aes(label = label),  # Label significant points
      color = "black",     # Ensure labels are in black
      size = 4, 
      max.overlaps = 10,  # Adjust this value to control label density
      box.padding = 0.5,
      point.padding = 0.3
    )
  
  if (!is.null(cell)){
    filename <- paste0(str_to_title(pathway),"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
  } else {
    filename <- paste0(str_to_title(pathway),"_Bulk_DEGs_for_",condition,".pdf") 
  }
  pdf(fs::path(dir.results,filename),width=10,height=7)
  plot(p2)
  dev.off()
  
  #Make sure gene nemaes are in printed file
  deg_results <- m 
  deg_results$Gene <- rownames(deg_results)
  
  #save results to excel file
  write_multiple_sheets <- function(output_file, sheet_data) {
    # Create a new workbook
    wb <- createWorkbook()
    
    # Loop over the list of data frames
    for (sheet_name in names(sheet_data)) {
      # Add a new sheet to the workbook
      addWorksheet(wb, sheet_name)
      
      # Write the data frame to the sheet
      writeData(wb, sheet_name, sheet_data[[sheet_name]])
    }
    # Save the workbook to the specified file
    saveWorkbook(wb, file = output_file, overwrite = TRUE)
    
    # Print a message
    message("Excel file with multiple sheets created: ", output_file)
  } 
  # Example usage
  df1 <- deg_results
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0(str_to_title(pathway),"_Results_",cell_name,"_cells_for_",condition,".xlsx"))  
  } else {
    output_file <- fs::path(dir.results,paste0(str_to_title(pathway),"_Bulk_Results_for_",condition,".xlsx"))
  }
  
  sheet_data <- list(
    "DEG" = df1
  )
  # Call the function
  write_multiple_sheets(output_file, sheet_data)
  return(p)
}



##c. Senescence Genes ----
degs_fxn_sens <- function(so,cell,exposure,gene_set,exp_group,ref_group,enrichment,top_gsea) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  condition <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group," vs. ",ref_group,")")
  
  #Differential Expression by Group
  if (!is.null(cell)) {
    Idents(so) <- so$celltype2
    cell_name <- str_replace_all(cell,"/","_")
  }
  # sens_genes <- c(sens_genes,"CDKN1A")
  # de.markers(so, gene_set, "group", id2 = "neither", id1 = "both", "PT", "")
  de.markers(so, gene_set, exposure, id2 = ref_group, id1 = exp_group, cell, "")
  colnames(m)[2] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group,")")
  colnames(m)[3] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",ref_group,")")
  
  # Filter for significant genes
  m_top <- m
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # Select the top 10 positive and 10 negative log2FC genes based on the largest magnitude of fold change
  top_positive_by_fc <- significant_genes %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
    head(10)  # Top 10 positive fold changes
  
  top_negative_by_fc <- significant_genes %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
    head(10)  # Top 10 negative fold changes
  
  # Select the top 10 positive and 10 negative log2FC genes based on significance (p_val_adj)
  top_positive_by_significance <- significant_genes %>% 
    filter(avg_log2FC > 0) %>% 
    arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
    head(10)  # Top 10 most significant positive fold changes
  
  top_negative_by_significance <- significant_genes %>% 
    filter(avg_log2FC < 0) %>% 
    arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
    head(10)  # Top 10 most significant negative fold changes
  
  # Combine top fold-change based and significance-based genes into a final list
  top_genes <- rbind(
    top_positive_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
    top_negative_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
    top_positive_by_significance %>% mutate(Selection = "Top 10 by Significance"),
    top_negative_by_significance %>% mutate(Selection = "Top 10 by Significance")
  )
  
  
  if (!is.null(cell)){
    title <- paste0("Senescence DEGs in ",cell_name," cells for ",condition)
  } else {
    title <- paste0("Senescence Bulk DEGs for ",condition)
  }
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       subtitle = paste0("Positive Log2 FC = Greater Expression in ", condition,"\n",
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
                       max.overlaps=60)
  if (!is.null(cell)){
    filename <- paste0("Senescence_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
  } else {
    filename <- paste0("Senescence_Bulk_DEGs_for_",condition,".pdf") 
  }
  pdf(fs::path(dir.results,filename),width=10,height=7)
  plot(p)
  dev.off()
  
  #GSEA
  if (enrichment=="Yes") {
    #Gene set enrichment analysis
    gc()
    sce_sn_hep <- as.SingleCellExperiment(so)
    ## C2 category is according to canonical pathways: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
    geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
    ### filter background to only include genes that we assessed
    geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
    geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_sn_hep),]
    m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
    stats <- m$p_val_adj
    names(stats) <- rownames(m)
    eaRes <- fgsea(pathways = m_list, stats = na.omit(stats))
    #ooEA <- order(eaRes$pval, decreasing = FALSE)
    #kable(head(eaRes[ooEA, 1:7], n = 20))
    # Convert the leadingEdge column to comma-separated strings
    eaRes$leadingEdge <- sapply(eaRes$leadingEdge, function(x) paste(x, collapse = ", "))
    gc()
    #Significant pathways
    sig <- eaRes %>% 
      filter(padj<0.05)
    
    #Plot top pathways
    # Subset top pathways for visualization
    top_pathways <- eaRes[order(-abs(eaRes$NES)), ][1:top_gsea, ]  # Top 10 pathways based on adjusted p-value
    
    # Define a significance threshold
    significance_threshold <- 0.05
    
    # Add a significance column for coloring
    top_pathways$significance <- ifelse(
      top_pathways$padj > significance_threshold, "Non-significant",
      ifelse(top_pathways$NES > 0, "Positive Significant", "Negative Significant")
    )
    
    if (!is.null(cell)) {
      title <- paste0("Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition)
    } else {
      title <- paste0("Top Enriched Pathways by NES (GSEA) for ",condition," (Bulk)")
    }
    # Create a bar plot
    gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
      geom_bar(stat = "identity") +
      coord_flip() +  # Flip coordinates for better readability
      scale_fill_manual(
        values = c(
          "Positive Significant" = "red",
          "Negative Significant" = "blue",
          "Non-significant" = "gray"
        ),
        name = "Significance"
      ) +
      labs(
        title = title,
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")
      )
    
    
    if (!is.null(cell)){
      filename <- paste0("Senescence_GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
    } else {
      filename <- paste0("Senescence_Bulk_GSEA_for_",condition,".pdf") 
    }
    pdf(fs::path(dir.results,filename),width=20,height=20)
    plot(gsea_plot)
    dev.off()
    
  }
  #Make sure gene nemaes are in printed file
  deg_results <- m 
  deg_results$Gene <- rownames(deg_results)
  
  #save results to excel file
  write_multiple_sheets <- function(output_file, sheet_data) {
    # Create a new workbook
    wb <- createWorkbook()
    
    # Loop over the list of data frames
    for (sheet_name in names(sheet_data)) {
      # Add a new sheet to the workbook
      addWorksheet(wb, sheet_name)
      
      # Write the data frame to the sheet
      writeData(wb, sheet_name, sheet_data[[sheet_name]])
    }
    # Save the workbook to the specified file
    saveWorkbook(wb, file = output_file, overwrite = TRUE)
    
    # Print a message
    message("Excel file with multiple sheets created: ", output_file)
  } 
  # Example usage
  df1 <- deg_results
  if (enrichment=="Yes") {
    df2 <- eaRes
  } else {
    df2 <- NULL
  }
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0("Senescence_Results_",cell_name,"_cells_for_",condition,".xlsx"))  
  } else {
    output_file <- fs::path(dir.results,paste0("Senescence_Bulk_Results_for_",condition,".xlsx"))
  }
  
  sheet_data <- list(
    "DEG" = df1,
    "Pathway_Results" = df2
  )
  
  # Call the function
  write_multiple_sheets(output_file, sheet_data)
}


#Visualize Function ----
##a. Bulk ----
visualize_function_bulk <- function(exposure) {
  
  # List all files in the results directory
  all_files <- dir_ls(path = dir.results, glob = "*.xlsx")  # List only .xlsx files
  
  # Filter the files that contain both the exposure and the word "Bulk" in the file name
  matching_files <- all_files[str_detect(all_files, exposure) & str_detect(all_files, "Bulk")]
  
  # Check if matching files exist for the current exposure
  if (length(matching_files) > 0) {
    # Ensure we have at least 3 files
    if (length(matching_files) >= 3) {
      # Read in the first, second, and third files as group1, group2, and group3
      group1 <- read.xlsx(matching_files[1])
      group2 <- read.xlsx(matching_files[2])
      group3 <- read.xlsx(matching_files[3])
    } else {
      # If there are not enough matching files, return a message
      print("Not enough matching files found.")
      return(NULL)
    }
  } else {
    # If no matching files, return a message
    print("No matching files found.")
    return(NULL)
  }
  
  # Create an empty list to store titles
  group_titles <- list()
  
  # Loop through each file in matching_files to extract titles
  for (i in 1:length(matching_files)) {
    file_path <- matching_files[i]
    
    # Extract the comparison group from the file name
    comparison_group <- sub(".*\\((.*?)\\).*", "\\1", basename(file_path))
    
    # Assign titles to group1_title, group2_title, etc.
    group_titles[[paste0("group", i, "_title")]] <- comparison_group
  }
  
  # Now you can access the group titles like this
  group1_title <- group_titles$group1_title
  group2_title <- group_titles$group2_title
  group3_title <- group_titles$group3_title
  
  # Define significance threshold
  significance_threshold <- 0.05
  
  # Filter significant genes and classify as upregulated or downregulated
  get_gene_sets <- function(data) {
    upregulated <- data %>%
      filter(p_val_adj < significance_threshold & avg_log2FC > 0) %>%
      pull(Gene) # Replace 'Gene' with the actual column name for gene IDs
    downregulated <- data %>%
      filter(p_val_adj < significance_threshold & avg_log2FC < 0) %>%
      pull(Gene) # Replace 'Gene' with the actual column name for gene IDs
    list(up = upregulated, down = downregulated)
  }
  
  # Calculate gene sets for each group
  group1_sets <- get_gene_sets(group1)
  group2_sets <- get_gene_sets(group2)
  group3_sets <- get_gene_sets(group3)
  
  # Extract upregulated and downregulated sets
  up1 <- group1_sets$up
  down1 <- group1_sets$down
  up2 <- group2_sets$up
  down2 <- group2_sets$down
  up3 <- group3_sets$up
  down3 <- group3_sets$down
  
  # Function to draw Venn diagrams for upregulated and downregulated genes
  draw_venn <- function(up1, up2, up3, down1, down2, down3) {
    
    # Plot Venn diagram for upregulated genes
    venn_up <- venn.diagram(
      x = list(
        Group1 = up1,
        Group2 = up2,
        Group3 = up3
      ),
      category.names = c("Group 1", "Group 2", "Group 3"),
      filename = NULL, # Don't save as file, display in plot window
      output = TRUE,
      main = "Upregulated Genes",
      col = "black",
      fill = c("lightblue", "lightgreen", "lightcoral"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.pos = c(0, 0, 0),  # Move Group 3 label to the bottom
      cat.dist = 0.05,
      margin = 0.1
    )
    
    # Plot Venn diagram for downregulated genes
    venn_down <- venn.diagram(
      x = list(
        Group1 = down1,
        Group2 = down2,
        Group3 = down3
      ),
      category.names = c("Group 1", "Group 2", "Group 3"),
      filename = NULL, # Don't save as file, display in plot window
      output = TRUE,
      main = "Downregulated Genes",
      col = "black",
      fill = c("lightblue", "lightgreen", "lightcoral"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.pos = c(0, 0, 0),  # Move Group 3 label to the bottom
      cat.dist = 0.05,
      margin = 0.1
    )
    
    # Draw the upregulated genes Venn diagram
    grid.draw(venn_up)
    
    # Draw a new page for downregulated Venn diagram
    grid.newpage() 
    
    # Draw the downregulated genes Venn diagram
    grid.draw(venn_down)
  }
  
  # Call the function to draw Venn diagrams
  draw_venn(up1, up2, up3, down1, down2, down3)
  
  # Print the titles for the groups
  print(paste0("Group 1 = ", group1_title))
  print(paste0("Group 2 = ", group2_title))
  print(paste0("Group 3 = ", group3_title))
}



##b. Single Cell ----
###i. Venn Diagrams ----
visualize_function <- function(exposure,cell) {
  
  # List all files in the results directory
  all_files <- dir_ls(path = dir.results, glob = "*.xlsx")  # List only .xlsx files
  
  # Filter the files that contain the current cell type and exposure in the file name
  matching_files <- all_files[str_detect(all_files, cell) & str_detect(all_files, exposure)]
  
  # Check if matching files exist for the current cell type and exposure
  if (length(matching_files) > 0) {
    # Ensure we have at least 3 files
    if (length(matching_files) >= 3) {
      # Read in the first, second, and third files as group1, group2, and group3
      group1 <- read.xlsx(matching_files[1])
      group2 <- read.xlsx(matching_files[2])
      group3 <- read.xlsx(matching_files[3])
      
    } 
  }
  # Create an empty list to store titles
  group_titles <- list()
  
  # Loop through each file in matching_files
  for (i in 1:length(matching_files)) {
    file_path <- matching_files[i]
    
    # Extract the comparison group from the file name
    comparison_group <- sub(".*\\((.*?)\\).*", "\\1", basename(file_path))
    
    # Assign titles to group1_title, group2_title, etc.
    group_titles[[paste0("group", i, "_title")]] <- comparison_group
  }
  
  # Now you can access the group titles like this
  group1_title <- group_titles$group1_title
  group2_title <- group_titles$group2_title
  group3_title <- group_titles$group3_title
  
  # Define significance threshold
  significance_threshold <- 0.05
  
  # Filter significant genes and classify as upregulated or downregulated
  get_gene_sets <- function(data) {
    upregulated <- data %>%
      filter(p_val_adj < significance_threshold & avg_log2FC > 0) %>%
      pull(Gene) # Replace 'gene' with the actual column name for gene IDs
    downregulated <- data %>%
      filter(p_val_adj < significance_threshold & avg_log2FC < 0) %>%
      pull(Gene) # Replace 'gene' with the actual column name for gene IDs
    list(up = upregulated, down = downregulated)
  }
  
  # Calculate gene sets for each group
  group1_sets <- get_gene_sets(group1)
  group2_sets <- get_gene_sets(group2)
  group3_sets <- get_gene_sets(group3)
  
  # Extract upregulated and downregulated sets
  up1 <- group1_sets$up
  down1 <- group1_sets$down
  up2 <- group2_sets$up
  down2 <- group2_sets$down
  up3 <- group3_sets$up
  down3 <- group3_sets$down
  
  # Find mutual upregulated genes between groups
  mutual_up12 <- intersect(up1, up2)  # Upregulated in both Group 1 and Group 2
  mutual_up13 <- intersect(up1, up3)  # Upregulated in both Group 1 and Group 3
  mutual_up23 <- intersect(up2, up3)  # Upregulated in both Group 2 and Group 3
  mutual_up_all <- Reduce(intersect, list(up1, up2, up3))  # Upregulated in all three groups
  
  # Find mutual downregulated genes between groups
  mutual_down12 <- intersect(down1, down2)  # Downregulated in both Group 1 and Group 2
  mutual_down13 <- intersect(down1, down3)  # Downregulated in both Group 1 and Group 3
  mutual_down23 <- intersect(down2, down3)  # Downregulated in both Group 2 and Group 3
  mutual_down_all <- Reduce(intersect, list(down1, down2, down3))  # Downregulated in all three groups
  
  # Find unique upregulated genes in each group
  unique_up1 <- setdiff(up1, union(up2, up3))  # Upregulated only in Group 1
  unique_up2 <- setdiff(up2, union(up1, up3))  # Upregulated only in Group 2
  unique_up3 <- setdiff(up3, union(up1, up2))  # Upregulated only in Group 3
  
  # Find unique downregulated genes in each group
  unique_down1 <- setdiff(down1, union(down2, down3))  # Downregulated only in Group 1
  unique_down2 <- setdiff(down2, union(down1, down3))  # Downregulated only in Group 2
  unique_down3 <- setdiff(down3, union(down1, down2))  # Downregulated only in Group 3
  
  
  
  # Function to draw Venn diagrams for upregulated and downregulated genes
  # Function to draw Venn diagrams for upregulated and downregulated genes
  draw_venn <- function(up1, up2, up3, down1, down2, down3) {
    
    # Plot Venn diagram for upregulated genes
    venn_up <- venn.diagram(
      x = list(
        Group1 = up1,
        Group2 = up2,
        Group3 = up3
      ),
      category.names = c("Group 1", "Group 2", "Group 3"),
      filename = NULL, # Don't save as file, display in plot window
      output = TRUE,
      main = "Upregulated Genes",
      col = "black",
      fill = c("lightblue", "lightgreen", "lightcoral"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.pos = c(0, 0, 0),  # Move Group 3 label to the bottom
      cat.dist = 0.05,
      margin = 0.1
    )
    
    # Plot Venn diagram for downregulated genes
    venn_down <- venn.diagram(
      x = list(
        Group1 = down1,
        Group2 = down2,
        Group3 = down3
      ),
      category.names = c("Group 1", "Group 2", "Group 3"),
      filename = NULL, # Don't save as file, display in plot window
      output = TRUE,
      main = "Downregulated Genes",
      col = "black",
      fill = c("lightblue", "lightgreen", "lightcoral"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.pos = c(0, 0, 0),  # Move Group 3 label to the bottom
      cat.dist = 0.05,
      margin = 0.1
    )
    
    # Draw the upregulated genes Venn diagram
    grid.draw(venn_up)
    
    # Draw a new page for downregulated Venn diagram
    grid.newpage() 
    
    # Draw the downregulated genes Venn diagram
    grid.draw(venn_down)
  }
  draw_venn(up1, up2, up3, down1, down2, down3)
  
  print(paste0("Group 1 = ",group1_title),paste0("Group 2 = ",group2_title),paste0("Group 3 = ",group3_title))
}

###ii. Bar/Dot plot ----
#exposure = "Group2"
#exp_group = "glp1_exclusive"
#ref_group = "neither"
#cell_types <- c("Hepatocyte","EC","Stellate","Cholang","NKC/NKT","Kup/MON","dHep","Kup/MAC","B/Plasma","Hep/Immune")
#cell_types <- str_replace_all(cell_type,"/","_")

# Function to extract cell type from file name using a list of cell types
extract_cell_type <- function(filename, cell_types) {
  for (cell_type in cell_types) {
    if (str_detect(filename, cell_type)) {
      return(cell_type)
    }
  }
  return(NA)  # Return NA if no cell type is found
}

visualize_cell <- function(exposure,exp_group,ref_group) {
  
  # List all files in the results directory
  all_files <- dir_ls(path = dir.results, glob = "*.xlsx")  # List only .xlsx files
  
  # Filter the files that contain both the exposure and the word "Bulk" in the file name
  matching_files <- all_files[str_detect(all_files, exposure) & str_detect(all_files, exp_group) & str_detect(all_files, ref_group) & !str_detect(all_files,"Bulk") & str_detect(all_files,"Senescence")]
  
  # Check if matching files exist for the current cell type and exposure
  if (length(matching_files) > 0) {
    # Loop through each matching file
    for (file in matching_files) {
      # Extract cell type from the filename
      cell_type <- extract_cell_type(basename(file), cell_types)
      
      # Dynamically assign the data to a variable named after the cell type
      if (!is.na(cell_type)) {
        assign(make.names(cell_type), read.xlsx(file), envir = .GlobalEnv)
      }
    }
  }
  
  # Initialize a list to store data
  all_data <- list()
  # Loop through each file and process
  for (file in matching_files) {
    # Read the file (assuming CSV format, modify if needed)
    data <- read.xlsx(file)
    
    # Extract cell type from the filename
    cell_type <- extract_cell_type(basename(file), cell_types)
    
    # Add the cell type as a column
    data$Cell_Type <- cell_type
    
    # Append to the list
    all_data[[length(all_data) + 1]] <- data
  }
  
  # Combine all data into a single dataframe
  combined_data <- bind_rows(all_data)
  combined_data$sig <- ifelse(combined_data$p_val_adj<0.05,"Significant","Non-Significant")
  combined_data$sig <- factor(combined_data$sig)
  
  # Filter significant data
  significant_data <- combined_data %>%
    filter(sig == "Significant")
  
  # Get top 3 positive and negative genes for each Cell_Type
  top_positive_per_cell_type <- significant_data %>%
    filter(avg_log2FC > 0) %>%  # Only positive values
    group_by(Cell_Type) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 3) %>%
    mutate(pos_or_neg = "Positive")
  
  top_negative_per_cell_type <- significant_data %>%
    filter(avg_log2FC < 0) %>%  # Only negative values
    group_by(Cell_Type) %>%
    arrange(avg_log2FC) %>%  # Sort in ascending order for negative values
    slice_head(n = 3) %>%
    mutate(pos_or_neg = "Negative")
  
  # Combine both positive and negative data
  top_genes_per_cell_type <- bind_rows(top_positive_per_cell_type, top_negative_per_cell_type)
  
  range_data <- combined_data %>%
    group_by(Cell_Type) %>%
    summarise(
      min_log2FC = min(avg_log2FC, na.rm = TRUE),
      max_log2FC = max(avg_log2FC, na.rm = TRUE)
    )
  exposure <- "Treatment"
  # Create the plot
  p <- ggplot(combined_data, aes(x = Cell_Type, y = avg_log2FC)) +
    # Add the bar chart (range from min to max)
    #geom_linerange(
    #  data = range_data,
    #  aes(x = Cell_Type, ymin = min_log2FC, ymax = max_log2FC),
    #  color = "lightblue",
    #  linewidth = 20,
    #  alpha = 0.5,
    #  inherit.aes = FALSE  # Disable inheriting global aesthetics
    #) +
    geom_rect(
      data = range_data,
      aes(
        xmin = as.numeric(as.factor(Cell_Type)) - 0.4,
        xmax = as.numeric(as.factor(Cell_Type)) + 0.4,
        ymin = min_log2FC,
        ymax = max_log2FC
      ),
      fill = "lightblue",
      color = "black", # Black border
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5)+
    # Overlay jittered dots
    geom_jitter(aes(color = sig), width = 0.3, alpha = 0.5, size = 0.3) +
    
    # Custom colors for dots
    scale_color_manual(values = c("Significant" = "red", "Non-Significant" = "darkblue")) +
    
    # Add non-overlapping text labels
    geom_text_repel(
      data = top_genes_per_cell_type,
      aes(label = Gene),
      color = "black",
      size = 2,
      box.padding = 0.35,      # Padding around text boxes
      point.padding = 0.2,    # Padding around data points
      segment.color = "gray", # Line color connecting labels to points
      max.overlaps = 10,       # Allow up to 10 overlaps (adjust as needed)
      angle = 45 
    ) +
    
    # Adjust theme and labels
    theme_minimal() +
    theme_minimal() +
    labs(
      title = paste0("Top Senescence Genes by ",exposure," Group (",exp_group," vs. ",ref_group,") by Cell Type"),
      x = "Cell Type",
      y = "Average Log2 Fold Change"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  pdf(fs::path(dir.results,paste0(exposure,"_",exp_group,"_",ref_group,".pdf")),width=10,height=10)
  plot(p)
  dev.off()
}

#so=so_sens_all
#cell=NULL
#exposure=exp
#covariate=NULL
#gene_set=rownames(so_sens_all)
#batch_size=length(gene_set)
#exp_group="Yes"
#ref_group="No"

#Mast Function----
##a. All Genes ----
mast_fxn <- function(so,cell,exposure,covariate,gene_set,batch_size,exp_group,ref_group) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  if (!is.null(exp_group)) {
    condition <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group," vs. ",ref_group,")")
  } else {
    condition <- paste0(str_to_title(exposure))
  }
  
  #Filter so to specified cell type
  if (!is.null(cell)){
    so <- subset(so,celltype2==cell)
    DefaultAssay(so) <- "RNA"
  } 
  
  #Association by Group
  combined_results <- data.frame()
  num_genes <- length(gene_set)
  for (start_batch in seq(1,num_genes,by=batch_size)) {
    
    end_batch <- min(start_batch + batch_size - 1, num_genes)
    batch <- gene_set[start_batch:end_batch]
    
    # Subset the Seurat object to include only genes in the gene set
    subset <- intersect(batch, rownames(so))
    so_sub <- subset(so_sub, features = subset)
    
    # Extract the expression data matrix from so (e.g., normalized counts)
    expression_matrix <- as.matrix(GetAssayData(so_sub, layer = "data"))
    
    # Extract metadata
    cell_metadata <- so_sub@meta.data
    #cell_metadata[[exposure]] <- factor(cell_metadata[[exposure]])
    
    # Create SingleCellAssay object
    sca <- FromMatrix(exprsArray = expression_matrix, cData = cell_metadata)
    
    # Assuming gene_set is a vector of gene names
    sca_gene_set <- sca[rownames(sca) %in% gene_set, ]
    
    #Define the formula
    model_formula <- as.formula(paste0("~", exposure))
    exp <- paste0(exposure,exp_group)
    
    # Run the linear model with zlm
    zlm_results <- zlm(formula = model_formula, sca = sca_gene_set)
    
    # Summarize results and perform likelihood ratio test (LRT) for outcome
    summary_zlm <- summary(zlm_results, doLRT = exp)
    
    # Convert the summary to a data table for easy manipulation
    summary_dt <- summary_zlm$datatable
    
    #Overall with Hurdle and logFC
    fcHurdle1 <- merge(summary_dt[contrast==exp & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                       summary_dt[contrast==exp & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle1[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdle1$component <- "FC_H"
    fcHurdle2 <- merge(summary_dt[contrast==exp & component=='C',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                       summary_dt[contrast==exp & component=='C', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle2[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdle2$component <- "C"
    fcHurdle3 <- merge(summary_dt[contrast==exp & component=='D',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                       summary_dt[contrast==exp & component=='D', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle3[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdle3$component <- "D"
    
    fcHurdle <- rbind(fcHurdle1,fcHurdle2)
    fcHurdle <- rbind(fcHurdle,fcHurdle3)
    
    result1 <- fcHurdle
    result1 <- result1 %>%
      filter(fdr<0.05) %>% 
      dplyr::rename(LogFoldChange=coef) %>% 
      dplyr::rename(Gene=primerid)
    
    # Append the batch results to the combined results using rbindlist
    # combined_results <- rbind(combined_results, batch_results)
    combined_results <- rbindlist(list(combined_results, result1), use.names = TRUE, fill = TRUE)
    
  }
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- paste0("Results_",cell_name,"_cells_for_",condition,".xlsx")
  } else {
    output_file <- paste0("Results_for_",condition,".xlsx")
  }
  
  write.csv(combined_results,fs::path(output_file))
  
  # Filter top 10 positive and negative log2FoldChange
  top_pos <- as.data.frame(combined_results) %>%
    filter(component=="FC_H") %>% 
    filter(fdr<0.05) %>%
    filter(LogFoldChange>0) %>% 
    arrange(desc(LogFoldChange)) %>%
    dplyr::slice(1:30) 
  # dplyr::rename(Gene=primerid,
  #               LogFC=coef)
  
  top_neg <- as.data.frame(combined_results) %>%
    filter(component=="FC_H") %>% 
    filter(fdr<0.05) %>%
    filter(LogFoldChange<0)  %>% 
    arrange(LogFoldChange) %>%
    dplyr::slice(1:30) 
  # dplyr::rename(Gene=primerid,
  #               LogFC=coef)
  if(dim(top_pos)[[1]] > 0 | dim(top_neg)[[1]] > 0) {
    top_pos$Direction <- "Positive"
    top_neg$Direction <- "Negative"
    
    # Combine and prepare for plotting
    top_genes <- bind_rows(top_pos, top_neg)
    
    # Set sorting order: all negative genes first, then all positive genes
    top_genes <- top_genes %>%
      mutate(Gene = factor(Gene, levels = Gene[order(Direction, LogFoldChange)]))  # Order by Direction and log2FC
    
    # # Bar chart with flipped x and y axes
    p <- ggplot(top_genes, aes(x = LogFoldChange, y = Gene, fill = Direction)) +
      geom_bar(stat = "identity") +
      labs(x = "logFC", y = "Gene",title=paste0(condition,"vs. ",reference," among", cell," cells in youth with Type 2 Diabetes")) +
      scale_fill_manual(values = c("blue", "red"), labels = c("Lower Expression", "Higher Expression")) +
      theme_minimal()+
      theme(element_text(family="Times"))+
      theme(legend.position="none")
    
    if (!is.null(cell)){
      figure_file<- paste0("MAST_Barchart_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      figure_file <- paste0("MAST_Barchart_for_",condition,".pdf")
    }
    
    pdf(fs::path(dir.results,figure_file),width=20,height=20)
    plot(p)
    dev.off()
  }
}

##b. Senescence Genes ----
mast_sens_fxn <- function(so,cell,exposure,covariate,gene_set,exp_group,ref_group) {
  DefaultAssay(so) <- "RNA"
  
  #Condition names
  if (!is.null(exp_group)) {
    condition <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group," vs. ",ref_group,")")
  } else {
    condition <- paste0(str_to_title(exposure))
  }
  
  #Filter so to specified cell type
  if (!is.null(cell)){
    so <- subset(so,celltype2==cell)
    DefaultAssay(so) <- "RNA"
  } 
  
  # Extract the expression data matrix from so (e.g., normalized counts)
  expression_matrix <- as.matrix(GetAssayData(so, layer = "data"))
  
  # Extract metadata
  cell_metadata <- so@meta.data
  #cell_metadata[[exposure]] <- factor(cell_metadata[[exposure]])
  
  # Create SingleCellAssay object
  sca <- FromMatrix(exprsArray = expression_matrix, cData = cell_metadata)
  
  # Assuming gene_set is a vector of gene names
  sca_gene_set <- sca[rownames(sca) %in% gene_set, ]
  
  #Check for missing data in exposure
  if (sum(is.na(colData(sca_gene_set)[exposure]))>0) {
    sca_filtered <- sca_gene_set[, !is.na(colData(sca_gene_set)[exposure])]
  } else {
    sca_filtered <- sca_gene_set
  }
  #sum(is.na(colData(sca_filtered)[exposure]))
  
  #Define the formula
  model_formula <- as.formula(paste0("~", exposure))
  exp <- paste0(exposure,exp_group)
  
  # Run the linear model with zlm
  zlm_results <- zlm(formula = model_formula, sca = sca_filtered)
  
  # Summarize results and perform likelihood ratio test (LRT) for outcome
  summary_zlm <- summary(zlm_results, doLRT = exp)
  
  # Convert the summary to a data table for easy manipulation
  summary_dt <- summary_zlm$datatable
  
  #Overall with Hurdle and logFC
  fcHurdle1 <- merge(summary_dt[contrast==exp & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                     summary_dt[contrast==exp & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle1[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle1$component <- "FC_H"
  fcHurdle2 <- merge(summary_dt[contrast==exp & component=='C',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                     summary_dt[contrast==exp & component=='C', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle2[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle2$component <- "C"
  fcHurdle3 <- merge(summary_dt[contrast==exp & component=='D',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                     summary_dt[contrast==exp & component=='D', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle3[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle3$component <- "D"
  
  #fcHurdle <- rbind(fcHurdle1,fcHurdle2)
  #fcHurdle <- rbind(fcHurdle,fcHurdle3)
  
  result1 <- fcHurdle1
  result1 <- result1 %>%
    filter(fdr<0.05) %>% 
    dplyr::rename(LogFoldChange=coef) %>% 
    dplyr::rename(Gene=primerid)
  
  # Specify the file name and data
  if (!is.null(cell)){
    cell_name <- str_replace_all(cell,"/","_")
    output_file <- paste0("Results_",cell_name,"_cells_for_",condition,".csv")
  } else {
    output_file <- paste0("Results_for_",condition,".csv")
  }
  
  write.csv(result1,fs::path(dir.results,output_file))
  
  # Filter top 10 positive and negative log2FoldChange
  top_pos <- as.data.frame(result1) %>%
    filter(component=="FC_H") %>% 
    filter(fdr<0.05) %>%
    filter(LogFoldChange>0) %>% 
    arrange(desc(LogFoldChange)) %>%
    dplyr::slice(1:30) 
  # dplyr::rename(Gene=primerid,
  #               LogFC=coef)
  
  top_neg <- as.data.frame(result1) %>%
    filter(component=="FC_H") %>% 
    filter(fdr<0.05) %>%
    filter(LogFoldChange<0)  %>% 
    arrange(LogFoldChange) %>%
    dplyr::slice(1:30) 
  # dplyr::rename(Gene=primerid,
  #               LogFC=coef)
  if(dim(top_pos)[[1]] > 0 | dim(top_neg)[[1]] > 0) {
    top_pos$Direction <- "Positive"
    top_neg$Direction <- "Negative"
    
    # Combine and prepare for plotting
    top_genes <- bind_rows(top_pos, top_neg)
    
    # Set sorting order: all negative genes first, then all positive genes
    top_genes <- top_genes %>%
      mutate(Gene = factor(Gene, levels = Gene[order(Direction, LogFoldChange)]))  # Order by Direction and log2FC
    
    if (!is.null(cell)) {
      title=paste0(condition,"in ",cell," cells")
    } else {
      title=paste0(condition)
    }
    
    # # Bar chart with flipped x and y axes
    p <- ggplot(top_genes, aes(x = LogFoldChange, y = Gene, fill = Direction)) +
      geom_bar(stat = "identity") +
      labs(x = "logFC", y = "Gene",title=title) +
      scale_fill_manual(values = c("blue", "red"), labels = c("Lower Expression", "Higher Expression")) +
      theme_minimal()+
      theme(element_text(family="Times"))+
      theme(legend.position="none")
    
    if (!is.null(cell)){
      figure_file<- paste0("MAST_Barchart_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      figure_file <- paste0("MAST_Barchart_for_",condition,".pdf")
    }
    
    pdf(fs::path(dir.results,figure_file),width=7,height=7)
    plot(p)
    dev.off()
    
  }
}



#Ideas Function ----
ideas_fxn <- function(so, exp_group, ref_group, covariates, var_type, id, method, exposure) {
  
  #Filter so to groups of interest
  so$filt_condition <- ifelse(so$group2==paste0(exp_group) | so$group2==paste0(ref_group),"Yes","No")
  so <- subset(so,filt_condition=="Yes")
  so$group3 <- ifelse(so$group2==paste0(exp_group),paste0(exp_group),paste0(ref_group))
  so$group3 <- factor(so$group3)
  
  # Create individual variable dynamically based on the input 'id'
  so$individual <- so[[id]]  # Assign the column specified by 'id' as the individual variable
  
  # Split genes into batches
  batch_size <- 1000  # Number of genes per batch
  count_matrix <- GetAssayData(so, assay = "RNA", layer = "counts")
  gene_batches <- split(rownames(count_matrix), ceiling(seq_along(rownames(count_matrix)) / batch_size))
  
  # Loop through each gene batch
  all_results <- list()  # To store results from all batches
  
  for (i in seq_along(gene_batches)) {
    cat("Processing batch", i, "of", length(gene_batches), "\n")
    
    # Get the list of genes for the current batch
    batch_genes <- gene_batches[[i]]
    
    # Extract the count matrix for this batch of genes
    batch_matrix <- count_matrix[batch_genes, , drop = FALSE]  # Drop = FALSE ensures it stays a matrix
    batch_matrix <- as.matrix(batch_matrix)
    # Round the count matrix to ensure integer values
    batch_matrix <- round(batch_matrix)  # Round values to ensure integer counts
    
    # Check if the matrix contains non-integer values
    if (any(batch_matrix != floor(batch_matrix))) {
      stop("Batch matrix contains non-integer values!")
    }
    
    # Extract the metadata for the cells corresponding to the batch
    meta_cell <- so@meta.data
    meta_cell$cell_id <- rownames(meta_cell)  # Add cell IDs if not already present
    
    # Ensure the metadata includes the required column (e.g., 'Kit_Lot' for individuals)
    meta_cell$cell_rd <- colSums(batch_matrix)  # Total counts per cell (read depth)
    
    # Add a small constant to ensure no zero or negative values
    meta_cell$cell_rd <- meta_cell$cell_rd + 1  # Adding 1 to avoid log(0) errors
    
    #    # Ensure the metadata includes the required column (e.g., 'Kit_Lot' for individuals)
    # meta_cell$cell_rd[meta_cell$cell_rd == 0] <- 1  # Ensure no zero read depth
    
    # Create individual-level metadata (e.g., group information, exposure)
    meta_ind <- unique(meta_cell[, c("individual", "group3")])  # Use dynamic 'id' for individual variable
    
    # Set variables for differential expression analysis
    var2test <- "group3"  # Group variable to test
    var2adjust <- NULL     # Adjusting variable (e.g., RIN, if applicable)
    var2test_type <- var_type  # "binary" or "continuous"
    var_per_cell <- "cell_rd"  # Cell-level adjustment variable (read depth)
    test_method <- paste0(method)
    
    # Perform IDEAS distance calculation for this batch
    dist_batch <- ideas_dist(
      count_input = batch_matrix,
      meta_cell = meta_cell,
      meta_ind = meta_ind,
      var_per_cell = var_per_cell,
      var2test = var2test,
      var2test_type = var2test_type,
      d_metric = "Was",  # Wasserstein distance
      fit_method = test_method  # Negative binomial fit
    )
    
    # Perform PERMANOVA for this batch
    pval_batch <- permanova(
      dist_array = dist_batch,
      meta_ind = meta_ind,
      var2test = var2test,
      var2adjust = NULL,
      var2test_type = var2test_type,
      n_perm = 999,
      r.seed = 903
    )
    
    # Combine results for this batch
    batch_results <- data.frame(
      gene = batch_genes,
      pvalue = pval_batch
    )
    
    # Store results
    all_results[[i]] <- batch_results
  }
  
  # Combine all batches into one final result
  final_results <- do.call(rbind, all_results)
  
  # Save or further analyze the results
  write.csv(final_results, "DEA_results_ideas.csv", row.names = FALSE)
  
  # Example mock data for log-fold changes (LFC); Replace with actual data
  set.seed(42)
  final_results$logFC <- rnorm(nrow(final_results), mean = 0, sd = 2)  # Mock log fold change values
  
  # Set thresholds for significance
  pval_threshold <- 0.05
  logFC_threshold <- 1
  
  # Classify genes as upregulated, downregulated, or not significant
  final_results <- final_results %>%
    mutate(
      significance = case_when(
        pvalue < pval_threshold & logFC > logFC_threshold ~ "Upregulated",
        pvalue < pval_threshold & logFC < -logFC_threshold ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Create the volcano plot
  volcano_plot <- ggplot(final_results, aes(x = logFC, y = -log10(pvalue), color = significance)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
    labs(
      title = "Volcano Plot",
      x = "Log2 Fold Change (logFC)",
      y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Save the plot
  ggsave("volcano_plot_ideas.png", volcano_plot, width = 8, height = 6)
  
  # Optionally, print the plot to the console
  print(volcano_plot)
}


# GeomSplitViolin <- ggproto(
#   "GeomSplitViolin", 
#   GeomViolin, 
#   draw_group = function(self, data, ..., draw_quantiles = NULL) {
#     data <- transform(data, 
#                       xminv = x - violinwidth * (x - xmin), 
#                       xmaxv = x + violinwidth * (xmax - x))
#     grp <- data[1,'group']
#     newdata <- plyr::arrange(
#       transform(data, x = if(grp%%2==1) xminv else xmaxv), 
#       if(grp%%2==1) y else -y
#     )
#     newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
#     newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
#     if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
#       stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
#       quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
#       aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
#       aesthetics$alpha <- rep(1, nrow(quantiles))
#       both <- cbind(quantiles, aesthetics)
#       quantile_grob <- GeomPath$draw_panel(both, ...)
#       ggplot2:::ggname("geom_split_violin", 
#                        grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
#     } else {
#       ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
#     }
#   }
# )
# 
# geom_split_violin <- function (mapping = NULL, 
#                                data = NULL, 
#                                stat = "ydensity", 
#                                position = "identity", ..., 
#                                draw_quantiles = NULL, 
#                                trim = TRUE, 
#                                scale = "area", 
#                                na.rm = FALSE, 
#                                show.legend = NA, 
#                                inherit.aes = TRUE) {
#   layer(data = data, 
#         mapping = mapping, 
#         stat = stat, 
#         geom = GeomSplitViolin, 
#         position = position, 
#         show.legend = show.legend, 
#         inherit.aes = inherit.aes, 
#         params = list(trim = trim, 
#                       scale = scale, 
#                       draw_quantiles = draw_quantiles, 
#                       na.rm = na.rm, ...)
#   )
# }
# 
# split.vp <- function(seurat_object, genes, filepath, color1 = "#6c9a8b", color2 = "#e8998d") {
#   for (i in 1:length(genes)){
#     cat("\n")
#     cat("###", genes[i])
#     cat("\n")  
#     d = VlnPlot(seurat_object, features = genes[i], split.by = "Group", idents = "PT", split.plot = F, pt.size = 0) 
#     d = d$data
#     p = ggplot(d,aes(x=ident, y = !!sym(genes[i]), fill=split))+
#       geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), 
#                   size = 0.2, alpha = 0.3, show.legend = T, aes(color = split)) +
#       geom_split_violin(trim = T) +
#       theme_bw()+
#       theme(legend.title = element_blank(),
#             axis.title.x = element_blank(), 
#             axis.text.x=element_blank(),
#             plot.title = element_text()) +
#       labs(title = (paste0(genes[i]," in  PT Cells")),
#            y = "Expression") +
#       scale_fill_manual(values=c(color1, color2)) + 
#       scale_color_manual(values=c(color1, color2))
#     print(p)
#     cat("\n")
#     # Save
#     ggsave(filename = paste0(filepath,"Violin_",genes[i],".jpeg"),plot = p,scale = 5,
#            width = 800,height = 600,units = "px")
#   }
# }
# 
# split.vp.combined <- function(seurat_object, genes, filepath, color1 = "#6c9a8b", color2 = "#e8998d", idents = NULL) {
#   compiled_d = data.frame()
#   for (i in 1:length(genes)){
#     d = VlnPlot(seurat_object, features = genes[i], split.by = "Group", idents = idents, split.plot = F, pt.size = 0) 
#     d = d$data
#     d$genename = colnames(d)[1]
#     colnames(d)[1] <- "expression"
#     compiled_d = rbind(compiled_d, d)
#   }
#   p = 
#     ggplot(compiled_d,aes(x=genename, y = expression, fill=split)) +
#     geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), 
#                 size = 0.2, alpha = 0.3, show.legend = T, aes(color = split)) +    
#     geom_split_violin(scale = "width", trim = T) +
#     theme_bw() +
#     theme(legend.title = element_blank(),
#           axis.title.x = element_blank(), 
#           plot.title = element_text()) +
#     labs(y = "Expression") +
#     scale_fill_manual(values=c(color1, color2)) + 
#     scale_color_manual(values=c(color1, color2))
#   print(p)
#   cat("\n")
#   # Save
#   ggsave(filename = paste0(filepath,"Violin_combined",".jpeg"),plot = p,scale = 5,
#          width = 1000,height = 600,units = "px")
#   
# }
# 
# # Formatted dot plot function
# # colorlow = "#8ecae6", colormid = "#fcbf49", colorhigh = "#d90429"
# dp.formatted <- function(seurat_object, genes, celltype, group.by, m,
#                          colorlow = "#83c5be", colormid = "#f4f1bb", colorhigh = "#d90429")
# {
#   pt.combined <- DotPlot(seurat_object,
#                          features = genes,idents = celltype, group.by = group.by,
#                          scale = F, cols = "RdYlBu"
#   )$data 
#   
#   pt.plot <- pt.combined %>% 
#     ggplot(aes(x=features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) + 
#     geom_point() +
#     theme_bw() +
#     scale_color_gradient2(low = colorlow, mid = colormid, high = colorhigh, midpoint = 2,
#                           guide = guide_colorbar(label.vjust = 0.8, ticks = F, draw.ulim = T, draw.llim = T),
#                           limits = c(0,4)) +
#     scale_size(range = c(0,4), 
#                limits = c(1,80)) +
#     theme(panel.grid = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#           legend.text = element_text(size = 8),
#           legend.title = element_text(size = 8, vjust = 0.5),
#           legend.spacing.x = unit(.1, "cm"),
#           legend.direction = "horizontal") +
#     guides(size = guide_legend(label.position = "bottom",
#                                title.position = "top"),
#            color = guide_colorbar(label.position = "bottom",
#                                   title.position = "top")) +
#     labs(color = "Scaled average expression",
#          size = "Expression (%) ") + 
#     scale_y_discrete(limits=rev)
#   
#   pt.table <- m %>%
#     filter(rownames(m) %in% genes) %>%
#     dplyr::mutate(p_val_rounded = round(p_val, 4),
#                   p_val = p_format(p_val_rounded, trailing.zero = T, accuracy = 0.001, digits = 3),
#                   p_val_adj_rounded = round(p_val_adj, 4),
#                   p_val_adj = p_format(p_val_adj_rounded, trailing.zero = T, accuracy = 0.001, digits = 3),
#                   pct.1 = sprintf("%.3f", pct.1),
#                   pct.2 = sprintf("%.3f", pct.2),
#                   avg_log2FC = sprintf("%.3f", avg_log2FC)) %>% 
#     dplyr::select(pct.1, pct.2, avg_log2FC, p_val, p_val_adj) 
#   gg.pt.table <- ggtexttable(pt.table,
#                              cols = c("T1D", "HC", "Log2FC", "p-value", "q-value"),
#                              theme = ttheme("blank")) %>%
#     tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1) %>%
#     tab_add_title("% Expressed in PT Cells")
#   
#   pt.plot_table <- ggarrange(pt.plot, NULL, gg.pt.table,
#                              nrow = 1, widths = c(1,-0.1,1), common.legend = F,
#                              legend = "top")
#   
# }
# 
# metabo_sc <- function(data, gene, transcripts, gene_name) {
#   plots <- list()
#   for (i in 1:length(transcripts)) {
#     plot <- subset(data, apply(!is.na(data[, c(transcripts[i], gene)]), 1, all)) %>%
#       ggplot(aes_string(y = gene, x = transcripts[i], color = "group")) +
#       geom_smooth(method = lm, se = F, aes(color=NULL), color = "black",
#                   linetype = "dashed") +
#       geom_point() +
#       geom_smooth(method = lm, se = F) +
#       labs(color = "Group",
#            y = gene_name,
#            x = if(is.null(label(data[[transcripts[i]]]))) (transcripts[i]) else label(data[[transcripts[i]]])) +
#       theme_bw()
#     plots[[i]] <- plot
#   }
#   ggarrange(plotlist = plots, common.legend = T)
# }
# 
# add_direction <- function(df) {
#   df <- df %>%
#     mutate(direction = case_when(
#       (avg_log2FC > 0 & p_val < 0.05) ~ "Upregulated", 
#       (avg_log2FC <= 0 & p_val < 0.05) ~ "Downregulated",
#       TRUE ~ "NS"
#     ))
#   df$direction <- factor(df$direction, levels = c("Downregulated", "Upregulated", "NS") )
#   return(df)
# }
# 

# # #GSEA Functions ----
# ##a. Set up functions ----
# # function for GSEA
# # Set relevant paths
# list.files()
# bg_path <- c()
# # bg_path <- c(fs::path(dir.dat,"/GSEA/"))
# list.files(bg_path)
# 
# # Functions
# ## Function: Adjacency matrix to list
# matrix_to_list <- function(pws){
#   pws.l <- list()
#   for (pw in colnames(pws)) {
#     pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
#   }
#   return(pws.l)
# }
# 
# ## Function: prepare_gmt
# prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
# 
#   # Read in gmt file
#   gmt <- gmtPathways(gmt_file)
#   hidden <- unique(unlist(gmt))
# 
#   # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
#   mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
#                 nrow = length(hidden), ncol = length(gmt))
#   for (i in 1:dim(mat)[2]){
#     mat[,i] <- as.numeric(hidden %in% gmt[[i]])
#   }
# 
#   #Subset to the genes that are present in our data to avoid bias
#   hidden1 <- intersect(genes_in_data, hidden)
#   mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
#   # And get the list again
#   final_list <- matrix_to_list(mat) # for this we use the function we previously defined
# 
#   if(savefile){
#     saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
#   }
# 
#   print('Wohoo! .gmt conversion successfull!:)')
#   return(final_list)
# }
# 
# # List of acronyms to uppercase
# acronyms <- c("rna", "dna", "mtor", "foxo", "ppar", "nmd", "fgfr", "robo",
#               "bhl", "cov", "jak", "stat", "wnt", "hiv", "bcl", "mapk",
#               "pt", "tal", "pc", "ic", "ec", "fibvsmcp")
# # unique(full_results$Gene))
# special_mixed <- c("rrna", "mrna", "trna", "gtpase", "atpase", "robos", "slits", "fibvsmcp")
# special_replacements <- c("rRNA", "mRNA", "tRNA", "GTPase", "ATPase", "ROBOs", "SLITs", "FIB/VSMC/P")
# 
# replace_mixed_case <- function(text, from, to) {
#   for (i in seq_along(from)) {
#     pattern <- paste0("\\b", from[i], "\\b")
#     text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), to[i])
#   }
#   return(text)
# }
# 
# capitalize_acronyms <- function(text, terms) {
#   for (term in terms) {
#     pattern <- paste0("\\b", term, "\\b")
#     replacement <- toupper(term)
#     text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), replacement)
#   }
#   return(text)
# }
# 
# ## b. Fgsea function ----
# plot_fgsea <- function(fgsea_res,
#                        top_n = 30,
#                        title = "Top Enriched Pathways",
#                        xlimit = 3,
#                        xnudge = xlimit/100,
#                        text1 = 6.5,
#                        text2 = 18,
#                        text3 = 20,
#                        face = "plain") {
# 
#   fgsea_res <- fgsea_res %>%
#     arrange(pval) %>%
#     head(top_n) %>%
#     mutate(
#       direction = case_when(NES < 0 ~ "Negative", NES > 0 ~ "Positive"),
#       pathway_clean = str_remove(pathway, "^KEGG_"),
#       pathway_clean = str_remove(pathway_clean, "^REACTOME_"),
#       pathway_clean = str_remove(pathway_clean, "^GOBP_"),
#       pathway_clean = str_remove(pathway_clean, "^GOMF_"),
#       pathway_clean = str_replace_all(pathway_clean, "_", " "),
#       pathway_clean = str_to_sentence(pathway_clean),
#       pathway_clean = str_replace_all(pathway_clean, "\\bi\\b", "I"),
#       pathway_clean = str_replace_all(pathway_clean, "\\bii\\b", "II"),
#       pathway_clean = str_replace_all(pathway_clean, "\\biii\\b", "III"),
#       pathway_clean = str_replace_all(pathway_clean, "\\biv\\b", "IV"),
#       pathway_clean = str_replace_all(pathway_clean, "\\bv\\b", "V"),
#       pathway_clean = str_replace_all(pathway_clean, regex("\\(immune\\)", ignore_case = TRUE), "(IMMUNE)"),
#       pathway_clean = capitalize_acronyms(pathway_clean, acronyms),
#       pathway_clean = replace_mixed_case(pathway_clean, special_mixed, special_replacements),
#       pathway_clean = paste0(pathway_clean, " (", size, ")")
#     ) %>%
#     arrange(pval)
# 
#   fgsea_res$pathway_clean <- reorder(fgsea_res$pathway_clean, fgsea_res$pval)
# 
#   fgsea_res %>%
#     ggplot(aes(x = -log10(pval), y = fct_rev(pathway_clean), label = pathway_clean)) +
#     geom_point(aes(size = abs(NES), color = direction, alpha = 0.8)) +
#     geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
#     geom_text(aes(group = pathway_clean, color = direction, fontface = face),
#               hjust = 0, size = text1, nudge_x = xnudge) +
#     scale_size_binned() +
#     scale_color_manual(values = c("Positive" = "#c75146", "Negative" = "#2c7da0")) +
#     scale_x_continuous(limits = c(0, xlimit), expand = expansion(mult = c(0, 0))) +
#     scale_y_discrete(expand = expansion(add = 1)) +
#     labs(
#       x = "-log(p-value)",
#       y = "Pathways",
#       color = "Direction",
#       size = "NES",
#       title = title
#     ) +
#     guides(alpha = "none") +
#     theme_bw() +
#     theme(
#       axis.text.y = element_blank(),
#       panel.grid = element_blank(),
#       axis.text.x = element_text(size = text3),
#       axis.title = element_text(size = text3),
#       axis.ticks.y = element_blank(),
#       legend.position = c(0.9, 0.2),
#       legend.background = element_blank(),
#       legend.box.background = element_rect(color = "black"),
#       legend.title = element_text(size = text2),
#       legend.text = element_text(size = text2),
#       title = element_text(size = text3)
#     )
# }
# ##c. Transposed fgsea function ----
# plot_fgsea_transpose <- function(fgsea_res,
#                                  top_n = 30,
#                                  title = "Top Enriched Pathways",
#                                  xlimit = 3,
#                                  xnudge = xlimit/100,
#                                  text1 = 6.5,
#                                  text2 = 18,
#                                  text3 = 20) {
# 
#   fgsea_res <- fgsea_res %>%
#     arrange(pval) %>%
#     head(top_n) %>%
#     mutate(
#       direction = case_when((NES < 0 & pval <= 0.05 ~ "Negative"),
#                             (NES > 0 & pval <= 0.05 ~ "Positive"),
#                             (NES < 0 & pval > 0.05 ~ "Negative p > 0.05"),
#                             (NES > 0 & pval > 0.05 ~ "Positive p > 0.05")),
#       face = case_when((NES < 0 & pval <= 0.05 ~ "plain"),
#                        (NES > 0 & pval <= 0.05 ~ "plain"),
#                        (NES < 0 & pval > 0.05 ~ "plain"),
#                        (NES > 0 & pval > 0.05 ~ "plain")),
#       pathway_clean = str_remove(pathway, "^KEGG_"),
#       pathway_clean = str_remove(pathway_clean, "^REACTOME_"),
#       pathway_clean = str_remove(pathway_clean, "^GOBP_"),
#       pathway_clean = str_remove(pathway_clean, "^GOMF_"),
#       pathway_clean = str_replace_all(pathway_clean, "_", " "),
#       pathway_clean = str_to_sentence(pathway_clean),
#       pathway_clean = str_replace_all(pathway_clean, "\\bi\\b", "I"),
#       pathway_clean = str_replace_all(pathway_clean, "\\bii\\b", "II"),
#       pathway_clean = str_replace_all(pathway_clean, "\\biii\\b", "III"),
#       pathway_clean = str_replace_all(pathway_clean, "\\biv\\b", "IV"),
#       pathway_clean = str_replace_all(pathway_clean, "\\bv\\b", "V"),
#       pathway_clean = str_replace_all(pathway_clean, regex("\\(immune\\)", ignore_case = TRUE), "(IMMUNE)"),
#       pathway_clean = capitalize_acronyms(pathway_clean, acronyms),
#       pathway_clean = replace_mixed_case(pathway_clean, special_mixed, special_replacements),
#       pathway_clean = paste0(pathway_clean, " (", size, ")")
#     ) %>%
#     arrange(pval)
# 
#   fgsea_res$pathway_clean <- reorder(fgsea_res$pathway_clean, -abs(fgsea_res$NES))
# 
#   fgsea_res %>%
#     ggplot(aes(x = abs(NES), y = fct_rev(pathway_clean), label = pathway_clean)) +
#     geom_point(aes(size = -log10(pval), color = direction, alpha = 0.8)) +
#     # geom_vline(xintercept = 2, linetype = "dashed") +
#     geom_text(aes(group = pathway_clean, color = direction, fontface = face),
#               hjust = 0, size = text1, nudge_x = xnudge) +
#     scale_size_binned() +
#     scale_color_manual(values = c("Positive" = "#963b32", "Negative" = "#1f566f",
#                                   "Positive p > 0.05" = "#f1bcb6", "Negative p > 0.05" = "#b3d6e5")) +
#     scale_x_continuous(limits = c(0, xlimit), expand = expansion(mult = c(0, 0))) +
#     scale_y_discrete(expand = expansion(add = 1)) +
#     labs(
#       x = "NES",
#       y = "Pathways",
#       color = "Direction",
#       size = "-log(p-value)",
#       title = title
#     ) +
#     guides(alpha = "none") +
#     theme_bw() +
#     theme(
#       axis.text.y = element_blank(),
#       panel.grid = element_blank(),
#       axis.text.x = element_text(size = text3),
#       axis.title = element_text(size = text3),
#       axis.ticks.y = element_blank(),
#       legend.position = c(0.15, 0.5),
#       legend.background = element_blank(),
#       legend.box.background = element_rect(color = "black"),
#       legend.title = element_text(size = text2),
#       legend.text = element_text(size = text2),
#       title = element_text(size = text3)
#     )
# }
# 
#Cell Chat Functions----
filter_by_celltype_size <- function(seurat_obj, celltype_col = "celltype", min_cells = 50) {
  # Get celltype
  celltypes <- seurat_obj@meta.data[[celltype_col]]

  # Count
  cell_counts <- table(celltypes)
  keep_celltypes <- names(cell_counts[cell_counts >= min_cells])

  # Subset
  keep_cells <- rownames(seurat_obj@meta.data)[celltypes %in% keep_celltypes]
  seurat_obj_filtered <- subset(seurat_obj, cells = keep_cells)

  return(seurat_obj_filtered)
}


run_cellchat_pipeline <- function(seurat_obj, group_name, group.by = "celltype", assay = "RNA", min_cells = 10) {
  # future::plan("multisession", workers = n_cores)

  cellchat <- createCellChat(object = seurat_obj, group.by = group.by, assay = assay)
  cellchat <- addMeta(cellchat, meta = seurat_obj@meta.data)
  cellchat <- setIdent(cellchat, ident.use = group.by)

  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  saveRDS(cellchat, file = paste0("cellchat_", group_name, ".rds"))
  return(cellchat)
}


get_group_size <- function(cellchat_obj) {
  return(as.numeric(table(cellchat_obj@idents)))
}

