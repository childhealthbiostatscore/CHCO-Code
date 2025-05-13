#Gene Sets ----
# Hemodynamics
# Sodium and fractional sodium excretion (sodium transporter family)
# Regulation of afferent/efferent dilation (KCNMA1, KCNN family, KCNJ family, RYR family)
# Auto-/paracrine modulation of vascular resistance (COX family, NOS family, PTGER family, PTGES3)
# Vascular remodeling/rarefaction: ICAM-1, VCAM1, NOS3, SEL family, ANGPT family, VWF, CD31, EMCN, KDR, PLVAP, VEGFA

# Inflammation
# NF-κβ: Il-6, IL-8, CXCL-1, CXCL2, CCL20, ICAM-1, TNFAIP3, NaF-1, TRAF-1 and TNF
# JAK/STAT: Fas, Bcl-2, Bcl-X, p21, p27
# Systemic inflammation: TNF and IL family
# Intra-renal inflammation: Complement family, BMP2, immune signatures CD45 CCL2

# Oxidative stress
# ROS regulation: SOD, NOX family, ERO1, ACOX, DAO, XO, SOD1, HMOX1, KEAP-1, NRF2, FOXO family
# Hypoxia: EPO, PHD2, HIF1A, VEGFA, glycolytic genes
# Lipid peroxidation
# Fibrosis
# Pro-fibrotic programs: TGFb family, PDGF family, CTGF and associated genes, MAPK, Wnt/B-catenin genes, STAT family, MYK genes, JAK family, Notch family, PI3K/AKT/PKB genes
# ECM genes: collagens, laminins, FN, VIM, ACTA2, thrombospondins, fibrillins, transglutaminases
# ECM degradation: MMPs, ADAMTS family
# Kidney injury
# Kidney function: LRP2, PPARA, PPARGC1A, SLCs, NPHS1, SYNPO, WT1
# Kidney injury: CSTC3, HAVCR1, LCN2

#DEGS ----
# de.markers(so_liver_sn, genes, "diagnosis_of_diabetes", id2 = "No", id1 = "Yes", NULL, "_top")
# seurat_object <- so_kidney_sc
# group.by <- "epic_glp1ra_1"
# id1 <- "No"
# id2 <- "Yes"
# celltype="PT"
# genes <- sens_genes
# extension="_top"
# function for de.markers
de.markers <- function(seurat_object, genes, group.by, id1, id2, celltype, extension,logfc.threshold,min.pct){
  m = FindMarkers(seurat_object, features = genes,group.by = group.by,ident.1 = id1,
                  ident.2 = id2, subset.ident = celltype, verbose = F)
  m$p_val_adj = p.adjust(m$p_val,method = "fdr")
  # m$fdr = p.adjust(m$p_val,method = "fdr")
  # m$p_val_adj = m$p_val
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

de.markers2 <- function(seurat_object, genes, group.by, id1, id2, celltype, extension, logfc.threshold, min.pct) {
  # Run FindMarkers to find differentially expressed genes
  m = FindMarkers(seurat_object, features = genes, group.by = group.by, ident.1 = id1, 
                  ident.2 = id2, subset.ident = celltype, verbose = F)
  
  # Adjust p-values using FDR
  m$fdr = p.adjust(m$p_val, method = "fdr")
  m$p_val_adj = m$p_val
  
  # Add columns for mean expression per group
  m$mean_expr_id1 = rowMeans(seurat_object@assays$RNA@data[genes, seurat_object$seurat_clusters == id1])
  m$mean_expr_id2 = rowMeans(seurat_object@assays$RNA@data[genes, seurat_object$seurat_clusters == id2])
  
  # Clean and select relevant columns
  m <- m %>%
    rownames_to_column('gene') %>%
    arrange(p_val) %>%
    column_to_rownames('gene') %>%
    dplyr::select(avg_log2FC, pct.1, pct.2, p_val, p_val_adj, mean_expr_id1, mean_expr_id2) %>%
    filter(!is.na(p_val))
  
  # Identify significant genes
  genes_subset <- rownames(m)[m$p_val <= 0.05]
  
  # Store significant genes in global environment if any
  if (length(genes_subset) > 0) {
    assign(paste0("genes_subset", extension), genes_subset, envir = .GlobalEnv)
  }
  
  # Assign results to global environment
  assign(paste0("m", extension), m, envir = .GlobalEnv)
  
  # Return the table of results
  return(knitr::kable(m, digits = 3))
}

#DEG & GSEA Function ----
#a. All Genes ----
degs_fxn <- function(so,cell,exposure,gene_set,additional_group,exp_group,ref_group,enrichment,top_gsea) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  if(is.null(additional_group)) {
    condition <- paste0(str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  } 
  if(!is.null(additional_group)) {
    condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  }
  #Differential Expression by Group
  if (!is.null(cell)) {
    cell_name <- str_replace_all(cell,"/","_")
  }
  
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
  
  if(is.null(additional_group)) {
    condition2 <- paste0(str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0("DEGs in ",cell_name," cells for ",condition2)
    } else {
      title <- paste0("Bulk DEGs for ",condition2)
    }
  }
  if(!is.null(additional_group)) {
    condition2 <- paste0(additional_group," ",str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0("DEGs in ",cell_name," cells for ",condition2)
    } else {
      title <- paste0("Bulk DEGs for ",condition2)
    }
  } 
  
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  total_cells <- ncol(so)  # Total number of cells
  total_genes <- nrow(m)  # Total number of genes in m_top
  # Create subtitle with additional information
  subtitle_text <- paste0(
    "Positive Log2 FC = Greater Expression in ", condition2, "\n",
    "(Significant at FDR-P<0.05, FC Threshold = 0.5, Min Gene Exp 5% Threshold) \n",
    "Total Cells = ", total_cells, " | Total Genes = ", total_genes
  )
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       subtitle = subtitle_text,
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
  
  
  if(!is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(additional_group,"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(additional_group,"_Bulk_DEGs_for_",condition,".pdf") 
    }
  }
  if(is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0("DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0("Bulk_DEGs_for_",condition,".pdf") 
    }
  }
  pdf(fs::path(dir.results,filename),width=20,height=10)
  plot(p)
  dev.off()
  
  # #GSEA
  # if (enrichment=="Yes") {
  #   #Gene set enrichment analysis
  #   gc()
  #   sce_sn_hep <- as.SingleCellExperiment(so)
  #   ## C2 category is according to canonical pathways: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
  #   geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  #   ### filter background to only include genes that we assessed
  #   geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
  #   geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_sn_hep),]
  #   m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
  #   stats <- m$p_val_adj
  #   names(stats) <- rownames(m)
  #   eaRes <- fgsea(pathways = m_list, stats = na.omit(stats))
  #   #ooEA <- order(eaRes$pval, decreasing = FALSE)
  #   #kable(head(eaRes[ooEA, 1:7], n = 20))
  #   # Convert the leadingEdge column to comma-separated strings
  #   eaRes$leadingEdge <- sapply(eaRes$leadingEdge, function(x) paste(x, collapse = ", "))
  #   gc()
  #   #Significant pathways
  #   sig <- eaRes %>% 
  #     filter(padj<0.05)
  #   
  #   #Plot top pathways
  #   # Subset top pathways for visualization
  #   top_pathways <- eaRes[order(-abs(eaRes$NES)), ][1:top_gsea, ]  # Top 10 pathways based on adjusted p-value
  #   
  #   # Define a significance threshold
  #   significance_threshold <- 0.05
  #   
  #   # Add a significance column for coloring
  #   top_pathways$significance <- ifelse(
  #     top_pathways$padj > significance_threshold, "Non-significant",
  #     ifelse(top_pathways$NES > 0, "Positive Significant", "Negative Significant")
  #   )
  #   
  #   if (!is.null(cell)) {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition2)
  #   } else {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) for ",condition2," (Bulk)")
  #   }
  #   # Create a bar plot
  #   gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
  #     geom_bar(stat = "identity") +
  #     coord_flip() +  # Flip coordinates for better readability
  #     scale_fill_manual(
  #       values = c(
  #         "Positive Significant" = "red",
  #         "Negative Significant" = "blue",
  #         "Non-significant" = "gray"
  #       ),
  #       name = "Significance"
  #     ) +
  #     labs(
  #       title = title,
  #       x = "Pathway",
  #       y = "Normalized Enrichment Score (NES)"
  #     ) +
  #     theme_minimal(base_size = 12) +
  #     theme(
  #       axis.text.y = element_text(size = 10),
  #       axis.title.y = element_text(size = 12),
  #       plot.title = element_text(size = 14, face = "bold")
  #     )
  #   
  #   
  #   if (!is.null(cell)){
  #     filename <- paste0(additional_group,"_GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
  #   } else {
  #     filename <- paste0(additional_group,"_Bulk_GSEA_for_",condition,".pdf") 
  #   }
  #   pdf(fs::path(dir.results,filename),width=20,height=20)
  #   plot(gsea_plot)
  #   dev.off()
  # }
  #Make sure gene nemaes are in printed file
  deg_results <- m 
  deg_results$Gene <- rownames(deg_results)
  
  # #save results to excel file
  # write_multiple_sheets <- function(output_file, sheet_data) {
  #   # Create a new workbook
  #   wb <- createWorkbook()
  #   
  #   # Loop over the list of data frames
  #   for (sheet_name in names(sheet_data)) {
  #     # Add a new sheet to the workbook
  #     addWorksheet(wb, sheet_name)
  #     
  #     # Write the data frame to the sheet
  #     writeData(wb, sheet_name, sheet_data[[sheet_name]])
  #   }
  #   # Save the workbook to the specified file
  #   saveWorkbook(wb, file = output_file, overwrite = TRUE)
  # } 
  
  # if (enrichment=="Yes") {
  #   # Example usage
  #   df1 <- deg_results
  #   df2 <- eaRes
  # } else {
  #   df1 <- deg_results
  #   df2 <- NULL
  # }
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
  } else {
    output_file <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".csv"))
  }
  write.csv(deg_results,output_file)
  
  # sheet_data <- list(
  #   "DEG" = df1,
  #   "Pathway_Results" = df2
  # )
  # 
  # # Call the function
  # write_multiple_sheets(output_file, sheet_data)
  # 
  #Save for IPA
  # Call the function
  ipa_results <- deg_results %>% 
    dplyr::select(all_of(c("Gene","avg_log2FC","p_val","p_val_adj"))) %>% 
    dplyr::rename(ID=Gene,
                  `Fold Change`=avg_log2FC,
                  `P-Value`=p_val,
                  `Adjusted P-Value`=p_val_adj)
  if (!is.null(cell)){
    output_file2 <- fs::path(dir.ipa,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
  } else {
    output_file2 <- fs::path(dir.ipa,paste0("Bulk_Results_for_",condition,".csv"))
  }
  write.csv(ipa_results, output_file2)
  return(p)
}

#b. Pseudobulk by Individual
degs_fxn2 <- function(so,cell,exposure,gene_set,additional_group,exp_group,ref_group,enrichment,top_gsea,cellname) {
  pathway=NULL
  DefaultAssay(so) <- "RNA"

  #Conidtion names
  if(is.null(additional_group)) {
    condition <- paste0(str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  }
  if(!is.null(additional_group)) {
    condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  }
  #Differential Expression by Group
  if (!is.null(cell)) {
    cell_name <- str_replace_all(cell,"/","_")
  }

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

  if(is.null(additional_group)) {
    condition2 <- paste0(str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0("DEGs in ",cell_name," cells for ",condition2)
    } else {
      title <- paste0(cellname,"Cells, Individual PseudoBulk DEGs for ",condition2)
    }
  }
  if(!is.null(additional_group)) {
    condition2 <- paste0(additional_group," ",str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0(pathway," DEGs in ",cell_name," cells for ",condition2)
    } else {
      title <- paste0(pathway," ",cellname,"Cells, Individual Pseudobulk DEGs for ",condition2)
    }
  }


  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  total_individuals <- ncol(so)  # Total number of cells
  total_genes <- nrow(m)  # Total number of genes in m_top
  # Create subtitle with additional information
  subtitle_text <- paste0(
    "Positive Log2 FC = Greater Expression in ", condition2, "\n",
    "(Significant at FDR-P<0.05, FC Threshold = 0.5, Min Gene Exp 5% Threshold) \n",
    "N = ", total_individuals, " | Total Genes = ", total_genes
  )

  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       subtitle = subtitle_text,
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


  if(!is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(pathway,"_",additional_group,"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(pathway,"_",additional_group,"_Bulk_DEGs_for_",condition,".pdf")
    }
  }
  if(is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(pathway,"_","DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(pathway,"_","Pseudobulk_DEGs_for_",condition,".pdf")
    }
  }
  pdf(fs::path(dir.results,filename),width=20,height=10)
  plot(p)
  dev.off()

  # #GSEA
  # if (enrichment=="Yes") {
  #   #Gene set enrichment analysis
  #   gc()
  #   sce_sn_hep <- as.SingleCellExperiment(so)
  #   ## C2 category is according to canonical pathways: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
  #   geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  #   ### filter background to only include genes that we assessed
  #   geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
  #   geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_sn_hep),]
  #   m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
  #   stats <- m$p_val_adj
  #   names(stats) <- rownames(m)
  #   eaRes <- fgsea(pathways = m_list, stats = na.omit(stats))
  #   #ooEA <- order(eaRes$pval, decreasing = FALSE)
  #   #kable(head(eaRes[ooEA, 1:7], n = 20))
  #   # Convert the leadingEdge column to comma-separated strings
  #   eaRes$leadingEdge <- sapply(eaRes$leadingEdge, function(x) paste(x, collapse = ", "))
  #   gc()
  #   #Significant pathways
  #   sig <- eaRes %>%
  #     filter(padj<0.05)
  #
  #   #Plot top pathways
  #   # Subset top pathways for visualization
  #   top_pathways <- eaRes[order(-abs(eaRes$NES)), ][1:top_gsea, ]  # Top 10 pathways based on adjusted p-value
  #
  #   # Define a significance threshold
  #   significance_threshold <- 0.05
  #
  #   # Add a significance column for coloring
  #   top_pathways$significance <- ifelse(
  #     top_pathways$padj > significance_threshold, "Non-significant",
  #     ifelse(top_pathways$NES > 0, "Positive Significant", "Negative Significant")
  #   )
  #
  #   if (!is.null(cell)) {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition2)
  #   } else {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) for ",condition2," (Bulk)")
  #   }
  #   # Create a bar plot
  #   gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
  #     geom_bar(stat = "identity") +
  #     coord_flip() +  # Flip coordinates for better readability
  #     scale_fill_manual(
  #       values = c(
  #         "Positive Significant" = "red",
  #         "Negative Significant" = "blue",
  #         "Non-significant" = "gray"
  #       ),
  #       name = "Significance"
  #     ) +
  #     labs(
  #       title = title,
  #       x = "Pathway",
  #       y = "Normalized Enrichment Score (NES)"
  #     ) +
  #     theme_minimal(base_size = 12) +
  #     theme(
  #       axis.text.y = element_text(size = 10),
  #       axis.title.y = element_text(size = 12),
  #       plot.title = element_text(size = 14, face = "bold")
  #     )
  #
  #
  #   if (!is.null(cell)){
  #     filename <- paste0(additional_group,"_GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
  #   } else {
  #     filename <- paste0(additional_group,"_Bulk_GSEA_for_",condition,".pdf")
  #   }
  #   pdf(fs::path(dir.results,filename),width=20,height=20)
  #   plot(gsea_plot)
  #   dev.off()
  # }
  #Make sure gene nemaes are in printed file
  deg_results <- m
  deg_results$Gene <- rownames(deg_results)

  # #save results to excel file
  # write_multiple_sheets <- function(output_file, sheet_data) {
  #   # Create a new workbook
  #   wb <- createWorkbook()
  #
  #   # Loop over the list of data frames
  #   for (sheet_name in names(sheet_data)) {
  #     # Add a new sheet to the workbook
  #     addWorksheet(wb, sheet_name)
  #
  #     # Write the data frame to the sheet
  #     writeData(wb, sheet_name, sheet_data[[sheet_name]])
  #   }
  #   # Save the workbook to the specified file
  #   saveWorkbook(wb, file = output_file, overwrite = TRUE)
  # }

  # if (enrichment=="Yes") {
  #   # Example usage
  #   df1 <- deg_results
  #   df2 <- eaRes
  # } else {
  #   df1 <- deg_results
  #   df2 <- NULL
  # }

  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))
  } else {
    output_file <- fs::path(dir.results,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
  }
  write.csv(deg_results,output_file)

  # sheet_data <- list(
  #   "DEG" = df1,
  #   "Pathway_Results" = df2
  # )
  #
  # # Call the function
  # write_multiple_sheets(output_file, sheet_data)
  #
  #Save for IPA
  # Call the function
  ipa_results <- deg_results %>%
    dplyr::select(all_of(c("Gene","avg_log2FC","p_val","p_val_adj"))) %>%
    dplyr::rename(ID=Gene,
                  `Fold Change`=avg_log2FC,
                  `P-Value`=p_val,
                  `Adjusted P-Value`=p_val_adj)
  if (!is.null(cell)){
    output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))
    # output_file2 <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))

  } else {
    output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
    # output_file2 <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".csv"))

  }
  write.csv(ipa_results, output_file2)
  return(p)
}


#c. Senescence Genes ----
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
    filename <- paste0("Senescence DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
  } else {
    filename <- paste0("Senescence Bulk_DEGs_for_",condition,".pdf") 
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
  df2 <- eaRes
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0("Senescence Results_",cell_name,"_cells_for_",condition,".xlsx"))  
  } else {
    output_file <- fs::path(dir.results,paste0("Senescence Bulk_Results_for_",condition,".xlsx"))
  }
  
  sheet_data <- list(
    "DEG" = df1,
    "Pathway_Results" = df2
  )
  
  # Call the function
  write_multiple_sheets(output_file, sheet_data)
}

#c. Pathway Genes ----
degs_fxn3 <- function(so,cell,exposure,gene_set,additional_group,exp_group,ref_group,enrichment,top_gsea,cellname,pathway) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  if(is.null(additional_group)) {
    condition <- paste0(str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  } 
  if(!is.null(additional_group)) {
    condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  }
  #Differential Expression by Group
  if (!is.null(cell)) {
    cell_name <- str_replace_all(cell,"/","_")
  }
  
  # # Scale the data
  # so <- NormalizeData(so)
  # so <- ScaleData(so)
  # 
  # # Perform PCA for dimensionality reduction
  # so <- RunPCA(so)
  # 
  # # Find clusters
  # so <- FindNeighbors(so, dims = 1:10)
  # so <- FindClusters(so, resolution = 0.5)
  
  de.markers(so, gene_set, exposure, id2 = ref_group, id1 = exp_group, cell,"")
  dp <- dp.formatted(so,genes=gene_set,group.by=exposure,m=m,celltype=NULL,exp_group=exp_group,ref_group=ref_group)
  colnames(m)[2] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group,")")
  colnames(m)[3] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",ref_group,")")
  
  
  # print(dp)
  
  # dp <- dp.formatted(so,genes=gene_set,group.by=exposure,m=m,celltype=NULL)
  # return(dp)
  # # Create the DotPlot
  # # p <- DotPlot(object = so, features = gene_set) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # # Assuming 'm' contains the results from de.markers() with p-values
  # # Filter for significant genes (adjusted p-value < 0.05)
  # significant_genes <- m[m$p_val_adj < 0.5, ]
  # 
  # # Extract gene names that are significant
  # significant_gene_names <- rownames(significant_genes)
  # 
  # # Create the DotPlot for the gene set
  # p <- DotPlot(object = so, features = gene_set) + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # 
  # # Add '*' annotations for significant genes
  # plot_new <-p + geom_text(
  #   data = subset(p$data, features.plot %in% significant_gene_names),  # Only annotate significant genes
  #   aes(label = "*"),  # Display '*' for significant genes
  #   size = 4,  # Adjust text size
  #   color = "black",  # Adjust text color if necessary
  #   vjust = 0.5  # Adjust vertical position of text
  # )
  
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
  
  if(is.null(additional_group)) {
    condition2 <- paste0(str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0("DEGs in ",cell_name," cells \n ",condition2)
    } else {
      title <- paste0(cellname,"Cells, \n Individual PseudoBulk DEGs for ",condition2)
    }
  }
  if(!is.null(additional_group)) {
    condition2 <- paste0(additional_group," ",str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0(pathway," DEGs in ",cell_name," cells \n ",condition2)
    } else {
      title <- paste0(pathway," ",cellname,"Cells, \n Individual Pseudobulk DEGs for ",condition2)
    }
  } 
  
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  total_individuals <- ncol(so)  # Total number of cells
  total_genes <- nrow(m)  # Total number of genes in m_top
  # Create subtitle with additional information
  subtitle_text <- paste0("Positive Log2 FC = Greater Expression \n in ", condition2)
  caption <- paste0("Significant FDR-P<0.05, N = ", total_individuals, " | Total Genes = ", total_genes)
  
  # Create a new column for significance (significant vs non-significant)
  # Create a new column for significance with direction (positive, negative, or non-significant)
  m_top$Significance <- ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "Significant Positive", 
                               ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "Significant Negative", 
                                      "Non-Significant"))
  m_top$Gene <- rownames(m_top)

  p_dot <- ggplot(m_top, aes(x = "Cell Type", y = Gene)) +
    # Color points based on significance
    geom_point(aes(color = Significance, size = abs(avg_log2FC)), alpha = 0.7) +
    scale_color_manual(values = c("Significant Positive" = "red",
                                  "Significant Negative" = "blue",
                                  "Non-Significant" = "gray")) +  # Set colors for positive, negative, and non-significant
    scale_size_continuous(range = c(2, 8), name = "Magnitude (|Log2 Fold Change|)") +
    labs(
      title = title,
      subtitle = subtitle_text,
      x = cellname,  # Display a generic x-axis label
      y = "Gene",
      caption = caption
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      axis.text.y = element_text(size = 6),  # To make the gene names readable
      axis.title.x = element_text(size=8),
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.ticks.x = element_blank()  # Remove x-axis ticks
    )
    # coord_fixed(ratio = 0.05)
  
  # Assuming 'm' is a matrix of gene expression data, where rows are genes and columns are cells.
  # Calculate mean expression for each gene in each group (ref_group and exp_group)
  # ref2 <- colnames(m)[3]
  # exp2 <- colnames(m)[2]
  # mean_exp_ref <- rowMeans(m[ref2])  # Mean expression for ref_group
  # mean_exp_exp <- rowMeans(m[exp2])  # Mean expression for exp_group
  # 
  # # Combine mean expressions into a data frame
  # m_top2 <- data.frame(
  #   Gene = rownames(m),
  #   Mean_Expression_Ref = mean_exp_ref,
  #   Mean_Expression_Exp = mean_exp_exp,
  #   Significance = ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "Significant Positive", 
  #                         ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "Significant Negative", 
  #                                "Non-Significant"))
  # )
  # 
  # # Create subtitle with additional information
  # subtitle_text <- paste0("Positive Log2 FC = Greater Expression \n in ", exp_group)
  # caption <- paste0("Significant FDR-P<0.05, N = ", total_individuals, " | Total Genes = ", total_genes)
  # 
  # # Reshape the data to make it suitable for ggplot (long format for groups)
  # m_top_long <- data.frame(
  #   Gene = rep(rownames(m_top2), 2),
  #   Group = rep(c(exp_group, ref_group), each = nrow(m_top2)),
  #   Mean_Expression = c(m_top2$Mean_Expression_Exp, m_top2$Mean_Expression_Ref),
  #   Significance = rep(m_top2$Significance, 2)
  # )
  # 
  # # p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
  # #   # Plot points for each gene's mean expression in both groups
  # #   geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7) + 
  # #   # Apply a more granular diverging color scale with fewer breaks in the legend
  # #   scale_color_gradient2(low = "white", mid = "red1", high = "darkred", 
  # #                         midpoint = 0.7, 
  # #                         limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
  # #                         breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 5),  # Reduce the number of breaks
  # #                         labels = scales::number_format(accuracy = 0.01),  # Round labels to 2 decimal places
  # #                         name = "% Expression") + 
  # #   scale_size_continuous(range = c(2, 8), name = "% Expression") +  # Size of points based on mean expression
  # #   # Draw boxes around Group 1 and Group 2
  # #   # geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = NA, color = "black", size = 1.5) + # Box around Group 1
  # #   # geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = NA, color = "black", size = 1.5) + # Box around Group 2
  # #   labs(
  # #     title = title,
  # #     subtitle = subtitle_text,
  # #     x = "Group",
  # #     y = "Gene",
  # #     caption = caption
  # #   ) + 
  # #   theme_classic() +
  # #   theme(
  # #     plot.title = element_text(size = 10, face = "bold"),
  # #     plot.subtitle = element_text(size = 8),
  # #     axis.text = element_text(size = 6),
  # #     axis.title = element_text(size = 8),
  # #     legend.title = element_text(size = 8),
  # #     legend.text = element_text(size = 6),
  # #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  # #     axis.title.x = element_text(size = 8),  
  # #     axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
  # #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  # #   )
  # # p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
  # #   # Plot points for each gene's mean expression in both groups (just for background grid, no shapes)
  # #   geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7, show.legend = FALSE) + 
  # #   
  # #   # Apply a more granular diverging color scale with fewer breaks in the legend
  # #   scale_color_gradient2(low = "white", mid = "red1", high = "darkred", 
  # #                         midpoint = 0.7, 
  # #                         limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
  # #                         breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 5),  # Reduce the number of breaks
  # #                         labels = scales::number_format(accuracy = 0.01),  # Round labels to 2 decimal places
  # #                         name = "% Expression") + 
  # #   
  # #   scale_size_continuous(range = c(2, 8), name = "% Expression") +  # Size of points based on mean expression
  # #   
  # #   # Add asterisks (*) for significant genes (both positive and negative)
  # #   geom_text(data = subset(m_top_long, Significance != "Non-Significant"),
  # #             aes(label = "*"),
  # #             size = 4, color = "black", vjust = 0.5) +  # Center the asterisk on the dot
  # #   
  # #   labs(
  # #     title = title,
  # #     subtitle = subtitle_text,
  # #     x = "Group",
  # #     y = "Gene",
  # #     caption = caption
  # #   ) + 
  # #   theme_classic() +
  # #   theme(
  # #     plot.title = element_text(size = 10, face = "bold"),
  # #     plot.subtitle = element_text(size = 8),
  # #     axis.text = element_text(size = 6),
  # #     axis.title = element_text(size = 8),
  # #     legend.title = element_text(size = 8),
  # #     legend.text = element_text(size = 6),
  # #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  # #     axis.title.x = element_text(size = 8),  
  # #     axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
  # #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  # #   )
  # 
  # # Scale the Mean_Expression values (e.g., z-score scaling)
  # m_top_long$Scaled_Mean_Expression <- scale(m_top_long$Mean_Expression)
  #                                            
  # p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
  #   # Plot points for each gene's mean expression in both groups (just for background grid, no shapes)
  #   geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7, show.legend = TRUE) +  # Keep the legend for size and color
  #   
  #   # Apply a more granular diverging color scale with fewer breaks in the legend
  #   scale_color_gradient2(low = "white", mid = "red1", high = "darkred", 
  #                         midpoint = 0.7, 
  #                         limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
  #                         breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 5),  # Reduce the number of breaks
  #                         labels = scales::number_format(accuracy = 0.01),  # Round labels to 2 decimal places
  #                         name = "% Expression") + 
  #   
  #   scale_size_continuous(range = c(2, 8), name = "% Expression", guide = guide_legend(title = "% Expression")) +  # Size of points based on mean expression, ensure size legend
  #   
  #   # Add asterisks (*) for significant genes (both positive and negative)
  #   geom_text(data = subset(m_top_long, Significance != "Non-Significant"),
  #             aes(label = "*"),
  #             size = 4, color = "black", vjust = 0.5) +  # Center the asterisk on the dot
  #   
  #   labs(
  #     title = title,
  #     subtitle = subtitle_text,
  #     x = "Group",
  #     y = "Gene",
  #     caption = caption
  #   ) + 
  #   theme_classic() +
  #   theme(
  #     plot.title = element_text(size = 10, face = "bold"),
  #     plot.subtitle = element_text(size = 8),
  #     axis.text = element_text(size = 6),
  #     axis.title = element_text(size = 8),
  #     legend.title = element_text(size = 8),
  #     legend.text = element_text(size = 6),
  #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  #     axis.title.x = element_text(size = 8),  
  #     axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
  #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  #   )
  # 
  # 
  # # Plot the mean expression for each gene in the two groups (side by side)
  # p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
  #   # Plot points for each gene's mean expression in both groups
  #   geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7) +
  #   # geom_point(aes(color = Significance, size = Mean_Expression), alpha = 0.7) + 
  #   # scale_color_manual(values = c("Significant Positive" = "red", 
  #   #                               "Significant Negative" = "blue", 
  #   #                               "Non-Significant" = "gray")) +  # Set colors for positive, negative, and non-significant
  #   # scale_color_gradient(low = "white", high = "darkred", name = "% Expression") +
  #   # scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "% Expression") +
  #   scale_color_gradient2(low = "white", mid = "orangered", high = "darkred", 
  #                         midpoint = 0.5, limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
  #                         breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 10),
  #                         name = "% Expression") + 
  #   scale_size_continuous(range = c(2, 8), name = "% Expression") +
  #     # Draw boxes around Group 1 and Group 2
  #   # geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = NA, color = "black", size = 1.5) + # Box around Group 1
  #   # geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = NA, color = "black", size = 1.5) + # Box around Group 2
  #   labs(
  #     title = title,
  #     subtitle = subtitle_text,
  #     x = "Group",
  #     y = "Gene",
  #     caption = caption
  #   ) + 
  #   theme_classic() +
  #   theme(
  #     plot.title = element_text(size = 10, face = "bold"),
  #     plot.subtitle = element_text(size = 8),
  #     axis.text = element_text(size = 6),
  #     axis.title = element_text(size = 8),
  #     legend.title = element_text(size = 8),
  #     legend.text = element_text(size = 6),
  #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  #     axis.title.x = element_text(size=8),  
  #     axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
  #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  #   )
  
  
  # p <- EnhancedVolcano(m_top,
  #                      lab = labels,
  #                      x = 'avg_log2FC',
  #                      y = 'p_val_adj',
  #                      title = title,
  #                      subtitle = subtitle_text,
  #                      pCutoff = 0.05,
  #                      FCcutoff = 0.5,
  #                      labFace = 'bold',
  #                      pointSize = 4,
  #                      labSize = 5,
  #                      drawConnectors = TRUE,
  #                      widthConnectors = 1.0,
  #                      colConnectors = 'black',
  #                      legendPosition=NULL,
  #                      boxedLabels = TRUE,
  #                      max.overlaps=60)
  # 
  
  if(!is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(pathway,"_",additional_group,"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(pathway,"_",additional_group,"_Bulk_DEGs_for_",condition,".pdf") 
    }
  }
  if(is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(pathway,"_","DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(pathway,"_","Pseudobulk_DEGs_for_",condition,".pdf") 
    }
  }
  pdf(fs::path(dir.results,filename),width=10,height=5)
  # plot(p_dot2)
  plot(dp)
  dev.off()
  
  # #GSEA
  # if (enrichment=="Yes") {
  #   #Gene set enrichment analysis
  #   gc()
  #   sce_sn_hep <- as.SingleCellExperiment(so)
  #   ## C2 category is according to canonical pathways: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
  #   geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  #   ### filter background to only include genes that we assessed
  #   geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
  #   geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_sn_hep),]
  #   m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
  #   stats <- m$p_val_adj
  #   names(stats) <- rownames(m)
  #   eaRes <- fgsea(pathways = m_list, stats = na.omit(stats))
  #   #ooEA <- order(eaRes$pval, decreasing = FALSE)
  #   #kable(head(eaRes[ooEA, 1:7], n = 20))
  #   # Convert the leadingEdge column to comma-separated strings
  #   eaRes$leadingEdge <- sapply(eaRes$leadingEdge, function(x) paste(x, collapse = ", "))
  #   gc()
  #   #Significant pathways
  #   sig <- eaRes %>% 
  #     filter(padj<0.05)
  #   
  #   #Plot top pathways
  #   # Subset top pathways for visualization
  #   top_pathways <- eaRes[order(-abs(eaRes$NES)), ][1:top_gsea, ]  # Top 10 pathways based on adjusted p-value
  #   
  #   # Define a significance threshold
  #   significance_threshold <- 0.05
  #   
  #   # Add a significance column for coloring
  #   top_pathways$significance <- ifelse(
  #     top_pathways$padj > significance_threshold, "Non-significant",
  #     ifelse(top_pathways$NES > 0, "Positive Significant", "Negative Significant")
  #   )
  #   
  #   if (!is.null(cell)) {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition2)
  #   } else {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) for ",condition2," (Bulk)")
  #   }
  #   # Create a bar plot
  #   gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
  #     geom_bar(stat = "identity") +
  #     coord_flip() +  # Flip coordinates for better readability
  #     scale_fill_manual(
  #       values = c(
  #         "Positive Significant" = "red",
  #         "Negative Significant" = "blue",
  #         "Non-significant" = "gray"
  #       ),
  #       name = "Significance"
  #     ) +
  #     labs(
  #       title = title,
  #       x = "Pathway",
  #       y = "Normalized Enrichment Score (NES)"
  #     ) +
  #     theme_minimal(base_size = 12) +
  #     theme(
  #       axis.text.y = element_text(size = 10),
  #       axis.title.y = element_text(size = 12),
  #       plot.title = element_text(size = 14, face = "bold")
  #     )
  #   
  #   
  #   if (!is.null(cell)){
  #     filename <- paste0(additional_group,"_GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
  #   } else {
  #     filename <- paste0(additional_group,"_Bulk_GSEA_for_",condition,".pdf") 
  #   }
  #   pdf(fs::path(dir.results,filename),width=20,height=20)
  #   plot(gsea_plot)
  #   dev.off()
  # }
  #Make sure gene nemaes are in printed file
  deg_results <- m 
  deg_results$Gene <- rownames(deg_results)
  
  # #save results to excel file
  # write_multiple_sheets <- function(output_file, sheet_data) {
  #   # Create a new workbook
  #   wb <- createWorkbook()
  #   
  #   # Loop over the list of data frames
  #   for (sheet_name in names(sheet_data)) {
  #     # Add a new sheet to the workbook
  #     addWorksheet(wb, sheet_name)
  #     
  #     # Write the data frame to the sheet
  #     writeData(wb, sheet_name, sheet_data[[sheet_name]])
  #   }
  #   # Save the workbook to the specified file
  #   saveWorkbook(wb, file = output_file, overwrite = TRUE)
  # } 
  
  # if (enrichment=="Yes") {
  #   # Example usage
  #   df1 <- deg_results
  #   df2 <- eaRes
  # } else {
  #   df1 <- deg_results
  #   df2 <- NULL
  # }
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))  
  } else {
    output_file <- fs::path(dir.results,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
  }
  write.csv(deg_results,output_file)
  
  # sheet_data <- list(
  #   "DEG" = df1,
  #   "Pathway_Results" = df2
  # )
  # 
  # # Call the function
  # write_multiple_sheets(output_file, sheet_data)
  # 
  #Save for IPA
  # Call the function
  ipa_results <- deg_results %>% 
    dplyr::select(all_of(c("Gene","avg_log2FC","p_val","p_val_adj"))) %>% 
    dplyr::rename(ID=Gene,
                  `Fold Change`=avg_log2FC,
                  `P-Value`=p_val,
                  `Adjusted P-Value`=p_val_adj)
  if (!is.null(cell)){
    output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))
    # output_file2 <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
    
  } else {
    output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
    # output_file2 <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".csv"))
    
  }
  write.csv(ipa_results, output_file2)
  # return(p)
  return(dp)
}

degs_fxn4 <- function(so,cell,exposure,gene_set,additional_group,exp_group,ref_group,enrichment,top_gsea,cellname,pathway) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  if(is.null(additional_group)) {
    condition <- paste0(str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  } 
  if(!is.null(additional_group)) {
    condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  }
  #Differential Expression by Group
  if (!is.null(cell)) {
    cell_name <- str_replace_all(cell,"/","_")
  }
  
  # # Scale the data
  # so <- NormalizeData(so)
  # so <- ScaleData(so)
  # 
  # # Perform PCA for dimensionality reduction
  # so <- RunPCA(so)
  # 
  # # Find clusters
  # so <- FindNeighbors(so, dims = 1:10)
  # so <- FindClusters(so, resolution = 0.5)
  
  de.markers(so, gene_set, exposure, id2 = ref_group, id1 = exp_group, cell,"")
  # colnames(m)[2] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group,")")
  # colnames(m)[3] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",ref_group,")")

  # dp <- dp.formatted(so,genes=gene_set,group.by=exposure,m=m,celltype=NULL)
  # print(dp)
  
  dp <- dp.formatted(so,genes=gene_set,group.by=exposure,m=m,celltype=NULL,exp_group=exp_group,ref_group=ref_group)
  
  # Create the DotPlot
  # p <- DotPlot(object = so, features = gene_set) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # Assuming 'm' contains the results from de.markers() with p-values
  # Filter for significant genes (adjusted p-value < 0.05)
  # significant_genes <- m[m$p_val_adj < 0.5, ]
  # 
  # # Extract gene names that are significant
  # significant_gene_names <- rownames(significant_genes)
  # 
  # Create the DotPlot for the gene set
  # p <- DotPlot(object = so, features = gene_set) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# return(p)
  # Add '*' annotations for significant genes
  # plot_new <-p + geom_text(
  #   data = subset(p$data, features.plot %in% significant_gene_names),  # Only annotate significant genes
  #   aes(label = "*"),  # Display '*' for significant genes
  #   size = 4,  # Adjust text size
  #   color = "black",  # Adjust text color if necessary
  #   vjust = 0.5  # Adjust vertical position of text
  # )
  # 
  # # Filter for significant genes
  # m_top <- m
  # 
  # significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  # 
  # # Select the top 10 positive and 10 negative log2FC genes based on the largest magnitude of fold change
  # top_positive_by_fc <- significant_genes %>% 
  #   filter(avg_log2FC > 0) %>% 
  #   arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
  #   head(10)  # Top 10 positive fold changes
  # 
  # top_negative_by_fc <- significant_genes %>% 
  #   filter(avg_log2FC < 0) %>% 
  #   arrange(desc(abs(avg_log2FC))) %>%  # Sort by absolute fold change (largest first)
  #   head(10)  # Top 10 negative fold changes
  # 
  # # Select the top 10 positive and 10 negative log2FC genes based on significance (p_val_adj)
  # top_positive_by_significance <- significant_genes %>% 
  #   filter(avg_log2FC > 0) %>% 
  #   arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
  #   head(10)  # Top 10 most significant positive fold changes
  # 
  # top_negative_by_significance <- significant_genes %>% 
  #   filter(avg_log2FC < 0) %>% 
  #   arrange(p_val_adj) %>%  # Sort by smallest p-value (most significant)
  #   head(10)  # Top 10 most significant negative fold changes
  # 
  # # Combine top fold-change based and significance-based genes into a final list
  # top_genes <- rbind(
  #   top_positive_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
  #   top_negative_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
  #   top_positive_by_significance %>% mutate(Selection = "Top 10 by Significance"),
  #   top_negative_by_significance %>% mutate(Selection = "Top 10 by Significance")
  # )
  # 
  # if(is.null(additional_group)) {
  #   condition2 <- paste0(str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
  #   if (!is.null(cell)){
  #     title <- paste0("DEGs in ",cell_name," cells \n ",condition2)
  #   } else {
  #     title <- paste0(cellname,"Cells, \n Individual PseudoBulk DEGs for ",condition2)
  #   }
  # }
  # if(!is.null(additional_group)) {
  #   condition2 <- paste0(additional_group," ",str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
  #   if (!is.null(cell)){
  #     title <- paste0(pathway," DEGs in ",cell_name," cells \n ",condition2)
  #   } else {
  #     title <- paste0(pathway," ",cellname,"Cells, \n Individual Pseudobulk DEGs for ",condition2)
  #   }
  # } 
  # 
  # 
  # labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  # total_individuals <- ncol(so)  # Total number of cells
  # total_genes <- nrow(m)  # Total number of genes in m_top
  # # Create subtitle with additional information
  # subtitle_text <- paste0("Positive Log2 FC = Greater Expression \n in ", condition2)
  # caption <- paste0("Significant FDR-P<0.05, N = ", total_individuals, " | Total Genes = ", total_genes)
  # 
  # # Create a new column for significance (significant vs non-significant)
  # # Create a new column for significance with direction (positive, negative, or non-significant)
  # m_top$Significance <- ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "Significant Positive", 
  #                              ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "Significant Negative", 
  #                                     "Non-Significant"))
  # m_top$Gene <- rownames(m_top)
  # 
  # p_dot <- ggplot(m_top, aes(x = "Cell Type", y = Gene)) +
  #   # Color points based on significance
  #   geom_point(aes(color = Significance, size = abs(avg_log2FC)), alpha = 0.7) +
  #   scale_color_manual(values = c("Significant Positive" = "red",
  #                                 "Significant Negative" = "blue",
  #                                 "Non-Significant" = "gray")) +  # Set colors for positive, negative, and non-significant
  #   scale_size_continuous(range = c(2, 8), name = "Magnitude (|Log2 Fold Change|)") +
  #   labs(
  #     title = title,
  #     subtitle = subtitle_text,
  #     x = cellname,  # Display a generic x-axis label
  #     y = "Gene",
  #     caption = caption
  #   ) +
  #   theme_classic() +
  #   theme(
  #     plot.title = element_text(size = 10, face = "bold"),
  #     plot.subtitle = element_text(size = 8),
  #     axis.text = element_text(size = 6),
  #     axis.title = element_text(size = 8),
  #     legend.title = element_text(size = 8),
  #     legend.text = element_text(size = 6),
  #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  #     axis.title.x = element_text(size=8),
  #     axis.text.x = element_blank(),  # Remove x-axis text
  #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  #   )
  # # coord_fixed(ratio = 0.05)
  # 
  # # Assuming 'm' is a matrix of gene expression data, where rows are genes and columns are cells.
  # # Calculate mean expression for each gene in each group (ref_group and exp_group)
  # ref2 <- colnames(m)[3]
  # exp2 <- colnames(m)[2]
  # mean_exp_ref <- rowMeans(m[ref2])  # Mean expression for ref_group
  # mean_exp_exp <- rowMeans(m[exp2])  # Mean expression for exp_group
  # 
  # # Combine mean expressions into a data frame
  # m_top2 <- data.frame(
  #   Gene = rownames(m),
  #   Mean_Expression_Ref = mean_exp_ref,
  #   Mean_Expression_Exp = mean_exp_exp,
  #   Significance = ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "Significant Positive", 
  #                         ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "Significant Negative", 
  #                                "Non-Significant"))
  # )
  # 
  # # Create subtitle with additional information
  # subtitle_text <- paste0("Positive Log2 FC = Greater Expression \n in ", exp_group)
  # caption <- paste0("Significant FDR-P<0.05, N = ", total_individuals, " | Total Genes = ", total_genes)
  # 
  # # Reshape the data to make it suitable for ggplot (long format for groups)
  # m_top_long <- data.frame(
  #   Gene = rep(rownames(m_top2), 2),
  #   Group = rep(c(exp_group, ref_group), each = nrow(m_top2)),
  #   Mean_Expression = c(m_top2$Mean_Expression_Exp, m_top2$Mean_Expression_Ref),
  #   Significance = rep(m_top2$Significance, 2)
  # )
  # 
  # # p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
  # #   # Plot points for each gene's mean expression in both groups
  # #   geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7) + 
  # #   # Apply a more granular diverging color scale with fewer breaks in the legend
  # #   scale_color_gradient2(low = "white", mid = "red1", high = "darkred", 
  # #                         midpoint = 0.7, 
  # #                         limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
  # #                         breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 5),  # Reduce the number of breaks
  # #                         labels = scales::number_format(accuracy = 0.01),  # Round labels to 2 decimal places
  # #                         name = "% Expression") + 
  # #   scale_size_continuous(range = c(2, 8), name = "% Expression") +  # Size of points based on mean expression
  # #   # Draw boxes around Group 1 and Group 2
  # #   # geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = NA, color = "black", size = 1.5) + # Box around Group 1
  # #   # geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = NA, color = "black", size = 1.5) + # Box around Group 2
  # #   labs(
  # #     title = title,
  # #     subtitle = subtitle_text,
  # #     x = "Group",
  # #     y = "Gene",
  # #     caption = caption
  # #   ) + 
  # #   theme_classic() +
  # #   theme(
  # #     plot.title = element_text(size = 10, face = "bold"),
  # #     plot.subtitle = element_text(size = 8),
  # #     axis.text = element_text(size = 6),
  # #     axis.title = element_text(size = 8),
  # #     legend.title = element_text(size = 8),
  # #     legend.text = element_text(size = 6),
  # #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  # #     axis.title.x = element_text(size = 8),  
  # #     axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
  # #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  # #   )
  # # p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
  # #   # Plot points for each gene's mean expression in both groups (just for background grid, no shapes)
  # #   geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7, show.legend = FALSE) + 
  # #   
  # #   # Apply a more granular diverging color scale with fewer breaks in the legend
  # #   scale_color_gradient2(low = "white", mid = "red1", high = "darkred", 
  # #                         midpoint = 0.7, 
  # #                         limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
  # #                         breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 5),  # Reduce the number of breaks
  # #                         labels = scales::number_format(accuracy = 0.01),  # Round labels to 2 decimal places
  # #                         name = "% Expression") + 
  # #   
  # #   scale_size_continuous(range = c(2, 8), name = "% Expression") +  # Size of points based on mean expression
  # #   
  # #   # Add asterisks (*) for significant genes (both positive and negative)
  # #   geom_text(data = subset(m_top_long, Significance != "Non-Significant"),
  # #             aes(label = "*"),
  # #             size = 4, color = "black", vjust = 0.5) +  # Center the asterisk on the dot
  # #   
  # #   labs(
  # #     title = title,
  # #     subtitle = subtitle_text,
  # #     x = "Group",
  # #     y = "Gene",
  # #     caption = caption
  # #   ) + 
  # #   theme_classic() +
  # #   theme(
  # #     plot.title = element_text(size = 10, face = "bold"),
  # #     plot.subtitle = element_text(size = 8),
  # #     axis.text = element_text(size = 6),
  # #     axis.title = element_text(size = 8),
  # #     legend.title = element_text(size = 8),
  # #     legend.text = element_text(size = 6),
  # #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  # #     axis.title.x = element_text(size = 8),  
  # #     axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
  # #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  # #   )
  # 
  # # Scale the Mean_Expression values (e.g., z-score scaling)
  # m_top_long$Scaled_Mean_Expression <- scale(m_top_long$Mean_Expression)
  # 
  # p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
  #   # Plot points for each gene's mean expression in both groups (just for background grid, no shapes)
  #   geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7, show.legend = TRUE) +  # Keep the legend for size and color
  #   
  #   # Apply a more granular diverging color scale with fewer breaks in the legend
  #   scale_color_gradient2(low = "white", mid = "red1", high = "darkred", 
  #                         midpoint = 0.7, 
  #                         limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
  #                         breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 5),  # Reduce the number of breaks
  #                         labels = scales::number_format(accuracy = 0.01),  # Round labels to 2 decimal places
  #                         name = "% Expression") + 
  #   
  #   scale_size_continuous(range = c(2, 8), name = "% Expression", guide = guide_legend(title = "% Expression")) +  # Size of points based on mean expression, ensure size legend
  #   
  #   # Add asterisks (*) for significant genes (both positive and negative)
  #   geom_text(data = subset(m_top_long, Significance != "Non-Significant"),
  #             aes(label = "*"),
  #             size = 4, color = "black", vjust = 0.5) +  # Center the asterisk on the dot
  #   
  #   labs(
  #     title = title,
  #     subtitle = subtitle_text,
  #     x = "Group",
  #     y = "Gene",
  #     caption = caption
  #   ) + 
  #   theme_classic() +
  #   theme(
  #     plot.title = element_text(size = 10, face = "bold"),
  #     plot.subtitle = element_text(size = 8),
  #     axis.text = element_text(size = 6),
  #     axis.title = element_text(size = 8),
  #     legend.title = element_text(size = 8),
  #     legend.text = element_text(size = 6),
  #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  #     axis.title.x = element_text(size = 8),  
  #     axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
  #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  #   )
  # 
  
  # # Plot the mean expression for each gene in the two groups (side by side)
  # p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
  #   # Plot points for each gene's mean expression in both groups
  #   geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7) +
  #   # geom_point(aes(color = Significance, size = Mean_Expression), alpha = 0.7) + 
  #   # scale_color_manual(values = c("Significant Positive" = "red", 
  #   #                               "Significant Negative" = "blue", 
  #   #                               "Non-Significant" = "gray")) +  # Set colors for positive, negative, and non-significant
  #   # scale_color_gradient(low = "white", high = "darkred", name = "% Expression") +
  #   # scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "% Expression") +
  #   scale_color_gradient2(low = "white", mid = "orangered", high = "darkred", 
  #                         midpoint = 0.5, limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
  #                         breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 10),
  #                         name = "% Expression") + 
  #   scale_size_continuous(range = c(2, 8), name = "% Expression") +
  #     # Draw boxes around Group 1 and Group 2
  #   # geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = NA, color = "black", size = 1.5) + # Box around Group 1
  #   # geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = NA, color = "black", size = 1.5) + # Box around Group 2
  #   labs(
  #     title = title,
  #     subtitle = subtitle_text,
  #     x = "Group",
  #     y = "Gene",
  #     caption = caption
  #   ) + 
  #   theme_classic() +
  #   theme(
  #     plot.title = element_text(size = 10, face = "bold"),
  #     plot.subtitle = element_text(size = 8),
  #     axis.text = element_text(size = 6),
  #     axis.title = element_text(size = 8),
  #     legend.title = element_text(size = 8),
  #     legend.text = element_text(size = 6),
  #     axis.text.y = element_text(size = 6),  # To make the gene names readable
  #     axis.title.x = element_text(size=8),  
  #     axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
  #     axis.ticks.x = element_blank()  # Remove x-axis ticks
  #   )
  
  
  # p <- EnhancedVolcano(m_top,
  #                      lab = labels,
  #                      x = 'avg_log2FC',
  #                      y = 'p_val_adj',
  #                      title = title,
  #                      subtitle = subtitle_text,
  #                      pCutoff = 0.05,
  #                      FCcutoff = 0.5,
  #                      labFace = 'bold',
  #                      pointSize = 4,
  #                      labSize = 5,
  #                      drawConnectors = TRUE,
  #                      widthConnectors = 1.0,
  #                      colConnectors = 'black',
  #                      legendPosition=NULL,
  #                      boxedLabels = TRUE,
  #                      max.overlaps=60)
  # 
  
  if(!is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(pathway,"_",additional_group,"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(pathway,"_",additional_group,"_Bulk_DEGs_for_",condition,".pdf")
    }
  }
  if(is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(pathway,"_","DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(pathway,"_","Pseudobulk_DEGs_for_",condition,".pdf")
    }
  }
  pdf(fs::path(dir.results,filename),width=10,height=5)
  # plot(p_dot2)
  plot(dp)
  dev.off()
  return(dp)
  # #GSEA
  # if (enrichment=="Yes") {
  #   #Gene set enrichment analysis
  #   gc()
  #   sce_sn_hep <- as.SingleCellExperiment(so)
  #   ## C2 category is according to canonical pathways: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
  #   geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  #   ### filter background to only include genes that we assessed
  #   geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
  #   geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_sn_hep),]
  #   m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
  #   stats <- m$p_val_adj
  #   names(stats) <- rownames(m)
  #   eaRes <- fgsea(pathways = m_list, stats = na.omit(stats))
  #   #ooEA <- order(eaRes$pval, decreasing = FALSE)
  #   #kable(head(eaRes[ooEA, 1:7], n = 20))
  #   # Convert the leadingEdge column to comma-separated strings
  #   eaRes$leadingEdge <- sapply(eaRes$leadingEdge, function(x) paste(x, collapse = ", "))
  #   gc()
  #   #Significant pathways
  #   sig <- eaRes %>% 
  #     filter(padj<0.05)
  #   
  #   #Plot top pathways
  #   # Subset top pathways for visualization
  #   top_pathways <- eaRes[order(-abs(eaRes$NES)), ][1:top_gsea, ]  # Top 10 pathways based on adjusted p-value
  #   
  #   # Define a significance threshold
  #   significance_threshold <- 0.05
  #   
  #   # Add a significance column for coloring
  #   top_pathways$significance <- ifelse(
  #     top_pathways$padj > significance_threshold, "Non-significant",
  #     ifelse(top_pathways$NES > 0, "Positive Significant", "Negative Significant")
  #   )
  #   
  #   if (!is.null(cell)) {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition2)
  #   } else {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) for ",condition2," (Bulk)")
  #   }
  #   # Create a bar plot
  #   gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
  #     geom_bar(stat = "identity") +
  #     coord_flip() +  # Flip coordinates for better readability
  #     scale_fill_manual(
  #       values = c(
  #         "Positive Significant" = "red",
  #         "Negative Significant" = "blue",
  #         "Non-significant" = "gray"
  #       ),
  #       name = "Significance"
  #     ) +
  #     labs(
  #       title = title,
  #       x = "Pathway",
  #       y = "Normalized Enrichment Score (NES)"
  #     ) +
  #     theme_minimal(base_size = 12) +
  #     theme(
  #       axis.text.y = element_text(size = 10),
  #       axis.title.y = element_text(size = 12),
  #       plot.title = element_text(size = 14, face = "bold")
  #     )
  #   
  #   
  #   if (!is.null(cell)){
  #     filename <- paste0(additional_group,"_GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
  #   } else {
  #     filename <- paste0(additional_group,"_Bulk_GSEA_for_",condition,".pdf") 
  #   }
  #   pdf(fs::path(dir.results,filename),width=20,height=20)
  #   plot(gsea_plot)
  #   dev.off()
  # }
  # #Make sure gene nemaes are in printed file
  # deg_results <- m 
  # deg_results$Gene <- rownames(deg_results)
  # 
  # #save results to excel file
  # write_multiple_sheets <- function(output_file, sheet_data) {
  #   # Create a new workbook
  #   wb <- createWorkbook()
  #   
  #   # Loop over the list of data frames
  #   for (sheet_name in names(sheet_data)) {
  #     # Add a new sheet to the workbook
  #     addWorksheet(wb, sheet_name)
  #     
  #     # Write the data frame to the sheet
  #     writeData(wb, sheet_name, sheet_data[[sheet_name]])
  #   }
  #   # Save the workbook to the specified file
  #   saveWorkbook(wb, file = output_file, overwrite = TRUE)
  # } 
  
  # if (enrichment=="Yes") {
  #   # Example usage
  #   df1 <- deg_results
  #   df2 <- eaRes
  # } else {
  #   df1 <- deg_results
  #   df2 <- NULL
  # }
  
  # # Specify the file name and data
  # if (!is.null(cell)){
  #   output_file <- fs::path(dir.results,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))  
  # } else {
  #   output_file <- fs::path(dir.results,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
  # }
  # write.csv(deg_results,output_file)
  # 
  # sheet_data <- list(
  #   "DEG" = df1,
  #   "Pathway_Results" = df2
  # )
  # 
  # # Call the function
  # write_multiple_sheets(output_file, sheet_data)
  # 
  #Save for IPA
  # Call the function
  # ipa_results <- deg_results %>% 
  #   dplyr::select(all_of(c("Gene","avg_log2FC","p_val","p_val_adj"))) %>% 
  #   dplyr::rename(ID=Gene,
  #                 `Fold Change`=avg_log2FC,
  #                 `P-Value`=p_val,
  #                 `Adjusted P-Value`=p_val_adj)
  # if (!is.null(cell)){
  #   output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))
  #   # output_file2 <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
  #   
  # } else {
  #   output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
  #   # output_file2 <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".csv"))
  #   
  # }
  # write.csv(ipa_results, output_file2)
  # return(p)
  # return(p_dot2)
}


# degs_fxn4 <- function(so,cell,exposure,gene_set,additional_group,exp_group,ref_group,enrichment,top_gsea,cellname,pathway) {
#   DefaultAssay(so) <- "RNA"
#   
#   #Conidtion names
#   if(is.null(additional_group)) {
#     condition <- paste0(str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
#     condition <- str_replace_all(condition," ","_")
#   } 
#   if(!is.null(additional_group)) {
#     condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
#     condition <- str_replace_all(condition," ","_")
#   }
#   #Differential Expression by Group
#   if (!is.null(cell)) {
#     cell_name <- str_replace_all(cell,"/","_")
#   }
#   
#   de.markers(so, gene_set, exposure, id2 = ref_group, id1 = exp_group, cell, "")
#   colnames(m)[2] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group,")")
#   colnames(m)[3] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",ref_group,")")
#   
#   # Filter for significant genes
#   m_top <- m
#   
#   
#   
#   if(is.null(additional_group)) {
#     condition2 <- paste0(str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
#     if (!is.null(cell)){
#       title <- paste0("DEGs in ",cell_name," cells \n ",condition2)
#     } else {
#       title <- paste0(cellname,"Cells, \n Individual PseudoBulk DEGs for ",condition2)
#     }
#   }
#   if(!is.null(additional_group)) {
#     condition2 <- paste0(additional_group," ",str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
#     if (!is.null(cell)){
#       title <- paste0(pathway," DEGs in ",cell_name," cells \n ",condition2)
#     } else {
#       title <- paste0(pathway," ",cellname,"Cells, \n Individual Pseudobulk DEGs for ",condition2)
#     }
#   } 
#   
#   
#   labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
#   total_individuals <- ncol(so)  # Total number of cells
#   total_genes <- nrow(m)  # Total number of genes in m_top
# 
#   # Assuming 'm' is a matrix of gene expression data, where rows are genes and columns are cells.
#   # Calculate mean expression for each gene in each group (ref_group and exp_group)
#   ref2 <- colnames(m)[3]
#   exp2 <- colnames(m)[2]
#   mean_exp_ref <- rowMeans(m[ref2])  # Mean expression for ref_group
#   mean_exp_exp <- rowMeans(m[exp2])  # Mean expression for exp_group
#   
#   # Combine mean expressions into a data frame
#   m_top2 <- data.frame(
#     Gene = rownames(m),
#     Mean_Expression_Ref = mean_exp_ref,
#     Mean_Expression_Exp = mean_exp_exp,
#     Significance = ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "Significant Positive", 
#                           ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "Significant Negative", 
#                                  "Non-Significant"))
#   )
#   
#   # Create subtitle with additional information
#   subtitle_text <- paste0("Positive Log2 FC = Greater Expression \n in ", exp_group)
#   caption <- paste0("Significant FDR-P<0.05, N = ", total_individuals, " | Total Genes = ", total_genes)
#   
#   # Reshape the data to make it suitable for ggplot (long format for groups)
#   m_top_long <- data.frame(
#     Gene = rep(rownames(m_top2), 2),
#     Group = rep(c(exp_group, ref_group), each = nrow(m_top2)),
#     Mean_Expression = c(m_top2$Mean_Expression_Exp, m_top2$Mean_Expression_Ref),
#     Significance = rep(m_top2$Significance, 2)
#   )
#   
#   
#   # Scale the Mean_Expression values (e.g., z-score scaling)
#   m_top_long$Scaled_Mean_Expression <- scale(m_top_long$Mean_Expression)
#   
#   p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
#     # Plot points for each gene's mean expression in both groups (just for background grid, no shapes)
#     geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7, show.legend = TRUE) +  # Keep the legend for size and color
#     
#     # Apply a more granular diverging color scale with fewer breaks in the legend
#     scale_color_gradient2(low = "white", mid = "red1", high = "darkred", 
#                           midpoint = 0.7, 
#                           limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
#                           breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 5),  # Reduce the number of breaks
#                           labels = scales::number_format(accuracy = 0.01),  # Round labels to 2 decimal places
#                           name = "% Expression") + 
#     
#     scale_size_continuous(range = c(2, 8), name = "% Expression", guide = guide_legend(title = "% Expression")) +  # Size of points based on mean expression, ensure size legend
#     
#     # Add asterisks (*) for significant genes (both positive and negative)
#     geom_text(data = subset(m_top_long, Significance != "Non-Significant"),
#               aes(label = "*"),
#               size = 4, color = "black", vjust = 0.5) +  # Center the asterisk on the dot
#     
#     labs(
#       title = title,
#       subtitle = subtitle_text,
#       x = "Group",
#       y = "Gene",
#       caption = caption
#     ) + 
#     theme_classic() +
#     theme(
#       plot.title = element_text(size = 10, face = "bold"),
#       plot.subtitle = element_text(size = 8),
#       axis.text = element_text(size = 6),
#       axis.title = element_text(size = 8),
#       legend.title = element_text(size = 8),
#       legend.text = element_text(size = 6),
#       axis.text.y = element_text(size = 6),  # To make the gene names readable
#       axis.title.x = element_text(size = 8),  
#       axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
#       axis.ticks.x = element_blank()  # Remove x-axis ticks
#     )
#   
# 
#   if(!is.null(additional_group)) {
#     if (!is.null(cell)){
#       filename <- paste0(pathway,"_",additional_group,"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
#     } else {
#       filename <- paste0(pathway,"_",additional_group,"_Bulk_DEGs_for_",condition,".pdf") 
#     }
#   }
#   if(is.null(additional_group)) {
#     if (!is.null(cell)){
#       filename <- paste0(pathway,"_","DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
#     } else {
#       filename <- paste0(pathway,"_","Pseudobulk_DEGs_for_",condition,".pdf") 
#     }
#   }
#   pdf(fs::path(dir.results,filename),width=5,height=10)
#   plot(p_dot2)
#   dev.off()
#   
#   #Make sure gene nemaes are in printed file
#   deg_results <- m 
#   deg_results$Gene <- rownames(deg_results)
#   
#   
#   
#   # Specify the file name and data
#   if (!is.null(cell)){
#     output_file <- fs::path(dir.results,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))  
#   } else {
#     output_file <- fs::path(dir.results,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
#   }
#   write.csv(deg_results,output_file)
#   
#  
#   #Save for IPA
#   # Call the function
#   ipa_results <- deg_results %>% 
#     dplyr::select(all_of(c("Gene","avg_log2FC","p_val","p_val_adj"))) %>% 
#     dplyr::rename(ID=Gene,
#                   `Fold Change`=avg_log2FC,
#                   `P-Value`=p_val,
#                   `Adjusted P-Value`=p_val_adj)
#   if (!is.null(cell)){
#     output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))
#     # output_file2 <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
#     
#   } else {
#     output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
#     # output_file2 <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".csv"))
#     
#   }
#   write.csv(ipa_results, output_file2)
#   # return(p)
#   return(p_dot2)
# }

# degs_fxn5 <- function(so,cell,exposure,gene_set,additional_group,exp_group,ref_group,enrichment,top_gsea,cellname,pathway) {
#   DefaultAssay(so) <- "RNA"
#   
#   #Conidtion names
#   if(is.null(additional_group)) {
#     condition <- paste0(str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
#     condition <- str_replace_all(condition," ","_")
#   } 
#   if(!is.null(additional_group)) {
#     condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
#     condition <- str_replace_all(condition," ","_")
#   }
#   #Differential Expression by Group
#   if (!is.null(cell)) {
#     cell_name <- str_replace_all(cell,"/","_")
#   }
#   
#   de.markers(so, gene_set, exposure, id2 = ref_group, id1 = exp_group, cell, "")
#   colnames(m)[2] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group,")")
#   colnames(m)[3] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",ref_group,")")
#   
#   # Filter for significant genes
#   m_top <- m
#   
#   
#   
#   if(is.null(additional_group)) {
#     condition2 <- paste0(str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
#     if (!is.null(cell)){
#       title <- paste0("DEGs in ",cell_name," cells \n ",condition2)
#     } else {
#       title <- paste0(cellname,"Cells, \n Individual PseudoBulk DEGs for ",condition2)
#     }
#   }
#   if(!is.null(additional_group)) {
#     condition2 <- paste0(additional_group," ",str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
#     if (!is.null(cell)){
#       title <- paste0(pathway," DEGs in ",cell_name," cells \n ",condition2)
#     } else {
#       title <- paste0(pathway," ",cellname,"Cells, \n Individual Pseudobulk DEGs for ",condition2)
#     }
#   } 
#   
#   
#   labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
#   total_individuals <- ncol(so)  # Total number of cells
#   total_genes <- nrow(m)  # Total number of genes in m_top
#   
#   # Assuming 'm' is a matrix of gene expression data, where rows are genes and columns are cells.
#   # Calculate mean expression for each gene in each group (ref_group and exp_group)
#   ref2 <- colnames(m)[3]
#   exp2 <- colnames(m)[2]
#   mean_exp_ref <- rowMeans(m[ref2])  # Mean expression for ref_group
#   mean_exp_exp <- rowMeans(m[exp2])  # Mean expression for exp_group
#   
#   # Combine mean expressions into a data frame
#   m_top2 <- data.frame(
#     Gene = rownames(m),
#     Mean_Expression_Ref = mean_exp_ref,
#     Mean_Expression_Exp = mean_exp_exp,
#     Significance = ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC > 0, "Significant Positive", 
#                           ifelse(m_top$p_val_adj < 0.05 & m_top$avg_log2FC < 0, "Significant Negative", 
#                                  "Non-Significant"))
#   )
#   
#   # Create subtitle with additional information
#   subtitle_text <- paste0("Positive Log2 FC = Greater Expression \n in ", exp_group)
#   caption <- paste0("Significant FDR-P<0.05, N = ", total_individuals, " | Total Genes = ", total_genes)
#   
#   # Reshape the data to make it suitable for ggplot (long format for groups)
#   m_top_long <- data.frame(
#     Gene = rep(rownames(m_top2), 2),
#     Group = rep(c(exp_group, ref_group), each = nrow(m_top2)),
#     Mean_Expression = c(m_top2$Mean_Expression_Exp, m_top2$Mean_Expression_Ref),
#     Significance = rep(m_top2$Significance, 2)
#   )
#   
#   
#   # Scale the Mean_Expression values (e.g., z-score scaling)
#   m_top_long$Scaled_Mean_Expression <- scale(m_top_long$Mean_Expression)
#   
#   p_dot2 <- ggplot(m_top_long, aes(x = Group, y = Gene)) + 
#     # Plot points for each gene's mean expression in both groups (just for background grid, no shapes)
#     geom_point(aes(color = Mean_Expression, size = Mean_Expression), alpha = 0.7, show.legend = TRUE) +  # Keep the legend for size and color
#     
#     # Apply a more granular diverging color scale with fewer breaks in the legend
#     scale_color_gradient2(low = "white", mid = "red1", high = "darkred", 
#                           midpoint = 0.7, 
#                           limits = c(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression)),
#                           breaks = seq(min(m_top_long$Mean_Expression), max(m_top_long$Mean_Expression), length.out = 5),  # Reduce the number of breaks
#                           labels = scales::number_format(accuracy = 0.01),  # Round labels to 2 decimal places
#                           name = "% Expression") + 
#     
#     scale_size_continuous(range = c(2, 8), name = "% Expression", guide = guide_legend(title = "% Expression")) +  # Size of points based on mean expression, ensure size legend
#     
#     # Add asterisks (*) for significant genes (both positive and negative)
#     geom_text(data = subset(m_top_long, Significance != "Non-Significant"),
#               aes(label = "*"),
#               size = 4, color = "black", vjust = 0.5) +  # Center the asterisk on the dot
#     
#     labs(
#       title = title,
#       subtitle = subtitle_text,
#       x = "Group",
#       y = "Gene",
#       caption = caption
#     ) + 
#     theme_classic() +
#     theme(
#       plot.title = element_text(size = 10, face = "bold"),
#       plot.subtitle = element_text(size = 8),
#       axis.text = element_text(size = 6),
#       axis.title = element_text(size = 8),
#       legend.title = element_text(size = 8),
#       legend.text = element_text(size = 6),
#       axis.text.y = element_text(size = 6),  # To make the gene names readable
#       axis.title.x = element_text(size = 8),  
#       axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text for better readability
#       axis.ticks.x = element_blank()  # Remove x-axis ticks
#     )
#   
#   
#   if(!is.null(additional_group)) {
#     if (!is.null(cell)){
#       filename <- paste0(pathway,"_",additional_group,"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
#     } else {
#       filename <- paste0(pathway,"_",additional_group,"_Bulk_DEGs_for_",condition,".pdf") 
#     }
#   }
#   if(is.null(additional_group)) {
#     if (!is.null(cell)){
#       filename <- paste0(pathway,"_","DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
#     } else {
#       filename <- paste0(pathway,"_","Pseudobulk_DEGs_for_",condition,".pdf") 
#     }
#   }
#   pdf(fs::path(dir.results,filename),width=5,height=10)
#   plot(p_dot2)
#   dev.off()
#   
#   #Make sure gene nemaes are in printed file
#   deg_results <- m 
#   deg_results$Gene <- rownames(deg_results)
#   
#   
#   
#   # Specify the file name and data
#   if (!is.null(cell)){
#     output_file <- fs::path(dir.results,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))  
#   } else {
#     output_file <- fs::path(dir.results,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
#   }
#   write.csv(deg_results,output_file)
#   
#   
#   #Save for IPA
#   # Call the function
#   ipa_results <- deg_results %>% 
#     dplyr::select(all_of(c("Gene","avg_log2FC","p_val","p_val_adj"))) %>% 
#     dplyr::rename(ID=Gene,
#                   `Fold Change`=avg_log2FC,
#                   `P-Value`=p_val,
#                   `Adjusted P-Value`=p_val_adj)
#   if (!is.null(cell)){
#     output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Results_",cell_name,"_cells_for_",condition,".csv"))
#     # output_file2 <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
#     
#   } else {
#     output_file2 <- fs::path(dir.ipa,paste0(pathway,"_","Pseudobulk_Results_for_",condition,".csv"))
#     # output_file2 <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".csv"))
#     
#   }
#   write.csv(ipa_results, output_file2)
#   # return(p)
#   return(p_dot2)
# }

degs_ind_fxn2 <- function(so,cell,exposure,gene_set,additional_group,exp_group,ref_group,enrichment,top_gsea,cellname) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  if(is.null(additional_group)) {
    condition <- paste0(str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  } 
  if(!is.null(additional_group)) {
    condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  }
  #Differential Expression by Group
  if (!is.null(cell)) {
    cell_name <- str_replace_all(cell,"/","_")
  }
  
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
  
  if(is.null(additional_group)) {
    condition2 <- paste0(str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0("DEGs in ",cell_name," cells for ",condition2)
    } else {
      title <- paste0(cellname,"Cells, Individual PseudoBulk DEGs for ",condition2)
    }
  }
  if(!is.null(additional_group)) {
    condition2 <- paste0(additional_group," ",str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0("DEGs in ",cell_name," cells for ",condition2)
    } else {
      title <- paste0(cellname,"Cells, Individual Pseudobulk DEGs for ",condition2)
    }
  } 
  
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  total_individuals <- ncol(so)  # Total number of cells
  total_genes <- nrow(m)  # Total number of genes in m_top
  # Create subtitle with additional information
  subtitle_text <- paste0(
    "Positive Log2 FC = Greater Expression in ", condition2, "\n",
    "(Significant at FDR-P<0.05, FC Threshold = 0.5, Min Gene Exp 5% Threshold) \n",
    "N = ", total_individuals, " | Total Genes = ", total_genes
  )
  # # Determine the max and min of your data for both x and y axes
  # x_max <- max(m_top$avg_log2FC, na.rm = TRUE)
  # x_min <- min(m_top$avg_log2FC, na.rm = TRUE)
  # 
  # y_max <- max(m_top$p_val_adj, na.rm = TRUE)
  # y_min <- min(m_top$p_val_adj, na.rm = TRUE)
  # 
  # # Add a little padding to the max and min values to extend the axes
  # x_padding <- 0.1 * (x_max - x_min)
  # y_padding <- 0.1 * (y_max - y_min)
  # 
  # # Create the plot with adjusted axis limits
  # p <- EnhancedVolcano(m_top,
  #                      lab = labels,
  #                      x = 'avg_log2FC',
  #                      y = 'p_val_adj',
  #                      title = title,
  #                      subtitle = subtitle_text,
  #                      pCutoff = 0.05,
  #                      FCcutoff = 0.5,
  #                      labFace = 'bold',
  #                      pointSize = 4,
  #                      labSize = 5,
  #                      drawConnectors = TRUE,
  #                      widthConnectors = 1.0,
  #                      colConnectors = 'black',
  #                      legendPosition=NULL,
  #                      boxedLabels = TRUE,
  #                      max.overlaps=60,
  #                      xlim = c(x_min - x_padding, x_max + x_padding),  # Extend x axis slightly beyond max
  #                      ylim = c(y_min - y_padding, y_max + y_padding))  # Extend y axis slightly beyond max
  # 
  
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       subtitle = subtitle_text,
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
  
  
  if(!is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(additional_group,"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(additional_group,"_Bulk_DEGs_for_",condition,".pdf") 
    }
  }
  if(is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0("DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0("Bulk_DEGs_for_",condition,".pdf") 
    }
  }
  pdf(fs::path(dir.results,filename),width=20,height=10)
  plot(p)
  dev.off()
  
  # #GSEA
  # if (enrichment=="Yes") {
  #   #Gene set enrichment analysis
  #   gc()
  #   sce_sn_hep <- as.SingleCellExperiment(so)
  #   ## C2 category is according to canonical pathways: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
  #   geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  #   ### filter background to only include genes that we assessed
  #   geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
  #   geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_sn_hep),]
  #   m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
  #   stats <- m$p_val_adj
  #   names(stats) <- rownames(m)
  #   eaRes <- fgsea(pathways = m_list, stats = na.omit(stats))
  #   #ooEA <- order(eaRes$pval, decreasing = FALSE)
  #   #kable(head(eaRes[ooEA, 1:7], n = 20))
  #   # Convert the leadingEdge column to comma-separated strings
  #   eaRes$leadingEdge <- sapply(eaRes$leadingEdge, function(x) paste(x, collapse = ", "))
  #   gc()
  #   #Significant pathways
  #   sig <- eaRes %>% 
  #     filter(padj<0.05)
  #   
  #   #Plot top pathways
  #   # Subset top pathways for visualization
  #   top_pathways <- eaRes[order(-abs(eaRes$NES)), ][1:top_gsea, ]  # Top 10 pathways based on adjusted p-value
  #   
  #   # Define a significance threshold
  #   significance_threshold <- 0.05
  #   
  #   # Add a significance column for coloring
  #   top_pathways$significance <- ifelse(
  #     top_pathways$padj > significance_threshold, "Non-significant",
  #     ifelse(top_pathways$NES > 0, "Positive Significant", "Negative Significant")
  #   )
  #   
  #   if (!is.null(cell)) {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition2)
  #   } else {
  #     title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) for ",condition2," (Bulk)")
  #   }
  #   # Create a bar plot
  #   gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
  #     geom_bar(stat = "identity") +
  #     coord_flip() +  # Flip coordinates for better readability
  #     scale_fill_manual(
  #       values = c(
  #         "Positive Significant" = "red",
  #         "Negative Significant" = "blue",
  #         "Non-significant" = "gray"
  #       ),
  #       name = "Significance"
  #     ) +
  #     labs(
  #       title = title,
  #       x = "Pathway",
  #       y = "Normalized Enrichment Score (NES)"
  #     ) +
  #     theme_minimal(base_size = 12) +
  #     theme(
  #       axis.text.y = element_text(size = 10),
  #       axis.title.y = element_text(size = 12),
  #       plot.title = element_text(size = 14, face = "bold")
  #     )
  #   
  #   
  #   if (!is.null(cell)){
  #     filename <- paste0(additional_group,"_GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
  #   } else {
  #     filename <- paste0(additional_group,"_Bulk_GSEA_for_",condition,".pdf") 
  #   }
  #   pdf(fs::path(dir.results,filename),width=20,height=20)
  #   plot(gsea_plot)
  #   dev.off()
  # }
  #Make sure gene nemaes are in printed file
  deg_results <- m 
  deg_results$Gene <- rownames(deg_results)
  
  # #save results to excel file
  # write_multiple_sheets <- function(output_file, sheet_data) {
  #   # Create a new workbook
  #   wb <- createWorkbook()
  #   
  #   # Loop over the list of data frames
  #   for (sheet_name in names(sheet_data)) {
  #     # Add a new sheet to the workbook
  #     addWorksheet(wb, sheet_name)
  #     
  #     # Write the data frame to the sheet
  #     writeData(wb, sheet_name, sheet_data[[sheet_name]])
  #   }
  #   # Save the workbook to the specified file
  #   saveWorkbook(wb, file = output_file, overwrite = TRUE)
  # } 
  
  # if (enrichment=="Yes") {
  #   # Example usage
  #   df1 <- deg_results
  #   df2 <- eaRes
  # } else {
  #   df1 <- deg_results
  #   df2 <- NULL
  # }
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
  } else {
    output_file <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".csv"))
  }
  write.csv(deg_results,output_file)
  
  # sheet_data <- list(
  #   "DEG" = df1,
  #   "Pathway_Results" = df2
  # )
  # 
  # # Call the function
  # write_multiple_sheets(output_file, sheet_data)
  # 
  #Save for IPA
  # Call the function
  ipa_results <- deg_results %>% 
    dplyr::select(all_of(c("Gene","avg_log2FC","p_val","p_val_adj"))) %>% 
    dplyr::rename(ID=Gene,
                  `Fold Change`=avg_log2FC,
                  `P-Value`=p_val,
                  `Adjusted P-Value`=p_val_adj)
  if (!is.null(cell)){
    #output_file2 <- fs::path(dir.ipa,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
    output_file2 <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".csv"))  
    
  } else {
    #output_file2 <- fs::path(dir.ipa,paste0("Bulk_Results_for_",condition,".csv"))
    output_file2 <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".csv"))
    
  }
  write.csv(ipa_results, output_file2)
  return(p)
}

degs_fxn_path <- function(so,cell,exposure,gene_set,additional_group,exp_group,ref_group,enrichment,top_gsea,pathway) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  if(is.null(additional_group)) {
    condition <- paste0(str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  } 
  if(!is.null(additional_group)) {
    condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," ")),"_",exp_group,"_vs_",ref_group)
    condition <- str_replace_all(condition," ","_")
  }
  #Differential Expression by Group
  if (!is.null(cell)) {
    cell_name <- str_replace_all(cell,"/","_")
  }
  
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
  
  if(is.null(additional_group)) {
    condition2 <- paste0(str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0("DEGs in ",cell_name," cells for ",condition2)
    } else {
      title <- paste0("Bulk DEGs for ",condition2)
    }
  }
  if(!is.null(additional_group)) {
    condition2 <- paste0(additional_group," ",str_to_title(exposure)," (",exp_group," vs. ",ref_group,")")
    if (!is.null(cell)){
      title <- paste0("DEGs in ",cell_name," cells for ",condition2)
    } else {
      title <- paste0("Bulk DEGs for ",condition2)
    }
  } 
  
  
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  total_cells <- ncol(so)  # Total number of cells
  total_genes <- nrow(m)  # Total number of genes in m_top
  # Create subtitle with additional information
  subtitle_text <- paste0(
    "Positive Log2 FC = Greater Expression in ", condition2, "\n",
    "(Significant at FDR-P<0.05, FC Threshold = 0.5, Min Gene Exp 5% Threshold) \n",
    "Total Cells = ", total_cells, " | Total Genes = ", total_genes
  )
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       subtitle = subtitle_text,
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
  
  
  if(!is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0(additional_group,"_DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0(additional_group,"_Bulk_DEGs_for_",condition,".pdf") 
    }
  }
  if(is.null(additional_group)) {
    if (!is.null(cell)){
      filename <- paste0("DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
    } else {
      filename <- paste0("Bulk_DEGs_for_",condition,".pdf") 
    }
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
      title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition2)
    } else {
      title <- paste0(additional_group,"_Top Enriched Pathways by NES (GSEA) for ",condition2," (Bulk)")
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
      filename <- paste0(additional_group,"_GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
    } else {
      filename <- paste0(additional_group,"_Bulk_GSEA_for_",condition,".pdf") 
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
  } 
  
  if (enrichment=="Yes") {
    # Example usage
    df1 <- deg_results
    df2 <- eaRes
  } else {
    df1 <- deg_results
    df2 <- NULL
  }
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".xlsx"))  
  } else {
    output_file <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".xlsx"))
  }
  
  sheet_data <- list(
    "DEG" = df1,
    "Pathway_Results" = df2
  )
  
  # Call the function
  write_multiple_sheets(output_file, sheet_data)
  return(p)
  
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

#Mixed Models ----
##a. Group ----
lmm_fxn <- function(exp,ref) {
  data_subgroup <- data %>% 
    # filter(group=="Type_2_Diabetes" | group=="Obese_Control")
    filter(group==ref | group==exp)
  data_subgroup$group <- factor(data_subgroup$group)
  data_subgroup$group <-relevel(data_subgroup$group, ref = ref)
  # full_results <- data.frame()
  # Prepare parallel execution of gene models
  full_results <- foreach(gene = gene_list[1:48], .combine = rbind, .packages = c("glmmTMB", "dplyr")) %dopar% {
    # for (gene in gene_list) {
    # for (gene in gene_list[1:20]) { #tester genes
    m0 <- as.formula(paste0(gene," ~ group + (1 | kit_id)"))
    model <- glmmTMB(m0,data=data_subgroup,family = nbinom2,  # Negative Binomial with a log link
                     ziformula = ~ 1)
    #If doesnt fit, try nb 
    if(summary(model)$coef$zi[1,4]=="NaN") {
      model2 <- glmmTMB(m0, data=data_subgroup,
                        family = nbinom2)
      foldchange <- round(exp(summary(model2)$coef$cond[2,1]) - 1,5)
      pval <- summary(model2)$coef$cond[2,4]
      model_name <- "Negative Binomial"
    } 
    if(summary(model)$coef$zi[1,4]!="NaN") {
      #Check if zero inflation component is non-significant, if non-sig, run nb
      if(summary(model)$coef$zi[1,4]>=0.05) {
        model2 <- glmmTMB(m0, data=data_subgroup,
                          family = nbinom2)
        foldchange <- round(exp(summary(model2)$coef$cond[2,1]) - 1,5)
        pval <- summary(model2)$coef$cond[2,4]
        model_name <- "Negative Binomial"
      } 
      #If sig, pull estimates from zinb 
      if(summary(model)$coef$zi[1,4]<0.05) {
        # foldchange <- round(exp(summary(model)$coef$cond[2,1]) - 1,5)
        foldchange <- round(exp(summary(model)$coef$cond[2,1]) - 1,5)
        pval <- summary(model)$coef$cond[2,4]
        model_name <- "Zero-Inflated Negative Binomial"
      }
    }
    # Collect results for each gene
    results <- data.frame(Gene = gene, FoldChange = foldchange, PValue = pval,ModelName=model_name)
    return(results)
    # results <- data.frame(Gene=gene,FoldChange=foldchange,PValue=pval)
    # full_results <- rbind(full_results,results)
    # # print(p)
    # # return(full_results)
  }
  #Make volcano plot of all gene results for group
  full_results <- full_results %>% 
    mutate(fdr=p.adjust(PValue,method="fdr"))
  
  m_top <- full_results
  # significant_genes <- m_top %>% filter(fdr < 0.05)
  # 
  # # Select the top 10 positive and 10 negative log2FC genes based on the largest magnitude of fold change
  # top_positive_by_fc <- significant_genes %>% 
  #   filter(FoldChange > 0) %>% 
  #   arrange(desc(abs(FoldChange))) %>%  # Sort by absolute fold change (largest first)
  #   head(10)  # Top 10 positive fold changes
  # 
  # top_negative_by_fc <- significant_genes %>% 
  #   filter(FoldChange < 0) %>% 
  #   arrange(desc(abs(FoldChange))) %>%  # Sort by absolute fold change (largest first)
  #   head(10)  # Top 10 negative fold changes
  # 
  # # Select the top 10 positive and 10 negative log2FC genes based on significance (fdr)
  # top_positive_by_significance <- significant_genes %>% 
  #   filter(FoldChange > 0) %>% 
  #   arrange(fdr) %>%  # Sort by smallest p-value (most significant)
  #   head(10)  # Top 10 most significant positive fold changes
  # 
  # top_negative_by_significance <- significant_genes %>% 
  #   filter(FoldChange < 0) %>% 
  #   arrange(fdr) %>%  # Sort by smallest p-value (most significant)
  #   head(10)  # Top 10 most significant negative fold changes
  # 
  # # Combine top fold-change based and significance-based genes into a final list
  # top_genes <- rbind(
  #   top_positive_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
  #   top_negative_by_fc %>% mutate(Selection = "Top 10 by Fold Change"),
  #   top_positive_by_significance %>% mutate(Selection = "Top 10 by Significance"),
  #   top_negative_by_significance %>% mutate(Selection = "Top 10 by Significance")
  # )
  # 
  # labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  # total_cells <-  length(unique(data_subgroup$barcode)) # Total number of cells
  # total_genes <-  length(unique(full_results$Gene)) # Total number of genes in m_top
  # # Create subtitle with additional information
  # 
  # title <- paste0("Differential Gene Expression for ",str_replace_all(exp,"_"," ")," vs. ",str_replace_all(ref,"_"," "))
  # subtitle_text <- paste0(
  #   "Fold Change = Greater Gene Expression in ",str_replace_all(exp,"_"," ")," compared to ",str_replace_all(ref,"_"," ")," \n(Significant at FDR-P<0.05, FC Threshold = 0.5, Min Gene Exp 5% Threshold) \n",
  #   "Total Cells = ", total_cells, " | Total Genes = ", total_genes
  # )
  # p <- EnhancedVolcano(m_top,
  #                      lab = labels,
  #                      x = 'FoldChange',
  #                      y = 'fdr',
  #                      title = title,
  #                      subtitle = subtitle_text,
  #                      pCutoff = 0.05,
  #                      FCcutoff = 0.5,
  #                      labFace = 'bold',
  #                      pointSize = 4,
  #                      labSize = 5,
  #                      drawConnectors = TRUE,
  #                      widthConnectors = 1.0,
  #                      colConnectors = 'black',
  #                      legendPosition=NULL,
  #                      boxedLabels = TRUE,
  #                      max.overlaps=60)
  
  # Assuming m_top contains FoldChange, PValue, and fdr columns
  # Transform PValue into -log10(p-value) for the y-axis
  m_top$logPValue <- -log10(m_top$PValue)
  
  # Calculate the log2 of the fold change (logFC)
  m_top$logFC <- log2(m_top$FoldChange)
  
  # Create a new column for significance based on the p-value cutoff and FC cutoff
  m_top$Significance <- ifelse(m_top$fdr<0.05, "Significant", "Not Significant")
  
  # Create a new column to combine both significance and direction (positive or negative)
  m_top$Color <- ifelse(m_top$Significance == "Significant" & m_top$logFC > 0, "Positive Significant", 
                        ifelse(m_top$Significance == "Significant" & m_top$logFC < 0, "Negative Significant", 
                               "Non-Significant"))
  
  # Select the genes that are significant to label
  m_top$Label <- ifelse(m_top$Significance == "Significant", as.character(m_top$Gene), NA)
  
  # Plot
  p <- ggplot(m_top, aes(x = logFC, y = logPValue, color = Color)) +
    geom_point(aes(shape = Color), size = 3, alpha = 0.7) +  # Set alpha for transparency
    scale_color_manual(values = c("Positive Significant" = "red", "Negative Significant" = "blue", "Non-Significant" = "gray")) +  # Custom color scale
    scale_shape_manual(values = c("Positive Significant" = 19, "Negative Significant" = 19, "Non-Significant" = 1)) +  # Same shape for all but Non-Sig
    geom_text(aes(label = Label), size = 3, color = "black", vjust = 1, hjust = 1) +  # Add labels for significant genes in black
    labs(x = "log2 Fold Change", y = "-log10(p-value)", title = "Volcano Plot") +
    theme_minimal() +
    theme(legend.position = "top") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
  
  # Print the plot
  print(p)
  
  return(p)
}



#so=so_sub
#cell=NULL
#exposure="Diabetes"
#gene_set=rownames(so_sub)
#batch_size=1000
#exp_group="Type 2 Diabetes"
#ref_group="Lean Control"

#Mast Function----
mast_fxn <- function(so,cell,exposure,covariate,gene_set,batch_size,exp_group,ref_group,additional_group) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  if (!is.null(exp_group)) {
    condition <- paste0(additional_group,"_",str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group," vs. ",ref_group,") ")
  } else {
    condition <- paste0(additional_group,"_",str_to_title(exposure))
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
    so_sub <- subset(so, features = subset)
    
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
      # filter(fdr<0.05) %>% 
      dplyr::rename(LogFoldChange=coef) %>% 
      dplyr::rename(Gene=primerid)
    
    # Append the batch results to the combined results using rbindlist
    #combined_results <- rbind(combined_results, result1)
    combined_results <- rbindlist(list(combined_results, result1), use.names = TRUE, fill = TRUE)
    
  }
  
  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- paste0("Results_MAST_for_",str_replace_all(" ", "_",additional_group),"_",cell_name,"_cells_for_",condition,".xlsx")
  } else {
    output_file <- paste0("Results_MAST_Bulk_",str_replace_all(" ", "_",additional_group),"_",condition,".xlsx")
  }
  
  write.csv(combined_results,fs::path(dir.results,output_file))
  
  # # Filter top positive and negative log2FoldChange, including significance based on fdr
  # top_pos <- as.data.frame(combined_results) %>%
  #   filter(component == "FC_H") %>%
  #   filter(LogFoldChange > 0 & fdr < 0.05)  # Positive fold change and significant (fdr < 0.05)
  # 
  # top_neg <- as.data.frame(combined_results) %>%
  #   filter(component == "FC_H") %>%
  #   filter(LogFoldChange < 0 & fdr < 0.05)  # Negative fold change and significant (fdr < 0.05)
  # 
  # # Non-significant data
  # top_non_sig <- as.data.frame(combined_results) %>%
  #   filter(component == "FC_H") %>%
  #   filter(fdr >= 0.05)  # Non-significant data
  # 
  # # Combine positive, negative, and non-significant results
  # top_genes <- bind_rows(top_pos, top_neg, top_non_sig)
  # 
  # # Add Direction column for positive and negative fold changes
  # top_genes$Direction <- ifelse(top_genes$LogFoldChange > 0, "Positive", 
  #                               ifelse(top_genes$LogFoldChange < 0, "Negative", "Non-significant"))
  # 
  # # Set sorting order: negative significant genes, positive significant genes, and non-significant
  # top_genes <- top_genes %>%
  #   mutate(Gene = factor(Gene, levels = Gene[order(Direction, LogFoldChange)]))  # Order by Direction and log2FC
  # 
  # # Assign colors based on significance and fold change
  # top_genes$Color <- case_when(
  #   top_genes$fdr < 0.05 & top_genes$LogFoldChange > 0 ~ "red",   # Positive and significant
  #   top_genes$fdr < 0.05 & top_genes$LogFoldChange < 0 ~ "blue",  # Negative and significant
  #   top_genes$fdr >= 0.05 ~ "grey"                                # Non-significant
  # )
  # 
  # # Title setup
  # if (!is.null(cell)) {
  #   title <- paste0(additional_group, " ", condition, " ", cell_name)
  # } else {
  #   title <- paste0("Bulk", additional_group, " ", condition)
  # }
  # 
  # # Volcano plot
  # p <- ggplot(top_genes, aes(x = LogFoldChange, y = -log10(fdr), color = Color)) +
  #   geom_point(size = 3) +  # Plot points
  #   labs(x = "logFC", y = "-log10(FDR)", title = title) +
  #   scale_color_identity() +  # Use predefined colors
  #   theme_minimal() +
  #   theme(text = element_text(family = "Times")) +
  #   theme(legend.position = "none") +
  #   theme(plot.title = element_text(hjust = 0.5))  # Center the plot title
  
  # # Filter top 10 positive and negative log2FoldChange
  # top_pos <- as.data.frame(combined_results) %>%
  #   filter(component == "FC_H") %>%
  #   # filter(fdr < 0.05) %>%
  #   filter(LogFoldChange > 0) 
  #   # arrange(desc(LogFoldChange)) %>%
  #   # dplyr::slice(1:30) 
  # 
  # top_neg <- as.data.frame(combined_results) %>%
  #   filter(component == "FC_H") %>%
  #   # filter(fdr < 0.05) %>%
  #   filter(LogFoldChange < 0) 
  #   # arrange(LogFoldChange) %>%
  #   # dplyr::slice(1:30) 
  # 
  # if (dim(top_pos)[1] > 0 | dim(top_neg)[1] > 0) {
  #   top_pos$Direction <- "Positive"
  #   top_neg$Direction <- "Negative"
  #   
  #   # Combine and prepare for plotting
  #   top_genes <- bind_rows(top_pos, top_neg)
  #   
  #   # Set sorting order: all negative genes first, then all positive genes
  #   top_genes <- top_genes %>%
  #     mutate(Gene = factor(Gene, levels = Gene[order(Direction, LogFoldChange)]))  # Order by Direction and log2FC
  #   
  #   if (!is.null(cell)) {
  #     title <- paste0(additional_group, " ", condition, " ", cell_name)
  #   } else {
  #     title <- paste0("Bulk", additional_group, " ", condition)
  #   }
  #   
  #   # Volcano plot
  #   p <- ggplot(top_genes, aes(x = LogFoldChange, y = -log10(fdr), color = Direction)) +
  #     geom_point(size = 3) +  # Plot points
  #     labs(x = "logFC", y = "-log10(FDR)", title = title) +
  #     scale_color_manual(values = c("blue", "red"), labels = c("Lower Expression", "Higher Expression")) +
  #     theme_minimal() +
  #     theme(text = element_text(family = "Times")) +
  #     theme(legend.position = "none") +
  #     theme(plot.title = element_text(hjust = 0.5))  # Center the plot title
  #   
  
  # Filter top positive and negative log2FoldChange, including significance based on fdr
  top_pos <- as.data.frame(combined_results) %>%
    filter(component == "FC_H") %>%
    filter(LogFoldChange > 0 & fdr < 0.05)  # Positive fold change and significant (fdr < 0.05)
  
  top_neg <- as.data.frame(combined_results) %>%
    filter(component == "FC_H") %>%
    filter(LogFoldChange < 0 & fdr < 0.05)  # Negative fold change and significant (fdr < 0.05)
  
  # Non-significant data
  top_non_sig <- as.data.frame(combined_results) %>%
    filter(component == "FC_H") %>%
    filter(fdr >= 0.05)  # Non-significant data
  
  # Combine positive, negative, and non-significant results
  top_genes <- bind_rows(top_pos, top_neg, top_non_sig)
  
  # Add Direction column for positive and negative fold changes
  top_genes$Direction <- ifelse(top_genes$LogFoldChange > 0, "Positive", 
                                ifelse(top_genes$LogFoldChange < 0, "Negative", "Non-significant"))
  
  # Set sorting order: negative significant genes, positive significant genes, and non-significant
  top_genes <- top_genes %>%
    mutate(Gene = factor(Gene, levels = Gene[order(Direction, LogFoldChange)]))  # Order by Direction and log2FC
  
  # Assign colors based on significance and fold change
  top_genes$Color <- case_when(
    top_genes$fdr < 0.05 & top_genes$LogFoldChange > 0 ~ "red",   # Positive and significant
    top_genes$fdr < 0.05 & top_genes$LogFoldChange < 0 ~ "blue",  # Negative and significant
    top_genes$fdr >= 0.05 ~ "grey"                                # Non-significant
  )
  
  # Identify top 10 most positive and top 10 most negative significant genes
  top_pos_10 <- top_pos %>%
    arrange(desc(LogFoldChange)) %>%
    slice_head(n = 10)  # Top 10 positive significant genes
  
  top_neg_10 <- top_neg %>%
    arrange(LogFoldChange) %>%
    slice_head(n = 10)  # Top 10 negative significant genes
  
  # Combine the top 10 positive and negative significant genes for labeling
  top_labels <- bind_rows(top_pos_10, top_neg_10)
  
  # Count the number of genes and cells included in the plot
  num_genes <- nrow(top_genes)
  num_cells <- length(unique(top_genes$Gene))  # Assuming each Gene represents a unique cell (adjust if you have specific cell data)
  
  # Title setup
  if (!is.null(cell)) {
    title <- paste0(additional_group, " ", condition, " ", cell_name)
  } else {
    title <- paste0("Bulk", additional_group, " ", condition)
  }
  subtitle <- paste0("Genes: ", num_genes, " | Cells: ", num_cells)
  
  # Volcano plot
  g <- ggplot(top_genes, aes(x = LogFoldChange, y = -log10(fdr), color = Color)) +
    geom_point(size = 3) +  # Plot points
    labs(x = "logFC", y = "-log10(FDR)", title = title,subtitle = subtitle) +
    scale_color_identity() +  # Use predefined colors
    theme_minimal() +
    theme(text = element_text(family = "Times")) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))  # Center the plot title
  
  # Add labels for the top 10 positive and negative significant genes
  p <- g + geom_text(data = top_labels, aes(x = LogFoldChange, y = -log10(fdr), label = Gene), 
                     vjust = 1, hjust = 1, size = 3, color = "black")
  
  # Generate the file name for the figure
  if (!is.null(cell)) {
    figure_file <- paste0("MAST_Volcano Plot_", cell_name, "_cells_for_", condition, "_", str_replace_all(" ", "_", additional_group), ".pdf")
  } else {
    figure_file <- paste0("MAST_Volcano_for_", condition, "_", str_replace_all(" ", "_", additional_group), ".pdf")
  }
  
  # Save the plot as a PDF
  pdf(fs::path(dir.results, figure_file), width = 20, height = 20)
  plot(p)
  dev.off()
  
  
  if (!is.null(cell)){
    figure_file<- paste0("MAST_Volcano Plot_",cell_name,"_cells_for_",condition,"_",str_replace_all(" ", "_",additional_group),".pdf")
  } else {
    figure_file <- paste0("MAST_Volcano_for_",condition,"_",str_replace_all(" ", "_",additional_group),".pdf")
  }
  
  pdf(fs::path(dir.results,figure_file),width=20,height=20)
  plot(p)
  dev.off()
  # }
  
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
    
    if (!is.null(cell)){
      title <- paste0(additional_group," ",condition," ",cell_name)
    } else {
      title <- paste0("Bulk",additional_group," ",condition)
    }
    
    # # Bar chart with flipped x and y axes
    p2 <- ggplot(top_genes, aes(x = LogFoldChange, y = Gene, fill = Direction)) +
      geom_bar(stat = "identity") +
      labs(x = "logFC", y = "Gene",title=title) +
      scale_fill_manual(values = c("blue", "red"), labels = c("Lower Expression", "Higher Expression")) +
      theme_minimal()+
      theme(element_text(family="Times"))+
      theme(legend.position="none")
    
    if (!is.null(cell)){
      figure_file<- paste0("MAST_Barchart_",cell_name,"_cells_for_",condition,"_",str_replace_all(" ", "_",additional_group),".pdf")
    } else {
      figure_file <- paste0("MAST_Barchart_for_",condition,"_",str_replace_all(" ", "_",additional_group),".pdf")
    }
    
    pdf(fs::path(dir.results,figure_file),width=20,height=20)
    plot(p2)
    dev.off()
  }
  return (p)
  # return(p2)
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

dp.formatted <- function(seurat_object, genes, celltype, group.by, m,
                         colorlow = "#83c5be", colormid = "#f4f1bb", colorhigh = "#d90429",exp_group,ref_group)
{
  pt.combined <- DotPlot(seurat_object,
                         features = genes,idents = celltype, group.by = group.by,
                         scale = F, cols = "RdYlBu"
  )$data 
  
  
  # m <- m %>%
  #   mutate(significant = ifelse(p_val_adj < 0.05, "*", NA))  # Add '*' for significant genes
  # 
  pt.plot <- pt.combined %>% 
    ggplot(aes(x=features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) + 
    geom_point() +
    # geom_text(aes(label = significant), size = 5, color = "black", vjust = -0.5, hjust = 0.5) +  # Add '*' for significant genes
    theme_bw() +
    scale_color_gradient2(low = colorlow, mid = colormid, high = colorhigh, midpoint = 2,
                          guide = guide_colorbar(label.vjust = 0.8, ticks = F, draw.ulim = T, draw.llim = T),
                          limits = c(0,6)) +
    scale_size(range = c(0,6), 
               limits = c(0,100)) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, vjust = 0.5),
          legend.spacing.x = unit(.1, "cm"),
          legend.direction = "horizontal") +
    guides(size = guide_legend(label.position = "bottom",
                               title.position = "top"),
           color = guide_colorbar(label.position = "bottom",
                                  title.position = "top")) +
    labs(color = "Scaled average expression",
         size = "Expression (%) ") + 
    scale_y_discrete(limits=rev)
  
  pt.table <- m %>%
    filter(rownames(m) %in% genes) %>%
    dplyr::mutate(p_val_rounded = round(p_val, 4),
                  p_val = p_format(p_val_rounded, trailing.zero = T, accuracy = 0.001, digits = 3),
                  p_val_adj_rounded = round(p_val_adj, 4),
                  p_val_adj = p_format(p_val_adj_rounded, trailing.zero = T, accuracy = 0.001, digits = 3),
                  pct.1 = sprintf("%.3f", pct.1),
                  pct.2 = sprintf("%.3f", pct.2),
                  avg_log2FC = sprintf("%.3f", avg_log2FC)) %>% 
    dplyr::select(pct.1, pct.2, avg_log2FC, p_val, p_val_adj) 
  gg.pt.table <- ggtexttable(pt.table,
                             cols = c(exp_group, ref_group, "Log2FC", "p-value", "q-value"),
                             theme = ttheme("blank")) %>%
    tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1) %>%
    tab_add_title(paste0("% Expressed in ",cellname," Cells"))
  
  pt.plot_table <- ggarrange(pt.plot, NULL, gg.pt.table,
                             nrow = 1, widths = c(1,-0.1,1), common.legend = F,
                             legend = "top")
  
}


#GSEA Functions ----
##a. Set up functions ----
# function for GSEA
# Set relevant paths
list.files()
bg_path <- c(fs::path(dir.dat,"/GSEA/"))
list.files(bg_path)

# Functions
## Function: Adjacency matrix to list 
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt 
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# List of acronyms to uppercase
acronyms <- c("rna", "dna", "mtor", "foxo", "ppar", "nmd", "fgfr", "robo", 
              "bhl", "cov", "jak", "stat", "wnt", "hiv", "bcl", "mapk",
              "pt", "tal", "pc", "ic", "ec", "fibvsmcp")
# unique(full_results$Gene))
special_mixed <- c("rrna", "mrna", "trna", "gtpase", "atpase", "robos", "slits", "fibvsmcp")
special_replacements <- c("rRNA", "mRNA", "tRNA", "GTPase", "ATPase", "ROBOs", "SLITs", "FIB/VSMC/P")

replace_mixed_case <- function(text, from, to) {
  for (i in seq_along(from)) {
    pattern <- paste0("\\b", from[i], "\\b")
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), to[i])
  }
  return(text)
}

capitalize_acronyms <- function(text, terms) {
  for (term in terms) {
    pattern <- paste0("\\b", term, "\\b")
    replacement <- toupper(term)
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), replacement)
  }
  return(text)
}

## b. Fgsea function ----
plot_fgsea <- function(fgsea_res,
                       top_n = 30,
                       title = "Top Enriched Pathways",
                       xlimit = 3,
                       xnudge = xlimit/100,
                       text1 = 6.5,
                       text2 = 18,
                       text3 = 20,
                       face = "plain") {
  
  fgsea_res <- fgsea_res %>%
    arrange(pval) %>%
    head(top_n) %>%
    mutate(
      direction = case_when(NES < 0 ~ "Negative", NES > 0 ~ "Positive"),
      pathway_clean = str_remove(pathway, "^KEGG_"), 
      pathway_clean = str_remove(pathway_clean, "^REACTOME_"), 
      pathway_clean = str_remove(pathway_clean, "^GOBP_"), 
      pathway_clean = str_remove(pathway_clean, "^GOMF_"), 
      pathway_clean = str_replace_all(pathway_clean, "_", " "),
      pathway_clean = str_to_sentence(pathway_clean),
      pathway_clean = str_replace_all(pathway_clean, "\\bi\\b", "I"),
      pathway_clean = str_replace_all(pathway_clean, "\\bii\\b", "II"),
      pathway_clean = str_replace_all(pathway_clean, "\\biii\\b", "III"),
      pathway_clean = str_replace_all(pathway_clean, "\\biv\\b", "IV"),
      pathway_clean = str_replace_all(pathway_clean, "\\bv\\b", "V"),
      pathway_clean = str_replace_all(pathway_clean, regex("\\(immune\\)", ignore_case = TRUE), "(IMMUNE)"),
      pathway_clean = capitalize_acronyms(pathway_clean, acronyms),
      pathway_clean = replace_mixed_case(pathway_clean, special_mixed, special_replacements),
      pathway_clean = paste0(pathway_clean, " (", size, ")")
    ) %>%
    arrange(pval)
  
  fgsea_res$pathway_clean <- reorder(fgsea_res$pathway_clean, fgsea_res$pval)
  
  fgsea_res %>%
    ggplot(aes(x = -log10(pval), y = fct_rev(pathway_clean), label = pathway_clean)) +
    geom_point(aes(size = abs(NES), color = direction, alpha = 0.8)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    geom_text(aes(group = pathway_clean, color = direction, fontface = face), 
              hjust = 0, size = text1, nudge_x = xnudge) +
    scale_size_binned() +
    scale_color_manual(values = c("Positive" = "#c75146", "Negative" = "#2c7da0")) +
    scale_x_continuous(limits = c(0, xlimit), expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(add = 1)) +
    labs(
      x = "-log(p-value)",
      y = "Pathways",
      color = "Direction",
      size = "NES",
      title = title
    ) +
    guides(alpha = "none") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = text3),
      axis.title = element_text(size = text3),
      axis.ticks.y = element_blank(), 
      legend.position = c(0.9, 0.2),
      legend.background = element_blank(),
      legend.box.background = element_rect(color = "black"),
      legend.title = element_text(size = text2),
      legend.text = element_text(size = text2),
      title = element_text(size = text3)
    )
}
##c. Transposed fgsea function ----
plot_fgsea_transpose <- function(fgsea_res,
                                 top_n = 30,
                                 title = "Top Enriched Pathways",
                                 xlimit = 3,
                                 xnudge = xlimit/100,
                                 text1 = 6.5,
                                 text2 = 18,
                                 text3 = 20) {
  
  fgsea_res <- fgsea_res %>%
    arrange(pval) %>%
    head(top_n) %>%
    mutate(
      direction = case_when((NES < 0 & pval <= 0.05 ~ "Negative (p < 0.05)"), 
                            (NES > 0 & pval <= 0.05 ~ "Positive (p < 0.05)"),
                            (NES < 0 & pval > 0.05 ~ "Negative (N.S.)"), 
                            (NES > 0 & pval > 0.05 ~ "Positive (N.S.)")),
      face = case_when((NES < 0 & pval <= 0.05 ~ "plain"), 
                       (NES > 0 & pval <= 0.05 ~ "plain"),
                       (NES < 0 & pval > 0.05 ~ "plain"), 
                       (NES > 0 & pval > 0.05 ~ "plain")),
      pathway_clean = str_remove(pathway, "^KEGG_"), 
      pathway_clean = str_remove(pathway_clean, "^REACTOME_"), 
      pathway_clean = str_remove(pathway_clean, "^GOBP_"), 
      pathway_clean = str_remove(pathway_clean, "^GOMF_"), 
      pathway_clean = str_replace_all(pathway_clean, "_", " "),
      pathway_clean = str_to_sentence(pathway_clean),
      pathway_clean = str_replace_all(pathway_clean, "\\bi\\b", "I"),
      pathway_clean = str_replace_all(pathway_clean, "\\bii\\b", "II"),
      pathway_clean = str_replace_all(pathway_clean, "\\biii\\b", "III"),
      pathway_clean = str_replace_all(pathway_clean, "\\biv\\b", "IV"),
      pathway_clean = str_replace_all(pathway_clean, "\\bv\\b", "V"),
      pathway_clean = str_replace_all(pathway_clean, regex("\\(immune\\)", ignore_case = TRUE), "(IMMUNE)"),
      pathway_clean = capitalize_acronyms(pathway_clean, acronyms),
      pathway_clean = replace_mixed_case(pathway_clean, special_mixed, special_replacements),
      pathway_clean = paste0(pathway_clean, " (", size, ")")
    ) %>%
    arrange(pval)
  
  fgsea_res$pathway_clean <- reorder(fgsea_res$pathway_clean, -abs(fgsea_res$NES))
  
  fgsea_res %>%
    ggplot(aes(x = abs(NES), y = fct_rev(pathway_clean), label = pathway_clean)) +
    geom_point(aes(size = -log10(pval), color = direction, alpha = 0.8)) +
    # geom_vline(xintercept = 2, linetype = "dashed") +
    geom_text(aes(group = pathway_clean, color = direction, fontface = face), 
              hjust = 0, size = text1, nudge_x = xnudge) +
    scale_size_binned() +
    scale_color_manual(values = c("Positive (p < 0.05)" = "#963b32", "Negative (p < 0.05)" = "#1f566f", 
                                  "Positive (N.S.)" = "#f1bcb6", "Negative (N.S.)" = "#b3d6e5")) +
    # scale_color_manual(values = c("Positive" = "#963b32", "Negative" = "#1f566f",
    #                               "Positive p > 0.05" = "#f1bcb6", "Negative p > 0.05" = "#b3d6e5")) +
    scale_x_continuous(limits = c(0, xlimit), expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(add = 1)) +
    labs(
      x = "NES",
      y = "Pathways",
      color = "Direction",
      size = "-log(p-value)",
      title = title
    ) +
    guides(alpha = "none") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = text3),
      axis.title = element_text(size = text3),
      axis.ticks.y = element_blank(), 
      legend.position = c(0.9, 0.2),
      legend.background = element_blank(),
      legend.box.background = element_rect(color = "black"),
      legend.title = element_text(size = text2),
      legend.text = element_text(size = text2),
      title = element_text(size = text3)
    )
}
