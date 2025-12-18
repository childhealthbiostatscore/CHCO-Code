library(tidyverse)
library(data.table)

# Set up environment for Kopah
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else {
  stop("Unknown user: please specify root path for this user.")
}

# Target genes to search for
target_genes <- c(
  'LSAMP', 'LMF1', 'KCNK6', 'EMC10', 'NES', 'FAT1',
  'FLT4', 'MAP3K7IP3', 'IRF3', 'CASR', 'ACSM2B', 'ACER1', 'CA12', 'GLS',
  'CDH5', "PECAM1"
)

celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  aPT = "aPT",
  `PT-S1/S2` = "PT-S1/S2",
  `PT-S3` = "PT-S3",
  `PT-1` = "PT-1",
  `PT-2` = "PT-2",
  `PT-3` = "PT-3",
  `PT-4` = "PT-4",
  `PT-5` = "PT-5",
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  `C-TAL-1` = "C-TAL-1",
  `C-TAL-2` = "C-TAL-2",
  aTAL = "aTAL",
  dTAL = "dTAL",
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A"),
  `EC-AVR` = "EC-AVR",
  `EC-GC`  = "EC-GC",
  `EC-PTC` = "EC-PTC",
  `EC-AEA` = "EC-AEA",
  `EC-LYM` = "EC-LYM",
  `EC/VSMC` = "EC/VSMC",
  Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  Immune_Myeloid = c("MAC", "MON", "cDC", "pDC"),
  Immune_Lymphoid = c("CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),
  DCT_CNT = c("DCT", "dDCT", "CNT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",
  PEC = "PEC",
  Schwann = "SchwannCells"
)

# Function to search for genes in a file
search_genes_in_file <- function(filename, cell_name, target_genes) {
  if (!file.exists(filename)) {
    return(NULL)
  }
  
  tryCatch({
    df <- fread(filename)
    
    # Search for target genes (case-insensitive)
    results <- df %>%
      filter(toupper(gene) %in% toupper(target_genes) | 
               toupper(Gene) %in% toupper(target_genes)) %>%
      mutate(
        Gene_Name = case_when(
          toupper(gene) %in% toupper(target_genes) ~ 
            target_genes[match(toupper(gene), toupper(target_genes))],
          toupper(Gene) %in% toupper(target_genes) ~ 
            target_genes[match(toupper(Gene), toupper(target_genes))],
          TRUE ~ gene
        ),
        Cell_Type = cell_name,
        Significant = p_glp_t2dobGLP_Y < 0.05,
        Direction = ifelse(logFC_glp_t2dobGLP_Y > 0, "Up", "Down")
      ) %>%
      dplyr::select(
        Gene = Gene_Name,
        Cell_Type,
        logFC = logFC_glp_t2dobGLP_Y,
        p_value = p_glp_t2dobGLP_Y,
        FDR = fdr,
        Significant,
        Direction,
        Description = description
      )
    
    return(results)
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", filename, e$message))
    return(NULL)
  })
}

# Search all files
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("COMPREHENSIVE GENE EXPRESSION ANALYSIS ACROSS KIDNEY CELL TYPES\n")
cat("Comparing DKD patients with/without GLP-1 receptor agonist treatment (KPMP)\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

all_results <- list()

for (file_prefix in names(celltype_groups)) {
  home_path <- file.path(root_path, "Renal HERITAGE/Results/nebula/csv/full")
  filename <- sprintf("%s_dkd_30_glpy_glpn_kpmp.csv", tolower(file_prefix))
  filepath <- file.path(home_path, filename)
  cell_name <- file_prefix
  
  results <- search_genes_in_file(filepath, cell_name, target_genes)
  if (!is.null(results) && nrow(results) > 0) {
    all_results[[file_prefix]] <- results
  }
}

# Combine all results
if (length(all_results) > 0) {
  results_df <- bind_rows(all_results)
  
  # Summary of findings
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("SUMMARY: GENES FOUND ACROSS ALL CELL TYPES\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")
  
  genes_found <- unique(results_df$Gene)
  cat(sprintf("Total genes found: %d out of %d\n", 
              length(genes_found), length(target_genes)))
  cat(sprintf("Genes found: %s\n", paste(sort(genes_found), collapse = ", ")))
  
  genes_not_found <- setdiff(target_genes, genes_found)
  if (length(genes_not_found) > 0) {
    cat(sprintf("\nGenes NOT found in any cell type: %s\n", 
                paste(sort(genes_not_found), collapse = ", ")))
  }
  
  # Significant findings
  significant_df <- results_df %>% filter(Significant)
  
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("SIGNIFICANT DIFFERENTIAL EXPRESSION (p_value < 0.05)\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")
  
  if (nrow(significant_df) > 0) {
    cat(sprintf("Total significant findings: %d\n\n", nrow(significant_df)))
    
    # Group by gene
    cat("--- By Gene ---\n")
    for (gene in sort(unique(significant_df$Gene))) {
      cat(sprintf("\n%s:\n", gene))
      gene_data <- significant_df %>% filter(Gene == gene)
      for (i in 1:nrow(gene_data)) {
        row <- gene_data[i,]
        direction <- ifelse(row$logFC > 0, "↑", "↓")
        cat(sprintf("  • %s: logFC=%.3f %s, p_value=%.4f\n", 
                    row$Cell_Type, row$logFC, direction, row$p_value))
      }
    }
    
    # Group by cell type
    cat("\n--- By Cell Type ---\n")
    cell_summary <- significant_df %>%
      group_by(Cell_Type) %>%
      summarise(Genes = paste(sort(unique(Gene)), collapse = ", "))
    
    for (i in 1:nrow(cell_summary)) {
      cat(sprintf("\n%s: %s\n", cell_summary$Cell_Type[i], cell_summary$Genes[i]))
    }
  } else {
    cat("No genes reached P < 0.05 significance threshold\n")
  }
  
  # All findings (including non-significant)
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("ALL DIFFERENTIAL EXPRESSION RESULTS\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")
  
  # Create detailed table grouped by gene
  for (gene in sort(unique(results_df$Gene))) {
    cat(sprintf("\n%s\n", gene))
    cat(paste0(rep("-", 60), collapse = ""), "\n")
    cat(sprintf("%-40s %8s %10s %8s %3s\n", 
                "Cell Type", "logFC", "p-value", "FDR", "Sig"))
    cat(paste0(rep("-", 60), collapse = ""), "\n")
    
    gene_data <- results_df %>% 
      filter(Gene == gene) %>%
      arrange(FDR)
    
    for (i in 1:nrow(gene_data)) {
      row <- gene_data[i,]
      sig_marker <- ifelse(row$Significant, "***", "")
      cat(sprintf("%-40s %8.3f %10.2e %8.4f %3s\n",
                  row$Cell_Type, row$logFC, row$p_value, row$FDR, sig_marker))
    }
  }
  
  # Pathway-relevant findings
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("PATHWAY-RELEVANT FINDINGS\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  
  # Vascular/Endothelial genes
  vascular_genes <- c('FLT4', 'NES', 'FAT1')
  vascular_patterns <- c('EC', 'VSMC', 'Glomerular', 'Peritubular', 'Arteriole', 'Endothelial')
  
  cat("\n--- Vascular Endothelial Dysfunction ---\n")
  for (gene in vascular_genes) {
    if (gene %in% results_df$Gene) {
      gene_data <- results_df %>% filter(Gene == gene)
      relevant <- gene_data %>%
        filter(grepl(paste(vascular_patterns, collapse = "|"), Cell_Type, ignore.case = TRUE))
      
      if (nrow(relevant) > 0) {
        cat(sprintf("\n%s:\n", gene))
        for (i in 1:nrow(relevant)) {
          row <- relevant[i,]
          sig <- ifelse(row$Significant, "**SIG**", "")
          cat(sprintf("  • %s: logFC=%.3f, FDR=%.4f %s\n",
                      row$Cell_Type, row$logFC, row$FDR, sig))
        }
      }
    }
  }
  
  # Metabolic genes
  metabolic_genes <- c('LMF1', 'ACSM2B', 'ACER1', 'GLS', 'CA12')
  metabolic_patterns <- c('PT', 'Proximal', 'Principal', 'Intercalated', 'Podocyte')
  
  cat("\n--- Metabolism ---\n")
  for (gene in metabolic_genes) {
    if (gene %in% results_df$Gene) {
      gene_data <- results_df %>% filter(Gene == gene)
      relevant <- gene_data %>%
        filter(grepl(paste(metabolic_patterns, collapse = "|"), Cell_Type, ignore.case = TRUE))
      
      if (nrow(relevant) > 0) {
        cat(sprintf("\n%s:\n", gene))
        for (i in 1:nrow(relevant)) {
          row <- relevant[i,]
          sig <- ifelse(row$Significant, "**SIG**", "")
          cat(sprintf("  • %s: logFC=%.3f, FDR=%.4f %s\n",
                      row$Cell_Type, row$logFC, row$FDR, sig))
        }
      }
    }
  }
  
  # Inflammatory genes
  inflammatory_genes <- c('MAP3K7IP3', 'IRF3')
  immune_patterns <- c('Immune', 'Mesangial', 'Fibroblast', 'Lymphoid', 'Myeloid')
  
  cat("\n--- Fibroinflammation ---\n")
  for (gene in inflammatory_genes) {
    if (gene %in% results_df$Gene) {
      gene_data <- results_df %>% filter(Gene == gene)
      relevant <- gene_data %>%
        filter(grepl(paste(immune_patterns, collapse = "|"), Cell_Type, ignore.case = TRUE))
      
      if (nrow(relevant) > 0) {
        cat(sprintf("\n%s:\n", gene))
        for (i in 1:nrow(relevant)) {
          row <- relevant[i,]
          sig <- ifelse(row$Significant, "**SIG**", "")
          cat(sprintf("  • %s: logFC=%.3f, FDR=%.4f %s\n",
                      row$Cell_Type, row$logFC, row$FDR, sig))
        }
      }
    }
  }
  
  # Export results
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("EXPORTING RESULTS\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")
  
  # Save to CSV
  output_file <- 'gene_expression_analysis_results.csv'
  write.csv(results_df %>% dplyr::select(-Description), file.path(root_path, "Renal HERITAGE/Results/Figures/", output_file), row.names = FALSE)
  cat(sprintf("Results exported to: %s\n", output_file))
  
  # Create pivot table for heatmap
  pivot_df <- results_df %>%
    dplyr::select(Gene, Cell_Type, logFC) %>%
    pivot_wider(names_from = Cell_Type, values_from = logFC)
  
  heatmap_file <- 'gene_expression_heatmap_data.csv'
  write.csv(pivot_df, file.path(root_path, "Renal HERITAGE/Results/Figures/", heatmap_file), row.names = FALSE)
  cat(sprintf("Heatmap data exported to: %s\n", heatmap_file))
  
  # Statistical summary
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("STATISTICAL SUMMARY\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")
  
  cat(sprintf("Total gene-cell type combinations tested: %d\n", nrow(results_df)))
  cat(sprintf("Significant findings (P < 0.05): %d\n", nrow(significant_df)))
  cat(sprintf("Percentage significant: %.1f%%\n", 
              nrow(significant_df)/nrow(results_df)*100))
  
  # Cell types with most differential expression
  cell_type_counts <- results_df %>%
    group_by(Cell_Type) %>%
    summarise(
      Total_Genes = n(),
      Significant_Genes = sum(Significant)
    ) %>%
    arrange(desc(Total_Genes))
  
  cat("\nCell types with most genes detected:\n")
  for (i in 1:min(10, nrow(cell_type_counts))) {
    row <- cell_type_counts[i,]
    cat(sprintf("  %s: %d genes (%d significant)\n",
                row$Cell_Type, row$Total_Genes, row$Significant_Genes))
  }
  
  # Create visualization-ready summary
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("CREATING VISUALIZATION DATA\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")
  
  # Summary for plotting
  plot_data <- results_df %>%
    mutate(
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      Label = sprintf("%.2f%s", logFC, Significance)
    )
  
  # Save plot-ready data
  plot_file <- 'gene_expression_plot_data.csv'
  write.csv(plot_data, file.path(root_path, "Renal HERITAGE/Results/Figures/", plot_file), row.names = FALSE)
  cat(sprintf("Plot-ready data exported to: %s\n", plot_file))
  
} else {
  cat("\nNo target genes found in any of the provided files.\n")
  cat("Please check that the gene names match exactly with those in your data.\n")
}

cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")

# Optional: Create a heatmap if you have the necessary packages
if (require(pheatmap, quietly = TRUE) && exists("pivot_df")) {
  cat("\nCreating heatmap visualization...\n")
  
  # Prepare matrix for heatmap
  heatmap_matrix <- as.matrix(pivot_df[,-1])
  rownames(heatmap_matrix) <- pivot_df$Gene
  
  # Replace NA with 0 for visualization
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Create heatmap
  pdf(file.path(root_path, "Renal HERITAGE/Results/Figures/gene_expression_heatmap.pdf"), width = 12, height = 8)
  pheatmap(
    heatmap_matrix,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "Gene Expression Across Kidney Cell Types\n(GLP-1 RA Treatment Effect)",
    fontsize_row = 10,
    fontsize_col = 8,
    angle_col = 45,
    border_color = NA
  )
  dev.off()
  
  cat("Heatmap saved as: gene_expression_heatmap.pdf\n")
}
