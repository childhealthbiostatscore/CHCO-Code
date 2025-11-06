# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Full list of 141 genes from bulk analysis
bulk_141_genes <- readxl::read_xlsx('Projects/Pierre/intercept_rna_NFS4_SIND_vs_MIND_protein_coding.xlsx',
                                    sheet = 1)
bulk_141_genes <- head(bulk_141_genes$gene_name, 142)

seurat_obj <- readRDS('../OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/scRNA/data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds')


seurat_obj$celltype1 <- case_when(grepl("PT-",seurat_obj$celltype_rpca)~"PT",
                                  grepl("TAL-",seurat_obj$celltype_rpca)~"TAL",
                                  grepl("EC-",seurat_obj$celltype_rpca)~"EC",
                                  grepl("POD",seurat_obj$celltype_rpca)~"POD",
                                  grepl("MAC",seurat_obj$celltype_rpca)~"MAC",
                                  grepl("MON",seurat_obj$celltype_rpca)~"MON",
                                  grepl("PC-",seurat_obj$celltype_rpca)~"PC",
                                  grepl("FIB",seurat_obj$celltype_rpca)~"FIB_MC_VSMC",
                                  grepl("DTL",seurat_obj$celltype_rpca)~"DTL",
                                  seurat_obj$celltype_rpca=="DCT"~"DCT",
                                  seurat_obj$celltype_rpca=="ATL"~"ATL",
                                  seurat_obj$celltype_rpca=="B"~"B",
                                  seurat_obj$celltype_rpca=="T"~"T")
seurat_obj$celltype1 <- as.character(seurat_obj$celltype1)

seurat_obj$KPMP_celltype2 <- as.character(seurat_obj$KPMP_celltype)
seurat_obj$celltype2 <- ifelse(seurat_obj$KPMP_celltype=="aPT" | 
                                 seurat_obj$KPMP_celltype=="PT-S1/S2" | 
                                 seurat_obj$KPMP_celltype == "PT-S3","PT",
                               ifelse(grepl("TAL",seurat_obj$KPMP_celltype),"TAL",
                                      ifelse(grepl("EC-",seurat_obj$KPMP_celltype),"EC",seurat_obj$KPMP_celltype2)))












# Check which genes are present in your Seurat object
genes_in_data <- bulk_141_genes[bulk_141_genes %in% rownames(seurat_obj)]
genes_missing <- bulk_141_genes[!bulk_141_genes %in% rownames(seurat_obj)]

cat("Genes found:", length(genes_in_data), "\n")
cat("Genes missing:", length(genes_missing), "\n")

# Subset to all cell types (not just monocytes)
# This answers "where they are found" across the entire dataset

# 1. Dot plot across ALL cell types

pdf('Projects/Pierre/GeneExpression_by_Celltype.pdf', width = 18, height = 20)
DotPlot(
  seurat_obj,
  features = genes_in_data[1:141],  # Show first 30 genes
  group.by = "celltype2"  # Main cell type
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Heatmap showing expression across both annotation levels
avg_exp_celltype2 <- AverageExpression(
  seurat_obj,
  features = genes_in_data,
  group.by = "celltype2"
)$RNA

avg_exp_KPMP <- AverageExpression(
  seurat_obj,
  features = genes_in_data,
  group.by = "KPMP_celltype"
)$RNA

# 3. Find which cell types express each gene
gene_celltype_table <- data.frame()

for (gene in genes_in_data) {
  # For each gene, find cell types with >10% expression
  exp_by_celltype <- FetchData(seurat_obj, vars = c(gene, "celltype2", "KPMP_celltype"))
  
  summary_broad <- exp_by_celltype %>%
    group_by(celltype2) %>%
    summarise(
      pct_expressing = mean(get(gene) > 0) * 100,
      mean_exp = mean(get(gene))
    ) %>%
    filter(pct_expressing > 10)
  
  if (nrow(summary_broad) > 0) {
    summary_broad$gene <- gene
    gene_celltype_table <- rbind(gene_celltype_table, summary_broad)
  }
}
dev.off()


top6_genes <- genes_in_data[1:min(10)]

p3 <- FeaturePlot(
  seurat_obj,
  features = top6_genes,
  ncol = 3,
  cols = c("lightgrey", "red")
) &
  theme_minimal() &
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave(paste0("Projects/Pierre/featureplot_ListOne.pdf"),
       p3, width = 15, height = 25)



# Save results


write.csv(gene_celltype_table, "Projects/Pierre/gene_expression_by_celltype.csv")





### Second file

gene_list <- readxl::read_xlsx('Projects/Pierre/scRNA_2025_mono_DEG_SIND_MIND.xlsx',
                               sheet = 2)


# Get all unique genes from the DEG analysis
genes_to_check <- unique(gene_list$gene)

# Check Gene Expression Across Cell Types - By Cluster
# Loop through each monocyte cluster and check gene expression

library(Seurat)
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(pheatmap)

# ============================================================================
# 1. DOT PLOT - Expression Across Cell Types (celltype2)
# ============================================================================
output_dir <- "Projects/Pierre/"

# Read DEG data from Excel (sheet 2 = p0.05)
deg_data <- read_excel("Projects/Pierre/scRNA_2025_mono_DEG_SIND_MIND.xlsx", sheet = "p0.05")

# Get unique clusters
clusters <- unique(deg_data$cluster)
cat("Found", length(clusters), "monocyte clusters in DEG data:\n")
print(clusters)

# ============================================================================
# LOOP THROUGH EACH CLUSTER
# ============================================================================

# Store all results
all_results <- list()

for (cluster_name in clusters) {
  
  cat("\n========================================\n")
  cat("Processing:", cluster_name, "\n")
  cat("========================================\n")
  
  # Get genes for this cluster
  cluster_genes <- deg_data %>%
    filter(cluster == cluster_name) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    pull(gene) %>%
    unique()
  
  # Filter to genes in Seurat object
  cluster_genes_in_data <- cluster_genes[cluster_genes %in% rownames(seurat_obj)]
  
  cat("Cluster has", length(cluster_genes), "DEGs\n")
  cat("Found", length(cluster_genes_in_data), "genes in Seurat object\n")
  
  # Get top 20 genes for visualization
  top20_genes <- cluster_genes_in_data[]
  
  # ============================================================================
  # 1. DOT PLOT - Top genes across all cell types (celltype2)
  # ============================================================================
  
  if (length(top20_genes) > 0) {
    
    p1 <- DotPlot(
      seurat_obj,
      features = top20_genes,
      group.by = "celltype2",
      dot.scale = 8
    ) +
      coord_flip() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
      labs(title = paste(cluster_name, "- Top DEGs Across All Cell Types"),
           x = "Gene", y = "Cell Type (celltype2)")
    
    ggsave(file.path(output_dir, paste0("dotplot_", gsub(" ", "_", cluster_name), "_celltype2.pdf")),
           p1, width = 30, height = 49)
    
    # ============================================================================
    # 2. DOT PLOT - Top genes across KPMP cell types
    # ============================================================================
    
    p2 <- DotPlot(
      seurat_obj,
      features = top20_genes,
      group.by = "KPMP_celltype",
      dot.scale = 8
    ) +
      coord_flip() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(title = paste(cluster_name, "- Top DEGs Across KPMP Cell Types"),
           x = "Gene", y = "Cell Type (KPMP_celltype)")
    
    ggsave(file.path(output_dir, paste0("dotplot_", gsub(" ", "_", cluster_name), "_KPMP.pdf")),
           p2, width = 30, height = 49)
    
    # ============================================================================
    # 3. EXPRESSION SUMMARY TABLE (OPTIMIZED)
    # ============================================================================
    
    # Fetch all gene data at once (MUCH faster than looping)
    cat("Fetching expression data for", length(cluster_genes_in_data), "genes...\n")
    exp_data_all <- FetchData(seurat_obj,
                              vars = c(cluster_genes_in_data, "celltype2", "KPMP_celltype"))
    
    cat("Calculating expression statistics...\n")
    expression_summary <- data.frame()
    
    # Process each gene using the pre-fetched data
    for (gene in cluster_genes_in_data) {
      
      # Summary by celltype2
      summary_ct2 <- exp_data_all %>%
        group_by(celltype2) %>%
        summarise(
          pct_expressing = mean(.data[[gene]] > 0) * 100,
          mean_expression = mean(.data[[gene]]),
          mean_expression_positive = mean(.data[[gene]][.data[[gene]] > 0], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        filter(pct_expressing > 5) %>%
        arrange(desc(pct_expressing)) %>%
        mutate(
          gene = gene,
          cluster = cluster_name,
          annotation = "celltype2"
        )
      
      expression_summary <- bind_rows(expression_summary, summary_ct2)
      
      # Summary by KPMP_celltype
      summary_kpmp <- exp_data_all %>%
        group_by(KPMP_celltype) %>%
        summarise(
          pct_expressing = mean(.data[[gene]] > 0) * 100,
          mean_expression = mean(.data[[gene]]),
          mean_expression_positive = mean(.data[[gene]][.data[[gene]] > 0], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        filter(pct_expressing > 5) %>%
        arrange(desc(pct_expressing)) %>%
        mutate(
          gene = gene,
          cluster = cluster_name,
          annotation = "KPMP_celltype"
        ) %>%
        rename(celltype2 = KPMP_celltype)  # Rename for consistency
      
      expression_summary <- bind_rows(expression_summary, summary_kpmp)
    }
    
    # Save table for this cluster
    write.csv(expression_summary,
              file.path(output_dir, paste0("expression_", gsub(" ", "_", cluster_name), ".csv")),
              row.names = FALSE)
    
    # Store in list
    all_results[[cluster_name]] <- expression_summary
    
    # ============================================================================
    # 4. PRINT SUMMARY TO CONSOLE
    # ============================================================================
    
    cat("\n--- Top 10 genes and where they're found ---\n")
    top10_location <- expression_summary %>%
      filter(annotation == "celltype2") %>%
      group_by(gene) %>%
      slice_head(n = 1) %>%
      arrange(desc(pct_expressing)) %>%
      slice_head(n = 10)
    
    for (i in 1:nrow(top10_location)) {
      cat(sprintf("%d. %s: %s (%.1f%% expressing)\n",
                  i,
                  top10_location$gene[i],
                  top10_location$celltype2[i],
                  top10_location$pct_expressing[i]))
    }
    
    # ============================================================================
    # 5. FEATURE PLOTS - Top 6 genes
    # ============================================================================
    
    top6_genes <- top20_genes[1:min(6, length(top20_genes))]
    
    p3 <- FeaturePlot(
      seurat_obj,
      features = top6_genes,
      ncol = 3,
      cols = c("lightgrey", "red")
    ) &
      theme_minimal() &
      theme(plot.title = element_text(size = 10, face = "bold"))
    
    ggsave(file.path(output_dir, paste0("featureplot_", gsub(" ", "_", cluster_name), ".pdf")),
           p3, width = 15, height = 10)
    
    # ============================================================================
    # 6. HEATMAP - Percentage expressing (celltype2)
    # ============================================================================
    
    pct_matrix <- expression_summary %>%
      filter(annotation == "celltype2") %>%
      select(gene, celltype2, pct_expressing) %>%
      pivot_wider(names_from = celltype2, values_from = pct_expressing, values_fill = 0)
    
    if (nrow(pct_matrix) > 1) {
      pct_mat <- as.matrix(pct_matrix[, -1])
      rownames(pct_mat) <- pct_matrix$gene
      
      pheatmap(
        pct_mat,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = colorRampPalette(c("white", "yellow", "red"))(100),
        main = paste(cluster_name, "- % Cells Expressing"),
        fontsize_row = 7,
        fontsize_col = 9,
        filename = file.path(output_dir, paste0("heatmap_", gsub(" ", "_", cluster_name), "_celltype2.pdf")),
        width = 30,
        height = 49
      )
    }
    
  } else {
    cat("No genes found in Seurat object for this cluster\n")
  }
}

# ============================================================================
# COMBINE ALL RESULTS
# ============================================================================

# Combine all clusters into one table
all_expression_data <- bind_rows(all_results)

write.csv(all_expression_data,
          file.path(output_dir, "all_clusters_expression_summary.csv"),
          row.names = FALSE)

# ============================================================================
# SUMMARY ACROSS ALL CLUSTERS
# ============================================================================

cat("\n\n========================================\n")
cat("SUMMARY ACROSS ALL CLUSTERS\n")
cat("========================================\n\n")

# Count genes per cluster
genes_per_cluster <- deg_data %>%
  group_by(cluster) %>%
  summarise(
    total_DEGs = n(),
    in_seurat = sum(gene %in% rownames(seurat_obj)),
    pct_in_seurat = round(in_seurat/total_DEGs * 100, 1)
  )

print(genes_per_cluster)

# Most commonly expressed cell types across all cluster DEGs
top_celltypes <- all_expression_data %>%
  filter(annotation == "celltype2") %>%
  group_by(celltype2) %>%
  summarise(
    n_genes = n_distinct(gene),
    avg_pct_expressing = mean(pct_expressing)
  ) %>%
  arrange(desc(n_genes)) %>%
  slice_head(n = 10)

cat("\nTop 10 cell types expressing cluster DEGs:\n")
print(top_celltypes)

# Genes found in multiple cell types
multi_celltype_genes <- all_expression_data %>%
  filter(annotation == "celltype2") %>%
  group_by(gene, cluster) %>%
  summarise(n_celltypes = n(), .groups = "drop") %>%
  filter(n_celltypes >= 3) %>%
  arrange(desc(n_celltypes))

cat("\nGenes expressed in 3+ cell types:\n")
print(head(multi_celltype_genes, 20))

write.csv(multi_celltype_genes,
          file.path(output_dir, "genes_multiple_celltypes.csv"),
          row.names = FALSE)




