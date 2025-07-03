library(Seurat)
library(dplyr)

# Function to extract top marker genes from an annotated Seurat object
extract_top_markers <- function(reference_seurat, 
                                cell_type_column = "celltype", 
                                top_n = 5,
                                min_pct = 0.25,
                                logfc_threshold = 0.25,
                                test_use = "wilcox") {
  
  cat("Extracting marker genes from reference dataset...\n")
  
  # Set the identity to the cell type column
  Idents(reference_seurat) <- reference_seurat@meta.data[[cell_type_column]]
  
  # Check available cell types
  cell_types <- levels(Idents(reference_seurat))
  cat("Found cell types:", paste(cell_types, collapse = ", "), "\n")
  
  # Find markers for all cell types
  cat("Running differential expression analysis...\n")
  all_markers <- FindAllMarkers(
    reference_seurat,
    only.pos = TRUE,  # Only positive markers
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    test.use = test_use,
    verbose = FALSE
  )
  
  # Filter significant markers
  significant_markers <- all_markers %>%
    filter(p_val_adj < 0.05) %>%
    arrange(cluster, desc(avg_log2FC))
  
  cat("Found", nrow(significant_markers), "significant markers\n")
  
  # Extract top N markers for each cell type
  top_markers_list <- list()
  
  for (celltype in cell_types) {
    celltype_markers <- significant_markers %>%
      filter(cluster == celltype) %>%
      slice_head(n = top_n)
    
    if (nrow(celltype_markers) > 0) {
      top_markers_list[[celltype]] <- celltype_markers$gene
      cat("Cell type", celltype, ":", length(celltype_markers$gene), "markers\n")
    } else {
      cat("Warning: No significant markers found for", celltype, "\n")
    }
  }
  
  # Return both the list and the full markers table
  return(list(
    marker_genes = top_markers_list,
    full_markers = significant_markers,
    summary = significant_markers %>%
      group_by(cluster) %>%
      summarise(
        n_markers = n(),
        avg_logfc = mean(avg_log2FC),
        avg_pct1 = mean(pct.1),
        avg_pct2 = mean(pct.2),
        .groups = 'drop'
      )
  ))
}

# Function to visualize marker quality
plot_marker_quality <- function(marker_results) {
  
  # Plot 1: Number of markers per cell type
  p1 <- ggplot(marker_results$summary, aes(x = cluster, y = n_markers)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Number of Significant Markers per Cell Type",
         x = "Cell Type", y = "Number of Markers")
  
  # Plot 2: Average log2FC distribution
  p2 <- ggplot(marker_results$full_markers, aes(x = cluster, y = avg_log2FC)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Log2FC Distribution of Markers",
         x = "Cell Type", y = "Average Log2FC")
  
  # Plot 3: Marker specificity (pct.1 vs pct.2)
  p3 <- ggplot(marker_results$full_markers, aes(x = pct.2, y = pct.1)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~cluster, scales = "free") +
    theme_minimal() +
    labs(title = "Marker Specificity (% cells expressing)",
         x = "% in other cell types", y = "% in target cell type")
  
  return(list(count_plot = p1, logfc_plot = p2, specificity_plot = p3))
}


# Additional utility function to compare marker lists
compare_marker_lists <- function(list1, list2, list1_name = "List1", list2_name = "List2") {
  
  cat("Comparing marker gene lists:\n")
  
  common_celltypes <- intersect(names(list1), names(list2))
  
  if (length(common_celltypes) == 0) {
    cat("No common cell types found between lists\n")
    return(NULL)
  }
  
  comparison <- data.frame(
    celltype = character(),
    list1_markers = integer(),
    list2_markers = integer(),
    common_markers = integer(),
    jaccard_index = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (celltype in common_celltypes) {
    markers1 <- list1[[celltype]]
    markers2 <- list2[[celltype]]
    
    common <- intersect(markers1, markers2)
    union <- union(markers1, markers2)
    jaccard <- length(common) / length(union)
    
    comparison <- rbind(comparison, data.frame(
      celltype = celltype,
      list1_markers = length(markers1),
      list2_markers = length(markers2),
      common_markers = length(common),
      jaccard_index = jaccard,
      stringsAsFactors = FALSE
    ))
    
    cat("Cell type:", celltype, "\n")
    cat("  Common markers:", paste(common, collapse = ", "), "\n")
  }
  
  colnames(comparison)[2:3] <- paste0(c(list1_name, list2_name), "_markers")
  
  return(comparison)
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(reshape2)

# Function to create comprehensive cell type visualizations
visualize_celltypes <- function(seurat_obj, 
                                celltype_column = "celltype_annotation",
                                reduction = "umap",
                                marker_genes = NULL,
                                sample_column = NULL) {
  
  cat("Creating cell type visualizations...\n")
  
  # Check if the cell type column exists
  if (!celltype_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", celltype_column, "not found in metadata"))
  }
  
  # Get cell types and create color palette
  cell_types <- unique(seurat_obj@meta.data[[celltype_column]])
  cell_types <- cell_types[!is.na(cell_types)]
  n_types <- length(cell_types)
  
  # Create a nice color palette
  if (n_types <= 8) {
    colors <- brewer.pal(max(3, n_types), "Set2")
  } else if (n_types <= 12) {
    colors <- brewer.pal(n_types, "Set3")
  } else {
    colors <- rainbow(n_types)
  }
  names(colors) <- cell_types
  
  plots <- list()
  
  # 1. Basic UMAP with cell type annotations
  plots$main_umap <- DimPlot(seurat_obj, 
                             group.by = celltype_column,
                             reduction = reduction,
                             label = TRUE,
                             label.size = 3,
                             cols = colors,
                             pt.size = 0.5) +
    ggtitle("Cell Type Annotations") +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  # 2. UMAP without labels for cleaner view
  plots$clean_umap <- DimPlot(seurat_obj, 
                              group.by = celltype_column,
                              reduction = reduction,
                              label = FALSE,
                              cols = colors,
                              pt.size = 0.3) +
    ggtitle("Cell Type Annotations (Clean)") +
    theme_minimal()
  
  # 3. Cell type proportions
  prop_data <- seurat_obj@meta.data %>%
    count(!!sym(celltype_column)) %>%
    mutate(proportion = n / sum(n) * 100) %>%
    arrange(desc(n))
  
  plots$proportions <- ggplot(prop_data, aes(x = reorder(!!sym(celltype_column), n), y = n)) +
    geom_col(aes(fill = !!sym(celltype_column)), alpha = 0.8) +
    scale_fill_manual(values = colors) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Cell Type Counts",
         x = "Cell Type",
         y = "Number of Cells") +
    theme(legend.position = "none") +
    geom_text(aes(label = paste0(n, " (", round(proportion, 1), "%)")), 
              hjust = -0.1, size = 3)
  
  # 4. If sample information is available, show composition
  if (!is.null(sample_column) && sample_column %in% colnames(seurat_obj@meta.data)) {
    
    composition_data <- seurat_obj@meta.data %>%
      group_by(!!sym(sample_column), !!sym(celltype_column)) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(!!sym(sample_column)) %>%
      mutate(proportion = count / sum(count) * 100)
    
    plots$sample_composition <- ggplot(composition_data, 
                                       aes(x = !!sym(sample_column), 
                                           y = proportion, 
                                           fill = !!sym(celltype_column))) +
      geom_col(position = "stack") +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Cell Type Composition by Sample",
           x = "Sample",
           y = "Percentage",
           fill = "Cell Type")
  }
  
  # 5. Annotation confidence metrics (if available)
  if ("score_difference" %in% colnames(seurat_obj@meta.data)) {
    plots$confidence <- ggplot(seurat_obj@meta.data, 
                               aes(x = !!sym(celltype_column), y = score_difference)) +
      geom_boxplot(aes(fill = !!sym(celltype_column)), alpha = 0.7) +
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "Annotation Confidence by Cell Type",
           subtitle = "Higher values = more confident annotations",
           x = "Cell Type",
           y = "Score Difference")
  }
  
  return(plots)
}

# Function to visualize marker gene expression
plot_marker_expression <- function(seurat_obj, marker_genes, 
                                   celltype_column = "celltype_annotation",
                                   reduction = "umap") {
  
  cat("Creating marker gene expression plots...\n")
  
  # Flatten marker genes list
  all_markers <- unique(unlist(marker_genes))
  
  # Check which markers are present
  available_markers <- intersect(all_markers, rownames(seurat_obj))
  missing_markers <- setdiff(all_markers, rownames(seurat_obj))
  
  if (length(missing_markers) > 0) {
    cat("Missing markers:", paste(missing_markers, collapse = ", "), "\n")
  }
  
  plots <- list()
  
  # 1. Feature plots for key markers (up to 8)
  if (length(available_markers) > 0) {
    key_markers <- head(available_markers, 8)
    
    plots$feature_plots <- FeaturePlot(seurat_obj, 
                                       features = key_markers,
                                       reduction = reduction,
                                       ncol = 4,
                                       pt.size = 0.3) &
      theme_minimal() &
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
  }
  
  # 2. Violin plots for marker expression by cell type
  if (length(available_markers) > 0) {
    # Select top markers for violin plot (max 6 for readability)
    violin_markers <- head(available_markers, 6)
    
    plots$violin_plots <- VlnPlot(seurat_obj, 
                                  features = violin_markers,
                                  group.by = celltype_column,
                                  ncol = 3,
                                  pt.size = 0) &
      theme_minimal() &
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # 3. Heatmap of average marker expression
  if (length(available_markers) > 0 && length(marker_genes) > 1) {
    
    # Calculate average expression per cell type
    avg_exp <- AverageExpression(seurat_obj, 
                                 features = available_markers,
                                 group.by = celltype_column,
                                 slot = "data")$RNA
    
    # Scale the data
    avg_exp_scaled <- t(scale(t(avg_exp)))
    
    # Create heatmap data
    heatmap_data <- melt(avg_exp_scaled)
    colnames(heatmap_data) <- c("Gene", "CellType", "Expression")
    
    plots$heatmap <- ggplot(heatmap_data, aes(x = CellType, y = Gene, fill = Expression)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, name = "Scaled\nExpression") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8)) +
      labs(title = "Average Marker Expression by Cell Type",
           x = "Cell Type", y = "Marker Gene")
  }
  
  # 4. Dot plot for marker expression
  if (length(available_markers) > 0) {
    plots$dot_plot <- DotPlot(seurat_obj, 
                              features = available_markers,
                              group.by = celltype_column) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Marker Gene Expression (Dot Plot)")
  }
  
  return(plots)
}

# Function to create summary statistics
celltype_summary <- function(seurat_obj, celltype_column = "celltype_annotation") {
  
  cat("Generating cell type summary statistics...\n")
  
  summary_stats <- seurat_obj@meta.data %>%
    group_by(!!sym(celltype_column)) %>%
    summarise(
      n_cells = n(),
      percentage = round(n() / nrow(seurat_obj@meta.data) * 100, 2),
      avg_nFeature = round(mean(nFeature_RNA, na.rm = TRUE), 0),
      avg_nCount = round(mean(nCount_RNA, na.rm = TRUE), 0),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_cells))
  
  # Add confidence metrics if available
  if ("score_difference" %in% colnames(seurat_obj@meta.data)) {
    confidence_stats <- seurat_obj@meta.data %>%
      group_by(!!sym(celltype_column)) %>%
      summarise(
        avg_confidence = round(mean(score_difference, na.rm = TRUE), 3),
        low_confidence_cells = sum(score_difference < 0.1, na.rm = TRUE),
        .groups = 'drop'
      )
    
    summary_stats <- left_join(summary_stats, confidence_stats, 
                               by = celltype_column)
  }
  
  print(summary_stats)
  return(summary_stats)
}

# Main execution function
create_all_visualizations <- function(seurat_obj, 
                                      marker_genes = NULL,
                                      celltype_column = "celltype_annotation",
                                      reduction = "umap",
                                      sample_column = NULL,
                                      save_plots = FALSE,
                                      output_prefix = "celltype_viz") {
  
  cat("=== Creating comprehensive cell type visualizations ===\n")
  
  # 1. Basic cell type visualizations
  basic_plots <- visualize_celltypes(seurat_obj, 
                                     celltype_column = celltype_column,
                                     reduction = reduction,
                                     sample_column = sample_column)
  
  # 2. Marker gene expression plots
  marker_plots <- NULL
  if (!is.null(marker_genes)) {
    marker_plots <- plot_marker_expression(seurat_obj, 
                                           marker_genes = marker_genes,
                                           celltype_column = celltype_column,
                                           reduction = reduction)
  }
  
  # 3. Summary statistics
  summary_stats <- celltype_summary(seurat_obj, celltype_column = celltype_column)
  
  # Display key plots
  print(basic_plots$main_umap)
  print(basic_plots$proportions)
  
  if (!is.null(basic_plots$confidence)) {
    print(basic_plots$confidence)
  }
  
  if (!is.null(marker_plots$feature_plots)) {
    print(marker_plots$feature_plots)
  }
  
  if (!is.null(marker_plots$heatmap)) {
    print(marker_plots$heatmap)
  }
  
  # Save plots if requested
  if (save_plots) {
    cat("Saving plots...\n")
    
    ggsave(paste0(output_prefix, "_main_umap.pdf"), 
           basic_plots$main_umap, width = 10, height = 8)
    ggsave(paste0(output_prefix, "_proportions.pdf"), 
           basic_plots$proportions, width = 8, height = 6)
    
    if (!is.null(marker_plots$heatmap)) {
      ggsave(paste0(output_prefix, "_marker_heatmap.pdf"), 
             marker_plots$heatmap, width = 10, height = 8)
    }
    
    # Save summary statistics
    write.csv(summary_stats, paste0(output_prefix, "_summary_stats.csv"), 
              row.names = FALSE)
  }
  
  # Return all plots and stats
  return(list(
    basic_plots = basic_plots,
    marker_plots = marker_plots,
    summary_stats = summary_stats
  ))
}

# Example usage:
# Assuming you have your annotated Seurat object 'sc_org' and marker genes

# Create all visualizations
# all_viz <- create_all_visualizations(
#   seurat_obj = sc_org,
#   marker_genes = marker_genes,  # from previous extraction or manual list
#   celltype_column = "celltype_annotation",
#   reduction = "umap",
#   sample_column = "orig.ident",  # if you have sample information
#   save_plots = TRUE,
#   output_prefix = "kidney_celltypes"
# )

# Access individual plot types:
# all_viz$basic_plots$main_umap
# all_viz$marker_plots$heatmap
# all_viz$summary_stats

cat("Visualization functions loaded!\n")
cat("Use create_all_visualizations() to generate comprehensive cell type plots\n")