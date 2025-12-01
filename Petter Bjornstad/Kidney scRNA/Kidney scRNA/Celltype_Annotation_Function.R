library(Seurat)
library(dplyr)

# Function to annotate cell types based on marker genes
annotate_celltypes <- function(seurat_obj, marker_list, min_score_diff = 0.1) {
  
  # Calculate module scores for each cell type
  cat("Calculating module scores for each cell type...\n")
  
  for (celltype in names(marker_list)) {
    # Check which markers are present in the dataset
    available_markers <- intersect(marker_list[[celltype]], rownames(seurat_obj))
    
    if (length(available_markers) == 0) {
      warning(paste("No markers found for", celltype, "in the dataset"))
      next
    }
    
    if (length(available_markers) < length(marker_list[[celltype]])) {
      missing_markers <- setdiff(marker_list[[celltype]], available_markers)
      cat(paste("Warning: Missing markers for", celltype, ":", paste(missing_markers, collapse = ", "), "\n"))
    }
    
    # Add module score
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(available_markers),
      name = paste0(celltype, "_score"),
      ctrl = 100
    )
  }
  
  # Get all score column names
  score_cols <- paste0(names(marker_list), "_score1")
  score_cols <- score_cols[score_cols %in% colnames(seurat_obj@meta.data)]
  
  if (length(score_cols) == 0) {
    stop("No module scores were calculated successfully")
  }
  
  # Create a matrix of scores
  score_matrix <- seurat_obj@meta.data[, score_cols, drop = FALSE]
  colnames(score_matrix) <- gsub("_score1", "", colnames(score_matrix))
  
  # Assign cell types based on highest score
  cat("Assigning cell types based on highest scores...\n")
  
  # Find the cell type with maximum score for each cell
  max_scores <- apply(score_matrix, 1, max)
  second_max_scores <- apply(score_matrix, 1, function(x) sort(x, decreasing = TRUE)[2])
  score_differences <- max_scores - second_max_scores
  
  # Assign cell types
  assigned_celltypes <- colnames(score_matrix)[apply(score_matrix, 1, which.max)]
  
  # Set cells with low confidence to "Unknown"
  low_confidence <- score_differences < min_score_diff
  assigned_celltypes[low_confidence] <- "Unknown"
  
  # Add annotations to metadata
  seurat_obj@meta.data$celltype_annotation <- assigned_celltypes
  seurat_obj@meta.data$max_score <- max_scores
  seurat_obj@meta.data$score_difference <- score_differences
  
  # Print summary
  cat("\nCell type annotation summary:\n")
  print(table(assigned_celltypes))
  
  cat(paste("\nCells with low confidence (score diff <", min_score_diff, "):", sum(low_confidence), "\n"))
  
  return(seurat_obj)
}

# Function to visualize cell type annotations
plot_celltype_results <- function(seurat_obj, reduction = "umap") {
  
  # Plot 1: UMAP with cell type annotations
  p1 <- DimPlot(seurat_obj, 
                group.by = "celltype_annotation", 
                reduction = reduction,
                label = TRUE, 
                label.size = 3) +
    ggtitle("Cell Type Annotations") +
    theme(legend.position = "bottom")
  
  # Plot 2: Score distribution
  score_cols <- grep("_score1$", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(score_cols) > 0) {
    score_data <- seurat_obj@meta.data[, c("celltype_annotation", score_cols)]
    colnames(score_data) <- c("celltype_annotation", gsub("_score1", "", score_cols))
    
    # Melt the data for plotting
    score_long <- reshape2::melt(score_data, id.vars = "celltype_annotation")
    
    p2 <- ggplot(score_long, aes(x = variable, y = value, fill = celltype_annotation)) +
      geom_boxplot() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Module Score Distribution by Cell Type",
           x = "Cell Type Markers", 
           y = "Module Score") +
      facet_wrap(~celltype_annotation, scales = "free")
  } else {
    p2 <- NULL
  }
  
  # Plot 3: Score difference distribution
  p3 <- ggplot(seurat_obj@meta.data, aes(x = celltype_annotation, y = score_difference)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Score Difference Distribution",
         subtitle = "Red line shows confidence threshold",
         x = "Assigned Cell Type", 
         y = "Score Difference")
  
  return(list(annotation_plot = p1, score_plot = p2, confidence_plot = p3))
}

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


library(Seurat)
library(dplyr)
library(stringr)

# Function to parse marker gene dataframe and create annotation list
parse_marker_dataframe <- function(marker_df, 
                                   celltype_col = "Celltype", 
                                   marker_col = "Pos_Marker_Genes") {
  
  cat("Parsing marker gene dataframe...\n")
  
  # Check if required columns exist
  if (!celltype_col %in% colnames(marker_df)) {
    stop(paste("Column", celltype_col, "not found in dataframe"))
  }
  
  if (!marker_col %in% colnames(marker_df)) {
    stop(paste("Column", marker_col, "not found in dataframe"))
  }
  
  # Create marker gene list
  marker_genes <- list()
  
  for (i in 1:nrow(marker_df)) {
    celltype <- marker_df[[celltype_col]][i]
    marker_string <- marker_df[[marker_col]][i]
    
    # Parse comma-separated genes and clean whitespace
    genes <- str_split(marker_string, ",")[[1]]
    genes <- str_trim(genes)  # Remove leading/trailing whitespace
    genes <- genes[genes != ""]  # Remove empty strings
    
    # Store in list
    marker_genes[[celltype]] <- genes
    
    cat("Cell type", celltype, ":", length(genes), "markers\n")
  }
  
  return(marker_genes)
}

# Enhanced annotation function that works with dataframe input
annotate_from_dataframe <- function(seurat_obj, 
                                    marker_df,
                                    celltype_col = "Celltype",
                                    marker_col = "Pos_Marker_Genes",
                                    min_score_diff = 0.1,
                                    ctrl_genes = 100) {
  
  cat("=== Annotating Seurat object using marker dataframe ===\n")
  
  # Parse the dataframe to create marker gene list
  marker_genes <- parse_marker_dataframe(marker_df, celltype_col, marker_col)
  
  # Print summary of marker genes
  cat("\nMarker gene summary:\n")
  for (celltype in names(marker_genes)) {
    cat(paste0(celltype, ": ", paste(marker_genes[[celltype]], collapse = ", "), "\n"))
  }
  
  # Calculate module scores for each cell type
  cat("\nCalculating module scores for each cell type...\n")
  
  for (celltype in names(marker_genes)) {
    # Check which markers are present in the dataset
    available_markers <- intersect(marker_genes[[celltype]], rownames(seurat_obj))
    
    if (length(available_markers) == 0) {
      warning(paste("No markers found for", celltype, "in the dataset"))
      next
    }
    
    if (length(available_markers) < length(marker_genes[[celltype]])) {
      missing_markers <- setdiff(marker_genes[[celltype]], available_markers)
      cat(paste("Warning: Missing markers for", celltype, ":", paste(missing_markers, collapse = ", "), "\n"))
    }
    
    cat(paste("Using", length(available_markers), "markers for", celltype, "\n"))
    
    # Add module score
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(available_markers),
      name = paste0(celltype, "_score"),
      ctrl = ctrl_genes
    )
  }
  
  # Get all score column names
  score_cols <- paste0(names(marker_genes), "_score1")
  score_cols <- score_cols[score_cols %in% colnames(seurat_obj@meta.data)]
  
  if (length(score_cols) == 0) {
    stop("No module scores were calculated successfully")
  }
  
  # Create a matrix of scores
  score_matrix <- seurat_obj@meta.data[, score_cols, drop = FALSE]
  colnames(score_matrix) <- gsub("_score1", "", colnames(score_matrix))
  
  # Assign cell types based on highest score
  cat("\nAssigning cell types based on highest scores...\n")
  
  # Find the cell type with maximum score for each cell
  max_scores <- apply(score_matrix, 1, max)
  second_max_scores <- apply(score_matrix, 1, function(x) {
    if(length(x) > 1) {
      return(sort(x, decreasing = TRUE)[2])
    } else {
      return(0)
    }
  })
  score_differences <- max_scores - second_max_scores
  
  # Assign cell types
  assigned_celltypes <- colnames(score_matrix)[apply(score_matrix, 1, which.max)]
  
  # Set cells with low confidence to "Unknown"
  low_confidence <- score_differences < min_score_diff
  assigned_celltypes[low_confidence] <- "Unknown"
  
  # Add annotations to metadata
  seurat_obj@meta.data$celltype_annotation <- assigned_celltypes
  seurat_obj@meta.data$max_score <- max_scores
  seurat_obj@meta.data$score_difference <- score_differences
  
  # Print summary
  cat("\nCell type annotation summary:\n")
  annotation_table <- table(assigned_celltypes)
  print(annotation_table)
  
  cat(paste("\nCells with low confidence (score diff <", min_score_diff, "):", sum(low_confidence), "\n"))
  cat(paste("Annotation rate:", round(100 * sum(assigned_celltypes != "Unknown") / length(assigned_celltypes), 2), "%\n"))
  
  # Create detailed summary
  summary_df <- data.frame(
    CellType = names(annotation_table),
    Count = as.numeric(annotation_table),
    Percentage = round(as.numeric(annotation_table) / length(assigned_celltypes) * 100, 2)
  )
  
  cat("\nDetailed annotation summary:\n")
  print(summary_df)
  
  return(list(
    seurat_object = seurat_obj,
    marker_genes = marker_genes,
    annotation_summary = summary_df,
    score_matrix = score_matrix
  ))
}

# Function to validate marker genes against Seurat object
validate_markers <- function(seurat_obj, marker_df, 
                             celltype_col = "Celltype", 
                             marker_col = "Pos_Marker_Genes") {
  
  cat("=== Validating marker genes against Seurat object ===\n")
  
  # Parse markers
  marker_genes <- parse_marker_dataframe(marker_df, celltype_col, marker_col)
  
  # Check availability
  validation_results <- data.frame(
    CellType = character(),
    Total_Markers = integer(),
    Available_Markers = integer(),
    Missing_Markers = integer(),
    Availability_Rate = numeric(),
    Available_Genes = character(),
    Missing_Genes = character(),
    stringsAsFactors = FALSE
  )
  
  for (celltype in names(marker_genes)) {
    markers <- marker_genes[[celltype]]
    available <- intersect(markers, rownames(seurat_obj))
    missing <- setdiff(markers, rownames(seurat_obj))
    
    validation_results <- rbind(validation_results, data.frame(
      CellType = celltype,
      Total_Markers = length(markers),
      Available_Markers = length(available),
      Missing_Markers = length(missing),
      Availability_Rate = round(length(available) / length(markers) * 100, 1),
      Available_Genes = paste(available, collapse = ", "),
      Missing_Genes = paste(missing, collapse = ", "),
      stringsAsFactors = FALSE
    ))
  }
  
  cat("\nMarker validation results:\n")
  print(validation_results[, 1:5])  # Print summary columns
  
  # Identify problematic cell types
  low_availability <- validation_results$Availability_Rate < 50
  if (any(low_availability)) {
    cat("\nWARNING: Cell types with <50% marker availability:\n")
    print(validation_results[low_availability, c("CellType", "Availability_Rate", "Missing_Genes")])
  }
  
  return(validation_results)
}

library(Seurat)
library(dplyr)
library(stringr)

# Function to parse marker gene dataframe with both positive and negative markers
parse_marker_dataframe_dual <- function(marker_df, 
                                        celltype_col = "Celltype", 
                                        pos_marker_col = "Pos_Marker_Genes",
                                        neg_marker_col = "Neg_Marker_Genes") {
  
  cat("Parsing marker gene dataframe with positive and negative markers...\n")
  
  # Check if required columns exist
  if (!celltype_col %in% colnames(marker_df)) {
    stop(paste("Column", celltype_col, "not found in dataframe"))
  }
  
  if (!pos_marker_col %in% colnames(marker_df)) {
    stop(paste("Column", pos_marker_col, "not found in dataframe"))
  }
  
  # Check if negative marker column exists
  has_neg_markers <- neg_marker_col %in% colnames(marker_df)
  
  if (!has_neg_markers) {
    cat("Warning: Negative marker column", neg_marker_col, "not found. Using positive markers only.\n")
  }
  
  # Create marker gene lists
  pos_markers <- list()
  neg_markers <- list()
  
  for (i in 1:nrow(marker_df)) {
    celltype <- marker_df[[celltype_col]][i]
    
    # Parse positive markers
    pos_marker_string <- marker_df[[pos_marker_col]][i]
    if (!is.na(pos_marker_string) && pos_marker_string != "") {
      pos_genes <- str_split(pos_marker_string, ",")[[1]]
      pos_genes <- str_trim(pos_genes)  # Remove leading/trailing whitespace
      pos_genes <- pos_genes[pos_genes != ""]  # Remove empty strings
      pos_markers[[celltype]] <- pos_genes
    } else {
      pos_markers[[celltype]] <- character(0)
    }
    
    # Parse negative markers if available
    if (has_neg_markers) {
      neg_marker_string <- marker_df[[neg_marker_col]][i]
      if (!is.na(neg_marker_string) && neg_marker_string != "") {
        neg_genes <- str_split(neg_marker_string, ",")[[1]]
        neg_genes <- str_trim(neg_genes)
        neg_genes <- neg_genes[neg_genes != ""]
        neg_markers[[celltype]] <- neg_genes
      } else {
        neg_markers[[celltype]] <- character(0)
      }
    }
    
    cat("Cell type", celltype, ":", length(pos_markers[[celltype]]), "positive markers")
    if (has_neg_markers) {
      cat(",", length(neg_markers[[celltype]]), "negative markers")
    }
    cat("\n")
  }
  
  return(list(
    positive_markers = pos_markers,
    negative_markers = if(has_neg_markers) neg_markers else NULL,
    has_negative_markers = has_neg_markers
  ))
}

# Enhanced annotation function using both positive and negative markers
annotate_from_dataframe_dual <- function(seurat_obj, 
                                         marker_df,
                                         celltype_col = "Celltype",
                                         pos_marker_col = "Pos_Marker_Genes",
                                         neg_marker_col = "Neg_Marker_Genes",
                                         min_score_diff = 0.1,
                                         neg_weight = 0.5,
                                         ctrl_genes = 100) {
  
  cat("=== Annotating Seurat object using positive and negative markers ===\n")
  
  # Parse the dataframe to create marker gene lists
  marker_data <- parse_marker_dataframe_dual(marker_df, celltype_col, pos_marker_col, neg_marker_col)
  
  pos_markers <- marker_data$positive_markers
  neg_markers <- marker_data$negative_markers
  has_neg_markers <- marker_data$has_negative_markers
  
  # Print summary of marker genes
  cat("\nMarker gene summary:\n")
  for (celltype in names(pos_markers)) {
    cat(paste0(celltype, " - Positive: ", paste(pos_markers[[celltype]], collapse = ", ")))
    if (has_neg_markers && length(neg_markers[[celltype]]) > 0) {
      cat(paste0(" | Negative: ", paste(neg_markers[[celltype]], collapse = ", ")))
    }
    cat("\n")
  }
  
  # Calculate positive module scores
  cat("\nCalculating positive module scores...\n")
  pos_score_cols <- c()
  
  for (celltype in names(pos_markers)) {
    if (length(pos_markers[[celltype]]) == 0) {
      cat("Warning: No positive markers for", celltype, "\n")
      next
    }
    
    # Check which positive markers are present
    available_pos_markers <- intersect(pos_markers[[celltype]], rownames(seurat_obj))
    
    if (length(available_pos_markers) == 0) {
      warning(paste("No positive markers found for", celltype, "in the dataset"))
      next
    }
    
    if (length(available_pos_markers) < length(pos_markers[[celltype]])) {
      missing_pos_markers <- setdiff(pos_markers[[celltype]], available_pos_markers)
      cat(paste("Warning: Missing positive markers for", celltype, ":", 
                paste(missing_pos_markers, collapse = ", "), "\n"))
    }
    
    cat(paste("Using", length(available_pos_markers), "positive markers for", celltype, "\n"))
    
    # Add positive module score
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(available_pos_markers),
      name = paste0(celltype, "_pos_score"),
      ctrl = ctrl_genes
    )
    
    pos_score_cols <- c(pos_score_cols, paste0(celltype, "_pos_score1"))
  }
  
  # Calculate negative module scores if available
  neg_score_cols <- c()
  
  if (has_neg_markers) {
    cat("\nCalculating negative module scores...\n")
    
    for (celltype in names(neg_markers)) {
      if (length(neg_markers[[celltype]]) == 0) {
        cat("No negative markers for", celltype, "\n")
        next
      }
      
      # Check which negative markers are present
      available_neg_markers <- intersect(neg_markers[[celltype]], rownames(seurat_obj))
      
      if (length(available_neg_markers) == 0) {
        cat("Warning: No negative markers found for", celltype, "in the dataset\n")
        next
      }
      
      if (length(available_neg_markers) < length(neg_markers[[celltype]])) {
        missing_neg_markers <- setdiff(neg_markers[[celltype]], available_neg_markers)
        cat(paste("Warning: Missing negative markers for", celltype, ":", 
                  paste(missing_neg_markers, collapse = ", "), "\n"))
      }
      
      cat(paste("Using", length(available_neg_markers), "negative markers for", celltype, "\n"))
      
      # Add negative module score
      seurat_obj <- AddModuleScore(
        object = seurat_obj,
        features = list(available_neg_markers),
        name = paste0(celltype, "_neg_score"),
        ctrl = ctrl_genes
      )
      
      neg_score_cols <- c(neg_score_cols, paste0(celltype, "_neg_score1"))
    }
  }
  
  # Filter available score columns
  pos_score_cols <- pos_score_cols[pos_score_cols %in% colnames(seurat_obj@meta.data)]
  neg_score_cols <- neg_score_cols[neg_score_cols %in% colnames(seurat_obj@meta.data)]
  
  if (length(pos_score_cols) == 0) {
    stop("No positive module scores were calculated successfully")
  }
  
  # Create combined scores
  cat("\nCombining positive and negative scores...\n")
  
  # Get positive scores
  pos_score_matrix <- seurat_obj@meta.data[, pos_score_cols, drop = FALSE]
  colnames(pos_score_matrix) <- gsub("_pos_score1", "", colnames(pos_score_matrix))
  
  # Get negative scores (if available)
  if (has_neg_markers && length(neg_score_cols) > 0) {
    neg_score_matrix <- seurat_obj@meta.data[, neg_score_cols, drop = FALSE]
    colnames(neg_score_matrix) <- gsub("_neg_score1", "", colnames(neg_score_matrix))
    
    # Ensure same cell types in both matrices
    common_celltypes <- intersect(colnames(pos_score_matrix), colnames(neg_score_matrix))
    
    if (length(common_celltypes) > 0) {
      # Calculate combined scores: positive - (weight * negative)
      combined_score_matrix <- pos_score_matrix[, common_celltypes, drop = FALSE] - 
        (neg_weight * neg_score_matrix[, common_celltypes, drop = FALSE])
      
      # Add cell types that only have positive scores
      pos_only_celltypes <- setdiff(colnames(pos_score_matrix), common_celltypes)
      if (length(pos_only_celltypes) > 0) {
        combined_score_matrix <- cbind(combined_score_matrix, 
                                       pos_score_matrix[, pos_only_celltypes, drop = FALSE])
      }
      
      cat("Combined scores calculated for", length(common_celltypes), "cell types with both positive and negative markers\n")
      cat("Positive-only scores for", length(pos_only_celltypes), "cell types\n")
      
    } else {
      cat("Warning: No cell types have both positive and negative markers. Using positive scores only.\n")
      combined_score_matrix <- pos_score_matrix
    }
  } else {
    cat("Using positive scores only (no negative markers available)\n")
    combined_score_matrix <- pos_score_matrix
  }
  
  # Assign cell types based on highest combined score
  cat("\nAssigning cell types based on highest combined scores...\n")
  
  # Find the cell type with maximum score for each cell
  max_scores <- apply(combined_score_matrix, 1, max)
  second_max_scores <- apply(combined_score_matrix, 1, function(x) {
    if(length(x) > 1) {
      return(sort(x, decreasing = TRUE)[2])
    } else {
      return(0)
    }
  })
  score_differences <- max_scores - second_max_scores
  
  # Assign cell types
  assigned_celltypes <- colnames(combined_score_matrix)[apply(combined_score_matrix, 1, which.max)]
  
  # Set cells with low confidence to "Unknown"
  low_confidence <- score_differences < min_score_diff
  assigned_celltypes[low_confidence] <- "Unknown"
  
  # Add annotations to metadata
  seurat_obj@meta.data$celltype_annotation_dual <- assigned_celltypes
  seurat_obj@meta.data$max_score_dual <- max_scores
  seurat_obj@meta.data$score_difference_dual <- score_differences
  
  # Print summary
  cat("\nDual marker cell type annotation summary:\n")
  annotation_table <- table(assigned_celltypes)
  print(annotation_table)
  
  cat(paste("\nCells with low confidence (score diff <", min_score_diff, "):", sum(low_confidence), "\n"))
  cat(paste("Annotation rate:", round(100 * sum(assigned_celltypes != "Unknown") / length(assigned_celltypes), 2), "%\n"))
  
  # Create detailed summary
  summary_df <- data.frame(
    CellType = names(annotation_table),
    Count = as.numeric(annotation_table),
    Percentage = round(as.numeric(annotation_table) / length(assigned_celltypes) * 100, 2)
  )
  
  cat("\nDetailed dual marker annotation summary:\n")
  print(summary_df)
  
  return(list(
    seurat_object = seurat_obj,
    positive_markers = pos_markers,
    negative_markers = neg_markers,
    annotation_summary = summary_df,
    positive_scores = pos_score_matrix,
    negative_scores = if(has_neg_markers && length(neg_score_cols) > 0) neg_score_matrix else NULL,
    combined_scores = combined_score_matrix
  ))
}

# Function to parse marker gene dataframe and create annotation list
parse_marker_dataframe <- function(marker_df, 
                                   celltype_col = "Celltype", 
                                   marker_col = "Pos_Marker_Genes") {
  
  cat("Parsing marker gene dataframe...\n")
  
  # Check if required columns exist
  if (!celltype_col %in% colnames(marker_df)) {
    stop(paste("Column", celltype_col, "not found in dataframe"))
  }
  
  if (!marker_col %in% colnames(marker_df)) {
    stop(paste("Column", marker_col, "not found in dataframe"))
  }
  
  # Create marker gene list
  marker_genes <- list()
  
  for (i in 1:nrow(marker_df)) {
    celltype <- marker_df[[celltype_col]][i]
    marker_string <- marker_df[[marker_col]][i]
    
    # Parse comma-separated genes and clean whitespace
    genes <- str_split(marker_string, ",")[[1]]
    genes <- str_trim(genes)  # Remove leading/trailing whitespace
    genes <- genes[genes != ""]  # Remove empty strings
    
    # Store in list
    marker_genes[[celltype]] <- genes
    
    cat("Cell type", celltype, ":", length(genes), "markers\n")
  }
  
  return(marker_genes)
}

# Enhanced annotation function that works with dataframe input
annotate_from_dataframe <- function(seurat_obj, 
                                    marker_df,
                                    celltype_col = "Celltype",
                                    marker_col = "Pos_Marker_Genes",
                                    min_score_diff = 0.1,
                                    ctrl_genes = 100) {
  
  cat("=== Annotating Seurat object using marker dataframe ===\n")
  
  # Parse the dataframe to create marker gene list
  marker_genes <- parse_marker_dataframe(marker_df, celltype_col, marker_col)
  
  # Print summary of marker genes
  cat("\nMarker gene summary:\n")
  for (celltype in names(marker_genes)) {
    cat(paste0(celltype, ": ", paste(marker_genes[[celltype]], collapse = ", "), "\n"))
  }
  
  # Calculate module scores for each cell type
  cat("\nCalculating module scores for each cell type...\n")
  
  for (celltype in names(marker_genes)) {
    # Check which markers are present in the dataset
    available_markers <- intersect(marker_genes[[celltype]], rownames(seurat_obj))
    
    if (length(available_markers) == 0) {
      warning(paste("No markers found for", celltype, "in the dataset"))
      next
    }
    
    if (length(available_markers) < length(marker_genes[[celltype]])) {
      missing_markers <- setdiff(marker_genes[[celltype]], available_markers)
      cat(paste("Warning: Missing markers for", celltype, ":", paste(missing_markers, collapse = ", "), "\n"))
    }
    
    cat(paste("Using", length(available_markers), "markers for", celltype, "\n"))
    
    # Add module score
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(available_markers),
      name = paste0(celltype, "_score"),
      ctrl = ctrl_genes
    )
  }
  
  # Get all score column names
  score_cols <- paste0(names(marker_genes), "_score1")
  score_cols <- score_cols[score_cols %in% colnames(seurat_obj@meta.data)]
  
  if (length(score_cols) == 0) {
    stop("No module scores were calculated successfully")
  }
  
  # Create a matrix of scores
  score_matrix <- seurat_obj@meta.data[, score_cols, drop = FALSE]
  colnames(score_matrix) <- gsub("_score1", "", colnames(score_matrix))
  
  # Assign cell types based on highest score
  cat("\nAssigning cell types based on highest scores...\n")
  
  # Find the cell type with maximum score for each cell
  max_scores <- apply(score_matrix, 1, max)
  second_max_scores <- apply(score_matrix, 1, function(x) {
    if(length(x) > 1) {
      return(sort(x, decreasing = TRUE)[2])
    } else {
      return(0)
    }
  })
  score_differences <- max_scores - second_max_scores
  
  # Assign cell types
  assigned_celltypes <- colnames(score_matrix)[apply(score_matrix, 1, which.max)]
  
  # Set cells with low confidence to "Unknown"
  low_confidence <- score_differences < min_score_diff
  assigned_celltypes[low_confidence] <- "Unknown"
  
  # Add annotations to metadata
  seurat_obj@meta.data$celltype_annotation <- assigned_celltypes
  seurat_obj@meta.data$max_score <- max_scores
  seurat_obj@meta.data$score_difference <- score_differences
  
  # Print summary
  cat("\nCell type annotation summary:\n")
  annotation_table <- table(assigned_celltypes)
  print(annotation_table)
  
  cat(paste("\nCells with low confidence (score diff <", min_score_diff, "):", sum(low_confidence), "\n"))
  cat(paste("Annotation rate:", round(100 * sum(assigned_celltypes != "Unknown") / length(assigned_celltypes), 2), "%\n"))
  
  # Create detailed summary
  summary_df <- data.frame(
    CellType = names(annotation_table),
    Count = as.numeric(annotation_table),
    Percentage = round(as.numeric(annotation_table) / length(assigned_celltypes) * 100, 2)
  )
  
  cat("\nDetailed annotation summary:\n")
  print(summary_df)
  
  return(list(
    seurat_object = seurat_obj,
    marker_genes = marker_genes,
    annotation_summary = summary_df,
    score_matrix = score_matrix
  ))
}

# Function to validate marker genes against Seurat object
validate_markers <- function(seurat_obj, marker_df, 
                             celltype_col = "Celltype", 
                             marker_col = "Pos_Marker_Genes") {
  
  cat("=== Validating marker genes against Seurat object ===\n")
  
  # Parse markers
  marker_genes <- parse_marker_dataframe(marker_df, celltype_col, marker_col)
  
  # Check availability
  validation_results <- data.frame(
    CellType = character(),
    Total_Markers = integer(),
    Available_Markers = integer(),
    Missing_Markers = integer(),
    Availability_Rate = numeric(),
    Available_Genes = character(),
    Missing_Genes = character(),
    stringsAsFactors = FALSE
  )
  
  for (celltype in names(marker_genes)) {
    markers <- marker_genes[[celltype]]
    available <- intersect(markers, rownames(seurat_obj))
    missing <- setdiff(markers, rownames(seurat_obj))
    
    validation_results <- rbind(validation_results, data.frame(
      CellType = celltype,
      Total_Markers = length(markers),
      Available_Markers = length(available),
      Missing_Markers = length(missing),
      Availability_Rate = round(length(available) / length(markers) * 100, 1),
      Available_Genes = paste(available, collapse = ", "),
      Missing_Genes = paste(missing, collapse = ", "),
      stringsAsFactors = FALSE
    ))
  }
  
  cat("\nMarker validation results:\n")
  print(validation_results[, 1:5])  # Print summary columns
  
  # Identify problematic cell types
  low_availability <- validation_results$Availability_Rate < 50
  if (any(low_availability)) {
    cat("\nWARNING: Cell types with <50% marker availability:\n")
    print(validation_results[low_availability, c("CellType", "Availability_Rate", "Missing_Genes")])
  }
  
  return(validation_results)
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)

# Function to debug celltype column issues
debug_celltype_column <- function(seurat_obj, celltype_column) {
  cat("=== Debugging celltype column ===\n")
  
  # Check if column exists
  if (!celltype_column %in% colnames(seurat_obj@meta.data)) {
    cat("ERROR: Column", celltype_column, "not found!\n")
    cat("Available columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
    return(NULL)
  }
  
  # Check data type and structure
  col_data <- seurat_obj@meta.data[[celltype_column]]
  cat("Column exists:", celltype_column, "\n")
  cat("Data type:", class(col_data), "\n")
  cat("Length:", length(col_data), "\n")
  
  # Check first few values
  if (is.list(col_data)) {
    cat("First 5 values (list):", str(head(col_data, 5)), "\n")
  } else {
    cat("First 5 values:", head(col_data, 5), "\n")
  }
  
  # Check for unique values
  if (is.list(col_data)) {
    unique_vals <- unique(unlist(col_data))
  } else {
    unique_vals <- unique(col_data)
  }
  
  cat("Number of unique values:", length(unique_vals), "\n")
  cat("Unique values:", paste(head(unique_vals, 10), collapse = ", "), "\n")
  
  return(col_data)
}

# Function to create comprehensive UMAP and PCA visualizations
visualize_celltype_umap_pca <- function(seurat_obj, 
                                        celltype_column = "celltype_annotation",
                                        sample_column = NULL,
                                        save_plots = FALSE,
                                        output_prefix = "celltype_viz") {
  
  cat("Creating UMAP and PCA visualizations for cell type annotations...\n")
  
  # Check if celltype column exists
  if (!celltype_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", celltype_column, "not found in metadata. Available columns:", 
               paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }
  
  # Check if required reductions exist
  if (!"umap" %in% names(seurat_obj@reductions)) {
    cat("UMAP not found. Running UMAP...\n")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  }
  
  if (!"pca" %in% names(seurat_obj@reductions)) {
    cat("PCA not found. Running PCA...\n")
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
  }
  
  # Get cell types and create color palette
  # Handle case where celltype_column might be a list or factor
  celltype_data <- seurat_obj@meta.data[[celltype_column]]
  
  # Convert to character if it's a factor or list
  if (is.list(celltype_data)) {
    celltype_data <- unlist(celltype_data)
  }
  if (is.factor(celltype_data)) {
    celltype_data <- as.character(celltype_data)
  }
  
  # Update the metadata with clean character vector
  seurat_obj@meta.data[[celltype_column]] <- celltype_data
  
  cell_types <- unique(celltype_data)
  cell_types <- cell_types[!is.na(cell_types)]
  cell_types <- sort(cell_types)  # Sort for consistent colors
  n_types <- length(cell_types)
  
  cat("Found", n_types, "cell types:", paste(cell_types, collapse = ", "), "\n")
  
  # Create a comprehensive color palette
  if (n_types <= 8) {
    colors <- brewer.pal(max(3, n_types), "Set2")
  } else if (n_types <= 12) {
    colors <- brewer.pal(n_types, "Set3")
  } else {
    # For more than 12 cell types, use a combination of palettes
    colors <- c(brewer.pal(12, "Set3"), 
                brewer.pal(min(n_types-12, 8), "Dark2"))
    if (n_types > 20) {
      colors <- c(colors, rainbow(n_types - 20))
    }
  }
  colors <- colors[1:n_types]
  names(colors) <- cell_types
  
  plots <- list()
  
  # 1. UMAP Plots
  cat("Creating UMAP visualizations...\n")
  
  # Main UMAP with labels
  plots$umap_labeled <- DimPlot(seurat_obj, 
                                reduction = "umap",
                                group.by = celltype_column,
                                label = TRUE,
                                label.size = 3,
                                cols = colors,
                                pt.size = 0.5) +
    ggtitle("UMAP: Cell Type Annotations") +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
  
  # Clean UMAP without labels
  plots$umap_clean <- DimPlot(seurat_obj, 
                              reduction = "umap",
                              group.by = celltype_column,
                              label = FALSE,
                              cols = colors,
                              pt.size = 0.3) +
    ggtitle("UMAP: Cell Type Annotations (Clean)") +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
  
  # 2. PCA Plots
  cat("Creating PCA visualizations...\n")
  
  # PCA with labels
  plots$pca_labeled <- DimPlot(seurat_obj, 
                               reduction = "pca",
                               group.by = celltype_column,
                               label = TRUE,
                               label.size = 3,
                               cols = colors,
                               pt.size = 0.5) +
    ggtitle("PCA: Cell Type Annotations") +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
  
  # Clean PCA without labels
  plots$pca_clean <- DimPlot(seurat_obj, 
                             reduction = "pca",
                             group.by = celltype_column,
                             label = FALSE,
                             cols = colors,
                             pt.size = 0.3) +
    ggtitle("PCA: Cell Type Annotations (Clean)") +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
  
  # 3. Split by cell type (UMAP and PCA)
  plots$umap_split <- DimPlot(seurat_obj, 
                              reduction = "umap",
                              group.by = celltype_column,
                              split.by = celltype_column,
                              ncol = 4,
                              cols = colors,
                              pt.size = 0.1) +
    ggtitle("UMAP: Cell Types (Split View)") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))
  
  plots$pca_split <- DimPlot(seurat_obj, 
                             reduction = "pca",
                             group.by = celltype_column,
                             split.by = celltype_column,
                             ncol = 4,
                             cols = colors,
                             pt.size = 0.1) +
    ggtitle("PCA: Cell Types (Split View)") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))
  
  # 4. Cell type proportions bar plot
  # Create a safe dataframe for counting
  count_df <- data.frame(
    celltype = seurat_obj@meta.data[[celltype_column]],
    stringsAsFactors = FALSE
  )
  
  prop_data <- count_df %>%
    count(celltype) %>%
    mutate(proportion = n / sum(n) * 100) %>%
    arrange(desc(n))
  
  # Rename column to match original celltype_column name
  colnames(prop_data)[1] <- celltype_column
  
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
  
  # 5. Sample composition (if sample column provided)
  if (!is.null(sample_column) && sample_column %in% colnames(seurat_obj@meta.data)) {
    
    # Create safe dataframe for composition analysis
    comp_df <- data.frame(
      sample = seurat_obj@meta.data[[sample_column]],
      celltype = seurat_obj@meta.data[[celltype_column]],
      stringsAsFactors = FALSE
    )
    
    composition_data <- comp_df %>%
      group_by(sample, celltype) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(sample) %>%
      mutate(proportion = count / sum(count) * 100)
    
    # Rename columns to match original names
    colnames(composition_data)[1] <- sample_column
    colnames(composition_data)[2] <- celltype_column
    
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
    
    # UMAP colored by sample
    plots$umap_by_sample <- DimPlot(seurat_obj, 
                                    reduction = "umap",
                                    group.by = sample_column,
                                    pt.size = 0.3) +
      ggtitle("UMAP: Colored by Sample") +
      theme_minimal()
  }
  
  # 6. Confidence plots (if confidence metrics available)
  if ("score_difference" %in% colnames(seurat_obj@meta.data)) {
    
    # Create safe dataframe for confidence analysis
    conf_df <- data.frame(
      celltype = seurat_obj@meta.data[[celltype_column]],
      score_difference = seurat_obj@meta.data$score_difference,
      stringsAsFactors = FALSE
    )
    
    # Confidence by cell type
    plots$confidence_boxplot <- ggplot(conf_df, 
                                       aes(x = celltype, 
                                           y = score_difference,
                                           fill = celltype)) +
      geom_boxplot(alpha = 0.7) +
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "Annotation Confidence by Cell Type",
           subtitle = "Red line shows confidence threshold (0.1)",
           x = "Cell Type",
           y = "Score Difference")
    
    # UMAP colored by confidence
    plots$umap_confidence <- FeaturePlot(seurat_obj, 
                                         features = "score_difference",
                                         reduction = "umap",
                                         pt.size = 0.3) +
      ggtitle("UMAP: Annotation Confidence") +
      theme_minimal()
  }
  
  # 7. Combined layout plots
  plots$combined_main <- plots$umap_labeled + plots$pca_labeled + 
    plot_layout(ncol = 2)
  
  plots$combined_clean <- plots$umap_clean + plots$pca_clean + 
    plot_layout(ncol = 2)
  
  # Display key plots
  cat("Displaying main visualizations...\n")
  print(plots$combined_main)
  print(plots$proportions)
  
  if (!is.null(plots$confidence_boxplot)) {
    print(plots$confidence_boxplot)
  }
  
  # Save plots if requested
  if (save_plots) {
    cat("Saving plots...\n")
    
    ggsave(paste0(output_prefix, "_umap_labeled.pdf"), 
           plots$umap_labeled, width = 12, height = 8)
    ggsave(paste0(output_prefix, "_pca_labeled.pdf"), 
           plots$pca_labeled, width = 12, height = 8)
    ggsave(paste0(output_prefix, "_umap_clean.pdf"), 
           plots$umap_clean, width = 10, height = 8)
    ggsave(paste0(output_prefix, "_pca_clean.pdf"), 
           plots$pca_clean, width = 10, height = 8)
    ggsave(paste0(output_prefix, "_combined_main.pdf"), 
           plots$combined_main, width = 20, height = 8)
    ggsave(paste0(output_prefix, "_proportions.pdf"), 
           plots$proportions, width = 10, height = 8)
    ggsave(paste0(output_prefix, "_umap_split.pdf"), 
           plots$umap_split, width = 16, height = 12)
    ggsave(paste0(output_prefix, "_pca_split.pdf"), 
           plots$pca_split, width = 16, height = 12)
    
    if (!is.null(plots$confidence_boxplot)) {
      ggsave(paste0(output_prefix, "_confidence.pdf"), 
             plots$confidence_boxplot, width = 12, height = 8)
    }
  }
  
  cat("Visualization complete!\n")
  
  return(plots)
}

# Function to create marker gene expression overlays
plot_marker_expression_overlay <- function(seurat_obj, 
                                           marker_genes,
                                           celltype_column = "celltype_annotation",
                                           max_markers = 8) {
  
  cat("Creating marker gene expression overlays...\n")
  
  # Flatten marker genes list and get available markers
  all_markers <- unique(unlist(marker_genes))
  available_markers <- intersect(all_markers, rownames(seurat_obj))
  
  if (length(available_markers) == 0) {
    cat("No marker genes found in the dataset!\n")
    return(NULL)
  }
  
  # Select top markers (limit for visualization)
  selected_markers <- head(available_markers, max_markers)
  
  plots <- list()
  
  # UMAP feature plots
  plots$umap_features <- FeaturePlot(seurat_obj, 
                                     features = selected_markers,
                                     reduction = "umap",
                                     ncol = 4,
                                     pt.size = 0.3) &
    theme_minimal() &
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # PCA feature plots
  plots$pca_features <- FeaturePlot(seurat_obj, 
                                    features = selected_markers,
                                    reduction = "pca",
                                    ncol = 4,
                                    pt.size = 0.3) &
    theme_minimal() &
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # Violin plots
  plots$violin_plots <- VlnPlot(seurat_obj, 
                                features = selected_markers,
                                group.by = celltype_column,
                                ncol = 4,
                                pt.size = 0) &
    theme_minimal() &
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plots)
}

# Main execution function
create_umap_pca_visualizations <- function(seurat_obj, 
                                           marker_genes = NULL,
                                           celltype_column = "celltype_annotation",
                                           sample_column = NULL,
                                           save_plots = FALSE,
                                           output_prefix = "celltype_umap_pca") {
  
  cat("=== Creating comprehensive UMAP and PCA visualizations ===\n")
  
  # Main visualizations
  main_plots <- visualize_celltype_umap_pca(seurat_obj,
                                            celltype_column = celltype_column,
                                            sample_column = sample_column,
                                            save_plots = save_plots,
                                            output_prefix = output_prefix)
  
  # Marker gene overlays
  marker_plots <- NULL
  if (!is.null(marker_genes)) {
    marker_plots <- plot_marker_expression_overlay(seurat_obj,
                                                   marker_genes,
                                                   celltype_column = celltype_column)
  }
  
  return(list(
    main_plots = main_plots,
    marker_plots = marker_plots
  ))
}

get_top_markers <- function(markers_df, top_n = 10) {
  
  top_markers <- markers_df %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
  
  return(top_markers)
}


create_marker_dotplot <- function(seurat_obj, markers_df, top_n = 5, 
                                  title = "Top Marker Genes by Cluster") {
  
  # Get top markers
  top_markers <- get_top_markers(markers_df, top_n = top_n)
  
  # Get unique marker genes
  marker_genes <- unique(top_markers$gene)
  
  # Create dot plot
  p <- DotPlot(seurat_obj, 
               features = marker_genes,
               group.by = "seurat_clusters") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    ) +
    labs(
      title = title,
      x = "Marker Genes",
      y = "Clusters"
    ) +
    scale_color_gradient(low = "lightgrey", high = "red", name = "Avg Expression") +
    guides(size = guide_legend(title = "Pct Expressed"))
  
  return(p)
}
create_marker_heatmap <- function(seurat_obj, markers_df, top_n = 5) {
  
  # Get top markers
  top_markers <- get_top_markers(markers_df, top_n = top_n)
  
  # Get unique marker genes
  marker_genes <- unique(top_markers$gene)
  
  # Create heatmap
  p <- DoHeatmap(seurat_obj, 
                 features = marker_genes,
                 group.by = "seurat_clusters",
                 size = 3,
                 angle = 45) +
    theme(axis.text.y = element_text(size = 8)) +
    ggtitle("Top Marker Genes Heatmap")
  
  return(p)
}

fix_scaling_and_create_heatmap <- function(seurat_obj, markers_df, top_n = 20) {
  
  cat("Fixing scaling issues and creating heatmap...\n")
  
  # Get top markers
  top_markers <- get_top_markers(markers_df, top_n = top_n)
  marker_genes <- unique(top_markers$gene)
  
  # Check which genes are available
  available_genes <- intersect(marker_genes, rownames(seurat_obj))
  missing_genes <- setdiff(marker_genes, rownames(seurat_obj))
  
  cat("Available marker genes:", length(available_genes), "\n")
  if (length(missing_genes) > 0) {
    cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  }
  
  # Scale all genes in the dataset to avoid this issue
  cat("Scaling all genes to ensure marker genes are included...\n")
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE)
  
  # Now create heatmap
  p <- DoHeatmap(seurat_obj, 
                 features = available_genes,
                 group.by = "seurat_clusters",
                 size = 3,
                 angle = 45) +
    theme(axis.text.y = element_text(size = 8)) +
    ggtitle(paste("Top", top_n, "Marker Genes Heatmap (Fixed)"))
  
  return(list(seurat_object = seurat_obj, heatmap = p))
}

create_marker_summary_table <- function(markers_df, top_n = 5) {
  
  summary_table <- markers_df %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
    mutate(
      avg_log2FC = round(avg_log2FC, 3),
      pct.1 = round(pct.1, 3),
      pct.2 = round(pct.2, 3),
      p_val_adj = format(p_val_adj, scientific = TRUE, digits = 3)
    )
  
  return(summary_table)
}

#Split marker plot
create_marker_dotplot_split <- function(seurat_obj, markers_df, top_n = 5, 
                                        title = "Top Marker Genes by Cluster",
                                        clusters_per_plot = 5) {
  
  # Get top markers
  top_markers <- get_top_markers(markers_df, top_n = top_n)
  
  # Get unique marker genes
  marker_genes <- unique(top_markers$gene)
  
  # Get all cluster levels
  all_clusters <- levels(seurat_obj$seurat_clusters)
  if (is.null(all_clusters)) {
    all_clusters <- sort(unique(seurat_obj$seurat_clusters))
  }
  
  # Split clusters into groups
  n_clusters <- length(all_clusters)
  cluster_groups <- split(all_clusters, ceiling(seq_along(all_clusters) / clusters_per_plot))
  
  # Create a plot for each group
  plot_list <- lapply(seq_along(cluster_groups), function(i) {
    
    # Subset to cells in these clusters
    cells_subset <- WhichCells(seurat_obj, 
                               expression = seurat_clusters %in% cluster_groups[[i]])
    so_subset <- subset(seurat_obj, cells = cells_subset)
    
    # Create dot plot
    p <- DotPlot(so_subset, 
                 features = marker_genes,
                 group.by = "seurat_clusters") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      labs(
        title = paste0(title, " (Part ", i, " of ", length(cluster_groups), ")"),
        x = "Marker Genes",
        y = "Clusters"
      ) +
      scale_color_gradient(low = "lightgrey", high = "red", name = "Avg Expression") +
      guides(size = guide_legend(title = "Pct Expressed"))
    
    return(p)
  })
  
  names(plot_list) <- paste0("plot_", seq_along(plot_list))
  return(plot_list)
}
