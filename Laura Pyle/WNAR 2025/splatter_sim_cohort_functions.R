library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)

########################################################
# FUNCTION TO ADD CORRELATION TO CELLS WITHIN A PERSON #
########################################################

add_cell_state_correlation <- function(sim, 
                                       n_states = 3,
                                       state_correlation = 0.5) {
  
  counts_original <- counts(sim)
  
  for(batch in unique(colData(sim)$Batch)) {
    batch_cells <- which(colData(sim)$Batch == batch)
    n_cells <- length(batch_cells)
    
    # Assign cells to states (some overlap creates correlation)
    if(state_correlation > 0) {
      # Soft clustering - cells can partially belong to multiple states
      state_weights <- matrix(runif(n_cells * n_states), 
                              nrow = n_cells, ncol = n_states)
      
      # Add correlation by making nearby cells similar
      for(i in 2:n_cells) {
        state_weights[i, ] <- (1 - state_correlation) * state_weights[i, ] + 
          state_correlation * state_weights[i-1, ]
      }
      
      # Normalize
      state_weights <- state_weights / rowSums(state_weights)
      
      # Create state-specific expression programs
      state_programs <- matrix(rnorm(nrow(counts_original) * n_states, 0, 0.2),
                               nrow = nrow(counts_original))
      
      # Apply weighted combination of programs
      for(i in 1:n_cells) {
        cell_program <- state_programs %*% state_weights[i, ]
        counts_original[, batch_cells[i]] <- 
          counts_original[, batch_cells[i]] * exp(cell_program)
      }
    }
  }
  
  assay(sim, "counts") <- round(counts_original)
  return(sim)
}

#sim_states <- add_cell_state_correlation(sim, 
#                                         n_states = 4,
 #                                        state_correlation = 0.6)

##################################
# BASE FUNCTION TO DO SIMULATION #
##################################

simulate_with_correlation <- function(n_individuals = 8,
                                          cells_per_ind = 150,
                                          cells_sd = NULL,
                                          cells_min = 30,
                                          n_genes = 3000,
                                          batch_facLoc = 0.1,
                                          batch_facScale = 0.25,
                                          within_ind_correlation = 0.4,
                                          n_cell_states = 5,
                                          verbose = TRUE) {
  
  # Handle cell count variability
  if(length(cells_per_ind) == 1) {
    # Single value provided - generate variable counts
    if(is.null(cells_sd)) {
      # Default to 30% coefficient of variation if sd not specified
      cells_sd <- cells_per_ind * 0.3
    }
    
    # Generate variable cell counts
    set.seed(123)  # For reproducibility
    cell_counts <- round(rnorm(n_individuals, 
                               mean = cells_per_ind, 
                               sd = cells_sd))
    
    # Apply minimum threshold
    cell_counts[cell_counts < cells_min] <- cells_min
    
  } else if(length(cells_per_ind) == n_individuals) {
    # Vector of counts provided - use directly
    cell_counts <- cells_per_ind
    
    # Still apply minimum threshold
    cell_counts[cell_counts < cells_min] <- cells_min
    
  } else {
    stop("cells_per_ind must be either a single value or a vector of length n_individuals")
  }
  
  # Print summary if verbose
  if(verbose) {
    cat("=== Simulation Parameters ===\n")
    cat("Genes:", n_genes, "\n")
    cat("Individuals:", n_individuals, "\n")
    cat("Batch effect location:", batch_facLoc, "\n")
    cat("Batch effect scale:", batch_facScale, "\n")
    cat("Within-individual correlation:", within_ind_correlation, "\n")
    cat("Cell states:", n_cell_states, "\n\n")
    
    cat("Cell distribution across individuals:\n")
    cat("  Mean:", round(mean(cell_counts), 1), "\n")
    cat("  SD:", round(sd(cell_counts), 1), "\n")
    cat("  Range:", min(cell_counts), "-", max(cell_counts), "\n")
    cat("  Total cells:", sum(cell_counts), "\n\n")
  }
  
  # Base simulation with variable parameters
  sim <- splatSimulate(
    nGenes = n_genes,           # Now customizable
    batchCells = cell_counts,
    batch.facLoc = batch_facLoc,    # Now customizable
    batch.facScale = batch_facScale, # Now customizable
    verbose = FALSE
  )
  
  # Add correlation through shared cell states
  sim <- add_cell_state_correlation(sim, 
                                    n_states = n_cell_states,
                                    state_correlation = within_ind_correlation)
  
  # Add metadata with actual cell counts
  colData(sim)$Individual <- factor(colData(sim)$Batch,
                                    labels = paste0("Patient_", 1:n_individuals))
  
  # Add cell count info to metadata
  ind_cell_counts <- table(colData(sim)$Individual)
  colData(sim)$CellsPerIndividual <- ind_cell_counts[colData(sim)$Individual]
  
  # Add simulation parameters as metadata
  metadata(sim) <- list(
    n_genes = n_genes,
    n_individuals = n_individuals,
    batch_facLoc = batch_facLoc,
    batch_facScale = batch_facScale,
    within_ind_correlation = within_ind_correlation,
    n_cell_states = n_cell_states,
    total_cells = sum(cell_counts)
  )
  
  return(sim)
}

###############################
# FUNCTION TO PLOT EXPRESSION #
###############################

#' Plot violin plots of average gene expression between two groups
#'
#' @param sce SingleCellExperiment object
#' @param group_var Character string specifying the column name in colData containing group information
#' @param assay_name Character string specifying which assay to use (default: "logcounts")
#' @param groups Optional character vector of length 2 specifying which groups to compare.
#'               If NULL, will use the first two unique values in group_var
#' @param plot_title Character string for the plot title
#' @param y_label Character string for y-axis label
#' @param colors Character vector of length 2 specifying colors for each group
#' @param add_stats Logical, whether to add statistical comparison (default: TRUE)
#' @param return_data Logical, whether to return the data frame used for plotting (default: FALSE)
#'
#' @return A ggplot object (or list with plot and data if return_data = TRUE)
#' 
#' @examples
#' # Basic usage
#' plot_expression_violin(sce, group_var = "cell_type")
#' 
#' # Specify groups and colors
#' plot_expression_violin(sce, 
#'                       group_var = "condition",
#'                       groups = c("Control", "Treatment"),
#'                       colors = c("#4DAF4A", "#E41A1C"))

plot_expression_violin <- function(sce,
                                   group_var,
                                   assay_name = NULL,
                                   groups = NULL,
                                   plot_title = "Average Gene Expression by Group",
                                   y_label = "Average Expression",
                                   colors = c("#377EB8", "#E41A1C"),
                                   add_stats = TRUE,
                                   return_data = FALSE) {
  
  # Input validation
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object")
  }
  
  if (!group_var %in% colnames(colData(sce))) {
    stop(paste("group_var '", group_var, "' not found in colData", sep = ""))
  }
  
  # Auto-select assay if not specified
  if (is.null(assay_name)) {
    available_assays <- assayNames(sce)
    # Priority order for common assay names
    preferred_assays <- c("logcounts", "normcounts", "counts", "CellMeans", "TrueCounts")
    assay_name <- preferred_assays[preferred_assays %in% available_assays][1]
    
    # If no preferred assay found, use the first available
    if (is.na(assay_name)) {
      assay_name <- available_assays[1]
    }
    message(paste("Using assay:", assay_name))
  }
  
  if (!assay_name %in% assayNames(sce)) {
    stop(paste("assay_name '", assay_name, "' not found in object. Available assays: ",
               paste(assayNames(sce), collapse = ", "), sep = ""))
  }
  
  # Extract expression matrix
  expr_matrix <- assay(sce, assay_name)
  
  # Calculate average expression per cell (mean across all genes)
  avg_expression <- colMeans(expr_matrix, na.rm = TRUE)
  
  # Create data frame for plotting
  plot_df <- data.frame(
    cell_id = colnames(sce),
    avg_expression = avg_expression,
    group = colData(sce)[[group_var]],
    stringsAsFactors = FALSE
  )
  
  # Filter to specified groups if provided
  if (!is.null(groups)) {
    if (length(groups) != 2) {
      stop("groups must be a vector of length 2")
    }
    if (!all(groups %in% plot_df$group)) {
      missing_groups <- groups[!groups %in% plot_df$group]
      stop(paste("Groups not found in data:", paste(missing_groups, collapse = ", ")))
    }
    plot_df <- plot_df[plot_df$group %in% groups, ]
    plot_df$group <- factor(plot_df$group, levels = groups)
  } else {
    # Use first two groups if not specified
    unique_groups <- unique(plot_df$group)
    if (length(unique_groups) < 2) {
      stop("Need at least 2 groups for comparison")
    }
    if (length(unique_groups) > 2) {
      warning(paste("More than 2 groups found. Using first two:", 
                    unique_groups[1], "and", unique_groups[2]))
      plot_df <- plot_df[plot_df$group %in% unique_groups[1:2], ]
    }
    plot_df$group <- factor(plot_df$group)
  }
  
  # Create violin plot
  p <- ggplot(plot_df, aes(x = group, y = avg_expression, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = colors) +
    labs(title = plot_title,
         x = "Group",
         y = y_label) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none"
    )
  
  # Add statistical test if requested
  if (add_stats) {
    # Perform Wilcoxon rank sum test
    group_levels <- levels(plot_df$group)
    group1_data <- plot_df$avg_expression[plot_df$group == group_levels[1]]
    group2_data <- plot_df$avg_expression[plot_df$group == group_levels[2]]
    
    test_result <- wilcox.test(group1_data, group2_data)
    p_value <- test_result$p.value
    
    # Format p-value for display
    if (p_value < 0.001) {
      p_label <- "p < 0.001"
    } else if (p_value < 0.01) {
      p_label <- sprintf("p = %.3f", p_value)
    } else {
      p_label <- sprintf("p = %.2f", p_value)
    }
    
    # Add p-value to plot
    y_max <- max(plot_df$avg_expression, na.rm = TRUE)
    y_range <- diff(range(plot_df$avg_expression, na.rm = TRUE))
    
    p <- p + 
      annotate("segment", x = 1, xend = 2, 
               y = y_max + y_range * 0.05, yend = y_max + y_range * 0.05) +
      annotate("text", x = 1.5, y = y_max + y_range * 0.08,
               label = p_label, size = 3.5)
  }
  
  # Add sample sizes to x-axis labels
  sample_sizes <- table(plot_df$group)
  x_labels <- paste0(names(sample_sizes), "\n(n=", sample_sizes, ")")
  p <- p + scale_x_discrete(labels = x_labels)
  
  # Return plot or plot with data
  if (return_data) {
    return(list(plot = p, data = plot_df))
  } else {
    return(p)
  }
}

# Alternative function for comparing expression of specific genes
plot_gene_expression_violin <- function(sce,
                                        genes,
                                        group_var,
                                        assay_name = NULL,
                                        groups = NULL,
                                        plot_title = NULL,
                                        y_label = "Expression",
                                        colors = c("#377EB8", "#E41A1C"),
                                        aggregate_fun = mean) {
  
  # Input validation
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object")
  }
  
  # Auto-select assay if not specified
  if (is.null(assay_name)) {
    available_assays <- assayNames(sce)
    # Priority order for common assay names
    preferred_assays <- c("logcounts", "normcounts", "counts", "CellMeans", "TrueCounts")
    assay_name <- preferred_assays[preferred_assays %in% available_assays][1]
    
    # If no preferred assay found, use the first available
    if (is.na(assay_name)) {
      assay_name <- available_assays[1]
    }
    message(paste("Using assay:", assay_name))
  }
  
  if (!assay_name %in% assayNames(sce)) {
    stop(paste("assay_name '", assay_name, "' not found in object. Available assays: ",
               paste(assayNames(sce), collapse = ", "), sep = ""))
  }
  
  # Check if genes exist
  available_genes <- rownames(sce)
  genes_found <- genes %in% available_genes
  if (!any(genes_found)) {
    stop("None of the specified genes found in the dataset")
  }
  if (!all(genes_found)) {
    warning(paste("Genes not found:", paste(genes[!genes_found], collapse = ", ")))
    genes <- genes[genes_found]
  }
  
  # Extract expression for specified genes
  expr_matrix <- assay(sce, assay_name)[genes, , drop = FALSE]
  
  # Aggregate expression across selected genes
  if (length(genes) > 1) {
    aggregated_expr <- apply(expr_matrix, 2, aggregate_fun, na.rm = TRUE)
  } else {
    aggregated_expr <- as.numeric(expr_matrix[1, ])
  }
  
  # Create plot data frame
  plot_df <- data.frame(
    cell_id = colnames(sce),
    expression = aggregated_expr,
    group = colData(sce)[[group_var]],
    stringsAsFactors = FALSE
  )
  
  # Set default title if not provided
  if (is.null(plot_title)) {
    if (length(genes) == 1) {
      plot_title <- paste(genes, "Expression")
    } else {
      plot_title <- paste("Average Expression of", length(genes), "Genes")
    }
  }
  
  # Filter groups if specified
  if (!is.null(groups)) {
    plot_df <- plot_df[plot_df$group %in% groups, ]
    plot_df$group <- factor(plot_df$group, levels = groups)
  } else {
    plot_df$group <- factor(plot_df$group)
  }
  
  # Create violin plot
  p <- ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = colors) +
    labs(title = plot_title,
         x = "Group",
         y = y_label) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none"
    )
  
  return(p)
}