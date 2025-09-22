
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

###############################
# FUNCTION TO DO SIMULATION.  #
###############################
simulate_patient_control_study <- function(
    n_controls = 4,
    n_patients = 4,
    cells_per_ind = 150,
    cells_sd = NULL,           # SD for cell count variation
    cells_min = 30,            # Minimum cells per individual
    cells_max = NULL,          # Maximum cells per individual
    n_genes = 3000,
    individual_variation = 0.1,
    disease_effect_size = 0.5,
    disease_gene_fraction = 0.2,
    verbose = TRUE,
    seed = 123) {              # For reproducibility
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Simulate base data
  n_individuals <- n_controls + n_patients
  
  # Handle cell count variability
  if(length(cells_per_ind) == 1) {
    # If single value provided, generate variable cell counts
    
    # Default SD is 30% of mean if not specified
    if(is.null(cells_sd)) {
      cells_sd <- cells_per_ind * 0.3
    }
    
    # Generate cell counts with variation
    if(cells_sd > 0) {
      cell_counts <- round(rnorm(n_individuals, 
                                 mean = cells_per_ind, 
                                 sd = cells_sd))
    } else {
      # No variation
      cell_counts <- rep(cells_per_ind, n_individuals)
    }
    
    # Apply minimum threshold
    cell_counts[cell_counts < cells_min] <- cells_min
    
    # Apply maximum threshold if specified
    if(!is.null(cells_max)) {
      cell_counts[cell_counts > cells_max] <- cells_max
    }
    
  } else if(length(cells_per_ind) == n_individuals) {
    # If vector provided, use directly (allows manual specification)
    cell_counts <- cells_per_ind
    
    # Still apply min/max thresholds
    cell_counts[cell_counts < cells_min] <- cells_min
    if(!is.null(cells_max)) {
      cell_counts[cell_counts > cells_max] <- cells_max
    }
    
  } else if(length(cells_per_ind) == 2) {
    # If two values provided, use different means for controls and patients
    control_cells <- round(rnorm(n_controls, 
                                 mean = cells_per_ind[1], 
                                 sd = ifelse(is.null(cells_sd), cells_per_ind[1] * 0.3, cells_sd)))
    patient_cells <- round(rnorm(n_patients, 
                                 mean = cells_per_ind[2], 
                                 sd = ifelse(is.null(cells_sd), cells_per_ind[2] * 0.3, cells_sd)))
    
    cell_counts <- c(control_cells, patient_cells)
    
    # Apply thresholds
    cell_counts[cell_counts < cells_min] <- cells_min
    if(!is.null(cells_max)) {
      cell_counts[cell_counts > cells_max] <- cells_max
    }
    
  } else {
    stop("cells_per_ind must be a single value, a vector of length 2 (control, patient means), or a vector of length n_individuals")
  }
  
  # Individual-specific expression variation
  ind_effects <- rnorm(n_individuals, 0, individual_variation)
  
  # Print summary if verbose
  if(verbose) {
    cat("=== Simulation Parameters ===\n")
    cat("Genes:", n_genes, "\n")
    cat("Individuals:", n_controls, "controls,", n_patients, "patients\n")
    cat("Individual variation:", individual_variation, "\n")
    cat("Disease effect size:", disease_effect_size, "\n")
    cat("Disease gene fraction:", disease_gene_fraction, "\n\n")
    
    cat("Cell distribution:\n")
    control_counts <- cell_counts[1:n_controls]
    patient_counts <- cell_counts[(n_controls+1):n_individuals]
    
    cat("  Controls:\n")
    cat("    Mean:", round(mean(control_counts), 1), "cells/individual\n")
    cat("    SD:", round(sd(control_counts), 1), "\n")
    cat("    Range:", min(control_counts), "-", max(control_counts), "\n")
    cat("    Total:", sum(control_counts), "cells\n")
    
    cat("  Patients:\n")
    cat("    Mean:", round(mean(patient_counts), 1), "cells/individual\n")
    cat("    SD:", round(sd(patient_counts), 1), "\n")
    cat("    Range:", min(patient_counts), "-", max(patient_counts), "\n")
    cat("    Total:", sum(patient_counts), "cells\n")
    
    cat("  Overall:\n")
    cat("    Total cells:", sum(cell_counts), "\n")
    cat("    CV of cell counts:", round(sd(cell_counts)/mean(cell_counts), 3), "\n\n")
  }
  
  # Create simulation with variable cell counts
  sim <- splatSimulate(
    nGenes = n_genes,
    batchCells = cell_counts,  # Variable counts per individual
    batch.facLoc = ind_effects,
    batch.facScale = rep(0.1, n_individuals),
    verbose = FALSE
  )
  
  # Create individual and group assignments
  unique_batches <- sort(unique(colData(sim)$Batch))
  batch_to_individual <- seq_along(unique_batches)
  names(batch_to_individual) <- as.character(unique_batches)
  
  # Map each cell to its individual number
  cell_individual_num <- batch_to_individual[as.character(colData(sim)$Batch)]
  
  # Create individual and group assignments
  colData(sim)$Individual <- factor(paste0("Ind_", cell_individual_num))
  
  # Assign groups
  individual_groups <- c(rep("Control", n_controls), rep("Patient", n_patients))
  colData(sim)$Group <- factor(individual_groups[cell_individual_num])
  
  # Add actual cell count to metadata for reference
  colData(sim)$CellsPerIndividual <- cell_counts[cell_individual_num]
  
  # Verify assignment if verbose - FIXED HERE
  if(verbose) {
    cat("Individual-Group Assignment:\n")
    
    # Create summary properly
    ind_table <- table(colData(sim)$Individual)
    
    assignment_check <- data.frame(
      Individual = names(ind_table),
      Group = individual_groups,  # Direct mapping since individuals are in order
      nCells = as.numeric(ind_table),
      stringsAsFactors = FALSE
    )
    
    print(assignment_check)
    cat("\n")
  }
  
  # Add disease effect to patient individuals
  if(disease_effect_size > 0 && disease_gene_fraction > 0) {
    counts_mat <- counts(sim)
    
    n_disease_genes <- max(1, round(n_genes * disease_gene_fraction))
    disease_genes <- sample(1:n_genes, size = n_disease_genes)
    
    patient_cells <- which(colData(sim)$Group == "Patient")
    
    if(length(patient_cells) > 0) {
      effect_multiplier <- exp(rnorm(length(disease_genes), 
                                     mean = log(disease_effect_size), 
                                     sd = 0.2))
      
      for(i in seq_along(disease_genes)) {
        gene_idx <- disease_genes[i]
        counts_mat[gene_idx, patient_cells] <- 
          round(counts_mat[gene_idx, patient_cells] * effect_multiplier[i])
      }
      
      assay(sim, "counts") <- counts_mat
      
      # Mark disease genes
      rowData(sim)$is_disease_gene <- FALSE
      rowData(sim)$is_disease_gene[disease_genes] <- TRUE
      rowData(sim)$disease_effect <- 1
      rowData(sim)$disease_effect[disease_genes] <- effect_multiplier
    }
  }
  
  # Add metadata about the simulation
  metadata(sim) <- list(
    n_controls = n_controls,
    n_patients = n_patients,
    n_genes = n_genes,
    individual_variation = individual_variation,
    disease_effect_size = disease_effect_size,
    disease_gene_fraction = disease_gene_fraction,
    cell_counts = cell_counts,
    total_cells = sum(cell_counts)
  )
  
  return(sim)
}

###########################
# PLOT DEGS               #
###########################

# If you have disease genes marked in rowData
plot_disease_genes <- function(sim, n_genes = 50) {
  # Get disease genes if they're marked
  if("is_disease_gene" %in% colnames(rowData(sim))) {
    disease_genes <- rownames(sim)[rowData(sim)$is_disease_gene]
    
    # Take first n_genes
    genes_to_plot <- disease_genes[1:min(n_genes, length(disease_genes))]
    
    p <- plot_gene_expression_violin(
      sim,
      genes = genes_to_plot,
      group_var = "Group",
      plot_title = sprintf("Disease Gene Expression\n(%d genes with designed effect)", 
                           length(genes_to_plot)),
      colors = c("Control" = "#377EB8", "Patient" = "#E41A1C")
    )
    
    return(p)
  } else {
    # Find differential genes empirically
    top_genes <- find_batch_effect_genes(sim, "Group", n_genes)
    
    p <- plot_gene_expression_violin(
      sim,
      genes = top_genes,
      group_var = "Group", 
      plot_title = sprintf("Top %d Differential Genes", n_genes),
      colors = c("#377EB8", "#E41A1C")
    )
    
    return(p)
  }
}