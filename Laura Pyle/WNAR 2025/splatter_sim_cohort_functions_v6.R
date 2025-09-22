
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
    within_ind_correlation = 0.4,  # NEW: correlation between cells within individual
    n_cell_states = 5,             # NEW: number of cell states per individual
    verbose = TRUE,
    seed = 123) {              
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Simulate base data
  n_individuals <- n_controls + n_patients
  
  # Handle cell count variability
  if(length(cells_per_ind) == 1) {
    if(is.null(cells_sd)) {
      cells_sd <- cells_per_ind * 0.3
    }
    
    if(cells_sd > 0) {
      cell_counts <- round(rnorm(n_individuals, 
                                 mean = cells_per_ind, 
                                 sd = cells_sd))
    } else {
      cell_counts <- rep(cells_per_ind, n_individuals)
    }
    
    cell_counts[cell_counts < cells_min] <- cells_min
    
    if(!is.null(cells_max)) {
      cell_counts[cell_counts > cells_max] <- cells_max
    }
    
  } else if(length(cells_per_ind) == n_individuals) {
    cell_counts <- cells_per_ind
    cell_counts[cell_counts < cells_min] <- cells_min
    if(!is.null(cells_max)) {
      cell_counts[cell_counts > cells_max] <- cells_max
    }
    
  } else if(length(cells_per_ind) == 2) {
    control_cells <- round(rnorm(n_controls, 
                                 mean = cells_per_ind[1], 
                                 sd = ifelse(is.null(cells_sd), cells_per_ind[1] * 0.3, cells_sd)))
    patient_cells <- round(rnorm(n_patients, 
                                 mean = cells_per_ind[2], 
                                 sd = ifelse(is.null(cells_sd), cells_per_ind[2] * 0.3, cells_sd)))
    
    cell_counts <- c(control_cells, patient_cells)
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
    cat("Within-individual correlation:", within_ind_correlation, "\n")  # NEW
    cat("Cell states per individual:", n_cell_states, "\n")              # NEW
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
    batchCells = cell_counts,
    batch.facLoc = ind_effects,
    batch.facScale = rep(0.1, n_individuals),
    verbose = FALSE
  )
  
  # ADD WITHIN-INDIVIDUAL CORRELATION
  if(within_ind_correlation > 0 && n_cell_states > 1) {
    sim <- add_cell_state_correlation(sim, 
                                      n_states = n_cell_states,
                                      state_correlation = within_ind_correlation)
  }
  
  # Create individual and group assignments
  unique_batches <- sort(unique(colData(sim)$Batch))
  batch_to_individual <- seq_along(unique_batches)
  names(batch_to_individual) <- as.character(unique_batches)
  
  cell_individual_num <- batch_to_individual[as.character(colData(sim)$Batch)]
  
  colData(sim)$Individual <- factor(paste0("Ind_", cell_individual_num))
  
  individual_groups <- c(rep("Control", n_controls), rep("Patient", n_patients))
  colData(sim)$Group <- factor(individual_groups[cell_individual_num])
  
  colData(sim)$CellsPerIndividual <- cell_counts[cell_individual_num]
  
  # Verify assignment if verbose
  if(verbose) {
    cat("Individual-Group Assignment:\n")
    
    ind_table <- table(colData(sim)$Individual)
    
    assignment_check <- data.frame(
      Individual = names(ind_table),
      Group = individual_groups,
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
    within_ind_correlation = within_ind_correlation,
    n_cell_states = n_cell_states,
    disease_effect_size = disease_effect_size,
    disease_gene_fraction = disease_gene_fraction,
    cell_counts = cell_counts,
    total_cells = sum(cell_counts)
  )
  
  return(sim)
}

# Include the correlation function (with fixes from earlier)
add_cell_state_correlation <- function(sim, 
                                       n_states = 3,
                                       state_correlation = 0.5) {
  
  # Validate parameters
  if(state_correlation < 0 || state_correlation > 1) {
    stop("state_correlation must be between 0 and 1")
  }
  
  counts_original <- counts(sim)
  
  for(batch in unique(colData(sim)$Batch)) {
    batch_cells <- which(colData(sim)$Batch == batch)
    n_cells <- length(batch_cells)
    
    # Assign cells to states (some overlap creates correlation)
    if(state_correlation > 0 && n_cells > 1) {  # Need at least 2 cells
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
      
      # Create state-specific expression programs (reduced variance)
      state_programs <- matrix(rnorm(nrow(counts_original) * n_states, 0, 0.1),
                               nrow = nrow(counts_original))
      
      # Apply weighted combination of programs
      for(i in 1:n_cells) {
        cell_program <- as.vector(state_programs %*% state_weights[i, ])
        counts_original[, batch_cells[i]] <- 
          counts_original[, batch_cells[i]] * exp(cell_program)
      }
    }
  }
  
  assay(sim, "counts") <- round(counts_original)
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

################################################################
# SIMULATE COHORT OVER A RANGE OF PARAMETERS FOR N_SIMULATIONS #
################################################################

simulate_parameter_sweep <- function(
    n_controls = 4,
    n_patients = 4,
    cells_per_ind = 150,
    cells_sd = NULL,
    cells_min = 30,
    cells_max = NULL,
    n_genes = 3000,
    individual_variation = 0.1,
    disease_effect_size = 0.5,
    disease_gene_fraction = 0.2,
    within_ind_correlation = 0.4,
    n_cell_states = 5,
    n_simulations = 1,
    verbose = FALSE,
    parallel = FALSE,
    n_cores = NULL,
    seed_start = 123,
    return_summary = TRUE,
    experiment_name = NULL) {  # NEW: optional experiment name
  
  # Generate experiment ID if not provided
  if(is.null(experiment_name)) {
    experiment_name <- paste0("exp_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # Convert single values to vectors for consistency
  make_vector <- function(x) {
    if(length(x) == 1) return(x)
    return(x)
  }
  
  # Create parameter grid
  param_grid <- expand.grid(
    n_controls = make_vector(n_controls),
    n_patients = make_vector(n_patients),
    cells_per_ind = make_vector(cells_per_ind),
    cells_sd = if(is.null(cells_sd)) NA else make_vector(cells_sd),
    cells_min = make_vector(cells_min),
    cells_max = if(is.null(cells_max)) NA else make_vector(cells_max),
    n_genes = make_vector(n_genes),
    individual_variation = make_vector(individual_variation),
    disease_effect_size = make_vector(disease_effect_size),
    disease_gene_fraction = make_vector(disease_gene_fraction),
    within_ind_correlation = make_vector(within_ind_correlation),
    n_cell_states = make_vector(n_cell_states),
    stringsAsFactors = FALSE
  )
  
  # Add simulation replicate column and unique IDs
  param_grid <- param_grid[rep(seq_len(nrow(param_grid)), each = n_simulations), ]
  param_grid$simulation_rep <- rep(1:n_simulations, nrow(param_grid) / n_simulations)
  param_grid$param_set_id <- rep(1:(nrow(param_grid) / n_simulations), each = n_simulations)
  param_grid$seed <- seed_start + seq_len(nrow(param_grid)) - 1
  
  # CREATE UNIQUE SIMULATION IDs
  param_grid$simulation_id <- paste0(
    experiment_name, "_",
    "ps", formatC(param_grid$param_set_id, width = 3, flag = "0"), "_",
    "rep", formatC(param_grid$simulation_rep, width = 2, flag = "0")
  )
  
  # Create short hash for each unique parameter combination
  param_grid$param_hash <- apply(param_grid[, 1:12], 1, function(row) {
    substr(digest::digest(row, algo = "md5"), 1, 8)
  })
  
  total_sims <- nrow(param_grid)
  
  cat("=== Parameter Sweep Setup ===\n")
  cat("Experiment name:", experiment_name, "\n")
  cat("Total parameter combinations:", nrow(param_grid) / n_simulations, "\n")
  cat("Simulations per combination:", n_simulations, "\n")
  cat("Total simulations to run:", total_sims, "\n\n")
  
  if(verbose && total_sims <= 10) {
    cat("Simulation IDs:\n")
    print(data.frame(
      sim_id = param_grid$simulation_id,
      param_set = param_grid$param_set_id,
      rep = param_grid$simulation_rep
    ))
    cat("\n")
  }
  
  # Function to run single simulation
  run_single_sim <- function(i) {
    params <- param_grid[i, ]
    
    if(verbose) {
      cat(sprintf("Running %s (param set %d, rep %d)...\n", 
                  params$simulation_id, params$param_set_id, params$simulation_rep))
    }
    
    # Handle NA values for optional parameters
    cells_sd_val <- if(is.na(params$cells_sd)) NULL else params$cells_sd
    cells_max_val <- if(is.na(params$cells_max)) NULL else params$cells_max
    
    # Run simulation
    sim <- simulate_patient_control_study(
      n_controls = params$n_controls,
      n_patients = params$n_patients,
      cells_per_ind = params$cells_per_ind,
      cells_sd = cells_sd_val,
      cells_min = params$cells_min,
      cells_max = cells_max_val,
      n_genes = params$n_genes,
      individual_variation = params$individual_variation,
      disease_effect_size = params$disease_effect_size,
      disease_gene_fraction = params$disease_gene_fraction,
      within_ind_correlation = params$within_ind_correlation,
      n_cell_states = params$n_cell_states,
      verbose = FALSE,
      seed = params$seed
    )
    
    # ADD UNIQUE IDs TO THE SIMULATION OBJECT ITSELF
    colData(sim)$simulation_id <- params$simulation_id
    colData(sim)$param_set_id <- params$param_set_id
    colData(sim)$simulation_rep <- params$simulation_rep
    colData(sim)$param_hash <- params$param_hash
    
    # Extract summary statistics
    if(return_summary) {
      # Calculate summary metrics
      disease_genes <- rownames(sim)[rowData(sim)$is_disease_gene]
      
      # Get expression differences for disease genes
      if(length(disease_genes) > 0) {
        disease_expr <- counts(sim)[disease_genes, , drop = FALSE]
        
        control_cells <- colData(sim)$Group == "Control"
        patient_cells <- colData(sim)$Group == "Patient"
        
        mean_control <- rowMeans(disease_expr[, control_cells, drop = FALSE])
        mean_patient <- rowMeans(disease_expr[, patient_cells, drop = FALSE])
        
        log_fc <- log2((mean_patient + 1) / (mean_control + 1))
        mean_log_fc <- mean(log_fc)
        median_log_fc <- median(log_fc)
      } else {
        mean_log_fc <- NA
        median_log_fc <- NA
      }
      
      # Cell count statistics
      cell_counts_by_ind <- table(colData(sim)$Individual)
      
      summary_stats <- data.frame(
        # UNIQUE IDENTIFIERS
        simulation_id = params$simulation_id,
        param_set_id = params$param_set_id,
        simulation_rep = params$simulation_rep,
        param_hash = params$param_hash,
        experiment_name = experiment_name,
        seed = params$seed,
        # Parameters
        n_controls = params$n_controls,
        n_patients = params$n_patients,
        cells_per_ind = params$cells_per_ind,
        cells_sd = params$cells_sd,
        cells_min = params$cells_min,
        cells_max = params$cells_max,
        n_genes = params$n_genes,
        individual_variation = params$individual_variation,
        disease_effect_size = params$disease_effect_size,
        disease_gene_fraction = params$disease_gene_fraction,
        within_ind_correlation = params$within_ind_correlation,
        n_cell_states = params$n_cell_states,
        # Results
        total_cells = ncol(sim),
        actual_disease_genes = sum(rowData(sim)$is_disease_gene),
        mean_cells_per_ind = mean(cell_counts_by_ind),
        sd_cells_per_ind = sd(cell_counts_by_ind),
        cv_cells_per_ind = sd(cell_counts_by_ind) / mean(cell_counts_by_ind),
        mean_log_fc_disease = mean_log_fc,
        median_log_fc_disease = median_log_fc,
        # Timestamp
        run_timestamp = Sys.time(),
        stringsAsFactors = FALSE
      )
      
      return(summary_stats)
    } else {
      # Return full simulation object with all identifiers in metadata
      metadata(sim)$simulation_id <- params$simulation_id
      metadata(sim)$param_set_id <- params$param_set_id
      metadata(sim)$simulation_rep <- params$simulation_rep
      metadata(sim)$param_hash <- params$param_hash
      metadata(sim)$experiment_name <- experiment_name
      metadata(sim)$parameters <- as.list(params)
      metadata(sim)$run_timestamp <- Sys.time()
      return(sim)
    }
  }
  
  # Run simulations (parallel/sequential code remains the same)
  if(parallel && total_sims > 1) {
    if(!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not available. Running sequentially.")
      parallel <- FALSE
    }
  }
  
  if(parallel && total_sims > 1) {
    # Parallel execution
    if(is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    n_cores <- min(n_cores, total_sims)
    
    cat("Running in parallel with", n_cores, "cores...\n")
    
    cl <- parallel::makeCluster(n_cores)
    
    parallel::clusterExport(cl, c("simulate_patient_control_study", 
                                  "add_cell_state_correlation",
                                  "param_grid",
                                  "experiment_name"), 
                            envir = environment())
    
    parallel::clusterEvalQ(cl, {
      library(splatter)
      library(SingleCellExperiment)
      library(digest)
    })
    
    results <- parallel::parLapply(cl, 1:total_sims, run_single_sim)
    parallel::stopCluster(cl)
    
  } else {
    # Sequential execution
    results <- lapply(1:total_sims, run_single_sim)
  }
  
  # Combine results
  if(return_summary) {
    # Combine summary statistics
    summary_df <- do.call(rbind, results)
    
    # Add experiment metadata as attributes
    attr(summary_df, "experiment_name") <- experiment_name
    attr(summary_df, "run_date") <- Sys.Date()
    attr(summary_df, "parameter_grid") <- param_grid
    
    # Add aggregate statistics across simulations if multiple reps
    if(n_simulations > 1) {
      cat("\n=== Aggregating results across simulations ===\n")
      
      library(dplyr)
      
      aggregate_stats <- summary_df %>%
        group_by(param_set_id, param_hash) %>%
        summarise(
          # Keep one copy of parameters
          n_controls = first(n_controls),
          n_patients = first(n_patients),
          cells_per_ind = first(cells_per_ind),
          individual_variation = first(individual_variation),
          disease_effect_size = first(disease_effect_size),
          disease_gene_fraction = first(disease_gene_fraction),
          within_ind_correlation = first(within_ind_correlation),
          n_cell_states = first(n_cell_states),
          # Aggregate statistics
          n_reps = n(),
          mean_total_cells = mean(total_cells),
          sd_total_cells = sd(total_cells),
          mean_log_fc = mean(mean_log_fc_disease, na.rm = TRUE),
          sd_log_fc = sd(mean_log_fc_disease, na.rm = TRUE),
          mean_cv_cells = mean(cv_cells_per_ind),
          .groups = "drop"
        )
      
      return(list(
        detailed = summary_df,
        aggregate = aggregate_stats,
        experiment_info = list(
          experiment_name = experiment_name,
          total_simulations = total_sims,
          parameter_combinations = nrow(aggregate_stats),
          replications = n_simulations,
          run_date = Sys.Date()
        )
      ))
    }
    
    return(summary_df)
    
  } else {
    # Return list of simulation objects with unique names
    names(results) <- param_grid$simulation_id
    
    # Add parameter summary as attribute
    attr(results, "parameters") <- param_grid
    attr(results, "experiment_name") <- experiment_name
    
    return(results)
  }
}

# Helper function to extract simulation by ID
get_simulation_by_id <- function(sweep_results, sim_id) {
  if(is.list(sweep_results) && sim_id %in% names(sweep_results)) {
    return(sweep_results[[sim_id]])
  } else {
    stop("Simulation ID not found: ", sim_id)
  }
}

# Function to compare specific simulations
compare_simulations <- function(sweep_results, sim_ids) {
  library(dplyr)
  
  if(is.data.frame(sweep_results) || "detailed" %in% names(sweep_results)) {
    # Working with summary data
    data <- if("detailed" %in% names(sweep_results)) {
      sweep_results$detailed
    } else {
      sweep_results
    }
    
    comparison <- data %>%
      filter(simulation_id %in% sim_ids) %>%
      select(simulation_id, param_set_id, simulation_rep,
             disease_effect_size, individual_variation,
             mean_log_fc_disease, total_cells)
    
    print(comparison)
    
  } else {
    # Working with full simulation objects
    for(id in sim_ids) {
      sim <- sweep_results[[id]]
      cat("\n", id, ":\n")
      cat("  Total cells:", ncol(sim), "\n")
      cat("  Disease genes:", sum(rowData(sim)$is_disease_gene), "\n")
      cat("  Parameters:\n")
      print(metadata(sim)$parameters[c("disease_effect_size", 
                                       "individual_variation")])
    }
  }
}