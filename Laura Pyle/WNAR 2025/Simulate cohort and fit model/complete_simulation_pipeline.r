# Complete Pipeline Script for Simulation and Model Fitting
# All required functions included

# Load required libraries
library(splatter)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(limma)
library(edgeR)
library(parallel)

#########################################
# SECTION 1: CORE SIMULATION FUNCTIONS #
#########################################

# Function to add correlation between cells within individuals
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
    if(state_correlation > 0 && n_cells > 1) {
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

# Main simulation function
simulate_patient_control_study <- function(
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
    stop("cells_per_ind must be a single value, a vector of length 2, or a vector of length n_individuals")
  }
  
  # Individual-specific expression variation
  ind_effects <- rnorm(n_individuals, 0, individual_variation)
  
  # Print summary if verbose
  if(verbose) {
    cat("=== Simulation Parameters ===\n")
    cat("Genes:", n_genes, "\n")
    cat("Individuals:", n_controls, "controls,", n_patients, "patients\n")
    cat("Total cells:", sum(cell_counts), "\n\n")
  }
  
  # Create simulation with variable cell counts
  sim <- splatSimulate(
    nGenes = n_genes,
    batchCells = cell_counts,
    batch.facLoc = ind_effects,
    batch.facScale = rep(0.1, n_individuals),
    verbose = FALSE
  )
  
  # Add within-individual correlation
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
    within_ind_correlation = within_ind_correlation,
    n_cell_states = n_cell_states,
    disease_effect_size = disease_effect_size,
    disease_gene_fraction = disease_gene_fraction,
    cell_counts = cell_counts,
    total_cells = sum(cell_counts)
  )
  
  return(sim)
}

#####################################
# SECTION 2: MODEL FITTING FUNCTIONS #
#####################################

# Default DE model fitting function
fit_default_de_model <- function(sim, params) {
  # Simple differential expression using limma
  counts_mat <- counts(sim)
  
  # Filter low-expressed genes
  keep <- rowSums(counts_mat > 1) >= 10
  counts_mat <- counts_mat[keep, ]
  
  # Create design matrix
  design <- model.matrix(~ colData(sim)$Group)
  
  # Normalize and fit
  dge <- DGEList(counts = counts_mat)
  dge <- calcNormFactors(dge)
  
  # Voom transformation
  v <- voom(dge, design, plot = FALSE)
  
  # Fit linear model
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Get results
  results <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  
  # Calculate summary statistics
  n_sig_005 <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
  n_sig_001 <- sum(results$adj.P.Val < 0.01, na.rm = TRUE)
  mean_logfc <- mean(abs(results$logFC), na.rm = TRUE)
  
  # If we know which are disease genes, calculate power
  power_005 <- NA
  power_001 <- NA
  
  if("is_disease_gene" %in% colnames(rowData(sim))) {
    disease_genes <- rownames(sim)[rowData(sim)$is_disease_gene]
    disease_genes_tested <- intersect(disease_genes, rownames(results))
    
    if(length(disease_genes_tested) > 0) {
      disease_results <- results[disease_genes_tested, ]
      power_005 <- mean(disease_results$adj.P.Val < 0.05, na.rm = TRUE)
      power_001 <- mean(disease_results$adj.P.Val < 0.01, na.rm = TRUE)
    }
  }
  
  return(data.frame(
    n_genes_tested = nrow(results),
    n_sig_005 = n_sig_005,
    n_sig_001 = n_sig_001,
    mean_abs_logfc = mean_logfc,
    power_005 = power_005,
    power_001 = power_001
  ))
}

# Pseudobulk model
fit_pseudobulk_model <- function(sim, params) {
  library(edgeR)
  
  # Aggregate to pseudobulk
  counts_mat <- counts(sim)
  
  pseudobulk <- aggregate(t(counts_mat), 
                          by = list(Individual = colData(sim)$Individual),
                          FUN = sum)
  
  # Reshape
  pb_counts <- t(pseudobulk[, -1])
  colnames(pb_counts) <- pseudobulk$Individual
  
  # Get group information
  ind_groups <- unique(data.frame(
    Individual = colData(sim)$Individual,
    Group = colData(sim)$Group
  ))
  
  # Ensure order matches
  ind_groups <- ind_groups[match(colnames(pb_counts), ind_groups$Individual), ]
  
  # DE analysis
  y <- DGEList(counts = pb_counts)
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ ind_groups$Group)
  y <- estimateDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  
  results <- topTags(lrt, n = Inf)$table
  
  return(data.frame(
    n_genes_pb = nrow(results),
    n_sig_pb_005 = sum(results$FDR < 0.05),
    n_sig_pb_001 = sum(results$FDR < 0.01),
    mean_logfc_pb = mean(abs(results$logFC))
  ))
}

##############################
# SECTION 3: MAIN PIPELINE #
##############################

run_simulation_analysis_pipeline <- function(
    # Simulation parameters
  param_grid = NULL,
  n_simulations = 10,
  experiment_name = NULL,
  # Model fitting function
  fit_function = NULL,
  # Storage options
  output_dir = "simulation_results",
  save_intermediate = FALSE,
  batch_size = 100,
  # Computation options
  parallel = TRUE,
  n_cores = NULL,
  verbose = TRUE,
  seed_start = 123) {
  
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate experiment name if not provided
  if(is.null(experiment_name)) {
    experiment_name <- paste0("exp_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # Create default parameter grid if not provided
  if(is.null(param_grid)) {
    param_grid <- expand.grid(
      n_controls = 10,
      n_patients = 10,
      cells_per_ind = 150,
      cells_sd = 45,
      cells_min = 30,
      cells_max = 300,
      n_genes = 3000,
      individual_variation = c(0.1, 0.3),
      disease_effect_size = c(2, 5, 10),
      disease_gene_fraction = 0.2,
      within_ind_correlation = 0.4,
      n_cell_states = 5,
      stringsAsFactors = FALSE
    )
  }
  
  # Ensure all required columns exist
  required_cols <- c("n_controls", "n_patients", "cells_per_ind", "n_genes",
                     "individual_variation", "disease_effect_size", 
                     "disease_gene_fraction")
  
  # Add default values for missing columns
  default_vals <- list(
    n_controls = 10, n_patients = 10, cells_per_ind = 150,
    cells_sd = 45, cells_min = 30, cells_max = 300,
    n_genes = 3000, individual_variation = 0.1,
    disease_effect_size = 2, disease_gene_fraction = 0.2,
    within_ind_correlation = 0.4, n_cell_states = 5
  )
  
  for(col in names(default_vals)) {
    if(!col %in% names(param_grid)) {
      param_grid[[col]] <- default_vals[[col]]
    }
  }
  
  # Add replicate information
  param_grid <- param_grid[rep(seq_len(nrow(param_grid)), each = n_simulations), ]
  param_grid$simulation_rep <- rep(1:n_simulations, nrow(param_grid) / n_simulations)
  param_grid$param_set_id <- rep(1:(nrow(param_grid) / n_simulations), each = n_simulations)
  param_grid$seed <- seed_start + seq_len(nrow(param_grid)) - 1
  param_grid$simulation_id <- paste0(
    experiment_name, "_",
    "ps", formatC(param_grid$param_set_id, width = 3, flag = "0"), "_",
    "rep", formatC(param_grid$simulation_rep, width = 2, flag = "0")
  )
  
  total_sims <- nrow(param_grid)
  n_batches <- ceiling(total_sims / batch_size)
  
  if(verbose) {
    cat("=== Simulation & Analysis Pipeline ===\n")
    cat("Experiment:", experiment_name, "\n")
    cat("Total simulations:", total_sims, "\n")
    cat("Batch size:", batch_size, "\n")
    cat("Number of batches:", n_batches, "\n")
    cat("Output directory:", output_dir, "\n\n")
  }
  
  # Use default fit function if none provided
  if(is.null(fit_function)) {
    fit_function <- fit_default_de_model
  }
  
  # Function to process single simulation
  process_single_simulation <- function(i) {
    params <- param_grid[i, ]
    
    tryCatch({
      # Generate simulation
      sim <- simulate_patient_control_study(
        n_controls = params$n_controls,
        n_patients = params$n_patients,
        cells_per_ind = params$cells_per_ind,
        cells_sd = if(is.na(params$cells_sd)) NULL else params$cells_sd,
        cells_min = params$cells_min,
        cells_max = if(is.na(params$cells_max)) NULL else params$cells_max,
        n_genes = params$n_genes,
        individual_variation = params$individual_variation,
        disease_effect_size = params$disease_effect_size,
        disease_gene_fraction = params$disease_gene_fraction,
        within_ind_correlation = params$within_ind_correlation,
        n_cell_states = params$n_cell_states,
        verbose = FALSE,
        seed = params$seed
      )
      
      # Add simulation ID to the object
      colData(sim)$simulation_id <- params$simulation_id
      
      # Fit model
      model_results <- fit_function(sim, params)
      
      # Combine parameters and results
      results <- cbind(params, model_results)
      
      # Optionally save intermediate simulation
      if(save_intermediate) {
        sim_file <- file.path(output_dir, "simulations", 
                              paste0(params$simulation_id, ".rds"))
        if(!dir.exists(dirname(sim_file))) {
          dir.create(dirname(sim_file), recursive = TRUE)
        }
        saveRDS(sim, sim_file)
      }
      
      # Clear memory
      rm(sim)
      gc(verbose = FALSE)
      
      return(results)
      
    }, error = function(e) {
      warning(sprintf("Error in simulation %s: %s", 
                      params$simulation_id, e$message))
      return(data.frame(simulation_id = params$simulation_id,
                        error = e$message))
    })
  }
  
  # Process in batches
  all_results <- list()
  
  for(batch in 1:n_batches) {
    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, total_sims)
    batch_indices <- batch_start:batch_end
    
    if(verbose) {
      cat(sprintf("\nProcessing batch %d/%d (simulations %d-%d)\n",
                  batch, n_batches, batch_start, batch_end))
    }
    
    # Process batch
    if(parallel && length(batch_indices) > 1) {
      # Parallel processing
      if(is.null(n_cores)) {
        n_cores <- parallel::detectCores() - 1
      }
      n_cores <- min(n_cores, length(batch_indices))
      
      cl <- parallel::makeCluster(n_cores)
      
      # Export required objects and functions
      parallel::clusterExport(cl, 
                              c("simulate_patient_control_study", 
                                "add_cell_state_correlation",
                                "fit_function",
                                "fit_default_de_model",
                                "param_grid",
                                "save_intermediate",
                                "output_dir"), 
                              envir = environment())
      
      # Load required packages on workers
      parallel::clusterEvalQ(cl, {
        library(splatter)
        library(SingleCellExperiment)
        library(limma)
        suppressPackageStartupMessages(library(edgeR))
      })
      
      batch_results <- parallel::parLapply(cl, batch_indices, 
                                           process_single_simulation)
      parallel::stopCluster(cl)
      
    } else {
      # Sequential processing
      batch_results <- lapply(batch_indices, process_single_simulation)
    }
    
    # Combine batch results
    batch_df <- do.call(rbind, batch_results)
    all_results[[batch]] <- batch_df
    
    # Save batch results immediately
    batch_file <- file.path(output_dir, 
                            sprintf("batch_%03d_results.rds", batch))
    saveRDS(batch_df, batch_file)
    
    if(verbose) {
      cat(sprintf("  Batch %d completed and saved to %s\n", 
                  batch, batch_file))
    }
    
    # Clear memory
    gc(verbose = FALSE)
  }
  
  # Combine all results
  final_results <- do.call(rbind, all_results)
  
  # Save combined results
  final_file <- file.path(output_dir, "combined_results.rds")
  saveRDS(final_results, final_file)
  
  # Save as CSV for easy viewing
  csv_file <- file.path(output_dir, "combined_results.csv")
  write.csv(final_results, csv_file, row.names = FALSE)
  
  # Save parameter grid
  param_file <- file.path(output_dir, "parameter_grid.rds")
  saveRDS(param_grid, param_file)
  
  if(verbose) {
    cat("\n=== Pipeline Complete ===\n")
    cat("Results saved to:", output_dir, "\n")
    cat("Total simulations processed:", nrow(final_results), "\n")
    if("error" %in% colnames(final_results)) {
      n_errors <- sum(!is.na(final_results$error))
      cat("Simulations with errors:", n_errors, "\n")
    }
  }
  
  return(final_results)
}

#################################
# SECTION 4: ANALYSIS FUNCTIONS #
#################################

# Function to analyze and visualize results
analyze_pipeline_results <- function(output_dir) {
  # Load combined results
  results <- readRDS(file.path(output_dir, "combined_results.rds"))
  
  library(ggplot2)
  library(dplyr)
  
  # Summarize across replicates
  summary_stats <- results %>%
    group_by(disease_effect_size, individual_variation) %>%
    summarise(
      mean_power = mean(power_005, na.rm = TRUE),
      sd_power = sd(power_005, na.rm = TRUE),
      mean_n_sig = mean(n_sig_005, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Plot power curves
  p1 <- ggplot(summary_stats, 
               aes(x = disease_effect_size, y = mean_power,
                   color = factor(individual_variation))) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_power - sd_power,
                      ymax = mean_power + sd_power),
                  width = 0.2) +
    theme_minimal() +
    labs(title = "Detection Power vs Effect Size",
         x = "Disease Effect Size",
         y = "Power (FDR < 0.05)",
         color = "Individual\nVariation") +
    scale_color_brewer(palette = "Set1")
  
  return(list(summary = summary_stats, plot = p1))
}

###########################
# SECTION 5: RUN EXAMPLE #
###########################

# Example usage
if(FALSE) {  # Set to TRUE to run
  
  # Small test run
  test_results <- run_simulation_analysis_pipeline(
    param_grid = expand.grid(
      n_controls = 10,
      n_patients = 10,
      cells_per_ind = 100,
      disease_effect_size = c(2, 5),
      individual_variation = c(0.1, 0.3),
      disease_gene_fraction = 0.2
    ),
    n_simulations = 3,
    experiment_name = "test_run",
    fit_function = fit_default_de_model,
    output_dir = "results/test",
    batch_size = 10,
    parallel = FALSE,  # Set to TRUE if you have multiple cores
    verbose = TRUE
  )
  
  # Analyze results
  analysis <- analyze_pipeline_results("results/test")
  print(analysis$plot)
}

cat("All functions loaded successfully!\n")
cat("To run a simulation pipeline, use: run_simulation_analysis_pipeline()\n")