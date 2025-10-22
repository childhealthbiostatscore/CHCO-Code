# complete_method_comparison_analysis.R
# Complete script for running simulation-based method comparison
# INCLUDING NEBULA, PSEUDOBULK, AND STANDARD DE

#############################################
# SECTION 1: LOAD LIBRARIES
#############################################

library(splatter)
library(SingleCellExperiment)
library(limma)
library(edgeR)
library(nebula)  # Add nebula
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(parallel)

#############################################
# SECTION 2: CORE SIMULATION FUNCTIONS
#############################################

# Function to add correlation between cells within individuals
add_cell_state_correlation <- function(sim, 
                                       n_states = 3,
                                       state_correlation = 0.5) {
  
  if(state_correlation < 0 || state_correlation > 1) {
    stop("state_correlation must be between 0 and 1")
  }
  
  counts_original <- counts(sim)
  
  for(batch in unique(colData(sim)$Batch)) {
    batch_cells <- which(colData(sim)$Batch == batch)
    n_cells <- length(batch_cells)
    
    if(state_correlation > 0 && n_cells > 1) {
      state_weights <- matrix(runif(n_cells * n_states), 
                              nrow = n_cells, ncol = n_states)
      
      for(i in 2:n_cells) {
        state_weights[i, ] <- (1 - state_correlation) * state_weights[i, ] + 
          state_correlation * state_weights[i-1, ]
      }
      
      state_weights <- state_weights / rowSums(state_weights)
      
      state_programs <- matrix(rnorm(nrow(counts_original) * n_states, 0, 0.1),
                               nrow = nrow(counts_original))
      
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

# Main simulation function with cell count variation
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
    verbose = FALSE,
    seed = 123) {
  
  set.seed(seed)
  n_individuals <- n_controls + n_patients
  
  # Handle cell count variability
  if(length(cells_per_ind) == 1) {
    if(is.null(cells_sd)) cells_sd <- cells_per_ind * 0.3
    
    if(cells_sd > 0) {
      cell_counts <- round(rnorm(n_individuals, mean = cells_per_ind, sd = cells_sd))
    } else {
      cell_counts <- rep(cells_per_ind, n_individuals)
    }
    
    cell_counts[cell_counts < cells_min] <- cells_min
    if(!is.null(cells_max)) cell_counts[cell_counts > cells_max] <- cells_max
    
  } else if(length(cells_per_ind) == 2) {
    # Different means for controls and patients
    control_cells <- round(rnorm(n_controls, 
                                 mean = cells_per_ind[1], 
                                 sd = ifelse(is.null(cells_sd), cells_per_ind[1] * 0.3, cells_sd)))
    patient_cells <- round(rnorm(n_patients, 
                                 mean = cells_per_ind[2], 
                                 sd = ifelse(is.null(cells_sd), cells_per_ind[2] * 0.3, cells_sd)))
    
    cell_counts <- c(control_cells, patient_cells)
    cell_counts[cell_counts < cells_min] <- cells_min
    if(!is.null(cells_max)) cell_counts[cell_counts > cells_max] <- cells_max
    
  } else {
    cell_counts <- cells_per_ind
  }
  
  # Individual-specific expression variation
  ind_effects <- rnorm(n_individuals, 0, individual_variation)
  
  if(verbose) {
    cat("Simulating", n_controls, "controls and", n_patients, "patients\n")
    cat("Cell counts per individual:", cell_counts, "\n")
    cat("Mean cells:", mean(cell_counts), "SD:", sd(cell_counts), "\n")
  }
  
  # Create simulation
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
  
  # Store actual cell counts
  colData(sim)$CellsPerIndividual <- cell_counts[cell_individual_num]
  
  # Add disease effect
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
    }
  }
  
  # Store simulation parameters in metadata
  metadata(sim) <- list(
    n_controls = n_controls,
    n_patients = n_patients,
    cells_per_ind = cells_per_ind,
    cells_sd = cells_sd,
    cells_min = cells_min,
    cells_max = cells_max,
    actual_cell_counts = cell_counts,
    disease_effect_size = disease_effect_size,
    disease_gene_fraction = disease_gene_fraction,
    individual_variation = individual_variation
  )
  
  return(sim)
}

#############################################
# SECTION 3: MODEL FITTING FUNCTIONS
#############################################

# Standard differential expression
fit_default_de_model <- function(sim, params = NULL) {
  counts_mat <- counts(sim)
  
  # Filter low-expressed genes
  keep <- rowSums(counts_mat > 1) >= 10
  counts_mat <- counts_mat[keep, ]
  
  # IMPORTANT: Also filter rowData to keep disease gene labels aligned
  filtered_rowdata <- rowData(sim)[keep, , drop = FALSE]
  
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
  
  # Calculate power if we know disease genes
  power_005 <- NA
  power_001 <- NA
  
  # FIX: Use filtered_rowdata instead of rowData(sim)
  if("is_disease_gene" %in% colnames(filtered_rowdata)) {
    disease_genes_in_filtered <- filtered_rowdata$is_disease_gene
    
    if(sum(disease_genes_in_filtered, na.rm = TRUE) > 0) {
      disease_results <- results[disease_genes_in_filtered, ]
      power_005 <- mean(disease_results$adj.P.Val < 0.05, na.rm = TRUE)
      power_001 <- mean(disease_results$adj.P.Val < 0.01, na.rm = TRUE)
    }
  }
  
  # Calculate cell count statistics
  cell_stats <- table(colData(sim)$Individual)
  
  return(data.frame(
    n_genes_tested = nrow(results),
    n_sig_005 = n_sig_005,
    n_sig_001 = n_sig_001,
    mean_abs_logfc = mean_logfc,
    power_005 = power_005,
    power_001 = power_001,
    mean_cells_per_ind = mean(cell_stats),
    sd_cells_per_ind = sd(cell_stats),
    min_cells = min(cell_stats),
    max_cells = max(cell_stats)
  ))
}

# Pseudobulk analysis
fit_pseudobulk_model <- function(sim, params = NULL) {
  # Aggregate to pseudobulk
  counts_mat <- counts(sim)
  
  pseudobulk <- aggregate(t(counts_mat), 
                          by = list(Individual = colData(sim)$Individual),
                          FUN = sum)
  
  pb_counts <- t(pseudobulk[, -1])
  colnames(pb_counts) <- pseudobulk$Individual
  
  # Get group information
  ind_groups <- unique(data.frame(
    Individual = colData(sim)$Individual,
    Group = colData(sim)$Group
  ))
  
  ind_groups <- ind_groups[match(colnames(pb_counts), ind_groups$Individual), ]
  
  # DE analysis - no filtering here initially to preserve all genes
  y <- DGEList(counts = pb_counts)
  
  # FILTER HERE and track which genes
  keep <- filterByExpr(y, min.count = 1)
  y <- y[keep, , keep.lib.sizes = FALSE]
  filtered_rowdata <- rowData(sim)[keep, , drop = FALSE]
  
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ ind_groups$Group)
  y <- estimateDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  
  results <- topTags(lrt, n = Inf)$table
  
  # Calculate power if we know disease genes
  power_pb_005 <- NA
  power_pb_001 <- NA
  
  # FIX: Use filtered_rowdata
  if("is_disease_gene" %in% colnames(filtered_rowdata)) {
    disease_genes_in_filtered <- filtered_rowdata$is_disease_gene
    
    if(sum(disease_genes_in_filtered, na.rm = TRUE) > 0) {
      disease_results <- results[disease_genes_in_filtered, ]
      power_pb_005 <- mean(disease_results$FDR < 0.05, na.rm = TRUE)
      power_pb_001 <- mean(disease_results$FDR < 0.01, na.rm = TRUE)
    }
  }
  
  return(data.frame(
    n_genes_pb = nrow(results),
    n_sig_pb_005 = sum(results$FDR < 0.05),
    n_sig_pb_001 = sum(results$FDR < 0.01),
    mean_logfc_pb = mean(abs(results$logFC)),
    power_pb_005 = power_pb_005,
    power_pb_001 = power_pb_001,
    n_individuals = ncol(pb_counts)
  ))
}
# NEBULA model fit function
fit_nebula_model <- function(sim, params = NULL) {
  
  tryCatch({
    # Prepare data for nebula
    counts_mat <- as.matrix(counts(sim))
    
    # Filter low-expressed genes (same as other methods)
    keep <- rowSums(counts_mat > 1) >= 10
    counts_mat <- counts_mat[keep, ]
    
    # IMPORTANT: Track filtered rowdata
    filtered_rowdata <- rowData(sim)[keep, , drop = FALSE]
    
    # Need at least some genes
    if(nrow(counts_mat) < 10) {
      warning("Too few genes pass filtering")
      return(data.frame(
        n_genes_nebula = NA,
        n_sig_nebula_005 = NA,
        n_sig_nebula_001 = NA,
        n_sig_nebula_raw_005 = NA,
        mean_logfc_nebula = NA,
        power_nebula_005 = NA,
        power_nebula_001 = NA,
        nebula_convergence = NA,
        nebula_method = "TOO_FEW_GENES"
      ))
    }
    
    # Create data frame with cell metadata
    cell_metadata <- data.frame(
      individual = as.character(colData(sim)$Individual),
      group = as.character(colData(sim)$Group),
      row.names = colnames(counts_mat),
      stringsAsFactors = FALSE
    )
    
    # Check we have both groups
    if(length(unique(cell_metadata$group)) < 2) {
      warning("Need at least 2 groups")
      return(data.frame(
        n_genes_nebula = NA,
        n_sig_nebula_005 = NA,
        n_sig_nebula_001 = NA,
        n_sig_nebula_raw_005 = NA,
        mean_logfc_nebula = NA,
        power_nebula_005 = NA,
        power_nebula_001 = NA,
        nebula_convergence = NA,
        nebula_method = "SINGLE_GROUP"
      ))
    }
    
    # Create design matrix
    pred_matrix <- model.matrix(~ group, data = cell_metadata)
    
    # Run NEBULA
    nebula_res <- nebula(
      count = counts_mat,
      id = cell_metadata$individual,
      pred = pred_matrix,
      method = "LN",
      cpc = 0.05,
      mincp = 0.01,
      ncore = 1
    )
    
    # Check if results are valid
    if(is.null(nebula_res) || is.null(nebula_res$summary)) {
      warning("NEBULA returned NULL results")
      return(data.frame(
        n_genes_nebula = NA,
        n_sig_nebula_005 = NA,
        n_sig_nebula_001 = NA,
        n_sig_nebula_raw_005 = NA,
        mean_logfc_nebula = NA,
        power_nebula_005 = NA,
        power_nebula_001 = NA,
        nebula_convergence = NA,
        nebula_method = "NULL_RESULTS"
      ))
    }
    
    # Extract results
    results <- nebula_res$summary
    
    # Get the coefficient name
    coef_names <- colnames(results)
    logfc_col <- grep("logFC_group", coef_names, value = TRUE)[1]
    p_col <- grep("p_group", coef_names, value = TRUE)[1]
    
    if(is.na(logfc_col) || is.na(p_col)) {
      logfc_col <- grep("logFC", coef_names, value = TRUE)[1]
      p_col <- grep("^p_", coef_names, value = TRUE)[1]
    }
    
    # Calculate summary statistics
    if(!is.na(p_col)) {
      p_values <- results[[p_col]]
      n_sig_005 <- sum(p_values < 0.05, na.rm = TRUE)
      n_sig_001 <- sum(p_values < 0.01, na.rm = TRUE)
      
      # Adjust p-values using BH method
      padj <- p.adjust(p_values, method = "BH")
      n_sig_adj_005 <- sum(padj < 0.05, na.rm = TRUE)
      n_sig_adj_001 <- sum(padj < 0.01, na.rm = TRUE)
    } else {
      n_sig_005 <- n_sig_001 <- n_sig_adj_005 <- n_sig_adj_001 <- NA
      padj <- rep(NA, nrow(results))
    }
    
    if(!is.na(logfc_col)) {
      mean_logfc <- mean(abs(results[[logfc_col]]), na.rm = TRUE)
    } else {
      mean_logfc <- NA
    }
    
    # Calculate power if we know disease genes
    power_nebula_005 <- NA
    power_nebula_001 <- NA
    
    # FIX: Use filtered_rowdata
    if("is_disease_gene" %in% colnames(filtered_rowdata)) {
      disease_genes_in_filtered <- filtered_rowdata$is_disease_gene
      
      if(sum(disease_genes_in_filtered, na.rm = TRUE) > 0) {
        disease_padj <- padj[disease_genes_in_filtered]
        power_nebula_005 <- mean(disease_padj < 0.05, na.rm = TRUE)
        power_nebula_001 <- mean(disease_padj < 0.01, na.rm = TRUE)
      }
    }
    
    # Get convergence info
    if(!is.null(nebula_res$convergence)) {
      convergence_rate <- mean(nebula_res$convergence == 0, na.rm = TRUE)
    } else {
      convergence_rate <- NA
    }
    
    return(data.frame(
      n_genes_nebula = nrow(results),
      n_sig_nebula_005 = n_sig_adj_005,
      n_sig_nebula_001 = n_sig_adj_001,
      n_sig_nebula_raw_005 = n_sig_005,
      mean_logfc_nebula = mean_logfc,
      power_nebula_005 = power_nebula_005,
      power_nebula_001 = power_nebula_001,
      nebula_convergence = convergence_rate,
      nebula_method = "LN"
    ))
    
  }, error = function(e) {
    warning("NEBULA failed: ", e$message)
    # ... keep your fallback code ...
  })
}

# Alternative: SIMPLER VERSION using MAST (another mixed model package)
fit_mast_model <- function(sim, params = NULL) {
  tryCatch({
    library(MAST)
    
    # Prepare data
    counts_mat <- as.matrix(counts(sim))
    
    # Filter
    keep <- rowSums(counts_mat > 1) >= 10
    counts_mat <- counts_mat[keep, ]
    
    # Log transform
    log_counts <- log1p(counts_mat)
    
    # Create MAST object
    cData <- data.frame(
      wellKey = colnames(counts_mat),
      group = colData(sim)$Group,
      individual = colData(sim)$Individual,
      row.names = colnames(counts_mat)
    )
    
    fData <- data.frame(
      primerid = rownames(counts_mat),
      row.names = rownames(counts_mat)
    )
    
    sca <- FromMatrix(log_counts, cData, fData)
    
    # Fit model
    zlm_output <- zlm(~ group + (1|individual), sca, method = "glmer")
    
    # Get results
    summary_zlm <- summary(zlm_output, doLRT = "groupPatient")
    
    # Extract p-values
    p_values <- summary_zlm$datatable[summary_zlm$datatable$component == "H",]$`Pr(>Chisq)`
    
    # Calculate statistics
    n_sig_005 <- sum(p_values < 0.05, na.rm = TRUE)
    n_sig_001 <- sum(p_values < 0.01, na.rm = TRUE)
    
    padj <- p.adjust(p_values, method = "BH")
    n_sig_adj_005 <- sum(padj < 0.05, na.rm = TRUE)
    n_sig_adj_001 <- sum(padj < 0.01, na.rm = TRUE)
    
    return(data.frame(
      n_genes_nebula = length(p_values),
      n_sig_nebula_005 = n_sig_adj_005,
      n_sig_nebula_001 = n_sig_adj_001,
      n_sig_nebula_raw_005 = n_sig_005,
      mean_logfc_nebula = NA,
      power_nebula_005 = NA,
      power_nebula_001 = NA,
      nebula_convergence = NA,
      nebula_method = "MAST"
    ))
    
  }, error = function(e) {
    warning("MAST failed: ", e$message)
    return(data.frame(
      n_genes_nebula = NA,
      n_sig_nebula_005 = NA,
      n_sig_nebula_001 = NA,
      n_sig_nebula_raw_005 = NA,
      mean_logfc_nebula = NA,
      power_nebula_005 = NA,
      power_nebula_001 = NA,
      nebula_convergence = NA,
      nebula_method = "MAST_FAILED"
    ))
  })
}

# Combined fitting function - all three methods
fit_all_methods <- function(sim, params = NULL) {
  de_results <- fit_default_de_model(sim, params)
  pb_results <- fit_pseudobulk_model(sim, params)
  nebula_results <- fit_nebula_model(sim, params)
  cbind(de_results, pb_results, nebula_results)
}

#############################################
# SECTION 4: MAIN PIPELINE
#############################################

run_method_comparison_pipeline <- function(
    param_grid = NULL,
    n_simulations = 10,
    experiment_name = NULL,
    output_dir = "results/method_comparison",
    batch_size = 50,
    parallel = TRUE,
    n_cores = 12,
    verbose = TRUE,
    seed_start = 123) {
  
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate experiment name if not provided
  if(is.null(experiment_name)) {
    experiment_name <- paste0("method_comp_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # Create default parameter grid if not provided
  if(is.null(param_grid)) {
    param_grid <- expand.grid(
      n_controls = 10,
      n_patients = 10,
      cells_per_ind = c(50, 150, 300),
      cells_sd = c(0, 45),
      cells_min = 30,
      cells_max = NA,
      disease_effect_size = c(2, 5, 10),
      individual_variation = c(0.1, 0.3),
      disease_gene_fraction = 0.2,
      stringsAsFactors = FALSE
    )
  }
  
  # Ensure all required columns exist
  if(!"cells_sd" %in% names(param_grid)) param_grid$cells_sd <- 0
  if(!"cells_min" %in% names(param_grid)) param_grid$cells_min <- 30
  if(!"cells_max" %in% names(param_grid)) param_grid$cells_max <- NA
  
  # Add simulation replicates
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
    cat("=== Method Comparison Pipeline (3 Methods) ===\n")
    cat("Methods: Standard DE, Pseudobulk, NEBULA\n")
    cat("Experiment:", experiment_name, "\n")
    cat("Total simulations:", total_sims, "\n")
    cat("Parameter combinations:", max(param_grid$param_set_id), "\n")
    cat("Replications per combination:", n_simulations, "\n")
    cat("Output directory:", output_dir, "\n")
    cat("Parallel processing:", parallel, "\n\n")
    
    # Show parameter ranges
    cat("Parameter ranges:\n")
    cat("  Cell counts:", paste(unique(param_grid$cells_per_ind), collapse=", "), "\n")
    cat("  Cell SD:", paste(unique(param_grid$cells_sd), collapse=", "), "\n")
    cat("  Disease effect sizes:", paste(unique(param_grid$disease_effect_size), collapse=", "), "\n")
    cat("  Individual variation:", paste(unique(param_grid$individual_variation), collapse=", "), "\n\n")
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
        cells_sd = params$cells_sd,
        cells_min = params$cells_min,
        cells_max = if(is.na(params$cells_max)) NULL else params$cells_max,
        disease_effect_size = params$disease_effect_size,
        individual_variation = params$individual_variation,
        disease_gene_fraction = params$disease_gene_fraction,
        verbose = FALSE,
        seed = params$seed
      )
      
      # Fit all three models
      model_results <- fit_all_methods(sim, params)
      
      # Combine parameters and results
      results <- cbind(params, model_results)
      
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
  
  # Process simulations
  all_results <- list()
  
  for(batch in 1:n_batches) {
    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, total_sims)
    batch_indices <- batch_start:batch_end
    
    if(verbose) {
      cat(sprintf("Processing batch %d/%d (simulations %d-%d)\n",
                  batch, n_batches, batch_start, batch_end))
    }
    
    if(parallel && length(batch_indices) > 1 && !is.null(n_cores)) {
      # Parallel processing
      cl <- parallel::makeCluster(min(n_cores, length(batch_indices)))
      
      # Export everything needed
      parallel::clusterExport(cl, 
                              c("simulate_patient_control_study", 
                                "add_cell_state_correlation",
                                "fit_default_de_model",
                                "fit_pseudobulk_model",
                                "fit_nebula_model",
                                "fit_all_methods",
                                "param_grid"), 
                              envir = environment())
      
      parallel::clusterEvalQ(cl, {
        library(splatter)
        library(SingleCellExperiment)
        library(limma)
        library(edgeR)
        library(nebula)
        library(dplyr)
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
    
    # Save batch results
    batch_file <- file.path(output_dir, sprintf("batch_%03d.rds", batch))
    saveRDS(batch_df, batch_file)
    
    gc(verbose = FALSE)
  }
  
  # Combine all results
  final_results <- do.call(rbind, all_results)
  
  # Save final results
  saveRDS(final_results, file.path(output_dir, "combined_results.rds"))
  write.csv(final_results, file.path(output_dir, "combined_results.csv"), 
            row.names = FALSE)
  
  if(verbose) {
    cat("\n=== Pipeline Complete ===\n")
    cat("Results saved to:", output_dir, "\n")
    
    # Check for errors
    if("error" %in% colnames(final_results)) {
      n_errors <- sum(!is.na(final_results$error))
      if(n_errors > 0) {
        cat("WARNING:", n_errors, "simulations had errors\n")
      }
    }
    
    # Check NEBULA convergence
    if("nebula_convergence" %in% colnames(final_results)) {
      mean_conv <- mean(final_results$nebula_convergence, na.rm = TRUE)
      cat("NEBULA mean convergence rate:", round(mean_conv, 3), "\n")
    }
  }
  
  return(final_results)
}
#############################################
# SECTION 5: ENHANCED ANALYSIS FUNCTIONS
#############################################

analyze_results <- function(results_df) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  
  # Load if file path
  if(is.character(results_df)) {
    if(dir.exists(results_df)) {
      results_df <- readRDS(file.path(results_df, "combined_results.rds"))
    } else {
      results_df <- readRDS(results_df)
    }
  }
  
  # Remove error rows if any
  if("error" %in% colnames(results_df)) {
    results_df <- results_df[is.na(results_df$error), ]
  }
  
  # Summarize across replicates - now for 3 methods
  summary_stats <- results_df %>%
    group_by(disease_effect_size, individual_variation, cells_per_ind, cells_sd) %>%
    summarise(
      # Standard DE
      mean_power_de = mean(power_005, na.rm = TRUE),
      se_power_de = sd(power_005, na.rm = TRUE) / sqrt(n()),
      mean_n_sig_de = mean(n_sig_005, na.rm = TRUE),
      
      # Pseudobulk
      mean_power_pb = mean(power_pb_005, na.rm = TRUE),
      se_power_pb = sd(power_pb_005, na.rm = TRUE) / sqrt(n()),
      mean_n_sig_pb = mean(n_sig_pb_005, na.rm = TRUE),
      
      # NEBULA
      mean_power_nebula = mean(power_nebula_005, na.rm = TRUE),
      se_power_nebula = sd(power_nebula_005, na.rm = TRUE) / sqrt(n()),
      mean_n_sig_nebula = mean(n_sig_nebula_005, na.rm = TRUE),
      mean_convergence = mean(nebula_convergence, na.rm = TRUE),
      
      # Cell count info
      actual_mean_cells = mean(mean_cells_per_ind, na.rm = TRUE),
      actual_sd_cells = mean(sd_cells_per_ind, na.rm = TRUE),
      
      n_reps = n(),
      .groups = "drop"
    )
  
  # Plot 1: Three-way power comparison
  power_long <- summary_stats %>%
    pivot_longer(
      cols = c(mean_power_de, mean_power_pb, mean_power_nebula),
      names_to = "method",
      values_to = "power",
      names_prefix = "mean_power_"
    ) %>%
    mutate(
      se = case_when(
        method == "de" ~ se_power_de,
        method == "pb" ~ se_power_pb,
        method == "nebula" ~ se_power_nebula
      ),
      method = case_when(
        method == "de" ~ "Standard DE",
        method == "pb" ~ "Pseudobulk",
        method == "nebula" ~ "NEBULA"
      ),
      cell_variation = ifelse(cells_sd == 0, "No variation", "With variation")
    )
  
  p1 <- ggplot(power_long, 
               aes(x = disease_effect_size, y = power, 
                   color = method, linetype = factor(individual_variation))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = power - se, ymax = power + se), 
                  width = 0.2, alpha = 0.5) +
    facet_grid(cell_variation ~ cells_per_ind,
               labeller = labeller(cells_per_ind = function(x) paste(x, "cells/person"))) +
    theme_minimal() +
    labs(
      title = "Detection Power: Three-Method Comparison",
      x = "Disease Effect Size",
      y = "Power (FDR < 0.05)",
      color = "Method",
      linetype = "Individual\nVariation"
    ) +
    scale_color_manual(values = c("Standard DE" = "#E41A1C", 
                                  "Pseudobulk" = "#377EB8",
                                  "NEBULA" = "#4DAF4A")) +
    scale_y_continuous(limits = c(0, 1))
  
  # Plot 2: Number of discoveries - all three methods
  discoveries_long <- summary_stats %>%
    pivot_longer(
      cols = c(mean_n_sig_de, mean_n_sig_pb, mean_n_sig_nebula),
      names_to = "method",
      values_to = "n_discoveries",
      names_prefix = "mean_n_sig_"
    ) %>%
    mutate(
      method = case_when(
        method == "de" ~ "Standard DE",
        method == "pb" ~ "Pseudobulk",
        method == "nebula" ~ "NEBULA"
      )
    )
  
  p2 <- ggplot(discoveries_long,
               aes(x = factor(cells_per_ind), y = n_discoveries, 
                   fill = method, alpha = factor(cells_sd))) +
    geom_boxplot(position = position_dodge(width = 0.9)) +
    facet_wrap(~disease_effect_size,
               labeller = labeller(disease_effect_size = function(x) paste("Effect size:", x))) +
    theme_minimal() +
    scale_alpha_manual(values = c("0" = 1, "45" = 0.6),
                       labels = c("0" = "Equal cells", "45" = "Variable cells")) +
    labs(
      title = "Number of Discoveries by Method",
      x = "Average Cells per Individual",
      y = "Number of Significant Genes",
      fill = "Method",
      alpha = "Cell Count"
    ) +
    scale_fill_manual(values = c("Standard DE" = "#E41A1C", 
                                 "Pseudobulk" = "#377EB8",
                                 "NEBULA" = "#4DAF4A"))
  
  # Plot 3: Method agreement heatmap
  if(all(c("n_sig_005", "n_sig_pb_005", "n_sig_nebula_005") %in% names(results_df))) {
    cor_matrix <- results_df %>%
      select(n_sig_005, n_sig_pb_005, n_sig_nebula_005) %>%
      cor(use = "complete.obs")
    
    colnames(cor_matrix) <- c("Standard DE", "Pseudobulk", "NEBULA")
    rownames(cor_matrix) <- c("Standard DE", "Pseudobulk", "NEBULA")
    
    # Convert to long format for ggplot
    cor_long <- as.data.frame(as.table(cor_matrix))
    names(cor_long) <- c("Method1", "Method2", "Correlation")
    
    p3 <- ggplot(cor_long, aes(x = Method1, y = Method2, fill = Correlation)) +
      geom_tile() +
      geom_text(aes(label = round(Correlation, 2)), color = "white", size = 5) +
      scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#A50026",
                           midpoint = 0.5, limits = c(0, 1)) +
      theme_minimal() +
      labs(title = "Method Correlation",
           subtitle = "Based on number of significant genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p3 <- NULL
  }
  
  # Plot 4: NEBULA-specific diagnostics
  p4 <- ggplot(summary_stats,
               aes(x = factor(cells_per_ind), y = mean_convergence,
                   fill = factor(individual_variation))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~cells_sd,
               labeller = labeller(cells_sd = function(x) paste("Cell SD:", x))) +
    theme_minimal() +
    labs(
      title = "NEBULA Convergence Rate",
      x = "Cells per Individual",
      y = "Convergence Rate",
      fill = "Individual\nVariation"
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(0, 1))
  
  # ENHANCED PARAMETER EFFECT PLOTS - INCLUDING NEBULA
  
  # Plot 5: Effect of Disease Strength - All 3 Methods
  effect_impact <- results_df %>%
    pivot_longer(cols = c(power_005, power_pb_005, power_nebula_005),
                 names_to = "method",
                 values_to = "power") %>%
    mutate(method = case_when(
      method == "power_005" ~ "Standard DE",
      method == "power_pb_005" ~ "Pseudobulk",
      method == "power_nebula_005" ~ "NEBULA"
    ))
  
  p5 <- ggplot(effect_impact, 
               aes(x = factor(disease_effect_size), y = power, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    theme_minimal() +
    labs(title = "Impact of Disease Effect Size on All Methods",
         x = "Disease Effect Size (Fold Change)",
         y = "Detection Power",
         fill = "Method") +
    scale_fill_manual(values = c("Standard DE" = "#E41A1C", 
                                 "Pseudobulk" = "#377EB8",
                                 "NEBULA" = "#4DAF4A")) +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, 1))
  
  # Plot 6: Effect of Individual Variation - All 3 Methods
  p6 <- ggplot(effect_impact, 
               aes(x = factor(individual_variation), y = power, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    theme_minimal() +
    labs(title = "Impact of Individual Variation on All Methods",
         x = "Individual Variation (SD)",
         y = "Detection Power",
         fill = "Method") +
    scale_fill_manual(values = c("Standard DE" = "#E41A1C", 
                                 "Pseudobulk" = "#377EB8",
                                 "NEBULA" = "#4DAF4A")) +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, 1))
  
  # Plot 7: Cell Count Effects - All 3 Methods
  cell_impact <- results_df %>%
    pivot_longer(cols = c(power_005, power_pb_005, power_nebula_005),
                 names_to = "method",
                 values_to = "power") %>%
    mutate(
      method = case_when(
        method == "power_005" ~ "Standard DE",
        method == "power_pb_005" ~ "Pseudobulk",
        method == "power_nebula_005" ~ "NEBULA"
      ),
      cell_variation = ifelse(cells_sd == 0, "Fixed", "Variable")
    )
  
  p7 <- ggplot(cell_impact, 
               aes(x = factor(cells_per_ind), y = power, 
                   fill = method, alpha = cell_variation)) +
    geom_boxplot(position = position_dodge(width = 0.9), outlier.alpha = 0.3) +
    facet_wrap(~disease_effect_size,
               labeller = labeller(disease_effect_size = function(x) paste("Effect:", x))) +
    scale_alpha_manual(values = c("Fixed" = 1, "Variable" = 0.6),
                       name = "Cell counts") +
    scale_fill_manual(values = c("Standard DE" = "#E41A1C", 
                                 "Pseudobulk" = "#377EB8",
                                 "NEBULA" = "#4DAF4A")) +
    labs(title = "Impact of Cell Count and Variability on Method Performance",
         x = "Average cells per individual",
         y = "Detection power",
         fill = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, 1))
  
  # TYPE I ERROR ANALYSIS
  # Identify null simulations (disease_effect_size == 0 if you have them, 
  # or use disease_gene_fraction == 0)
  
  # For Type I error, we need to look at non-disease genes
  # Create a Type I error dataset
  type1_data <- results_df %>%
    mutate(
      # Calculate false positives for each method
      # Assuming disease_gene_fraction tells us what proportion are disease genes
      n_null_genes = n_genes_tested * (1 - disease_gene_fraction),
      n_null_genes_pb = n_genes_pb * (1 - disease_gene_fraction),
      n_null_genes_nebula = n_genes_nebula * (1 - disease_gene_fraction),
      
      # Expected false positives at alpha = 0.05
      expected_fp = n_null_genes * 0.05,
      expected_fp_pb = n_null_genes_pb * 0.05,
      expected_fp_nebula = n_null_genes_nebula * 0.05,
      
      # Observed false positive rate (approximation)
      # This is the rate among all genes minus the true positives
      fpr_de = pmax(0, (n_sig_005 - (disease_gene_fraction * n_genes_tested * power_005)) / n_null_genes),
      fpr_pb = pmax(0, (n_sig_pb_005 - (disease_gene_fraction * n_genes_pb * power_pb_005)) / n_null_genes_pb),
      fpr_nebula = pmax(0, (n_sig_nebula_005 - (disease_gene_fraction * n_genes_nebula * power_nebula_005)) / n_null_genes_nebula)
    )
  
  # Plot 8: Type I Error Rate by Method
  type1_long <- type1_data %>%
    pivot_longer(cols = c(fpr_de, fpr_pb, fpr_nebula),
                 names_to = "method",
                 values_to = "fpr",
                 names_prefix = "fpr_") %>%
    mutate(method = case_when(
      method == "de" ~ "Standard DE",
      method == "pb" ~ "Pseudobulk",
      method == "nebula" ~ "NEBULA"
    ))
  
  p8 <- ggplot(type1_long, 
               aes(x = method, y = fpr, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = 1.5, y = 0.055, label = "Nominal α = 0.05", 
             color = "red", size = 3.5) +
    theme_minimal() +
    labs(title = "Type I Error Rate by Method",
         subtitle = "False positive rate among non-disease genes",
         x = "Method",
         y = "False Positive Rate",
         fill = "Method") +
    scale_fill_manual(values = c("Standard DE" = "#E41A1C", 
                                 "Pseudobulk" = "#377EB8",
                                 "NEBULA" = "#4DAF4A")) +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, max(0.1, max(type1_long$fpr, na.rm = TRUE))))
  
  # Plot 9: Type I Error by Cell Count
  p9 <- ggplot(type1_long, 
               aes(x = factor(cells_per_ind), y = fpr, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    facet_wrap(~cells_sd,
               labeller = labeller(cells_sd = function(x) 
                 paste("Cell SD:", x))) +
    theme_minimal() +
    labs(title = "Type I Error Rate by Cell Count",
         x = "Cells per Individual",
         y = "False Positive Rate",
         fill = "Method") +
    scale_fill_manual(values = c("Standard DE" = "#E41A1C", 
                                 "Pseudobulk" = "#377EB8",
                                 "NEBULA" = "#4DAF4A")) +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, max(0.1, max(type1_long$fpr, na.rm = TRUE))))
  
  # Plot 10: Type I Error by Individual Variation
  p10 <- ggplot(type1_long, 
                aes(x = factor(individual_variation), y = fpr, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Type I Error Rate by Individual Variation",
         x = "Individual Variation (SD)",
         y = "False Positive Rate",
         fill = "Method") +
    scale_fill_manual(values = c("Standard DE" = "#E41A1C", 
                                 "Pseudobulk" = "#377EB8",
                                 "NEBULA" = "#4DAF4A")) +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, max(0.1, max(type1_long$fpr, na.rm = TRUE))))
  
  # Summary statistics for Type I error
  type1_summary <- type1_long %>%
    group_by(method) %>%
    summarise(
      mean_fpr = mean(fpr, na.rm = TRUE),
      sd_fpr = sd(fpr, na.rm = TRUE),
      median_fpr = median(fpr, na.rm = TRUE),
      q25_fpr = quantile(fpr, 0.25, na.rm = TRUE),
      q75_fpr = quantile(fpr, 0.75, na.rm = TRUE),
      prop_above_nominal = mean(fpr > 0.05, na.rm = TRUE),
      max_fpr = max(fpr, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Combine parameter effect plots
  param_effects <- (p5 | p6) / p7
  
  # Combine Type I error plots
  type1_plots <- (p8 | p9) / p10
  
  # Overall combined plot
  combined <- (p1 / p2) | (p3 / p4)
  
  return(list(
    summary = summary_stats,
    type1_summary = type1_summary,
    type1_data = type1_data,
    plots = list(
      power = p1,
      discoveries = p2,
      correlation = p3,
      nebula_diagnostics = p4,
      effect_size = p5,
      individual_variation = p6,
      cell_effects = p7,
      type1_overall = p8,
      type1_by_cells = p9,
      type1_by_variation = p10,
      combined = combined,
      param_effects = param_effects,
      type1_combined = type1_plots
    ),
    raw_data = results_df
  ))
}

#############################################
# SECTION 6: RUN ANALYSIS WITH NEW PLOTS
#############################################

# Set parameters including cell count variation
param_grid <- expand.grid(
  n_controls = 10,
  n_patients = 10,
  cells_per_ind = c(50, 150, 300),      # Low, medium, high cell counts
  cells_sd = c(0, 45),                   # No variation vs ~30% CV
  cells_min = 30,                        # Minimum 30 cells
  cells_max = NA,                        # No maximum
  disease_effect_size = c(2, 5, 10),    # Weak, moderate, strong effects
  individual_variation = c(0.1, 0.3),    # Low vs high individual variation
  disease_gene_fraction = c(0.05, 0.01, 0.2),           # 20% of genes affected
  stringsAsFactors = FALSE
)

cat("=== Three-Method Comparison Analysis ===\n")
cat("Methods: Standard DE, Pseudobulk, NEBULA\n")
cat("Total parameter combinations:", nrow(param_grid), "\n")
cat("Cell count scenarios:", length(unique(param_grid$cells_per_ind)), "x", 
    length(unique(param_grid$cells_sd)), "\n")
cat("Biological scenarios:", length(unique(param_grid$disease_effect_size)), "x",
    length(unique(param_grid$individual_variation)), "\n\n")

# Run pipeline
cat("Starting three-method comparison analysis...\n\n")

results <- run_method_comparison_pipeline(
  param_grid = param_grid,
  n_simulations = 10,  # 10 replicates per parameter combination
  experiment_name = "three_method_comparison",
  output_dir = "results/three_method_comparison",
  batch_size = 50,
  parallel = TRUE,  # Set to TRUE if you have multiple cores
  n_cores = 12,       # Adjust based on your system
  verbose = TRUE
)


# Analyze results
cat("\nAnalyzing results...\n")
analysis <- analyze_results(results)

# Display summaries
cat("\n=== Power Summary ===\n")
print(analysis$summary)

cat("\n=== Type I Error Summary ===\n")
print(analysis$type1_summary)

# Show plots
print(analysis$plots$combined)
print(analysis$plots$param_effects)
print(analysis$plots$type1_combined)

# Save all plots
ggsave("results/three_method_comparison/combined_plot.png", 
       analysis$plots$combined, 
       width = 16, height = 12, dpi = 300)

ggsave("results/three_method_comparison/parameter_effects.png", 
       analysis$plots$param_effects, 
       width = 16, height = 12, dpi = 300)

ggsave("results/three_method_comparison/type1_error_analysis.png", 
       analysis$plots$type1_combined, 
       width = 16, height = 12, dpi = 300)

# Individual parameter effect plots
ggsave("results/three_method_comparison/effect_size_impact.png", 
       analysis$plots$effect_size, 
       width = 10, height = 6, dpi = 300)

ggsave("results/three_method_comparison/individual_variation_impact.png", 
       analysis$plots$individual_variation, 
       width = 10, height = 6, dpi = 300)

ggsave("results/three_method_comparison/cell_count_impact.png", 
       analysis$plots$cell_effects, 
       width = 12, height = 8, dpi = 300)

cat("\n=== All plots saved ===\n")
#############################################
# SECTION 5B: REPORT GENERATION FUNCTIONS
#############################################

# Function to create HTML report
create_html_report <- function(analysis, results_df, output_dir, 
                               report_title = "Method Comparison Report") {
  
  library(rmarkdown)
  library(knitr)
  library(DT)
  
  # Create report content
  report_content <- '---
title: "`r report_title`"
author: "Method Comparison Analysis"
date: "`r Sys.Date()`"
output: 
  html_document:
    toc: true
    toc_float: true
    toc_depth: 3
    theme: united
    highlight: tango
    code_folding: hide
    df_print: paged
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.width = 10, fig.height = 6)
library(ggplot2)
library(dplyr)
library(knitr)
library(DT)
library(patchwork)
```

# Executive Summary

This report compares **Standard Differential Expression** and **Pseudobulk** methods across different simulation parameters.
```{r summary_stats}
# Overall performance summary
overall_summary <- results_df %>%
  filter(!is.na(power_005)) %>%
  summarise(
    `Total Simulations` = n(),
    `Parameter Combinations` = n_distinct(param_set_id),
    `Mean DE Power` = round(mean(power_005, na.rm = TRUE), 3),
    `Mean PB Power` = round(mean(power_pb_005, na.rm = TRUE), 3),
    `Mean DE Discoveries` = round(mean(n_sig_005, na.rm = TRUE), 0),
    `Mean PB Discoveries` = round(mean(n_sig_pb_005, na.rm = TRUE), 0)
  )

kable(t(overall_summary), col.names = "Value", caption = "Overall Summary Statistics")
```

## Key Findings
```{r key_findings}
# Which method performs better?
power_comparison <- results_df %>%
  group_by(disease_effect_size, individual_variation) %>%
  summarise(
    de_wins = sum(power_005 > power_pb_005, na.rm = TRUE),
    pb_wins = sum(power_pb_005 > power_005, na.rm = TRUE),
    ties = sum(power_005 == power_pb_005, na.rm = TRUE),
    .groups = "drop"
  )

# Best performing conditions
best_conditions <- analysis$summary %>%
  mutate(
    best_de = paste0(round(mean_power_de * 100, 1), "%"),
    best_pb = paste0(round(mean_power_pb * 100, 1), "%"),
    winner = ifelse(mean_power_de > mean_power_pb, "Standard DE", "Pseudobulk")
  )

kable(best_conditions %>% 
      select(disease_effect_size, individual_variation, best_de, best_pb, winner),
      col.names = c("Effect Size", "Individual Variation", "DE Power", "PB Power", "Better Method"),
      caption = "Power Comparison Across Conditions")
```

# Detailed Results

## Power Analysis

### Detection Power by Effect Size
```{r power_plot, fig.width=12, fig.height=7}
print(analysis$plots$power + 
      theme(legend.position = "bottom") +
      labs(subtitle = "Higher values indicate better detection of true disease genes"))
```

### Statistical Comparison
```{r statistical_tests}
# Paired t-tests for each condition
stat_tests <- results_df %>%
  group_by(disease_effect_size, individual_variation) %>%
  summarise(
    n = n(),
    mean_diff = mean(power_005 - power_pb_005, na.rm = TRUE),
    t_test_p = ifelse(n() > 1, 
                      t.test(power_005, power_pb_005, paired = TRUE)$p.value, 
                      NA),
    .groups = "drop"
  ) %>%
  mutate(
    p_adjusted = p.adjust(t_test_p, method = "fdr"),
    significance = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01 ~ "**",
      p_adjusted < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

kable(stat_tests,
      col.names = c("Effect Size", "Ind. Variation", "N", "Mean Difference", "P-value", "Adj. P-value", "Significance"),
      digits = 4,
      caption = "Paired t-tests comparing methods (positive difference favors Standard DE)")
```

## Discovery Analysis

### Number of Significant Genes
```{r discoveries_plot, fig.width=12, fig.height=6}
print(analysis$plots$discoveries +
      theme(legend.position = "bottom") +
      labs(subtitle = "Average number of genes detected as significant (FDR < 0.05)"))
```

## Method Agreement

### Correlation Between Methods
```{r correlation_plot, fig.width=12, fig.height=8}
print(analysis$plots$correlation +
      labs(subtitle = "Each point represents one simulation"))
```

### Agreement Statistics
```{r agreement_stats}
agreement <- results_df %>%
  mutate(
    both_detect = (n_sig_005 > 0) & (n_sig_pb_005 > 0),
    only_de = (n_sig_005 > 0) & (n_sig_pb_005 == 0),
    only_pb = (n_sig_pb_005 > 0) & (n_sig_005 == 0),
    neither = (n_sig_005 == 0) & (n_sig_pb_005 == 0)
  ) %>%
  group_by(disease_effect_size, individual_variation) %>%
  summarise(
    `Both Detect` = paste0(round(mean(both_detect) * 100, 1), "%"),
    `Only DE` = paste0(round(mean(only_de) * 100, 1), "%"),
    `Only PB` = paste0(round(mean(only_pb) * 100, 1), "%"),
    `Neither` = paste0(round(mean(neither) * 100, 1), "%"),
    Correlation = round(cor(n_sig_005, n_sig_pb_005), 3),
    .groups = "drop"
  )

kable(agreement,
      caption = "Method Agreement Across Conditions")
```

# Parameter Effects

## Effect of Disease Strength
```{r effect_size_impact, fig.width=10, fig.height=6}
effect_impact <- results_df %>%
  pivot_longer(cols = c(power_005, power_pb_005),
               names_to = "method",
               values_to = "power") %>%
  mutate(method = ifelse(method == "power_005", "Standard DE", "Pseudobulk"))

p_effect <- ggplot(effect_impact, 
                   aes(x = factor(disease_effect_size), y = power, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Impact of Disease Effect Size",
       x = "Disease Effect Size (Fold Change)",
       y = "Detection Power",
       fill = "Method") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

print(p_effect)
```

## Effect of Individual Variation
```{r individual_variation_impact, fig.width=10, fig.height=6}
p_var <- ggplot(effect_impact, 
                aes(x = factor(individual_variation), y = power, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Impact of Individual Variation",
       x = "Individual Variation (SD)",
       y = "Detection Power",
       fill = "Method") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

print(p_var)
```

## Cell variation analysis
```{r cell_table}
# Analyze impact of cell counts
cell_impact <- results %>%
  group_by(cells_per_ind, cells_sd, disease_effect_size) %>%
  summarise(
    de_power = mean(power_005, na.rm = TRUE),
    pb_power = mean(power_pb_005, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(de_power, pb_power),
               names_to = "method",
               values_to = "power")

p_cell <- ggplot(cell_impact, 
            aes(x = factor(cells_per_ind), y = power, 
                fill = method, alpha = factor(cells_sd))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~disease_effect_size) +
  scale_alpha_manual(values = c("0" = 1, "45" = 0.6),
                     labels = c("0" = "No variation", "45" = "With variation")) +
  labs(title = "Impact of Cell Count on Method Performance",
       x = "Average cells per individual",
       y = "Detection power",
       fill = "Method",
       alpha = "Cell count variation") +
  theme_minimal()

print(p_cell)
```

# Detailed Data Tables

## Summary Statistics by Condition
```{r summary_table}
datatable(analysis$summary %>%
            select(-starts_with("se_")) %>%
            mutate(across(where(is.numeric), ~round(., 3))),
          options = list(pageLength = 15),
          caption = "Complete summary statistics for all parameter combinations")
```

## Simulation Parameters
```{r param_summary}
param_summary <- results_df %>%
  group_by(param_set_id) %>%
  slice(1) %>%
  select(param_set_id, n_controls, n_patients, cells_per_ind, 
         disease_effect_size, individual_variation, disease_gene_fraction) %>%
  ungroup()

datatable(param_summary,
          options = list(pageLength = 10),
          caption = "Parameter values for each parameter set")
```

# Recommendations
```{r recommendations}
# Find best conditions for each method
best_de <- analysis$summary %>%
  filter(mean_power_de == max(mean_power_de)) %>%
  slice(1)

best_pb <- analysis$summary %>%
  filter(mean_power_pb == max(mean_power_pb)) %>%
  slice(1)

best_overall <- analysis$summary %>%
  mutate(avg_power = (mean_power_de + mean_power_pb) / 2) %>%
  filter(avg_power == max(avg_power)) %>%
  slice(1)
```

Based on the analysis:

1. **Best conditions for Standard DE**: Effect size = `r best_de$disease_effect_size`, Individual variation = `r best_de$individual_variation` (Power: `r round(best_de$mean_power_de * 100, 1)`%)

2. **Best conditions for Pseudobulk**: Effect size = `r best_pb$disease_effect_size`, Individual variation = `r best_pb$individual_variation` (Power: `r round(best_pb$mean_power_pb * 100, 1)`%)

3. **Overall recommendation**: `r ifelse(mean(results_df$power_005, na.rm=TRUE) > mean(results_df$power_pb_005, na.rm=TRUE), "Standard DE performs better on average", "Pseudobulk performs better on average")` across the tested conditions.

4. **Method agreement**: The methods show `r ifelse(cor(results_df$n_sig_005, results_df$n_sig_pb_005, use="complete.obs") > 0.7, "high", ifelse(cor(results_df$n_sig_005, results_df$n_sig_pb_005, use="complete.obs") > 0.4, "moderate", "low"))` correlation (r = `r round(cor(results_df$n_sig_005, results_df$n_sig_pb_005, use="complete.obs"), 3)`).

# Session Information
```{r session_info}
sessionInfo()
```
'
  
  # Write the Rmd file
  rmd_file <- file.path(output_dir, "method_comparison_report.Rmd")
  writeLines(report_content, rmd_file)
  
  # Render the report
  output_file <- file.path(output_dir, paste0(report_title, "_", Sys.Date(), ".html"))
  
  render(rmd_file, 
         output_file = basename(output_file),
         output_dir = output_dir,
         params = list(
           analysis = analysis,
           results_df = results_df,
           report_title = report_title
         ),
         envir = new.env())
  
  cat("HTML report saved to:", output_file, "\n")
  
  # Clean up Rmd file
  file.remove(rmd_file)
  
  return(output_file)
}

#############################################
# SECTION 6: RUN ANALYSIS (UPDATED)
#############################################

# Analyze results
cat("\nAnalyzing results...\n")
analysis <- analyze_results(results)

# Display summary
print(analysis$summary)

# Show plots
print(analysis$plots$combined)

# Save plots
ggsave("results/method_comparison_final/combined_plot.png", 
       analysis$plots$combined, 
       width = 14, height = 10, dpi = 300)

# CREATE HTML REPORT
cat("\nGenerating HTML report...\n")
report_file <- create_html_report(
  analysis = analysis,
  results_df = results,
  output_dir = "results/method_comparison_final",
  report_title = "Method Comparison Analysis Report"
)

cat("\n=== Analysis Complete ===\n")
cat("Results saved to: results/method_comparison_final/\n")
cat("HTML report: ", report_file, "\n")

# Optionally open the report
if(interactive()) {
  browseURL(report_file)
}