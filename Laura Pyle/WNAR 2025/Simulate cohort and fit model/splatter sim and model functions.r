# Main function that combines simulation and model fitting
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
      disease_effect_size = c(2, 5, 10),
      individual_variation = c(0.1, 0.3),
      disease_gene_fraction = 0.2,
      stringsAsFactors = FALSE
    )
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
  
  # Function to process single simulation
  process_single_simulation <- function(i) {
    params <- param_grid[i, ]
    
    tryCatch({
      # Generate simulation
      sim <- simulate_patient_control_study(
        n_controls = params$n_controls,
        n_patients = params$n_patients,
        cells_per_ind = params$cells_per_ind,
        disease_effect_size = params$disease_effect_size,
        individual_variation = params$individual_variation,
        disease_gene_fraction = params$disease_gene_fraction,
        verbose = FALSE,
        seed = params$seed
      )
      
      # Add simulation ID to the object
      colData(sim)$simulation_id <- params$simulation_id
      
      # Fit model
      if(!is.null(fit_function)) {
        model_results <- fit_function(sim, params)
      } else {
        # Default: simple differential expression
        model_results <- fit_default_de_model(sim, params)
      }
      
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

# Example 1: Mixed effects model accounting for individual
fit_mixed_model <- function(sim, params) {
  library(lme4)
  library(lmerTest)
  
  # Get top variable genes
  counts_mat <- counts(sim)
  gene_vars <- apply(counts_mat, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:100])
  
  # Fit mixed model to top genes
  results_list <- list()
  
  for(gene in top_genes[1:10]) {  # Just top 10 for speed
    expr_data <- data.frame(
      expression = counts_mat[gene, ],
      group = colData(sim)$Group,
      individual = colData(sim)$Individual
    )
    
    # Fit mixed model
    model <- lmer(expression ~ group + (1|individual), data = expr_data)
    
    # Extract p-value
    sum_model <- summary(model)
    p_value <- sum_model$coefficients["groupPatient", "Pr(>|t|)"]
    
    results_list[[gene]] <- p_value
  }
  
  return(data.frame(
    n_genes_mixed = 10,
    mean_pvalue = mean(unlist(results_list), na.rm = TRUE),
    n_sig_mixed = sum(unlist(results_list) < 0.05, na.rm = TRUE)
  ))
}

# Example 2: Pseudobulk aggregation approach
fit_pseudobulk_model <- function(sim, params) {
  library(edgeR)
  library(dplyr)
  
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

# Example 3: Machine learning classifier
fit_ml_model <- function(sim, params) {
  library(ranger)  # Random forest
  
  # Prepare data
  expr_mat <- t(log1p(counts(sim)))  # Cells x genes
  
  # Select top variable genes
  gene_vars <- apply(expr_mat, 2, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:500])
  
  expr_mat <- expr_mat[, top_genes]
  
  # Create training data
  train_data <- data.frame(
    expr_mat,
    Group = colData(sim)$Group
  )
  
  # Fit random forest
  rf_model <- ranger(Group ~ ., 
                     data = train_data,
                     num.trees = 100,
                     importance = "impurity")
  
  # Out-of-bag accuracy
  oob_accuracy <- 1 - rf_model$prediction.error
  
  # Variable importance
  var_importance <- importance(rf_model)
  top_important <- sort(var_importance, decreasing = TRUE)[1:10]
  
  return(data.frame(
    rf_oob_accuracy = oob_accuracy,
    rf_num_trees = 100,
    top_gene_importance = mean(top_important)
  ))
}

# Load and analyze results
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
    geom_line(linewidth = 1.2) +
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