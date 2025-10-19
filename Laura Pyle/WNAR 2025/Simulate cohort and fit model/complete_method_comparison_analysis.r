# complete_method_comparison_analysis.R
# Complete script for running simulation-based method comparison
# All functions included and properly organized

#############################################
# SECTION 1: LOAD LIBRARIES
#############################################

library(splatter)
library(SingleCellExperiment)
library(limma)
library(edgeR)
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
  } else {
    cell_counts <- cells_per_ind
  }
  
  # Individual-specific expression variation
  ind_effects <- rnorm(n_individuals, 0, individual_variation)
  
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
  
  # DE analysis
  y <- DGEList(counts = pb_counts)
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ ind_groups$Group)
  y <- estimateDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  
  results <- topTags(lrt, n = Inf)$table
  
  # Calculate power if we know disease genes
  power_pb_005 <- NA
  power_pb_001 <- NA
  
  if("is_disease_gene" %in% colnames(rowData(sim))) {
    disease_genes <- rownames(sim)[rowData(sim)$is_disease_gene]
    disease_genes_tested <- intersect(disease_genes, rownames(results))
    
    if(length(disease_genes_tested) > 0) {
      disease_results <- results[disease_genes_tested, ]
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
    power_pb_001 = power_pb_001
  ))
}

# Combined fitting function
fit_multiple_methods <- function(sim, params = NULL) {
  de_results <- fit_default_de_model(sim, params)
  pb_results <- fit_pseudobulk_model(sim, params)
  cbind(de_results, pb_results)
}

#############################################
# SECTION 4: MAIN PIPELINE
#############################################

run_method_comparison_pipeline <- function(
    param_grid = NULL,
    n_simulations = 10,
    experiment_name = NULL,
    output_dir = "/Users/pylell/Library/CloudStorage/OneDrive-UW/Pyle/scRNAseq simulations/results/method_comparison",
    batch_size = 50,
    parallel = FALSE,
    n_cores = NULL,
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
      cells_per_ind = 150,
      disease_effect_size = c(2, 5, 10),
      individual_variation = c(0.1, 0.3),
      disease_gene_fraction = 0.2,
      stringsAsFactors = FALSE
    )
  }
  
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
    cat("=== Method Comparison Pipeline ===\n")
    cat("Experiment:", experiment_name, "\n")
    cat("Total simulations:", total_sims, "\n")
    cat("Parameter combinations:", max(param_grid$param_set_id), "\n")
    cat("Replications per combination:", n_simulations, "\n")
    cat("Output directory:", output_dir, "\n")
    cat("Parallel processing:", parallel, "\n\n")
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
      
      # Fit both models
      model_results <- fit_multiple_methods(sim, params)
      
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
                                "fit_multiple_methods",
                                "param_grid"), 
                              envir = environment())
      
      parallel::clusterEvalQ(cl, {
        library(splatter)
        library(SingleCellExperiment)
        library(limma)
        library(edgeR)
        library(dplyr)
      })
      
      batch_results <- parallel::parLapply(cl, batch_indices, 
                                           process_single_simulation)
      parallel::stopCluster(cl)
      
    } else {
      # Sequential processing (more reliable)
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
  }
  
  return(final_results)
}

#############################################
# SECTION 5: ANALYSIS FUNCTIONS
#############################################

analyze_results <- function(results_df) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  
  results_df <- "/Users/pylell/Library/CloudStorage/OneDrive-UW/Pyle/scRNAseq simulations/results/method_comparison"
  
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
  
  # Summarize across replicates
  summary_stats <- results_df %>%
    group_by(disease_effect_size, individual_variation) %>%
    summarise(
      # Standard DE
      mean_power_de = mean(power_005, na.rm = TRUE),
      se_power_de = sd(power_005, na.rm = TRUE) / sqrt(n()),
      mean_n_sig_de = mean(n_sig_005, na.rm = TRUE),
      
      # Pseudobulk
      mean_power_pb = mean(power_pb_005, na.rm = TRUE),
      se_power_pb = sd(power_pb_005, na.rm = TRUE) / sqrt(n()),
      mean_n_sig_pb = mean(n_sig_pb_005, na.rm = TRUE),
      
      n_reps = n(),
      .groups = "drop"
    )
  
  # Reshape for plotting
  power_long <- summary_stats %>%
    pivot_longer(
      cols = c(mean_power_de, mean_power_pb),
      names_to = "method",
      values_to = "power",
      names_prefix = "mean_power_"
    ) %>%
    mutate(
      se = ifelse(method == "de", se_power_de, se_power_pb),
      method = ifelse(method == "de", "Standard DE", "Pseudobulk")
    )
  
  # Plot 1: Power comparison
  p1 <- ggplot(power_long, 
               aes(x = disease_effect_size, y = power, 
                   color = method, linetype = factor(individual_variation))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = power - se, ymax = power + se), 
                  width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "Detection Power Comparison",
      x = "Disease Effect Size",
      y = "Power (FDR < 0.05)",
      color = "Method",
      linetype = "Individual\nVariation"
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(limits = c(0, 1))
  
  # Plot 2: Number of discoveries
  discoveries_long <- summary_stats %>%
    pivot_longer(
      cols = c(mean_n_sig_de, mean_n_sig_pb),
      names_to = "method",
      values_to = "n_discoveries",
      names_prefix = "mean_n_sig_"
    ) %>%
    mutate(
      method = ifelse(method == "de", "Standard DE", "Pseudobulk")
    )
  
  p2 <- ggplot(discoveries_long,
               aes(x = factor(disease_effect_size), y = n_discoveries,
                   fill = method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~individual_variation, 
               labeller = labeller(individual_variation = function(x) paste("Ind. Var. =", x))) +
    theme_minimal() +
    labs(
      title = "Number of Discoveries",
      x = "Disease Effect Size",
      y = "Mean Significant Genes",
      fill = "Method"
    ) +
    scale_fill_brewer(palette = "Set2")
  
  # Plot 3: Method correlation
  p3 <- ggplot(results_df,
               aes(x = n_sig_005, y = n_sig_pb_005)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", color = "red") +
    facet_grid(individual_variation ~ disease_effect_size,
               labeller = labeller(
                 individual_variation = function(x) paste("Var:", x),
                 disease_effect_size = function(x) paste("Effect:", x)
               )) +
    theme_minimal() +
    labs(
      title = "Method Correlation",
      x = "Standard DE (# significant)",
      y = "Pseudobulk (# significant)"
    ) +
    coord_fixed()
  
  # Combine plots
  combined <- (p1 / p2) | p3
  
  return(list(
    summary = summary_stats,
    plots = list(
      power = p1,
      discoveries = p2,
      correlation = p3,
      combined = combined
    ),
    raw_data = results_df
  ))
}

#############################################
# SECTION 6: RUN ANALYSIS
#############################################

# Set parameters
param_grid <- expand.grid(
  n_controls = 10,
  n_patients = 10,
  cells_per_ind = 150,
  disease_effect_size = c(1, 2, 5, 10),
  individual_variation = c(0.1, 0.2, 0.3),
  disease_gene_fraction = 0.2,
  stringsAsFactors = FALSE
)

# Run pipeline
cat("Starting method comparison analysis...\n\n")

results <- run_method_comparison_pipeline(
  param_grid = param_grid,
  n_simulations = 20,  # 20 replicates per parameter combination
  experiment_name = "method_comparison_final",
  output_dir = "results/method_comparison_final",
  batch_size = 50,
  parallel = FALSE,  # Set to TRUE if you have multiple cores
  n_cores = 4,       # Adjust based on your system
  verbose = TRUE
)

# Analyze results
cat("\nAnalyzing results...\n")
analysis <- analyze_results(results)

# Display summary
print(analysis$summary)

# Show plots
print(analysis$plots$combined)

# Save plots
ggsave("/Users/pylell/Library/CloudStorage/OneDrive-UW/Pyle/scRNAseq simulations/results/method_comparison/combined_plot.png", 
       analysis$plots$combined, 
       width = 14, height = 10, dpi = 300)

cat("\n=== Analysis Complete ===\n")
cat("Results saved to: /Users/pylell/Library/CloudStorage/OneDrive-UW/Pyle/scRNAseq simulations/results/method_comparison/\n")