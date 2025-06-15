# Single Cell RNA-seq Data Simulation Using Splatter
# Comparison of pseudobulk vs mixed effects methods with cell imbalance analysis

# Load required libraries
library(splatter)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(Matrix)

# Check for optional libraries (will load when needed)
check_optional_packages <- function() {
  required_packages <- c("DESeq2", "nebula")
  missing_packages <- c()
  
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if(length(missing_packages) > 0) {
    cat("Optional packages needed for DE analysis:\n")
    for(pkg in missing_packages) {
      if(pkg == "DESeq2") {
        cat("  install DESeq2: BiocManager::install('DESeq2')\n")
      } else {
        cat("  install", pkg, ": install.packages('", pkg, "')\n")
      }
    }
    cat("Note: Simulations will work without these, but DE comparisons require them.\n\n")
  }
}

# Check packages on load
check_optional_packages()

# Function to simulate scRNA-seq data with cell type imbalance
simulate_scrna_imbalanced <- function(
    n_cells_per_condition = 500,
    cell_imbalance_ratio = c(1, 1),
    n_genes = 2000,
    de_prob = 0.1,
    de_facLoc = 1.0,
    de_facScale = 0.4,
    group_effect_size = 1.5,
    time_effect_size = 1.2,
    interaction_effect_size = 1.3,
    dropout_type = "experiment",
    seed = 123
) {
  
  set.seed(seed)
  
  # Calculate actual cell numbers based on imbalance ratio
  if(length(cell_imbalance_ratio) != 2) {
    stop("cell_imbalance_ratio must be a vector of length 2, e.g., c(3, 1)")
  }
  
  total_ratio <- sum(cell_imbalance_ratio)
  control_cells_t1 <- max(10, round(n_cells_per_condition * cell_imbalance_ratio[1] / total_ratio))
  control_cells_t2 <- control_cells_t1
  treatment_cells_t1 <- max(10, round(n_cells_per_condition * cell_imbalance_ratio[2] / total_ratio))
  treatment_cells_t2 <- treatment_cells_t1
  
  total_cells <- control_cells_t1 + control_cells_t2 + treatment_cells_t1 + treatment_cells_t2
  
  # Create group probabilities based on actual cell numbers
  group_prob <- c(
    control_cells_t1 / total_cells,
    control_cells_t2 / total_cells,
    treatment_cells_t1 / total_cells,
    treatment_cells_t2 / total_cells
  )
  
  # Estimate parameters from reference
  params <- newSplatParams()
  
  # Set parameters
  params <- setParams(params, 
                      nGenes = n_genes,
                      batchCells = total_cells,
                      group.prob = group_prob,
                      de.prob = de_prob,
                      de.facLoc = de_facLoc,
                      de.facScale = de_facScale,
                      dropout.type = dropout_type,
                      dropout.mid = 3,
                      dropout.shape = -1,
                      seed = seed)
  
  # Simulate the base data
  sim <- splatSimulate(params, method = "groups", verbose = FALSE)
  
  # Extract count matrix
  counts <- as.matrix(counts(sim))
  
  # Create condition labels
  condition_labels <- c(
    rep("Control_T1", control_cells_t1),
    rep("Control_T2", control_cells_t2),
    rep("Treatment_T1", treatment_cells_t1),
    rep("Treatment_T2", treatment_cells_t2)
  )
  
  # Create metadata
  cell_metadata <- data.frame(
    CellID = colnames(counts),
    Group = c(
      rep("Control", control_cells_t1 + control_cells_t2),
      rep("Treatment", treatment_cells_t1 + treatment_cells_t2)
    ),
    Time = c(
      rep("T1", control_cells_t1),
      rep("T2", control_cells_t2),
      rep("T1", treatment_cells_t1),
      rep("T2", treatment_cells_t2)
    ),
    Condition = condition_labels,
    stringsAsFactors = FALSE
  )
  
  # Apply custom effect sizes
  gene_metadata <- as.data.frame(rowData(sim))
  de_genes <- which(gene_metadata$DEFacGroup1 != 1 | 
                      gene_metadata$DEFacGroup2 != 1 |
                      gene_metadata$DEFacGroup3 != 1 |
                      gene_metadata$DEFacGroup4 != 1)
  
  modified_counts <- counts
  
  for(gene_idx in de_genes) {
    control_t1_idx <- which(cell_metadata$Condition == "Control_T1")
    control_t2_idx <- which(cell_metadata$Condition == "Control_T2") 
    treatment_t1_idx <- which(cell_metadata$Condition == "Treatment_T1")
    treatment_t2_idx <- which(cell_metadata$Condition == "Treatment_T2")
    
    # Apply effects
    modified_counts[gene_idx, control_t2_idx] <- 
      round(modified_counts[gene_idx, control_t2_idx] * time_effect_size)
    modified_counts[gene_idx, treatment_t2_idx] <- 
      round(modified_counts[gene_idx, treatment_t2_idx] * time_effect_size)
    
    modified_counts[gene_idx, treatment_t1_idx] <- 
      round(modified_counts[gene_idx, treatment_t1_idx] * group_effect_size)
    modified_counts[gene_idx, treatment_t2_idx] <- 
      round(modified_counts[gene_idx, treatment_t2_idx] * group_effect_size)
    
    modified_counts[gene_idx, treatment_t2_idx] <- 
      round(modified_counts[gene_idx, treatment_t2_idx] * interaction_effect_size)
  }
  
  # Calculate QC metrics
  cell_metadata$Total_UMI <- colSums(modified_counts)
  cell_metadata$Genes_Detected <- colSums(modified_counts > 0)
  
  # Enhanced gene metadata
  gene_metadata$Cells_Expressing <- rowSums(modified_counts > 0)
  gene_metadata$Mean_Expression <- rowMeans(modified_counts)
  gene_metadata$Is_DE <- rownames(gene_metadata) %in% rownames(counts)[de_genes]
  
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(counts = modified_counts),
    colData = cell_metadata,
    rowData = gene_metadata
  )
  
  # Add normalized data
  sce <- logNormCounts(sce)
  
  return(list(
    sce = sce,
    counts = modified_counts,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata,
    original_sim = sim,
    de_genes = de_genes,
    simulation_params = list(
      cell_imbalance_ratio = cell_imbalance_ratio,
      actual_cell_numbers = c(
        Control_T1 = control_cells_t1,
        Control_T2 = control_cells_t2,  
        Treatment_T1 = treatment_cells_t1,
        Treatment_T2 = treatment_cells_t2
      ),
      n_genes = n_genes,
      de_prob = de_prob,
      group_effect_size = group_effect_size,
      time_effect_size = time_effect_size,
      interaction_effect_size = interaction_effect_size
    )
  ))
}

# Wrapper function for backward compatibility
simulate_scrna_splatter <- function(...) {
  return(simulate_scrna_imbalanced(...))
}

# Function to visualize the simulated data
plot_splatter_simulation <- function(sim_data) {
  
  # Plot 1: UMI distribution
  p1 <- ggplot(sim_data$cell_metadata, aes(x = Condition, y = Total_UMI, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Total UMI per Cell", x = "Condition", y = "Total UMI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_summary(fun = length, geom = "text", aes(label = paste("n =", after_stat(y))), 
                 vjust = -0.5, size = 3)
  
  # Plot 2: Genes detected
  p2 <- ggplot(sim_data$cell_metadata, aes(x = Condition, y = Genes_Detected, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Genes Detected per Cell", x = "Condition", y = "Genes Detected") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_summary(fun = length, geom = "text", aes(label = paste("n =", after_stat(y))), 
                 vjust = -0.5, size = 3)
  
  # Plot 3: Expression of top DE genes
  if(length(sim_data$de_genes) >= 6) {
    top_de_genes <- rownames(sim_data$counts)[sim_data$de_genes[1:6]]
    
    expr_data <- data.frame()
    for(gene in top_de_genes) {
      temp_data <- data.frame(
        Gene = gene,
        Expression = as.numeric(sim_data$counts[gene, ]),
        Group = sim_data$cell_metadata$Group,
        Time = sim_data$cell_metadata$Time,
        Condition = sim_data$cell_metadata$Condition
      )
      expr_data <- rbind(expr_data, temp_data)
    }
    
    # Calculate mean expression per condition
    expr_summary <- expr_data %>%
      group_by(Gene, Group, Time, Condition) %>%
      summarise(Mean_Expression = mean(Expression), .groups = 'drop')
    
    p3 <- ggplot(expr_summary, aes(x = Time, y = Mean_Expression, 
                                   color = Group, group = Group)) +
      geom_point(size = 2) +
      geom_line(linewidth = 1) +
      facet_wrap(~Gene, scales = "free_y", nrow = 2) +
      theme_minimal() +
      labs(title = "Mean Expression of Top DE Genes", 
           x = "Time", y = "Mean Expression")
  } else {
    p3 <- ggplot() + theme_void() + labs(title = "No DE genes to plot")
  }
  
  # Plot 4: PCA
  sce_temp <- sim_data$sce
  sce_temp <- runPCA(sce_temp, ncomponents = 10)
  
  pca_data <- data.frame(
    PC1 = reducedDim(sce_temp, "PCA")[, 1],
    PC2 = reducedDim(sce_temp, "PCA")[, 2],
    Group = sim_data$cell_metadata$Group,
    Time = sim_data$cell_metadata$Time,
    Condition = sim_data$cell_metadata$Condition
  )
  
  p4 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(alpha = 0.6, size = 0.8) +
    theme_minimal() +
    labs(title = "PCA of Simulated Data", x = "PC1", y = "PC2")
  
  return(list(umi_plot = p1, genes_plot = p2, expression_plot = p3, pca_plot = p4))
}

# Function to calculate effect sizes achieved
calculate_achieved_effects <- function(sim_data) {
  
  if(length(sim_data$de_genes) == 0) {
    cat("No DE genes found!\n")
    return(NULL)
  }
  
  # Take first DE gene as example
  example_gene <- rownames(sim_data$counts)[sim_data$de_genes[1]]
  gene_expr <- as.numeric(sim_data$counts[example_gene, ])
  
  # Calculate mean expression per condition
  conditions <- unique(sim_data$cell_metadata$Condition)
  mean_expr <- sapply(conditions, function(cond) {
    mean(gene_expr[sim_data$cell_metadata$Condition == cond])
  })
  
  # Calculate fold changes
  control_t1 <- mean_expr["Control_T1"]
  control_t2 <- mean_expr["Control_T2"] 
  treatment_t1 <- mean_expr["Treatment_T1"]
  treatment_t2 <- mean_expr["Treatment_T2"]
  
  time_effect_control <- control_t2 / control_t1
  time_effect_treatment <- treatment_t2 / treatment_t1
  group_effect_t1 <- treatment_t1 / control_t1
  group_effect_t2 <- treatment_t2 / control_t2
  
  cat("=== Achieved Effect Sizes (Example Gene:", example_gene, ") ===\n")
  cat("Mean expression levels:\n")
  cat("  Control T1:", round(control_t1, 2), "\n")
  cat("  Control T2:", round(control_t2, 2), "\n")
  cat("  Treatment T1:", round(treatment_t1, 2), "\n") 
  cat("  Treatment T2:", round(treatment_t2, 2), "\n")
  cat("\nFold changes:\n")
  cat("  Time effect (Control):", round(time_effect_control, 2), "\n")
  cat("  Time effect (Treatment):", round(time_effect_treatment, 2), "\n")
  cat("  Group effect (T1):", round(group_effect_t1, 2), "\n")
  cat("  Group effect (T2):", round(group_effect_t2, 2), "\n")
  cat("  Interaction (differential time response):", 
      round(time_effect_treatment / time_effect_control, 2), "\n")
}

# Enhanced DE testing function with pseudobulk vs NEBULA comparison
run_de_test_comparison <- function(sim_data, methods = c("pseudobulk", "nebula"), alpha = 0.05) {
  
  results <- list()
  
  for(method in methods) {
    if(method == "pseudobulk") {
      results[[method]] <- run_pseudobulk_deseq2(sim_data, alpha = alpha)
    } else if(method == "nebula") {
      results[[method]] <- run_nebula_mixed_effects(sim_data, alpha = alpha)
    } else {
      stop("Unknown method: ", method)
    }
  }
  
  return(results)
}

# Pseudobulk analysis using DESeq2
run_pseudobulk_deseq2 <- function(sim_data, alpha = 0.05) {
  
  if(!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 package required. Install with: BiocManager::install('DESeq2')")
  }
  
  tryCatch({
    
    # Create pseudobulk by summing counts within each condition
    cell_meta <- sim_data$cell_metadata
    counts_matrix <- sim_data$counts
    
    # Filter to T1 timepoint for simplicity
    t1_cells <- cell_meta$Time == "T1"
    t1_meta <- cell_meta[t1_cells, ]
    t1_counts <- counts_matrix[, t1_cells]
    
    # Create pseudobulk samples
    n_reps <- 6
    
    pseudobulk_counts <- matrix(0, nrow = nrow(t1_counts), ncol = n_reps)
    pseudobulk_meta <- data.frame(
      sample_id = paste0("Rep_", 1:n_reps),
      group = rep(c("Control", "Treatment"), each = n_reps/2),
      stringsAsFactors = FALSE
    )
    rownames(pseudobulk_counts) <- rownames(t1_counts)
    colnames(pseudobulk_counts) <- pseudobulk_meta$sample_id
    
    # Randomly assign cells to replicates
    control_cells <- which(t1_meta$Group == "Control")
    treatment_cells <- which(t1_meta$Group == "Treatment")
    
    control_reps <- split(sample(control_cells), rep(1:(n_reps/2), length.out = length(control_cells)))
    treatment_reps <- split(sample(treatment_cells), rep(1:(n_reps/2), length.out = length(treatment_cells)))
    
    # Sum counts for each replicate
    for(i in 1:(n_reps/2)) {
      pseudobulk_counts[, i] <- rowSums(t1_counts[, control_reps[[i]], drop = FALSE])
      pseudobulk_counts[, i + n_reps/2] <- rowSums(t1_counts[, treatment_reps[[i]], drop = FALSE])
    }
    
    # Remove genes with very low counts
    keep_genes <- rowSums(pseudobulk_counts >= 10) >= 2
    pseudobulk_counts <- pseudobulk_counts[keep_genes, ]
    
    # Create DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = pseudobulk_counts,
      colData = pseudobulk_meta,
      design = ~ group
    )
    
    # Run DESeq2
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c("group", "Treatment", "Control"))
    
    # Extract results
    p_values <- rep(1, nrow(sim_data$counts))
    names(p_values) <- rownames(sim_data$counts)
    
    # Fill in p-values for tested genes
    tested_genes <- rownames(res)[!is.na(res$pvalue)]
    p_values[tested_genes] <- res$pvalue[!is.na(res$pvalue)]
    
    # Multiple testing correction
    p_adj <- p.adjust(p_values, method = "BH")
    significant_genes <- which(p_adj < alpha)
    
  }, error = function(e) {
    warning("Pseudobulk analysis failed: ", e$message)
    p_values <- rep(1, nrow(sim_data$counts))
    names(p_values) <- rownames(sim_data$counts)
    p_adj <- rep(1, nrow(sim_data$counts))
    names(p_adj) <- rownames(sim_data$counts)
    significant_genes <- integer(0)
  })
  
  # Ensure significant_genes is always defined
  if(!exists("significant_genes")) {
    significant_genes <- integer(0)
  }
  
  # Calculate performance metrics
  true_de_genes <- sim_data$de_genes
  true_positives <- sum(significant_genes %in% true_de_genes)
  false_positives <- length(significant_genes) - true_positives
  
  return(list(
    method = "pseudobulk",
    p_values = p_values,
    p_adj = p_adj,
    significant_genes = significant_genes,
    n_significant = length(significant_genes),
    true_positives = true_positives,
    false_positives = false_positives,
    power = true_positives / length(true_de_genes),
    fdr = false_positives / max(length(significant_genes), 1)
  ))
}

# NEBULA mixed effects analysis
run_nebula_mixed_effects <- function(sim_data, alpha = 0.05) {
  
  if(!requireNamespace("nebula", quietly = TRUE)) {
    stop("nebula package required. Install with: install.packages('nebula')")
  }
  
  tryCatch({
    
    # Prepare data for NEBULA
    cell_meta <- sim_data$cell_metadata
    counts_matrix <- sim_data$counts
    
    # Filter to T1 timepoint for simplicity
    t1_cells <- cell_meta$Time == "T1"
    t1_meta <- cell_meta[t1_cells, ]
    t1_counts <- counts_matrix[, t1_cells]
    
    # Create subject IDs
    n_subjects_per_group <- 6
    
    control_indices <- which(t1_meta$Group == "Control")
    treatment_indices <- which(t1_meta$Group == "Treatment")
    
    subject_assignment <- character(nrow(t1_meta))
    subject_assignment[control_indices] <- sample(
      paste0("Control_S", 1:n_subjects_per_group), 
      length(control_indices), replace = TRUE
    )
    subject_assignment[treatment_indices] <- sample(
      paste0("Treatment_S", 1:n_subjects_per_group), 
      length(treatment_indices), replace = TRUE
    )
    
    t1_meta$subject_id <- subject_assignment
    t1_meta$group_numeric <- ifelse(t1_meta$Group == "Treatment", 1, 0)
    
    # Convert to NEBULA format
    gene_means <- rowMeans(t1_counts)
    keep_genes <- gene_means > 0.1 & rowSums(t1_counts > 0) >= 10
    
    if(sum(keep_genes) < 10) {
      keep_genes <- gene_means > 0.01 & rowSums(t1_counts > 0) >= 5
    }
    
    filtered_counts <- t1_counts[keep_genes, ]
    
    # Create NEBULA data object
    data_obj <- nebula::group_cell(
      count = filtered_counts,
      id = t1_meta$subject_id,
      pred = data.frame(
        group = t1_meta$group_numeric,
        intercept = 1
      )
    )
    
    # Run NEBULA-HL
    res <- nebula::nebula(
      data_obj$count,
      data_obj$id, 
      pred = data_obj$pred,
      offset = data_obj$offset,
      model = "NBLMM"
    )
    
    # Extract p-values
    p_values <- rep(1, nrow(sim_data$counts))
    names(p_values) <- rownames(sim_data$counts)
    
    # Fill in p-values for tested genes
    if(nrow(res$summary) > 0) {
      tested_genes <- rownames(res$summary)
      p_values[tested_genes] <- res$summary$p_group
    }
    
    # Multiple testing correction
    p_adj <- p.adjust(p_values, method = "BH")
    significant_genes <- which(p_adj < alpha)
    
  }, error = function(e) {
    warning("NEBULA analysis failed: ", e$message)
    p_values <- rep(1, nrow(sim_data$counts))
    names(p_values) <- rownames(sim_data$counts)
    p_adj <- rep(1, nrow(sim_data$counts))
    names(p_adj) <- rownames(sim_data$counts)
    significant_genes <- integer(0)
  })
  
  # Ensure significant_genes is always defined
  if(!exists("significant_genes")) {
    significant_genes <- integer(0)
  }
  
  # Calculate performance metrics
  true_de_genes <- sim_data$de_genes
  true_positives <- sum(significant_genes %in% true_de_genes)
  false_positives <- length(significant_genes) - true_positives
  
  return(list(
    method = "nebula",
    p_values = p_values,
    p_adj = p_adj,
    significant_genes = significant_genes,
    n_significant = length(significant_genes),
    true_positives = true_positives,
    false_positives = false_positives,
    power = true_positives / length(true_de_genes),
    fdr = false_positives / max(length(significant_genes), 1)
  ))
}

# Method comparison function
run_method_comparison <- function(
    n_simulations = 50,
    methods = c("pseudobulk", "nebula"),
    n_cells_per_condition = 500,
    cell_imbalance_ratio = c(1, 1),
    effect_size = 1.5,
    alpha = 0.05,
    ...
) {
  
  cat("Comparing methods:", paste(methods, collapse = " vs "), "\n")
  cat("Cell imbalance ratio:", paste(cell_imbalance_ratio, collapse = ":"), "\n")
  cat("Running", n_simulations, "simulations...\n")
  
  # Initialize results storage
  results <- data.frame(
    simulation = rep(1:n_simulations, length(methods)),
    method = rep(methods, each = n_simulations),
    power = NA,
    fdr = NA,
    n_significant = NA,
    true_positives = NA,
    false_positives = NA
  )
  
  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
  
  for(i in 1:n_simulations) {
    
    # Generate simulation with imbalance
    sim_data <- simulate_scrna_imbalanced(
      n_cells_per_condition = n_cells_per_condition,
      cell_imbalance_ratio = cell_imbalance_ratio,
      group_effect_size = effect_size,
      seed = i,
      ...
    )
    
    # Run both methods
    method_results <- run_de_test_comparison(sim_data, methods = methods, alpha = alpha)
    
    # Store results
    for(j in seq_along(methods)) {
      row_idx <- (i-1) * length(methods) + j
      method_res <- method_results[[methods[j]]]
      
      results$power[row_idx] <- method_res$power
      results$fdr[row_idx] <- method_res$fdr
      results$n_significant[row_idx] <- method_res$n_significant
      results$true_positives[row_idx] <- method_res$true_positives
      results$false_positives[row_idx] <- method_res$false_positives
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Calculate summary statistics
  summary_stats <- results %>%
    group_by(method) %>%
    summarise(
      mean_power = mean(power, na.rm = TRUE),
      se_power = sd(power, na.rm = TRUE) / sqrt(n()),
      mean_fdr = mean(fdr, na.rm = TRUE),
      se_fdr = sd(fdr, na.rm = TRUE) / sqrt(n()),
      mean_n_significant = mean(n_significant, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\n=== Method Comparison Results ===\n")
  print(summary_stats)
  
  # Create comparison plots
  p1 <- ggplot(results, aes(x = method, y = power, fill = method)) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
    labs(title = "Statistical Power Comparison", 
         x = "Method", y = "Power",
         subtitle = paste("Effect size:", effect_size, "| Cells:", n_cells_per_condition, 
                          "| Imbalance:", paste(cell_imbalance_ratio, collapse = ":"))) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p2 <- ggplot(results, aes(x = method, y = fdr, fill = method)) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(title = "False Discovery Rate Comparison", 
         x = "Method", y = "FDR",
         subtitle = paste("Target FDR =", alpha)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p1)
  print(p2)
  
  return(list(
    individual_results = results,
    summary = summary_stats,
    plots = list(power_plot = p1, fdr_plot = p2),
    parameters = list(
      n_simulations = n_simulations,
      methods = methods,
      effect_size = effect_size,
      n_cells_per_condition = n_cells_per_condition,
      cell_imbalance_ratio = cell_imbalance_ratio
    )
  ))
}

# Multi-effect size comparison
run_multi_effect_comparison <- function(
    effect_sizes = c(1.2, 1.3, 1.8, 2.5),
    effect_labels = c("Very Small", "Small", "Medium", "Large"),
    n_simulations = 50,
    methods = c("pseudobulk", "nebula"),
    n_cells_per_condition = 400,
    alpha = 0.05,
    ...
) {
  
  cat("Running multi-effect size comparison...\n")
  cat("Effect sizes:", paste(effect_sizes, collapse = ", "), "\n")
  cat("Methods:", paste(methods, collapse = " vs "), "\n")
  
  # Initialize results storage
  all_results <- data.frame()
  
  for(i in seq_along(effect_sizes)) {
    
    effect_size <- effect_sizes[i]
    effect_label <- effect_labels[i]
    
    cat("\n--- Testing", effect_label, "effect size (", effect_size, ") ---\n")
    
    # Run comparison for this effect size
    comparison_result <- run_method_comparison(
      n_simulations = n_simulations,
      methods = methods,
      n_cells_per_condition = n_cells_per_condition,
      effect_size = effect_size,
      alpha = alpha,
      ...
    )
    
    # Add effect size information to results
    temp_results <- comparison_result$individual_results
    temp_results$effect_size <- effect_size
    temp_results$effect_label <- effect_label
    
    # Combine with previous results
    all_results <- rbind(all_results, temp_results)
  }
  
  # Calculate summary statistics across all effect sizes
  summary_stats <- all_results %>%
    group_by(effect_label, effect_size, method) %>%
    summarise(
      mean_power = mean(power, na.rm = TRUE),
      se_power = sd(power, na.rm = TRUE) / sqrt(n()),
      mean_fdr = mean(fdr, na.rm = TRUE),
      se_fdr = sd(fdr, na.rm = TRUE) / sqrt(n()),
      mean_n_significant = mean(n_significant, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\n=== Multi-Effect Size Comparison Results ===\n")
  print(summary_stats)
  
  # Create summary table for easy interpretation
  summary_table <- summary_stats %>%
    select(effect_label, method, mean_power, mean_fdr) %>%
    pivot_wider(names_from = method, values_from = c(mean_power, mean_fdr), 
                names_sep = "_") %>%
    mutate(
      power_difference = abs(mean_power_pseudobulk - mean_power_nebula),
      better_power = ifelse(mean_power_pseudobulk > mean_power_nebula, "Pseudobulk", "NEBULA"),
      fdr_difference = abs(mean_fdr_pseudobulk - mean_fdr_nebula),
      better_fdr = ifelse(mean_fdr_pseudobulk < mean_fdr_nebula, "Pseudobulk", "NEBULA")
    )
  
  cat("\n=== Summary Table ===\n")
  print(summary_table)
  
  # Create plots
  p1 <- ggplot(all_results, aes(x = effect_label, y = power, fill = method)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, 
                 position = position_dodge(width = 0.8), color = "white") +
    labs(title = "Statistical Power vs Effect Size", 
         x = "Effect Size", y = "Power") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  print(p1)
  
  return(list(
    all_results = all_results,
    summary_stats = summary_stats,
    summary_table = summary_table,
    plots = list(power_boxplot = p1),
    parameters = list(
      effect_sizes = effect_sizes,
      effect_labels = effect_labels,
      n_simulations = n_simulations,
      methods = methods,
      n_cells_per_condition = n_cells_per_condition
    )
  ))
}

# Imbalanced comparison function
run_imbalanced_comparison <- function(
    effect_sizes = c(1.2, 1.3, 1.8, 2.5),
    effect_labels = c("Very Small", "Small", "Medium", "Large"),
    imbalance_ratios = list(c(1, 1), c(2, 1), c(3, 1), c(5, 1), c(10,1)),
    imbalance_labels = c("Balanced (1:1)", "Low (2:1)", "Moderate (3:1)", "High (5:1)", "Extreme (10:1)"),
    n_simulations = 25,
    methods = c("pseudobulk", "nebula"),
    base_cells = 1000,
    alpha = 0.05,
    ...
) {
  
  cat("Running imbalanced comparison across effect sizes and cell imbalances...\n")
  
  # Initialize results storage
  all_results <- data.frame()
  
  # Progress tracking
  total_combinations <- length(effect_sizes) * length(imbalance_ratios)
  current_combo <- 0
  
  for(i in seq_along(effect_sizes)) {
    for(j in seq_along(imbalance_ratios)) {
      
      current_combo <- current_combo + 1
      effect_size <- effect_sizes[i]
      effect_label <- effect_labels[i]
      imbalance_ratio <- imbalance_ratios[[j]]
      imbalance_label <- imbalance_labels[j]
      
      cat(sprintf("\n--- Combination %d/%d: %s effect + %s imbalance ---\n", 
                  current_combo, total_combinations, effect_label, imbalance_label))
      
      # Run comparison for this combination
      comparison_result <- run_method_comparison(
        n_simulations = n_simulations,
        methods = methods,
        n_cells_per_condition = base_cells,
        effect_size = effect_size,
        alpha = alpha,
        cell_imbalance_ratio = imbalance_ratio,
        ...
      )
      
      # Add experimental condition information
      temp_results <- comparison_result$individual_results
      temp_results$effect_size <- effect_size
      temp_results$effect_label <- effect_label
      temp_results$imbalance_ratio_text <- imbalance_label
      temp_results$imbalance_ratio_1 <- imbalance_ratio[1]
      temp_results$imbalance_ratio_2 <- imbalance_ratio[2]
      temp_results$imbalance_severity <- imbalance_ratio[1] / imbalance_ratio[2]
      
      # Combine with previous results
      all_results <- rbind(all_results, temp_results)
    }
  }
  
  # Calculate comprehensive summary statistics
  summary_stats <- all_results %>%
    group_by(effect_label, effect_size, imbalance_ratio_text, imbalance_severity, method) %>%
    summarise(
      mean_power = mean(power, na.rm = TRUE),
      se_power = sd(power, na.rm = TRUE) / sqrt(n()),
      mean_fdr = mean(fdr, na.rm = TRUE),
      se_fdr = sd(fdr, na.rm = TRUE) / sqrt(n()),
      mean_n_significant = mean(n_significant, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\n=== Comprehensive Imbalanced Comparison Results ===\n")
  print(summary_stats)
  
  # Create power difference heatmap
  power_wide <- summary_stats %>%
    select(effect_label, imbalance_ratio_text, method, mean_power) %>%
    pivot_wider(names_from = method, values_from = mean_power) %>%
    mutate(power_difference = pseudobulk - nebula,
           better_method = ifelse(power_difference > 0, "Pseudobulk", "NEBULA"))
  
  cat("\n=== Power Difference Summary ===\n")
  print(power_wide)
  
  # Create heatmap plot
  p1 <- ggplot(power_wide, aes(x = effect_label, y = imbalance_ratio_text, fill = power_difference)) +
    geom_tile() +
    geom_text(aes(label = round(power_difference, 3)), color = "white", fontface = "bold") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Power\nDifference\n(Pseudo - NEBULA)") +
    labs(title = "Power Difference Heatmap", 
         x = "Effect Size", y = "Cell Imbalance",
         subtitle = "Positive values = Pseudobulk better, Negative = NEBULA better") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p1)
  
  return(list(
    all_results = all_results,
    summary_stats = summary_stats,
    power_differences = power_wide,
    plots = list(power_difference_heatmap = p1),
    parameters = list(
      effect_sizes = effect_sizes,
      imbalance_ratios = imbalance_ratios,
      n_simulations = n_simulations,
      methods = methods,
      base_cells = base_cells
    )
  ))
}

# Enhanced comparison examples
enhanced_method_comparison_example <- function() {
  
  cat("=== Enhanced Pseudobulk vs NEBULA Comparison ===\n")
  cat("Testing 4 effect sizes: Very Small (1.2x), Small (1.3x), Medium (1.8x), Large (2.5x)\n")
  
  # Run multi-effect size comparison
  multi_results <- run_multi_effect_comparison(
    effect_sizes = c(1.2, 1.3, 1.8, 2.5),
    effect_labels = c("Very Small (1.2x)", "Small (1.3x)", "Medium (1.8x)", "Large (2.5x)"),
    n_simulations = 30,
    methods = c("pseudobulk", "nebula"),
    n_cells_per_condition = 350,
    n_genes = 600,
    de_prob = 0.12,
    alpha = 0.05
  )
  
  return(multi_results)
}

# imbalanced_method_comparison_example <- function() {
#   
#   cat("=== Imbalanced Cell Type Analysis ===\n")
#   cat("Testing how cell imbalance affects DE method performance\n")
#   
#  # Run comprehensive imbalanced comparison
#   imbalanced_results <- run_imbalanced_comparison(
#     effect_sizes = c(1.2, 1.8, 2.5),
#     effect_labels = c("Small (1.2x)", "Medium (1.8x)", "Large (2.5x)"),
#     imbalance_ratios = list(c(1, 1), c(2, 1), c(5, 1)),
#     imbalance_labels = c("Balanced (1:1)", "Moderate (2:1)", "Extreme (5:1)"),
#     n_simulations = 50,
#     methods = c("pseudobulk", "nebula"),
#     base_cells = 5000,
#     n_genes = 500,
#     de_prob = 0.1,
#     alpha = 0.05
#   )
#   
#   return(imbalanced_results)
# }

imbalanced_method_comparison_example <- function() {
  cat("=== Imbalanced Cell Type Analysis ===\n")
  cat("Testing how cell imbalance affects DE method performance\n")
  
  # Use ALL combinations (this will take longer)
  imbalanced_results <- run_imbalanced_comparison(
    # Remove the parameters to use defaults, or specify all:
    effect_sizes = c(1.2, 1.3, 1.8, 2.5),
    effect_labels = c("Very Small (1.2x)", "Small (1.3x)", "Medium (1.8x)", "Large (2.5x)"),
    imbalance_ratios = list(c(1, 1), c(2, 1), c(3, 1), c(5, 1), c(10,1)),
    imbalance_labels = c("Balanced (1:1)", "Low (2:1)", "Moderate (3:1)", "High (5:1)", "Extreme (10:1)"),
    n_simulations = 20,
    methods = c("pseudobulk", "nebula"),
    base_cells = 5000,
    n_genes = 500,
    de_prob = 0.1,
    alpha = 0.05
  )
  
  return(imbalanced_results)
}

# Example usage and testing
cat("\n=== SINGLE CELL RNA-SEQ SIMULATION READY ===\n")
cat("Running a test simulation...\n")

# Test simulation
test_sim <- simulate_scrna_imbalanced(
  n_cells_per_condition = 100,
  cell_imbalance_ratio = c(2, 1),
  n_genes = 300,
  group_effect_size = 1.5,
  seed = 123
)

cat("Test simulation completed!\n")
cat("Cell counts per condition:\n")
print(table(test_sim$cell_metadata$Condition))

# Generate plots
test_plots <- plot_splatter_simulation(test_sim)
print(test_plots$umi_plot)

# Show achieved effects
calculate_achieved_effects(test_sim)

cat("\n=== SAMPLE USAGE EXAMPLES ===\n")
cat("1. Basic simulation:\n")
cat("   sim <- simulate_scrna_splatter(n_cells_per_condition=200, n_genes=500)\n")
cat("2. Imbalanced simulation:\n")
cat("   sim <- simulate_scrna_imbalanced(cell_imbalance_ratio=c(3,1))\n")
cat("3. Method comparison:\n")
cat("   results <- run_method_comparison(n_simulations=20)\n")
cat("4. Multi-effect size analysis (WITH SUMMARY TABLE):\n")
cat("   multi_results <- enhanced_method_comparison_example()\n")
cat("   print(multi_results$summary_table)  # <-- THE SUMMARY TABLE!\n")
cat("5. Full imbalance analysis (WITH HEATMAP):\n")
cat("   imbalance_results <- imbalanced_method_comparison_example()\n")
cat("   print(imbalance_results$power_differences)  # <-- POWER COMPARISON TABLE!\n")

cat("\n=== RUNNING COMPREHENSIVE ANALYSIS ===\n")
cat("Uncomment the sections below to run the full analyses:\n")

cat("\n# ===== MULTI-EFFECT SIZE ANALYSIS ===== #\n")
cat(" Uncomment this block to run:\n")
cat(" multi_results <- enhanced_method_comparison_example()\n")
cat(" print(multi_results$summary_table)\n")

# To actually run it, uncomment these lines:
# cat("Running multi-effect size analysis...\n")
# multi_results <- enhanced_method_comparison_example()
# cat("\nSUMMARY TABLE:\n")
# print(multi_results$summary_table)

cat("\n# ===== IMBALANCE ANALYSIS ===== #\n") 
cat(" Uncomment this block to run:\n")
cat(" imbalance_results <- imbalanced_method_comparison_example()\n")
cat(" print(imbalance_results$power_differences)\n")

# To actually run it, uncomment these lines:
# cat("Running imbalance analysis...\n")
# imbalance_results <- imbalanced_method_comparison_example()
# cat("\nPOWER DIFFERENCES TABLE:\n")
# print(imbalance_results$power_differences)
# cat("\nDETAILED SUMMARY STATS:\n")
# print(imbalance_results$summary_stats)

cat("\n=== TO RUN ANALYSES: UNCOMMENT THE BLOCKS ABOVE ===\n")
cat("The # symbols at the start of lines prevent execution.\n")
cat("Remove the # to actually run the comprehensive analyses.\n")

cat("\n=== QUICK DEMO AVAILABLE ===\n")
cat("For a quick demo with fewer simulations, uncomment this:\n")
cat("demo_results <- run_method_comparison(n_simulations=5, n_cells_per_condition=100)\n")

# Quick demo (uncomment to run)
# cat("Running quick demo...\n")
# demo_results <- run_method_comparison(n_simulations=5, n_cells_per_condition=100)
# cat("Demo completed!\n")

cat("\n=== READY TO USE! ===\n")

imbalance_results <- imbalanced_method_comparison_example()
print(imbalance_results$power_differences)