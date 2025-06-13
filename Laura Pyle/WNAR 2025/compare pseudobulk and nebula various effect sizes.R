# Single Cell RNA-seq Data Simulation Using Splatter
# Two groups, two timepoints, with prespecified effect sizes

library(splatter)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

# Function to simulate scRNA-seq data using Splatter
simulate_scrna_splatter <- function(
    n_cells_per_condition = 500,   # cells per group-timepoint combination
    n_genes = 2000,                # total number of genes
    de_prob = 0.1,                 # probability a gene is DE (10% = 200 DE genes if n_genes=2000)
    de_facLoc = 1.0,              # location parameter for DE factor (controls effect size)
    de_facScale = 0.4,            # scale parameter for DE factor (controls variability)
    group_effect_size = 1.5,      # multiplier for group effects
    time_effect_size = 1.2,       # multiplier for time effects
    interaction_effect_size = 1.3, # multiplier for interaction effects
    dropout_type = "experiment",   # dropout type
    seed = 123
) {
  
  set.seed(seed)
  
  # Set up parameters for 4 conditions (2 groups × 2 timepoints)
  n_conditions <- 4
  total_cells <- n_cells_per_condition * n_conditions
  
  # Create group probabilities (equal cells per condition)
  group_prob <- rep(1/n_conditions, n_conditions)
  
  # Estimate parameters from a reference (using default Splat parameters)
  params <- newSplatParams()
  
  # Set parameters using the updated API (nGroups is determined by length of group.prob)
  params <- setParams(params, 
                      nGenes = n_genes,
                      batchCells = total_cells,
                      group.prob = group_prob,  # This determines nGroups automatically
                      de.prob = de_prob,
                      de.facLoc = de_facLoc,
                      de.facScale = de_facScale,
                      dropout.type = dropout_type,
                      dropout.mid = 3,
                      dropout.shape = -1,
                      seed = seed)
  
  # Simulate the base data
  sim <- splatSimulate(params, method = "groups", verbose = FALSE)
  
  # Extract count matrix and convert to regular matrix for easier manipulation
  counts <- as.matrix(counts(sim))
  
  # Create condition labels
  n_cells_total <- ncol(counts)
  cells_per_condition <- n_cells_total / 4
  
  condition_labels <- rep(c("Control_T1", "Control_T2", "Treatment_T1", "Treatment_T2"), 
                          each = cells_per_condition)
  
  # Create metadata
  cell_metadata <- data.frame(
    CellID = colnames(counts),
    Group = rep(c("Control", "Control", "Treatment", "Treatment"), each = cells_per_condition),
    Time = rep(c("T1", "T2", "T1", "T2"), each = cells_per_condition),
    Condition = condition_labels,
    stringsAsFactors = FALSE
  )
  
  # Add custom effects to DE genes based on our experimental design
  gene_metadata <- as.data.frame(rowData(sim))
  de_genes <- which(gene_metadata$DEFacGroup1 != 1 | 
                      gene_metadata$DEFacGroup2 != 1 |
                      gene_metadata$DEFacGroup3 != 1 |
                      gene_metadata$DEFacGroup4 != 1)
  
  cat("Number of DE genes identified:", length(de_genes), "\n")
  
  # Apply custom effect sizes to DE genes
  modified_counts <- counts
  
  for(gene_idx in de_genes) {
    base_expression <- counts[gene_idx, ]
    
    # Get indices for each condition
    control_t1_idx <- which(cell_metadata$Condition == "Control_T1")
    control_t2_idx <- which(cell_metadata$Condition == "Control_T2") 
    treatment_t1_idx <- which(cell_metadata$Condition == "Treatment_T1")
    treatment_t2_idx <- which(cell_metadata$Condition == "Treatment_T2")
    
    # Apply effects (multiplicative on count scale)
    # Time effect: T2 vs T1 (applies to both groups)
    modified_counts[gene_idx, control_t2_idx] <- 
      round(modified_counts[gene_idx, control_t2_idx] * time_effect_size)
    modified_counts[gene_idx, treatment_t2_idx] <- 
      round(modified_counts[gene_idx, treatment_t2_idx] * time_effect_size)
    
    # Group effect: Treatment vs Control (applies to both timepoints)  
    modified_counts[gene_idx, treatment_t1_idx] <- 
      round(modified_counts[gene_idx, treatment_t1_idx] * group_effect_size)
    modified_counts[gene_idx, treatment_t2_idx] <- 
      round(modified_counts[gene_idx, treatment_t2_idx] * group_effect_size)
    
    # Interaction effect: additional effect for Treatment at T2
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
      n_cells_per_condition = n_cells_per_condition,
      n_genes = n_genes,
      de_prob = de_prob,
      group_effect_size = group_effect_size,
      time_effect_size = time_effect_size,
      interaction_effect_size = interaction_effect_size
    )
  ))
}

# Function to plot simulation results
plot_splatter_simulation <- function(sim_data) {
  
  # Plot 1: UMI distribution
  p1 <- ggplot(sim_data$cell_metadata, aes(x = Condition, y = Total_UMI, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Total UMI per Cell", x = "Condition", y = "Total UMI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 2: Genes detected
  p2 <- ggplot(sim_data$cell_metadata, aes(x = Condition, y = Genes_Detected, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Genes Detected per Cell", x = "Condition", y = "Genes Detected") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
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

# Function for power analysis
run_power_analysis <- function(
    n_simulations = 100,               # Number of simulated datasets
    n_cells_per_condition = 250,
    effect_size = 1.5,                 # Effect size to test
    alpha = 0.05,                      # Significance threshold
    test_type = "wilcox",              # Test type: "wilcox", "t.test", etc.
    ...                                # Other parameters passed to simulation
) {
  
  cat("Running power analysis with", n_simulations, "simulations...\n")
  
  # Store results
  results <- data.frame(
    simulation = 1:n_simulations,
    n_de_detected = NA,
    n_true_de = NA,
    power = NA,
    fdr = NA
  )
  
  # Progress tracking
  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
  
  for(i in 1:n_simulations) {
    
    # Generate simulation with different seed each time
    sim_data <- simulate_scrna_splatter(
      n_cells_per_condition = n_cells_per_condition,
      group_effect_size = effect_size,
      seed = i,  # Different seed for each simulation
      ...
    )
    
    # Run differential expression test (example with simple t-test)
    de_results <- run_de_test(sim_data, test_type = test_type, alpha = alpha)
    
    # Calculate metrics
    results$n_de_detected[i] <- de_results$n_significant
    results$n_true_de[i] <- length(sim_data$de_genes)
    results$power[i] <- de_results$true_positives / length(sim_data$de_genes)
    results$fdr[i] <- de_results$false_positives / max(de_results$n_significant, 1)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Summary statistics
  summary_stats <- data.frame(
    mean_power = mean(results$power, na.rm = TRUE),
    se_power = sd(results$power, na.rm = TRUE) / sqrt(n_simulations),
    mean_fdr = mean(results$fdr, na.rm = TRUE),
    se_fdr = sd(results$fdr, na.rm = TRUE) / sqrt(n_simulations)
  )
  
  cat("\n=== Power Analysis Results ===\n")
  cat("Effect size tested:", effect_size, "\n")
  cat("Mean power:", round(summary_stats$mean_power, 3), "±", round(summary_stats$se_power, 3), "\n")
  cat("Mean FDR:", round(summary_stats$mean_fdr, 3), "±", round(summary_stats$se_fdr, 3), "\n")
  
  return(list(
    individual_results = results,
    summary = summary_stats,
    parameters = list(
      n_simulations = n_simulations,
      effect_size = effect_size,
      n_cells_per_condition = n_cells_per_condition
    )
  ))
}

# Simple DE testing function (you'd replace this with your preferred method)
run_de_test <- function(sim_data, test_type = "wilcox", alpha = 0.05) {
  
  # Example: test group effect at T1 timepoint
  t1_cells <- sim_data$cell_metadata$Time == "T1"
  control_cells <- sim_data$cell_metadata$Group == "Control" & t1_cells
  treatment_cells <- sim_data$cell_metadata$Group == "Treatment" & t1_cells
  
  n_genes <- nrow(sim_data$counts)
  p_values <- numeric(n_genes)
  
  # Run statistical test for each gene
  for(g in 1:n_genes) {
    control_expr <- sim_data$counts[g, control_cells]
    treatment_expr <- sim_data$counts[g, treatment_cells]
    
    if(test_type == "wilcox") {
      test_result <- wilcox.test(treatment_expr, control_expr)
    } else if(test_type == "t.test") {
      test_result <- t.test(treatment_expr, control_expr)
    }
    
    p_values[g] <- test_result$p.value
  }
  
  # Multiple testing correction
  p_adj <- p.adjust(p_values, method = "BH")
  significant_genes <- which(p_adj < alpha)
  
  # Calculate true/false positives
  true_de_genes <- sim_data$de_genes
  true_positives <- sum(significant_genes %in% true_de_genes)
  false_positives <- length(significant_genes) - true_positives
  
  return(list(
    p_values = p_values,
    p_adj = p_adj,
    significant_genes = significant_genes,
    n_significant = length(significant_genes),
    true_positives = true_positives,
    false_positives = false_positives
  ))
}

# Example: Power analysis across different effect sizes
power_analysis_example <- function() {
  effect_sizes <- c(1.2, 1.5, 2.0, 2.5)
  power_results <- list()
  
  for(i in seq_along(effect_sizes)) {
    cat("\n--- Testing effect size:", effect_sizes[i], "---\n")
    power_results[[i]] <- run_power_analysis(
      n_simulations = 500,  # Increase this for real analysis (500-1000)
      effect_size = effect_sizes[i],
      n_cells_per_condition = 200,
      n_genes = 500,
      de_prob = 0.1
    )
  }
  
  # Combine results
  power_summary <- data.frame(
    effect_size = effect_sizes,
    power = sapply(power_results, function(x) x$summary$mean_power),
    power_se = sapply(power_results, function(x) x$summary$se_power)
  )
  
  # Plot power curve
  library(ggplot2)
  p <- ggplot(power_summary, aes(x = effect_size, y = power)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = power - power_se, ymax = power + power_se), width = 0.05) +
    ylim(0, 1) +
    labs(title = "Power Analysis", x = "Effect Size (Fold Change)", y = "Statistical Power") +
    theme_minimal()
  
  print(p)
  
  return(power_summary)
}

# Example usage - Single simulation
cat("=== Single Simulation Example ===\n")
sim_data <- simulate_scrna_splatter(
  n_cells_per_condition = 250,
  n_genes = 1000,
  de_prob = 0.15,
  group_effect_size = 1.8,
  time_effect_size = 1.5,
  interaction_effect_size = 1.3,
  seed = 123
)

# Generate and display plots
cat("\n=== Generating Plots ===\n")
plots <- plot_splatter_simulation(sim_data)

# Display plots
print(plots$umi_plot + plots$genes_plot)
print(plots$expression_plot)
print(plots$pca_plot)

# Calculate and display effect sizes
calculate_achieved_effects(sim_data)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total cells:", ncol(sim_data$counts), "\n")
cat("Total genes:", nrow(sim_data$counts), "\n") 
cat("DE genes:", length(sim_data$de_genes), "\n")
cat("Median UMI per cell:", median(sim_data$cell_metadata$Total_UMI), "\n")
cat("Median genes per cell:", median(sim_data$cell_metadata$Genes_Detected), "\n")

cat("\n=== Usage ===\n")
cat("Access the SingleCellExperiment object: sim_data$sce\n")
cat("Access raw counts: sim_data$counts\n")
cat("Access metadata: sim_data$cell_metadata\n")
cat("DE gene indices: sim_data$de_genes\n")

# Function to generate and save multiple simulations for custom DE methods
generate_simulation_batch <- function(
    n_simulations = 100,
    save_directory = "simulations",
    file_prefix = "sim_",
    ...  # Parameters passed to simulate_scrna_splatter
) {
  
  # Create directory if it doesn't exist
  if(!dir.exists(save_directory)) {
    dir.create(save_directory)
  }
  
  cat("Generating", n_simulations, "simulations...\n")
  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
  
  simulation_info <- data.frame(
    simulation_id = 1:n_simulations,
    file_path = NA,
    n_cells = NA,
    n_genes = NA,
    n_de_genes = NA,
    seed = NA
  )
  
  for(i in 1:n_simulations) {
    
    # Generate simulation
    sim_data <- simulate_scrna_splatter(seed = i, ...)
    
    # Save as RDS file (preserves all R objects)
    file_path <- file.path(save_directory, paste0(file_prefix, sprintf("%03d", i), ".rds"))
    saveRDS(sim_data, file_path)
    
    # Store metadata
    simulation_info$file_path[i] <- file_path
    simulation_info$n_cells[i] <- ncol(sim_data$counts)
    simulation_info$n_genes[i] <- nrow(sim_data$counts)
    simulation_info$n_de_genes[i] <- length(sim_data$de_genes)
    simulation_info$seed[i] <- i
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Save simulation metadata
  write.csv(simulation_info, file.path(save_directory, "simulation_info.csv"), row.names = FALSE)
  
  cat("\nSimulations saved to:", save_directory, "\n")
  cat("Load individual simulations with: readRDS('path/to/simulation.rds')\n")
  
  return(simulation_info)
}

# Function to apply custom DE method to saved simulations
apply_custom_de_method <- function(
    simulation_info,  # Output from generate_simulation_batch
    custom_de_function,  # Your DE method function
    alpha = 0.05,
    ...  # Additional parameters for your DE method
) {
  
  n_sims <- nrow(simulation_info)
  results <- data.frame(
    simulation_id = simulation_info$simulation_id,
    n_de_detected = NA,
    n_true_de = NA,
    power = NA,
    fdr = NA
  )
  
  cat("Applying custom DE method to", n_sims, "simulations...\n")
  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
  
  for(i in 1:n_sims) {
    
    # Load simulation
    sim_data <- readRDS(simulation_info$file_path[i])
    
    # Apply your custom DE method
    de_results <- custom_de_function(
      counts = sim_data$counts,
      cell_metadata = sim_data$cell_metadata,
      gene_metadata = sim_data$gene_metadata,
      alpha = alpha,
      ...
    )
    
    # Extract results (assumes your function returns p-values or similar)
    significant_genes <- de_results$significant_genes  # Adjust based on your output
    true_de_genes <- sim_data$de_genes
    
    # Calculate metrics
    true_positives <- sum(significant_genes %in% true_de_genes)
    false_positives <- length(significant_genes) - true_positives
    
    results$n_de_detected[i] <- length(significant_genes)
    results$n_true_de[i] <- length(true_de_genes)
    results$power[i] <- true_positives / length(true_de_genes)
    results$fdr[i] <- false_positives / max(length(significant_genes), 1)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Calculate summary statistics
  summary_stats <- data.frame(
    mean_power = mean(results$power, na.rm = TRUE),
    se_power = sd(results$power, na.rm = TRUE) / sqrt(n_sims),
    mean_fdr = mean(results$fdr, na.rm = TRUE),
    se_fdr = sd(results$fdr, na.rm = TRUE) / sqrt(n_sims)
  )
  
  cat("\n=== Custom DE Method Results ===\n")
  cat("Mean power:", round(summary_stats$mean_power, 3), "±", round(summary_stats$se_power, 3), "\n")
  cat("Mean FDR:", round(summary_stats$mean_fdr, 3), "±", round(summary_stats$se_fdr, 3), "\n")
  
  return(list(
    individual_results = results,
    summary = summary_stats
  ))
}

# Example workflow for custom DE methods
custom_de_workflow_example <- function() {
  
  cat("=== Example: Custom DE Method Workflow ===\n")
  
  # Step 1: Generate and save simulations
  sim_info <- generate_simulation_batch(
    n_simulations = 20,  # Use more for real analysis
    n_cells_per_condition = 200,
    n_genes = 500,
    group_effect_size = 1.8,
    save_directory = "my_simulations"
  )
  
  # Step 2: Define your custom DE function (example template)
  my_custom_de_method <- function(counts, cell_metadata, gene_metadata, alpha = 0.05) {
    
    # Your novel DE method goes here
    # This is just a placeholder example
    
    # Example: Simple fold change + t-test (replace with your method)
    group1_cells <- cell_metadata$Group == "Control" & cell_metadata$Time == "T1"
    group2_cells <- cell_metadata$Group == "Treatment" & cell_metadata$Time == "T1"
    
    p_values <- numeric(nrow(counts))
    for(g in 1:nrow(counts)) {
      test_result <- t.test(counts[g, group2_cells], counts[g, group1_cells])
      p_values[g] <- test_result$p.value
    }
    
    # Multiple testing correction
    p_adj <- p.adjust(p_values, method = "BH")
    significant_genes <- which(p_adj < alpha)
    
    return(list(
      p_values = p_values,
      p_adj = p_adj,
      significant_genes = significant_genes
    ))
  }
  
  # Step 3: Apply your method to all simulations
  power_results <- apply_custom_de_method(
    simulation_info = sim_info,
    custom_de_function = my_custom_de_method,
    alpha = 0.05
  )
  
  return(power_results)
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
  
  library(DESeq2)
  library(dplyr)
  
  tryCatch({
    
    # Create pseudobulk by summing counts within each condition
    # For this example, we'll compare Treatment vs Control at T1
    cell_meta <- sim_data$cell_metadata
    counts_matrix <- sim_data$counts
    
    # Filter to T1 timepoint for simplicity
    t1_cells <- cell_meta$Time == "T1"
    t1_meta <- cell_meta[t1_cells, ]
    t1_counts <- counts_matrix[, t1_cells]
    
    # Create pseudobulk samples (one per group)
    # In real analysis, you'd have multiple biological replicates
    # Here we simulate by randomly splitting cells into "replicates"
    n_reps <- 6  # 3 reps per group
    
    pseudobulk_counts <- matrix(0, nrow = nrow(t1_counts), ncol = n_reps)
    pseudobulk_meta <- data.frame(
      sample_id = paste0("Rep_", 1:n_reps),
      group = rep(c("Control", "Treatment"), each = n_reps/2),
      stringsAsFactors = FALSE
    )
    rownames(pseudobulk_counts) <- rownames(t1_counts)
    colnames(pseudobulk_counts) <- pseudobulk_meta$sample_id
    
    # Randomly assign cells to replicates within each group
    control_cells <- which(t1_meta$Group == "Control")
    treatment_cells <- which(t1_meta$Group == "Treatment")
    
    # Split cells into replicates
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
    dds <- DESeqDataSetFromMatrix(
      countData = pseudobulk_counts,
      colData = pseudobulk_meta,
      design = ~ group
    )
    
    # Run DESeq2
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("group", "Treatment", "Control"))
    
    # Extract results
    p_values <- rep(1, nrow(sim_data$counts))  # Initialize with 1s
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
    p_adj <- rep(1, nrow(sim_data$counts))
    significant_genes <- integer(0)
  })
  
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
  
  library(nebula)
  
  tryCatch({
    
    # Prepare data for NEBULA
    cell_meta <- sim_data$cell_metadata
    counts_matrix <- sim_data$counts
    
    # Filter to T1 timepoint for simplicity
    t1_cells <- cell_meta$Time == "T1"
    t1_meta <- cell_meta[t1_cells, ]
    t1_counts <- counts_matrix[, t1_cells]
    
    # Create subject IDs (simulating multiple subjects per group)
    # In real data, these would be actual subject/patient IDs
    n_subjects_per_group <- 6
    subject_ids <- c(
      rep(paste0("Control_S", 1:n_subjects_per_group), 
          length.out = sum(t1_meta$Group == "Control")),
      rep(paste0("Treatment_S", 1:n_subjects_per_group), 
          length.out = sum(t1_meta$Group == "Treatment"))
    )
    
    # Randomly assign subjects within groups
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
    # Remove genes with very low expression
    gene_means <- rowMeans(t1_counts)
    keep_genes <- gene_means > 0.1 & rowSums(t1_counts > 0) >= 10
    
    if(sum(keep_genes) < 10) {
      # If too few genes pass filter, relax criteria
      keep_genes <- gene_means > 0.01 & rowSums(t1_counts > 0) >= 5
    }
    
    filtered_counts <- t1_counts[keep_genes, ]
    
    # Create NEBULA data object
    data_obj <- group_cell(
      count = filtered_counts,
      id = t1_meta$subject_id,
      pred = data.frame(
        group = t1_meta$group_numeric,
        intercept = 1
      )
    )
    
    # Run NEBULA-HL (negative binomial mixed model)
    res <- nebula(
      data_obj$count,
      data_obj$id, 
      pred = data_obj$pred,
      offset = data_obj$offset,
      model = "NBLMM"  # Negative binomial linear mixed model
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
    p_adj <- rep(1, nrow(sim_data$counts))
    significant_genes <- integer(0)
  })
  
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

# Power analysis comparing methods
run_method_comparison <- function(
    n_simulations = 100,
    methods = c("pseudobulk", "nebula"),
    n_cells_per_condition = 500,
    effect_size = 1.5,
    alpha = 0.05,
    ...
) {
  
  cat("Comparing methods:", paste(methods, collapse = " vs "), "\n")
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
    
    # Generate simulation
    sim_data <- simulate_scrna_splatter(
      n_cells_per_condition = n_cells_per_condition,
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
  library(ggplot2)
  
  # Power comparison
  p1 <- ggplot(results, aes(x = method, y = power, fill = method)) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
    labs(title = "Statistical Power Comparison", 
         x = "Method", y = "Power",
         subtitle = paste("Effect size:", effect_size, "| Cells per condition:", n_cells_per_condition)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # FDR comparison
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
      n_cells_per_condition = n_cells_per_condition
    )
  ))
}

# Example method comparison
method_comparison_example <- function() {
  
  cat("=== Pseudobulk vs NEBULA Method Comparison ===\n")
  cat("Note: Install required packages with:\n")
  cat("# install.packages('DESeq2') # from Bioconductor\n") 
  cat("# install.packages('nebula')\n")
  
  # Run comparison with moderate effect size
  comparison_results <- run_method_comparison(
    n_simulations = 50,  # Increase for real analysis
    methods = c("pseudobulk", "nebula"),
    n_cells_per_condition = 400,
    n_genes = 800,
    de_prob = 0.1,
    effect_size = 1.8,  # 1.8-fold change
    alpha = 0.05
  )
  
  return(comparison_results)
}

# Multi-effect size comparison
run_multi_effect_comparison <- function(
    effect_sizes = c(1.3, 1.8, 2.5),  # small, medium, large
    effect_labels = c("Small", "Medium", "Large"),
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
  
  # Create comprehensive plots
  library(ggplot2)
  
  # Power across effect sizes
  p1 <- ggplot(all_results, aes(x = effect_label, y = power, fill = method)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, 
                 position = position_dodge(width = 0.8), color = "white") +
    labs(title = "Statistical Power vs Effect Size", 
         x = "Effect Size", y = "Power",
         subtitle = paste("Cells per condition:", n_cells_per_condition, "| Simulations per effect:", n_simulations)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  # FDR across effect sizes
  p2 <- ggplot(all_results, aes(x = effect_label, y = fdr, fill = method)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, 
                 position = position_dodge(width = 0.8), color = "white") +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(title = "False Discovery Rate vs Effect Size", 
         x = "Effect Size", y = "FDR",
         subtitle = paste("Target FDR =", alpha)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  # Power curves (line plot)
  power_summary <- summary_stats %>%
    select(effect_size, effect_label, method, mean_power, se_power)
  
  p3 <- ggplot(power_summary, aes(x = effect_size, y = mean_power, color = method)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_power - se_power, ymax = mean_power + se_power), 
                  width = 0.05) +
    labs(title = "Power Curves", 
         x = "Effect Size (Fold Change)", y = "Mean Power",
         subtitle = "Error bars show standard error") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_brewer(type = "qual", palette = "Set1") +
    ylim(0, 1)
  
  # Number of significant genes detected
  p4 <- ggplot(all_results, aes(x = effect_label, y = n_significant, fill = method)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, 
                 position = position_dodge(width = 0.8), color = "white") +
    labs(title = "Number of Significant Genes Detected", 
         x = "Effect Size", y = "Number of Significant Genes") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  # Print plots
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
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
  
  return(list(
    all_results = all_results,
    summary_stats = summary_stats,
    summary_table = summary_table,
    plots = list(
      power_boxplot = p1,
      fdr_boxplot = p2, 
      power_curves = p3,
      n_significant = p4
    ),
    parameters = list(
      effect_sizes = effect_sizes,
      effect_labels = effect_labels,
      n_simulations = n_simulations,
      methods = methods,
      n_cells_per_condition = n_cells_per_condition
    )
  ))
}

# Enhanced method comparison example with multiple effect sizes
enhanced_method_comparison_example <- function() {
  
  cat("=== Enhanced Pseudobulk vs NEBULA Comparison ===\n")
  cat("Testing 3 effect sizes: Small (1.3x), Medium (1.8x), Large (2.5x)\n")
  cat("Note: Install required packages with:\n")
  cat("# BiocManager::install('DESeq2')\n") 
  cat("# install.packages('nebula')\n")
  
  # Run multi-effect size comparison
  multi_results <- run_multi_effect_comparison(
    effect_sizes = c(1.3, 1.8, 2.5),  # small, medium, large fold changes
    effect_labels = c("Small (1.3x)", "Medium (1.8x)", "Large (2.5x)"),
    n_simulations = 30,  # Reduced for faster testing, increase to 100+ for publication
    methods = c("pseudobulk", "nebula"),
    n_cells_per_condition = 350,
    n_genes = 600,
    de_prob = 0.12,  # 12% of genes are DE
    alpha = 0.05
  )
  
  # Print key insights
  cat("\n=== Key Insights ===\n")
  
  # Find which method performs better for each effect size
  insights <- multi_results$summary_table
  for(i in 1:nrow(insights)) {
    row <- insights[i, ]
    cat("Effect size", row$effect_label, ":\n")
    cat("  Power: ", row$better_power, " is better (difference: ", 
        round(row$power_difference, 3), ")\n")
    cat("  FDR control: ", row$better_fdr, " is better (difference: ", 
        round(row$fdr_difference, 3), ")\n")
  }
  
  return(multi_results)
}

cat("\n=== Multi-Effect Size Comparison Added ===\n")
cat("Tests small (1.3x), medium (1.8x), and large (2.5x) effect sizes\n")
cat("Usage: enhanced_method_comparison_example()\n")

cat("\n=== Running Enhanced Method Comparison ===\n")
cat("Uncomment the line below to run the multi-effect size comparison:\n")
cat("# multi_results <- enhanced_method_comparison_example()\n")

# To run automatically, uncomment this line:
# multi_results <- enhanced_method_comparison_example()

# Generate plots
plots <- plot_splatter_simulation(sim_data)

# Display plots
print(plots$umi_plot + plots$genes_plot)
print(plots$expression_plot)
print(plots$pca_plot)

# Calculate and display effect sizes
calculate_achieved_effects(sim_data)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total cells:", ncol(sim_data$counts), "\n")
cat("Total genes:", nrow(sim_data$counts), "\n") 
cat("DE genes:", length(sim_data$de_genes), "\n")
cat("Median UMI per cell:", median(sim_data$cell_metadata$Total_UMI), "\n")
cat("Median genes per cell:", median(sim_data$cell_metadata$Genes_Detected), "\n")

cat("\n=== Usage ===\n")
cat("Access the SingleCellExperiment object: sim_data$sce\n")
cat("Access raw counts: sim_data$counts\n")
cat("Access metadata: sim_data$cell_metadata\n")
cat("DE gene indices: sim_data$de_genes\n")

# Run the enhanced comparison
multi_results <- enhanced_method_comparison_example()

# View the summary table
print(multi_results$summary_table)

# Access specific results
print(multi_results$summary_stats)