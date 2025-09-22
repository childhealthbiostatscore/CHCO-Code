# Single Cell RNA-seq Data Simulation Using Splatter
# Two groups, two timepoints, with prespecified effect sizes

library(splatter)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(ggplot2)
library(patchwork)

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

cat("\n=== Custom DE Method Support Added ===\n")
cat("Use generate_simulation_batch() to save simulations\n")
cat("Use apply_custom_de_method() to test your DE method\n")
cat("Example: custom_de_workflow_example()\n")

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