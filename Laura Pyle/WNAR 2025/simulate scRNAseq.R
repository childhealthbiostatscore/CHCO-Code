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
  
  # Set basic parameters
  params <- setParams(params, 
                      nGenes = n_genes,
                      nCells = total_cells,
                      seed = seed)
  
  # Set group parameters
  params <- setParams(params,
                      group.prob = group_prob,
                      de.prob = de_prob,
                      de.facLoc = de_facLoc,
                      de.facScale = de_facScale)
  
  # Set dropout parameters
  params <- setParams(params,
                      dropout.type = dropout_type,
                      dropout.mid = 3,    # midpoint of dropout logistic
                      dropout.shape = -1) # shape parameter
  
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

# Example usage
cat("=== Simulating scRNA-seq data with Splatter ===\n")

# Simulation with moderate effects
sim_data <- simulate_scrna_splatter(
  n_cells_per_condition = 250,
  n_genes = 1000,
  de_prob = 0.15,                    # 15% of genes are DE
  group_effect_size = 1.8,           # 1.8x fold change for group effect
  time_effect_size = 1.5,            # 1.5x fold change for time effect  
  interaction_effect_size = 1.3,     # 1.3x additional effect for interaction
  seed = 123
)

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