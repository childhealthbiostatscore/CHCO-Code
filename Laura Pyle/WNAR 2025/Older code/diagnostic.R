# DIAGNOSTIC VERSION: Debug why NEBULA vs Pseudobulk results flipped
# This version adds extensive debugging to identify the issue

# Diagnostic function to compare a simple case
diagnostic_comparison <- function() {
  
  cat("=== DIAGNOSTIC: Investigating NEBULA vs Pseudobulk Issue ===\n")
  
  # Simple test simulation
  sim_data <- simulate_scrna_imbalanced(
    n_cells_per_condition = 200,
    cell_imbalance_ratio = c(1, 1),  # Balanced
    n_genes = 100,
    group_effect_size = 2.0,  # Large effect
    de_prob = 0.3,  # Many DE genes
    seed = 42
  )
  
  cat("True DE genes:", length(sim_data$de_genes), "\n")
  
  # Test both methods with debugging
  cat("\n--- Testing Pseudobulk ---\n")
  pseudo_result <- debug_pseudobulk_deseq2(sim_data)
  
  cat("\n--- Testing NEBULA ---\n")
  nebula_result <- debug_nebula_mixed_effects(sim_data)
  
  # Compare results
  cat("\n=== COMPARISON RESULTS ===\n")
  cat("Pseudobulk:\n")
  cat("  Significant genes:", pseudo_result$n_significant, "\n")
  cat("  True positives:", pseudo_result$true_positives, "\n")
  cat("  Power:", round(pseudo_result$power, 3), "\n")
  cat("  FDR:", round(pseudo_result$fdr, 3), "\n")
  
  cat("NEBULA:\n")
  cat("  Significant genes:", nebula_result$n_significant, "\n")
  cat("  True positives:", nebula_result$true_positives, "\n")
  cat("  Power:", round(nebula_result$power, 3), "\n")
  cat("  FDR:", round(nebula_result$fdr, 3), "\n")
  
  # Check p-value distributions
  cat("\n=== P-VALUE DIAGNOSTICS ===\n")
  cat("Pseudobulk p-values summary:\n")
  print(summary(pseudo_result$p_values))
  cat("NEBULA p-values summary:\n")
  print(summary(nebula_result$p_values))
  
  # Count how many p-values are exactly 1 (indicating failure)
  cat("Pseudobulk p-values = 1:", sum(pseudo_result$p_values == 1), "out of", length(pseudo_result$p_values), "\n")
  cat("NEBULA p-values = 1:", sum(nebula_result$p_values == 1), "out of", length(nebula_result$p_values), "\n")
  
  return(list(
    sim_data = sim_data,
    pseudobulk = pseudo_result,
    nebula = nebula_result
  ))
}

# Debug version of pseudobulk with extensive logging
debug_pseudobulk_deseq2 <- function(sim_data, alpha = 0.05) {
  
  cat("Starting pseudobulk analysis...\n")
  
  if(!requireNamespace("DESeq2", quietly = TRUE)) {
    cat("DESeq2 not available - simulating results\n")
    return(create_dummy_result("pseudobulk", sim_data, success = FALSE))
  }
  
  tryCatch({
    
    cell_meta <- sim_data$cell_metadata
    counts_matrix <- sim_data$counts
    
    # Filter to T1 timepoint
    t1_cells <- cell_meta$Time == "T1"
    t1_meta <- cell_meta[t1_cells, ]
    t1_counts <- counts_matrix[, t1_cells]
    
    cat("T1 cells: Control =", sum(t1_meta$Group == "Control"), 
        ", Treatment =", sum(t1_meta$Group == "Treatment"), "\n")
    
    # Create pseudobulk samples
    n_reps <- 6
    cat("Creating", n_reps, "pseudobulk replicates...\n")
    
    pseudobulk_counts <- matrix(0, nrow = nrow(t1_counts), ncol = n_reps)
    pseudobulk_meta <- data.frame(
      sample_id = paste0("Rep_", 1:n_reps),
      group = rep(c("Control", "Treatment"), each = n_reps/2),
      stringsAsFactors = FALSE
    )
    rownames(pseudobulk_counts) <- rownames(t1_counts)
    colnames(pseudobulk_counts) <- pseudobulk_meta$sample_id
    
    # Assign cells to replicates
    control_cells <- which(t1_meta$Group == "Control")
    treatment_cells <- which(t1_meta$Group == "Treatment")
    
    if(length(control_cells) < 3 || length(treatment_cells) < 3) {
      cat("ERROR: Not enough cells per group\n")
      return(create_dummy_result("pseudobulk", sim_data, success = FALSE))
    }
    
    control_reps <- split(sample(control_cells), rep(1:(n_reps/2), length.out = length(control_cells)))
    treatment_reps <- split(sample(treatment_cells), rep(1:(n_reps/2), length.out = length(treatment_cells)))
    
    # Sum counts for each replicate
    for(i in 1:(n_reps/2)) {
      if(length(control_reps[[i]]) > 0) {
        pseudobulk_counts[, i] <- rowSums(t1_counts[, control_reps[[i]], drop = FALSE])
      }
      if(length(treatment_reps[[i]]) > 0) {
        pseudobulk_counts[, i + n_reps/2] <- rowSums(t1_counts[, treatment_reps[[i]], drop = FALSE])
      }
    }
    
    cat("Pseudobulk counts summary:\n")
    print(summary(as.vector(pseudobulk_counts)))
    
    # Filter genes
    keep_genes <- rowSums(pseudobulk_counts >= 5) >= 2
    cat("Keeping", sum(keep_genes), "out of", length(keep_genes), "genes\n")
    
    if(sum(keep_genes) < 10) {
      cat("ERROR: Too few genes passing filter\n")
      return(create_dummy_result("pseudobulk", sim_data, success = FALSE))
    }
    
    pseudobulk_counts <- pseudobulk_counts[keep_genes, ]
    
    # Run DESeq2
    cat("Running DESeq2...\n")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = pseudobulk_counts,
      colData = pseudobulk_meta,
      design = ~ group
    )
    
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c("group", "Treatment", "Control"))
    
    cat("DESeq2 completed. Results summary:\n")
    print(summary(res))
    
    # Extract results
    p_values <- rep(1, nrow(sim_data$counts))
    names(p_values) <- rownames(sim_data$counts)
    
    # Fill in p-values for tested genes
    tested_genes <- rownames(res)[!is.na(res$pvalue)]
    p_values[tested_genes] <- res$pvalue[!is.na(res$pvalue)]
    
    cat("Tested", length(tested_genes), "genes\n")
    cat("P-value range:", range(p_values[tested_genes]), "\n")
    
    # Multiple testing correction
    p_adj <- p.adjust(p_values, method = "BH")
    significant_genes <- which(p_adj < alpha)
    
    cat("Significant genes at alpha =", alpha, ":", length(significant_genes), "\n")
    
  }, error = function(e) {
    cat("Pseudobulk ERROR:", e$message, "\n")
    return(create_dummy_result("pseudobulk", sim_data, success = FALSE))
  })
  
  # Calculate performance
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
    fdr = false_positives / max(length(significant_genes), 1),
    success = TRUE
  ))
}

# Debug version of NEBULA with extensive logging
debug_nebula_mixed_effects <- function(sim_data, alpha = 0.05) {
  
  cat("Starting NEBULA analysis...\n")
  
  if(!requireNamespace("nebula", quietly = TRUE)) {
    cat("NEBULA not available - simulating results\n")
    return(create_dummy_result("nebula", sim_data, success = FALSE))
  }
  
  tryCatch({
    
    cell_meta <- sim_data$cell_metadata
    counts_matrix <- sim_data$counts
    
    # Filter to T1 timepoint
    t1_cells <- cell_meta$Time == "T1"
    t1_meta <- cell_meta[t1_cells, ]
    t1_counts <- counts_matrix[, t1_cells]
    
    cat("T1 cells: Control =", sum(t1_meta$Group == "Control"), 
        ", Treatment =", sum(t1_meta$Group == "Treatment"), "\n")
    
    # Create subject IDs - THIS MIGHT BE THE BUG
    n_subjects_per_group <- 8  # Increased from 6
    
    control_indices <- which(t1_meta$Group == "Control")
    treatment_indices <- which(t1_meta$Group == "Treatment")
    
    cat("Creating subject assignments...\n")
    cat("Control cells:", length(control_indices), "\n")
    cat("Treatment cells:", length(treatment_indices), "\n")
    
    # Better subject assignment strategy
    subject_assignment <- character(nrow(t1_meta))
    
    # Assign control cells to subjects more evenly
    if(length(control_indices) > 0) {
      control_subjects <- paste0("Control_S", 1:n_subjects_per_group)
      subject_assignment[control_indices] <- rep(control_subjects, 
                                                 length.out = length(control_indices))
    }
    
    # Assign treatment cells to subjects more evenly  
    if(length(treatment_indices) > 0) {
      treatment_subjects <- paste0("Treatment_S", 1:n_subjects_per_group)
      subject_assignment[treatment_indices] <- rep(treatment_subjects, 
                                                   length.out = length(treatment_indices))
    }
    
    t1_meta$subject_id <- subject_assignment
    t1_meta$group_numeric <- ifelse(t1_meta$Group == "Treatment", 1, 0)
    
    # Check subject distribution
    cat("Subject distribution:\n")
    print(table(t1_meta$subject_id, t1_meta$Group))
    
    # Filter genes more conservatively for NEBULA
    gene_means <- rowMeans(t1_counts)
    gene_detection <- rowSums(t1_counts > 0)
    
    # More lenient filtering for NEBULA
    keep_genes <- gene_means > 0.01 & gene_detection >= 5
    
    cat("Gene filtering: mean >", 0.01, ", detected in >=", 5, "cells\n")
    cat("Keeping", sum(keep_genes), "out of", length(keep_genes), "genes\n")
    
    if(sum(keep_genes) < 10) {
      cat("ERROR: Too few genes for NEBULA\n")
      return(create_dummy_result("nebula", sim_data, success = FALSE))
    }
    
    filtered_counts <- t1_counts[keep_genes, ]
    
    cat("Filtered counts summary:\n")
    print(summary(as.vector(filtered_counts)))
    
    # Create NEBULA data object
    cat("Creating NEBULA data object...\n")
    data_obj <- nebula::group_cell(
      count = filtered_counts,
      id = t1_meta$subject_id,
      pred = data.frame(
        group = t1_meta$group_numeric,
        intercept = 1
      )
    )
    
    cat("NEBULA data object created successfully\n")
    cat("Data dimensions:", dim(data_obj$count), "\n")
    cat("Number of subjects:", length(unique(data_obj$id)), "\n")
    
    # Run NEBULA-HL
    cat("Running NEBULA-HL...\n")
    res <- nebula::nebula(
      data_obj$count,
      data_obj$id, 
      pred = data_obj$pred,
      offset = data_obj$offset,
      model = "NBLMM"
    )
    
    cat("NEBULA completed successfully\n")
    cat("Results summary:\n")
    if(!is.null(res$summary) && nrow(res$summary) > 0) {
      print(summary(res$summary))
      cat("Number of results:", nrow(res$summary), "\n")
      cat("P-value range:", range(res$summary$p_group, na.rm = TRUE), "\n")
    } else {
      cat("ERROR: No results from NEBULA\n")
      return(create_dummy_result("nebula", sim_data, success = FALSE))
    }
    
    # Extract p-values
    p_values <- rep(1, nrow(sim_data$counts))
    names(p_values) <- rownames(sim_data$counts)
    
    # Fill in p-values for tested genes
    if(nrow(res$summary) > 0) {
      tested_genes <- rownames(res$summary)
      valid_pvals <- !is.na(res$summary$p_group)
      
      if(sum(valid_pvals) > 0) {
        p_values[tested_genes[valid_pvals]] <- res$summary$p_group[valid_pvals]
        cat("Successfully extracted", sum(valid_pvals), "p-values from NEBULA\n")
      } else {
        cat("ERROR: All NEBULA p-values are NA\n")
        return(create_dummy_result("nebula", sim_data, success = FALSE))
      }
    }
    
    # Multiple testing correction
    p_adj <- p.adjust(p_values, method = "BH")
    significant_genes <- which(p_adj < alpha)
    
    cat("Significant genes at alpha =", alpha, ":", length(significant_genes), "\n")
    
  }, error = function(e) {
    cat("NEBULA ERROR:", e$message, "\n")
    traceback()
    return(create_dummy_result("nebula", sim_data, success = FALSE))
  })
  
  # Calculate performance
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
    fdr = false_positives / max(length(significant_genes), 1),
    success = TRUE
  ))
}

# Helper function to create dummy results when methods fail
create_dummy_result <- function(method_name, sim_data, success = FALSE) {
  p_values <- rep(1, nrow(sim_data$counts))
  names(p_values) <- rownames(sim_data$counts)
  p_adj <- rep(1, nrow(sim_data$counts))
  names(p_adj) <- rownames(sim_data$counts)
  
  return(list(
    method = method_name,
    p_values = p_values,
    p_adj = p_adj,
    significant_genes = integer(0),
    n_significant = 0,
    true_positives = 0,
    false_positives = 0,
    power = 0,
    fdr = 0,
    success = success
  ))
}

# Test if packages are available
test_package_availability <- function() {
  cat("=== PACKAGE AVAILABILITY TEST ===\n")
  
  deseq2_available <- requireNamespace("DESeq2", quietly = TRUE)
  nebula_available <- requireNamespace("nebula", quietly = TRUE)
  
  cat("DESeq2 available:", deseq2_available, "\n")
  cat("NEBULA available:", nebula_available, "\n")
  
  if(!deseq2_available) {
    cat("Install DESeq2 with: BiocManager::install('DESeq2')\n")
  }
  
  if(!nebula_available) {
    cat("Install NEBULA with: install.packages('nebula')\n")
  }
  
  return(list(deseq2 = deseq2_available, nebula = nebula_available))
}

cat("=== DIAGNOSTIC VERSION LOADED ===\n")
cat("Run these functions to debug:\n")
cat("1. test_package_availability()  # Check if packages are installed\n")
cat("2. diagnostic_comparison()      # Run detailed comparison with debugging\n")

# Run the diagnostic
cat("\n=== RUNNING DIAGNOSTICS ===\n")
pkg_status <- test_package_availability()

if(pkg_status$deseq2 && pkg_status$nebula) {
  cat("\nBoth packages available - running diagnostic comparison...\n")
  diagnostic_results <- diagnostic_comparison()
} else {
  cat("\nMissing packages - cannot run full diagnostic\n")
  if(!pkg_status$deseq2) cat("Missing: DESeq2\n")
  if(!pkg_status$nebula) cat("Missing: NEBULA\n")
}