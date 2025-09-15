# FUNCTION TO ADD CORRELATION TO CELLS WITHIN A PERSON
add_cell_state_correlation <- function(sim, 
                                       n_states = 3,
                                       state_correlation = 0.5) {
  
  counts_original <- counts(sim)
  
  for(batch in unique(colData(sim)$Batch)) {
    batch_cells <- which(colData(sim)$Batch == batch)
    n_cells <- length(batch_cells)
    
    # Assign cells to states (some overlap creates correlation)
    if(state_correlation > 0) {
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
      state_programs <- matrix(rnorm(nrow(counts_original) * n_states, 0, 0.2),
                               nrow = nrow(counts_original))
      
      # Apply weighted combination of programs
      for(i in 1:n_cells) {
        cell_program <- state_programs %*% state_weights[i, ]
        counts_original[, batch_cells[i]] <- 
          counts_original[, batch_cells[i]] * exp(cell_program)
      }
    }
  }
  
  assay(sim, "counts") <- round(counts_original)
  return(sim)
}

#sim_states <- add_cell_state_correlation(sim, 
#                                         n_states = 4,
 #                                        state_correlation = 0.6)

# BASE FUNCTION TO DO SIMULATION
simulate_with_correlation <- function(n_individuals = 8,
                                          cells_per_ind = 150,
                                          cells_sd = NULL,
                                          cells_min = 30,
                                          n_genes = 3000,
                                          batch_facLoc = 0.1,
                                          batch_facScale = 0.25,
                                          within_ind_correlation = 0.4,
                                          n_cell_states = 5,
                                          verbose = TRUE) {
  
  # Handle cell count variability
  if(length(cells_per_ind) == 1) {
    # Single value provided - generate variable counts
    if(is.null(cells_sd)) {
      # Default to 30% coefficient of variation if sd not specified
      cells_sd <- cells_per_ind * 0.3
    }
    
    # Generate variable cell counts
    set.seed(123)  # For reproducibility
    cell_counts <- round(rnorm(n_individuals, 
                               mean = cells_per_ind, 
                               sd = cells_sd))
    
    # Apply minimum threshold
    cell_counts[cell_counts < cells_min] <- cells_min
    
  } else if(length(cells_per_ind) == n_individuals) {
    # Vector of counts provided - use directly
    cell_counts <- cells_per_ind
    
    # Still apply minimum threshold
    cell_counts[cell_counts < cells_min] <- cells_min
    
  } else {
    stop("cells_per_ind must be either a single value or a vector of length n_individuals")
  }
  
  # Print summary if verbose
  if(verbose) {
    cat("=== Simulation Parameters ===\n")
    cat("Genes:", n_genes, "\n")
    cat("Individuals:", n_individuals, "\n")
    cat("Batch effect location:", batch_facLoc, "\n")
    cat("Batch effect scale:", batch_facScale, "\n")
    cat("Within-individual correlation:", within_ind_correlation, "\n")
    cat("Cell states:", n_cell_states, "\n\n")
    
    cat("Cell distribution across individuals:\n")
    cat("  Mean:", round(mean(cell_counts), 1), "\n")
    cat("  SD:", round(sd(cell_counts), 1), "\n")
    cat("  Range:", min(cell_counts), "-", max(cell_counts), "\n")
    cat("  Total cells:", sum(cell_counts), "\n\n")
  }
  
  # Base simulation with variable parameters
  sim <- splatSimulate(
    nGenes = n_genes,           # Now customizable
    batchCells = cell_counts,
    batch.facLoc = batch_facLoc,    # Now customizable
    batch.facScale = batch_facScale, # Now customizable
    verbose = FALSE
  )
  
  # Add correlation through shared cell states
  sim <- add_cell_state_correlation(sim, 
                                    n_states = n_cell_states,
                                    state_correlation = within_ind_correlation)
  
  # Add metadata with actual cell counts
  colData(sim)$Individual <- factor(colData(sim)$Batch,
                                    labels = paste0("Patient_", 1:n_individuals))
  
  # Add cell count info to metadata
  ind_cell_counts <- table(colData(sim)$Individual)
  colData(sim)$CellsPerIndividual <- ind_cell_counts[colData(sim)$Individual]
  
  # Add simulation parameters as metadata
  metadata(sim) <- list(
    n_genes = n_genes,
    n_individuals = n_individuals,
    batch_facLoc = batch_facLoc,
    batch_facScale = batch_facScale,
    within_ind_correlation = within_ind_correlation,
    n_cell_states = n_cell_states,
    total_cells = sum(cell_counts)
  )
  
  return(sim)
}
