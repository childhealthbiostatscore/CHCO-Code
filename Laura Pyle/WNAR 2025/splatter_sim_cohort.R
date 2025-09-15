# Load package
library(splatter)
library(scater)

source("/Users/laurapyle/Documents/GitHub/CHCO-Code/Laura Pyle/WNAR 2025/splatter_sim_cohort_functions.R")

# Different patterns for different cohorts
simulate_cohort_specific <- function(n_controls = 5, 
                                     n_patients = 8,
                                     control_cells_mean = 180,
                                     patient_cells_mean = 120,
                                     control_facLoc = 0,
                                     control_facScale = 0.15,
                                     patient_facLoc = 0.1,
                                     patient_facScale = 0.3,
                                     correlation = 0.4,
                                     n_genes = 3000,
                                     verbose = TRUE) {
  
  set.seed(42)
  
  # Generate cell counts
  control_counts <- round(rnorm(n_controls, control_cells_mean, 30))
  control_counts[control_counts < 30] <- 30
  
  patient_counts <- round(rnorm(n_patients, patient_cells_mean, 50))
  patient_counts[patient_counts < 30] <- 30
  
  # Combine counts
  all_counts <- c(control_counts, patient_counts)
  
  # Create vectors of batch parameters
  all_facLoc <- c(rep(control_facLoc, n_controls), 
                  rep(patient_facLoc, n_patients))
  
  all_facScale <- c(rep(control_facScale, n_controls),
                    rep(patient_facScale, n_patients))
  
  if(verbose) {
    cat("=== Cohort Design ===\n")
    cat("Controls: n =", n_controls, "\n")
    cat("  Cells (mean):", control_cells_mean, "\n")
    cat("  batch.facLoc:", control_facLoc, "\n")
    cat("  batch.facScale:", control_facScale, "\n\n")
    
    cat("Patients: n =", n_patients, "\n")
    cat("  Cells (mean):", patient_cells_mean, "\n")
    cat("  batch.facLoc:", patient_facLoc, "\n")
    cat("  batch.facScale:", patient_facScale, "\n\n")
  }
  
  # Generate simulation with group-specific parameters
  sim <- simulate_with_correlation(
    n_individuals = n_controls + n_patients,
    cells_per_ind = all_counts,
    n_genes = n_genes,
    batch_facLoc = all_facLoc,      # Different for each group
    batch_facScale = all_facScale,  # Different for each group
    within_ind_correlation = correlation,
    verbose = verbose
  )
  
  # Add disease status
  disease_status <- c(rep("Control", n_controls), 
                      rep("CKM", n_patients))
  
  # Map to cells
  ind_number <- as.integer(factor(colData(sim)$Individual))
  colData(sim)$DiseaseStatus <- factor(disease_status[ind_number])
  
  # Update individual labels
  individual_labels <- c(paste0("Control_", 1:n_controls),
                         paste0("Patient_", 1:n_patients))
  colData(sim)$Individual <- factor(individual_labels[ind_number])
  
  if(verbose) {
    cat("\n=== Final Cohort Summary ===\n")
    print(table(colData(sim)$DiseaseStatus, colData(sim)$Individual))
  }
  
  return(sim)
}

# this doesn't seem to allow for different parameter estimates in the groups
sim_cohort <- simulate_cohort_specific(
  n_controls = 4,
  n_patients = 6,
  control_cells_mean = 180,
  patient_cells_mean = 120,
  control_facLoc = 0,
  control_facScale = 0.15,
  patient_facLoc = 10,
  patient_facScale = 10,
  correlation = 0.4,
  n_genes = 3000,
  verbose = TRUE)

# Check distribution by disease status
table(colData(sim_cohort)$DiseaseStatus, colData(sim_cohort)$Individual)
plot_expression_violin(sim_cohort, group_var = "DiseaseStatus")



# remaining workflow
# visualize simulation
# simulate remaining data and save
