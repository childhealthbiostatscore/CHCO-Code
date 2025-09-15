# Load package
library(splatter)
library(scater)

source("/Users/laurapyle/Documents/GitHub/CHCO-Code/Laura Pyle/WNAR 2025/splatter_sim_cohort_functions.R")

# Different patterns for different cohorts
simulate_cohort_specific <- function(n_controls = 5, 
                                     n_patients = 8,
                                     control_cells_mean = 180,
                                     patient_cells_mean = 120,
                                     correlation = 0.4) {
  
  # Controls often have better cell yield
  control_counts <- round(rnorm(n_controls, control_cells_mean, 30))
  
  # Patients might have more variable/lower yields
  patient_counts <- round(rnorm(n_patients, patient_cells_mean, 50))
  
  # Combine
  all_counts <- c(control_counts, patient_counts)
  
  # Generate simulation
  sim <- simulate_with_correlation(
    n_individuals = n_controls + n_patients,
    cells_per_ind = all_counts,
    within_ind_correlation = correlation
  )
  
  # Add disease status
  disease_status <- c(rep("Control", n_controls), 
                      rep("Treatment", n_patients))
  
  # Fix: Extract batch number from "Batch1", "Batch2", etc.
  batch_numbers <- as.numeric(gsub("Batch", "", colData(sim)$Batch))
  batch_to_status <- disease_status[batch_numbers]
  colData(sim)$DiseaseStatus <- factor(batch_to_status)
  
  return(sim)
}

sim_cohort <- simulate_cohort_specific(
  n_controls = 4,
  n_patients = 6
)

# Check distribution by disease status
table(colData(sim_cohort)$DiseaseStatus, colData(sim_cohort)$Individual)
control_cells <-  colData(sim_cohort)$DiseaseStatus == "Control"
tx_cells <- colData(sim_cohort)$DiseaseStatus == "Treatment"
sce_control <- sim_cohort[, control_cells]
sce_tx <- sim_cohort[, tx_cells]
sce_list <- list()
sce_list[["Control"]] <- sce_control
sce_list[["CKM"]] <- tx_cells

comparison <- compareSCEs(sce_list)
names(comparison$Plots)
comparison$Plots$Means