# Test function to debug NEBULA
test_nebula_on_single_sim <- function() {
  # Create a simple test simulation
  sim <- simulate_patient_control_study(
    n_controls = 5,
    n_patients = 5,
    cells_per_ind = 100,
    disease_effect_size = 5,
    verbose = TRUE
  )
  
  # Try to fit nebula
  cat("\nTesting NEBULA...\n")
  result <- fit_nebula_model(sim)
  
  print(result)
  
  # Check what went wrong
  if(is.na(result$n_genes_nebula)) {
    cat("\nNEBULA failed. Checking data structure:\n")
    cat("Unique individuals:", length(unique(colData(sim)$Individual)), "\n")
    cat("Unique groups:", unique(colData(sim)$Group), "\n")
    cat("Cells per individual:\n")
    print(table(colData(sim)$Individual))
  }
  
  return(result)
}

# Run the test
test_result <- test_nebula_on_single_sim()

# test the fixed functions
# Test with one simulation
test_sim <- simulate_patient_control_study(
  n_controls = 5,
  n_patients = 5,
  cells_per_ind = 100,
  disease_effect_size = 5
)

cat("Disease genes in simulation:", sum(rowData(test_sim)$is_disease_gene), "\n")

# Test each method
de_result <- fit_default_de_model(test_sim)
pb_result <- fit_pseudobulk_model(test_sim)
nebula_result <- fit_nebula_model(test_sim)

cat("\nPower results:\n")
cat("DE power:", de_result$power_005, "\n")
cat("PB power:", pb_result$power_pb_005, "\n")
cat("NEBULA power:", nebula_result$power_nebula_005, "\n")
