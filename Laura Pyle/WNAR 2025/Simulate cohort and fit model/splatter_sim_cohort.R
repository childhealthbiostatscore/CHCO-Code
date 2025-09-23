# Load package
library(splatter)
library(scater)

source('/Users/laurapyle/Documents/GitHub/CHCO-Code/Laura Pyle/WNAR 2025/Simulate cohort and fit model/splatter sim and model functions.r')

# Example 1: Basic pipeline with default DE
results_basic <- run_simulation_analysis_pipeline(
  param_grid = expand.grid(
    n_controls = 10,
    n_patients = 10,
    cells_per_ind = c(100, 200),
    disease_effect_size = c(2, 5, 10),
    individual_variation = c(0.1, 0.3),
    disease_gene_fraction = 0.2
  ),
  n_simulations = 10,
  experiment_name = "de_power_analysis",
  fit_function = fit_default_de_model,
  output_dir = "results/de_power",
  batch_size = 50,
  parallel = TRUE,
  n_cores = 8
)

# Example 2: Compare multiple methods
fit_multiple_methods <- function(sim, params) {
  # Run multiple analyses
  de_results <- fit_default_de_model(sim, params)
  pb_results <- fit_pseudobulk_model(sim, params)
  
  # Combine
  cbind(de_results, pb_results)
}

results_comparison <- run_simulation_analysis_pipeline(
  param_grid = expand.grid(
    disease_effect_size = c(1, 2, 5, 10),
    individual_variation = c(0.1, 0.2, 0.3)
  ),
  n_simulations = 20,
  fit_function = fit_multiple_methods,
  output_dir = "results/method_comparison"
)

# Example 3: Large-scale sweep with minimal storage
results_large <- run_simulation_analysis_pipeline(
  param_grid = expand.grid(
    cells_per_ind = seq(50, 300, by = 50),
    disease_effect_size = seq(1, 10, by = 1),
    individual_variation = seq(0.1, 0.5, by = 0.1),
    disease_gene_fraction = c(0.05, 0.1, 0.2, 0.3)
  ),
  n_simulations = 5,
  experiment_name = "large_parameter_sweep",
  save_intermediate = FALSE,  # Don't save simulations
  batch_size = 200,
  parallel = TRUE
)

analysis <- analyze_pipeline_results("results/de_power")
print(analysis$plot)
