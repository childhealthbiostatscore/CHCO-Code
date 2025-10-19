# Load package
library(splatter)
library(scater)

# specify user for paths
user <- Sys.info()[["user"]]
if (user == "laurapyle") {
  data_path <- "/Users/laurapyle/Library/CloudStorage/OneDrive-UW/Pyle/scRNAseq simulations"
  github_path <- "/Users/laurapyle/Documents/GitHub/CHCO-Code/Laura Pyle/WNAR 2025"
} else if (user == "lpyle") {
  data_path <- "/Users/lpyle/Library/CloudStorage/OneDrive-UW/Pyle/scRNAseq simulations"
  github_path <- "/Users/lpyle/Documents/GitHub/CHCO-Code/Laura Pyle/WNAR 2025"
} else if (user == "pylell") {
  data_path <- "/Users/pylell/Library/CloudStorage/OneDrive-UW/Pyle/scRNAseq simulations"
  github_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Laura Pyle/WNAR 2025"
} else {
  stop("Unknown user: please specify root path for this user.")
}

filename <- "Simulate cohort and fit model/complete_simulation_pipeline.r"
full_path <- file.path(github_path, filename)
source(full_path)

#############################################
# Example 1: Basic pipeline with default DE #
#############################################

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
  output_dir = paste0(data_path,"/results/de_power"),
  batch_size = 50,
  parallel = TRUE,
  n_cores = 8
)
analyze_pipeline_results(output_dir = paste0(data_path,"/results/de_power"))
print(analysis$plot)

#######################################
# Example 2: Compare multiple methods #
#######################################

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
  output_dir = paste0(data_path,"/results/method_comparison")
)

# STEP 4: NOW use the analysis functions from the second script
filename <- "Simulate cohort and fit model/method_comparison_analysis.R"
full_path <- file.path(github_path, filename)
source(full_path)

# Check structure first
check_results_structure(paste0(data_path,"/results/method_comparison")) 

# Then analyze
analysis <- analyze_method_comparison(paste0(data_path,"/results/method_comparison"))
print(analysis$plots$combined)

# Get statistics
stats <- compare_methods_statistically(results)
print(stats$overall_metrics)

#####################################################
# Example 3: Large-scale sweep with minimal storage #
#####################################################

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



