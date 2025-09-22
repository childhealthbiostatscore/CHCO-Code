# Load package
library(splatter)
library(scater)

source("/Users/laurapyle/Documents/GitHub/CHCO-Code/Laura Pyle/WNAR 2025/splatter_sim_cohort_functions_v6.R")

# ONE SIMULATION WITH SET PARAMETERS
sim_cohort <- simulate_patient_control_study(
    n_controls = 10,
    n_patients = 10,
    cells_per_ind = 150,
    cells_sd = 150,           # SD for cell count variation
    cells_min = 30,            # Minimum cells per individual
    cells_max = NULL,          # Maximum cells per individual
    n_genes = 3000,
    individual_variation = 0.1,
    disease_effect_size = 2,
    disease_gene_fraction = 0.2,
    within_ind_correlation = 0.1,
    n_cell_states = 5,
    verbose = TRUE,
    seed = 123)

# Check distribution by disease status
table(colData(sim_cohort)$Group, colData(sim_cohort)$Individual)
plot_expression_violin(sim_cohort, group_var = "Group")

p1 <- plot_disease_genes(sim_cohort, n_genes = 50)
print(p1)

# SIMULATE OVER A RANGE OF PARAMETERS
# PSXXX in the simulation ID variable indicates a unique set of parameter combinations
sweep_exp1 <- simulate_parameter_sweep(
  experiment_name = "dose_response_study",
  n_controls = 5,
  n_patients = 5,
  disease_effect_size = c(1, 2, 5),
  individual_variation = c(0.1, 0.3),
  n_simulations = 3,  # 3 replicates
  verbose = TRUE,
  return_summary = TRUE
)

# View the unique IDs
head(sweep_exp1$detailed[, c("simulation_id", "param_set_id", 
                             "simulation_rep", "param_hash")])

# Example 2: Return full simulations with IDs
sweep_sims <- simulate_parameter_sweep(
  experiment_name = "cell_variation_test",
  cells_per_ind = c(100, 200),
  disease_effect_size = c(2, 5),
  n_simulations = 2,
  return_summary = FALSE
)

# Get specific simulation
sim1 <- get_simulation_by_id(sweep_sims, "cell_variation_test_ps002_rep01")
metadata(sim1)$simulation_id
table(colData(sim1)$simulation_id)

# Example 3: Track and compare specific simulations
sweep_test <- simulate_parameter_sweep(
  experiment_name = "parameter_test",
  disease_effect_size = c(1, 5, 10),
  individual_variation = c(0.1, 0.5),
  n_simulations = 2,
  return_summary = TRUE
)

# Find best performing simulations
library(dplyr)
best_sims <- sweep_test$detailed %>%
  arrange(desc(mean_log_fc_disease)) %>%
  head(3) %>%
  select(simulation_id, disease_effect_size, 
         individual_variation, mean_log_fc_disease)

print(best_sims)

# Compare specific simulations
compare_simulations(sweep_test, 
                    c("parameter_test_ps001_rep01", 
                      "parameter_test_ps006_rep01"))

# Example 4: Save and reload with preserved IDs
# Save results
saveRDS(sweep_test, "parameter_sweep_results.rds")
loaded_sweep <- readRDS("parameter_sweep_results.rds")

# Plot one simulation
sim_to_plot <- sweep_sims[[5]]
p3 <- plot_expression_violin(
  sim_to_plot, 
  group_var = "Group",
  plot_title = paste("Simulation:", metadata(sim_to_plot)$simulation_id)
)
print(p3)
