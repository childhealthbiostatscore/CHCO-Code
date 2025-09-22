# Load package
library(splatter)
library(scater)

source("/Users/laurapyle/Documents/GitHub/CHCO-Code/Laura Pyle/WNAR 2025/splatter_sim_cohort_functions_v4.R")

sim_cohort <- simulate_patient_control_study(n_controls = 10,
                               n_patients = 10,
                               cells_per_ind = 100,
                               n_genes = 3000,
                               individual_variation = 0.1,
                               disease_effect_size = 10,
                               disease_gene_fraction = 0.1,
                               verbose = TRUE)

# Check distribution by disease status
table(colData(sim_cohort)$Group, colData(sim_cohort)$Individual)
plot_expression_violin(sim_cohort, group_var = "Group")

p1 <- plot_disease_genes(sim_cohort, n_genes = 50)
print(p1)

# remaining workflow
# visualize simulation
# simulate remaining data and save
