# Load package
library(splatter)
library(scater)

# Create mock data
set.seed(1)
sce <- mockSCE()

# Estimate parameters from mock data
params <- splatEstimate(sce)
# Simulate data using estimated parameters
sim <- splatSimulate(params)

# Population level simulation

# TEST SIM 1
# Define individual-level parameters
n_individuals <- 12
mean_cells_per_ind <- 150
sd_cells_per_ind <- 40
# Generate realistic cell counts per individual
set.seed(42)
cells_per_individual <- round(abs(rnorm(n_individuals, 
                                        mean = mean_cells_per_ind, 
                                        sd = sd_cells_per_ind)))
cells_per_individual[cells_per_individual < 30] <- 30  # Minimum threshold
# Create parameters with individual-level variation
params <- newSplatParams(
  nGenes = 3000,
  batchCells = cells_per_individual,
  batch.facLoc = 0.1,      # Log-normal mean for individual effects
  batch.facScale = 0.25,   # Log-normal SD for individual effects
  mean.shape = 0.6,
  mean.rate = 0.3,
  dropout.mid = 2,
  dropout.shape = -1
)
# Simulate
sim_individuals1 <- splatSimulate(params, verbose = FALSE)
# Add clear individual labels
colData(sim_individuals1)$Individual <- factor(
  colData(sim_individuals1)$Batch,
  labels = paste0("Patient_", 1:n_individuals)
)

# TEST SIM 2
# Define individual-level parameters
n_individuals <- 12
mean_cells_per_ind <- 150
sd_cells_per_ind <- 40
# Generate realistic cell counts per individual
set.seed(42)
cells_per_individual <- round(abs(rnorm(n_individuals, 
                                        mean = mean_cells_per_ind, 
                                        sd = sd_cells_per_ind)))
cells_per_individual[cells_per_individual < 30] <- 30  # Minimum threshold
# Create parameters with individual-level variation
params <- newSplatParams(
  nGenes = 3000,
  batchCells = cells_per_individual,
  batch.facLoc = 5,      # Log-normal mean for individual effects
  batch.facScale = 0.25,   # Log-normal SD for individual effects
  mean.shape = 0.6,
  mean.rate = 0.3,
  dropout.mid = 2,
  dropout.shape = -1
)
# Simulate
sim_individuals2 <- splatSimulate(params, verbose = FALSE)
# Add clear individual labels
colData(sim_individuals2)$Individual <- factor(
  colData(sim_individuals2)$Batch,
  labels = paste0("Patient_", 1:n_individuals)
)
comparison <- compareSCEs(list(Splat = sim_individuals1, Simple = sim_individuals2))
names(comparison$Plots)
comparison$Plots$Means
