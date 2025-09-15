# Load package
suppressPackageStartupMessages({
  library(splatter)
  library(scater)
})

# Create mock data
set.seed(1)
sce <- mockSCE()

# Estimate parameters from mock data
params <- splatEstimate(sce)
# Simulate data using estimated parameters
sim <- splatSimulate(params)