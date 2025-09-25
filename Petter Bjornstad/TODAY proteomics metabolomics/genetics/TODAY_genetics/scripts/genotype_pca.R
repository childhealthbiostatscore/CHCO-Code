suppressPackageStartupMessages(library(reticulate    ))
reticulate::use_python("/mmfs1/gscratch/scrubbed/mdanto/envs/bioinfo/bin/python") 
use_condaenv("bioinfo", required = TRUE)

Sys.setenv(PATH = paste("/mmfs1/gscratch/scrubbed/mdanto/envs/bioinfo/bin", Sys.getenv("PATH"), sep = ":"))

source("R/functions.R")
source("/mmfs1/gscratch/togo/matteo/projects/software/R_scripts/plink_pca_wrapper.r")

# Run PCA on imputed dataset
ref_prefix    = "/mmfs1/gscratch/togo/matteo/projects/data/1kgp_reference_imputation/1kgp_plink2"
target_prefix = here("data/genotyping/imputation/plink2/merged_imputed-merge")

results <- fast_genotype_pca(
  ref_prefix    = ref_prefix,
  target_prefix = target_prefix,
  output_dir    = here("output/genotype_pca"),
  n_pcs = 20,
  n_snps = 50000,
  threads = 32,       
  memory = 256000 
)


