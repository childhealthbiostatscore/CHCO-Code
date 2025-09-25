source("R/functions.R")

################################################################################
# Get input files

suppressWarnings(file.symlink("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/data/genotyping/imputation/plink2/merged_imputed-merge.pgen", here("data/to_qtl/genotypes.pgen")))
suppressWarnings(file.symlink("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/data/genotyping/imputation/plink2/merged_imputed-merge.pvar", here("data/to_qtl/genotypes.pvar")))
suppressWarnings(file.symlink("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/data/genotyping/imputation/plink2/merged_imputed-merge.psam", here("data/to_qtl/genotypes.psam")))

phenotype_files    = here("data/to_qtl", paste0("phenotypes_", c("metabolomics_plasma", "metabolomics_urine", "proteomics"), ".txt"))
outcomes_file      = here("data/to_qtl/phenotypes_outcomes.txt")
genotype_prefix    = here("data/to_qtl/genotypes"        )
id_conversion_file = here("data/to_qtl/id_conversion.txt")
pca_file           = here("data/to_qtl/genotype_pcs.txt" )
covariate_file     = here("data/to_qtl/covariates.txt"   )
threads            = 100

################################################################################
# Run test
source(here("R/qtl_analysis_wrapper.r"))

# Check basic file formats
check_input_files(
  phenotype_file = phenotype_files[[1]],
  genotype_prefix = genotype_prefix,
  id_conversion_file = id_conversion_file,
  pca_file = pca_file,
  covariate_file = covariate_file
)

# Diagnose ID matching issues
diagnose_id_matching(
  phenotype_file = phenotype_files[[1]],
  genotype_prefix = genotype_prefix,
  id_conversion_file = id_conversion_file,
  pca_file = pca_file,
  covariate_file = covariate_file
)

################################################################################
# Run analysis: TEST
source(here("R/qtl_analysis_wrapper.r"))

results <- run_qtl_analysis(
  phenotype_file       = phenotype_files[[1]],
  genotype_prefix      = genotype_prefix,
  id_conversion_file   = id_conversion_file,
  pca_file             = pca_file,
  covariate_file       = covariate_file,
  output_prefix        = here("output", "test"),
  temp_dir             = here("output", "test", "processing"),
  phenotype_list       = c("_s_-3-Hydroxyisobutyric.Acid.in.uM", "3-Methylcrotonylglycine_Tiglyglycine.in.uM"),  # Only these phenotypes
  n_pcs                = 10,
  covariate_columns    = c("sex", "age", paste0("race", 1:3)),
  #normalization        = "rank_normal", # Phenotypes: rank-normal for proteomics
  normalization        = "log_scale",    # phenotypes log-transform + scaling for metabolomics
  normalize_covariates = TRUE,           # Covariates: smart scaling
  normalize_pcs        = TRUE,           # PCs: z-score scaling
  outlier_removal      = 4,
  maf_threshold        = 0.01,           # 1% MAF threshold
  hwe_threshold        = 1e-6,           # HWE p-value threshold  
  threads              = threads,
  memory               = 1000000
)

################################################################################


################################################################################
# Run analysis: metabolomics plasma
source(here("R/qtl_analysis_wrapper.r"))

results <- run_qtl_analysis(
  phenotype_file       = phenotype_files[[1]],
  genotype_prefix      = genotype_prefix,
  id_conversion_file   = id_conversion_file,
  pca_file             = pca_file,
  covariate_file       = covariate_file,
  output_prefix        = here("output", "metabolomics_plasma.standard_covariates"),
  temp_dir             = here("output", "metabolomics_plasma.standard_covariates", "processing"),
  #phenotype_list       = c("Creatinine.in.uM", "Lysine.in.uM", "Ornithine.in.uM", "Arginine.in.uM"),  # Only these phenotypes
  n_pcs                = 10,
  covariate_columns    = c("sex", "age", paste0("race", 1:3)),
  #normalization        = "rank_normal", # Phenotypes: rank-normal for proteomics
  normalization        = "log_scale",    # phenotypes log-transform + scaling for metabolomics
  normalize_covariates = TRUE,           # Covariates: smart scaling
  normalize_pcs        = TRUE,           # PCs: z-score scaling
  outlier_removal      = 4,
  maf_threshold        = 0.05,           # 1% MAF threshold
  hwe_threshold        = 1e-6,           # HWE p-value threshold  
  threads              = threads,
  memory               = 1000000,
  parallel_phenotypes  = 12,      # Run 16 phenotypes simultaneously
  batch_size           = 36 
)

################################################################################

################################################################################
# Run analysis: metabolomics urine
source(here("R/qtl_analysis_wrapper.r"))

results <- run_qtl_analysis(
  phenotype_file       = phenotype_files[[2]],
  genotype_prefix      = genotype_prefix,
  id_conversion_file   = id_conversion_file,
  pca_file             = pca_file,
  covariate_file       = covariate_file,
  output_prefix        = here("output", "metabolomics_urine.standard_covariates"),
  temp_dir             = here("output", "metabolomics_urine.standard_covariates", "processing"),
  #phenotype_list       = c("Creatinine.in.uM", "Lysine.in.uM", "Ornithine.in.uM", "Arginine.in.uM"),  # Only these phenotypes
  n_pcs                = 10,
  covariate_columns    = c("sex", "age", paste0("race", 1:3)),
  #normalization        = "rank_normal", # Phenotypes: rank-normal for proteomics
  normalization        = "log_scale",    # phenotypes log-transform + scaling for metabolomics
  normalize_covariates = TRUE,           # Covariates: smart scaling
  normalize_pcs        = TRUE,           # PCs: z-score scaling
  outlier_removal      = 4,
  maf_threshold        = 0.05,           # 1% MAF threshold
  hwe_threshold        = 1e-6,           # HWE p-value threshold  
  threads              = threads,
  memory               = 1000000,
  parallel_phenotypes  = 12,      # Run 16 phenotypes simultaneously
  batch_size           = 36  
)

################################################################################

################################################################################
# Run analysis: proteomics
source(here("R/qtl_analysis_wrapper.r"))

results <- run_qtl_analysis(
  phenotype_file       = phenotype_files[[3]],
  genotype_prefix      = genotype_prefix,
  id_conversion_file   = id_conversion_file,
  pca_file             = pca_file,
  covariate_file       = covariate_file,
  output_prefix        = here("output", "proteomics.standard_covariates"),
  temp_dir             = here("output", "proteomics.standard_covariates", "processing"),
  #phenotype_list       = c("Creatinine.in.uM", "Lysine.in.uM", "Ornithine.in.uM", "Arginine.in.uM"),  # Only these phenotypes
  n_pcs                = 10,
  covariate_columns    = c("sex", "age", paste0("race", 1:3)),
  normalization        = "rank_normal",  # Phenotypes: rank-normal for proteomics
  #normalization        = "log_scale",   # phenotypes log-transform + scaling for metabolomics
  normalize_covariates = TRUE,           # Covariates: smart scaling
  normalize_pcs        = TRUE,           # PCs: z-score scaling 
  outlier_removal      = 4,
  maf_threshold        = 0.05,           # 1% MAF threshold
  hwe_threshold        = 1e-6,           # HWE p-value threshold  
  threads              = threads,
  memory               = 1000000,
  parallel_phenotypes  = 25,      # Run 16 phenotypes simultaneously
  batch_size           = 100   
)

################################################################################


################################################################################
# Run analysis: outcomes
source(here("R/qtl_analysis_wrapper.r"))

results <- run_qtl_analysis(
  phenotype_file       = outcomes_file,
  genotype_prefix      = genotype_prefix,
  id_conversion_file   = id_conversion_file,
  pca_file             = pca_file,
  covariate_file       = covariate_file,
  output_prefix        = here("output", "outcomes.standard_covariates"),
  temp_dir             = here("output", "outcomes.standard_covariates", "processing"),
  #phenotype_list       = colnames(fread(outcomes_file, sep = "\t", header = TRUE, data.table = FALSE) %>% select(-sample_id)),
  #n_pcs                = 10,
  #covariate_columns    = c("sex", "age", paste0("race", 1:3)),
  n_pcs                = 3,
  covariate_columns    = c("sex", "age"),
  normalization        = "none",  
  normalize_covariates = TRUE,           # Covariates: smart scaling
  normalize_pcs        = TRUE,           # PCs: z-score scaling 
  outlier_removal      = 0,
  maf_threshold        = 0.05,           # 1% MAF threshold
  hwe_threshold        = 1e-6,           # HWE p-value threshold  
  threads              = threads,
  memory               = 1000000,
  parallel_phenotypes  = 40,      # Run 16 phenotypes simultaneously
  batch_size           = 40   
)

################################################################################

