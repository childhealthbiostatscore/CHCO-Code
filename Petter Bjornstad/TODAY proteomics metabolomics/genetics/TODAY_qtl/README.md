# prepare_input_data.R
- Get all the required data (ID conversions, phenotypes and covariates) and return text files ready for running QTLs

# run_qtls.R
- Run the QTL analysis for each phenotype

# qtl_analysis_wrapper.r
- V01 (2025/08/29): works, not well optimized for 7k phenotypes
- V02: parallelized, works in batches