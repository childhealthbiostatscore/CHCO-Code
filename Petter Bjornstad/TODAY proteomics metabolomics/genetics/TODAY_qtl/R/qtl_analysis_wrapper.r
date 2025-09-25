# QTL Analysis Wrapper
# Author: Research Assistant
# Description: Comprehensive R wrapper for QTL analysis using PLINK2

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(qqman)
library(parallel)

#' Run QTL Analysis using PLINK2
#' 
#' @param phenotype_file Path to phenotype file (TAB-separated)
#' @param genotype_prefix Prefix for PLINK2 genotype files (.pgen, .pvar, .psam)
#' @param id_conversion_file Path to ID conversion file (genotype_id, phenotype_id)
#' @param pca_file Path to genotype PCA file
#' @param covariate_file Path to covariate file (sex, age, etc.)
#' @param output_prefix Prefix for output files
#' @param phenotype_list Vector of phenotype column names to analyze (default: NULL = analyze all)
#' @param n_pcs Number of PCs to include as covariates (default: 10)
#' @param covariate_columns Vector of column names from covariate file to include
#' @param threads Number of threads for PLINK2 (default: 4)
#' @param memory Memory limit in MB for PLINK2 (default: 8000)
#' @param plink2_path Path to PLINK2 executable (default: "plink2")
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @param normalization Method for phenotype normalization (default: "rank_normal")
#'   Options: "none", "log", "log_scale", "scale", "rank_normal"
#' @param normalize_covariates Whether to normalize continuous covariates (default: TRUE)
#' @param normalize_pcs Whether to normalize genotype PCs (default: TRUE)
#' @param maf_threshold Minor allele frequency threshold (default: 0.01 = 1%)
#' @param hwe_threshold Hardy-Weinberg equilibrium p-value threshold (default: 1e-6)
#' @param skip_completed Whether to skip phenotypes with existing summary stats files (default: TRUE)
#' @param parallel_phenotypes Number of phenotypes to run in parallel (default: 1, max: threads/2)
#' @param batch_size Number of phenotypes to process in each batch for memory management (default: 50)
#' @param outlier_removal Remove outliers beyond this many SDs from mean (default: NULL)
#' @param temp_dir Temporary directory for intermediate files (default: "processing")
#'   If "processing", creates a processing subfolder in the output directory
#' 
#' @return List containing summary statistics and file paths to plots
#' 
run_qtl_analysis <- function(
    phenotype_file,
    genotype_prefix,
    id_conversion_file,
    pca_file,
    covariate_file,
    output_prefix,
    phenotype_list = NULL,
    n_pcs = 10,
    covariate_columns = c("sex", "age"),
    threads = 4,
    memory = 8000,
    plink2_path = "plink2",
    fdr_threshold = 0.05,
    normalization = "rank_normal",
    normalize_covariates = TRUE,
    normalize_pcs = TRUE,
    maf_threshold = 0.01,
    hwe_threshold = 1e-6,
    skip_completed = TRUE,
    parallel_phenotypes = 1,
    batch_size = 50,
    outlier_removal = NULL,
    temp_dir = "processing"
) {
    
    # Create output directory structure
    output_dir <- output_prefix
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create subfolders
    stats_dir <- file.path(output_dir, "summary_statistics")
    manhattan_dir <- file.path(output_dir, "manhattan_plots")
    qq_dir <- file.path(output_dir, "qq_plots")
    
    dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(manhattan_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(qq_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Handle temp directory - if "processing", create in output folder
    if(temp_dir == "processing") {
        temp_dir <- file.path(output_dir, "processing")
    }
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
    
    cat("Starting QTL Analysis Pipeline...\n")
    cat("=================================\n")
    cat(sprintf("Output directory: %s\n", output_dir))
    cat(sprintf("Summary statistics: %s\n", stats_dir))
    cat(sprintf("Manhattan plots: %s\n", manhattan_dir))
    cat(sprintf("QQ plots: %s\n", qq_dir))
    cat(sprintf("Processing directory: %s\n", temp_dir))
    
    # Validate and adjust parallelization settings
    max_parallel <- max(1, threads %/% 2)  # Don't use more than half available threads for parallel phenotypes
    if(parallel_phenotypes > max_parallel) {
        warning(sprintf("Reducing parallel_phenotypes from %d to %d (max: threads/2)", parallel_phenotypes, max_parallel))
        parallel_phenotypes <- max_parallel
    }
    
    # Calculate resources per worker for parallel processing
    if(parallel_phenotypes > 1) {
        threads_per_worker <- max(1, threads %/% parallel_phenotypes)
        memory_per_worker <- max(1000, memory %/% parallel_phenotypes)
        cat(sprintf("Parallel processing: %d workers, %d threads and %d MB per worker\n", 
                    parallel_phenotypes, threads_per_worker, memory_per_worker))
    } else {
        threads_per_worker <- threads
        memory_per_worker <- memory
        cat(sprintf("Sequential processing: %d threads, %d MB total\n", threads, memory))
    }
    
    # Save analysis parameters
    params_file <- file.path(output_dir, "analysis_parameters.txt")
    params_info <- data.frame(
        parameter = c(
            "phenotype_file", "genotype_prefix", "id_conversion_file", "pca_file", "covariate_file",
            "output_prefix", "phenotype_list", "n_pcs", "covariate_columns", "threads", "memory",
            "plink2_path", "fdr_threshold", "normalization", "normalize_covariates", "normalize_pcs",
            "maf_threshold", "hwe_threshold", "skip_completed", "parallel_phenotypes", 
            "batch_size", "outlier_removal", "temp_dir", "analysis_date", "analysis_time"
        ),
        value = c(
            phenotype_file, genotype_prefix, id_conversion_file, pca_file, covariate_file,
            output_prefix, 
            if(is.null(phenotype_list)) "NULL (all phenotypes)" else paste(phenotype_list, collapse = ", "),
            n_pcs, paste(covariate_columns, collapse = ", "), threads, memory,
            plink2_path, fdr_threshold, normalization, normalize_covariates, normalize_pcs,
            maf_threshold, hwe_threshold, skip_completed, parallel_phenotypes, 
            batch_size, if(is.null(outlier_removal)) "NULL" else outlier_removal,
            temp_dir, Sys.Date(), format(Sys.time(), "%H:%M:%S")
        ),
        stringsAsFactors = FALSE
    )
    fwrite(params_info, params_file, sep = "\t")
    cat(sprintf("Analysis parameters saved: %s\n", basename(params_file)))
    
    # Step 1: Load and validate input files
    cat("Step 1: Loading input files...\n")
    
    # Load ID conversion file
    id_conv <- fread(id_conversion_file, header = TRUE)
    colnames(id_conv)[1:2] <- c("genotype_id", "phenotype_id")
    cat(sprintf("  - Loaded %d ID mappings\n", nrow(id_conv)))
    
    # Load phenotype file
    phenotypes <- fread(phenotype_file, header = TRUE)
    cat(sprintf("  - Loaded phenotypes: %d samples x %d phenotypes\n", 
                nrow(phenotypes), ncol(phenotypes)-1))
    
    # Load PCA file
    pca_data <- fread(pca_file, header = TRUE)
    cat(sprintf("  - Loaded PCA data: %d samples x %d PCs\n", 
                nrow(pca_data), ncol(pca_data)-1))
    
    # Load covariate file
    covariates <- fread(covariate_file, header = TRUE)
    cat(sprintf("  - Loaded covariates: %d samples x %d variables\n", 
                nrow(covariates), ncol(covariates)-1))
    
    # Step 2: Merge datasets
    cat("Step 2: Merging datasets...\n")
    
    # Debug: Show column names and sample data
    cat(sprintf("  - ID conversion columns: %s\n", paste(colnames(id_conv), collapse = ", ")))
    cat(sprintf("  - Phenotype file first column: %s\n", colnames(phenotypes)[1]))
    cat(sprintf("  - PCA file first column: %s\n", colnames(pca_data)[1]))
    cat(sprintf("  - Covariate file first column: %s\n", colnames(covariates)[1]))
    
    # Merge phenotypes with ID conversion (phenotypes have sample_id = phenotype_id)
    pheno_merged <- merge(phenotypes, id_conv, 
                         by.x = colnames(phenotypes)[1], 
                         by.y = "phenotype_id", all = FALSE)
    
    cat(sprintf("  - After phenotype + ID conversion merge: %d samples\n", nrow(pheno_merged)))
    
    # Select PCs to include
    pc_cols <- paste0("PC", 1:n_pcs)
    available_pcs <- intersect(pc_cols, colnames(pca_data))
    if(length(available_pcs) < n_pcs) {
        warning(sprintf("Only %d PCs available, using all", length(available_pcs)))
        pc_cols <- available_pcs
    }
    
    # Merge with PCA data (PCA has genotype_id)
    pca_subset <- pca_data[, c(colnames(pca_data)[1], pc_cols), with = FALSE]
    pheno_merged <- merge(pheno_merged, pca_subset, 
                         by.x = "genotype_id", 
                         by.y = colnames(pca_data)[1], all = FALSE)
    
    cat(sprintf("  - After PCA merge: %d samples\n", nrow(pheno_merged)))
    
    # Merge with covariates (covariates have sample_id = phenotype_id)
    available_covs <- intersect(covariate_columns, colnames(covariates))
    if(length(available_covs) < length(covariate_columns)) {
        missing_covs <- setdiff(covariate_columns, available_covs)
        warning(sprintf("Missing covariates: %s", paste(missing_covs, collapse = ", ")))
    }
    
    cov_subset <- covariates[, c(colnames(covariates)[1], available_covs), with = FALSE]
    pheno_merged <- merge(pheno_merged, cov_subset, 
                         by.x = colnames(phenotypes)[1], 
                         by.y = colnames(covariates)[1], all = FALSE)
    
    cat(sprintf("  - After covariate merge: %d samples\n", nrow(pheno_merged)))
    cat(sprintf("  - Final merged dataset: %d samples\n", nrow(pheno_merged)))
    
    # Step 2b: Normalize covariates and PCs
    cat("Step 2b: Normalizing covariates and PCs...\n")
    
    # Normalize genotype PCs
    if(normalize_pcs && length(pc_cols) > 0) {
        cat(sprintf("  - Normalizing %d genotype PCs...\n", length(pc_cols)))
        for(pc in pc_cols) {
            original_data <- pheno_merged[[pc]]
            if(sum(!is.na(original_data)) > 5) {
                # Always use scaling (mean=0, sd=1) for PCs
                pheno_merged[, (pc) := scale(original_data)[,1]]
                cat(sprintf("    - Normalized %s\n", pc))
            } else {
                warning(sprintf("Skipping normalization for %s: too few values", pc))
            }
        }
    } else if(!normalize_pcs) {
        cat("  - Genotype PCs: using raw values (normalize_pcs = FALSE)\n")
    }
    
    # Normalize continuous covariates
    if(normalize_covariates && length(available_covs) > 0) {
        cat(sprintf("  - Normalizing continuous covariates...\n"))
        for(cov in available_covs) {
            original_data <- pheno_merged[[cov]]
            
            # Check if covariate appears to be categorical/binary
            unique_vals <- unique(original_data[!is.na(original_data)])
            n_unique <- length(unique_vals)
            
            if(n_unique <= 2) {
                # Binary/categorical - don't normalize
                cat(sprintf("    - %s: detected as categorical/binary (%d unique values), not normalized\n", 
                           cov, n_unique))
            } else if(n_unique < 10 && all(unique_vals == round(unique_vals))) {
                # Appears to be categorical integers
                cat(sprintf("    - %s: detected as categorical integer (%d unique values), not normalized\n", 
                           cov, n_unique))
            } else if(sum(!is.na(original_data)) > 5) {
                # Continuous variable - normalize using scaling
                pheno_merged[, (cov) := scale(original_data)[,1]]
                cat(sprintf("    - %s: normalized as continuous variable\n", cov))
            } else {
                warning(sprintf("Skipping normalization for %s: too few values", cov))
            }
        }
    } else if(!normalize_covariates) {
        cat("  - Covariates: using raw values (normalize_covariates = FALSE)\n")
    }
    
    # Check for covariate issues that cause PLINK2 VIF_INFINITE errors
    cat("Step 2c: Checking covariate correlations...\n")
    covar_cols <- c(pc_cols, available_covs)
    if(length(covar_cols) > 1) {
        # Get covariate data for correlation check
        covar_data <- pheno_merged[, covar_cols, with = FALSE]
        covar_data_complete <- covar_data[complete.cases(covar_data), ]
        
        if(nrow(covar_data_complete) > 0) {
            # Check for constant variables
            constant_vars <- sapply(covar_data_complete, function(x) length(unique(x)) <= 1)
            if(any(constant_vars)) {
                constant_var_names <- names(constant_vars)[constant_vars]
                warning(sprintf("Removing constant covariates: %s", paste(constant_var_names, collapse = ", ")))
                covar_cols <- setdiff(covar_cols, constant_var_names)
                available_covs <- setdiff(available_covs, constant_var_names)
                pc_cols <- setdiff(pc_cols, constant_var_names)
            }
            
            # Check correlations if we still have multiple variables
            if(length(covar_cols) > 1) {
                covar_data_clean <- pheno_merged[, covar_cols, with = FALSE]
                covar_data_clean <- covar_data_clean[complete.cases(covar_data_clean), ]
                
                if(nrow(covar_data_clean) > 0) {
                    cor_matrix <- cor(covar_data_clean, use = "complete.obs")
                    
                    # Check for perfect correlations (excluding diagonal)
                    cor_matrix_abs <- abs(cor_matrix)
                    diag(cor_matrix_abs) <- 0
                    perfect_cor <- which(cor_matrix_abs > 0.99, arr.ind = TRUE)
                    
                    if(nrow(perfect_cor) > 0) {
                        # Get names of perfectly correlated variables
                        cor_pairs <- unique(apply(perfect_cor, 1, function(x) {
                            paste(sort(c(rownames(cor_matrix)[x[1]], colnames(cor_matrix)[x[2]])), collapse = " & ")
                        }))
                        warning(sprintf("High correlations detected between covariates: %s", 
                               paste(cor_pairs, collapse = ", ")))
                        cat(sprintf("    - Consider removing one variable from each highly correlated pair\n"))
                    }
                    
                    cat(sprintf("    - Covariate correlation check completed for %d variables\n", length(covar_cols)))
                } else {
                    warning("No complete covariate data available for correlation check")
                }
            }
        } else {
            warning("No complete covariate data available")
        }
    }

    # Step 4: Prepare PLINK2 input files
    cat("Step 4: Preparing PLINK2 input files...\n")
    
    # Identify phenotype columns (exclude ID and covariate columns)
    # We need to exclude both the genotype_id and the original sample_id column
    exclude_cols <- c("genotype_id", "sample_id", pc_cols, available_covs)
    all_pheno_cols <- setdiff(colnames(pheno_merged), exclude_cols)
    
    cat(sprintf("  - Excluding columns: %s\n", paste(exclude_cols, collapse = ", ")))
    cat(sprintf("  - Available phenotype columns: %s\n", paste(head(all_pheno_cols, 5), collapse = ", ")))
    if(length(all_pheno_cols) > 5) {
        cat(sprintf("    ... and %d more\n", length(all_pheno_cols) - 5))
    }
    
    # Handle phenotype selection
    if(is.null(phenotype_list)) {
        # Analyze all phenotypes
        pheno_cols <- all_pheno_cols
        cat(sprintf("  - No phenotype list provided, analyzing all %d phenotypes\n", length(pheno_cols)))
    } else {
        # Validate requested phenotypes exist
        missing_phenos <- setdiff(phenotype_list, all_pheno_cols)
        available_phenos <- intersect(phenotype_list, all_pheno_cols)
        
        if(length(missing_phenos) > 0) {
            warning(sprintf("The following phenotypes were not found in the data: %s", 
                          paste(missing_phenos, collapse = ", ")))
        }
        
        if(length(available_phenos) == 0) {
            stop("None of the requested phenotypes were found in the data!")
        }
        
        pheno_cols <- available_phenos
        cat(sprintf("  - Analyzing %d out of %d requested phenotypes\n", 
                    length(pheno_cols), length(phenotype_list)))
        
        if(length(missing_phenos) > 0) {
            cat(sprintf("  - Available phenotypes: %s\n", paste(pheno_cols, collapse = ", ")))
        }
    }
    
    cat(sprintf("  - Total phenotypes to analyze: %d\n", length(pheno_cols)))
    
    # IMMEDIATELY clean all phenotype names upfront to avoid any file system issues
    cat("Step 3c: Cleaning phenotype names for file system compatibility...\n")
    clean_pheno_names <- gsub("[/\\\\:*?\"<>|()]", "_", pheno_cols)  # Replace problematic characters
    clean_pheno_names <- gsub("_{2,}", "_", clean_pheno_names)       # Replace multiple underscores with single
    clean_pheno_names <- gsub("^_|_$", "", clean_pheno_names)        # Remove leading/trailing underscores
    
    # Check for duplicate clean names and handle them
    duplicate_clean_names <- duplicated(clean_pheno_names) | duplicated(clean_pheno_names, fromLast = TRUE)
    if(any(duplicate_clean_names)) {
        # Add suffixes to make names unique
        clean_name_counts <- table(clean_pheno_names)
        for(clean_name in names(clean_name_counts)[clean_name_counts > 1]) {
            indices <- which(clean_pheno_names == clean_name)
            for(j in seq_along(indices)) {
                if(j > 1) {  # Keep first one as-is, add suffix to others
                    clean_pheno_names[indices[j]] <- paste0(clean_name, "_", j)
                }
            }
        }
        cat(sprintf("  - Fixed %d duplicate clean names\n", sum(duplicate_clean_names)))
    }
    
    # Create mapping between original and clean names
    names(clean_pheno_names) <- pheno_cols
    
    # Show cleaning results
    for(i in seq_along(pheno_cols)) {
        orig_name <- pheno_cols[i]
        clean_name <- clean_pheno_names[i]
        if(orig_name != clean_name) {
            cat(sprintf("  - '%s' -> '%s'\n", orig_name, clean_name))
        }
    }
    cat(sprintf("  - %d phenotype names cleaned\n", sum(pheno_cols != clean_pheno_names)))
    
    # Step 3b: Normalize phenotype data
    cat("Step 3b: Normalizing phenotype data...\n")
    cat(sprintf("  - Normalization method: %s\n", normalization))
    
    if(normalization != "none") {
        for(pheno_name in pheno_cols) {
            original_data <- pheno_merged[[pheno_name]]
            
            # Remove missing values for normalization
            non_missing <- !is.na(original_data)
            if(sum(non_missing) < 10) {
                warning(sprintf("Skipping normalization for %s: too few non-missing values", pheno_name))
                next
            }
            
            # Remove outliers if specified
            if(!is.null(outlier_removal) && outlier_removal > 0) {
                mean_val <- mean(original_data, na.rm = TRUE)
                sd_val <- sd(original_data, na.rm = TRUE)
                outlier_threshold_low <- mean_val - outlier_removal * sd_val
                outlier_threshold_high <- mean_val + outlier_removal * sd_val
                
                outliers <- original_data < outlier_threshold_low | original_data > outlier_threshold_high
                outliers[is.na(outliers)] <- FALSE
                
                if(sum(outliers) > 0) {
                    cat(sprintf("    - Removed %d outliers from %s\n", sum(outliers), pheno_name))
                    pheno_merged[outliers, (pheno_name) := NA]
                    original_data[outliers] <- NA
                    non_missing <- !is.na(original_data)
                }
            }
            
            # Apply normalization
            normalized_data <- original_data
            
            if(normalization == "log") {
                # Natural log transformation (add small constant if needed)
                min_val <- min(original_data, na.rm = TRUE)
                if(min_val <= 0) {
                    shift_val <- abs(min_val) + 1e-6
                    normalized_data[non_missing] <- log(original_data[non_missing] + shift_val)
                    cat(sprintf("    - Applied log transformation with shift %.6f to %s\n", shift_val, pheno_name))
                } else {
                    normalized_data[non_missing] <- log(original_data[non_missing])
                    cat(sprintf("    - Applied log transformation to %s\n", pheno_name))
                }
                
            } else if(normalization == "log_scale") {
                # Log transform then scale (mean=0, sd=1)
                min_val <- min(original_data, na.rm = TRUE)
                if(min_val <= 0) {
                    shift_val <- abs(min_val) + 1e-6
                    log_data <- log(original_data[non_missing] + shift_val)
                } else {
                    log_data <- log(original_data[non_missing])
                }
                normalized_data[non_missing] <- scale(log_data)[,1]
                cat(sprintf("    - Applied log + scaling transformation to %s\n", pheno_name))
                
            } else if(normalization == "scale") {
                # Z-score normalization (mean=0, sd=1)
                normalized_data[non_missing] <- scale(original_data[non_missing])[,1]
                cat(sprintf("    - Applied scaling transformation to %s\n", pheno_name))
                
            } else if(normalization == "rank_normal") {
                # Rank-based inverse normal transformation
                ranks <- rank(original_data[non_missing], ties.method = "average")
                # Convert to quantiles (0,1)
                quantiles <- (ranks - 0.5) / length(ranks)
                # Apply inverse normal transformation
                normalized_data[non_missing] <- qnorm(quantiles)
                cat(sprintf("    - Applied rank-based inverse normal transformation to %s\n", pheno_name))
            }
            
            # Update the merged data
            pheno_merged[, (pheno_name) := normalized_data]
        }
        
        cat(sprintf("  - Normalization completed for %d phenotypes\n", length(pheno_cols)))
    } else {
        cat("  - No normalization applied (normalization = 'none')\n")
    }

    # Initialize results storage
    all_results <- list()
    plot_files <- list()
    summary_stats <- data.frame(
        phenotype = character(),
        n_variants = integer(),
        n_significant = integer(),
        min_p = numeric(),
        top_variant = character(),
        top_variant_p = numeric(),
        top_variant_effect = numeric(),
        stringsAsFactors = FALSE
    )

    # Step 5: Run analysis for each phenotype (with optional parallelization)
    cat("Step 5: Running PLINK2 analysis...\n")
    
    if(parallel_phenotypes > 1) {
        cat(sprintf("Setting up parallel processing with %d workers...\n", parallel_phenotypes))
        
        # Process phenotypes in batches to manage memory
        n_batches <- ceiling(length(pheno_cols) / batch_size)
        cat(sprintf("Processing %d phenotypes in %d batches of up to %d phenotypes each\n", 
                    length(pheno_cols), n_batches, batch_size))
        cat(sprintf("Each worker will use %d threads and %d MB memory\n", threads_per_worker, memory_per_worker))
        
        for(batch_i in 1:n_batches) {
            start_idx <- (batch_i - 1) * batch_size + 1
            end_idx <- min(batch_i * batch_size, length(pheno_cols))
            batch_phenos <- pheno_cols[start_idx:end_idx]
            
            cat(sprintf("\nProcessing batch %d/%d: phenotypes %d-%d\n", 
                        batch_i, n_batches, start_idx, end_idx))
            
            # Split batch into chunks for parallel processing
            batch_chunks <- split(batch_phenos, ceiling(seq_along(batch_phenos) / ceiling(length(batch_phenos) / parallel_phenotypes)))
            
            # Process chunks in parallel using mclapply (simpler than foreach)
            batch_results <- mclapply(batch_chunks, function(chunk_phenos) {
                chunk_results <- list()
                for(pheno_name in chunk_phenos) {
                    result <- process_single_phenotype(
                        pheno_name = pheno_name, 
                        pheno_merged = pheno_merged, 
                        clean_pheno_names = clean_pheno_names, 
                        stats_dir = stats_dir, 
                        temp_dir = temp_dir, 
                        genotype_prefix = genotype_prefix, 
                        pc_cols = pc_cols, 
                        available_covs = available_covs, 
                        covar_cols = covar_cols, 
                        maf_threshold = maf_threshold, 
                        hwe_threshold = hwe_threshold, 
                        threads = threads_per_worker, 
                        memory = memory_per_worker, 
                        plink2_path = plink2_path, 
                        fdr_threshold = fdr_threshold, 
                        skip_completed = skip_completed
                    )
                    if(!is.null(result)) {
                        chunk_results[[pheno_name]] <- result
                    }
                }
                return(chunk_results)
            }, mc.cores = parallel_phenotypes, mc.silent = FALSE)
            
            # Collect results from this batch
            batch_count <- 0
            for(chunk_result in batch_results) {
                if(!is.null(chunk_result)) {
                    for(result in chunk_result) {
                        if(!is.null(result)) {
                            all_results[[result$phenotype]] <- result$data
                            summary_stats <- rbind(summary_stats, result$summary)
                            batch_count <- batch_count + 1
                        }
                    }
                }
            }
            
            cat(sprintf("Batch %d completed: %d results collected\n", batch_i, batch_count))
        }
        
    } else {
        # Sequential processing (original code)
        for(i in seq_along(pheno_cols)) {
            pheno_name <- pheno_cols[i]
            
            result <- process_single_phenotype(
                pheno_name = pheno_name, 
                pheno_merged = pheno_merged, 
                clean_pheno_names = clean_pheno_names, 
                stats_dir = stats_dir, 
                temp_dir = temp_dir, 
                genotype_prefix = genotype_prefix, 
                pc_cols = pc_cols, 
                available_covs = available_covs, 
                covar_cols = covar_cols, 
                maf_threshold = maf_threshold, 
                hwe_threshold = hwe_threshold, 
                threads = threads, 
                memory = memory, 
                plink2_path = plink2_path, 
                fdr_threshold = fdr_threshold, 
                skip_completed = skip_completed
            )
            
            if(!is.null(result)) {
                all_results[[result$phenotype]] <- result$data
                summary_stats <- rbind(summary_stats, result$summary)
            }
        }
    }

    # Step 6: Generate plots and final output
    cat("Step 6: Generating output files...\n")
    
    # Write overall summary statistics (main output folder)
    summary_file <- file.path(output_dir, "analysis_summary.txt")
    fwrite(summary_stats, summary_file, sep = "\t")
    
    # Note: Not creating combined_results.txt to avoid large files with many phenotypes

    # Generate plots for each phenotype
    for(pheno_name in names(all_results)) {
        results <- all_results[[pheno_name]]
        
        if(nrow(results) == 0) next
        
        # Get the correct individual stats file path for this phenotype (RDS format)
        # Find the clean name for this phenotype
        clean_pheno_name <- clean_pheno_names[pheno_name]
        current_individual_stats_file <- file.path(stats_dir, paste0(clean_pheno_name, "_summary_stats.rds"))
        
        # Manhattan plot - derive filename from the correct stats file (RDS)
        manhattan_file <- file.path(manhattan_dir, gsub("_summary_stats\\.rds$", "_manhattan.png", basename(current_individual_stats_file)))
        
        cat(sprintf("    - Creating Manhattan plot for %s: %s\n", pheno_name, basename(manhattan_file)))
        
        tryCatch({
            png(manhattan_file, width = 1200, height = 600, res = 150)
            
            # Prepare data for manhattan plot with optimization
            plot_data <- results[!is.na(results$P), ]
            
            if(nrow(plot_data) > 0) {
                # Speed optimization: Remove variants with very high p-values (close to 1)
                # Keep all significant variants + sample of non-significant variants
                significant_variants <- plot_data[plot_data$FDR < fdr_threshold, ]
                non_significant <- plot_data[plot_data$FDR >= fdr_threshold, ]
                
                # For non-significant variants, keep only those with p < 0.5 to speed up plotting
                if(nrow(non_significant) > 0) {
                    non_significant_filtered <- non_significant[non_significant$P < 0.5, ]
                    
                    # Further thin if still too many points
                    if(nrow(non_significant_filtered) > 500000) {
                        thin_factor <- ceiling(nrow(non_significant_filtered) / 500000)
                        keep_indices <- seq(1, nrow(non_significant_filtered), by = thin_factor)
                        non_significant_filtered <- non_significant_filtered[keep_indices, ]
                        cat(sprintf("    - Thinned non-significant variants for plotting: keeping 1 in %d variants\n", thin_factor))
                    }
                    
                    plot_data <- rbind(significant_variants, non_significant_filtered)
                } else {
                    plot_data <- significant_variants
                }
                
                cat(sprintf("    - Manhattan plot data: %d total variants (%d significant, %d non-significant)\n", 
                            nrow(plot_data), nrow(significant_variants), nrow(plot_data) - nrow(significant_variants)))
                
                # Check required columns exist
                required_cols <- c("CHROM", "POS", "P", "ID")
                missing_cols <- setdiff(required_cols, colnames(plot_data))
                
                if(length(missing_cols) > 0) {
                    cat(sprintf("    - Warning: Missing columns for Manhattan plot: %s\n", 
                               paste(missing_cols, collapse = ", ")))
                    
                    # Simple plot without fancy formatting
                    plot(1:nrow(plot_data), -log10(plot_data$P), 
                         main = paste("Manhattan Plot -", pheno_name),  # Use original name in title
                         xlab = "Variant Index", ylab = "-log10(P-value)",
                         pch = 20, cex = 0.5)
                } else {
                    # Standard manhattan plot - clean and simple  
                    manhattan(plot_data, chr = "CHROM", bp = "POS", p = "P", snp = "ID",
                             main = paste("Manhattan Plot -", pheno_name))  # Use original name in title
                }
            } else {
                # Create empty plot if no data
                plot(1, 1, type = "n", main = paste("Manhattan Plot -", pheno_name),
                     xlab = "Position", ylab = "-log10(P-value)")
                text(1, 1, "No data available", cex = 1.5)
            }
            
            dev.off()
            
            # Verify file was created
            if(file.exists(manhattan_file)) {
                cat(sprintf("    - Manhattan plot created: %s\n", basename(manhattan_file)))
            } else {
                warning(sprintf("Manhattan plot file not created: %s", manhattan_file))
                manhattan_file <- NA  # Mark as failed
            }
            
        }, error = function(e) {
            dev.off()  # Make sure device is closed
            warning(sprintf("Manhattan plot failed for %s: %s", pheno_name, e$message))
            manhattan_file <<- NA  # Mark as failed
        })
        
        # QQ plot - derive filename from the correct stats file (RDS)
        qq_file <- file.path(qq_dir, gsub("_summary_stats\\.rds$", "_qq.png", basename(current_individual_stats_file)))
        
        cat(sprintf("    - Creating QQ plot for %s: %s\n", pheno_name, basename(qq_file)))
        
        tryCatch({
            png(qq_file, width = 800, height = 600, res = 150)
            
            if(nrow(plot_data) > 0) {
                qq(plot_data$P, main = paste("QQ Plot -", pheno_name))  # Use original name in title
            } else {
                # Create empty plot if no data
                plot(1, 1, type = "n", main = paste("QQ Plot -", pheno_name),
                     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
                text(1, 1, "No data available", cex = 1.5)
            }
            
            dev.off()
            
            # Verify file was created
            if(file.exists(qq_file)) {
                cat(sprintf("    - QQ plot created: %s\n", basename(qq_file)))
            } else {
                warning(sprintf("QQ plot file not created: %s", qq_file))
                qq_file <- NA  # Mark as failed
            }
            
        }, error = function(e) {
            dev.off()  # Make sure device is closed
            warning(sprintf("QQ plot failed for %s: %s", pheno_name, e$message))
            qq_file <<- NA  # Mark as failed
        })
        
        plot_files[[pheno_name]] <- list(
            manhattan = manhattan_file, 
            qq = qq_file,
            stats = current_individual_stats_file
        )
    }

    # Clean up processing directory only if analysis completed successfully
    success <- length(all_results) > 0
    if(success) {
        cat(sprintf("Analysis completed successfully! Cleaning up processing directory: %s\n", temp_dir))
        #unlink(temp_dir, recursive = TRUE)
    } else {
        warning(sprintf("Analysis failed or no results generated. Processing files retained in: %s", temp_dir))
    }

    cat("Analysis pipeline finished!\n")
    cat("=================================\n")
    if(success) {
        cat(sprintf("Overall summary: %s\n", basename(summary_file)))
        cat(sprintf("Analysis parameters: %s\n", basename(params_file)))
        cat(sprintf("Individual statistics: %d files in %s/\n", length(all_results), basename(stats_dir)))
        
        # Count actual plot files that exist
        manhattan_count <- sum(sapply(plot_files, function(x) !is.na(x$manhattan) && file.exists(x$manhattan)))
        qq_count <- sum(sapply(plot_files, function(x) !is.na(x$qq) && file.exists(x$qq)))
        
        cat(sprintf("Manhattan plots: %d files in %s/\n", manhattan_count, basename(manhattan_dir)))
        cat(sprintf("QQ plots: %d files in %s/\n", qq_count, basename(qq_dir)))
    } else {
        cat("No results generated - check processing files for debugging\n")
    }

    # Return summary
    return(list(
        summary_stats = summary_stats,
        summary_file = summary_file,
        params_file = params_file,
        plot_files = plot_files,
        n_phenotypes_analyzed = length(all_results),
        processing_dir = if(!success) temp_dir else NULL
    ))
}

#' Process a single phenotype (helper function for parallel processing)
#' 
process_single_phenotype <- function(pheno_name, pheno_merged, clean_pheno_names, stats_dir, temp_dir, 
                                    genotype_prefix, pc_cols, available_covs, covar_cols, maf_threshold, 
                                    hwe_threshold, threads, memory, plink2_path, fdr_threshold, skip_completed) {
    
    clean_pheno_name <- clean_pheno_names[pheno_name]
    
    # Check if analysis already completed
    if(skip_completed) {
        existing_stats_file <- file.path(stats_dir, paste0(clean_pheno_name, "_summary_stats.rds"))
        
        if(file.exists(existing_stats_file)) {
            # Load existing results
            existing_results <- readRDS(existing_stats_file)
            existing_results$phenotype <- pheno_name
            
            # Create summary statistics
            n_significant <- sum(existing_results$FDR < fdr_threshold, na.rm = TRUE)
            min_p_idx <- which.min(existing_results$P)
            min_p_value <- min(existing_results$P, na.rm = TRUE)
            
            if(length(min_p_idx) > 0 && !is.na(min_p_value)) {
                top_variant <- if("ID" %in% colnames(existing_results)) existing_results$ID[min_p_idx] else paste0("chr", existing_results$CHROM[min_p_idx], ":", existing_results$POS[min_p_idx])
                top_variant_p <- min_p_value
                top_variant_effect <- if("BETA" %in% colnames(existing_results)) existing_results$BETA[min_p_idx] else if("OR" %in% colnames(existing_results)) existing_results$OR[min_p_idx] else NA
            } else {
                top_variant <- NA
                top_variant_p <- NA
                top_variant_effect <- NA
            }
            
            summary_row <- data.frame(
                phenotype = pheno_name,
                n_variants = nrow(existing_results),
                n_significant = n_significant,
                min_p = min_p_value,
                top_variant = top_variant,
                top_variant_p = top_variant_p,
                top_variant_effect = top_variant_effect
            )
            
            return(list(phenotype = pheno_name, data = existing_results, summary = summary_row))
        }
    }
    
    # Determine .psam format (cache this for efficiency)
    psam_check <- fread(paste0(genotype_prefix, ".psam"), nrows = 1)
    if(grepl("#?IID", colnames(psam_check)[1], ignore.case = TRUE)) {
        psam_format <- "IID_only"
    } else {
        psam_format <- "FID_IID"
    }
    
    # Prepare phenotype file - format to match .psam
    if(psam_format == "IID_only") {
        current_pheno <- pheno_merged[, c("genotype_id", pheno_name, pc_cols, available_covs), with = FALSE]
        colnames(current_pheno)[1:2] <- c("IID", "PHENO")
    } else {
        current_pheno <- pheno_merged[, c("genotype_id", "genotype_id", pheno_name, pc_cols, available_covs), with = FALSE]
        colnames(current_pheno)[1:3] <- c("FID", "IID", "PHENO")
    }
    
    # Remove missing values
    current_pheno <- current_pheno[complete.cases(current_pheno), ]
    
    if(nrow(current_pheno) < 10) {
        warning(sprintf("Skipping %s: insufficient samples (%d)", pheno_name, nrow(current_pheno)))
        return(NULL)
    }
    
    # Create unique temporary files for this worker
    worker_id <- Sys.getpid()  # Use process ID to ensure uniqueness
    pheno_file <- file.path(temp_dir, paste0("pheno_", clean_pheno_name, "_", worker_id, ".txt"))
    keep_file <- file.path(temp_dir, paste0("keep_samples_", clean_pheno_name, "_", worker_id, ".txt"))
    output_file <- file.path(temp_dir, paste0("results_", clean_pheno_name, "_", worker_id))
    
    # Write phenotype file
    fwrite(current_pheno, pheno_file, sep = "\t")
    
    # Create keep file
    if(psam_format == "IID_only") {
        keep_data <- data.frame(IID = pheno_merged$genotype_id)
    } else {
        keep_data <- data.frame(FID = pheno_merged$genotype_id, IID = pheno_merged$genotype_id)
    }
    fwrite(keep_data, keep_file, sep = "\t", col.names = FALSE)
    
    # Prepare covariate string
    covar_string <- paste(covar_cols, collapse = ",")
    
    # Run PLINK2
    plink2_cmd <- sprintf(
        "%s --pfile %s --keep %s --pheno %s --pheno-name PHENO --covar %s --covar-name %s --maf %f --hwe %e --glm hide-covar --no-input-missing-phenotype --threads %d --memory %d --out %s",
        plink2_path, genotype_prefix, keep_file, pheno_file, pheno_file, 
        covar_string, maf_threshold, hwe_threshold, threads, memory, output_file
    )
    
    # Execute PLINK2 (suppress output for cleaner parallel logs)
    system_result <- system(plink2_cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if(system_result != 0) {
        # Clean up files
        unlink(c(pheno_file, keep_file, paste0(output_file, "*")))
        warning(sprintf("PLINK2 failed for phenotype %s", pheno_name))
        return(NULL)
    }
    
    # Read results
    result_file <- paste0(output_file, ".PHENO.glm.linear")
    if(!file.exists(result_file)) {
        unlink(c(pheno_file, keep_file, paste0(output_file, "*")))
        warning(sprintf("Output file not found for %s", pheno_name))
        return(NULL)
    }
    
    results <- fread(result_file)
    
    # Standardize column names
    if(!"CHROM" %in% colnames(results) && "#CHROM" %in% colnames(results)) {
        setnames(results, "#CHROM", "CHROM")
    }
    if(!"POS" %in% colnames(results) && "BP" %in% colnames(results)) {
        setnames(results, "BP", "POS")
    }
    if(!"ID" %in% colnames(results) && "SNP" %in% colnames(results)) {
        setnames(results, "SNP", "ID")
    }
    
    # Apply FDR correction
    results$FDR <- p.adjust(results$P, method = "BH")
    results$phenotype <- pheno_name
    
    # Optimize results for storage
    essential_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "BETA", "SE", "P", "FDR", "phenotype")
    available_essential_cols <- intersect(essential_cols, colnames(results))
    results_optimized <- results[, available_essential_cols, with = FALSE]
    
    # Reduce precision
    numeric_cols <- c("BETA", "SE", "P", "FDR")
    available_numeric_cols <- intersect(numeric_cols, colnames(results_optimized))
    for(col in available_numeric_cols) {
        results_optimized[[col]] <- signif(results_optimized[[col]], 4)
    }
    
    # Create summary statistics
    n_significant <- sum(results_optimized$FDR < fdr_threshold, na.rm = TRUE)
    min_p_idx <- which.min(results_optimized$P)
    min_p_value <- min(results_optimized$P, na.rm = TRUE)
    
    if(length(min_p_idx) > 0 && !is.na(min_p_value)) {
        top_variant <- if("ID" %in% colnames(results_optimized)) results_optimized$ID[min_p_idx] else paste0("chr", results_optimized$CHROM[min_p_idx], ":", results_optimized$POS[min_p_idx])
        top_variant_p <- min_p_value
        top_variant_effect <- if("BETA" %in% colnames(results_optimized)) results_optimized$BETA[min_p_idx] else if("OR" %in% colnames(results_optimized)) results_optimized$OR[min_p_idx] else NA
    } else {
        top_variant <- NA
        top_variant_p <- NA
        top_variant_effect <- NA
    }
    
    summary_row <- data.frame(
        phenotype = pheno_name,
        n_variants = nrow(results_optimized),
        n_significant = n_significant,
        min_p = min_p_value,
        top_variant = top_variant,
        top_variant_p = top_variant_p,
        top_variant_effect = top_variant_effect
    )
    
    # Save results
    individual_stats_file <- file.path(stats_dir, paste0(clean_pheno_name, "_summary_stats.rds"))
    saveRDS(results_optimized, individual_stats_file, compress = TRUE)
    
    # Clean up temporary files
    unlink(c(pheno_file, keep_file, paste0(output_file, "*")))
    
    return(list(phenotype = pheno_name, data = results_optimized, summary = summary_row))
}

#' Example usage function
#' 
example_usage <- function() {
    # Example function call
    results <- run_qtl_analysis(
        phenotype_file = "path/to/phenotypes.txt",
        genotype_prefix = "path/to/genotypes",
        id_conversion_file = "path/to/id_conversion.txt",
        pca_file = "path/to/pca.txt",
        covariate_file = "path/to/covariates.txt",
        output_prefix = "output/qtl_analysis_folder",  # Now a directory
        phenotype_list = c("protein1", "protein2", "metabolite1"),  # Specific phenotypes
        n_pcs = 10,
        covariate_columns = c("sex", "age", "batch"),
        threads = 8,
        memory = 16000,
        fdr_threshold = 0.05,
        normalization = "rank_normal",        # Recommended for QTL analysis
        normalize_covariates = TRUE,          # Normalize continuous covariates
        normalize_pcs = TRUE,                 # Normalize genotype PCs
        maf_threshold = 0.01,                 # MAF > 1%
        hwe_threshold = 1e-6,                 # HWE p > 1e-6
        skip_completed = TRUE,                # Skip already completed analyses
        parallel_phenotypes = 4,              # Run 4 phenotypes simultaneously
        batch_size = 50,                      # Process 50 phenotypes per batch
        outlier_removal = 4                   # Remove outliers beyond 4 SDs
        # temp_dir = "processing"             # Default: creates processing/ in output folder
    )
    
    # Print summary
    print(results$summary_stats)
}

#' Utility function to check file formats and requirements
#' 
check_input_files <- function(phenotype_file, genotype_prefix, id_conversion_file, 
                             pca_file, covariate_file) {
    
    cat("Checking input files...\n")
    
    # Check if files exist
    files_to_check <- c(phenotype_file, id_conversion_file, pca_file, covariate_file,
                        paste0(genotype_prefix, c(".pgen", ".pvar", ".psam")))
    
    for(file in files_to_check) {
        if(!file.exists(file)) {
            stop(sprintf("File not found: %s", file))
        } else {
            cat(sprintf("✓ Found: %s\n", file))
        }
    }
    
    # Check file formats
    cat("\nChecking file formats...\n")
    
    # Check phenotype file
    pheno <- fread(phenotype_file, nrows = 5)
    cat(sprintf("✓ Phenotype file: %d columns, first few: %s\n", 
                ncol(pheno), paste(colnames(pheno)[1:min(5, ncol(pheno))], collapse = ", ")))
    
    # Check ID conversion file
    id_conv <- fread(id_conversion_file, nrows = 5)
    cat(sprintf("✓ ID conversion file: %d columns, headers: %s\n", 
                ncol(id_conv), paste(colnames(id_conv), collapse = ", ")))
    
    # Check PCA file
    pca <- fread(pca_file, nrows = 5)
    cat(sprintf("✓ PCA file: %d columns, first few: %s\n", 
                ncol(pca), paste(colnames(pca)[1:min(5, ncol(pca))], collapse = ", ")))
    
    # Check covariate file
    cov <- fread(covariate_file, nrows = 5)
    cat(sprintf("✓ Covariate file: %d columns, headers: %s\n", 
                ncol(cov), paste(colnames(cov), collapse = ", ")))
    
    # Check PLINK2 .psam file (sample information)
    psam_file <- paste0(genotype_prefix, ".psam")
    if(file.exists(psam_file)) {
        psam <- fread(psam_file, nrows = 10)
        cat(sprintf("✓ PLINK2 .psam file: %d columns, headers: %s\n", 
                    ncol(psam), paste(colnames(psam), collapse = ", ")))
        
        # PLINK2 .psam files should have FID, IID as first two columns
        if(ncol(psam) >= 2) {
            cat(sprintf("  Sample IDs in genotype file (FID): %s\n", 
                        paste(head(psam[[1]], 3), collapse = ", ")))
            cat(sprintf("  Sample IDs in genotype file (IID): %s\n", 
                        paste(head(psam[[2]], 3), collapse = ", ")))
        }
    }
    
    cat("\nAll files passed basic checks!\n")
}