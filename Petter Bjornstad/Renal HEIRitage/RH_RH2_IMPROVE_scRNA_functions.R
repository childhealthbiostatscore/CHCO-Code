library(Seurat)
library(nebula)
library(Matrix)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(purrr)
library(aws.s3)
library(jsonlite)

# Set up environment for Kopah
keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

pb90_subset_clean <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset.rds", bucket = "scrna", region = "")

# Set up environment for Kopah
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else {
  stop("Unknown user: please specify root path for this user.")
}

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)


# ===========================================================================
# Function: run_nebula_parallel
# ===========================================================================

run_nebula_parallel <- function(seurat_obj,
                                n_cores = 10,
                                layer = "counts",
                                subject_id_col = "record_id",
                                offset_col = "pooled_offset",
                                formula = ~ group,
                                model = "NBLMM",
                                reml = 1,
                                output_re = TRUE,
                                covariance = TRUE,
                                s3_bucket = "scrna",
                                s3_key = NULL,
                                verbose = TRUE, 
                                group = F) {
  
  # Extract counts and gene list
  counts_mat <- round(GetAssayData(seurat_obj, layer = layer))
  genes_list <- rownames(counts_mat)
  
  # Set up parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Ensure cleanup on exit
  on.exit({
    stopCluster(cl)
  }, add = TRUE)
  
  start_time <- Sys.time()
  
  # Run nebula in parallel
  nebula_results_list <- foreach(g = genes_list, 
                                 .packages = c("nebula", "Matrix"),
                                 .errorhandling = "pass") %dopar% {
                                   warn <- NULL
                                   err <- NULL
                                   res <- NULL
                                   
                                   tryCatch({
                                     # Subset data for single gene
                                     count_gene <- counts_mat[g, , drop = FALSE]
                                     meta_gene <- subset(seurat_obj, features = g)@meta.data
                                     
                                     # Create model matrix
                                     pred_gene <- model.matrix(formula, data = meta_gene)
                                     
                                     # Group cells
                                     if (group) {
                                       data_g_gene <- group_cell(count = count_gene, 
                                                                 id = meta_gene[[subject_id_col]], 
                                                                 pred = pred_gene)
                                     } else {
                                       data_g_gene <- list(count = count_gene, 
                                                           id = meta_gene[[subject_id_col]], 
                                                           pred = pred_gene)
                                     }
                                     
                                     # Run nebula with warning handling
                                     res <- withCallingHandlers({
                                       nebula(count = data_g_gene$count, 
                                              id = data_g_gene$id, 
                                              pred = data_g_gene$pred,
                                              ncore = 1, 
                                              output_re = output_re, 
                                              covariance = covariance, 
                                              reml = reml, 
                                              model = model, 
                                              offset = meta_gene[[offset_col]])
                                     }, warning = function(w) {
                                       warn <<- conditionMessage(w)
                                       invokeRestart("muffleWarning")
                                     })
                                   }, error = function(e) {
                                     err <<- conditionMessage(e)
                                   })
                                   
                                   list(gene = g, result = res, warning = warn, error = err)
                                 }
  
  # Report warnings and errors if verbose
  if (verbose) {
    for (res in nebula_results_list) {
      if (!is.null(res$warning)) {
        cat(sprintf("Warning for gene %s: %s\n", res$gene, res$warning))
      }
      if (!is.null(res$error)) {
        cat(sprintf("Error for gene %s: %s\n", res$gene, res$error))
      }
    }
  }
  
  # Clean up results
  names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)
  nebula_results_list <- lapply(nebula_results_list, function(x) x$result)
  nebula_results_list <- Filter(Negate(is.null), nebula_results_list)
  
  # Calculate runtime and non-convergence rate
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  nebula_nonconverged_percent <- (length(genes_list) - length(nebula_results_list)) / length(genes_list)
  
  if (verbose) {
    cat(sprintf("\nRuntime: %.2f minutes\n", runtime))
    cat(sprintf("%.2f%% of genes filtered due to low expression or convergence issues\n", 
                nebula_nonconverged_percent * 100))
  }
  
  # Save to S3 if requested
  if (!is.null(s3_key) && !is.null(s3_bucket)) {
    s3saveRDS(nebula_results_list, bucket = s3_bucket, object = s3_key, region = "")
    
    if (verbose) {
      cat(sprintf("Results uploaded to s3://%s/%s\n", s3_bucket, s3_key))
    }
  }
  
  # Return results and metadata
  return(list(
    results = nebula_results_list,
    n_genes_tested = length(genes_list),
    n_genes_converged = length(nebula_results_list),
    nonconverged_percent = nebula_nonconverged_percent,
    runtime_minutes = as.numeric(runtime)
  ))
}

# ===========================================================================
# Function: process_nebula_results
# ===========================================================================

process_nebula_results <- function(nebula_list, 
                                   pval_col = "p_groupType_1_Diabetes", 
                                   logfc_col = "logFC_groupType_1_Diabetes",
                                   convergence_cut = -10,
                                   logfc_cut = 10) {
  # Extract convergence codes
  convergence_df <- purrr::map_dfr(names(nebula_list), function(gene_name) {
    convergence_code <- nebula_list[[gene_name]]$convergence
    data.frame(Gene = gene_name, Convergence_Code = convergence_code)
  })
  
  # Filter to converged models
  converged_genes <- convergence_df %>%
    filter(Convergence_Code >= convergence_cut) %>%
    pull(Gene)
  
  # Combine model summary results
  summary_df <- purrr::map_dfr(converged_genes, function(gene_name) {
    nebula_list[[gene_name]]$summary %>%
      dplyr::mutate(Gene = gene_name)
  }) %>%
    filter(abs(.data[[logfc_col]]) < logfc_cut)
  
  # Add FDR adjustment
  if (pval_col %in% names(summary_df)) {
    summary_df <- summary_df %>%
      dplyr::mutate(fdr = p.adjust(.data[[pval_col]], method = "fdr")) 
  } else {
    warning(paste("Column", pval_col, "not found in summary data. FDR not computed."))
    summary_df$fdr <- NA
  }
  
  # Extract overdispersion estimates
  overdisp_df <- purrr::map_dfr(names(nebula_list), function(gene_name) {
    od <- nebula_list[[gene_name]]$overdispersion
    od$Gene <- gene_name
    od
  })
  
  return(list(
    convergence     = convergence_df,
    results         = summary_df,
    overdispersion  = overdisp_df
  ))
}
