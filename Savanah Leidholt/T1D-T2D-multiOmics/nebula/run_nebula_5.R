#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
cell_name <- args[1]
contrast_name <- args[2]

if (is.na(cell_name) || is.na(contrast_name)) {
  stop("Usage: Rscript run_nebula_5.R <cell_name> <contrast_name>")
}

.libPaths("/mmfs1/gscratch/togo/leidholt/R_SLL_Seurat/")

suppressPackageStartupMessages({
  library(Matrix)
  library(nebula)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(aws.s3)
  library(jsonlite)
})

ncore <-  min(8, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
chunk_size <- 2000

keys <- jsonlite::fromJSON("/mmfs1/home/leidholt/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

#importing cleaned celltypes subsets
in_object <- paste0(
  "data/nebula_prepped_inputs/",
  "nebula_input_",
  cell_name,
  ".rds"
)

cat("Reading prepped input from S3:", in_object, "\n")

prep <- s3readRDS(
  object = in_object,
  bucket = "triad",
  region = "",
  base_url = "s3.kopah.uw.edu",
  use_https = TRUE
)

counts_nebula <- prep$count
meta_nebula <- prep$meta

#running nebula pipeline
run_nebula_chunk <- function(count_chunk, chunk_id) {
  
  warn_msg <- NULL
  err_msg <- NULL
  
  fit <- tryCatch(
    withCallingHandlers(
      nebula(
        count = round(count_chunk),
        id = meta_sub$record_id,
        pred = pred,
        offset = offset,
        model = "NBLMM",
        reml = 1,
        ncore = ncore
      ),
      warning = function(w) {
        warn_msg <<- paste(c(warn_msg, conditionMessage(w)), collapse = " | ")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err_msg <<- conditionMessage(e)
      NULL
    }
  )
  
  list(
    chunk_id = chunk_id,
    genes = rownames(count_chunk),
    fit = fit,
    warning = warn_msg,
    error = err_msg
  )
}

gene_chunks <- split(
  rownames(counts_sub),
  ceiling(seq_along(rownames(counts_sub)) / chunk_size)
)

chunk_results <- map2(
  gene_chunks,
  seq_along(gene_chunks),
  function(gene_set, chunk_id) {
    cat("Running chunk", chunk_id, "of", length(gene_chunks), "\n")
    run_nebula_chunk(
      count_chunk = counts_sub[gene_set, , drop = FALSE],
      chunk_id = chunk_id
    )
  }
)

diagnostics <- map_dfr(
  chunk_results,
  function(x) {
    tibble(
      chunk_id = x$chunk_id,
      gene = x$genes,
      chunk_success = !is.null(x$fit),
      warning = x$warning,
      error = x$error
    )
  }
)

successful_chunks <- keep(chunk_results, ~ !is.null(.x$fit))

if (length(successful_chunks) == 0) {
  nebula_out <- list(
    results = NULL,
    diagnostics = diagnostics,
    donor_table = donor_tab
  )
} else {
  
  tt <- map_dfr(
    successful_chunks,
    function(x) {
      
      df <- as.data.frame(x$fit$summary)
      
      if (!"gene" %in% colnames(df)) {
        df <- rownames_to_column(df, "gene")
      }
      
      if (!is.null(x$fit$convergence)) {
        conv_df <- tibble(
          gene = rownames(x$fit$summary),
          convergence = x$fit$convergence
        )
        
        df <- left_join(df, conv_df, by = "gene")
      } else {
        df$convergence <- NA_real_
      }
      
      mutate(df, chunk_id = x$chunk_id)
    }
  )
  
  tt <- tt %>%
    mutate(
      stable_fit = is.na(convergence) | convergence >= -10,
      comparison = contrast_name,
      celltype = cell_name,
      contrast_level_0 = pair[1],
      contrast_level_1 = pair[2],
      logFC_direction = paste0(pair[2], "_vs_", pair[1]),
      n_cells = nrow(meta_sub),
      n_genes_tested = nrow(counts_sub),
      n_genes_successful = n_distinct(gene),
      n_donors_group1 = donor_tab$n[match(pair[1], donor_tab$group)],
      n_donors_group2 = donor_tab$n[match(pair[2], donor_tab$group)],
      p_adj = p.adjust(p_group_contrast, method = "fdr"),
      neg_log10_p = -log10(pmax(p_group_contrast, 1e-300)),
      neg_log10_fdr = -log10(pmax(p_adj, 1e-300))
    )
  
  nebula_out <- list(
    results = tt,
    diagnostics = diagnostics,
    donor_table = donor_tab,
    gene_filter_summary = list(
      n_genes_after_light_filter = nrow(counts_sub)
    )
  )
}

tmp_rds <- tempfile(fileext = ".rds")
saveRDS(nebula_out, tmp_rds)

out_object <- paste0(
  "results/nebula/nebula_5/",
  "nebula_",
  cell_name,
  "_",
  contrast_name,
  ".rds"
)

put_object(
  file = tmp_rds,
  object = out_object,
  bucket = "triad",
  base_url = "s3.kopah.uw.edu",
  region = "",
  use_https = TRUE
)

unlink(tmp_rds)

cat("Saved to:", out_object, "\n")
cat("Done\n")
