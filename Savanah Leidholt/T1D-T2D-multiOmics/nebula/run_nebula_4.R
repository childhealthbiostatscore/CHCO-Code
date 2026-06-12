#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

cell_name <- args[1]
contrast_name <- args[2]

if (is.na(cell_name) || is.na(contrast_name)) {
  stop("Usage: Rscript run_nebula_4.R <cell_name> <contrast_name>")
}

.libPaths("/mmfs1/gscratch/togo/leidholt/R_SLL_Seurat/")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(nebula)
  library(purrr)
  library(aws.s3)
  library(jsonlite)
  library(tibble)
})

cat("Cell type:", cell_name, "\n")
cat("Contrast:", contrast_name, "\n")

# Hard-coded NEBULA settings
ncore <- 8
chunk_size <- 2000

cat("Using ncore:", ncore, "\n")
cat("Using chunk size:", chunk_size, "\n")

root_path <- "/mmfs1/gscratch/togo/leidholt/"
keys <- fromJSON("/mmfs1/home/leidholt/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)


#Load Seurat object

pb90_multiomics_subset <- s3readRDS(
  object = "data/pb90_multiomics_SLL_subset_20260527.rds",
  bucket = "triad",
  region = ""
)


#Metadata
meta_nebula <- pb90_multiomics_subset@meta.data %>%
  mutate(
    cell_barcode = rownames(pb90_multiomics_subset@meta.data),
    record_id = as.character(record_id),
    group = as.character(group),
    sex = as.character(sex),
    study = as.character(study),
    celltype = as.character(KPMP_celltype_general)
  ) %>%
  filter(
    !is.na(record_id),
    !is.na(group),
    !is.na(age),
    !is.na(sex),
    !is.na(celltype),
    sex %in% c("Female", "Male")
  ) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    sex_num = ifelse(sex == "Male", 1, 0)
  )


#Counts
counts_nebula <- GetAssayData(
  pb90_multiomics_subset,
  assay = "RNA",
  layer = "counts"
)

if (any(duplicated(rownames(counts_nebula)))) {
  cat("Making duplicated gene names unique\n")
  rownames(counts_nebula) <- make.unique(rownames(counts_nebula))
}

counts_nebula <- counts_nebula[, meta_nebula$cell_barcode, drop = FALSE]

rm(pb90_multiomics_subset)
gc()


#contrast setup
contrast_pairs <- list(
  T1D_vs_LC  = c("Lean Control", "Type 1 Diabetes"),
  T2D_vs_LC  = c("Lean Control", "Type 2 Diabetes"),
  T1D_vs_T2D = c("Type 2 Diabetes", "Type 1 Diabetes"),
  OC_vs_LC   = c("Lean Control", "Obese Control"),
  T1D_vs_OC  = c("Obese Control", "Type 1 Diabetes"),
  T2D_vs_OC  = c("Obese Control", "Type 2 Diabetes")
)

pair <- contrast_pairs[[contrast_name]]

if (is.null(pair)) {
  stop("Unknown contrast_name: ", contrast_name)
}


#Filtering parameters
min_cells_per_donor_celltype <- 20
min_donors_per_group <- 3

# Lighter gene filtering only
min_cells_gene <- 100
min_mean_gene <- 0.01


#run one NEBULA chunk
run_nebula_chunk <- function(count_chunk, meta_sub, pred, offset, chunk_id) {
  
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


#Main NEBULA function
run_nebula <- function(cell_name, pair, contrast_name) {
  
  cat("\n=============================================================\n")
  cat("NEBULA:", cell_name, contrast_name, "\n")
  cat("=============================================================\n")
  
  meta_sub <- meta_nebula %>%
    filter(
      celltype == cell_name,
      group %in% pair,
      sex %in% c("Female", "Male")
    )
  
  if (nrow(meta_sub) == 0) {
    stop("No cells for this cell type / contrast")
  }
  
  counts_sub <- counts_nebula[, meta_sub$cell_barcode, drop = FALSE]
  
  # donor-cell filter
  donor_cell_counts <- meta_sub %>%
    count(record_id, group, name = "n_cells_donor")
  
  keep_donors <- donor_cell_counts %>%
    filter(n_cells_donor >= min_cells_per_donor_celltype) %>%
    pull(record_id)
  
  meta_sub <- meta_sub %>%
    filter(record_id %in% keep_donors)
  
  counts_sub <- counts_sub[, meta_sub$cell_barcode, drop = FALSE]
  
  donor_tab <- meta_sub %>%
    distinct(record_id, group) %>%
    count(group)
  
  cat("Donor table:\n")
  print(donor_tab)
  
  if (nrow(donor_tab) < 2 || any(donor_tab$n < min_donors_per_group)) {
    stop("Too few donors per group after donor-cell filtering")
  }
  
  # Raw library-size offset.
  # Do NOT log-transform; NEBULA handles this internally.
  lib_size <- Matrix::colSums(counts_sub)
  
  meta_sub <- meta_sub %>%
    mutate(
      lib_size = lib_size
    ) %>%
    filter(is.finite(lib_size), lib_size > 0)
  
  counts_sub <- counts_sub[, meta_sub$cell_barcode, drop = FALSE]
  offset <- meta_sub$lib_size
  
  # lighter gene filtering
  keep_genes <- Matrix::rowSums(counts_sub > 0) >= min_cells_gene &
    Matrix::rowMeans(counts_sub) > min_mean_gene
  
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  
  if (nrow(counts_sub) == 0) {
    stop("No genes after light gene filtering")
  }
  
  meta_sub <- meta_sub %>%
    mutate(
      group = factor(group, levels = pair),
      group_contrast = ifelse(group == pair[2], 1, 0),
      sex = factor(sex, levels = c("Female", "Male")),
      sex_num = ifelse(sex == "Male", 1, 0),
      record_id = factor(record_id)
    )
  
  pred <- model.matrix(
    ~ group_contrast + age + sex_num,
    data = meta_sub
  )
  
  cat("Cells:", nrow(meta_sub), "\n")
  cat("Donors:", length(unique(meta_sub$record_id)), "\n")
  cat("Genes tested:", nrow(counts_sub), "\n")
  cat("Groups:", pair[1], "vs", pair[2], "\n")
  
  gene_chunks <- split(
    rownames(counts_sub),
    ceiling(seq_along(rownames(counts_sub)) / chunk_size)
  )
  
  cat("Number of chunks:", length(gene_chunks), "\n")
  
  start_time <- Sys.time()
  
  chunk_results <- map2(
    gene_chunks,
    seq_along(gene_chunks),
    function(gene_set, chunk_id) {
      cat("Running chunk", chunk_id, "of", length(gene_chunks),
          "with", length(gene_set), "genes\n")
      
      run_nebula_chunk(
        count_chunk = counts_sub[gene_set, , drop = FALSE],
        meta_sub = meta_sub,
        pred = pred,
        offset = offset,
        chunk_id = chunk_id
      )
    }
  )
  
  end_time <- Sys.time()
  cat("NEBULA runtime:\n")
  print(end_time - start_time)
  
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
    return(list(
      results = NULL,
      diagnostics = diagnostics,
      donor_table = donor_tab
    ))
  }
  
  tt <- map_dfr(
    successful_chunks,
    function(x) {
      
      df <- as.data.frame(x$fit$summary)
      
      if (!"gene" %in% colnames(df)) {
        df <- df %>%
          rownames_to_column("gene")
      }
      
      if (!is.null(x$fit$convergence)) {
        conv_df <- tibble(
          gene = rownames(x$fit$summary),
          convergence = x$fit$convergence
        )
        
        df <- df %>%
          left_join(conv_df, by = "gene")
      } else {
        df$convergence <- NA_real_
      }
      
      df %>%
        mutate(chunk_id = x$chunk_id)
    }
  )
  
  # Post-fit convergence filter
  if ("convergence" %in% colnames(tt)) {
    tt <- tt %>%
      mutate(
        stable_fit = is.na(convergence) | convergence >= -10
      )
  } else {
    tt <- tt %>%
      mutate(stable_fit = NA)
  }
  
  if (!"p_group_contrast" %in% colnames(tt)) {
    stop(
      "Could not find p_group_contrast. Available columns: ",
      paste(colnames(tt), collapse = ", ")
    )
  }
  
  tt <- tt %>%
    mutate(
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
  
  list(
    results = tt,
    diagnostics = diagnostics,
    donor_table = donor_tab,
    gene_filter_summary = list(
      n_genes_after_light_filter = nrow(counts_sub)
    )
  )
}


# Run
nebula_out <- run_nebula(
  cell_name = cell_name,
  pair = pair,
  contrast_name = contrast_name
)


# Save to S3
tmp_rds <- tempfile(fileext = ".rds")

saveRDS(nebula_out, file = tmp_rds)

out_object <- paste0(
  "results/nebula/nebula_4/",
  "nebula_",
  cell_name,
  "_",
  contrast_name,
  ".rds"
)

cat("Saving to:", out_object, "\n")

put_object(
  file = tmp_rds,
  object = out_object,
  bucket = "triad",
  base_url = "s3.kopah.uw.edu",
  region = "",
  use_https = TRUE
)

unlink(tmp_rds)

cat("Done\n")