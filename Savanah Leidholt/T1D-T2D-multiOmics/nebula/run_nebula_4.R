#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

cell_name <- args[1]
contrast_name <- args[2]

.libPaths("/mmfs1/gscratch/togo/leidholt/R_SLL_Seurat/")

library(Seurat)
library(dplyr)
library(Matrix)
library(nebula)
library(future)
library(purrr)
library(aws.s3)
library(jsonlite)

future::plan(sequential)

user <- Sys.info()[["user"]]

root_path <- "/mmfs1/gscratch/togo/leidholt/"
keys <- fromJSON("/mmfs1/home/leidholt/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

pb90_multiomics_subset <- s3readRDS(
  object = "data/pb90_multiomics_SLL_subset_20260527.rds",
  bucket = "triad",
  region = ""
)

#Clean metadata
meta_nebula <- pb90_multiomics_subset@meta.data %>%
  dplyr::mutate(
    cell_barcode = rownames(pb90_multiomics_subset@meta.data),
    record_id = as.character(record_id),
    group = as.character(group),
    sex = as.character(sex),
    study = as.character(study),
    celltype = as.character(KPMP_celltype_general)
  ) %>%
  dplyr::filter(
    !is.na(record_id),
    !is.na(group),
    !is.na(age),
    !is.na(sex),
    !is.na(celltype),
    sex %in% c("Female", "Male")
  ) %>%
  dplyr::mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    sex_num = ifelse(sex == "Male", 1, 0)
  )

#Extract full raw count matrix
counts_nebula <- GetAssayData(
  pb90_multiomics_subset,
  assay = "RNA",
  layer = "counts"
)

if (any(duplicated(rownames(counts_nebula)))) {
  message("Making duplicated gene names unique")
  rownames(counts_nebula) <- make.unique(rownames(counts_nebula))
}

counts_nebula <- counts_nebula[, meta_nebula$cell_barcode, drop = FALSE]

contrast_pairs <- list(
  T1D_vs_LC = c("Lean Control", "Type 1 Diabetes"),
  T2D_vs_LC = c("Lean Control", "Type 2 Diabetes"),
  T1D_vs_T2D = c("Type 2 Diabetes", "Type 1 Diabetes"),
  OC_vs_LC = c("Lean Control", "Obese Control"),
  T1D_vs_OC = c("Obese Control", "Type 1 Diabetes"),
  T2D_vs_OC = c("Obese Control", "Type 2 Diabetes")
)

pair <- contrast_pairs[[contrast_name]]

if (is.null(pair)) {
  stop("Unknown contrast_name: ", contrast_name)
}

min_cells_per_donor_celltype <- 20
min_donors_per_group <- 3
min_cells_gene <- 50
min_mean_gene <- 0.01

#------------------------------------
# NEBULA PIPELINE
#------------------------------------
run_nebula <- function(cell_name, pair, contrast_name) {
  
  cat("\n=== NEBULA:", cell_name, contrast_name, "===\n")
  
  #Subset metadata to target cell type and contrast groups
  meta_sub <- meta_nebula %>%
    dplyr::filter(
      celltype == cell_name,
      group %in% pair,
      sex %in% c("Female", "Male")
    )
  
  if (nrow(meta_sub) == 0) {
    message("Skipping: no cells")
    return(NULL)
  }
  
  #Subset counts using al genes
  counts_sub <- counts_nebula[, meta_sub$cell_barcode, drop = FALSE]
  
  #Donor cell count filter
  donor_cell_counts <- meta_sub %>%
    dplyr::count(record_id, group, name = "n_cells_donor")
  
  keep_donors <- donor_cell_counts %>%
    dplyr::filter(n_cells_donor >= min_cells_per_donor_celltype) %>%
    dplyr::pull(record_id)
  
  meta_sub <- meta_sub %>%
    dplyr::filter(record_id %in% keep_donors)
  
  counts_sub <- counts_sub[, meta_sub$cell_barcode, drop = FALSE]
  
  donor_tab <- meta_sub %>%
    dplyr::distinct(record_id, group) %>%
    dplyr::count(group)
  
  if (nrow(donor_tab) < 2 || any(donor_tab$n < min_donors_per_group)) {
    message("Skipping: too few donors per group")
    print(donor_tab)
    return(NULL)
  }
  
  #Moved this step from version run_nebula_3 before filtering so that the entire library was used
  meta_sub <- meta_sub %>%
    dplyr::mutate(
      lib_size = Matrix::colSums(counts_sub),
      pooled_offset = log(lib_size)
    ) %>%
    dplyr::filter(is.finite(pooled_offset), lib_size > 0)
  
  counts_sub <- counts_sub[, meta_sub$cell_barcode, drop = FALSE]
  
  #lowly expressed gene filtering
  keep_genes <- Matrix::rowSums(counts_sub > 0) >= min_cells_gene &
    Matrix::rowMeans(counts_sub) > min_mean_gene
  
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  
  if (nrow(counts_sub) == 0) {
    message("Skipping: no genes after filtering")
    return(NULL)
  }
  
  #metadata factor setting
  meta_sub <- meta_sub %>%
    dplyr::mutate(
      group = factor(group, levels = pair),
      group_contrast = ifelse(group == pair[2], 1, 0),
      sex = factor(sex, levels = c("Female", "Male")),
      sex_num = ifelse(sex == "Male", 1, 0),
      record_id = factor(record_id)
    )
  
  gene_index <- seq_len(nrow(counts_sub))
  genes_list <- rownames(counts_sub)
  
  cat(
    "Testing",
    length(genes_list),
    "genes with",
    length(unique(meta_sub$record_id)),
    "subjects and",
    nrow(meta_sub),
    "cells.\n"
  )
  
  start_time <- Sys.time()
  
  nebula_gene_results <- purrr::map(
    gene_index,
    function(i) {
      g <- genes_list[i]
      warn <- NULL
      err <- NULL
      res <- NULL
      
      tryCatch({
        count_gene <- counts_sub[i, , drop = FALSE]
        rownames(count_gene) <- g
        
        pred_gene <- model.matrix(
          ~ group_contrast + age + sex_num,
          data = meta_sub
        )
        
        res <- withCallingHandlers({
          nebula(
            count = count_gene,
            id = meta_sub$record_id,
            pred = pred_gene,
            offset = meta_sub$pooled_offset,
            model = "NBLMM",
            covariance = TRUE,
            reml = 1,
            output_re = TRUE,
            ncore = 1
          )
        }, warning = function(w) {
          warn <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        })
        
      }, error = function(e) {
        err <<- conditionMessage(e)
      })
      
      list(
        gene = g,
        result = res,
        warning = warn,
        error = err
      )
    }
  )
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # Diagnostics
  diagnostics <- purrr::map_dfr(
    nebula_gene_results,
    function(x) {
      tibble::tibble(
        gene = x$gene,
        has_result = !is.null(x$result),
        warning = x$warning,
        error = x$error
      )
    }
  )
  
  print(diagnostics %>% dplyr::count(has_result))
  
  #results table
  successful_results <- nebula_gene_results %>%
    purrr::keep(~ !is.null(.x$result))
  
  if (length(successful_results) == 0) {
    message("No successful genes")
    return(list(results = NULL, diagnostics = diagnostics))
  }
  
  tt <- purrr::map_dfr(
    successful_results,
    function(x) {
      df <- x$result$summary %>% as.data.frame()
      df$gene_tested <- x$gene
      df
    }
  )
  
  if (!"p_group_contrast" %in% colnames(tt)) {
    stop("Could not find p_group_contrast. Run colnames(tt) to inspect NEBULA output.")
  }
  
  tt <- tt %>%
    dplyr::mutate(
      comparison = contrast_name,
      celltype = cell_name,
      contrast_level_0 = pair[1],
      contrast_level_1 = pair[2],
      logFC_direction = paste0(pair[2], "_vs_", pair[1]),
      n_cells = nrow(meta_sub),
      n_genes_tested = length(genes_list),
      n_genes_successful = nrow(tt),
      n_donors_group1 = donor_tab$n[match(pair[1], donor_tab$group)],
      n_donors_group2 = donor_tab$n[match(pair[2], donor_tab$group)],
      p_adj = p.adjust(p_group_contrast, method = "fdr"),
      neg_log10_p = -log10(pmax(p_group_contrast, 1e-300)),
      neg_log10_fdr = -log10(pmax(p_adj, 1e-300))
    )
  
  return(list(results = tt, diagnostics = diagnostics))
}

nebula_out <- run_nebula(
  cell_name = cell_name,
  pair = pair,
  contrast_name = contrast_name
)

if (is.null(nebula_out)) {
  stop("NEBULA returned NULL for ", cell_name, " ", contrast_name)
}

tmp_rds <- tempfile(fileext = ".rds")
saveRDS(nebula_out, file = tmp_rds)

put_object(
  file = tmp_rds,
  object = paste0("results/nebula/nebula_4/","nebula_", cell_name,"_",contrast_name,".rds"),
  bucket = "triad",
  base_url = "s3.kopah.uw.edu",
  region = "",
  use_https = TRUE
)

unlink(tmp_rds)
