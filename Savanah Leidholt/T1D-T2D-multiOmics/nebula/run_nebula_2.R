#!/usr/bin/env Rscript

# NEBULA single-cell DE by cell type and contrast
# Uses entire count matrix, donor-aware gene filtering, and chunked NEBULA fitting

args <- commandArgs(trailingOnly = TRUE)

cell_name <- args[1]
contrast_name <- args[2]

.libPaths("/mmfs1/gscratch/togo/leidholt/R_SLL_Seurat/")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(nebula)
  library(purrr)
  library(aws.s3)
  library(jsonlite)
  library(tidyr)
  library(tibble)
})


#Paths
root_path <- "/mmfs1/gscratch/togo/leidholt/"
keys <- fromJSON("/mmfs1/home/leidholt/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)


#Settings
min_cells_per_donor_celltype <- 20
min_donors_per_group <- 3

min_cells_gene <- 50
min_mean_gene <- 0.01
min_donors_detected_per_group <- 2

max_abs_logFC_for_stable <- 20
nebula_chunk_size <- 500

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


#Read data
pb90_multiomics_subset <- s3readRDS(
  object = "data/pb90_multiomics_SLL_subset_20260527.rds",
  bucket = "triad",
  region = ""
)

counts_all <- GetAssayData(
  pb90_multiomics_subset,
  assay = "RNA",
  layer = "counts"
)

if (any(duplicated(rownames(counts_all)))) {
  message("Making duplicated gene names unique")
  rownames(counts_all) <- make.unique(rownames(counts_all))
}

full_lib_size <- Matrix::colSums(counts_all)

meta_all <- pb90_multiomics_subset@meta.data %>%
  mutate(
    cell_barcode = rownames(pb90_multiomics_subset@meta.data),
    record_id = as.character(record_id),
    group = as.character(group),
    sex = as.character(sex),
    study = as.character(study),
    celltype = as.character(KPMP_celltype_general),
    full_lib_size = full_lib_size[cell_barcode]
  ) %>%
  filter(
    !is.na(record_id),
    !is.na(group),
    !is.na(age),
    !is.na(sex),
    !is.na(celltype),
    sex %in% c("Female", "Male"),
    is.finite(full_lib_size),
    full_lib_size > 0
  )

counts_all <- counts_all[, meta_all$cell_barcode, drop = FALSE]


#Donor-aware gene filter
filter_genes_donor_aware <- function(counts_sub, meta_sub, pair) {
  
  basic_keep <- Matrix::rowSums(counts_sub > 0) >= min_cells_gene &
    Matrix::rowMeans(counts_sub) > min_mean_gene
  
  counts_basic <- counts_sub[basic_keep, , drop = FALSE]
  
  if (nrow(counts_basic) == 0) {
    return(character(0))
  }
  
  detected_mat <- counts_basic > 0
  
  donor_detect <- map_dfr(
    rownames(counts_basic),
    function(g) {
      tibble(
        gene = g,
        record_id = meta_sub$record_id,
        group = meta_sub$group,
        detected = as.logical(detected_mat[g, ])
      ) %>%
        group_by(gene, record_id, group) %>%
        summarise(detected_donor = any(detected), .groups = "drop")
    }
  )
  
  gene_prev <- donor_detect %>%
    group_by(gene, group) %>%
    summarise(
      n_donors_detected = sum(detected_donor),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = group,
      values_from = n_donors_detected,
      values_fill = 0
    )
  
  keep_genes <- gene_prev %>%
    filter(
      .data[[pair[1]]] >= min_donors_detected_per_group,
      .data[[pair[2]]] >= min_donors_detected_per_group
    ) %>%
    pull(gene)
  
  keep_genes
}


#nebula function to run in chunks and drop one gene until it passes at a time
run_nebula_keep_going <- function(counts_sub, meta_sub, pred, chunk_size = 500) {
  
  fit_genes <- function(genes) {
    
    fit <- tryCatch(
      nebula(
        count = counts_sub[genes, , drop = FALSE],
        id = meta_sub$record_id,
        pred = pred,
        offset = meta_sub$pooled_offset,
        model = "NBLMM",
        covariance = TRUE,
        reml = 1,
        output_re = TRUE,
        ncore = 1
      ),
      error = function(e) e
    )
    
    if (!inherits(fit, "error")) {
      return(list(
        results = as.data.frame(fit$summary) %>%
          tibble::rownames_to_column("gene"),
        failed_genes = character(0)
      ))
    }
    
    if (length(genes) == 1) {
      message("Dropping failed gene: ", genes, " | ", fit$message)
      return(list(
        results = NULL,
        failed_genes = genes
      ))
    }
    
    mid <- floor(length(genes) / 2)
    
    left <- fit_genes(genes[seq_len(mid)])
    right <- fit_genes(genes[(mid + 1):length(genes)])
    
    list(
      results = dplyr::bind_rows(left$results, right$results),
      failed_genes = c(left$failed_genes, right$failed_genes)
    )
  }
  
  genes <- rownames(counts_sub)
  chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
  
  out <- purrr::map(chunks, fit_genes)
  
  list(
    results = dplyr::bind_rows(purrr::map(out, "results")),
    failed_genes = unlist(purrr::map(out, "failed_genes"))
  )
}


#NEBULA function
run_nebula <- function(cell_name, pair, contrast_name) {
  
  cat("\n=== NEBULA:", cell_name, contrast_name, "===\n")
  cat("Reference:", pair[1], "\n")
  cat("Contrast:", pair[2], "\n\n")
  
  meta_sub <- meta_all %>%
    filter(
      celltype == cell_name,
      group %in% pair
    )
  
  if (nrow(meta_sub) == 0) {
    stop("No cells found for ", cell_name, " and ", contrast_name)
  }
  
  counts_sub <- counts_all[, meta_sub$cell_barcode, drop = FALSE]
  
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
  
  print(donor_tab)
  
  if (nrow(donor_tab) < 2 || any(donor_tab$n < min_donors_per_group)) {
    stop("Too few donors per group after filtering")
  }
  
  meta_sub <- meta_sub %>%
    mutate(
      group = factor(group, levels = pair),
      group_contrast = ifelse(group == pair[2], 1, 0),
      sex = factor(sex, levels = c("Female", "Male")),
      sex_num = ifelse(sex == "Male", 1, 0),
      record_id = factor(record_id),
      age_scaled = as.numeric(scale(age)),
      pooled_offset = log(full_lib_size)
    ) %>%
    filter(
      is.finite(pooled_offset),
      is.finite(age_scaled)
    )
  
  counts_sub <- counts_sub[, meta_sub$cell_barcode, drop = FALSE]
  
  keep_genes <- filter_genes_donor_aware(
    counts_sub = counts_sub,
    meta_sub = meta_sub,
    pair = pair
  )
  
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  
  if (nrow(counts_sub) == 0) {
    stop("No genes retained after donor-aware filtering")
  }
  
  genes_list <- rownames(counts_sub)
  
  cat(
    "Testing", length(genes_list), "genes with",
    length(unique(meta_sub$record_id)), "donors and",
    nrow(meta_sub), "cells.\n"
  )
  
  pred <- model.matrix(
    ~ group_contrast + age_scaled + sex_num,
    data = meta_sub
  )
  
  nebula_safe <- run_nebula_keep_going(
    counts_sub = counts_sub,
    meta_sub = meta_sub,
    pred = pred,
    chunk_size = nebula_chunk_size
  )
  
  tt <- nebula_safe$results
  failed_genes <- nebula_safe$failed_genes
  
  if (is.null(tt) || nrow(tt) == 0) {
    stop("No genes successfully fit")
  }
  
  cat("Successfully fit genes:", nrow(tt), "\n")
  cat("Failed genes dropped:", length(failed_genes), "\n")
  
  diagnostics <- tibble(
    gene = c(tt$gene, failed_genes),
    has_result = c(rep(TRUE, nrow(tt)), rep(FALSE, length(failed_genes))),
    error = c(rep(NA_character_, nrow(tt)), rep("NEBULA failed; dropped by recursive chunking", length(failed_genes)))
  )
  
  print(diagnostics %>% count(has_result))
  
  if (!"p_group_contrast" %in% colnames(tt)) {
    stop("Could not find p_group_contrast. Inspect colnames(tt).")
  }
  
  tt <- tt %>%
    mutate(
      comparison = contrast_name,
      celltype = cell_name,
      contrast_level_0 = pair[1],
      contrast_level_1 = pair[2],
      logFC_direction = paste0(pair[2], "_vs_", pair[1]),
      n_cells = nrow(meta_sub),
      n_genes_tested = length(genes_list),
      n_genes_successful = nrow(tt),
      n_genes_failed = length(failed_genes),
      n_donors_group1 = donor_tab$n[match(pair[1], donor_tab$group)],
      n_donors_group2 = donor_tab$n[match(pair[2], donor_tab$group)],
      p_adj = p.adjust(p_group_contrast, method = "fdr"),
      neg_log10_p = -log10(pmax(p_group_contrast, 1e-300)),
      neg_log10_fdr = -log10(pmax(p_adj, 1e-300)),
      stable_fit = case_when(
        !is.finite(logFC_group_contrast) ~ FALSE,
        !is.finite(se_group_contrast) ~ FALSE,
        is.na(p_adj) ~ FALSE,
        abs(logFC_group_contrast) > max_abs_logFC_for_stable ~ FALSE,
        TRUE ~ TRUE
      ),
      result_class = case_when(
        stable_fit & p_adj < 0.05 & logFC_group_contrast > 0 ~ paste0("Increased in ", pair[2]),
        stable_fit & p_adj < 0.05 & logFC_group_contrast < 0 ~ paste0("Decreased in ", pair[2]),
        stable_fit ~ "Not significant",
        TRUE ~ "Unstable fit"
      )
    )
  
  list(
    results = tt,
    diagnostics = diagnostics,
    failed_genes = failed_genes,
    settings = list(
      cell_name = cell_name,
      contrast_name = contrast_name,
      pair = pair,
      min_cells_per_donor_celltype = min_cells_per_donor_celltype,
      min_donors_per_group = min_donors_per_group,
      min_cells_gene = min_cells_gene,
      min_mean_gene = min_mean_gene,
      min_donors_detected_per_group = min_donors_detected_per_group,
      max_abs_logFC_for_stable = max_abs_logFC_for_stable,
      nebula_chunk_size = nebula_chunk_size,
      offset = "full RNA library size",
      age = "scaled",
      gene_set = "all genes after donor-aware filtering, not HVGs"
    )
  )
}


nebula_out <- run_nebula(
  cell_name = cell_name,
  pair = pair,
  contrast_name = contrast_name
)



tmp_rds <- tempfile(fileext = ".rds")

saveRDS(nebula_out, file = tmp_rds)

put_object(
  file = tmp_rds,
  object = paste0(
    "results/nebula/nebula_3/",
    "nebula_",
    cell_name,
    "_",
    contrast_name,
    ".rds"
  ),
  bucket = "triad",
  base_url = "s3.kopah.uw.edu",
  region = "",
  use_https = TRUE
)

unlink(tmp_rds)

cat("\nSaved NEBULA results for", cell_name, contrast_name, "\n")