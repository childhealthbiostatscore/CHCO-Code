#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
cell_name <- args[1]

if (is.na(cell_name)) {
  stop("Usage: Rscript prep_nebula_celltype.R <cell_name>")
}

.libPaths("/mmfs1/gscratch/togo/leidholt/R_SLL_Seurat/")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(aws.s3)
  library(jsonlite)
})

keys <- jsonlite::fromJSON("/mmfs1/home/leidholt/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

cat("Preparing NEBULA input for cell type:", cell_name, "\n")

cat("Loading Seurat object from S3\n")

obj <- s3readRDS(
  object = "data/pb90_multiomics_SLL_subset_20260527.rds",
  bucket = "triad",
  region = "",
  base_url = "s3.kopah.uw.edu",
  use_https = TRUE
)

meta <- obj@meta.data %>%
  mutate(
    cell_barcode = rownames(obj@meta.data),
    record_id = as.character(record_id),
    group = as.character(group),
    sex = as.character(sex),
    study = as.character(study),
    celltype = as.character(KPMP_celltype_general)
  ) %>%
  filter(
    celltype == cell_name,
    !is.na(record_id),
    !is.na(group),
    !is.na(age),
    !is.na(sex),
    !is.na(celltype),
    sex %in% c("Female", "Male")
  )

if (nrow(meta) == 0) {
  stop("No cells found for cell type: ", cell_name)
}

cat("Cells after metadata filtering:", nrow(meta), "\n")

counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

if (any(duplicated(rownames(counts)))) {
  cat("Making duplicated gene names unique\n")
  rownames(counts) <- make.unique(rownames(counts))
}

counts <- counts[, meta$cell_barcode, drop = FALSE]

stopifnot(identical(colnames(counts), meta$cell_barcode))

rm(obj)
gc()

prep <- list(
  count = counts,
  meta = meta,
  cell_name = cell_name,
  n_cells = nrow(meta),
  n_genes = nrow(counts),
  created = Sys.time()
)

tmp_rds <- tempfile(fileext = ".rds")
saveRDS(prep, tmp_rds)

out_object <- paste0(
  "data/nebula_prepped_inputs/",
  "nebula_input_",
  cell_name,
  ".rds"
)

cat("Saving prepped input to S3:", out_object, "\n")

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