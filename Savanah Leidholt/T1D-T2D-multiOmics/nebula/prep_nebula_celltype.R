#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
cell_name <- args[1]
contrast_name <- args[2]

if (is.na(cell_name) || is.na(contrast_name)) {
  stop("Usage: Rscript prep_nebula_inputs.R <cell_name> <contrast_name>")
}

.libPaths("/mmfs1/gscratch/togo/leidholt/R_SLL_Seurat/")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(aws.s3)
  library(jsonlite)
})

root_path <- "/mmfs1/gscratch/togo/leidholt/"
out_dir <- file.path(root_path, "nebula_prepped_inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keys <- fromJSON("/mmfs1/home/leidholt/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

contrast_pairs <- list(
  T1D_vs_LC  = c("Lean Control", "Type 1 Diabetes"),
  T2D_vs_LC  = c("Lean Control", "Type 2 Diabetes"),
  T1D_vs_T2D = c("Type 2 Diabetes", "Type 1 Diabetes"),
  OC_vs_LC   = c("Lean Control", "Obese Control"),
  T1D_vs_OC  = c("Obese Control", "Type 1 Diabetes"),
  T2D_vs_OC  = c("Obese Control", "Type 2 Diabetes")
)

pair <- contrast_pairs[[contrast_name]]
if (is.null(pair)) stop("Unknown contrast_name: ", contrast_name)

min_cells_per_donor_celltype <- 20
min_donors_per_group <- 3
min_cells_gene <- 100
min_mean_gene <- 0.01

cat("Loading Seurat object\n")

obj <- s3readRDS(
  object = "data/pb90_multiomics_SLL_subset_20260527.rds",
  bucket = "triad",
  region = ""
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
    group %in% pair,
    !is.na(record_id),
    !is.na(group),
    !is.na(age),
    !is.na(sex),
    sex %in% c("Female", "Male")
  )

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

if (any(duplicated(rownames(counts)))) {
  rownames(counts) <- make.unique(rownames(counts))
}

counts <- counts[, meta$cell_barcode, drop = FALSE]

rm(obj)
gc()

donor_cell_counts <- meta %>%
  count(record_id, group, name = "n_cells_donor")

keep_donors <- donor_cell_counts %>%
  filter(n_cells_donor >= min_cells_per_donor_celltype) %>%
  pull(record_id)

meta <- meta %>%
  filter(record_id %in% keep_donors)

counts <- counts[, meta$cell_barcode, drop = FALSE]

donor_tab <- meta %>%
  distinct(record_id, group) %>%
  count(group)

print(donor_tab)

if (nrow(donor_tab) < 2 || any(donor_tab$n < min_donors_per_group)) {
  stop("Too few donors per group after filtering")
}

lib_size <- Matrix::colSums(counts)

meta <- meta %>%
  mutate(lib_size = lib_size) %>%
  filter(is.finite(lib_size), lib_size > 0)

counts <- counts[, meta$cell_barcode, drop = FALSE]
offset <- meta$lib_size

keep_genes <- Matrix::rowSums(counts > 0) >= min_cells_gene &
  Matrix::rowMeans(counts) > min_mean_gene

counts <- counts[keep_genes, , drop = FALSE]

if (nrow(counts) == 0) {
  stop("No genes left after filtering")
}

meta <- meta %>%
  mutate(
    group = factor(group, levels = pair),
    group_contrast = ifelse(group == pair[2], 1, 0),
    sex = factor(sex, levels = c("Female", "Male")),
    sex_num = ifelse(sex == "Male", 1, 0),
    record_id = factor(record_id)
  )

pred <- model.matrix(
  ~ group_contrast + age + sex_num,
  data = meta
)

prep <- list(
  count = counts_celltype,
  meta = meta_celltype,
  cell_name = cell_name
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
