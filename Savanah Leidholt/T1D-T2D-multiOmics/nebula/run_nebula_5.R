#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
cell_name <- args[1]
contrast_name <- args[2]

if (is.na(cell_name) || is.na(contrast_name)) {
  stop("Usage: Rscript run_nebula_5_hvg.R <cell_name> <contrast_name>")
}

.libPaths(c("/mmfs1/gscratch/togo/leidholt/R_SLL_Seurat", .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(nebula)
  library(dplyr)
  library(tibble)
  library(aws.s3)
  library(jsonlite)
})

ncore <- min(8, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
chunk_size <- 1000

min_cells_per_donor_celltype <- 20
min_donors_per_group <- 3
min_cells_gene <- 50
min_mean_gene <- 0.005

keys <- jsonlite::fromJSON("/mmfs1/home/leidholt/keys.json")

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

if (is.null(pair)) {
  stop("Unknown contrast_name: ", contrast_name)
}

in_object <- paste0(
  "data/nebula_prepped_inputs_hvg/",
  "nebula_input_hvg_",
  cell_name,
  ".rds"
)

cat("Reading HVG prepped input from S3:", in_object, "\n")

prep <- s3readRDS(
  object = in_object,
  bucket = "triad",
  region = "",
  base_url = "s3.kopah.uw.edu",
  use_https = TRUE
)

counts_nebula <- prep$count
meta_nebula <- prep$meta

meta_sub <- meta_nebula %>%
  filter(
    group %in% pair,
    !is.na(record_id),
    !is.na(group),
    !is.na(age),
    !is.na(sex),
    sex %in% c("Female", "Male")
  )

counts_sub <- counts_nebula[, meta_sub$cell_barcode, drop = FALSE]

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

lib_size <- Matrix::colSums(counts_sub)

meta_sub <- meta_sub %>%
  mutate(lib_size = lib_size) %>%
  filter(is.finite(lib_size), lib_size > 0)

counts_sub <- counts_sub[, meta_sub$cell_barcode, drop = FALSE]
offset <- meta_sub$lib_size

keep_genes <- Matrix::rowSums(counts_sub > 0) >= min_cells_gene &
  Matrix::rowMeans(counts_sub) > min_mean_gene

counts_sub <- counts_sub[keep_genes, , drop = FALSE]

if (nrow(counts_sub) == 0) {
  stop("No HVGs left after filtering")
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

cat("Running HVG NEBULA\n")
cat("Cell type:", cell_name, "\n")
cat("Contrast:", contrast_name, "\n")
cat("Cells:", nrow(meta_sub), "\n")
cat("Donors:", length(unique(meta_sub$record_id)), "\n")
cat("HVGs tested:", nrow(counts_sub), "\n")
cat("ncore:", ncore, "\n")

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

chunk_results <- vector("list", length(gene_chunks))

for (chunk_id in seq_along(gene_chunks)) {
  gene_set <- gene_chunks[[chunk_id]]
  cat("Running chunk", chunk_id, "of", length(gene_chunks), "\n")
  
  chunk_results[[chunk_id]] <- run_nebula_chunk(
    count_chunk = counts_sub[gene_set, , drop = FALSE],
    chunk_id = chunk_id
  )
  
  gc()
}

diagnostics <- bind_rows(
  lapply(chunk_results, function(x) {
    tibble(
      chunk_id = x$chunk_id,
      gene = x$genes,
      chunk_success = !is.null(x$fit),
      warning = x$warning,
      error = x$error
    )
  })
)

successful_chunks <- Filter(function(x) !is.null(x$fit), chunk_results)

if (length(successful_chunks) == 0) {
  
  tt <- NULL
  
} else {
  
  tt <- bind_rows(
    lapply(successful_chunks, function(x) {
      
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
      
      df %>%
        mutate(chunk_id = x$chunk_id)
    })
  )
  
  if (!"p_group_contrast" %in% colnames(tt)) {
    stop(
      "Could not find p_group_contrast. Available columns: ",
      paste(colnames(tt), collapse = ", ")
    )
  }
  
  tt <- tt %>%
    mutate(
      stable_fit = is.na(convergence) | convergence >= -10,
      comparison = contrast_name,
      celltype = cell_name,
      hvg_only = TRUE,
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
}

convergence_summary <- tibble(
  celltype = cell_name,
  comparison = contrast_name,
  hvg_only = TRUE,
  n_cells = nrow(meta_sub),
  n_donors = length(unique(meta_sub$record_id)),
  n_genes_tested = nrow(counts_sub),
  n_chunks = length(gene_chunks),
  n_chunks_successful = length(successful_chunks),
  n_genes_returned = ifelse(is.null(tt), 0, n_distinct(tt$gene)),
  n_genes_stable = ifelse(is.null(tt), 0, sum(tt$stable_fit, na.rm = TRUE)),
  convergence_rate = ifelse(
    is.null(tt) || nrow(tt) == 0,
    0,
    mean(tt$stable_fit, na.rm = TRUE)
  )
)

nebula_out <- list(
  results = tt,
  diagnostics = diagnostics,
  donor_table = donor_tab,
  convergence_summary = convergence_summary,
  gene_filter_summary = list(
    n_hvg_original = length(prep$hvg),
    n_hvg_after_filter = nrow(counts_sub)
  )
)

tmp_rds <- tempfile(fileext = ".rds")
saveRDS(nebula_out, tmp_rds)

out_object <- paste0(
  "results/nebula/nebula_5_hvg/",
  "nebula_hvg_",
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