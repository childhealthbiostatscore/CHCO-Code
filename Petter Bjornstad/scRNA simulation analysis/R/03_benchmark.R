################################################################################
# 03_benchmark.R
#
# PURPOSE:
#   Aggregate results from all completed simulation runs (output of
#   02_simulate_analyze.R) and compute per-scenario x method summaries:
#     - Power        (sensitivity at BH FDR 5%)
#     - Observed FDR (FP / called significant)
#     - Type I error (FP rate when prop_de == 0)
#     - AUC-ROC      (gene-level p-value as score)
#     - P-value calibration (KS statistic vs Uniform under null)
#     - LogFC bias   (median bias & MAE of estimated vs true log2FC, DE genes)
#     - Compute time (wall seconds per method)
#
# INPUT (from S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/
#     param_grid/param_grid.rds
#     stats/array_NNNNN/{nebula_stats,edger_stats,limma_stats,
#                        wilcox_stats,mast_stats,timing,params}.rds
#
# OUTPUT (to S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/benchmark/
#     benchmark_raw.rds      -- one row per array_id x method (all reps)
#     benchmark_avg.rds      -- averaged/SD across reps per scenario x method
#     benchmark_avg.csv      -- same, human-readable
#
# USAGE:
#   Rscript 03_benchmark.R \
#       --fdr_thr  0.05 \
#       --n_cores  32
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(optparse)
  library(BiocParallel)
  library(ROCR)
  library(aws.s3)
})

# ── S3 / Multi-user setup ─────────────────────────────────────────────────────
setup_s3 <- function() {
  user <- Sys.info()[["user"]]
  
  if (user == "choiyej") {
    # --- local (choiyej Mac) ---
    keys_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json"
  } else if (user %in% c("rameshsh", "yejichoi", "pylell")) {
    # --- Hyak HPC ---
    keys_path <- "/mmfs1/home/yejichoi/keys.json"
  } else {
    stop("Unknown user '", user, "'. Add credentials path to setup_s3().")
  }
  
  keys <- jsonlite::fromJSON(keys_path)
  Sys.setenv(
    AWS_ACCESS_KEY_ID     = keys$MY_ACCESS_KEY,
    AWS_SECRET_ACCESS_KEY = keys$MY_SECRET_KEY,
    AWS_S3_ENDPOINT       = "s3.kopah.uw.edu"
  )
  message(sprintf("S3 configured for user '%s'", user))
}

setup_s3()

# S3 paths
S3_BUCKET    <- "scrna"
S3_BASE      <- "Projects/Paired scRNA simulation analysis/results/"
S3_GRID_PFX  <- paste0(S3_BASE, "param_grid/")
S3_STATS_PFX <- paste0(S3_BASE, "stats/")
S3_BENCH_PFX <- paste0(S3_BASE, "benchmark/")

# ── S3 helper functions ───────────────────────────────────────────────────────
s3write_using_region <- function(FUN, ..., object, bucket,
                                 region = NULL, opts = NULL, filename = NULL) {
  ext  <- if (!is.null(filename)) tools::file_ext(filename) else tools::file_ext(object)
  tmp  <- tempfile(fileext = if (nchar(ext) > 0) paste0(".", ext) else "")
  on.exit(unlink(tmp))
  FUN(..., file = tmp)
  args <- list(file = tmp, object = object, bucket = bucket)
  if (!is.null(region)) args$region <- region
  if (!is.null(opts))   args        <- c(args, opts)
  do.call(aws.s3::put_object, args)
}

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--fdr_thr", type = "double",  default = 0.05),
  make_option("--n_cores", type = "integer", default = 32L)
)
opt <- parse_args(OptionParser(option_list = option_list))
BiocParallel::register(BiocParallel::MulticoreParam(opt$n_cores))

# ── Load param grid from S3 ───────────────────────────────────────────────────
message("── [03] Loading param_grid from S3 ──")
param_grid  <- s3readRDS(
  object = paste0(S3_GRID_PFX, "param_grid.rds"),
  bucket = S3_BUCKET,
  region = ""
)
n_scenarios <- nrow(param_grid)
message(sprintf("── [03] Aggregating %d tasks ──", n_scenarios))

# ── Metric computation ────────────────────────────────────────────────────────
# Methods that test the interaction coefficient (logFC comparable to log2FC_true)
INTERACTION_METHODS <- c("nebula", "edger", "limma")

compute_metrics <- function(stats_df, fdr_thr, method_name) {
  if (is.null(stats_df) || nrow(stats_df) == 0) return(NULL)
  
  df <- stats_df %>%
    mutate(
      padj_int = ifelse(is.na(padj_int), 1, padj_int),
      is_de    = ifelse(is.na(is_de),    FALSE, is_de),
      is_sig   = padj_int < fdr_thr,
      TP = is_sig  &  is_de,
      FP = is_sig  & !is_de,
      FN = !is_sig &  is_de,
      TN = !is_sig & !is_de
    )
  
  n_de    <- sum(df$is_de)
  n_nonde <- sum(!df$is_de)
  n_sig   <- sum(df$is_sig)
  n_tp    <- sum(df$TP)
  n_fp    <- sum(df$FP)
  n_genes_tested <- sum(!is.na(stats_df$pval_int))
  
  # AUC-ROC
  auc_roc <- tryCatch({
    pred_obj <- ROCR::prediction(1 - df$pval_int, df$is_de)
    as.numeric(ROCR::performance(pred_obj, "auc")@y.values)
  }, error = function(e) NA_real_)
  
  # P-value calibration: KS test of non-DE p-values against Uniform(0,1)
  # Only meaningful when there are non-DE genes with valid p-values
  pvals_nonde <- df$pval_int[!df$is_de & !is.na(df$pval_int)]
  ks_stat <- if (length(pvals_nonde) >= 10) {
    ks.test(pvals_nonde, "punif")$statistic
  } else {
    NA_real_
  }
  
  # LogFC bias: only for methods that estimate the interaction coefficient
  # (Wilcoxon/MAST test a different contrast, so bias is not meaningful)
  logfc_bias    <- NA_real_
  logfc_mae     <- NA_real_
  if (method_name %in% INTERACTION_METHODS && n_de > 0 &&
      "log2FC_true" %in% colnames(df)) {
    de_genes <- df %>% filter(is_de, !is.na(logFC_int))
    if (nrow(de_genes) > 0) {
      residuals   <- de_genes$logFC_int - de_genes$log2FC_true
      logfc_bias  <- median(residuals, na.rm = TRUE)
      logfc_mae   <- median(abs(residuals), na.rm = TRUE)
    }
  }
  
  data.frame(
    n_de_true      = n_de,
    n_genes_tested = n_genes_tested,
    n_sig          = n_sig,
    power          = if (n_de    > 0) n_tp / n_de    else NA_real_,
    fdr_obs        = if (n_sig   > 0) n_fp / n_sig   else 0,
    t1e            = if (n_nonde > 0) n_fp / n_nonde else NA_real_,
    auc_roc        = auc_roc,
    ks_stat        = as.numeric(ks_stat),
    logfc_bias     = logfc_bias,
    logfc_mae      = logfc_mae,
    stringsAsFactors = FALSE
  )
}

# ── Helper: safely read one RDS file from S3 (returns NULL on missing/error) ──
s3_read_safe <- function(object, bucket, region = "") {
  tryCatch(
    s3readRDS(object = object, bucket = bucket, region = region),
    error = function(e) NULL
  )
}

# ── Process one array task ────────────────────────────────────────────────────
process_task <- function(i) {
  arr_id  <- param_grid$array_id[i]
  arr_str <- sprintf("array_%05d", arr_id)
  pfx_i   <- paste0(S3_STATS_PFX, arr_str, "/")
  
  # Timing (shared across methods, split per-method)
  timing <- tryCatch(
    s3readRDS(object = paste0(pfx_i, "timing.rds"),
              bucket = S3_BUCKET, region = ""),
    error = function(e) c(sim    = NA_real_, nebula = NA_real_,
                          edger  = NA_real_, limma  = NA_real_,
                          wilcox = NA_real_, mast   = NA_real_)
  )
  
  # Scenario parameters to attach
  p_cols <- c("array_id","scenario_id","sim_rep",
              "prop_de","lfc_label","indiv_var_label",
              "cells_label","corr_cells",
              "n_subjects_per_arm",
              "indiv_var_factor","cells_mean","interaction_lfc")
  prow <- param_grid[i, intersect(p_cols, colnames(param_grid))]
  
  # Methods — all five DE methods produced by 02_simulate_analyze.R
  methods <- list(
    nebula = list(file = "nebula_stats.rds",
                  time = as.numeric(timing["nebula.elapsed"])),
    edger  = list(file = "edger_stats.rds",
                  time = as.numeric(timing["edger.elapsed"])),
    limma  = list(file = "limma_stats.rds",
                  time = as.numeric(timing["limma.elapsed"])),
    wilcox = list(file = "wilcox_stats.rds",
                  time = as.numeric(timing["wilcox.elapsed"])),
    mast   = list(file = "mast_stats.rds",
                  time = as.numeric(timing["mast.elapsed"]))
  )
  
  rows <- lapply(names(methods), function(meth) {
    s3_obj <- paste0(pfx_i, methods[[meth]]$file)
    st <- s3_read_safe(s3_obj, S3_BUCKET)
    if (is.null(st)) return(NULL)
    m  <- compute_metrics(st, opt$fdr_thr, meth)
    if (is.null(m)) return(NULL)
    m$method       <- meth
    m$compute_time <- methods[[meth]]$time
    cbind(prow, m)
  })
  
  bind_rows(Filter(Negate(is.null), rows))
}

message("  Processing tasks in parallel...")
raw_list <- bplapply(seq_len(n_scenarios), process_task,
                     BPPARAM = BiocParallel::MulticoreParam(opt$n_cores))

benchmark_raw <- bind_rows(Filter(Negate(is.null), raw_list))
message(sprintf("  Collected %d rows from %d tasks",
                nrow(benchmark_raw), n_scenarios))

# ── Summarise across replicates ───────────────────────────────────────────────
group_vars <- c("method","prop_de","lfc_label","indiv_var_label",
                "cells_label","corr_cells","n_subjects_per_arm")

benchmark_avg <- benchmark_raw %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    n_reps_done      = n(),
    mean_power       = mean(power,        na.rm = TRUE),
    sd_power         = sd(power,          na.rm = TRUE),
    mean_fdr         = mean(fdr_obs,      na.rm = TRUE),
    sd_fdr           = sd(fdr_obs,        na.rm = TRUE),
    mean_t1e         = mean(t1e,          na.rm = TRUE),
    sd_t1e           = sd(t1e,            na.rm = TRUE),
    mean_auc_roc     = mean(auc_roc,      na.rm = TRUE),
    mean_ks_stat     = mean(ks_stat,      na.rm = TRUE),
    sd_ks_stat       = sd(ks_stat,        na.rm = TRUE),
    mean_logfc_bias  = mean(logfc_bias,   na.rm = TRUE),
    sd_logfc_bias    = sd(logfc_bias,     na.rm = TRUE),
    mean_logfc_mae   = mean(logfc_mae,    na.rm = TRUE),
    sd_logfc_mae     = sd(logfc_mae,      na.rm = TRUE),
    mean_compute_s   = mean(compute_time, na.rm = TRUE),
    sd_compute_s     = sd(compute_time,   na.rm = TRUE),
    .groups = "drop"
  )

# ── Save to S3 ────────────────────────────────────────────────────────────────
message("── [03] Saving benchmark outputs to S3 ──")

s3saveRDS(benchmark_raw,
          object = paste0(S3_BENCH_PFX, "benchmark_raw.rds"),
          bucket = S3_BUCKET,
          region = "")
message("  benchmark_raw.rds saved to S3")

s3saveRDS(benchmark_avg,
          object = paste0(S3_BENCH_PFX, "benchmark_avg.rds"),
          bucket = S3_BUCKET,
          region = "")
message("  benchmark_avg.rds saved to S3")

s3write_using_region(
  write.csv,
  benchmark_avg,
  row.names = FALSE,
  object = paste0(S3_BENCH_PFX, "benchmark_avg.csv"),
  bucket = S3_BUCKET,
  region = ""
)
message("  benchmark_avg.csv saved to S3")

message("── [03] Done ──")
message(sprintf("  Unique scenario x method rows: %d", nrow(benchmark_avg)))
print(head(benchmark_avg, 6))