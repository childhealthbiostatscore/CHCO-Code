################################################################################
# 05_benchmark.R
#
# PURPOSE:
#   Aggregate results from all completed simulation runs (output of
#   02_simulate_analyze.R) and compute per-scenario x method summaries:
#     - Power        (sensitivity at BH FDR 5%)
#     - Observed FDR (FP / called significant)
#     - Type I error (FP rate when prop_de == 0)
#     - AUC-ROC      (gene-level p-value as score)
#     - Compute time (wall seconds per method)
#
# INPUT:
#   results/stats/array_NNNNN/{nebula_stats,deseq2_stats,edger_stats,timing,params}.rds
#
# OUTPUT:
#   results/benchmark/
#     benchmark_raw.rds      -- one row per array_id x method (all reps)
#     benchmark_avg.rds      -- averaged/SD across 50 reps per scenario x method
#     benchmark_avg.csv      -- same, human-readable
#
# USAGE:
#   Rscript 05_benchmark.R \
#       --param_grid  results/param_grid/param_grid.rds \
#       --stats_root  results/stats \
#       --out_dir     results/benchmark \
#       --fdr_thr     0.05 \
#       --n_cores     32
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(optparse)
  library(BiocParallel)
  library(ROCR)
})

option_list <- list(
  make_option("--param_grid",  type = "character",
              default = "results/param_grid/param_grid.rds"),
  make_option("--stats_root",  type = "character",
              default = "results/stats"),
  make_option("--out_dir",     type = "character",
              default = "results/benchmark"),
  make_option("--fdr_thr",     type = "double",    default = 0.05),
  make_option("--n_cores",     type = "integer",   default = 32L)
)
opt <- parse_args(OptionParser(option_list = option_list))
register(MulticoreParam(opt$n_cores))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

param_grid  <- readRDS(opt$param_grid)
n_scenarios <- nrow(param_grid)
message(sprintf("── [05] Aggregating %d tasks ──", n_scenarios))

# ── Metric computation ────────────────────────────────────────────────────────
compute_metrics <- function(stats_df, fdr_thr) {
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

  auc_roc <- tryCatch({
    pred_obj <- ROCR::prediction(1 - df$pval_int, df$is_de)
    as.numeric(ROCR::performance(pred_obj, "auc")@y.values)
  }, error = function(e) NA_real_)

  data.frame(
    n_de_true  = n_de,
    n_sig      = n_sig,
    power      = if (n_de    > 0) n_tp / n_de    else NA_real_,
    fdr_obs    = if (n_sig   > 0) n_fp / n_sig   else 0,
    t1e        = if (n_nonde > 0) n_fp / n_nonde else NA_real_,
    auc_roc    = auc_roc,
    stringsAsFactors = FALSE
  )
}

# ── Process one array task ────────────────────────────────────────────────────
process_task <- function(i) {
  arr_id  <- param_grid$array_id[i]
  arr_str <- sprintf("array_%05d", arr_id)
  dir_i   <- file.path(opt$stats_root, arr_str)

  if (!dir.exists(dir_i)) return(NULL)

  # Timing (shared across methods, split per-method)
  timing <- tryCatch(readRDS(file.path(dir_i, "timing.rds")),
                     error = function(e) c(sim=NA, nebula=NA, deseq2=NA, edger=NA))

  # Scenario parameters to attach
  p_cols <- c("array_id","scenario_id","sim_rep",
               "prop_de","effect_size_label","indiv_var_label",
               "cells_label","corr_cells","interaction_label",
               "n_subjects_per_arm",
               "effect_size","indiv_var_factor","cells_mean","interaction_lfc")
  prow <- param_grid[i, p_cols]

  # Methods
  methods <- list(
    nebula = list(file = "nebula_stats.rds",
                  time = as.numeric(timing["nebula"])),
    deseq2 = list(file = "deseq2_stats.rds",
                  time = as.numeric(timing["deseq2"])),
    edger  = list(file = "edger_stats.rds",
                  time = as.numeric(timing["edger"]))
  )

  rows <- lapply(names(methods), function(meth) {
    fpath <- file.path(dir_i, methods[[meth]]$file)
    if (!file.exists(fpath)) return(NULL)
    st <- tryCatch(readRDS(fpath), error = function(e) NULL)
    m  <- compute_metrics(st, opt$fdr_thr)
    if (is.null(m)) return(NULL)
    m$method       <- meth
    m$compute_time <- methods[[meth]]$time
    cbind(prow, m)
  })

  bind_rows(Filter(Negate(is.null), rows))
}

message("  Processing tasks in parallel...")
raw_list <- bplapply(seq_len(n_scenarios), process_task,
                     BPPARAM = MulticoreParam(opt$n_cores))

benchmark_raw <- bind_rows(Filter(Negate(is.null), raw_list))
message(sprintf("  Collected %d rows from %d tasks",
                nrow(benchmark_raw), n_scenarios))

# ── Summarise across replicates ───────────────────────────────────────────────
group_vars <- c("method","prop_de","effect_size_label","indiv_var_label",
                "cells_label","corr_cells","interaction_label","n_subjects_per_arm")

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
    mean_compute_s   = mean(compute_time, na.rm = TRUE),
    sd_compute_s     = sd(compute_time,   na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(benchmark_raw, file.path(opt$out_dir, "benchmark_raw.rds"))
saveRDS(benchmark_avg, file.path(opt$out_dir, "benchmark_avg.rds"))
write.csv(benchmark_avg, file.path(opt$out_dir, "benchmark_avg.csv"),
          row.names = FALSE)

message("── [05] Done ──")
message(sprintf("  Unique scenario x method rows: %d", nrow(benchmark_avg)))
print(head(benchmark_avg, 6))
