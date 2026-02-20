################################################################################
# 01_parameter_grid.R
#
# PURPOSE:
#   Define and save the trimmed simulation parameter grid.  Each row is one
#   simulation replicate; SLURM array jobs index into this grid.
#
# DESIGN CHOICES (see README for rationale):
#   - 432 unique scenarios × 50 reps = 21,600 total runs
#   - n_subjects_per_arm = subjects per treatment arm (Dapagliflozin or
#     Placebo).  Each subject contributes BOTH a PRE and POST visit, so the
#     full design per arm is:
#         n_subjects_per_arm subjects × 2 visits = 2 × n pseudobulk samples
#     Total pseudobulk samples = 2 arms × n_subjects_per_arm × 2 visits
#     Values: 5, 10, 20  (stress-tests pseudobulk at very low N)
#
# PARAMETERS:
#
#   prop_de              : 0, 0.05, 0.15
#                          (0 = pure null for T1E; 0.05/0.15 = power scenarios)
#
#   effect_size_label    : "med", "high"
#                          (derived from ATTEMPT empirical quantiles)
#                          "no" dropped -- redundant with prop_de = 0
#
#   indiv_var_label      : "no", "med", "high"
#                          (sigma inflation factor: 1.0x, 1.5x, 2.5x)
#                          KEY parameter for NEBULA vs pseudobulk comparison
#
#   cells_label          : "low" (~500 cells/subj), "high" (~6000 cells/subj)
#
#   corr_cells           : 0.1, 0.7
#                          (within-subject cell-cell Gaussian copula correlation)
#
#   interaction_label    : "med", "high"
#                          (0.5x or 1.0x of effect_size; "no" dropped -- use
#                           prop_de = 0 for the null instead)
#
#   n_subjects_per_arm   : 5, 10, 20
#                          (subjects per treatment arm; each has PRE + POST)
#
#   sim_rep              : 1 ... n_reps  (default 50)
#
# TOTAL: 3 x 2 x 3 x 2 x 2 x 2 x 3 x 50 = 21,600 rows
#
# INPUT (from S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/reference/
#     effect_size_summary.rds
#
# OUTPUT (to S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/param_grid/
#     param_grid.rds
#     param_grid.csv
#
# USAGE:
#   Rscript 01_parameter_grid.R \
#       --n_reps 50
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
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
S3_REF_PFX   <- "Projects/Paired scRNA simulation analysis/results/reference/"
S3_GRID_PFX  <- "Projects/Paired scRNA simulation analysis/results/param_grid/"

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

s3read_using_region <- function(FUN, ..., object, bucket,
                                 region = NULL, opts = NULL, filename = NULL) {
  ext  <- if (!is.null(filename)) tools::file_ext(filename) else tools::file_ext(object)
  tmp  <- tempfile(fileext = if (nchar(ext) > 0) paste0(".", ext) else "")
  on.exit(unlink(tmp))
  args <- list(object = object, file = tmp, bucket = bucket)
  if (!is.null(region)) args$region <- region
  if (!is.null(opts))   args        <- c(args, opts)
  do.call(aws.s3::save_object, args)
  FUN(tmp, ...)
}

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--n_reps", type = "integer", default = 50L,
              help = "Simulation replicates per scenario [default: 50]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── Load empirical effect sizes from S3 ──────────────────────────────────────
message("── [01] Loading effect size summary from S3 ──")
eff <- s3readRDS(
  object = paste0(S3_REF_PFX, "effect_size_summary.rds"),
  bucket = S3_BUCKET,
  region = ""
)
es_med  <- as.numeric(eff$effect_sizes["med_effect"])
es_high <- as.numeric(eff$effect_sizes["high_effect"])
message(sprintf("Effect sizes from ATTEMPT:  med=%.3f | high=%.3f", es_med, es_high))

# ── Build parameter grid ──────────────────────────────────────────────────────
param_grid <- expand.grid(

  # -- Gene level --------------------------------------------------------------
  prop_de           = c(0, 0.05, 0.15),
  effect_size_label = c("med", "high"),         # "no" dropped

  # -- Individual level --------------------------------------------------------
  indiv_var_label   = c("no", "med", "high"),   # most important for the question
  cells_label       = c("low", "high"),
  corr_cells        = c(0.1, 0.7),
  corr_time         = 0.7,                      # fixed per protocol

  # -- Group/Time level --------------------------------------------------------
  interaction_label = c("med", "high"),         # "no" dropped

  # -- Sample size -------------------------------------------------------------
  # n_subjects_per_arm = subjects per treatment arm (Dapa or Placebo).
  # Each subject contributes PRE + POST visits, so the pseudobulk design has:
  #   2 visits x n_subjects_per_arm samples per arm
  #   = 4 x n_subjects_per_arm total pseudobulk samples
  # e.g. n=5  -> 20 pseudobulk samples total (10 per arm)
  #      n=10 -> 40 pseudobulk samples total
  #      n=20 -> 80 pseudobulk samples total
  n_subjects_per_arm = c(5L, 10L, 20L),

  # -- Replicates --------------------------------------------------------------
  sim_rep           = seq_len(opt$n_reps),

  stringsAsFactors  = FALSE
)

# ── Map labels to numeric values ──────────────────────────────────────────────
es_map <- c(med  = es_med,  high = es_high)
iv_map <- c(no   = 1.0,     med  = 1.5,    high = 2.5)

# cells_mean / cells_sd: mean and SD of cells per subject
# (realistic range from ATTEMPT; truncated to minimum 50 in simulation script)
cl_mean_map <- c(low = 500L,  high = 6000L)
cl_sd_map   <- c(low = 150L,  high = 1500L)

# interaction_multiplier: fraction of effect_size assigned to interaction log-FC
int_map <- c(med = 0.5, high = 1.0)

param_grid <- param_grid %>%
  mutate(
    effect_size            = es_map[effect_size_label],
    indiv_var_factor       = iv_map[indiv_var_label],
    cells_mean             = cl_mean_map[cells_label],
    cells_sd               = cl_sd_map[cells_label],
    interaction_multiplier = int_map[interaction_label],
    interaction_lfc        = interaction_multiplier * effect_size,

    # Human-readable scenario ID (excludes rep index)
    scenario_id = paste(
      sprintf("pde%.2f",  prop_de),
      paste0("es_",       effect_size_label),
      paste0("iv_",       indiv_var_label),
      paste0("nc_",       cells_label),
      sprintf("cc%.1f",   corr_cells),
      paste0("int_",      interaction_label),
      sprintf("n%02d",    n_subjects_per_arm),
      sep = "_"
    ),

    # 1-based linear index for SLURM array
    array_id = seq_len(n())
  )

n_unique <- nrow(distinct(param_grid, scenario_id))
n_total  <- nrow(param_grid)
message(sprintf("Unique scenarios : %d",   n_unique))
message(sprintf("Total runs       : %d  (%d unique x %d reps)",
                n_total, n_unique, opt$n_reps))

# Sanity check: 3 x 2 x 3 x 2 x 2 x 2 x 3 = 432
stopifnot(n_unique == 432L)

# ── Save to S3 ───────────────────────────────────────────────────────────────
message("── [01] Saving param_grid to S3 ──")

s3writeRDS(param_grid,
           object = paste0(S3_GRID_PFX, "param_grid.rds"),
           bucket = S3_BUCKET,
           region = "")
message("  param_grid.rds saved to S3")

s3write_using_region(
  write.csv,
  param_grid,
  row.names = FALSE,
  object = paste0(S3_GRID_PFX, "param_grid.csv"),
  bucket = S3_BUCKET,
  region = ""
)
message("  param_grid.csv saved to S3")

message(sprintf("── [01] Done. Grid saved to s3://%s/%s ──", S3_BUCKET, S3_GRID_PFX))
message("First 3 rows:")
print(head(param_grid, 3))
