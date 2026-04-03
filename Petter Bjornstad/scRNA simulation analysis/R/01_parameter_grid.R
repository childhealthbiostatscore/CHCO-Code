################################################################################
# 01_parameter_grid.R
#
# PURPOSE:
#   Define and save the simulation parameter grid.  Each row is one
#   simulation replicate; SLURM array jobs index into this grid.
#
# DESIGN CHOICES (see README for rationale):
#   - n_subjects_per_arm = subjects per treatment arm (Group A or
#     Group B).  Each subject contributes BOTH a timepoint1 and timepoint2
#     visit, so the full design per arm is:
#         n_subjects_per_arm subjects × 2 visits = 2 × n pseudobulk samples
#     Total pseudobulk samples = 2 arms × n_subjects_per_arm × 2 visits
#     Values: 5, 10, 20  (stress-tests pseudobulk at very low N)
#
# PARAMETERS:
#
#   prop_de              : 0, 0.05, 0.15
#                          (0 = pure null for T1E; 0.05/0.15 = power scenarios)
#
#   lfc_label            : "vlow", "low", "med", "high"
#                          Realistic condition-level log2FC tiers informed by
#                          scRNA-seq DE benchmarking literature (Zimmerman et al.
#                          2021, Murphy & Skene 2022, Gagnon et al.):
#                            vlow = 0.10 log2FC -> 0.069 natural log  (subtle)
#                            low  = 0.25 log2FC -> 0.173 natural log  (Seurat default)
#                            med  = 0.50 log2FC -> 0.347 natural log  (Zimmerman/Gagnon)
#                            high = 1.00 log2FC -> 0.693 natural log  (upper realistic)
#                          These are applied as the interaction coefficient
#                          (visittimepoint2:treatmentTreatment) in the scDesign3
#                          NB GLM, which uses a log link (natural log scale).
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
#   n_subjects_per_arm   : 5, 10, 20
#                          (subjects per treatment arm; each has PRE + POST)
#
#   sim_rep              : 1 ... n_reps  (default 20)
#
# TOTAL:
#   DE scenarios:   2 x 4 x 3 x 2 x 2 x 3 x 20 = 5,760 rows (288 unique)
#   Null scenarios: 1 x 1 x 3 x 2 x 2 x 3 x 20 =   720 rows ( 36 unique)
#   Grand total:                                    6,480 rows (324 unique)
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

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--n_reps", type = "integer", default = 50L,
              help = "Simulation replicates per scenario [default: 50]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── Realistic condition-level effect sizes ─────────────────────────────────────
# log2FC tiers informed by scRNA-seq DE benchmarking literature
# (Zimmerman et al. 2021, Murphy & Skene 2022, Gagnon et al.).
# Most real condition-level effects (disease/treatment within a cell type)
# fall in log2FC 0.1–1.0; larger effects (>1.0) are typical of cell-type
# marker genes, not condition contrasts.
# scDesign3 uses a NB GLM with log link (natural log), so we convert:
#   interaction_lfc = log2FC * log(2)
lfc_tiers <- data.frame(
  lfc_label       = c("vlow", "low", "med", "high"),
  log2FC          = c(0.10,   0.25,  0.50,  1.00),
  interaction_lfc = c(0.10,   0.25,  0.50,  1.00) * log(2),
  stringsAsFactors = FALSE
)
message("Effect size tiers (realistic condition-level log2FCs):")
message(sprintf("  vlow : log2FC=0.10 -> interaction_lfc=%.4f", lfc_tiers$interaction_lfc[1]))
message(sprintf("  low  : log2FC=0.25 -> interaction_lfc=%.4f", lfc_tiers$interaction_lfc[2]))
message(sprintf("  med  : log2FC=0.50 -> interaction_lfc=%.4f", lfc_tiers$interaction_lfc[3]))
message(sprintf("  high : log2FC=1.00 -> interaction_lfc=%.4f", lfc_tiers$interaction_lfc[4]))

# ── Build parameter grid ──────────────────────────────────────────────────────
# Two sub-grids:
#   1. DE scenarios (prop_de > 0): vary lfc_label across vlow/low/med/high
#   2. Null scenarios (prop_de = 0): lfc_label fixed to "vlow" (value is irrelevant
#      since no DE genes are injected; avoids 4× redundant null runs)

shared_vars <- list(
  indiv_var_label    = c("no", "med", "high"),
  cells_label        = c("low", "high"),
  corr_cells         = c(0.1, 0.7),
  n_subjects_per_arm = c(5L, 10L, 20L),
  sim_rep            = seq_len(opt$n_reps)
)

# DE grid: 2 prop_de × 4 lfc × (shared)
grid_de <- do.call(expand.grid, c(
  list(prop_de   = c(0.05, 0.15),
       lfc_label = c("vlow", "low", "med", "high")),
  shared_vars,
  list(stringsAsFactors = FALSE)
))

# Null grid: 1 prop_de × 1 lfc × (shared)
grid_null <- do.call(expand.grid, c(
  list(prop_de   = 0,
       lfc_label = "vlow"),
  shared_vars,
  list(stringsAsFactors = FALSE)
))

param_grid <- bind_rows(grid_de, grid_null)

# ── Map labels to numeric values ──────────────────────────────────────────────
iv_map <- c(no = 1.0, med = 1.5, high = 2.5)

# cells_mean / cells_sd: mean and SD of cells per subject
# (realistic range from ATTEMPT; truncated to minimum 50 in simulation script)
cl_mean_map <- c(low = 500L,  high = 6000L)
cl_sd_map   <- c(low = 150L,  high = 1500L)

# interaction_lfc: merge from lfc_tiers
param_grid <- param_grid %>%
  left_join(lfc_tiers, by = "lfc_label") %>%
  mutate(
    indiv_var_factor = iv_map[indiv_var_label],
    cells_mean       = cl_mean_map[cells_label],
    cells_sd         = cl_sd_map[cells_label],
    
    # Human-readable scenario ID (excludes rep index)
    scenario_id = paste(
      sprintf("pde%.2f",  prop_de),
      paste0("lfc_",      lfc_label),
      paste0("iv_",       indiv_var_label),
      paste0("nc_",       cells_label),
      sprintf("cc%.1f",   corr_cells),
      sprintf("n%02d",    n_subjects_per_arm),
      sep = "_"
    )
  )

# ── Sort: cells_label="low" first, then "high" ───────────────────────────────
# This groups array IDs so you can submit separate SLURM jobs:
#   Low-cell  tasks: array IDs 1 .. n_low
#   High-cell tasks: array IDs (n_low+1) .. n_total
param_grid <- param_grid %>%
  arrange(desc(cells_label == "low")) %>%   # low=TRUE first
  mutate(array_id = seq_len(n()))

n_low  <- sum(param_grid$cells_label == "low")
n_high <- sum(param_grid$cells_label == "high")
message(sprintf("  Array ID layout:  low-cell = 1-%d  |  high-cell = %d-%d",
                n_low, n_low + 1, n_low + n_high))

n_unique <- nrow(distinct(param_grid, scenario_id))
n_total  <- nrow(param_grid)
message(sprintf("Unique scenarios : %d",   n_unique))
message(sprintf("Total runs       : %d  (%d unique x %d reps)",
                n_total, n_unique, opt$n_reps))

# Sanity check:
#   DE:   2 prop_de × 4 lfc × 3 iv × 2 cells × 2 corr × 3 n_subj = 288
#   Null: 1 prop_de × 1 lfc × 3 iv × 2 cells × 2 corr × 3 n_subj =  36
#   Total unique scenarios = 324
stopifnot(n_unique == 324L)

# ── Save to S3 ───────────────────────────────────────────────────────────────
message("── [01] Saving param_grid to S3 ──")

s3saveRDS(param_grid,
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