################################################################################
# 00_fit_reference.R
#
# PURPOSE:
#   Fit a scDesign3 reference model on the top 2,000 HVGs from the ATTEMPT
#   Seurat object (attempt_so).  This script should be run ONCE; the fitted
#   model object is saved to S3 and re-used by all subsequent simulation
#   scripts.
#
# DESIGN:
#   - Paired study: participants have a PRE and POST visit
#   - Two treatment groups: Dapagliflozin (Group A) and Placebo (Group B)
#   - Cell type: user specifies CELL_TYPE via --cell_type argument (or edit
#     CELL_TYPE below).  Use KPMP_celltype_general for broad categories or
#     KPMP_celltype for granular ones.
#   - Subject-level random effects are captured via scDesign3's sigma/corr
#     components.
#
# INPUT (from S3):
#   bucket: attempt
#     cleaned_data/attempt_clean_so.rds   -- ATTEMPT Seurat object
#
# OUTPUT (to S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/reference/
#     hvg_genes.rds           -- character vector of top 2k HVG gene names
#     sce_ref.rds             -- SingleCellExperiment used for fitting
#     scdesign3_fit.rds       -- fitted scDesign3 model object
#     effect_size_summary.rds -- realistic mean/SD effect sizes from data
#
# USAGE (command line):
#   Rscript 00_fit_reference.R \
#       --cell_type "PT"   \  # KPMP_celltype_general value
#       --n_cores   32
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDesign3)
  library(dplyr)
  library(Matrix)
  library(BiocParallel)
  library(optparse)
  library(tibble)
  library(aws.s3)
})

# ── S3 / Multi-user setup ─────────────────────────────────────────────────────
# Detect user and set up AWS credentials + endpoint
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

# S3 output prefix
S3_BUCKET <- "scrna"
S3_PREFIX <- "Projects/Paired scRNA simulation analysis/results/reference/"

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
  make_option("--cell_type",   type = "character", default = "PT",
              help = "Value in KPMP_celltype_general to subset [default: PT]"),
  make_option("--n_hvg",       type = "integer",   default = 2000L,
              help = "Number of HVGs to retain [default: 2000]"),
  make_option("--n_cores",     type = "integer",   default = 32L,
              help = "Number of parallel cores [default: 32]"),
  make_option("--seed",        type = "integer",   default = 42L,
              help = "Random seed [default: 42]")
)
opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)
BiocParallel::register(BiocParallel::MulticoreParam(opt$n_cores))

# ── Load attempt_so from S3 ───────────────────────────────────────────────────
message("── [00] Loading Seurat object from S3 (bucket: attempt) ──")
so <- s3readRDS(
  object = "cleaned_data/attempt_clean_so.rds",
  bucket = "attempt",
  region = ""
)
message(sprintf("  Loaded: %d genes x %d cells", nrow(so), ncol(so)))

celltype_groups <- list(
  PT              = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL             = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC              = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  EC              = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC"),
  IC              = c("IC-A", "IC-B", "aIC"),
  Immune          = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK"),
  # Immune_myeloid  = c("MAC", "MON", "cDC", "pDC"),
  # Immune_lymphoid = c("CD4+ T", "CD8+ T", "B", "NK"),
  VSMC_P_FIB      = c("VSMC/P", "FIB"),
  POD             = "POD"
)

lookup <- unlist(lapply(names(celltype_groups), function(group) {
  setNames(rep(group, length(celltype_groups[[group]])),
           celltype_groups[[group]])
}))

so$KPMP_celltype_general <- ifelse(
  so$KPMP_celltype %in% names(lookup),
  lookup[so$KPMP_celltype],
  so$KPMP_celltype
)

# ── 1. Subset to cell type ────────────────────────────────────────────────────
message(sprintf("── [00] Subsetting to cell type: %s ──", opt$cell_type))
cells_keep <- which(so$KPMP_celltype_general == opt$cell_type)
if (length(cells_keep) == 0) {
  stop(sprintf("No cells found with KPMP_celltype_general == '%s'. ",
               opt$cell_type),
       "Available values: ", paste(unique(so$KPMP_celltype_general), collapse = ", "))
}
so_sub <- so[, cells_keep]
rm(so); gc()

# ── 2. Select top 2k HVGs ─────────────────────────────────────────────────────
message("── [00] Identifying top HVGs ──")
so_sub <- FindVariableFeatures(so_sub, selection.method = "vst",
                               nfeatures = opt$n_hvg, verbose = FALSE)
hvg_genes <- VariableFeatures(so_sub)[seq_len(opt$n_hvg)]
message(sprintf("  Selected %d HVGs", length(hvg_genes)))

s3saveRDS(hvg_genes,
          object = paste0(S3_PREFIX, "hvg_genes.rds"),
          bucket = S3_BUCKET,
          region = "")
message("  hvg_genes.rds saved to S3")

# ── 3. Build SCE with required colData ───────────────────────────────────────
# Recode clinical labels to generic simulation labels:
#   visit:     PRE           -> timepoint1
#              POST          -> timepoint2
#   treatment: Placebo       -> groupB  (reference group)
#              Dapagliflozin -> groupA  (treatment group)
#
# Factor levels set so reference = first level:
#   visit:     c("timepoint1", "timepoint2")
#   treatment: c("groupB", "groupA")
#
message("── [00] Building SingleCellExperiment ──")
counts_mat <- GetAssayData(so_sub, slot = "counts")[hvg_genes, ]

col_df <- data.frame(
  cell        = colnames(so_sub),
  subject_id  = factor(as.character(so_sub$subject_id)),
  visit       = factor(
    ifelse(so_sub$visit == "PRE", "timepoint1", "timepoint2"),
    levels = c("timepoint1", "timepoint2")
  ),
  treatment   = factor(
    ifelse(so_sub$treatment == "Placebo", "groupB", "groupA"),
    levels = c("groupB", "groupA")
  ),
  celltype    = so_sub$KPMP_celltype_general,
  stringsAsFactors = FALSE
)

sce_ref <- SingleCellExperiment(
  assays   = list(counts = counts_mat),
  colData  = col_df
)

s3saveRDS(sce_ref,
          object = paste0(S3_PREFIX, "sce_ref.rds"),
          bucket = S3_BUCKET,
          region = "")
message("  sce_ref.rds saved to S3")

# ── 4. Compute realistic effect sizes from data ───────────────────────────────
# We compute log2FC(timepoint2/timepoint1) per treatment group and summarise
# the distribution so the simulation parameter grid can draw from realistic ranges.
message("── [00] Computing empirical effect size distribution ──")

norm_counts <- log1p(as.matrix(counts_mat))

# Per-gene, per-group mean across visits
compute_group_lfc <- function(grp) {
  idx_pre  <- col_df$treatment == grp & col_df$visit == "timepoint1"
  idx_post <- col_df$treatment == grp & col_df$visit == "timepoint2"
  mu_pre   <- rowMeans(norm_counts[, idx_pre,  drop = FALSE])
  mu_post  <- rowMeans(norm_counts[, idx_post, drop = FALSE])
  mu_post - mu_pre  # log-scale difference (roughly log2FC)
}

lfc_groupA     <- compute_group_lfc("groupA")
lfc_groupB     <- compute_group_lfc("groupB")
interaction_lfc <- lfc_groupA - lfc_groupB  # delta-delta

effect_summary <- list(
  lfc_groupA              = summary(lfc_groupA),
  lfc_groupB              = summary(lfc_groupB),
  interaction_lfc         = summary(interaction_lfc),
  interaction_sd          = sd(interaction_lfc),
  interaction_median_abs  = median(abs(interaction_lfc[interaction_lfc != 0])),
  # Suggested simulation levels
  effect_sizes = list(
    no_effect   = 0,
    med_effect  = quantile(abs(interaction_lfc), 0.75, na.rm = TRUE),
    high_effect = quantile(abs(interaction_lfc), 0.90, na.rm = TRUE)
  )
)

s3saveRDS(effect_summary,
          object = paste0(S3_PREFIX, "effect_size_summary.rds"),
          bucket = S3_BUCKET,
          region = "")
message("  effect_size_summary.rds saved to S3")
message("  Effect size summary:")
print(effect_summary$effect_sizes)

# ── 5. Fit scDesign3 reference model ─────────────────────────────────────────
# The formula captures:
#   - visit main effect
#   - treatment main effect
#   - visit × treatment interaction  (the quantity of interest)
#   - subject_id as a grouping variable for the random-effect copula
#
# Factor levels:
#   visit:     reference = "timepoint1"
#   treatment: reference = "groupB"
# -> interaction coefficient name = "visittimepoint2:treatmentgroupA"
message("── [00] Fitting scDesign3 reference model (this may take a while) ──")

set.seed(opt$seed)
fit_obj <- scdesign3(
  sce              = sce_ref,
  assay_use        = "counts",
  celltype         = "celltype",           # colData column for cell type label
  pseudotime       = NULL,                 # not a trajectory experiment
  spatial          = NULL,
  other_covariates = c("subject_id", "visit", "treatment"),
  mu_formula       = "visit * treatment",
  sigma_formula    = "1",                  # global dispersion
  family_use       = "nb",                 # negative binomial
  n_cores          = opt$n_cores,
  usebam           = FALSE,
  corr_formula     = "subject_id",         # cell-cell correlation within subject
  copula           = "gaussian",
  DT               = TRUE,
  pseudo_obs       = TRUE,
  return_model     = TRUE,
  parallelization  = "pbmcapply",
  BPPARAM          = BiocParallel::MulticoreParam(opt$n_cores)
)

s3saveRDS(fit_obj,
          object = paste0(S3_PREFIX, "scdesign3_fit.rds"),
          bucket = S3_BUCKET,
          region = "")
message("  scdesign3_fit.rds saved to S3")

# ── 6. Save construct_data output for use in 02_simulate_analyze.R ────────────
# scdesign3() does not expose the construct_data() output directly.
# We re-run construct_data() to obtain $dat and $filtered_gene, which are
# required by extract_para() and simu_new() in step 02.
message("── [00] Running construct_data() to save dat + filtered_gene ──")

# Build the new_covariate (= covariate of the reference data itself)
ref_covariate <- as.data.frame(colData(sce_ref)[, c("subject_id", "visit",
                                                    "treatment", "celltype")])

construct_data_result <- construct_data(
  sce              = sce_ref,
  assay_use        = "counts",
  celltype         = "celltype",
  pseudotime       = NULL,
  spatial          = NULL,
  other_covariates = c("subject_id", "visit", "treatment"),
  corr_by          = "subject_id"
)

# Save only the components needed downstream (not the full SCE again)
ref_data <- list(
  dat           = construct_data_result$dat,
  filtered_gene = construct_data_result$filtered_gene
)

s3saveRDS(ref_data,
          object = paste0(S3_PREFIX, "construct_data.rds"),
          bucket = S3_BUCKET,
          region = "")
message("  construct_data.rds saved to S3")

message(sprintf("── [00] Done. All reference outputs saved to s3://%s/%s ──",
                S3_BUCKET, S3_PREFIX))