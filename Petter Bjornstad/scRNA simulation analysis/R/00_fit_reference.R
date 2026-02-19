################################################################################
# 00_fit_reference.R
#
# PURPOSE:
#   Fit a scDesign3 reference model on the top 2,000 HVGs from the ATTEMPT
#   Seurat object (attempt_so).  This script should be run ONCE; the fitted
#   model object is saved to disk and re-used by all subsequent simulation
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
# OUTPUT:
#   results/reference/
#     ├── hvg_genes.rds           -- character vector of top 2k HVG gene names
#     ├── sce_ref.rds             -- SingleCellExperiment used for fitting
#     ├── scdesign3_fit.rds       -- fitted scDesign3 model object
#     └── effect_size_summary.rds -- realistic mean/SD effect sizes from data
#
# USAGE (command line):
#   Rscript 00_fit_reference.R \
#       --seurat_path /path/to/attempt_so.rds \
#       --cell_type "PT"                        \  # KPMP_celltype_general value
#       --out_dir    results/reference           \
#       --n_cores    32
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
})

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--seurat_path", type = "character", default = NULL,
              help = "Path to attempt_so.rds [required]"),
  make_option("--cell_type",   type = "character", default = "PT",
              help = "Value in KPMP_celltype_general to subset [default: PT]"),
  make_option("--out_dir",     type = "character", default = "results/reference",
              help = "Output directory [default: results/reference]"),
  make_option("--n_hvg",       type = "integer",   default = 2000L,
              help = "Number of HVGs to retain [default: 2000]"),
  make_option("--n_cores",     type = "integer",   default = 32L,
              help = "Number of parallel cores [default: 32]"),
  make_option("--seed",        type = "integer",   default = 42L,
              help = "Random seed [default: 42]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$seurat_path)) stop("--seurat_path is required")

set.seed(opt$seed)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

register(MulticoreParam(opt$n_cores))

message("── [00] Loading Seurat object ──")
so <- readRDS(opt$seurat_path)

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
# Use scran-based HVG selection via Seurat FindVariableFeatures
so_sub <- FindVariableFeatures(so_sub, selection.method = "vst",
                               nfeatures = opt$n_hvg, verbose = FALSE)
hvg_genes <- VariableFeatures(so_sub)[seq_len(opt$n_hvg)]
message(sprintf("  Selected %d HVGs", length(hvg_genes)))
saveRDS(hvg_genes, file.path(opt$out_dir, "hvg_genes.rds"))

# ── 3. Build SCE with required colData ───────────────────────────────────────
message("── [00] Building SingleCellExperiment ──")
counts_mat <- GetAssayData(so_sub, slot = "counts")[hvg_genes, ]

col_df <- data.frame(
  cell        = colnames(so_sub),
  subject_id  = so_sub$subject_id,
  visit       = factor(so_sub$visit, levels = c("PRE", "POST")),
  treatment   = factor(so_sub$treatment,
                        levels = c("Placebo", "Dapagliflozin")),
  celltype    = so_sub$KPMP_celltype_general,
  stringsAsFactors = FALSE
)

sce_ref <- SingleCellExperiment(
  assays   = list(counts = counts_mat),
  colData  = col_df
)
saveRDS(sce_ref, file.path(opt$out_dir, "sce_ref.rds"))

# ── 4. Compute realistic effect sizes from data ───────────────────────────────
# We compute log2FC(POST/PRE) per treatment group and summarise the distribution
# so the simulation parameter grid can draw from realistic ranges.
message("── [00] Computing empirical effect size distribution ──")

norm_counts <- log1p(as.matrix(counts_mat))

# Per-gene, per-group mean across visits
compute_group_lfc <- function(grp) {
  idx_pre  <- col_df$treatment == grp & col_df$visit == "PRE"
  idx_post <- col_df$treatment == grp & col_df$visit == "POST"
  mu_pre   <- rowMeans(norm_counts[, idx_pre,  drop = FALSE])
  mu_post  <- rowMeans(norm_counts[, idx_post, drop = FALSE])
  mu_post - mu_pre  # log-scale difference (roughly log2FC)
}

lfc_dapagl  <- compute_group_lfc("Dapagliflozin")
lfc_placebo <- compute_group_lfc("Placebo")
interaction_lfc <- lfc_dapagl - lfc_placebo  # delta-delta

effect_summary <- list(
  lfc_dapagliflozin = summary(lfc_dapagl),
  lfc_placebo       = summary(lfc_placebo),
  interaction_lfc   = summary(interaction_lfc),
  interaction_sd    = sd(interaction_lfc),
  interaction_median_abs = median(abs(interaction_lfc[interaction_lfc != 0])),
  # Suggested simulation levels
  effect_sizes = list(
    no_effect  = 0,
    med_effect = quantile(abs(interaction_lfc), 0.75, na.rm = TRUE),
    high_effect = quantile(abs(interaction_lfc), 0.90, na.rm = TRUE)
  )
)
saveRDS(effect_summary, file.path(opt$out_dir, "effect_size_summary.rds"))
message("  Effect size summary:")
print(effect_summary$effect_sizes)

# ── 5. Fit scDesign3 reference model ─────────────────────────────────────────
# The formula captures:
#   - visit main effect
#   - treatment main effect
#   - visit × treatment interaction  (the quantity of interest)
#   - subject_id as a grouping variable for the random-effect copula
message("── [00] Fitting scDesign3 reference model (this may take a while) ──")

# scDesign3 requires specific colData format
# 'mu_formula' drives the marginal mean; 'sigma_formula' drives the dispersion
# 'corr_formula' drives the copula correlation structure

set.seed(opt$seed)
fit_obj <- scdesign3(
  sce             = sce_ref,
  assay_use       = "counts",
  celltype        = "celltype",        # colData column for cell type label
  pseudotime      = NULL,              # not a trajectory experiment
  spatial         = NULL,
  other_covariates = c("subject_id", "visit", "treatment"),
  mu_formula      = "visit * treatment",
  sigma_formula   = "1",              # global dispersion
  family_use      = "nb",             # negative binomial
  n_cores         = opt$n_cores,
  usebam          = FALSE,
  corr_formula    = "subject_id",      # cell-cell correlation within subject
  copula          = "gaussian",
  DT              = TRUE,
  pseudo_obs      = TRUE,
  return_model    = TRUE,
  parallelization = "pbmcapply",
  BPPARAM         = MulticoreParam(opt$n_cores)
)

saveRDS(fit_obj, file.path(opt$out_dir, "scdesign3_fit.rds"))
message(sprintf("── [00] Done. Reference model saved to %s ──", opt$out_dir))
