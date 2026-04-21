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
#     construct_data.rds      -- list(dat, filtered_gene) from construct_data()
#                                required by extract_para() and simu_new() in step 02
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
  library(scuttle)
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
  make_option("--n_hvg",       type = "integer",   default = 2000,
              help = "Number of HVGs to retain [default: 2000]"),
  make_option("--n_cores",     type = "integer",   default = 32L,
              help = "Number of parallel cores [default: 32]"),
  make_option("--seed",        type = "integer",   default = 42L,
              help = "Random seed [default: 42]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# manual mods for debugging
# opt$cell_type <- "PT"
# opt$n_hvg <- 10

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
#   treatment: Placebo       -> Placebo  (reference group)
#              Dapagliflozin -> Treatment  (treatment group)
#
# Factor levels set so reference = first level:
#   visit:     c("timepoint1", "timepoint2")
#   treatment: c("Placebo", "Treatment")
#
message("── [00] Building SingleCellExperiment ──")
counts_mat <- GetAssayData(so_sub, layer = "counts")[hvg_genes, ]

# Create a clean subject ID mapping
orig_ids <- sort(unique(as.character(so_sub$subject_id)))

# figure out which are Placebo vs Treatment from the original data
orig_trt <- sapply(orig_ids, function(id) {
  unique(as.character(so_sub$treatment[so_sub$subject_id == id]))
})

# assign grpPlacebo_S01, grpTreatment_S01, etc.
placebo_ids   <- orig_ids[orig_trt == "Placebo"]
treatment_ids <- orig_ids[orig_trt != "Placebo"]

clean_ids <- character(length(orig_ids))
names(clean_ids) <- orig_ids

clean_ids[placebo_ids]   <- sprintf("grpPlacebo_S%02d", seq_along(placebo_ids))
clean_ids[treatment_ids] <- sprintf("grpTreatment_S%02d", seq_along(treatment_ids))

col_df <- data.frame(
  cell        = colnames(so_sub),
  subject_id  = factor(clean_ids[as.character(so_sub$subject_id)]),
  visit       = factor(
    ifelse(so_sub$visit == "PRE", "timepoint1", "timepoint2"),
    levels = c("timepoint1", "timepoint2")
  ),
  treatment   = factor(
    ifelse(so_sub$treatment == "Placebo", "Placebo", "Treatment"),
    levels = c("Placebo", "Treatment")
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


sce_ref <- s3readRDS(object = paste0(S3_PREFIX, "sce_ref.rds"),
                     bucket = S3_BUCKET,
                     region = "")

# 1. Pseudobulk
pb <- aggregateAcrossCells(sce_ref,
                           ids = DataFrame(subject = colData(sce_ref)$subject_id,
                                           visit   = colData(sce_ref)$visit))

logcounts <- log1p(sweep(counts(pb), 2, colSums(counts(pb)), "/") * 1e4)

# 2. Variance decomposition using one-way ANOVA per gene
#    ICC = var_between / (var_between + var_within)
genes <- rownames(logcounts)
var_decomp <- do.call(rbind, lapply(genes, function(g) {
  dat <- data.frame(
    y       = logcounts[g, ],
    subject = factor(pb$subject)
  )
  # One-way ANOVA: y ~ subject
  fit <- tryCatch(aov(y ~ subject, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  ms <- summary(fit)[[1]]
  ms_between <- ms["subject",   "Mean Sq"]
  ms_within  <- ms["Residuals", "Mean Sq"]
  n_per_subj <- nrow(dat) / nlevels(dat$subject)  # ~2 (visits)
  
  # Variance components from expected mean squares
  var_within  <- ms_within
  var_between <- max((ms_between - ms_within) / n_per_subj, 0)
  icc         <- var_between / (var_between + var_within)
  
  data.frame(gene = g, var_between = var_between,
             var_within = var_within, icc = icc)
}))

# 3. Summary
cat("ICC summary:\n")
print(summary(var_decomp$icc))
cat("\nBetween-subject variance:\n")
print(summary(var_decomp$var_between))
cat("\nWithin-subject variance:\n")
print(summary(var_decomp$var_within))
cat("\nRatio (between/within):\n")
print(summary(var_decomp$var_between / var_decomp$var_within))

# 4. Quick histogram
hist(var_decomp$icc, breaks = 50, main = "ICC Distribution",
     xlab = "ICC (fraction of variance due to subject)")
abline(v = median(var_decomp$icc), col = "red", lty = 2)

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

lfc_Treatment <- compute_group_lfc("Treatment")
lfc_Placebo <- compute_group_lfc("Placebo")
interaction_lfc <- lfc_Treatment - lfc_Placebo  # delta-delta

effect_summary <- list(
  lfc_Treatment              = summary(lfc_Treatment),
  lfc_Placebo              = summary(lfc_Placebo),
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
#   treatment: reference = "Placebo"
# -> interaction coefficient name = "visittimepoint2:treatmentTreatment"
message("── [00] Fitting scDesign3 reference model (this may take a while) ──")

set.seed(opt$seed)
ref_construct_data <- construct_data(
  sce              = sce_ref,
  assay_use        = "counts",
  celltype         = "celltype",           # colData column for cell type label
  pseudotime       = NULL,                 # not a trajectory experiment
  spatial          = NULL,
  other_covariates = c("subject_id", "visit", "treatment"),
  corr_by = "subject_id")

# fit the NB models per gene (NB and ZINB look similar but NB is more efficient)
ref_fit_marginal <- fit_marginal(
  data = ref_construct_data,
  mu_formula       = "visit * treatment",
  sigma_formula    = "1",                  # global dispersion
  family_use = "nb",
  n_cores = 4,
  usebam = FALSE
)

# fit the copula correlation structure
ref_fit_copula <- fit_copula(
  sce = sce_ref,
  assay_use = "counts",
  input_data = ref_construct_data$dat,
  marginal_list = ref_fit_marginal,
  family_use = "nb",
  copula           = "gaussian",
  DT               = TRUE,
  pseudo_obs       = FALSE,
  important_feature = "all",
  if_sparse = FALSE,
  n_cores = 4,
  parallelization = "mcmapply",
  BPPARAM = NULL
)

# extract parameters from the fitted marginals + copula
ref_extract_para <- extract_para(
  sce              = sce_ref,
  assay_use        = "counts",
  marginal_list = ref_fit_marginal,
  n_cores = 4,
  family_use = "nb",
  new_covariate = ref_construct_data$new_covariate,
  parallelization = "mcmapply",
  BPPARAM = NULL,
  data = ref_construct_data$dat)

# Save only the components needed downstream (not the full SCE again)
ref_data <- list(
  ref_construct_data = ref_construct_data,
  ref_fit_marginal = ref_fit_marginal,
  ref_fit_copula = ref_fit_copula,
  ref_extract_para = ref_extract_para
)

s3saveRDS(ref_data,
          object = paste0(S3_PREFIX, "construct_data.rds"),
          bucket = S3_BUCKET,
          region = "")
message("  construct_data.rds saved to S3")

message(sprintf("── [00] Done. All reference outputs saved to s3://%s/%s ──",
                S3_BUCKET, S3_PREFIX))

# test_res <- simu_new(
#   sce = sce_ref,
#   assay_use = "counts",
#   mean_mat = ref_extract_para$mean_mat,
#   sigma_mat = ref_extract_para$sigma_mat,
#   zero_mat = ref_extract_para$zero_mat,
#   quantile_mat = NULL,
#   copula_list = ref_fit_copula$copula_list,
#   n_cores = 1,
#   family_use = "gaussian",
#   input_data = ref_construct_data$dat,
#   new_covariate = ref_construct_data$new_covariate,
#   important_feature = ref_fit_copula$important_feature,
#   filtered_gene = ref_construct_data$filtered_gene
# )
# 
# 
# ref_cov <- ref_data$ref_construct_data$newCovariate
# cells_per_group <- ref_cov %>%
#   group_by(subject_id, celltype, visit, treatment) %>%
#   summarise(n_cells = n(), .groups = "drop")

################################################################################
# 6. REFERENCE DATA DIAGNOSTICS
#
# Characterize individual-level variance, zero inflation, and cell count
# distributions in the reference data. These diagnostics help interpret
# simulation results — particularly why indiv_var_factor had negligible
# effect on power and why NEBULA vs pseudobulk power gap is smaller in
# simulation than in real data.
################################################################################
message("── [00] Running reference data diagnostics ──")

library(ggplot2)
counts <- assay(sce_ref, "counts")
cd     <- as.data.frame(colData(sce_ref))

# ─────────────────────────────────────────────────────────────────────────────
# 6a. Variance decomposition: between-subject ICC per gene
#     (using pseudobulk log-CPM, separate from the ANOVA above)
# ─────────────────────────────────────────────────────────────────────────────
message("  Computing ICC distribution across genes...")

pb_diag <- aggregate(t(as.matrix(counts)),
                     by = list(subject_id = cd$subject_id, visit = cd$visit),
                     FUN = sum)
pb_mat_diag <- t(as.matrix(pb_diag[, -(1:2)]))
lib_diag    <- colSums(pb_mat_diag)
log_cpm     <- log2(t(t(pb_mat_diag) / lib_diag * 1e6) + 1)

icc_df <- do.call(rbind, lapply(seq_len(nrow(log_cpm)), function(g) {
  y <- log_cpm[g, ]
  s <- pb_diag$subject_id
  total_var <- var(y)
  if (total_var == 0) return(data.frame(gene = rownames(log_cpm)[g],
                                        between = 0, within = 0, icc = 0))
  subj_means  <- tapply(y, s, mean)
  between_var <- var(subj_means)
  within_var  <- mean(tapply(y, s, var), na.rm = TRUE)
  icc <- between_var / (between_var + within_var)
  data.frame(gene = rownames(log_cpm)[g],
             between = between_var, within = within_var, icc = icc)
}))

cat("\n── ICC distribution (pseudobulk log-CPM) ──\n")
cat(sprintf("  Median ICC:           %.3f\n", median(icc_df$icc, na.rm = TRUE)))
cat(sprintf("  Mean ICC:             %.3f\n", mean(icc_df$icc, na.rm = TRUE)))
cat(sprintf("  Genes with ICC > 0.3: %d / %d (%.1f%%)\n",
            sum(icc_df$icc > 0.3, na.rm = TRUE), nrow(icc_df),
            100 * mean(icc_df$icc > 0.3, na.rm = TRUE)))
cat(sprintf("  Genes with ICC > 0.5: %d / %d (%.1f%%)\n",
            sum(icc_df$icc > 0.5, na.rm = TRUE), nrow(icc_df),
            100 * mean(icc_df$icc > 0.5, na.rm = TRUE)))

# ─────────────────────────────────────────────────────────────────────────────
# 6b. scDesign3 fitted dispersion (sigma) — what does inflating it do?
# ─────────────────────────────────────────────────────────────────────────────
message("  Extracting fitted dispersion parameters...")

sigmas <- sapply(ref_data$ref_fit_marginal, function(gm) {
  if (!is.null(gm$fit)) {
    log_theta <- gm$fit$family$getTheta()
    if (!is.null(log_theta) && is.finite(log_theta)) {
      theta <- exp(log_theta)
      return(1 / theta)  # sigma = 1/theta
    }
  }
  NA
})
sigmas <- sigmas[!is.na(sigmas)]

cat("\n── Fitted NB dispersion (sigma) ──\n")
cat(sprintf("  N genes:  %d\n", length(sigmas)))
cat(sprintf("  Median:   %.4f\n", median(sigmas)))
cat(sprintf("  Mean:     %.4f\n", mean(sigmas)))
cat(sprintf("  SD:       %.4f\n", sd(sigmas)))
cat(sprintf("  Range:    [%.4f, %.4f]\n", min(sigmas), max(sigmas)))

# NB variance = mu + sigma * mu^2
# So sigma controls how much variance exceeds Poisson.
# If sigma is already small, multiplying by 1.5x or 2.5x barely changes counts.
# Show what the NB CV looks like at a typical mean expression level:
typical_mu <- 1.0  # ~ median expression for HVGs in raw counts
for (fct in c(1.0, 1.5, 2.5)) {
  nb_var <- typical_mu + fct * median(sigmas) * typical_mu^2
  nb_cv  <- sqrt(nb_var) / typical_mu
  cat(sprintf("  At mu=%.1f, sigma_factor=%.1fx: NB_var=%.3f, CV=%.3f\n",
              typical_mu, fct, nb_var, nb_cv))
}

# ─────────────────────────────────────────────────────────────────────────────
# 6c. Zero inflation: observed vs NB-predicted zero rates
#
# P(Y=0 | NB) = (theta / (theta + mu))^theta where theta = 1/sigma
# If observed >> predicted, the NB model underestimates sparsity.
# This artificially helps pseudobulk in simulation (cleaner aggregates)
# while in real data, sparsity hurts pseudobulk more than NEBULA.
# ─────────────────────────────────────────────────────────────────────────────
message("  Computing observed vs NB-predicted zero rates...")

obs_zero_rate <- rowMeans(counts == 0)

pred_zero_rate <- sapply(seq_along(ref_data$ref_fit_marginal), function(g) {
  gm <- ref_data$ref_fit_marginal[[g]]
  if (is.null(gm$fit)) return(NA)
  log_theta <- gm$fit$family$getTheta()
  if (is.null(log_theta) || !is.finite(log_theta)) return(NA)
  mu    <- mean(fitted(gm$fit))
  theta <- exp(log_theta)
  (theta / (theta + mu))^theta
})
names(pred_zero_rate) <- names(ref_data$ref_fit_marginal)

zero_df <- data.frame(
  gene         = rownames(counts),
  observed     = obs_zero_rate,
  predicted_NB = pred_zero_rate[rownames(counts)]
) %>% filter(!is.na(predicted_NB))

excess <- zero_df$observed - zero_df$predicted_NB

cat("\n── Zero inflation diagnostic ──\n")
cat(sprintf("  Median observed zero rate:     %.3f\n", median(zero_df$observed)))
cat(sprintf("  Median NB-predicted zero rate: %.3f\n", median(zero_df$predicted_NB)))
cat(sprintf("  Median excess zeros:           %.3f\n", median(excess)))
cat(sprintf("  Genes with >10%% excess zeros: %d / %d (%.1f%%)\n",
            sum(excess > 0.10), nrow(zero_df),
            100 * mean(excess > 0.10)))
cat(sprintf("  Genes with >20%% excess zeros: %d / %d (%.1f%%)\n",
            sum(excess > 0.20), nrow(zero_df),
            100 * mean(excess > 0.20)))

# ─────────────────────────────────────────────────────────────────────────────
# 6d. Cell count distribution per subject x visit
#
# If highly uneven in real data but symmetric in simulation (clipped normal),
# pseudobulk gets artificially uniform precision in simulation.
# ─────────────────────────────────────────────────────────────────────────────
message("  Summarizing cell counts per subject x visit...")

cell_counts <- cd %>%
  group_by(subject_id, visit, treatment) %>%
  summarise(n_cells = n(), .groups = "drop")

cat("\n── Cell counts per subject x visit ──\n")
cat(sprintf("  Min:    %d\n",   min(cell_counts$n_cells)))
cat(sprintf("  Median: %d\n",   median(cell_counts$n_cells)))
cat(sprintf("  Mean:   %.0f\n", mean(cell_counts$n_cells)))
cat(sprintf("  Max:    %d\n",   max(cell_counts$n_cells)))
cat(sprintf("  SD:     %.0f\n", sd(cell_counts$n_cells)))
cat(sprintf("  CV:     %.2f\n", sd(cell_counts$n_cells) / mean(cell_counts$n_cells)))
print(as.data.frame(cell_counts))

# ─────────────────────────────────────────────────────────────────────────────
# 6e. Save diagnostics and plots
# ─────────────────────────────────────────────────────────────────────────────
diag_output <- list(
  icc_df      = icc_df,
  sigmas      = sigmas,
  zero_df     = zero_df,
  cell_counts = cell_counts
)
s3saveRDS(diag_output,
          object = paste0(S3_PREFIX, "ref_diagnostics.rds"),
          bucket = S3_BUCKET, region = "")
message("  ref_diagnostics.rds saved to S3")

# plots
pdf_tmp <- tempfile(fileext = ".pdf")
pdf(pdf_tmp, width = 10, height = 8)

# ICC histogram
hist(icc_df$icc, breaks = 50, col = "steelblue",
     main = "ICC Distribution (between-subject / total variance)",
     xlab = "ICC", ylab = "Number of genes")
abline(v = median(icc_df$icc, na.rm = TRUE), col = "red", lty = 2)

# sigma distribution with inflation overlay
sigma_all <- data.frame(
  sigma  = c(sigmas, sigmas * 1.5, sigmas * 2.5),
  factor = rep(c("1.0x (ref)", "1.5x", "2.5x"), each = length(sigmas))
)
print(ggplot(sigma_all, aes(x = sigma, fill = factor)) +
        geom_density(alpha = 0.4) +
        labs(title = "NB dispersion: reference vs inflated",
             x = "sigma", y = "Density") +
        theme_minimal())

# observed vs predicted zeros
print(ggplot(zero_df, aes(x = predicted_NB, y = observed)) +
        geom_point(alpha = 0.3, size = 0.8) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        labs(title = "Observed vs NB-predicted zero rates",
             subtitle = "Above diagonal = excess zeros not captured by NB",
             x = "Predicted P(Y=0) under NB", y = "Observed zero rate") +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        theme_minimal())

# cell counts bar chart
print(ggplot(cell_counts, aes(x = interaction(subject_id, visit),
                              y = n_cells, fill = treatment)) +
        geom_col() +
        labs(title = "Cell counts per subject x visit (reference)",
             x = NULL, y = "Number of cells") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)))

dev.off()
put_object(file = pdf_tmp,
           object = paste0(S3_PREFIX, "ref_diagnostics.pdf"),
           bucket = S3_BUCKET, region = "")
unlink(pdf_tmp)

message("  ref_diagnostics.pdf saved to S3")
message("── [00] Diagnostics complete ──")