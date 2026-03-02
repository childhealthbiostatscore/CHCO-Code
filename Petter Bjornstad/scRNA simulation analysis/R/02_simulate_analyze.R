################################################################################
# 02_simulate_analyze.R
#
# PURPOSE:
#   COMBINED script: simulate one scRNA-seq dataset, run all three DE methods
#   (NEBULA-LN, DESeq2 pseudobulk, edgeR pseudobulk), save only the small
#   stats output files to S3, and discard the raw count matrix.
#
#   This avoids writing 5-15 MB sparse matrices for each of 21,600 tasks,
#   reducing total storage from ~300-800 GB down to ~1-2 GB.
#
# DESIGN:
#   Each SLURM array task corresponds to one row in param_grid.rds.
#   The task:
#     1. Draws the simulated design (subjects x visits x cells)
#     2. Modifies the scDesign3 fitted model to inject DE signal + indiv variance
#     3. Calls extract_para() + simu_new() -- count matrix lives only in memory
#     4. Immediately runs NEBULA, DESeq2, edgeR, Wilcoxon, MAST on the in-memory counts
#     5. Saves three small stats data.frames + timing + truth to S3
#     6. Frees the count matrix (gc())
#
# STUDY DESIGN CLARIFICATION:
#   n_subjects_per_arm = number of SUBJECTS in each treatment arm.
#   Each subject has BOTH a timepoint1 and a timepoint2 visit (paired design).
#   So the full cell-level design is:
#     2 arms x n_subjects_per_arm subjects x 2 visits x ~cells_mean cells
#   And the pseudobulk design is:
#     2 arms x n_subjects_per_arm subjects x 2 visits
#     = 4 x n_subjects_per_arm pseudobulk samples total
#   Example: n_subjects_per_arm=5  -> 20 pseudobulk samples (10 per arm)
#            n_subjects_per_arm=10 -> 40 pseudobulk samples
#            n_subjects_per_arm=20 -> 80 pseudobulk samples
#
# INPUT (from S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/
#     param_grid/param_grid.rds
#     reference/scdesign3_fit.rds
#     reference/hvg_genes.rds
#     reference/sce_ref.rds
#
# OUTPUT (to S3, per array task, ~1-2 MB total):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/stats/array_NNNNN/
#     nebula_stats.rds    -- gene x {logFC, pval, padj, is_de, beta_int_true}
#     deseq2_stats.rds    -- same structure
#     edger_stats.rds     -- same structure
#     wilcox_stats.rds    -- same structure (contrast: groupA POST vs groupB POST)
#     mast_stats.rds      -- same structure (MAST + subject_id latent covariate)
#     timing.rds          -- named vector: sim, nebula, deseq2, edger, wilcox, mast (s)
#     params.rds          -- single-row data.frame of this task's parameters
#     truth.rds           -- data.frame: gene, is_de, beta_int_true
#
# NOTE on Wilcoxon / MAST contrast:
#   These methods compare groupA POST cells vs groupB POST cells, which tests
#   (b_trt + b_int), not the pure interaction b_int. Under prop_de = 0 (b_int=0),
#   detections reflect the main treatment effect from the reference model,
#   demonstrating type I error inflation from ignoring paired structure.
#
# USAGE:
#   Rscript 02_simulate_analyze.R \
#       --array_id      1 \
#       --n_cores       4 \
#       --mast_max_cells 5000
################################################################################

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(nebula)
  library(DESeq2)
  library(edgeR)
  library(Matrix)
  library(dplyr)
  library(BiocParallel)
  library(optparse)
  library(aws.s3)
  library(Seurat)
  library(MAST)
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
S3_BUCKET   <- "scrna"
S3_BASE     <- "Projects/Paired scRNA simulation analysis/results/"
S3_REF_PFX  <- paste0(S3_BASE, "reference/")
S3_GRID_PFX <- paste0(S3_BASE, "param_grid/")

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--array_id",      type = "integer",   default = NULL,
              help = "Row index in param_grid (1-based) [required]"),
  make_option("--n_cores",       type = "integer",   default = 4L),
  make_option("--nebula_method",   type = "character", default = "LN",
              help = "NEBULA method: LN or HL"),
  make_option("--pb_min_count",    type = "integer",   default = 10L,
              help = "Min rowSum for pseudobulk filtering"),
  make_option("--pb_min_samples",  type = "integer",   default = 2L,
              help = "Min samples meeting pb_min_count threshold"),
  make_option("--mast_max_cells",  type = "integer",   default = 5000L,
              help = "Max POST cells to pass to MAST (subsample if exceeded) [default: 5000]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$array_id)) stop("--array_id is required")

# Seed is array_id so each task is reproducible and unique
set.seed(opt$array_id)
message(sprintf("  Seed set to: %d (= array_id)", opt$array_id))
BiocParallel::register(BiocParallel::MulticoreParam(opt$n_cores))

# Per-task S3 output prefix
arr_str <- sprintf("array_%05d", opt$array_id)
S3_OUT  <- paste0(S3_BASE, "stats/", arr_str, "/")

message(sprintf("══ [02] array_id=%d  pid=%d ══", opt$array_id, Sys.getpid()))

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD INPUTS FROM S3
# ─────────────────────────────────────────────────────────────────────────────
message("  Loading param_grid from S3...")
param_grid <- s3readRDS(
  object = paste0(S3_GRID_PFX, "param_grid.rds"),
  bucket = S3_BUCKET,
  region = ""
)
params <- param_grid[opt$array_id, ]

message("  Loading reference model from S3...")
fit_obj      <- s3readRDS(object = paste0(S3_REF_PFX, "scdesign3_fit.rds"),
                          bucket = S3_BUCKET, region = "")
hvg_genes    <- s3readRDS(object = paste0(S3_REF_PFX, "hvg_genes.rds"),
                          bucket = S3_BUCKET, region = "")
sce_ref      <- s3readRDS(object = paste0(S3_REF_PFX, "sce_ref.rds"),
                          bucket = S3_BUCKET, region = "")
ref_data     <- s3readRDS(object = paste0(S3_REF_PFX, "construct_data.rds"),
                          bucket = S3_BUCKET, region = "")
# ref_data$dat           -- input data.frame for extract_para() and simu_new()
# ref_data$filtered_gene -- genes filtered in construct_data(), for simu_new()

n_genes <- length(hvg_genes)  # 2000

message(sprintf("  Scenario: %s  rep=%d", params$scenario_id, params$sim_rep))

# ─────────────────────────────────────────────────────────────────────────────
# 2. BUILD CELL-LEVEL DESIGN
#
# n_subjects_per_arm subjects in each of 2 arms (groupA / groupB),
# each with timepoint1 + timepoint2 visit (paired design).
# Cells per subject drawn from N(cells_mean, cells_sd^2), clipped to >= 50.
# ─────────────────────────────────────────────────────────────────────────────
treatments <- c("groupB", "groupA")
visits     <- c("timepoint1", "timepoint2")

make_subject_cells <- function(subj_id, trt, cells_mean, cells_sd) {
  # Each subject has cells from BOTH visits; split roughly 50/50 between visits
  n_total <- max(50L, round(rnorm(1, cells_mean, cells_sd)))
  # Assign visits proportionally (slightly random split)
  n_pre   <- round(n_total * runif(1, 0.4, 0.6))
  n_post  <- n_total - n_pre
  data.frame(
    subject_id = rep(subj_id, n_total),
    treatment  = rep(trt,     n_total),
    visit      = c(rep("timepoint1", n_pre), rep("timepoint2", n_post)),
    stringsAsFactors = FALSE
  )
}

cell_list <- lapply(treatments, function(trt) {
  lapply(seq_len(params$n_subjects_per_arm), function(j) {
    # Subject ID encodes arm and subject number
    subj <- sprintf("%s_S%02d",
                    ifelse(trt == "groupA", "grpA", "grpB"), j)
    make_subject_cells(subj, trt, params$cells_mean, params$cells_sd)
  })
})
new_covariate <- do.call(rbind, do.call(c, cell_list))
new_covariate$visit     <- factor(new_covariate$visit,
                                  levels = c("timepoint1", "timepoint2"))
new_covariate$treatment <- factor(new_covariate$treatment,
                                  levels = c("groupB", "groupA"))
new_covariate$celltype  <- unique(sce_ref$celltype)[1]

n_cells_total <- nrow(new_covariate)
n_subjects    <- length(unique(new_covariate$subject_id))
message(sprintf("  Cells: %d  |  Subjects: %d (%d per arm x 2 arms)",
                n_cells_total, n_subjects, params$n_subjects_per_arm))

# ─────────────────────────────────────────────────────────────────────────────
# 3. GROUND TRUTH: ASSIGN DE GENES AND TRUE INTERACTION EFFECTS
# ─────────────────────────────────────────────────────────────────────────────
n_de   <- round(params$prop_de * n_genes)
de_idx <- if (n_de > 0) sample(seq_len(n_genes), n_de) else integer(0)

beta_int_true <- numeric(n_genes)
if (n_de > 0) {
  signs <- sample(c(-1L, 1L), n_de, replace = TRUE)
  beta_int_true[de_idx] <- signs * params$interaction_lfc
}

truth_df <- data.frame(
  gene          = hvg_genes,
  is_de         = seq_len(n_genes) %in% de_idx,
  beta_int_true = beta_int_true,
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────────────────
# 4. MODIFY scDesign3 MODEL: inject DE signal + individual-level variance
#
# Factor levels in sce_ref colData are:
#   visit:     c("timepoint1", "timepoint2")  -> reference = "timepoint1"
#   treatment: c("groupB", "groupA")          -> reference = "groupB"
# So the interaction coefficient is: "visittimepoint2:treatmentgroupA"
# ─────────────────────────────────────────────────────────────────────────────
COEF_INT <- "visittimepoint2:treatmentgroupA"

modify_gene_model <- function(gm, beta_int, iv_factor) {
  if (!is.null(gm$fit)) {
    cf <- coef(gm$fit)
    if (COEF_INT %in% names(cf)) {
      cf[COEF_INT] <- beta_int
    } else {
      cf <- c(cf, setNames(beta_int, COEF_INT))
    }
    gm$fit$coefficients <- cf
  }
  if (!is.null(gm$sigma) && iv_factor > 1.0) {
    gm$sigma <- gm$sigma * iv_factor
  }
  gm
}

message("  Modifying scDesign3 model (DE injection + individual variance)...")
modified_marginal <- bplapply(
  seq_along(fit_obj$marginal_list),
  function(g) modify_gene_model(fit_obj$marginal_list[[g]],
                                beta_int  = beta_int_true[g],
                                iv_factor = params$indiv_var_factor),
  BPPARAM = BiocParallel::MulticoreParam(opt$n_cores)
)
names(modified_marginal) <- names(fit_obj$marginal_list)

# copula_list is set to NULL: simulate gene marginals independently.
# The fitted copula from 00_fit_reference.R is keyed on reference subject IDs
# and cannot be directly applied to new simulated subjects. For benchmarking
# DE methods, marginal accuracy (per-gene mean/dispersion) matters more than
# between-gene correlation. Between-gene correlation structure is dropped here.
fit_mod <- fit_obj
fit_mod$marginal_list <- modified_marginal

# ─────────────────────────────────────────────────────────────────────────────
# 5. SIMULATE COUNTS (held in memory only)
#
# scDesign3 modular API:
#   extract_para()  -- converts marginal_list -> parameter matrices
#   simu_new()      -- uses parameter matrices + copula to generate counts
# ─────────────────────────────────────────────────────────────────────────────
t_sim <- proc.time()
message("  Extracting parameters from modified model...")

sim_para <- tryCatch(
  extract_para(
    sce           = sce_ref,
    marginal_list = fit_mod$marginal_list,
    n_cores       = opt$n_cores,
    family_use    = "nb",
    new_covariate = new_covariate,
    data          = ref_data$dat
  ),
  error = function(e) {
    message("  ERROR in extract_para: ", conditionMessage(e))
    NULL
  }
)

if (is.null(sim_para)) {
  s3saveRDS(list(error    = "extract_para failed",
                 params   = params,
                 array_id = opt$array_id),
            object = paste0(S3_OUT, "error.rds"),
            bucket = S3_BUCKET,
            region = "")
  quit(status = 1)
}

message("  Simulating counts with simu_new()...")
# simu_new() in scDesign3 v1.x does not handle copula_list = NULL gracefully —
# it enters the copula branch regardless and crashes on an empty list.
# The workaround: supply quantile_mat directly (independent U(0,1) draws per
# gene per cell). simu_new() skips copula sampling when quantile_mat is provided
# and uses these values to invert the fitted NB marginals instead.
quantile_mat_ind <- matrix(
  runif(prod(dim(sim_para$mean_mat))),
  nrow     = nrow(sim_para$mean_mat),
  ncol     = ncol(sim_para$mean_mat),
  dimnames = dimnames(sim_para$mean_mat)
)

sim_result <- tryCatch(
  simu_new(
    sce               = sce_ref,
    mean_mat          = sim_para$mean_mat,
    sigma_mat         = sim_para$sigma_mat,
    zero_mat          = sim_para$zero_mat,
    quantile_mat      = quantile_mat_ind,  # bypasses copula; genes are independent
    copula_list       = NULL,
    n_cores           = opt$n_cores,
    family_use        = "nb",
    input_data        = ref_data$dat,
    new_covariate     = new_covariate,
    important_feature = fit_obj$important_feature,
    filtered_gene     = ref_data$filtered_gene
  ),
  error = function(e) {
    message("  ERROR in simu_new: ", conditionMessage(e))
    NULL
  }
)
rm(quantile_mat_ind); gc(verbose = FALSE)

t_sim_elapsed <- (proc.time() - t_sim)["elapsed"]

if (is.null(sim_result)) {
  s3saveRDS(list(error    = "simu_new failed",
                 params   = params,
                 array_id = opt$array_id),
            object = paste0(S3_OUT, "error.rds"),
            bucket = S3_BUCKET,
            region = "")
  quit(status = 1)
}

counts <- sim_result               # simu_new() returns dgCMatrix directly in this build
if (is.null(rownames(counts))) rownames(counts) <- hvg_genes

# ── Validate simulated counts ──────────────────────────────────────────────
# qnbinom() returns NA when mu ~ 0 and the uniform draw is near 1 (extreme
# tail of a near-degenerate NB). Inf can appear with tiny size + large mu.
# Both poison every downstream DE method's likelihood. Replace with 0 (treat
# as structural dropout — equivalent to the cell not expressing that gene).
if (inherits(counts, "dgCMatrix")) {
  bad <- !is.finite(counts@x)   # catches NA, NaN, Inf, -Inf
  if (any(bad)) {
    message(sprintf("  WARNING: %d non-finite values in counts (NA/Inf from qnbinom) — set to 0",
                    sum(bad)))
    counts@x[bad] <- 0L
    counts <- drop0(counts)      # remove structural zeros introduced above
  }
} else {
  bad <- !is.finite(counts)
  if (any(bad)) {
    message(sprintf("  WARNING: %d non-finite values in counts — set to 0", sum(bad)))
    counts[bad] <- 0L
  }
}
ls_vec <- colSums(counts)
message(sprintf("  lib_size: min=%.0f  median=%.0f  max=%.0f  (%.1f%% zero-count cells)",
                min(ls_vec), median(ls_vec), max(ls_vec), 100 * mean(ls_vec == 0)))


message(sprintf("  Simulated: %d genes x %d cells  (%.1f s)",
                nrow(counts), ncol(counts), t_sim_elapsed))

# Free heavy objects no longer needed (ref_data retained only as long as needed above)
rm(fit_mod, modified_marginal, sim_result, sim_para, ref_data)
gc(verbose = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 6. HELPER: tidy results and join ground truth
# ─────────────────────────────────────────────────────────────────────────────
tidy_join <- function(gene, logfc, pval) {
  df <- data.frame(gene      = gene,
                   logFC_int = logfc,
                   pval_int  = pval,
                   padj_int  = p.adjust(pval, method = "BH"),
                   stringsAsFactors = FALSE)
  left_join(df, truth_df, by = "gene")
}

# ─────────────────────────────────────────────────────────────────────────────
# 7. NEBULA
# ─────────────────────────────────────────────────────────────────────────────
t_neb <- proc.time()
message(sprintf("  Running NEBULA-%s...", opt$nebula_method))

run_nebula <- function(counts, new_covariate, method, n_cores,
                       offset_vec = NULL) {
  # Build offset if not supplied.
  # Independent simulation (no copula) creates extreme library-size variability:
  # top HVGs dominate some cells, making ratio-based offsets (lib/median) take
  # values like log(100) = 4.6 that push the NB log-likelihood into NaN.
  # Fix: mean-center on the log scale and clip to [-5, 5] (~150-fold range).
  # As last resort the caller may pass rep(0, ncol(counts)) (zero offset).
  if (is.null(offset_vec)) {
    log_ls     <- log(pmax(colSums(counts), 1L))
    offset_vec <- pmin(pmax(log_ls - mean(log_ls), -5), 5)
  }

  nc <- new_covariate
  nc$visit     <- factor(nc$visit,      levels = c("timepoint1", "timepoint2"))
  nc$treatment <- factor(nc$treatment,  levels = c("groupB",     "groupA"))

  X   <- model.matrix(~ visit * treatment, data = nc)
  neb <- list(count  = counts,
                    id     = nc$subject_id,
                    pred   = X,
                    offset = offset_vec)

  # Identify subjects per arm dynamically from new_covariate (not hardcoded)
  subjects_grpA <- unique(nc$subject_id[nc$treatment == "groupA"])
  subjects_grpB <- unique(nc$subject_id[nc$treatment == "groupB"])

  keep <- sapply(seq_len(nrow(neb$count)), function(g) {
    expr_vec <- neb$count[g, ]
    subA <- length(unique(neb$id[expr_vec > 0 & neb$id %in% subjects_grpA]))
    subB <- length(unique(neb$id[expr_vec > 0 & neb$id %in% subjects_grpB]))
    (subA >= 2) & (subB >= 2)
  })
  
  neb$count <- neb$count[keep, ]
  
  nebula(count   = neb$count,
         id      = neb$id,
         pred    = neb$pred,
         offset  = neb$offset,
         method  = method,
         ncore   = n_cores,
         verbose = T)
}


extract_nebula_stats <- function(fit_neb) {
  summ     <- fit_neb$summary
  # NEBULA stores gene names in a $gene column; rownames are just integers 1..n
  gene_vec <- if ("gene" %in% colnames(summ)) summ$gene else rownames(summ)
  lfc_col  <- grep("logFC.*visittimepoint2.*groupA|logFC.*groupA.*visittimepoint2",
                   colnames(summ), value = TRUE)
  if (!length(lfc_col)) lfc_col <- grep(":", colnames(summ), value = TRUE)[1]
  pval_col <- sub("^logFC_", "p_", lfc_col[1])
  tidy_join(gene  = gene_vec,
            logfc = summ[[lfc_col[1]]],
            pval  = summ[[pval_col]])
}

nebula_stats <- tryCatch({

  fit_neb <- run_nebula(counts, new_covariate, opt$nebula_method, opt$n_cores)
  extract_nebula_stats(fit_neb)

}, error = function(e) {
  # Fallback 1: switch method (LN <-> HL), serial to avoid parallel noise
  fb1 <- if (opt$nebula_method == "LN") "HL" else "LN"
  message(sprintf("  NEBULA-%s failed; retrying %s ncore=1  [%s]",
                  opt$nebula_method, fb1, conditionMessage(e)))
  tryCatch({
    r <- extract_nebula_stats(run_nebula(counts, new_covariate, fb1, 1L))
    attr(r, "nebula_fallback") <- fb1
    r
  }, error = function(e2) {
    # Fallback 2: same alternate method but zero offset.
    # Independent simulation produces extreme library-size variance; dropping
    # the offset removes that source of NaN gradients entirely.
    message(sprintf("  NEBULA-%s also failed; retrying zero offset  [%s]",
                    fb1, conditionMessage(e2)))
    tryCatch({
      r2 <- extract_nebula_stats(
              run_nebula(counts, new_covariate, fb1, 1L,
                         offset_vec = rep(0, ncol(counts))))
      attr(r2, "nebula_fallback") <- paste0(fb1, "_zero_offset")
      r2
    }, error = function(e3) {
      message("  NEBULA all fallbacks failed: ", conditionMessage(e3))
      NULL
    })
  })
})

t_neb_elapsed <- (proc.time() - t_neb)["elapsed"]
message(sprintf("  NEBULA done: %.1f s", t_neb_elapsed))

# ─────────────────────────────────────────────────────────────────────────────
# 8. PSEUDOBULK AGGREGATION (shared between DESeq2 and edgeR)
# ─────────────────────────────────────────────────────────────────────────────
message("  Aggregating pseudobulk...")

# One pseudobulk sample per subject x visit
new_covariate$pb_id <- paste(new_covariate$subject_id,
                             new_covariate$visit, sep = "__")
unique_pbs <- unique(new_covariate$pb_id)

pb_counts <- vapply(unique_pbs, function(pb) {
  idx <- which(new_covariate$pb_id == pb)
  rowSums(counts[, idx, drop = FALSE])
}, numeric(nrow(counts)))

# Pseudobulk sample metadata
pb_meta <- new_covariate[match(unique_pbs, new_covariate$pb_id),
                         c("pb_id", "subject_id", "visit", "treatment")]
rownames(pb_meta) <- unique_pbs
pb_meta$visit     <- factor(pb_meta$visit,
                            levels = c("timepoint1", "timepoint2"))
pb_meta$treatment <- factor(pb_meta$treatment,
                            levels = c("groupB", "groupA"))

# Filter low-count genes
keep           <- rowSums(pb_counts >= opt$pb_min_count) >= opt$pb_min_samples
pb_counts_filt <- pb_counts[keep, ]
message(sprintf("  Pseudobulk: %d samples | %d/%d genes pass filter",
                ncol(pb_counts), sum(keep), nrow(pb_counts)))

# Design matrix: paired pseudobulk
#
# CORRECT FORMULA FOR A PAIRED 2x2 DESIGN:
#   ~ subject_id + visit + visit:treatment
#
# Why NOT ~ subject_id + visit * treatment:
#   Each subject belongs to exactly one treatment arm for their entire study
#   duration, so treatment is a linear combination of the subject_id dummies.
#   Adding a treatment main effect makes the matrix rank-deficient and
#   DESeq2/edgeR will refuse to fit.
#
# What the terms mean:
#   subject_id       absorbs all between-subject variance, INCLUDING the
#                    arm-level mean difference (i.e. the treatment main effect)
#   visit            the within-subject visit effect averaged across arms
#   visit:treatment  the INTERACTION -- how much the visit effect differs
#                    between groupA and groupB -- this is the estimand of interest
pb_design <- model.matrix(~ subject_id + visit + visit:treatment, data = pb_meta)

# ─────────────────────────────────────────────────────────────────────────────
# 9. DESeq2
# ─────────────────────────────────────────────────────────────────────────────
t_d2 <- proc.time()
message("  Running DESeq2...")

deseq2_stats <- tryCatch({
  dds <- DESeqDataSetFromMatrix(countData = pb_counts_filt,
                                colData   = pb_meta,
                                design    = ~ subject_id + visit + visit:treatment)
  
  dds <- DESeq(dds, parallel = TRUE,
               BPPARAM  = BiocParallel::MulticoreParam(opt$n_cores),
               quiet    = TRUE)
  
  # Find interaction result name
  rn       <- resultsNames(dds)
  int_name <- grep("visittimepoint2.*groupA|groupA.*visittimepoint2",
                   rn, value = TRUE)
  if (!length(int_name)) {
    message("  WARNING: DESeq2 interaction term not found; using last coef: ",
            tail(rn, 1))
    int_name <- tail(rn, 1)
  }
  
  res <- results(dds, name = int_name[1], independentFiltering = FALSE)
  
  # Return stats for ALL genes (NAs for filtered-out genes)
  full <- data.frame(gene      = hvg_genes,
                     logFC_int = NA_real_,
                     pval_int  = NA_real_,
                     stringsAsFactors = FALSE)
  hit                 <- match(rownames(res), full$gene)
  full$logFC_int[hit] <- res$log2FoldChange
  full$pval_int[hit]  <- res$pvalue
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")
  
}, error = function(e) {
  message("  DESeq2 error: ", conditionMessage(e))
  NULL
})

t_d2_elapsed <- (proc.time() - t_d2)["elapsed"]
message(sprintf("  DESeq2 done: %.1f s", t_d2_elapsed))

# ─────────────────────────────────────────────────────────────────────────────
# 10. edgeR (quasi-likelihood)
# ─────────────────────────────────────────────────────────────────────────────
t_er <- proc.time()
message("  Running edgeR...")

edger_stats <- tryCatch({
  dge <- DGEList(counts = pb_counts_filt)
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design = pb_design)
  
  fit_er  <- glmQLFit(dge, design = pb_design, robust = TRUE)
  
  # Identify interaction column
  int_col <- grep("visittimepoint2.*groupA|groupA.*visittimepoint2",
                  colnames(pb_design))
  if (!length(int_col)) {
    message("  WARNING: edgeR interaction column not found; using last column.")
    int_col <- ncol(pb_design)
  }
  
  qlf    <- glmQLFTest(fit_er, coef = int_col[1])
  res_er <- topTags(qlf, n = Inf, sort.by = "none")$table
  
  full <- data.frame(gene      = hvg_genes,
                     logFC_int = NA_real_,
                     pval_int  = NA_real_,
                     stringsAsFactors = FALSE)
  hit                 <- match(rownames(res_er), full$gene)
  full$logFC_int[hit] <- res_er$logFC
  full$pval_int[hit]  <- res_er$PValue
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")
  
}, error = function(e) {
  message("  edgeR error: ", conditionMessage(e))
  NULL
})

t_er_elapsed <- (proc.time() - t_er)["elapsed"]
message(sprintf("  edgeR done: %.1f s", t_er_elapsed))

# ─────────────────────────────────────────────────────────────────────────────
# 11. BUILD SEURAT OBJECT FOR CELL-LEVEL METHODS (Wilcoxon, MAST)
#
# Both FindMarkers methods compare groupA POST cells vs groupB POST cells.
# NOTE: This contrast tests (b_trt + b_int), NOT the pure interaction b_int.
# For null scenarios (prop_de = 0, b_int = 0), any detections reflect the
# main treatment effect in the reference model (illustrating type I error
# inflation from ignoring the paired structure).
# ─────────────────────────────────────────────────────────────────────────────
message("  Building Seurat object for cell-level methods...")
so_sim <- CreateSeuratObject(
  counts    = counts,
  meta.data = new_covariate,
  min.cells = 0, min.features = 0
)
so_sim <- NormalizeData(so_sim, verbose = FALSE)

# Subset to POST cells only; set identity to treatment group
so_post <- subset(so_sim, subset = visit == "timepoint2")
Idents(so_post) <- "treatment"
n_post <- ncol(so_post)
message(sprintf("  POST cells: %d  (groupA=%d, groupB=%d)",
                n_post,
                sum(so_post$treatment == "groupA"),
                sum(so_post$treatment == "groupB")))

# ─────────────────────────────────────────────────────────────────────────────
# 12. FINDMARKERS — Wilcoxon (naive, treats each cell as independent)
# ─────────────────────────────────────────────────────────────────────────────
t_wilcox <- proc.time()
message("  Running FindMarkers (Wilcoxon)...")

wilcox_stats <- tryCatch({
  fm_w <- FindMarkers(
    so_post,
    ident.1         = "groupA",
    ident.2         = "groupB",
    test.use        = "wilcox",
    logfc.threshold = 0,    # no pre-filtering: we evaluate all 2k genes
    min.pct         = 0,
    verbose         = FALSE
  )

  full <- data.frame(gene = hvg_genes,
                     logFC_int = NA_real_, pval_int = NA_real_,
                     stringsAsFactors = FALSE)
  hit                 <- match(rownames(fm_w), full$gene)
  full$logFC_int[hit] <- fm_w$avg_log2FC
  full$pval_int[hit]  <- fm_w$p_val
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")

}, error = function(e) {
  message("  Wilcoxon error: ", conditionMessage(e))
  NULL
})

t_wilcox_elapsed <- (proc.time() - t_wilcox)["elapsed"]
message(sprintf("  Wilcoxon done: %.1f s", t_wilcox_elapsed))

# ─────────────────────────────────────────────────────────────────────────────
# 13. FINDMARKERS — MAST (hurdle model; subject_id as fixed latent covariate)
#
# Seurat's latent.vars adds subject_id as a fixed covariate in the continuous
# component of the MAST hurdle model — partial control for subject effects,
# NOT a true random effect.
#
# Cell subsampling: MAST is slow at high cell counts. If POST cells exceed
# mast_max_cells, randomly subsample to that limit (stratified by treatment).
# ─────────────────────────────────────────────────────────────────────────────
t_mast <- proc.time()
message(sprintf("  Running FindMarkers (MAST, max_cells=%d)...", opt$mast_max_cells))

mast_stats <- tryCatch({
  so_mast <- so_post  # may be subsampled below

  if (n_post > opt$mast_max_cells) {
    n_per_group <- floor(opt$mast_max_cells / 2)
    cells_A <- WhichCells(so_post, idents = "groupA")
    cells_B <- WhichCells(so_post, idents = "groupB")
    keep <- c(
      sample(cells_A, min(n_per_group, length(cells_A))),
      sample(cells_B, min(n_per_group, length(cells_B)))
    )
    so_mast <- subset(so_post, cells = keep)
    message(sprintf("  MAST: subsampled to %d POST cells (%d/group)",
                    ncol(so_mast), n_per_group))
  }

  fm_m <- FindMarkers(
    so_mast,
    ident.1         = "groupA",
    ident.2         = "groupB",
    test.use        = "MAST",
    latent.vars     = "subject_id",
    logfc.threshold = 0,
    min.pct         = 0,
    verbose         = FALSE
  )

  full <- data.frame(gene = hvg_genes,
                     logFC_int = NA_real_, pval_int = NA_real_,
                     stringsAsFactors = FALSE)
  hit                 <- match(rownames(fm_m), full$gene)
  full$logFC_int[hit] <- fm_m$avg_log2FC
  full$pval_int[hit]  <- fm_m$p_val
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")

}, error = function(e) {
  message("  MAST error: ", conditionMessage(e))
  NULL
})

t_mast_elapsed <- (proc.time() - t_mast)["elapsed"]
message(sprintf("  MAST done: %.1f s", t_mast_elapsed))

# Free Seurat objects
rm(so_sim, so_post)
if (exists("so_mast")) rm(so_mast)
gc(verbose = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 14. SAVE RESULTS TO S3 (small files only; counts already freed)
# ─────────────────────────────────────────────────────────────────────────────
timing <- c(sim    = t_sim_elapsed,
            nebula = t_neb_elapsed,
            deseq2 = t_d2_elapsed,
            edger  = t_er_elapsed,
            wilcox = t_wilcox_elapsed,
            mast   = t_mast_elapsed)

message(sprintf("  Saving results to s3://%s/%s ...", S3_BUCKET, S3_OUT))

s3saveRDS(nebula_stats,
          object = paste0(S3_OUT, "nebula_stats.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(deseq2_stats,
          object = paste0(S3_OUT, "deseq2_stats.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(edger_stats,
          object = paste0(S3_OUT, "edger_stats.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(wilcox_stats,
          object = paste0(S3_OUT, "wilcox_stats.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(mast_stats,
          object = paste0(S3_OUT, "mast_stats.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(timing,
          object = paste0(S3_OUT, "timing.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(truth_df,
          object = paste0(S3_OUT, "truth.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(params,
          object = paste0(S3_OUT, "params.rds"),
          bucket = S3_BUCKET, region = "")

total_elapsed <- sum(timing)
message(sprintf(
  "══ [02] DONE  array=%d  total=%.1f s  (sim=%.1f | neb=%.1f | d2=%.1f | er=%.1f | wilcox=%.1f | mast=%.1f) ══",
  opt$array_id, total_elapsed,
  timing["sim"], timing["nebula"], timing["deseq2"], timing["edger"],
  timing["wilcox"], timing["mast"]
))
