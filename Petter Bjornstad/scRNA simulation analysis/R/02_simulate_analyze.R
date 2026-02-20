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
#     3. Calls simulate_new_data() -- count matrix lives only in memory
#     4. Immediately runs NEBULA, DESeq2, edgeR on the in-memory counts
#     5. Saves three small stats data.frames + timing + truth to S3
#     6. Frees the count matrix (gc())
#
# STUDY DESIGN CLARIFICATION:
#   n_subjects_per_arm = number of SUBJECTS in each treatment arm.
#   Each subject has BOTH a PRE and a POST visit (paired design).
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
# OUTPUT (to S3, per array task, ~1 MB total):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/stats/array_NNNNN/
#     nebula_stats.rds    -- gene x {logFC, pval, padj, is_de, beta_int_true}
#     deseq2_stats.rds    -- same structure
#     edger_stats.rds     -- same structure
#     timing.rds          -- named vector: sim, nebula, deseq2, edger (seconds)
#     params.rds          -- single-row data.frame of this task's parameters
#     truth.rds           -- data.frame: gene, is_de, beta_int_true
#
# USAGE:
#   Rscript 02_simulate_analyze.R \
#       --array_id    1 \
#       --n_cores     4
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
S3_BUCKET  <- "scrna"
S3_BASE    <- "Projects/Paired scRNA simulation analysis/results/"
S3_REF_PFX <- paste0(S3_BASE, "reference/")
S3_GRID_PFX <- paste0(S3_BASE, "param_grid/")

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--array_id",      type = "integer",   default = NULL,
              help = "Row index in param_grid (1-based) [required]"),
  make_option("--n_cores",       type = "integer",   default = 4L),
  make_option("--nebula_method", type = "character", default = "LN",
              help = "NEBULA method: LN or HL"),
  make_option("--pb_min_count",  type = "integer",   default = 10L,
              help = "Min rowSum for pseudobulk filtering"),
  make_option("--pb_min_samples",type = "integer",   default = 2L,
              help = "Min samples meeting pb_min_count threshold")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$array_id)) stop("--array_id is required")

# Seed is array_id so each task is reproducible and unique
set.seed(opt$array_id)
message(sprintf("  Seed set to: %d (= array_id)", opt$array_id))
register(MulticoreParam(opt$n_cores))

# Per-task S3 output prefix
arr_str  <- sprintf("array_%05d", opt$array_id)
S3_OUT   <- paste0(S3_BASE, "stats/", arr_str, "/")

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
fit_obj   <- s3readRDS(object = paste0(S3_REF_PFX, "scdesign3_fit.rds"),
                       bucket = S3_BUCKET, region = "")
hvg_genes <- s3readRDS(object = paste0(S3_REF_PFX, "hvg_genes.rds"),
                       bucket = S3_BUCKET, region = "")
sce_ref   <- s3readRDS(object = paste0(S3_REF_PFX, "sce_ref.rds"),
                       bucket = S3_BUCKET, region = "")

n_genes <- length(hvg_genes)  # 2000

message(sprintf("  Scenario: %s  rep=%d", params$scenario_id, params$sim_rep))

# ─────────────────────────────────────────────────────────────────────────────
# 2. BUILD CELL-LEVEL DESIGN
#
# n_subjects_per_arm subjects in each of 2 arms, each with PRE + POST visit.
# Cells per subject drawn from N(cells_mean, cells_sd^2), clipped to >= 50.
# ─────────────────────────────────────────────────────────────────────────────
treatments <- c("Placebo", "Dapagliflozin")
visits     <- c("PRE", "POST")

make_subject_cells <- function(subj_id, trt, cells_mean, cells_sd) {
  # Each subject has cells from BOTH visits; split roughly 50/50 between visits
  n_total <- max(50L, round(rnorm(1, cells_mean, cells_sd)))
  # Assign visits proportionally (slightly random split)
  n_pre   <- round(n_total * runif(1, 0.4, 0.6))
  n_post  <- n_total - n_pre
  data.frame(
    subject_id = rep(subj_id, n_total),
    treatment  = rep(trt,     n_total),
    visit      = c(rep("PRE", n_pre), rep("POST", n_post)),
    stringsAsFactors = FALSE
  )
}

cell_list <- lapply(treatments, function(trt) {
  lapply(seq_len(params$n_subjects_per_arm), function(j) {
    # Subject ID encodes arm and subject number
    subj <- sprintf("%s_S%02d",
                    ifelse(trt == "Dapagliflozin", "Dapa", "Plac"), j)
    make_subject_cells(subj, trt, params$cells_mean, params$cells_sd)
  })
})
new_covariate <- do.call(rbind, do.call(c, cell_list))
new_covariate$visit     <- factor(new_covariate$visit,
                                  levels = c("PRE", "POST"))
new_covariate$treatment <- factor(new_covariate$treatment,
                                  levels = c("Placebo", "Dapagliflozin"))
new_covariate$celltype  <- unique(sce_ref$celltype)[1]

n_cells_total <- nrow(new_covariate)
n_subjects    <- length(unique(new_covariate$subject_id))
message(sprintf("  Cells: %d  |  Subjects: %d (%d per arm x 2 arms)",
                n_cells_total, n_subjects, params$n_subjects_per_arm))

# ─────────────────────────────────────────────────────────────────────────────
# 3. GROUND TRUTH: ASSIGN DE GENES AND TRUE INTERACTION EFFECTS
# ─────────────────────────────────────────────────────────────────────────────
n_de    <- round(params$prop_de * n_genes)
de_idx  <- if (n_de > 0) sample(seq_len(n_genes), n_de) else integer(0)

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
# ─────────────────────────────────────────────────────────────────────────────
# The interaction coefficient name in scDesign3/mgcv output:
COEF_INT <- "visitPOST:treatmentDapagliflozin"

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

build_cs_corr <- function(n, rho) {
  S      <- matrix(rho, n, n)
  diag(S) <- 1.0
  S
}

message("  Modifying scDesign3 model (DE injection + individual variance)...")
modified_marginal <- bplapply(
  seq_along(fit_obj$marginal_list),
  function(g) modify_gene_model(fit_obj$marginal_list[[g]],
                                beta_int  = beta_int_true[g],
                                iv_factor = params$indiv_var_factor),
  BPPARAM = MulticoreParam(opt$n_cores)
)
names(modified_marginal) <- names(fit_obj$marginal_list)

# Rebuild copula with requested within-subject cell correlation
modified_copula <- fit_obj$copula_list
if (!is.null(modified_copula)) {
  for (subj in names(modified_copula)) {
    if (!is.null(modified_copula[[subj]]$corr)) {
      nc <- nrow(modified_copula[[subj]]$corr)
      modified_copula[[subj]]$corr <- build_cs_corr(nc, params$corr_cells)
    }
  }
}

fit_mod <- fit_obj
fit_mod$marginal_list <- modified_marginal
if (!is.null(modified_copula)) fit_mod$copula_list <- modified_copula

# ─────────────────────────────────────────────────────────────────────────────
# 5. SIMULATE COUNTS (held in memory only)
# ─────────────────────────────────────────────────────────────────────────────
t_sim <- proc.time()
message("  Simulating counts...")

sim_result <- tryCatch(
  simulate_new_data(
    fitted_model      = fit_mod,
    n_cores           = opt$n_cores,
    parallelization   = "pbmcapply",
    BPPARAM           = MulticoreParam(opt$n_cores),
    new_covariate     = new_covariate,
    family_use        = "nb",
    important_feature = "all",
    nonnegative       = TRUE,
    set_zero_prob     = TRUE
  ),
  error = function(e) {
    message("  ERROR in simulate_new_data: ", conditionMessage(e))
    NULL
  }
)

t_sim_elapsed <- (proc.time() - t_sim)["elapsed"]

if (is.null(sim_result)) {
  s3saveRDS(list(error = "simulate_new_data failed",
                  params = params,
                  array_id = opt$array_id),
             object = paste0(S3_OUT, "error.rds"),
             bucket = S3_BUCKET,
             region = "")
  quit(status = 1)
}

counts <- sim_result$new_count   # genes x cells, in memory
if (is.null(rownames(counts))) rownames(counts) <- hvg_genes
message(sprintf("  Simulated: %d genes x %d cells  (%.1f s)",
                nrow(counts), ncol(counts), t_sim_elapsed))

# Free heavy objects no longer needed
rm(fit_mod, modified_marginal, modified_copula, sim_result)
gc(verbose = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 6. HELPER: tidy results and join ground truth
# ─────────────────────────────────────────────────────────────────────────────
tidy_join <- function(gene, logfc, pval) {
  df <- data.frame(gene     = gene,
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

nebula_stats <- tryCatch({
  lib_size    <- colSums(counts)
  size_factor <- lib_size / median(lib_size)

  new_covariate$visit     <- factor(new_covariate$visit,     levels = c("PRE","POST"))
  new_covariate$treatment <- factor(new_covariate$treatment, levels = c("Placebo","Dapagliflozin"))

  X   <- model.matrix(~ visit * treatment, data = new_covariate)
  neb <- group_cell(count  = counts,
                    id     = new_covariate$subject_id,
                    pred   = X,
                    offset = log(size_factor))

  fit_neb <- nebula(count   = neb$count,
                    id      = neb$id,
                    pred    = neb$pred,
                    offset  = neb$offset,
                    method  = opt$nebula_method,
                    ncore   = opt$n_cores,
                    verbose = FALSE)

  summ     <- fit_neb$summary
  # Robustly find interaction column (visitPOST:treatmentDapagliflozin)
  lfc_col  <- grep("logFC.*visitPOST.*Dapa|logFC.*Dapa.*visitPOST",
                   colnames(summ), value = TRUE)
  if (!length(lfc_col)) lfc_col <- grep(":", colnames(summ), value = TRUE)[1]
  pval_col <- sub("^logFC_", "p_", lfc_col[1])

  tidy_join(gene  = rownames(summ),
            logfc = summ[[lfc_col[1]]],
            pval  = summ[[pval_col]])
}, error = function(e) {
  message("  NEBULA error: ", conditionMessage(e))
  NULL
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
                         c("pb_id","subject_id","visit","treatment")]
rownames(pb_meta) <- unique_pbs
pb_meta$visit     <- factor(pb_meta$visit,     levels = c("PRE","POST"))
pb_meta$treatment <- factor(pb_meta$treatment, levels = c("Placebo","Dapagliflozin"))

# Filter low-count genes
keep <- rowSums(pb_counts >= opt$pb_min_count) >= opt$pb_min_samples
pb_counts_filt <- pb_counts[keep, ]
message(sprintf("  Pseudobulk: %d samples | %d/%d genes pass filter",
                ncol(pb_counts), sum(keep), nrow(pb_counts)))

# Design matrix: paired pseudobulk (subject_id controls within-subject pairing)
pb_design <- model.matrix(~ subject_id + visit * treatment, data = pb_meta)

# ─────────────────────────────────────────────────────────────────────────────
# 9. DESeq2
# ─────────────────────────────────────────────────────────────────────────────
t_d2 <- proc.time()
message("  Running DESeq2...")

deseq2_stats <- tryCatch({
  dds <- DESeqDataSetFromMatrix(countData = pb_counts_filt,
                                colData   = pb_meta,
                                design    = ~ subject_id + visit * treatment)

  dds <- DESeq(dds, parallel = TRUE,
               BPPARAM  = MulticoreParam(opt$n_cores),
               quiet    = TRUE)

  # Find interaction result name
  rn       <- resultsNames(dds)
  int_name <- grep("visitPOST.*Dapa|Dapa.*visitPOST", rn, value = TRUE)
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
  hit                <- match(rownames(res), full$gene)
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
  int_col <- grep("visitPOST.*Dapa|Dapa.*visitPOST", colnames(pb_design))
  if (!length(int_col)) {
    message("  WARNING: edgeR interaction column not found; using last column.")
    int_col <- ncol(pb_design)
  }

  qlf     <- glmQLFTest(fit_er, coef = int_col[1])
  res_er  <- topTags(qlf, n = Inf, sort.by = "none")$table

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
# 11. SAVE RESULTS TO S3 (small files only; counts already freed)
# ─────────────────────────────────────────────────────────────────────────────
timing <- c(sim    = t_sim_elapsed,
            nebula = t_neb_elapsed,
            deseq2 = t_d2_elapsed,
            edger  = t_er_elapsed)

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
  "══ [02] DONE  array=%d  total=%.1f s  (sim=%.1f | neb=%.1f | d2=%.1f | er=%.1f) ══",
  opt$array_id, total_elapsed,
  timing["sim"], timing["nebula"], timing["deseq2"], timing["edger"]
))
