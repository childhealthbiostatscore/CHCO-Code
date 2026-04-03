################################################################################
# 02_simulate_analyze.R
#
# PURPOSE:
#   COMBINED script: simulate one scRNA-seq dataset, run all six DE methods
#   (NEBULA-NBLMM, edgeR pseudobulk, limma-voom + duplicateCorrelation,
#   Wilcoxon, MAST), save only the small stats output files to S3, and
#   discard the raw count matrix.
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
#     4. Immediately runs NEBULA (gene-by-gene NBLMM+REML), edgeR,
#        limma-voom+duplicateCorrelation, Wilcoxon, MAST on the in-memory counts
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
#     edger_stats.rds     -- same structure
#     limma_stats.rds     -- same structure (limma-voom + duplicateCorrelation)
#     wilcox_stats.rds    -- same structure (contrast: Treatment POST vs Placebo POST)
#     mast_stats.rds      -- same structure (MAST + subject_id latent covariate)
#     agreement.rds       -- pairwise Jaccard / overlap coeff between methods
#     sig_membership.rds  -- binary gene x method matrix of sig hits (for UpSet)
#     timing.rds          -- named vector: sim, nebula, edger, limma, wilcox, mast (s)
#     params.rds          -- single-row data.frame of this task's parameters
#     truth.rds           -- data.frame: gene, is_de, beta_int_true
#
# NOTE on Wilcoxon / MAST contrast:
#   These methods compare Treatment POST cells vs Placebo POST cells, which tests
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
  library(edgeR)
  library(limma)
  library(scran)
  library(Matrix)
  library(dplyr)
  library(BiocParallel)
  library(optparse)
  library(aws.s3)
  library(Seurat)
  library(MAST)
  library(foreach)
  library(doParallel)
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
  make_option("--nebula_method",   type = "character", default = "NBLMM",
              help = "NEBULA model (now uses gene-by-gene NBLMM+REML; kept for backward compat)"),
  make_option("--pb_min_count",    type = "integer",   default = 10L,
              help = "Min rowSum for pseudobulk filtering"),
  make_option("--pb_min_samples",  type = "integer",   default = 2L,
              help = "Min samples meeting pb_min_count threshold"),
  make_option("--mast_max_cells",  type = "integer",   default = 5000L,
              help = "Max POST cells to pass to MAST (subsample if exceeded) [default: 5000]")
)
opt <- parse_args(OptionParser(option_list = option_list))
# manual settings for debugging
# opt$array_id = 1
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
# n_subjects_per_arm subjects in each of 2 arms (Placebo / Treatment),
# each with timepoint1 + timepoint2 visit (paired design).
# Cells per subject drawn from N(cells_mean, cells_sd^2), clipped to >= 50.
# ─────────────────────────────────────────────────────────────────────────────
treatments <- c("Placebo", "Treatment")
visits     <- c("timepoint1", "timepoint2")

# grab the exact subject IDs the copula knows about
ref_subjects <- levels(ref_data$ref_construct_data$newCovariate$subject_id)
# e.g. c("grpPlacebo_S01", "grpPlacebo_S02", ..., "grpTreatment_S01", ...)

n_ref_per_arm <- sum(grepl("^grpPlacebo", ref_subjects))

make_subject_cells <- function(subj_id, trt, cells_mean, cells_sd) {
  n_total <- max(50L, round(rnorm(1, cells_mean, cells_sd)))
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
  prefix <- ifelse(trt == "Treatment", "grpTreatment", "grpPlacebo")
  lapply(seq_len(params$n_subjects_per_arm), function(j) {
    subj <- sprintf("%s_S%02d", prefix, j)
    make_subject_cells(subj, trt, params$cells_mean, params$cells_sd)
  })
})

new_covariate <- do.call(rbind, do.call(c, cell_list))
new_covariate$visit     <- factor(new_covariate$visit,
                                  levels = c("timepoint1", "timepoint2"))
new_covariate$treatment <- factor(new_covariate$treatment,
                                  levels = c("Placebo", "Treatment"))
new_covariate$celltype  <- unique(sce_ref$celltype)[1]

n_cells_total <- nrow(new_covariate)
n_subjects    <- length(unique(new_covariate$subject_id))
message(sprintf("  Cells: %d  |  Subjects: %d (%d per arm x 2 arms)",
                n_cells_total, n_subjects, params$n_subjects_per_arm))

# Per-subject cell counts (attach to params for downstream diagnostics)
cells_per_subject <- table(new_covariate$subject_id)
params$n_cells_total    <- n_cells_total
params$cells_per_subject <- list(setNames(as.integer(cells_per_subject),
                                          names(cells_per_subject)))
params$cells_min        <- min(cells_per_subject)
params$cells_median     <- median(as.integer(cells_per_subject))
params$cells_max        <- max(cells_per_subject)

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
  log2FC_true   = beta_int_true / log(2),   # convert natural log -> log2 for bias
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────────────────
# 4. MODIFY scDesign3 MODEL: inject DE signal + individual-level variance
#
# Factor levels in sce_ref colData are:
#   visit:     c("timepoint1", "timepoint2")  -> reference = "timepoint1"
#   treatment: c("Placebo", "Treatment")          -> reference = "Placebo"
# So the interaction coefficient is: "visittimepoint2:treatmentTreatment"
# ─────────────────────────────────────────────────────────────────────────────
COEF_INT <- "visittimepoint2:treatmentTreatment"

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
  seq_along(ref_data$ref_fit_marginal),
  function(g) modify_gene_model(ref_data$ref_fit_marginal[[g]],
                                beta_int  = beta_int_true[g],
                                iv_factor = params$indiv_var_factor),
  BPPARAM = BiocParallel::MulticoreParam(opt$n_cores)
)
names(modified_marginal) <- names(ref_data$ref_fit_marginal)

# The fitted copula from 00_fit_reference.R is keyed on reference subject IDs.
# When the simulation has more subjects than the reference, we cycle through
# reference copulas and clone input_data rows for the extra subjects.
# This preserves gene-gene correlation structure while allowing flexible sample sizes.

# ─────────────────────────────────────────────────────────────────────────────
# 5. SIMULATE COUNTS (held in memory only)
#
# scDesign3 modular API:
#   extract_para()  -- converts marginal_list -> parameter matrices
#   simu_new()      -- uses parameter matrices + copula to generate counts
# ─────────────────────────────────────────────────────────────────────────────
t_sim <- proc.time()
message("  Extracting parameters from modified model...")
new_covariate_for_sim <- new_covariate[, colnames(ref_data$ref_construct_data$newCovariate)]
new_covariate_for_sim <- new_covariate_for_sim %>%
  dplyr::mutate(corr_group = subject_id)

sim_para <- tryCatch(
  extract_para(
    sce           = sce_ref,
    marginal_list = modified_marginal,
    n_cores       = opt$n_cores,
    family_use    = "nb",
    new_covariate = new_covariate_for_sim,
    data          = ref_data$ref_construct_data$dat
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

message("  Remap copula list and input data for new subject design...")
ref_copula   <- ref_data$ref_fit_copula$copula_list
ref_names    <- names(ref_copula)
ref_dat_full <- ref_data$ref_construct_data$dat %>%
  dplyr::mutate(corr_group = as.character(corr_group))

# New subjects from the simulation design
new_subjects <- as.character(unique(new_covariate$subject_id))

# Map each new subject to a reference subject by cycling through available refs.
# Works whether n_new < n_ref (subset), n_new == n_ref (1:1), or n_new > n_ref (recycle).
ref_assignment <- setNames(
  ref_names[((seq_along(new_subjects) - 1) %% length(ref_names)) + 1],
  new_subjects
)
message(sprintf("  Mapping %d new subjects -> %d reference copulas (cycling)",
                length(new_subjects), length(ref_names)))

# Build mapped_copula: one entry per new subject, copied from assigned ref
mapped_copula <- setNames(
  lapply(new_subjects, function(s) ref_copula[[ ref_assignment[[s]] ]]),
  new_subjects
)

# Build input_data: for each new subject, we need rows in input_data so that
# simu_new() can iterate over unique(input_data$corr_group).
#   - If the new subject exists in the reference data, use its real rows.
#   - If it doesn't (n_new > n_ref), clone rows from the assigned donor and relabel.
ref_subject_ids <- unique(as.character(ref_dat_full$subject_id))
input_data_list <- lapply(new_subjects, function(s) {
  if (s %in% ref_subject_ids) {
    # Subject exists in reference — use real rows
    ref_dat_full %>% dplyr::filter(subject_id == s)
  } else {
    # Subject not in reference — clone from the assigned donor
    donor <- ref_assignment[[s]]
    ref_dat_full %>%
      dplyr::filter(subject_id == donor) %>%
      dplyr::mutate(subject_id = s,
                    corr_group = s)
  }
})
ref_data_dat <- do.call(rbind, input_data_list)

message("  Simulating counts with simu_new()...")

sim_result <- simu_new(
  sce           = sce_ref,
  mean_mat      = sim_para$mean_mat,
  sigma_mat     = sim_para$sigma_mat,
  zero_mat      = sim_para$zero_mat,
  copula_list   = mapped_copula,
  n_cores       = opt$n_cores,
  family_use    = "nb",
  important_feature = ref_data$ref_fit_copula$important_feature,
  input_data    = ref_data_dat,
  filtered_gene = ref_data$ref_construct_data$filtered_gene,
  new_covariate = new_covariate_for_sim
)

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
rm(modified_marginal, sim_result, sim_para, ref_data, ref_data_dat, mapped_copula)
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
# 7. NEBULA  (gene-by-gene NBLMM + REML — matches ATTEMPT reference analysis)
# ─────────────────────────────────────────────────────────────────────────────

# Fixed version of nebula::group_cell that properly handles factor conversion.
# The original group_cell() has a bug: is.unsorted() on a factor can return
# FALSE even when cells aren't sorted, causing all subjects to collapse to 1.
group_cell_mod <- function(count, id, pred = NULL, offset = NULL) {
  ng <- nrow(count)
  nc <- ncol(count)
  if (nc != length(id)) {
    stop("The length of id is not equal to the number of columns of the count matrix.")
  }
  id <- as.character(id)
  levels <- unique(id)
  id <- factor(id, levels = levels)
  if (is.unsorted(id) == FALSE) {
    cat("The cells are already grouped.")
    return(NULL)
  }
  k <- length(levels)
  o <- order(id)
  count <- count[, o]
  id <- id[o]
  if (is.null(pred) == FALSE) {
    if (nc != nrow(as.matrix(pred))) {
      stop("The number of rows of the design matrix is not equal to the number of columns of the count matrix")
    }
    pred <- as.matrix(pred)[o, , drop = FALSE]
  }
  if (is.null(offset) == FALSE) {
    if (nc != length(offset)) {
      stop("The length of offset is not equal to the number of columns of the count matrix")
    }
    if (sum(offset <= 0) > 0) {
      stop("Some elements in the scaling factor are not positive.")
    }
    offset <- offset[o]
  }
  grouped <- list(count = count, id = id, pred = pred, offset = offset)
  return(grouped)
}

t_neb <- proc.time()
message("  Running NEBULA (gene-by-gene NBLMM + REML)...")

nebula_stats <- tryCatch({
  nc <- new_covariate
  nc$visit      <- factor(nc$visit,      levels = c("timepoint1", "timepoint2"))
  nc$treatment  <- factor(nc$treatment,  levels = c("Placebo",     "Treatment"))
  nc$subject_id <- as.character(nc$subject_id)

  counts_rounded <- round(counts)
  genes_list     <- rownames(counts_rounded)

  # Pooled offset via scran (matches ATTEMPT QMD approach)
  sce_offset <- SingleCellExperiment(assays = list(counts = counts_rounded))
  sce_offset <- computeSumFactors(sce_offset,
                                  BPPARAM = BiocParallel::MulticoreParam(opt$n_cores))
  pooled_offset <- sizeFactors(sce_offset)
  rm(sce_offset); gc(verbose = FALSE)
  message(sprintf("  Pooled offset: min=%.4f  median=%.4f  max=%.4f",
                  min(pooled_offset), median(pooled_offset), max(pooled_offset)))

  message(sprintf("  NEBULA: fitting %d genes across %d cells, %d subjects",
                  length(genes_list), ncol(counts_rounded),
                  length(unique(nc$subject_id))))

  # Gene-by-gene parallel loop (matches ATTEMPT QMD / 05_volcano_reference.R)
  cl <- makeCluster(opt$n_cores)
  registerDoParallel(cl)

  nebula_res_list <- foreach(
    g = genes_list,
    .packages = c("nebula", "Matrix"),
    .export   = c("counts_rounded", "nc", "pooled_offset", "group_cell_mod"),
    .errorhandling = "pass"
  ) %dopar% {
    warn <- err <- res <- NULL
    tryCatch({
      count_gene <- counts_rounded[g, , drop = FALSE]
      pred_gene  <- model.matrix(~ visit * treatment, data = nc)
      data_gene  <- group_cell_mod(count  = count_gene,
                                   id     = nc$subject_id,
                                   pred   = pred_gene,
                                   offset = pooled_offset)

      # group_cell returns NULL when cells appear sorted; use raw inputs
      if (is.null(data_gene)) {
        data_gene <- list(count = count_gene, id = nc$subject_id,
                          pred = pred_gene, offset = pooled_offset)
      }

      res <- withCallingHandlers(
        nebula(
          count      = data_gene$count,
          id         = data_gene$id,
          pred       = data_gene$pred,
          offset     = data_gene$offset,
          ncore      = 1,
          output_re  = TRUE,
          covariance = TRUE,
          reml       = 1,
          model      = "NBLMM"
        ),
        warning = function(w) {
          warn <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      )
    }, error = function(e) {
      err <<- conditionMessage(e)
    })
    list(gene = g, result = res, warning = warn, error = err)
  }

  stopCluster(cl)

  # Extract results — filter out errors (NULL) and non-converged genes
  names(nebula_res_list) <- vapply(nebula_res_list, `[[`, "", "gene")
  result_list <- lapply(nebula_res_list, `[[`, "result")
  result_list <- Filter(Negate(is.null), result_list)

  # Check NEBULA's internal convergence flag (conv == 1 means converged)
  converged <- sapply(result_list, function(r) {
    if (!is.null(r$convergence)) all(r$convergence == 1) else TRUE
  })
  n_returned  <- length(result_list)
  n_converged <- sum(converged)
  n_total     <- length(genes_list)
  result_list <- result_list[converged]

  message(sprintf("  NEBULA: %d/%d returned results, %d converged (%.1f%% failed)",
                  n_returned, n_total, n_converged,
                  100 * (n_total - n_converged) / n_total))

  # Parse each converged gene's NEBULA result into a row
  parsed <- lapply(names(result_list), function(g) {
    res  <- result_list[[g]]
    summ <- res$summary
    # Find the interaction column (visit * treatment formula order)
    lfc_col <- grep("logFC.*visit.*treatment|logFC.*Treatment.*timepoint|logFC.*visittimepoint2.*Treatment",
                    colnames(summ), value = TRUE, ignore.case = TRUE)
    if (!length(lfc_col)) {
      # Fallback: any logFC column with a colon (interaction)
      lfc_col <- grep("^logFC_.*:", colnames(summ), value = TRUE)
    }
    if (!length(lfc_col)) return(NULL)

    pval_col <- sub("^logFC_", "p_", lfc_col[1])
    if (!pval_col %in% colnames(summ)) return(NULL)

    # Additional sanity check: skip genes with extreme logFC (non-converged artifacts)
    lfc_val <- summ[[lfc_col[1]]] / log(2)
    if (!is.finite(lfc_val) || abs(lfc_val) > 20) return(NULL)

    data.frame(
      gene      = g,
      logFC_int = lfc_val,
      pval_int  = summ[[pval_col]],
      stringsAsFactors = FALSE
    )
  })
  parsed <- Filter(Negate(is.null), parsed)

  if (length(parsed) == 0) {
    message("  NEBULA: no genes converged")
    NULL
  } else {
    df <- do.call(rbind, parsed)
    df$padj_int <- p.adjust(df$pval_int, method = "BH")
    # Join with ground truth
    left_join(df, truth_df, by = "gene")
  }

}, error = function(e) {
  message("  NEBULA error: ", conditionMessage(e))
  NULL
})

t_neb_elapsed <- (proc.time() - t_neb)["elapsed"]
message(sprintf("  NEBULA done: %.1f s  (%s genes)",
                t_neb_elapsed, if (!is.null(nebula_stats)) nrow(nebula_stats) else "0"))

# ─────────────────────────────────────────────────────────────────────────────
# 8. PSEUDOBULK AGGREGATION (shared between edgeR and limma)
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
# Sparse-matrix rowSums can produce floating-point artifacts (e.g. 1547.0000000000002);
# round once here so all downstream methods (edgeR, limma) receive clean integers.
pb_counts <- round(pb_counts)

# Pseudobulk sample metadata
pb_meta <- new_covariate[match(unique_pbs, new_covariate$pb_id),
                         c("pb_id", "subject_id", "visit", "treatment")]
rownames(pb_meta) <- unique_pbs
pb_meta$visit     <- factor(pb_meta$visit,
                            levels = c("timepoint1", "timepoint2"))
pb_meta$treatment <- factor(pb_meta$treatment,
                            levels = c("Placebo", "Treatment"))

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
#   edgeR will refuse to fit.
#
# What the terms mean:
#   visit            the within-subject visit effect averaged across arms
#   visit:treatment  the INTERACTION -- how much the visit effect differs
#                    between Treatment and Placebo -- this is the estimand of interest
pb_design <- model.matrix(~ visit + treatment + visit:treatment, data = pb_meta)

# ─────────────────────────────────────────────────────────────────────────────
# 9. edgeR (quasi-likelihood)
# ─────────────────────────────────────────────────────────────────────────────
t_er <- proc.time()
message("  Running edgeR...")

edger_stats <- tryCatch({
  dge <- DGEList(counts = pb_counts_filt)
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design = pb_design)
  
  fit_er  <- glmQLFit(dge, design = pb_design, robust = TRUE)
  
  # Identify interaction column
  int_col <- grep("visittimepoint2.*Treatment|Treatment.*visittimepoint2",
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
# 10. limma-voom + duplicateCorrelation
#
# Uses voom to transform pseudobulk counts, then estimates the within-subject
# correlation via duplicateCorrelation() with subject_id as the blocking factor.
# The design matrix is ~ visit * treatment (full interaction, no subject_id
# dummies), and the estimated correlation is passed to both lmFit calls.
#
# This is the standard limma approach for repeated-measures designs where
# subjects are measured at multiple timepoints.
# ─────────────────────────────────────────────────────────────────────────────
t_limma <- proc.time()
message("  Running limma-voom + duplicateCorrelation...")

limma_stats <- tryCatch({
  # DGEList + normalization (reuse filtered pseudobulk counts)
  dge_limma <- DGEList(counts = pb_counts_filt)
  dge_limma <- calcNormFactors(dge_limma)

  # Design: full interaction (no subject dummies — correlation handles pairing)
  limma_design <- model.matrix(~ visit * treatment, data = pb_meta)

  # voom transform
  v <- limma::voom(dge_limma, design = limma_design, plot = FALSE)

  # Estimate within-subject correlation
  corfit <- limma::duplicateCorrelation(v, limma_design,
                                        block = pb_meta$subject_id)
  message(sprintf("  duplicateCorrelation consensus = %.4f", corfit$consensus))

  # Re-run voom with the correlation estimate for better precision weights
  v <- limma::voom(dge_limma, design = limma_design, plot = FALSE,
                   block = pb_meta$subject_id,
                   correlation = corfit$consensus)

  # Fit with correlation structure
  fit_limma <- limma::lmFit(v, limma_design,
                            block = pb_meta$subject_id,
                            correlation = corfit$consensus)
  fit_limma <- limma::eBayes(fit_limma)

  # Extract interaction term
  int_col <- grep("visittimepoint2.*Treatment|Treatment.*visittimepoint2",
                  colnames(limma_design), value = TRUE)
  if (!length(int_col)) {
    message("  WARNING: limma interaction column not found; using last coef.")
    int_col <- colnames(limma_design)[ncol(limma_design)]
  }

  res_limma <- limma::topTable(fit_limma, coef = int_col[1],
                               number = Inf, sort.by = "none")

  # Build full-gene output (NAs for filtered-out genes)
  full <- data.frame(gene      = hvg_genes,
                     logFC_int = NA_real_,
                     pval_int  = NA_real_,
                     stringsAsFactors = FALSE)
  hit                 <- match(rownames(res_limma), full$gene)
  full$logFC_int[hit] <- res_limma$logFC
  full$pval_int[hit]  <- res_limma$P.Value
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")

}, error = function(e) {
  message("  limma error: ", conditionMessage(e))
  NULL
})

t_limma_elapsed <- (proc.time() - t_limma)["elapsed"]
message(sprintf("  limma done: %.1f s", t_limma_elapsed))

# ─────────────────────────────────────────────────────────────────────────────
# 11. BUILD SEURAT OBJECT FOR CELL-LEVEL METHODS (Wilcoxon, MAST)
#
# Both FindMarkers methods compare Treatment POST cells vs Placebo POST cells.
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
# Increase future.globals.maxSize to avoid error when NormalizeData exports
# large objects to parallel workers (default 500 MiB is too small for big sims)
options(future.globals.maxSize = 2 * 1024^3)  # 2 GiB
so_sim <- NormalizeData(so_sim, verbose = FALSE)

# Subset to POST cells only; set identity to treatment group
so_post <- subset(so_sim, subset = visit == "timepoint2")
Idents(so_post) <- "treatment"
n_post <- ncol(so_post)
message(sprintf("  POST cells: %d  (Treatment=%d, Placebo=%d)",
                n_post,
                sum(so_post$treatment == "Treatment"),
                sum(so_post$treatment == "Placebo")))

# ─────────────────────────────────────────────────────────────────────────────
# 12. FINDMARKERS — Wilcoxon (naive, treats each cell as independent)
# ─────────────────────────────────────────────────────────────────────────────
t_wilcox <- proc.time()
message("  Running FindMarkers (Wilcoxon)...")

wilcox_stats <- tryCatch({
  fm_w <- FindMarkers(
    so_post,
    ident.1         = "Treatment",
    ident.2         = "Placebo",
    test.use        = "wilcox",
    logfc.threshold = 0,    # no pre-filtering: we evaluate all genes
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
    cells_A <- WhichCells(so_post, idents = "Treatment")
    cells_B <- WhichCells(so_post, idents = "Placebo")
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
    ident.1         = "Treatment",
    ident.2         = "Placebo",
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
# 13. METHOD AGREEMENT: pairwise overlap of significant gene sets
#
# For each method pair, compute:
#   - Jaccard index  = |A ∩ B| / |A ∪ B|
#   - Overlap coeff  = |A ∩ B| / min(|A|, |B|)
#   - Counts: n_A, n_B, n_intersect, n_union
# Also store the per-method significant gene lists so the downstream
# aggregation script can build UpSet plots across scenarios.
# ─────────────────────────────────────────────────────────────────────────────
message("  Computing method agreement...")

method_results <- list(
  NEBULA   = nebula_stats,
  edgeR    = edger_stats,
  limma    = limma_stats,
  Wilcoxon = wilcox_stats,
  MAST     = mast_stats
)

# Significant genes per method (padj < 0.05)
sig_genes <- lapply(method_results, function(df) {
  if (is.null(df)) return(character(0))
  df$gene[!is.na(df$padj_int) & df$padj_int < 0.05]
})

# Pairwise agreement
method_names <- names(sig_genes)
pairs <- combn(method_names, 2, simplify = FALSE)

agreement <- do.call(rbind, lapply(pairs, function(p) {
  a <- sig_genes[[p[1]]]
  b <- sig_genes[[p[2]]]
  n_a   <- length(a)
  n_b   <- length(b)
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  data.frame(
    method_A    = p[1],
    method_B    = p[2],
    n_A         = n_a,
    n_B         = n_b,
    n_intersect = inter,
    n_union     = uni,
    jaccard     = if (uni > 0) inter / uni else NA_real_,
    overlap_coef = if (min(n_a, n_b) > 0) inter / min(n_a, n_b) else NA_real_,
    stringsAsFactors = FALSE
  )
}))

# Binary membership matrix (genes x methods) for UpSet-style analysis downstream
all_genes <- unique(unlist(sig_genes))
if (length(all_genes) > 0) {
  sig_membership <- as.data.frame(
    sapply(method_names, function(m) as.integer(all_genes %in% sig_genes[[m]])),
    stringsAsFactors = FALSE
  )
  sig_membership$gene <- all_genes
} else {
  sig_membership <- data.frame(gene = character(0), stringsAsFactors = FALSE)
}

n_any_sig <- length(all_genes)
n_all_sig <- sum(rowSums(sig_membership[, method_names, drop = FALSE]) == length(method_names))
message(sprintf("  Agreement: %d genes sig in any method, %d in all %d methods",
                n_any_sig, n_all_sig, length(method_names)))
if (nrow(agreement) > 0) {
  message(sprintf("  Jaccard range: %.3f – %.3f",
                  min(agreement$jaccard, na.rm = TRUE),
                  max(agreement$jaccard, na.rm = TRUE)))
}

# 14. PERFORMANCE METRICS: power, FDR, type I error at both p<0.05 and FDR<0.05
#
# Since we know ground truth (is_de, beta_int_true), compute per method:
#   - power_nominal:   sensitivity at raw p < 0.05 among true DE genes
#   - power_fdr:       sensitivity at FDR < 0.05 among true DE genes
#   - typeI_nominal:   type I error rate at raw p < 0.05 among true null genes
#   - typeI_fdr:       type I error rate at FDR < 0.05 among true null genes
#   - FDP:             false discovery proportion among genes called sig at FDR < 0.05
#   - n_sig_nominal:   total genes significant at p < 0.05
#   - n_sig_fdr:       total genes significant at FDR < 0.05
#   - n_tested:        genes with non-NA p-values
# ─────────────────────────────────────────────────────────────────────────────
message("  Computing performance metrics...")

compute_metrics <- function(stats_df, method_name) {
  if (is.null(stats_df) || nrow(stats_df) == 0) {
    return(data.frame(
      method = method_name,
      n_tested = 0L, n_de_tested = 0L, n_null_tested = 0L,
      n_sig_nominal = 0L, n_sig_fdr = 0L,
      power_nominal = NA_real_, power_fdr = NA_real_,
      typeI_nominal = NA_real_, typeI_fdr = NA_real_,
      FDP = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  df <- stats_df[!is.na(stats_df$pval_int), ]
  n_tested <- nrow(df)
  
  # Split by ground truth
  de   <- df[!is.na(df$is_de) & df$is_de == TRUE, ]
  null <- df[!is.na(df$is_de) & df$is_de == FALSE, ]
  n_de   <- nrow(de)
  n_null <- nrow(null)
  
  # Counts
  n_sig_nominal <- sum(df$pval_int < 0.05, na.rm = TRUE)
  n_sig_fdr     <- sum(df$padj_int < 0.05, na.rm = TRUE)
  
  # Power (sensitivity among true DE)
  power_nominal <- if (n_de > 0) mean(de$pval_int < 0.05, na.rm = TRUE) else NA_real_
  power_fdr     <- if (n_de > 0) mean(de$padj_int < 0.05, na.rm = TRUE) else NA_real_
  
  # Type I error (among true nulls)
  typeI_nominal <- if (n_null > 0) mean(null$pval_int < 0.05, na.rm = TRUE) else NA_real_
  typeI_fdr     <- if (n_null > 0) mean(null$padj_int < 0.05, na.rm = TRUE) else NA_real_
  
  # False discovery proportion (among called positives at FDR < 0.05)
  sig_fdr_genes <- df[!is.na(df$padj_int) & df$padj_int < 0.05, ]
  FDP <- if (nrow(sig_fdr_genes) > 0) {
    mean(!sig_fdr_genes$is_de | is.na(sig_fdr_genes$is_de), na.rm = TRUE)
  } else {
    NA_real_
  }
  
  data.frame(
    method        = method_name,
    n_tested      = n_tested,
    n_de_tested   = n_de,
    n_null_tested = n_null,
    n_sig_nominal = n_sig_nominal,
    n_sig_fdr     = n_sig_fdr,
    power_nominal = round(power_nominal, 4),
    power_fdr     = round(power_fdr, 4),
    typeI_nominal = round(typeI_nominal, 4),
    typeI_fdr     = round(typeI_fdr, 4),
    FDP           = round(FDP, 4),
    stringsAsFactors = FALSE
  )
}

performance <- do.call(rbind, list(
  compute_metrics(nebula_stats, "NEBULA"),
  compute_metrics(edger_stats,  "edgeR"),
  compute_metrics(limma_stats,  "limma"),
  compute_metrics(wilcox_stats, "Wilcoxon"),
  compute_metrics(mast_stats,   "MAST")
))

message("  Performance metrics:")
print(performance)

# ─────────────────────────────────────────────────────────────────────────────
# 15. SAVE RESULTS TO S3 (small files only; counts already freed)
# ─────────────────────────────────────────────────────────────────────────────
timing <- c(sim    = t_sim_elapsed,
            nebula = t_neb_elapsed,
            edger  = t_er_elapsed,
            limma  = t_limma_elapsed,
            wilcox = t_wilcox_elapsed,
            mast   = t_mast_elapsed)

message(sprintf("  Saving results to s3://%s/%s ...", S3_BUCKET, S3_OUT))

s3saveRDS(nebula_stats,
          object = paste0(S3_OUT, "nebula_stats.rds"),
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
s3saveRDS(limma_stats,
          object = paste0(S3_OUT, "limma_stats.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(agreement,
          object = paste0(S3_OUT, "agreement.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(performance,
          object = paste0(S3_OUT, "performance.rds"),
          bucket = S3_BUCKET, region = "")
s3saveRDS(sig_membership,
          object = paste0(S3_OUT, "sig_membership.rds"),
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
  "══ [02] DONE  array=%d  total=%.1f s  (sim=%.1fs | neb=%.1fs | er=%.1fs | limma=%.1fs | wilcox=%.1fs | mast=%.1fs) ══",
  opt$array_id, total_elapsed,
  timing["sim.elapsed"], timing["nebula.elapsed"], timing["edger.elapsed"],
  timing["limma.elapsed"], timing["wilcox.elapsed"], timing["mast.elapsed"]
))
