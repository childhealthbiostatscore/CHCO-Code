################################################################################
# 05_volcano_single_sim.R
#
# PURPOSE:
#   Run ONE simulation using the reference model (with a specific parameter
#   combination) and produce volcano plots of the interaction term for every
#   DE method (NEBULA, edgeR, limma-voom+dupCor, Wilcoxon, MAST).
#
#   Two volcano plots per method:
#     1. Raw p-value on the y-axis
#     2. FDR-adjusted p-value on the y-axis
#
#   Output: a single multi-page PDF saved to S3.
#
# INPUT (from S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/
#     reference/hvg_genes.rds
#     reference/sce_ref.rds
#     reference/construct_data.rds
#
# OUTPUT (to S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/plots/
#     volcano_single_sim.pdf
#
# USAGE:
#   Rscript 05_volcano_single_sim.R \
#       --prop_de       0.15 \
#       --lfc_label     "med" \
#       --indiv_var     "no"  \
#       --cells_label   "low" \
#       --corr_cells    0.1   \
#       --n_subjects    10    \
#       --seed          42    \
#       --n_cores       4     \
#       --mast_max_cells 5000
################################################################################

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(nebula)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(BiocParallel)
  library(optparse)
  library(aws.s3)
  library(Seurat)
  library(MAST)
})

# =============================================================================
# S3 SETUP
# =============================================================================
setup_s3 <- function() {
  user <- Sys.info()[["user"]]

  if (user == "choiyej") {
    keys_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json"
  } else if (user %in% c("rameshsh", "yejichoi", "pylell")) {
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

S3_BUCKET   <- "scrna"
S3_BASE     <- "Projects/Paired scRNA simulation analysis/results/"
S3_REF_PFX  <- paste0(S3_BASE, "reference/")

# =============================================================================
# CLI ARGS
# =============================================================================
option_list <- list(
  make_option("--prop_de",        type = "double",    default = 0.15,
              help = "Proportion of DE genes [default: 0.15]"),
  make_option("--lfc_label",      type = "character", default = "med",
              help = "Effect size tier: vlow/low/med/high [default: med]"),
  make_option("--indiv_var",      type = "character", default = "no",
              help = "Individual variance: no/med/high [default: no]"),
  make_option("--cells_label",    type = "character", default = "low",
              help = "Cells per subject: low/high [default: low]"),
  make_option("--corr_cells",     type = "double",    default = 0.1,
              help = "Within-subject cell correlation [default: 0.1]"),
  make_option("--n_subjects",     type = "integer",   default = 10L,
              help = "Subjects per treatment arm [default: 10]"),
  make_option("--seed",           type = "integer",   default = 42L,
              help = "Random seed [default: 42]"),
  make_option("--n_cores",        type = "integer",   default = 4L,
              help = "Number of cores [default: 4]"),
  make_option("--nebula_method",  type = "character", default = "LN",
              help = "NEBULA method: LN or HL [default: LN]"),
  make_option("--pb_min_count",   type = "integer",   default = 10L,
              help = "Min rowSum for pseudobulk filtering [default: 10]"),
  make_option("--pb_min_samples", type = "integer",   default = 2L,
              help = "Min samples meeting pb_min_count [default: 2]"),
  make_option("--mast_max_cells", type = "integer",   default = 5000L,
              help = "Max POST cells for MAST [default: 5000]"),
  make_option("--output_name",    type = "character", default = "volcano_single_sim.pdf",
              help = "Output PDF filename [default: volcano_single_sim.pdf]")
)
opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)
BiocParallel::register(BiocParallel::MulticoreParam(opt$n_cores))

# Map labels to numeric values
lfc_tiers <- data.frame(
  lfc_label       = c("vlow", "low", "med", "high"),
  log2FC          = c(0.10,   0.25,  0.50,  1.00),
  interaction_lfc = c(0.10,   0.25,  0.50,  1.00) * log(2),
  stringsAsFactors = FALSE
)
iv_map      <- c(no = 1.0, med = 1.5, high = 2.5)
cl_mean_map <- c(low = 500L,  high = 6000L)
cl_sd_map   <- c(low = 150L,  high = 1500L)

interaction_lfc <- lfc_tiers$interaction_lfc[lfc_tiers$lfc_label == opt$lfc_label]
indiv_var_factor <- iv_map[opt$indiv_var]
cells_mean <- cl_mean_map[opt$cells_label]
cells_sd   <- cl_sd_map[opt$cells_label]

scenario_label <- paste0(
  "prop_de=", opt$prop_de,
  " | log2FC=", lfc_tiers$log2FC[lfc_tiers$lfc_label == opt$lfc_label],
  " (", opt$lfc_label, ")",
  " | indiv_var=", opt$indiv_var,
  " | cells=", opt$cells_label,
  " | corr=", opt$corr_cells,
  " | n_subj/arm=", opt$n_subjects
)
message(sprintf("== [05] Volcano plot: %s ==", scenario_label))

# =============================================================================
# 1. LOAD REFERENCE FROM S3
# =============================================================================
message("  Loading reference model from S3...")
hvg_genes <- s3readRDS(object = paste0(S3_REF_PFX, "hvg_genes.rds"),
                        bucket = S3_BUCKET, region = "")
sce_ref   <- s3readRDS(object = paste0(S3_REF_PFX, "sce_ref.rds"),
                        bucket = S3_BUCKET, region = "")
ref_data  <- s3readRDS(object = paste0(S3_REF_PFX, "construct_data.rds"),
                        bucket = S3_BUCKET, region = "")

n_genes <- length(hvg_genes)

# =============================================================================
# 2. BUILD CELL-LEVEL DESIGN
# =============================================================================
treatments <- c("Placebo", "Treatment")
visits     <- c("timepoint1", "timepoint2")

ref_subjects <- levels(ref_data$ref_construct_data$newCovariate$subject_id)

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
  lapply(seq_len(opt$n_subjects), function(j) {
    subj <- sprintf("%s_S%02d", prefix, j)
    make_subject_cells(subj, trt, cells_mean, cells_sd)
  })
})

new_covariate <- do.call(rbind, do.call(c, cell_list))
new_covariate$visit     <- factor(new_covariate$visit,
                                  levels = c("timepoint1", "timepoint2"))
new_covariate$treatment <- factor(new_covariate$treatment,
                                  levels = c("Placebo", "Treatment"))
new_covariate$celltype  <- unique(sce_ref$celltype)[1]

n_cells_total <- nrow(new_covariate)
message(sprintf("  Cells: %d  |  Subjects: %d (%d per arm x 2 arms)",
                n_cells_total, length(unique(new_covariate$subject_id)),
                opt$n_subjects))

# =============================================================================
# 3. GROUND TRUTH
# =============================================================================
n_de   <- round(opt$prop_de * n_genes)
de_idx <- if (n_de > 0) sample(seq_len(n_genes), n_de) else integer(0)

beta_int_true <- numeric(n_genes)
if (n_de > 0) {
  signs <- sample(c(-1L, 1L), n_de, replace = TRUE)
  beta_int_true[de_idx] <- signs * interaction_lfc
}

truth_df <- data.frame(
  gene          = hvg_genes,
  is_de         = seq_len(n_genes) %in% de_idx,
  beta_int_true = beta_int_true,
  log2FC_true   = beta_int_true / log(2),
  stringsAsFactors = FALSE
)

# =============================================================================
# 4. MODIFY scDesign3 MODEL
# =============================================================================
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

message("  Modifying scDesign3 model...")
modified_marginal <- bplapply(
  seq_along(ref_data$ref_fit_marginal),
  function(g) modify_gene_model(ref_data$ref_fit_marginal[[g]],
                                beta_int  = beta_int_true[g],
                                iv_factor = indiv_var_factor),
  BPPARAM = BiocParallel::MulticoreParam(opt$n_cores)
)
names(modified_marginal) <- names(ref_data$ref_fit_marginal)

# =============================================================================
# 5. SIMULATE COUNTS
# =============================================================================
message("  Extracting parameters...")
new_covariate_for_sim <- new_covariate[, colnames(ref_data$ref_construct_data$newCovariate)]
new_covariate_for_sim <- new_covariate_for_sim %>%
  dplyr::mutate(corr_group = subject_id)

sim_para <- extract_para(
  sce           = sce_ref,
  marginal_list = modified_marginal,
  n_cores       = opt$n_cores,
  family_use    = "nb",
  new_covariate = new_covariate_for_sim,
  data          = ref_data$ref_construct_data$dat
)

# Remap copula
ref_copula   <- ref_data$ref_fit_copula$copula_list
ref_names    <- names(ref_copula)
ref_dat_full <- ref_data$ref_construct_data$dat %>%
  dplyr::mutate(corr_group = as.character(corr_group))

new_subjects <- as.character(unique(new_covariate$subject_id))
ref_assignment <- setNames(
  ref_names[((seq_along(new_subjects) - 1) %% length(ref_names)) + 1],
  new_subjects
)

mapped_copula <- setNames(
  lapply(new_subjects, function(s) ref_copula[[ ref_assignment[[s]] ]]),
  new_subjects
)

ref_subject_ids <- unique(as.character(ref_dat_full$subject_id))
input_data_list <- lapply(new_subjects, function(s) {
  if (s %in% ref_subject_ids) {
    ref_dat_full %>% dplyr::filter(subject_id == s)
  } else {
    donor <- ref_assignment[[s]]
    ref_dat_full %>%
      dplyr::filter(subject_id == donor) %>%
      dplyr::mutate(subject_id = s, corr_group = s)
  }
})
ref_data_dat <- do.call(rbind, input_data_list)

message("  Simulating counts...")
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

counts <- sim_result
if (is.null(rownames(counts))) rownames(counts) <- hvg_genes

# Clean non-finite values
if (inherits(counts, "dgCMatrix")) {
  bad <- !is.finite(counts@x)
  if (any(bad)) {
    message(sprintf("  WARNING: %d non-finite values set to 0", sum(bad)))
    counts@x[bad] <- 0L
    counts <- drop0(counts)
  }
} else {
  bad <- !is.finite(counts)
  if (any(bad)) {
    message(sprintf("  WARNING: %d non-finite values set to 0", sum(bad)))
    counts[bad] <- 0L
  }
}

message(sprintf("  Simulated: %d genes x %d cells", nrow(counts), ncol(counts)))

rm(modified_marginal, sim_result, sim_para, ref_data, ref_data_dat, mapped_copula)
gc(verbose = FALSE)

# =============================================================================
# 6. HELPER: tidy results + join truth
# =============================================================================
tidy_join <- function(gene, logfc, pval) {
  df <- data.frame(gene      = gene,
                   logFC_int = logfc,
                   pval_int  = pval,
                   padj_int  = p.adjust(pval, method = "BH"),
                   stringsAsFactors = FALSE)
  left_join(df, truth_df, by = "gene")
}

# =============================================================================
# 7. RUN ALL 5 DE METHODS
# =============================================================================

# --- NEBULA ---
message("  Running NEBULA...")
run_nebula <- function(counts, new_covariate, method, n_cores) {
  nc <- new_covariate
  nc$visit     <- factor(nc$visit,      levels = c("timepoint1", "timepoint2"))
  nc$treatment <- factor(nc$treatment,  levels = c("Placebo",     "Treatment"))
  nc$subject_id <- as.character(nc$subject_id)

  X   <- model.matrix(~ visit * treatment, data = nc)
  neb <- group_cell(count = counts, id = nc$subject_id, pred = X)

  if (is.null(neb)) {
    neb <- list(count = counts, id = nc$subject_id, pred = X)
  }

  subjects_Treatment <- unique(nc$subject_id[nc$treatment == "Treatment"])
  subjects_Placebo <- unique(nc$subject_id[nc$treatment == "Placebo"])
  keep <- sapply(seq_len(nrow(neb$count)), function(g) {
    expr_vec <- neb$count[g, ]
    subA <- length(unique(neb$id[expr_vec > 0 & neb$id %in% subjects_Treatment]))
    subB <- length(unique(neb$id[expr_vec > 0 & neb$id %in% subjects_Placebo]))
    (subA >= 2) & (subB >= 2)
  })
  neb$count <- neb$count[keep, ]

  nebula(count = neb$count, id = neb$id, pred = neb$pred,
         method = method, ncore = n_cores, verbose = TRUE)
}

nebula_stats <- tryCatch({
  fit_neb <- run_nebula(counts, new_covariate, opt$nebula_method, opt$n_cores)
  summ     <- fit_neb$summary
  gene_vec <- if ("gene" %in% colnames(summ)) summ$gene else rownames(summ)
  lfc_col  <- grep("logFC.*visittimepoint2.*Treatment|logFC.*Treatment.*visittimepoint2",
                   colnames(summ), value = TRUE)
  if (!length(lfc_col)) lfc_col <- grep(":", colnames(summ), value = TRUE)[1]
  pval_col <- sub("^logFC_", "p_", lfc_col[1])
  tidy_join(gene  = gene_vec,
            logfc = summ[[lfc_col[1]]] / log(2),
            pval  = summ[[pval_col]])
}, error = function(e) {
  message("  NEBULA error: ", conditionMessage(e))
  NULL
})

# --- Pseudobulk aggregation ---
message("  Aggregating pseudobulk...")
new_covariate$pb_id <- paste(new_covariate$subject_id,
                             new_covariate$visit, sep = "__")
unique_pbs <- unique(new_covariate$pb_id)

pb_counts <- vapply(unique_pbs, function(pb) {
  idx <- which(new_covariate$pb_id == pb)
  rowSums(counts[, idx, drop = FALSE])
}, numeric(nrow(counts)))

pb_meta <- new_covariate[match(unique_pbs, new_covariate$pb_id),
                         c("pb_id", "subject_id", "visit", "treatment")]
rownames(pb_meta) <- unique_pbs
pb_meta$visit     <- factor(pb_meta$visit,
                            levels = c("timepoint1", "timepoint2"))
pb_meta$treatment <- factor(pb_meta$treatment,
                            levels = c("Placebo", "Treatment"))

keep           <- rowSums(pb_counts >= opt$pb_min_count) >= opt$pb_min_samples
pb_counts_filt <- pb_counts[keep, ]
pb_design <- model.matrix(~ visit + treatment + visit:treatment, data = pb_meta)

# --- limma-voom + duplicateCorrelation ---
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

  full <- data.frame(gene = hvg_genes, logFC_int = NA_real_,
                     pval_int = NA_real_, stringsAsFactors = FALSE)
  hit                 <- match(rownames(res_limma), full$gene)
  full$logFC_int[hit] <- res_limma$logFC
  full$pval_int[hit]  <- res_limma$P.Value
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")
}, error = function(e) {
  message("  limma error: ", conditionMessage(e))
  NULL
})

# --- edgeR ---
message("  Running edgeR...")
edger_stats <- tryCatch({
  dge <- DGEList(counts = pb_counts_filt)
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design = pb_design)
  fit_er  <- glmQLFit(dge, design = pb_design, robust = TRUE)
  int_col <- grep("visittimepoint2.*Treatment|Treatment.*visittimepoint2",
                  colnames(pb_design))
  if (!length(int_col)) int_col <- ncol(pb_design)
  qlf    <- glmQLFTest(fit_er, coef = int_col[1])
  res_er <- topTags(qlf, n = Inf, sort.by = "none")$table

  full <- data.frame(gene = hvg_genes, logFC_int = NA_real_,
                     pval_int = NA_real_, stringsAsFactors = FALSE)
  hit                 <- match(rownames(res_er), full$gene)
  full$logFC_int[hit] <- res_er$logFC
  full$pval_int[hit]  <- res_er$PValue
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")
}, error = function(e) {
  message("  edgeR error: ", conditionMessage(e))
  NULL
})

# --- Seurat object for cell-level methods ---
message("  Building Seurat object for cell-level methods...")
so_sim <- CreateSeuratObject(counts = counts, meta.data = new_covariate,
                             min.cells = 0, min.features = 0)
options(future.globals.maxSize = 2 * 1024^3)
so_sim <- NormalizeData(so_sim, verbose = FALSE)
so_post <- subset(so_sim, subset = visit == "timepoint2")
Idents(so_post) <- "treatment"

# --- Wilcoxon ---
message("  Running FindMarkers (Wilcoxon)...")
wilcox_stats <- tryCatch({
  fm_w <- FindMarkers(so_post, ident.1 = "Treatment", ident.2 = "Placebo",
                      test.use = "wilcox", logfc.threshold = 0,
                      min.pct = 0, verbose = FALSE)
  full <- data.frame(gene = hvg_genes, logFC_int = NA_real_,
                     pval_int = NA_real_, stringsAsFactors = FALSE)
  hit                 <- match(rownames(fm_w), full$gene)
  full$logFC_int[hit] <- fm_w$avg_log2FC
  full$pval_int[hit]  <- fm_w$p_val
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")
}, error = function(e) {
  message("  Wilcoxon error: ", conditionMessage(e))
  NULL
})

# --- MAST ---
message("  Running FindMarkers (MAST)...")
mast_stats <- tryCatch({
  so_mast <- so_post
  n_post <- ncol(so_post)
  if (n_post > opt$mast_max_cells) {
    n_per_group <- floor(opt$mast_max_cells / 2)
    cells_A <- WhichCells(so_post, idents = "Treatment")
    cells_B <- WhichCells(so_post, idents = "Placebo")
    keep_cells <- c(
      sample(cells_A, min(n_per_group, length(cells_A))),
      sample(cells_B, min(n_per_group, length(cells_B)))
    )
    so_mast <- subset(so_post, cells = keep_cells)
  }
  fm_m <- FindMarkers(so_mast, ident.1 = "Treatment", ident.2 = "Placebo",
                      test.use = "MAST", latent.vars = "subject_id",
                      logfc.threshold = 0, min.pct = 0, verbose = FALSE)
  full <- data.frame(gene = hvg_genes, logFC_int = NA_real_,
                     pval_int = NA_real_, stringsAsFactors = FALSE)
  hit                 <- match(rownames(fm_m), full$gene)
  full$logFC_int[hit] <- fm_m$avg_log2FC
  full$pval_int[hit]  <- fm_m$p_val
  full$padj_int       <- p.adjust(full$pval_int, method = "BH")
  left_join(full, truth_df, by = "gene")
}, error = function(e) {
  message("  MAST error: ", conditionMessage(e))
  NULL
})

rm(so_sim, so_post, counts)
if (exists("so_mast")) rm(so_mast)
gc(verbose = FALSE)

# =============================================================================
# 8. VOLCANO PLOT FUNCTION
# =============================================================================
make_volcano <- function(data, p_col, fc, title = NULL,
                         test_name,
                         sig_type = "pval", fdr_col = NULL,
                         p_thresh = 0.05,
                         cell_type = "",
                         formula_text = NULL,
                         cohort_text = NULL,
                         positive_text = "Positive",
                         negative_text = "Negative",
                         text_size = 20,
                         geom_text_size = 4,
                         caption_size = 8.5,
                         legend_text_size = 10,
                         volcano_force = 6,
                         volcano_box_padding = 0,
                         off_chart_threshold = 0.95,
                         off_chart_y_position = 0.85,
                         off_chart_arrow_length = 0.02) {

  if (!p_col %in% names(data) || !fc %in% names(data)) return(NULL)

  # Rename gene column if needed
  if (!"Gene" %in% names(data) && "gene" %in% names(data)) {
    data$Gene <- data$gene
  }

  data <- data %>% dplyr::filter(!is.na(.data[[p_col]]) & !is.na(.data[[fc]]))
  if (nrow(data) == 0) return(NULL)

  # Apply significance filter based on sig_type
  if (sig_type == "fdr") {
    if (is.null(fdr_col) || !fdr_col %in% names(data)) return(NULL)
    data <- data %>% dplyr::filter(!is.na(.data[[fdr_col]]))
    data$is_sig <- data[[fdr_col]] < p_thresh
  } else {
    data$is_sig <- data[[p_col]] < p_thresh
  }

  set.seed(1)
  epsilon <- 1e-300

  # For FDR plots, use adjusted p-values on the y-axis; otherwise use raw p-values
  y_col <- if (sig_type == "fdr") fdr_col else p_col
  data <- data %>%
    dplyr::mutate(neg_log_p = -log10(.data[[y_col]] + epsilon))

  y_max <- max(data$neg_log_p, na.rm = TRUE) * 1.1
  y_cutoff <- y_max * off_chart_threshold

  top_pos <- data %>%
    dplyr::filter(.data[[fc]] > 0 & is_sig) %>%
    dplyr::arrange(.data[[y_col]])
  n_pos <- nrow(top_pos)
  top_pos_n <- top_pos %>% dplyr::slice_head(n = 20)

  top_neg <- data %>%
    dplyr::filter(.data[[fc]] < 0 & is_sig) %>%
    dplyr::arrange(.data[[y_col]])
  n_neg <- nrow(top_neg)
  top_neg_n <- top_neg %>% dplyr::slice_head(n = 20)

  # Identify off-chart genes
  off_chart_genes <- data %>%
    dplyr::filter(Gene %in% c(top_pos_n$Gene, top_neg_n$Gene) & neg_log_p > y_cutoff) %>%
    dplyr::mutate(
      is_positive = .data[[fc]] > 0,
      x_position = ifelse(is_positive,
                          .data[[fc]] + seq(from = 0.1, by = 0.2, length.out = n()),
                          .data[[fc]] - seq(from = 0.1, by = 0.2, length.out = n())),
      y_position = y_max * off_chart_y_position
    )

  on_chart_genes <- c(top_pos_n$Gene, top_neg_n$Gene)[!c(top_pos_n$Gene, top_neg_n$Gene) %in% off_chart_genes$Gene]

  data <- data %>%
    dplyr::mutate(
      top_color = case_when(
        Gene %in% top_pos$Gene ~ "#f28482",
        Gene %in% top_neg$Gene ~ "#457b9d",
        TRUE ~ "#ced4da"
      ),
      top_size = if_else(Gene %in% c(top_pos$Gene, top_neg$Gene), 1.3, 1),
      top_lab  = if_else(Gene %in% on_chart_genes, Gene, ""),
      display_neg_log_p = pmin(neg_log_p, y_cutoff)
    ) %>%
    dplyr::filter(abs(.data[[fc]]) < 10)

  max_fc <- max(data[[fc]], na.rm = TRUE)
  min_fc <- min(data[[fc]], na.rm = TRUE)

  # Build caption
  has_formula <- !is.null(formula_text) && !is.na(formula_text)
  has_cohort  <- !is.null(cohort_text) && !is.na(cohort_text)

  caption_text <- paste0(
    if (has_formula && has_cohort) {
      paste0("Formula: ", formula_text, " + (1|subject) | Cohort: ", cohort_text, "\n\n")
    } else if (has_formula) {
      paste0("Formula: ", formula_text, " + (1|subject)\n\n")
    } else if (has_cohort) {
      paste0("Cohort: ", cohort_text, "\n\n")
    } else { "" },
    "Cell type: ", cell_type,
    if (cell_type != "") " | " else "",
    "Positive n = ", n_pos, " | Negative n = ", n_neg,
    if (nrow(off_chart_genes) > 0) paste0("\n", nrow(off_chart_genes), " gene(s) with p-values near zero (arrows indicate off-scale values)") else ""
  )

  # Y-axis label
  y_label <- ifelse(sig_type == "fdr", "-log10(adj. p-value)", "-log10(p-value)")

  p <- ggplot(data, aes(x = .data[[fc]], y = display_neg_log_p)) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.5, aes(color = top_color), size = 3) +
    geom_label_repel(seed = 1,
                     data = dplyr::filter(data, top_lab != ""),
                     aes(label = top_lab, color = top_color),
                     fontface = "bold",
                     size = geom_text_size, max.overlaps = Inf,
                     force = volcano_force, segment.alpha = 0.3, segment.size = 0.3,
                     box.padding = volcano_box_padding,
                     fill = fill_alpha("white", 0.7),
                     label.size = 0
    )

  # Add arrows for off-chart genes
  if (nrow(off_chart_genes) > 0) {
    p <- p +
      geom_segment(
        data = off_chart_genes,
        aes(x = .data[[fc]], y = y_cutoff * 0.98,
            xend = .data[[fc]], yend = y_cutoff - (y_max * off_chart_arrow_length)),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        color = "black", size = 0.6
      ) +
      geom_text(
        data = off_chart_genes,
        aes(x = x_position, y = y_position, label = Gene),
        size = geom_text_size * 0.9,
        hjust = ifelse(off_chart_genes$is_positive, 0, 1),
        color = "black", fontface = "italic"
      ) +
      geom_text(
        data = off_chart_genes,
        aes(x = x_position, y = y_position - (y_max * 0.03),
            label = paste0("p=", format(.data[[y_col]], scientific = TRUE, digits = 2))),
        size = geom_text_size * 0.7,
        hjust = ifelse(off_chart_genes$is_positive, 0, 1),
        color = "darkgrey"
      )
  }

  p <- p +
    labs(title = title,
         x = paste0("log2 FC ", test_name),
         y = y_label,
         caption = caption_text) +
    scale_size_continuous(range = c(1, 1.3)) +
    scale_color_manual(values = c("#457b9d" = "#457b9d", "#ced4da" = "#ced4da", "#f28482" = "#f28482")) +
    guides(color = "none", size = "none") +
    annotate("segment",
             x = max_fc / 8, xend = (max_fc * 7) / 8,
             y = -y_max * 0.13,
             col = "darkgrey", arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text",
             x = mean(c(max_fc / 8, (max_fc * 7) / 8)),
             y = -y_max * 0.18,
             label = positive_text,
             size = geom_text_size, color = "#343a40") +
    annotate("segment",
             x = min_fc / 8, xend = (min_fc * 7) / 8,
             y = -y_max * 0.13,
             col = "darkgrey", arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text",
             x = mean(c(min_fc / 8, (min_fc * 7) / 8)),
             y = -y_max * 0.18,
             label = negative_text,
             size = geom_text_size, color = "#343a40") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, y_max), clip = "off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text_size),
          title = element_text(size = legend_text_size),
          plot.margin = margin(t = 10, r = 20, b = 25, l = 20),
          axis.title.x = element_text(margin = margin(t = 32)),
          plot.caption = element_text(size = caption_size, hjust = 0.5, margin = margin(t = 15)),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background  = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", color = NA))

  return(p)
}

# =============================================================================
# 9. GENERATE VOLCANO PLOTS
# =============================================================================
message("  Generating volcano plots...")

methods <- list(
  list(name = "NEBULA",    data = nebula_stats, test_name = "(Interaction)",
       formula_text = "~ visit * treatment", note = "Cell-level mixed model"),
  list(name = "limma",     data = limma_stats,  test_name = "(Interaction)",
       formula_text = "~ visit * treatment + dupCor(subject)", note = "Pseudobulk limma-voom + dupCor"),
  list(name = "edgeR",     data = edger_stats,  test_name = "(Interaction)",
       formula_text = "~ visit + treatment + visit:treatment", note = "Pseudobulk QL F-test"),
  list(name = "Wilcoxon",  data = wilcox_stats, test_name = "(Trt POST vs Plc POST)",
       formula_text = NULL, note = "Cell-level naive; tests b_trt + b_int"),
  list(name = "MAST",      data = mast_stats,   test_name = "(Trt POST vs Plc POST)",
       formula_text = NULL, note = "Cell-level hurdle; subject as fixed covariate")
)

plot_list <- list()

for (m in methods) {
  if (is.null(m$data)) {
    message(sprintf("  Skipping %s (NULL result)", m$name))
    next
  }

  # Raw p-value volcano
  p_pval <- make_volcano(
    data      = m$data,
    p_col     = "pval_int",
    fc        = "logFC_int",
    title     = paste0(m$name, " — Raw p-value"),
    test_name = m$test_name,
    sig_type  = "pval",
    fdr_col   = NULL,
    p_thresh  = 0.05,
    cell_type = unique(sce_ref$celltype)[1],
    formula_text = m$formula_text,
    cohort_text  = NULL,
    positive_text = "",
    negative_text = ""
  )

  # FDR volcano
  p_fdr <- make_volcano(
    data      = m$data,
    p_col     = "pval_int",
    fc        = "logFC_int",
    title     = paste0(m$name, " — FDR-adjusted p-value"),
    test_name = m$test_name,
    sig_type  = "fdr",
    fdr_col   = "padj_int",
    p_thresh  = 0.05,
    cell_type = unique(sce_ref$celltype)[1],
    formula_text = m$formula_text,
    cohort_text  = NULL,
    positive_text = "",
    negative_text = ""
  )

  if (!is.null(p_pval)) plot_list[[paste0(m$name, "_pval")]] <- p_pval
  if (!is.null(p_fdr))  plot_list[[paste0(m$name, "_fdr")]]  <- p_fdr
}

# =============================================================================
# 10. SAVE TO PDF AND UPLOAD TO S3
# =============================================================================
message("  Saving volcano plots to PDF...")

tmp_pdf <- tempfile(fileext = ".pdf")
pdf(tmp_pdf, width = 14, height = 10)

for (pname in names(plot_list)) {
  # Add a super-title with the scenario parameters
  p <- plot_list[[pname]] +
    labs(subtitle = scenario_label)
  print(p)
}

# Also create a combined side-by-side page per method (pval | fdr)
for (m in methods) {
  pval_key <- paste0(m$name, "_pval")
  fdr_key  <- paste0(m$name, "_fdr")
  if (pval_key %in% names(plot_list) && fdr_key %in% names(plot_list)) {
    combined <- plot_list[[pval_key]] + plot_list[[fdr_key]] +
      plot_annotation(
        title    = paste0(m$name, " — Interaction Volcano Plots"),
        subtitle = scenario_label,
        theme    = theme(
          plot.title    = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey40")
        )
      )
    print(combined)
  }
}

dev.off()

# Upload to S3
s3_output_path <- paste0(S3_BASE, "plots/", opt$output_name)
aws.s3::put_object(
  file   = tmp_pdf,
  object = s3_output_path,
  bucket = S3_BUCKET,
  region = ""
)
unlink(tmp_pdf)

message(sprintf("== [05] Done. Volcano plots saved to s3://%s/%s ==",
                S3_BUCKET, s3_output_path))
