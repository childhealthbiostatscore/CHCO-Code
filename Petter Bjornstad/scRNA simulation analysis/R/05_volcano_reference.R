################################################################################
# 05_volcano_reference.R
#
# PURPOSE:
#   Run all 5 DE methods (NEBULA, edgeR, limma-voom+dupCor, Wilcoxon, MAST) on the
#   ACTUAL ATTEMPT reference dataset (not simulated data) and produce volcano
#   plots of the interaction term for each method.
#
#   Two volcano plots per method:
#     1. Raw p-value on the y-axis
#     2. FDR-adjusted p-value on the y-axis
#
#   Since this is real data, there is no ground truth — the volcano plots show
#   observed differential expression from the actual clinical trial.
#
# INPUT (from S3):
#   bucket: attempt
#     cleaned_data/attempt_clean_so.rds   -- ATTEMPT Seurat object
#
# OUTPUT (to S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/plots/
#     volcano_reference.pdf
#
# USAGE:
#   Rscript 05_volcano_reference.R \
#       --cell_type     "PT"   \
#       --n_hvg         2000   \
#       --n_cores       4      \
#       --mast_max_cells 5000  \
#       --output_name   "volcano_reference.pdf"
################################################################################

suppressPackageStartupMessages({
  library(nebula)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(gridExtra)
  library(UpSetR)
  library(BiocParallel)
  library(optparse)
  library(aws.s3)
  library(Seurat)
  library(MAST)
  library(scran)
  library(SingleCellExperiment)
  library(foreach)
  library(doParallel)
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

S3_BUCKET <- "scrna"
S3_BASE   <- "Projects/Paired scRNA simulation analysis/results/"

# =============================================================================
# CLI ARGS
# =============================================================================
option_list <- list(
  make_option("--cell_type",      type = "character", default = "PT",
              help = "KPMP_celltype_general to subset [default: PT]"),
  make_option("--n_hvg",          type = "integer",   default = 2000L,
              help = "Number of HVGs [default: 2000]"),
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
  make_option("--seed",           type = "integer",   default = 42L,
              help = "Random seed [default: 42]"),
  make_option("--output_name",    type = "character", default = "volcano_reference.pdf",
              help = "Output PDF filename [default: volcano_reference.pdf]")
)
opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)
BiocParallel::register(BiocParallel::MulticoreParam(opt$n_cores))

message(sprintf("== [05-ref] Volcano plots on ATTEMPT reference data | cell_type=%s ==",
                opt$cell_type))

# =============================================================================
# 1. LOAD ATTEMPT SEURAT OBJECT FROM S3
# =============================================================================
message("  Loading ATTEMPT Seurat object from S3 (bucket: attempt)...")
so <- s3readRDS(
  object = "cleaned_data/attempt_clean_so.rds",
  bucket = "attempt",
  region = ""
)
message(sprintf("  Loaded: %d genes x %d cells", nrow(so), ncol(so)))

# =============================================================================
# 2. SUBSET TO CELL TYPE
# =============================================================================
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

cells_keep <- which(so$KPMP_celltype_general == opt$cell_type)
if (length(cells_keep) == 0) {
  stop(sprintf("No cells found with KPMP_celltype_general == '%s'. ",
               opt$cell_type),
       "Available: ", paste(unique(so$KPMP_celltype_general), collapse = ", "))
}
so_sub <- so[, cells_keep]
rm(so); gc(verbose = FALSE)
message(sprintf("  Cell type '%s': %d cells", opt$cell_type, ncol(so_sub)))

# =============================================================================
# 3. SELECT HVGs
# =============================================================================
message("  Selecting top HVGs...")
so_sub <- FindVariableFeatures(so_sub, selection.method = "vst",
                               nfeatures = opt$n_hvg, verbose = FALSE)
hvg_genes <- VariableFeatures(so_sub)[seq_len(opt$n_hvg)]
message(sprintf("  Selected %d HVGs", length(hvg_genes)))

# =============================================================================
# 4. RECODE METADATA
# =============================================================================
# Recode to generic labels matching simulation pipeline
#   visit:     PRE -> timepoint1, POST -> timepoint2
#   treatment: Placebo -> Placebo, Dapagliflozin -> Treatment
# Clean subject IDs: grpPlacebo_S01, grpTreatment_S01, etc.

orig_ids <- sort(unique(as.character(so_sub$subject_id)))
orig_trt <- sapply(orig_ids, function(id) {
  unique(as.character(so_sub$treatment[so_sub$subject_id == id]))
})

placebo_ids   <- orig_ids[orig_trt == "Placebo"]
treatment_ids <- orig_ids[orig_trt != "Placebo"]

clean_ids <- character(length(orig_ids))
names(clean_ids) <- orig_ids
clean_ids[placebo_ids]   <- sprintf("grpPlacebo_S%02d", seq_along(placebo_ids))
clean_ids[treatment_ids] <- sprintf("grpTreatment_S%02d", seq_along(treatment_ids))

so_sub$subject_id_clean <- unname(clean_ids[as.character(so_sub$subject_id)])
so_sub$visit_recode     <- factor(
  ifelse(so_sub$visit == "PRE", "timepoint1", "timepoint2"),
  levels = c("timepoint1", "timepoint2")
)
so_sub$treatment_recode <- factor(
  ifelse(so_sub$treatment == "Placebo", "Placebo", "Treatment"),
  levels = c("Placebo", "Treatment")
)

n_subjects <- length(unique(so_sub$subject_id_clean))
n_per_arm  <- table(orig_trt)
message(sprintf("  Subjects: %d total (Placebo=%d, Treatment=%d)",
                n_subjects, n_per_arm["Placebo"],
                n_per_arm[names(n_per_arm) != "Placebo"]))

# Build a metadata data.frame for DE methods
meta <- data.frame(
  cell        = colnames(so_sub),
  subject_id  = so_sub$subject_id_clean,
  visit       = so_sub$visit_recode,
  treatment   = so_sub$treatment_recode,
  celltype    = so_sub$KPMP_celltype_general,
  stringsAsFactors = FALSE
)

# Get raw counts for the HVGs
counts <- GetAssayData(so_sub, layer = "counts")[hvg_genes, ]

# =============================================================================
# 5. HELPER: tidy results (no ground truth for real data)
# =============================================================================
tidy_results <- function(gene, logfc, pval) {
  data.frame(
    gene      = gene,
    logFC_int = logfc,
    pval_int  = pval,
    padj_int  = p.adjust(pval, method = "BH"),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 6. NEBULA
# =============================================================================
t_neb <- proc.time()
message("  Running NEBULA...")

group_cell_mod <- function (count, id, pred = NULL, offset = NULL) 
{
  ng = nrow(count)
  nc = ncol(count)
  if (nc != length(id)) {
    stop("The length of id is not equal to the number of columns of the count matrix.")
  }
  id = as.character(id)
  levels = unique(id)
  id = factor(id, levels = levels)
  if (is.unsorted(id) == FALSE) {
    cat("The cells are already grouped.")
    return(NULL)
  }
  k = length(levels)
  o = order(id)
  count = count[, o]
  id = id[o]
  if (is.null(pred) == FALSE) {
    if (nc != nrow(as.matrix(pred))) {
      stop("The number of rows of the design matrix is not equal to the number of columns of the count matrix")
    }
    pred = as.matrix(pred)[o, , drop = FALSE]
  }
  if (is.null(offset) == FALSE) {
    if (nc != length(offset)) {
      stop("The length of offset is not equal to the number of columns of the count matrix")
    }
    if (sum(offset <= 0) > 0) {
      stop("Some elements in the scaling factor are not positive.")
    }
    offset = offset[o]
  }
  grouped = list(count = count, id = id, pred = pred, offset = offset)
  return(grouped)
}

nebula_stats <- tryCatch({
  nc <- meta
  nc$visit     <- factor(nc$visit,      levels = c("timepoint1", "timepoint2"))
  nc$treatment <- factor(nc$treatment,  levels = c("Placebo",     "Treatment"))
  nc$subject_id <- as.character(nc$subject_id)
  counts_rounded <- round(counts)
  genes_list <- rownames(counts_rounded)

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
  
  # Gene-by-gene parallel loop (matches ATTEMPT QMD approach)
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
      pred_gene  <- model.matrix(~ treatment * visit, data = nc)
      data_gene  <- group_cell_mod(count  = count_gene,
                                   id     = nc$subject_id,
                                   pred   = pred_gene,
                                   offset = pooled_offset)

      # If group_cell_mod returns NULL (cells already sorted), use raw inputs
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
          model      = "NBLMM",
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
    # Find the interaction column (treatment * visit formula order)
    lfc_col <- grep("logFC.*treatment.*visit|logFC.*Treatment.*timepoint",
                    colnames(summ), value = TRUE, ignore.case = TRUE)
    if (!length(lfc_col)) {
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
  df <- do.call(rbind, parsed)
  df$padj_int <- p.adjust(df$pval_int, method = "BH")
  df

}, error = function(e) {
  message("  NEBULA error: ", conditionMessage(e))
  NULL
})

t_neb_elapsed <- (proc.time() - t_neb)["elapsed"]
message(sprintf("  NEBULA done: %.1f s  (%s genes)",
                t_neb_elapsed, if (!is.null(nebula_stats)) nrow(nebula_stats) else "0"))

# =============================================================================
# 7. PSEUDOBULK AGGREGATION (shared by edgeR and limma)
# =============================================================================
message("  Aggregating pseudobulk...")

meta$pb_id <- paste(meta$subject_id, meta$visit, sep = "__")
unique_pbs <- unique(meta$pb_id)

pb_counts <- vapply(unique_pbs, function(pb) {
  idx <- which(meta$pb_id == pb)
  rowSums(counts[, idx, drop = FALSE])
}, numeric(nrow(counts)))

pb_meta <- meta[match(unique_pbs, meta$pb_id),
                c("pb_id", "subject_id", "visit", "treatment")]
rownames(pb_meta) <- unique_pbs
pb_meta$visit     <- factor(pb_meta$visit,
                            levels = c("timepoint1", "timepoint2"))
pb_meta$treatment <- factor(pb_meta$treatment,
                            levels = c("Placebo", "Treatment"))

pb_counts      <- round(pb_counts)
keep           <- rowSums(pb_counts >= opt$pb_min_count) >= opt$pb_min_samples
pb_counts_filt <- pb_counts[keep, ]
pb_design      <- model.matrix(~ visit + treatment + visit:treatment, data = pb_meta)

message(sprintf("  Pseudobulk: %d samples | %d/%d genes pass filter",
                ncol(pb_counts), sum(keep), nrow(pb_counts)))

# =============================================================================
# 8. edgeR
# =============================================================================
t_er <- proc.time()
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
  full
}, error = function(e) {
  message("  edgeR error: ", conditionMessage(e))
  NULL
})

t_er_elapsed <- (proc.time() - t_er)["elapsed"]
message(sprintf("  edgeR done: %.1f s", t_er_elapsed))

# =============================================================================
# 9. limma-voom + duplicateCorrelation
#
# Uses voom to transform pseudobulk counts, then estimates the within-subject
# correlation via duplicateCorrelation() with subject_id as the blocking factor.
# Design: ~ visit * treatment (no subject dummies — correlation handles pairing).
# =============================================================================
t_limma <- proc.time()
message("  Running limma-voom + duplicateCorrelation...")

limma_stats <- tryCatch({
  dge_limma <- DGEList(counts = pb_counts_filt)
  dge_limma <- calcNormFactors(dge_limma)

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
  full
}, error = function(e) {
  message("  limma error: ", conditionMessage(e))
  NULL
})

t_limma_elapsed <- (proc.time() - t_limma)["elapsed"]
message(sprintf("  limma done: %.1f s", t_limma_elapsed))

# =============================================================================
# 10. BUILD SEURAT OBJECT FOR CELL-LEVEL METHODS
# =============================================================================
message("  Building Seurat object for cell-level methods...")

# Need a Seurat object with the recoded metadata
so_de <- CreateSeuratObject(
  counts    = counts,
  min.cells = 0, min.features = 0
)
so_de$subject_id <- meta$subject_id
so_de$visit      <- meta$visit
so_de$treatment  <- meta$treatment

options(future.globals.maxSize = 2 * 1024^3)
so_de <- NormalizeData(so_de, verbose = FALSE)

# Subset to POST cells; set identity to treatment
so_post <- subset(so_de, subset = visit == "timepoint2")
Idents(so_post) <- "treatment"
n_post <- ncol(so_post)
message(sprintf("  POST cells: %d  (Treatment=%d, Placebo=%d)",
                n_post,
                sum(so_post$treatment == "Treatment"),
                sum(so_post$treatment == "Placebo")))

# =============================================================================
# 11. WILCOXON (FindMarkers — naive, ignores paired structure)
# =============================================================================
t_wilcox <- proc.time()
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
  full
}, error = function(e) {
  message("  Wilcoxon error: ", conditionMessage(e))
  NULL
})

t_wilcox_elapsed <- (proc.time() - t_wilcox)["elapsed"]
message(sprintf("  Wilcoxon done: %.1f s", t_wilcox_elapsed))

# =============================================================================
# 12. MAST (hurdle model; subject_id as fixed latent covariate)
# =============================================================================
t_mast <- proc.time()
message(sprintf("  Running FindMarkers (MAST, max_cells=%d)...", opt$mast_max_cells))

mast_stats <- tryCatch({
  so_mast <- so_post
  
  if (n_post > opt$mast_max_cells) {
    n_per_group <- floor(opt$mast_max_cells / 2)
    cells_A <- WhichCells(so_post, idents = "Treatment")
    cells_B <- WhichCells(so_post, idents = "Placebo")
    keep_cells <- c(
      sample(cells_A, min(n_per_group, length(cells_A))),
      sample(cells_B, min(n_per_group, length(cells_B)))
    )
    so_mast <- subset(so_post, cells = keep_cells)
    message(sprintf("  MAST: subsampled to %d POST cells (%d/group)",
                    ncol(so_mast), n_per_group))
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
  full
}, error = function(e) {
  message("  MAST error: ", conditionMessage(e))
  NULL
})

t_mast_elapsed <- (proc.time() - t_mast)["elapsed"]
message(sprintf("  MAST done: %.1f s", t_mast_elapsed))

rm(so_de, so_post)
if (exists("so_mast")) rm(so_mast)
gc(verbose = FALSE)

# =============================================================================
# 13. TIMING SUMMARY
# =============================================================================
timing <- c(nebula = t_neb_elapsed,
            edger  = t_er_elapsed,
            limma  = t_limma_elapsed,
            wilcox = t_wilcox_elapsed,
            mast   = t_mast_elapsed)
message(sprintf("  Timing: %s",
                paste(sprintf("%s=%.1fs", names(timing), timing), collapse = " | ")))

# =============================================================================
# 14. VOLCANO PLOT FUNCTION
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
  
  on_chart_genes <- c(top_pos_n$Gene, top_neg_n$Gene)[
    !c(top_pos_n$Gene, top_neg_n$Gene) %in% off_chart_genes$Gene
  ]
  
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
    if (nrow(off_chart_genes) > 0) {
      paste0("\n", nrow(off_chart_genes),
             " gene(s) with p-values near zero (arrows indicate off-scale values)")
    } else ""
  )
  
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
    scale_color_manual(values = c("#457b9d" = "#457b9d", "#ced4da" = "#ced4da",
                                  "#f28482" = "#f28482")) +
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
          plot.caption = element_text(size = caption_size, hjust = 0.5,
                                      margin = margin(t = 15)),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background  = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", color = NA))
  
  return(p)
}

# =============================================================================
# 15. METHOD LIST + SHARED THEME
# =============================================================================
message("  Preparing analyses...")

methods <- list(
  list(name = "NEBULA",    data = nebula_stats, test_name = "(Interaction)",
       formula_text = "~ visit * treatment", note = "Cell-level NBLMM + REML + pooled offset"),
  list(name = "edgeR",     data = edger_stats,  test_name = "(Interaction)",
       formula_text = "~ visit + treatment + visit:treatment", note = "Pseudobulk QL F-test"),
  list(name = "limma",     data = limma_stats,  test_name = "(Interaction)",
       formula_text = "~ visit * treatment + duplicateCorrelation", note = "Pseudobulk voom + dupCor"),
  list(name = "Wilcoxon",  data = wilcox_stats, test_name = "(Trt POST vs Plc POST)",
       formula_text = NULL, note = "Cell-level naive; tests b_trt + b_int"),
  list(name = "MAST",      data = mast_stats,   test_name = "(Trt POST vs Plc POST)",
       formula_text = NULL, note = "Cell-level hurdle; subject as fixed covariate")
)

# Keep only methods that ran successfully
methods_ok <- Filter(function(m) !is.null(m$data), methods)
method_names <- sapply(methods_ok, `[[`, "name")

method_colors <- c(NEBULA = "#e63946", edgeR = "#2a9d8f", limma = "#457b9d",
                   Wilcoxon = "#e9c46a", MAST = "#264653")

subtitle_text <- sprintf("ATTEMPT reference | Cell type: %s | %d HVGs | %d subjects",
                         opt$cell_type, opt$n_hvg, n_subjects)

shared_theme <- theme_minimal(base_size = 14) +
  theme(plot.title    = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        panel.grid.minor = element_blank())

# =============================================================================
# 16. VOLCANO PLOTS
# =============================================================================
message("  Generating volcano plots...")

volcano_list <- list()

for (m in methods_ok) {
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
    cell_type = opt$cell_type,
    formula_text = m$formula_text,
    cohort_text  = paste0("ATTEMPT (real data) | ", m$note),
    positive_text = "Up in Treatment",
    negative_text = "Down in Treatment"
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
    cell_type = opt$cell_type,
    formula_text = m$formula_text,
    cohort_text  = paste0("ATTEMPT (real data) | ", m$note),
    positive_text = "Up in Treatment",
    negative_text = "Down in Treatment"
  )
  
  if (!is.null(p_pval)) volcano_list[[paste0(m$name, "_pval")]] <- p_pval
  if (!is.null(p_fdr))  volcano_list[[paste0(m$name, "_fdr")]]  <- p_fdr
}

# =============================================================================
# 17. P-VALUE DISTRIBUTION HISTOGRAMS
# =============================================================================
message("  Generating p-value distribution histograms...")

pval_hist_list <- list()

for (m in methods_ok) {
  df <- m$data %>% dplyr::filter(!is.na(pval_int))
  if (nrow(df) == 0) next
  
  p <- ggplot(df, aes(x = pval_int)) +
    geom_histogram(bins = 50, fill = method_colors[m$name],
                   color = "white", alpha = 0.8) +
    geom_hline(yintercept = nrow(df) / 50, linetype = "dashed",
               color = "red", linewidth = 0.6) +
    labs(title = paste0(m$name, " — Raw p-value distribution"),
         subtitle = paste0(m$note, "\n",
                           "Dashed red line = expected count under uniform (null)"),
         x = "p-value", y = "Count") +
    shared_theme +
    theme(plot.subtitle = element_text(size = 9, color = "grey50"))
  pval_hist_list[[m$name]] <- p
}

# Combined overlay histogram
pval_combined_df <- do.call(rbind, lapply(methods_ok, function(m) {
  m$data %>%
    dplyr::filter(!is.na(pval_int)) %>%
    dplyr::mutate(method = m$name) %>%
    dplyr::select(gene, pval_int, method)
}))

pval_overlay <- ggplot(pval_combined_df, aes(x = pval_int, fill = method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity", color = NA) +
  scale_fill_manual(values = method_colors) +
  labs(title = "P-value distributions — All methods overlaid",
       subtitle = subtitle_text,
       x = "p-value", y = "Count", fill = "Method") +
  shared_theme

# Faceted version (cleaner when methods overlap a lot)
pval_faceted <- ggplot(pval_combined_df, aes(x = pval_int, fill = method)) +
  geom_histogram(bins = 50, alpha = 0.8, color = "white") +
  facet_wrap(~ method, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = method_colors) +
  guides(fill = "none") +
  labs(title = "P-value distributions — By method",
       subtitle = subtitle_text,
       x = "p-value", y = "Count") +
  shared_theme

# =============================================================================
# 18. METHOD AGREEMENT: UPSET PLOT + RANK CORRELATIONS
# =============================================================================
message("  Generating method agreement analyses...")

# --- Significant gene sets (FDR < 0.05) ---
sig_genes_list <- lapply(methods_ok, function(m) {
  m$data %>%
    dplyr::filter(!is.na(padj_int) & padj_int < 0.05) %>%
    dplyr::pull(gene)
})
names(sig_genes_list) <- method_names

# UpSet plot data: binary membership matrix
all_genes_tested <- unique(unlist(lapply(methods_ok, function(m) {
  m$data %>% dplyr::filter(!is.na(pval_int)) %>% dplyr::pull(gene)
})))
upset_df <- as.data.frame(sapply(sig_genes_list, function(gs) {
  as.integer(all_genes_tested %in% gs)
}))
rownames(upset_df) <- all_genes_tested

# --- Pairwise rank correlations (p-value ranks and logFC) ---
# Build a merged data.frame of all method results
merged_all <- NULL
for (m in methods_ok) {
  tmp <- m$data %>%
    dplyr::filter(!is.na(pval_int) & !is.na(logFC_int)) %>%
    dplyr::select(gene, logFC_int, pval_int) %>%
    dplyr::rename(
      !!paste0("logFC_", m$name)  := logFC_int,
      !!paste0("pval_", m$name)   := pval_int
    )
  if (is.null(merged_all)) {
    merged_all <- tmp
  } else {
    merged_all <- dplyr::inner_join(merged_all, tmp, by = "gene")
  }
}

# Compute Spearman correlations for logFC and p-value ranks
n_meth <- length(method_names)
cor_logfc <- matrix(NA, n_meth, n_meth, dimnames = list(method_names, method_names))
cor_pval  <- matrix(NA, n_meth, n_meth, dimnames = list(method_names, method_names))

for (i in seq_len(n_meth)) {
  for (j in seq_len(n_meth)) {
    fc_i <- merged_all[[paste0("logFC_", method_names[i])]]
    fc_j <- merged_all[[paste0("logFC_", method_names[j])]]
    pv_i <- merged_all[[paste0("pval_", method_names[i])]]
    pv_j <- merged_all[[paste0("pval_", method_names[j])]]
    cor_logfc[i, j] <- cor(fc_i, fc_j, method = "spearman", use = "complete.obs")
    cor_pval[i, j]  <- cor(pv_i, pv_j, method = "spearman", use = "complete.obs")
  }
}

# Heatmap helper
make_cor_heatmap <- function(cor_mat, title_text) {
  df <- as.data.frame(as.table(cor_mat))
  names(df) <- c("Method1", "Method2", "Correlation")
  ggplot(df, aes(x = Method1, y = Method2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4, color = "black") +
    scale_fill_gradient2(low = "#457b9d", mid = "white", high = "#e63946",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = title_text, x = NULL, y = NULL) +
    shared_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

cor_logfc_plot <- make_cor_heatmap(cor_logfc,
                                   "Spearman correlation of log2FC between methods")
cor_pval_plot <- make_cor_heatmap(cor_pval,
                                  "Spearman correlation of p-value ranks between methods")

# =============================================================================
# 19. LOG2FC SCATTER PLOTS BETWEEN METHOD PAIRS
# =============================================================================
message("  Generating log2FC comparison scatter plots...")

fc_scatter_list <- list()

method_pairs <- combn(method_names, 2, simplify = FALSE)

for (pair in method_pairs) {
  m1 <- pair[1]; m2 <- pair[2]
  fc1_col <- paste0("logFC_", m1)
  fc2_col <- paste0("logFC_", m2)
  
  if (!fc1_col %in% names(merged_all) || !fc2_col %in% names(merged_all)) next
  
  df_pair <- merged_all %>%
    dplyr::select(gene, all_of(c(fc1_col, fc2_col))) %>%
    dplyr::filter(!is.na(.data[[fc1_col]]) & !is.na(.data[[fc2_col]]))
  
  r_val <- cor(df_pair[[fc1_col]], df_pair[[fc2_col]],
               method = "spearman", use = "complete.obs")
  
  p <- ggplot(df_pair, aes(x = .data[[fc1_col]], y = .data[[fc2_col]])) +
    geom_point(alpha = 0.3, size = 1.5, color = "grey40") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = FALSE, color = "#457b9d", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("rho = %.3f", r_val),
             hjust = -0.2, vjust = 1.5, size = 5, fontface = "italic",
             color = "#e63946") +
    labs(title = paste0("log2FC: ", m1, " vs ", m2),
         x = paste0("log2FC (", m1, ")"),
         y = paste0("log2FC (", m2, ")")) +
    shared_theme +
    coord_fixed(ratio = 1)
  
  fc_scatter_list[[paste0(m1, "_vs_", m2)]] <- p
}

# =============================================================================
# 20. SUMMARY STATISTICS TABLE
# =============================================================================
message("  Generating summary statistics table...")

summary_rows <- lapply(methods_ok, function(m) {
  df <- m$data %>% dplyr::filter(!is.na(pval_int))
  n_tested  <- nrow(df)
  n_sig_p   <- sum(df$pval_int < 0.05, na.rm = TRUE)
  n_sig_fdr <- sum(df$padj_int < 0.05, na.rm = TRUE)
  
  sig_genes <- df %>% dplyr::filter(padj_int < 0.05)
  median_fc     <- if (nrow(sig_genes) > 0) median(abs(sig_genes$logFC_int), na.rm = TRUE) else NA
  median_fc_all <- median(abs(df$logFC_int), na.rm = TRUE)
  
  # Direction breakdown
  n_up   <- sum(sig_genes$logFC_int > 0, na.rm = TRUE)
  n_down <- sum(sig_genes$logFC_int < 0, na.rm = TRUE)
  
  data.frame(
    Method             = m$name,
    Contrast           = m$test_name,
    Genes_Tested       = n_tested,
    Sig_p05            = n_sig_p,
    Sig_FDR05          = n_sig_fdr,
    Pct_Sig_FDR        = sprintf("%.1f%%", 100 * n_sig_fdr / n_tested),
    Up_FDR05           = n_up,
    Down_FDR05         = n_down,
    Median_absFC_Sig   = round(median_fc, 3),
    Median_absFC_All   = round(median_fc_all, 3),
    Time_sec           = round(timing[paste0(tolower(gsub("Wilcoxon", "wilcox", m$name)), ".elapsed")], 1),
    stringsAsFactors   = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)

message("  Summary table:")
print(summary_table)

# Create a table grob for the PDF
summary_grob <- tableGrob(
  summary_table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 10,
    core    = list(bg_params = list(fill = c("grey95", "white"), col = NA)),
    colhead = list(fg_params = list(fontface = "bold"),
                   bg_params = list(fill = "#457b9d", col = NA),
                   fg_params = list(col = "white", fontface = "bold"))
  )
)

# =============================================================================
# 21. SAVE EVERYTHING TO PDF AND UPLOAD TO S3
# =============================================================================
message("  Saving all plots to PDF...")

tmp_pdf <- tempfile(fileext = ".pdf")
pdf(tmp_pdf, width = 14, height = 10)

# ── Page 1: Summary table ──
grid::grid.newpage()
grid::grid.draw(
  gridExtra::arrangeGrob(
    grid::textGrob(
      paste0("DE Method Comparison — ATTEMPT Reference (", opt$cell_type, ")"),
      gp = grid::gpar(fontsize = 18, fontface = "bold")
    ),
    grid::textGrob(subtitle_text,
                   gp = grid::gpar(fontsize = 12, col = "grey40")
    ),
    summary_grob,
    ncol = 1,
    heights = grid::unit(c(0.08, 0.05, 0.87), "npc")
  )
)

# ── Volcano plots: individual pages ──
for (pname in names(volcano_list)) {
  p <- volcano_list[[pname]] + labs(subtitle = subtitle_text)
  print(p)
}

# ── Volcano plots: combined side-by-side (pval | fdr) per method ──
for (m in methods_ok) {
  pval_key <- paste0(m$name, "_pval")
  fdr_key  <- paste0(m$name, "_fdr")
  if (pval_key %in% names(volcano_list) && fdr_key %in% names(volcano_list)) {
    combined <- volcano_list[[pval_key]] + volcano_list[[fdr_key]] +
      plot_annotation(
        title    = paste0(m$name, " — Interaction Volcano Plots (ATTEMPT Reference)"),
        subtitle = subtitle_text,
        theme    = theme(
          plot.title    = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey40")
        )
      )
    print(combined)
  }
}

# ── P-value distribution: faceted ──
print(pval_faceted)

# ── P-value distribution: overlay ──
print(pval_overlay)

# ── P-value distribution: individual histograms ──
for (pname in names(pval_hist_list)) {
  print(pval_hist_list[[pname]])
}

# ── Rank correlation heatmaps ──
cor_combined <- cor_logfc_plot + cor_pval_plot +
  plot_annotation(
    title    = "Method Agreement — Spearman Rank Correlations",
    subtitle = subtitle_text,
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )
print(cor_combined)

# ── UpSet plot: overlap of FDR-significant gene sets ──
# Only include methods that have at least 1 significant gene AND appear as columns
upset_sets <- intersect(
  names(sig_genes_list)[sapply(sig_genes_list, length) > 0],
  colnames(upset_df)
)
if (length(upset_sets) >= 2) {
  # UpSetR uses base graphics, so we print directly
  tryCatch({
    print(
      upset(upset_df,
            sets        = rev(upset_sets),
            nsets       = length(upset_sets),
            keep.order  = TRUE,
            order.by    = "freq",
            decreasing  = TRUE,
            mb.ratio    = c(0.6, 0.4),
            text.scale  = c(1.5, 1.2, 1.2, 1.0, 1.5, 1.2),
            main.bar.color  = "#457b9d",
            sets.bar.color  = method_colors[rev(upset_sets)],
            mainbar.y.label = "Intersection size",
            sets.x.label    = "Significant genes (FDR < 0.05)")
    )
  }, error = function(e) {
    message("  UpSet plot failed (non-fatal): ", conditionMessage(e))
  })
} else {
  message("  Skipping UpSet plot: fewer than 2 methods have FDR-significant genes")
}

# ── log2FC scatter plots between method pairs ──
# Print in pairs of 2 per page if possible
scatter_names <- names(fc_scatter_list)
if (length(scatter_names) > 0) {
  pairs_per_page <- 2
  for (i in seq(1, length(scatter_names), by = pairs_per_page)) {
    idx <- i:min(i + pairs_per_page - 1, length(scatter_names))
    if (length(idx) == 2) {
      combined_fc <- fc_scatter_list[[scatter_names[idx[1]]]] +
        fc_scatter_list[[scatter_names[idx[2]]]] +
        plot_annotation(
          title    = "log2FC Agreement Between Methods",
          subtitle = subtitle_text,
          theme    = theme(
            plot.title    = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 10, color = "grey40")
          )
        )
      print(combined_fc)
    } else {
      print(fc_scatter_list[[scatter_names[idx[1]]]] +
              labs(subtitle = subtitle_text))
    }
  }
}

# ── Pairwise p-value scatter (log10-log10) for all method pairs ──
pval_scatter_list <- list()
for (pair in method_pairs) {
  m1 <- pair[1]; m2 <- pair[2]
  pv1_col <- paste0("pval_", m1)
  pv2_col <- paste0("pval_", m2)
  if (!pv1_col %in% names(merged_all) || !pv2_col %in% names(merged_all)) next
  
  df_pair <- merged_all %>%
    dplyr::select(gene, all_of(c(pv1_col, pv2_col))) %>%
    dplyr::filter(!is.na(.data[[pv1_col]]) & !is.na(.data[[pv2_col]]))
  
  p <- ggplot(df_pair, aes(x = -log10(.data[[pv1_col]] + 1e-300),
                           y = -log10(.data[[pv2_col]] + 1e-300))) +
    geom_point(alpha = 0.3, size = 1.5, color = "grey40") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste0("-log10(p): ", m1, " vs ", m2),
         x = paste0("-log10(p) ", m1),
         y = paste0("-log10(p) ", m2)) +
    shared_theme
  
  pval_scatter_list[[paste0(m1, "_vs_", m2)]] <- p
}

pv_scatter_names <- names(pval_scatter_list)
if (length(pv_scatter_names) > 0) {
  for (i in seq(1, length(pv_scatter_names), by = 2)) {
    idx <- i:min(i + 1, length(pv_scatter_names))
    if (length(idx) == 2) {
      combined_pv <- pval_scatter_list[[pv_scatter_names[idx[1]]]] +
        pval_scatter_list[[pv_scatter_names[idx[2]]]] +
        plot_annotation(
          title    = "P-value Agreement Between Methods (-log10 scale)",
          subtitle = subtitle_text,
          theme    = theme(
            plot.title    = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 10, color = "grey40")
          )
        )
      print(combined_pv)
    } else {
      print(pval_scatter_list[[pv_scatter_names[idx[1]]]] +
              labs(subtitle = subtitle_text))
    }
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

# Also save stats results + summary table to S3 for downstream use
ref_results <- list(
  nebula_stats  = nebula_stats,
  edger_stats   = edger_stats,
  limma_stats   = limma_stats,
  wilcox_stats  = wilcox_stats,
  mast_stats    = mast_stats,
  summary_table = summary_table,
  timing        = timing,
  cor_logfc     = cor_logfc,
  cor_pval      = cor_pval,
  sig_genes_list = sig_genes_list,
  cell_type     = opt$cell_type,
  n_hvg         = opt$n_hvg,
  n_subjects    = n_subjects
)
s3saveRDS(ref_results,
          object = paste0(S3_BASE, "reference/de_results_", opt$cell_type, ".rds"),
          bucket = S3_BUCKET,
          region = "")

# =============================================================================
# 22. SAVE CSV SUMMARIES TO S3
# =============================================================================
message("  Saving CSV files to S3...")

# Helper to write a CSV to S3
s3write_csv <- function(df, s3_path) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(df, file = tmp, row.names = FALSE)
  aws.s3::put_object(file = tmp, object = s3_path, bucket = S3_BUCKET, region = "")
}

S3_CSV_PFX <- paste0(S3_BASE, "reference/csv/")

# 1. Summary statistics table
s3write_csv(summary_table,
            paste0(S3_CSV_PFX, "summary_stats_", opt$cell_type, ".csv"))
message("  Saved: summary_stats.csv")

# 2. Per-method full results (gene, logFC, pval, padj)
for (m in methods_ok) {
  df <- m$data %>%
    dplyr::filter(!is.na(pval_int)) %>%
    dplyr::arrange(pval_int)
  s3write_csv(df,
              paste0(S3_CSV_PFX, tolower(m$name), "_results_", opt$cell_type, ".csv"))
  message(sprintf("  Saved: %s_results.csv (%d genes)", tolower(m$name), nrow(df)))
}

# 3. Significant genes per method (FDR < 0.05)
for (m in methods_ok) {
  df_sig <- m$data %>%
    dplyr::filter(!is.na(padj_int) & padj_int < 0.05) %>%
    dplyr::arrange(padj_int)
  s3write_csv(df_sig,
              paste0(S3_CSV_PFX, tolower(m$name), "_sig_genes_", opt$cell_type, ".csv"))
  message(sprintf("  Saved: %s_sig_genes.csv (%d genes)", tolower(m$name), nrow(df_sig)))
}

# 4. Merged results across all methods (wide format)
merged_wide <- NULL
for (m in methods_ok) {
  tmp <- m$data %>%
    dplyr::filter(!is.na(pval_int)) %>%
    dplyr::select(gene, logFC_int, pval_int, padj_int) %>%
    dplyr::rename(
      !!paste0("logFC_", m$name)  := logFC_int,
      !!paste0("pval_", m$name)   := pval_int,
      !!paste0("padj_", m$name)   := padj_int
    )
  if (is.null(merged_wide)) {
    merged_wide <- tmp
  } else {
    merged_wide <- dplyr::full_join(merged_wide, tmp, by = "gene")
  }
}
s3write_csv(merged_wide,
            paste0(S3_CSV_PFX, "all_methods_merged_", opt$cell_type, ".csv"))
message(sprintf("  Saved: all_methods_merged.csv (%d genes x %d cols)",
                nrow(merged_wide), ncol(merged_wide)))

# 5. Correlation matrices
cor_logfc_df <- as.data.frame(as.table(cor_logfc))
names(cor_logfc_df) <- c("Method1", "Method2", "Spearman_rho_logFC")
s3write_csv(cor_logfc_df,
            paste0(S3_CSV_PFX, "correlation_logfc_", opt$cell_type, ".csv"))

cor_pval_df <- as.data.frame(as.table(cor_pval))
names(cor_pval_df) <- c("Method1", "Method2", "Spearman_rho_pval")
s3write_csv(cor_pval_df,
            paste0(S3_CSV_PFX, "correlation_pval_", opt$cell_type, ".csv"))
message("  Saved: correlation_logfc.csv, correlation_pval.csv")

# 6. Gene overlap matrix (which methods call each gene significant)
overlap_df <- data.frame(gene = all_genes_tested, upset_df, stringsAsFactors = FALSE) %>%
  dplyr::mutate(n_methods_sig = rowSums(dplyr::across(all_of(method_names)))) %>%
  dplyr::arrange(desc(n_methods_sig))
s3write_csv(overlap_df,
            paste0(S3_CSV_PFX, "gene_sig_overlap_", opt$cell_type, ".csv"))
message(sprintf("  Saved: gene_sig_overlap.csv (%d genes)", nrow(overlap_df)))

message(sprintf("== [05-ref] Done. All outputs saved to s3://%s/%s ==",
                S3_BUCKET, S3_BASE))
message(sprintf("  PDF:  s3://%s/%s", S3_BUCKET, s3_output_path))
message(sprintf("  RDS:  s3://%s/%sde_results_%s.rds", S3_BUCKET,
                paste0(S3_BASE, "reference/"), opt$cell_type))
message(sprintf("  CSVs: s3://%s/%s", S3_BUCKET, S3_CSV_PFX))

