# Core
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(splatter)
  library(Matrix)
  library(nebula)
  library(Seurat)
  library(scater)
  library(progressr)
  library(aws.s3)
  library(jsonlite)
  library(edgeR)
  library(MAST)
})

handlers(global = TRUE)
handlers("txtprogressbar")   # console-friendly

keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

s3write_using_region <- function(FUN, ..., object, bucket, region = NULL, opts = NULL, filename = NULL) {
  if (missing(bucket)) {
    bucket <- get_bucketname(object)
  }
  object <- get_objectkey(object)
  
  tmp <- if (is.character(filename)) {
    file.path(tempdir(TRUE), filename)
  } else {
    # if object has an extension, keep it; otherwise make a generic tmp
    ext <- tools::file_ext(object)
    if (nzchar(ext)) tempfile(fileext = paste0(".", ext)) else tempfile()
  }
  
  on.exit(unlink(tmp), add = TRUE)
  
  # Add region to opts if provided
  if (!is.null(region)) {
    if (is.null(opts)) {
      opts <- list(region = region)
    } else {
      opts$region <- region
    }
  }
  
  FUN(tmp, ...)
  
  if (is.null(opts)) {
    r <- put_object(file = tmp, bucket = bucket, object = object)
  } else {
    r <- do.call("put_object", c(list(file = tmp, bucket = bucket, object = object), opts))
  }
  
  return(invisible(r))
}

write_csv_s3 <- function(file, x) {
  write.csv(x, file = file, row.names = FALSE)
}

# --- fixed simulation settings ---
n_genes <- 2000
n_simulations <- 500

# --- variable factor simulation settings ---
study_design <- c("paired", "unpaired") #nGroups
n_subjects <- c(5, 10 , 20) # nBatches
cells_per_subject <- c(1000, 3000, 6000) #batchCells
btw_subject_var <- c(0.001, 0.01, 0.1, 0.2) # batch.facScale
gene_wise_disp <- c(0.1, 0.2, 0.4, 0.6) # bcv.common
effect_size <- c(0.01, 0.05, 0.1, 0.2, 0.3) # de.prob

# ---- parameter grid (all combinations) ----
grid <- expand.grid(
  study_design      = study_design,
  n_subjects        = n_subjects,
  cells_per_subject = cells_per_subject,
  btw_subject_var   = btw_subject_var,
  gene_wise_disp    = gene_wise_disp,
  effect_size       = effect_size,
  KEEP.OUT.ATTRS    = FALSE,
  stringsAsFactors  = FALSE
)

grid$combo_id <- sprintf("combo_%04d", seq_len(nrow(grid)))

s3write_using_region(FUN = write_csv_s3,
                     grid,
                     object = "scRNA simulation analysis/simulation_grid.csv", 
                     bucket = "simulation", 
                     region = "")

# Splatter params
# ---- base params ----
base_params <- newSplatParams()

# ---- helper to simulate one dataset for a given grid row ----
fmt_time <- function(t0) {
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  hrs <- floor(secs / 3600)
  mins <- floor((secs %% 3600) / 60)
  ss <- round(secs %% 60)
  sprintf("%02dh:%02dm:%02ds", hrs, mins, ss)
}

fast_auc <- function(score, truth01) {
  ok <- is.finite(score) & !is.na(score) & !is.na(truth01)
  score <- score[ok]; truth01 <- truth01[ok]
  if (length(unique(truth01)) < 2) return(NA_real_)
  r <- rank(score, ties.method = "average")
  n1 <- sum(truth01 == 1)
  n0 <- sum(truth01 == 0)
  as.numeric((sum(r[truth01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0))
}

null_pval_diag <- function(pnull) {
  pnull <- pnull[is.finite(pnull) & !is.na(pnull)]
  if (length(pnull) < 10) return(list(ks_p = NA_real_, hist = rep(NA_real_, 10)))
  ks_p <- suppressWarnings(stats::ks.test(pnull, "punif")$p.value)
  hist <- hist(pnull, breaks = seq(0, 1, by = 0.1), plot = FALSE)$counts
  list(ks_p = as.numeric(ks_p), hist = as.numeric(hist))
}

mem_mb <- function() sum(gc()[, 2])  # MB (R heap used). Not OS peak RSS.

safe_s3_exists <- function(bucket, object) {
  # aws.s3::object_exists() returns TRUE/FALSE; wrap for safety
  out <- tryCatch(aws.s3::object_exists(object = object, bucket = bucket),
                  error = function(e) FALSE)
  isTRUE(out)
}

safe_s3saveRDS <- function(x, bucket, object, region = "") {
  tryCatch({
    s3saveRDS(x, bucket = bucket, object = object, region = region)
    TRUE
  }, error = function(e) {
    message("S3 save failed for: ", object, " | ", conditionMessage(e))
    FALSE
  })
}

# helper for compiling metrics
calc_metrics <- function(res, truth_is_de, true_lfc, alpha = 0.05) {
  genes <- unique(res$gene)
  truth01 <- as.integer(truth_is_de[match(genes, names(truth_is_de))])
  truth01[is.na(truth01)] <- 0L
  
  res2 <- res[match(genes, res$gene), ]
  pval <- res2$pval
  sig <- is.finite(pval) & !is.na(pval) & (pval < alpha)
  
  nde <- sum(truth01 == 1)
  nnull <- sum(truth01 == 0)
  
  power <- if (nde > 0) mean(sig[truth01 == 1], na.rm = TRUE) else NA_real_
  type1 <- if (nnull > 0) mean(sig[truth01 == 0], na.rm = TRUE) else NA_real_
  
  # bias/mse for DE genes only (where truth defined)
  est <- res2$logFC
  tru <- true_lfc[match(genes, names(true_lfc))]
  ok_eff <- (truth01 == 1) & is.finite(est) & !is.na(est) & is.finite(tru) & !is.na(tru)
  bias <- if (sum(ok_eff) > 0) mean(est[ok_eff] - tru[ok_eff]) else NA_real_
  mse  <- if (sum(ok_eff) > 0) mean((est[ok_eff] - tru[ok_eff])^2) else NA_real_
  
  score <- res2$stat
  if (all(!is.finite(score) | is.na(score))) score <- -log10(pval)
  auc <- fast_auc(score, truth01)
  
  null_diag <- null_pval_diag(pval[truth01 == 0])
  
  list(
    power = power,
    type1 = type1,
    bias = bias,
    mse = mse,
    auc = auc,
    null_ks_p = null_diag$ks_p,
    null_hist = null_diag$hist
  )
}

# helper for each simulation
simulate_one <- function(row, base_params, n_genes) {
  design <- row[["study_design"]]
  nsubj  <- row[["n_subjects"]]
  cps    <- row[["cells_per_subject"]]
  
  batch_cells <- rep(cps, nsubj)
  
  # IMPORTANT: we set de.prob = 0 in Splatter and inject DE ourselves
  if (design == "unpaired") {
    params <- setParams(
      base_params,
      nGenes         = n_genes,
      batchCells     = batch_cells,
      batch.facScale = row[["btw_subject_var"]],
      bcv.common     = row[["gene_wise_disp"]],
      group.prob     = c(0.5, 0.5),
      de.prob        = 0
    )
    sce <- splatSimulate(params, method = "groups", verbose = FALSE)
    
  } else {
    pop <- setParams(
      newSplatPopParams(),
      nGenes         = n_genes,
      batchCells     = batch_cells,
      batch.facScale = row[["btw_subject_var"]],
      bcv.common     = row[["gene_wise_disp"]]
    )
    sce <- splatSimulate(pop, method = "single", verbose = FALSE)
  }
  
  counts <- SummarizedExperiment::assay(sce, "counts")
  ng <- nrow(counts)
  nc <- ncol(counts)
  
  # make gene ids stable
  if (is.null(rownames(counts))) rownames(counts) <- paste0("gene_", seq_len(ng))
  genes <- rownames(counts)
  
  # ----- build cell-level metadata -----
  # subject (batch) from Splatter is in colData(sce)$Batch usually
  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!("Batch" %in% names(meta))) {
    stop("Expected Splatter to provide Batch in colData; not found.")
  }
  meta$subject <- factor(meta$Batch)
  
  # subject -> group assignment (balanced)
  subj_levels <- levels(meta$subject)
  nA <- floor(length(subj_levels) / 2)
  group_by_subj <- setNames(c(rep("A", nA), rep("B", length(subj_levels) - nA)), subj_levels)
  meta$group <- factor(group_by_subj[as.character(meta$subject)], levels = c("A","B"))
  
  if (design == "paired") {
    # within each subject: split cells into t1/t2
    meta$time <- NA_character_
    for (sj in subj_levels) {
      idx <- which(meta$subject == sj)
      # deterministic split within subject
      n1 <- floor(length(idx) / 2)
      meta$time[idx[seq_len(n1)]] <- "t1"
      meta$time[idx[-seq_len(n1)]] <- "t2"
    }
    meta$time <- factor(meta$time, levels = c("t1","t2"))
  } else {
    meta$time <- factor(NA_character_, levels = c("t1","t2"))
  }
  
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(meta)
  
  # ----- inject DE and store truth in rowData -----
  de_prob <- row[["effect_size"]]  # you originally used this as de.prob; now treated as LFC magnitude?
  # If you want *both* de.prob and lfc, split them; for now:
  #   - fraction DE genes is controlled by row$effect_size? (not good)
  # Better: interpret row$effect_size as LFC, and use a fixed DE fraction.
  # Here: assume 10% DE genes unless you add de_frac to grid.
  de_frac <- 0.10
  lfc <- row[["effect_size"]]  # natural-log fold-change
  
  is_de <- stats::runif(ng) < de_frac
  true_lfc <- rep(0, ng)
  true_lfc[is_de] <- lfc
  names(is_de) <- genes
  names(true_lfc) <- genes
  
  # Apply effects:
  if (design == "unpaired") {
    # group effect: B vs A
    b_cells <- which(meta$group == "B")
    if (length(b_cells) > 0) {
      counts[is_de, b_cells] <- round(counts[is_de, b_cells] * exp(lfc))
    }
  } else {
    # interaction: only B at t2 gets effect
    b_t2 <- which(meta$group == "B" & meta$time == "t2")
    if (length(b_t2) > 0) {
      counts[is_de, b_t2] <- round(counts[is_de, b_t2] * exp(lfc))
    }
  }
  
  SummarizedExperiment::assay(sce, "counts") <- counts
  SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(
    is_de = is_de,
    true_lfc = true_lfc
  )
  
  sce
}

# nebula function
run_nebula <- function(sce) {
  counts <- SummarizedExperiment::assay(sce, "counts")
  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  
  # nebula uses id + predictors
  id <- meta$subject
  
  if (all(is.na(meta$time))) {
    # unpaired: group only
    pred <- data.frame(group = meta$group)
    fit <- nebula::nebula(counts = counts, id = id, pred = pred)
    sm <- as.data.frame(fit$summary)
    sm$gene <- rownames(sm)
    out <- data.frame(
      gene  = sm$gene,
      pval  = sm$p,          # adjust if your nebula summary differs
      logFC = sm$beta,       # group effect
      stat  = sm$wald,
      stringsAsFactors = FALSE
    )
    return(out)
  }
  
  # paired: group*time interaction
  pred <- data.frame(group = meta$group, time = meta$time)
  pred$gB <- as.integer(pred$group == "B")
  pred$t2 <- as.integer(pred$time == "t2")
  pred$gBxt2 <- pred$gB * pred$t2
  
  fit <- nebula::nebula(counts = counts, id = id, pred = pred[, c("gB","t2","gBxt2")])
  sm <- as.data.frame(fit$summary)
  sm$gene <- rownames(sm)
  
  # We want the interaction coefficient. nebula summary column naming can differ.
  # If beta is a matrix, you need to extract the gBxt2 column. Handle both.
  if (is.matrix(sm$beta)) {
    logfc <- sm$beta[, "gBxt2"]
    pval  <- sm$p[, "gBxt2"]
    stat  <- sm$wald[, "gBxt2"]
  } else {
    # fallback: assume returned for last term (not ideal but avoids crashing)
    logfc <- sm$beta
    pval  <- sm$p
    stat  <- sm$wald
  }
  
  data.frame(
    gene  = sm$gene,
    pval  = as.numeric(pval),
    logFC = as.numeric(logfc),
    stat  = as.numeric(stat),
    stringsAsFactors = FALSE
  )
}

# pseudobulk edgeR
run_pseudobulk_edger <- function(sce) {
  counts <- SummarizedExperiment::assay(sce, "counts")
  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  
  if (all(is.na(meta$time))) {
    # aggregate per subject
    key <- meta$subject
    pb <- rowsum(t(counts), group = key, reorder = FALSE)  # subj x genes
    pb <- t(pb)                                            # genes x subj
    
    info <- data.frame(subject = colnames(pb), group = tapply(as.character(meta$group), key, `[`, 1),
                       stringsAsFactors = FALSE)
    info$group <- factor(info$group, levels = c("A","B"))
    
    y <- edgeR::DGEList(pb)
    y <- edgeR::calcNormFactors(y)
    design <- model.matrix(~ group, data = info)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)  # groupB
    
    tt <- edgeR::topTags(qlf, n = Inf)$table
    return(data.frame(gene = rownames(tt), pval = tt$PValue, logFC = tt$logFC, stat = tt$F,
                      stringsAsFactors = FALSE))
  }
  
  # paired: aggregate per subject x time
  key <- paste(meta$subject, meta$time, sep = "||")
  pb <- rowsum(t(counts), group = key, reorder = FALSE)  # (subj,time) x genes
  pb <- t(pb)                                            # genes x samples
  
  parts <- do.call(rbind, strsplit(colnames(pb), "\\|\\|"))
  info <- data.frame(subject = parts[,1], time = parts[,2], stringsAsFactors = FALSE)
  info$time <- factor(info$time, levels = c("t1","t2"))
  
  # group per subject (stable)
  subj_group <- tapply(as.character(meta$group), meta$subject, `[`, 1)
  info$group <- factor(subj_group[info$subject], levels = c("A","B"))
  
  y <- edgeR::DGEList(pb)
  y <- edgeR::calcNormFactors(y)
  
  design <- model.matrix(~ group * time, data = info)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  
  # interaction term is groupB:timet2 (usually last column)
  cn <- colnames(design)
  coef_int <- which(grepl("groupB:timet2", cn, fixed = TRUE))
  if (length(coef_int) != 1) coef_int <- ncol(design)
  
  qlf <- edgeR::glmQLFTest(fit, coef = coef_int)
  tt <- edgeR::topTags(qlf, n = Inf)$table
  
  data.frame(gene = rownames(tt), pval = tt$PValue, logFC = tt$logFC, stat = tt$F,
             stringsAsFactors = FALSE)
}

# MAST
run_mast_zlm <- function(sce) {
  counts <- SummarizedExperiment::assay(sce, "counts")
  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  
  # MAST expects log-expression (typically log2 CPM-ish); use log1p counts here (simple, stable)
  log_expr <- log1p(counts)
  
  # Create SingleCellAssay
  fdat <- data.frame(primerid = rownames(log_expr), stringsAsFactors = FALSE)
  cdat <- meta
  cdat$cngeneson <- colSums(log_expr > 0)
  
  sca <- MAST::FromMatrix(exprsArray = as.matrix(log_expr),
                          cData = cdat, fData = fdat)
  
  if (all(is.na(meta$time))) {
    # unpaired: ~ group + cngeneson
    sca$group <- relevel(factor(sca$group), ref = "A")
    z <- MAST::zlm(~ group + cngeneson, sca, method = "glm", ebayes = FALSE)
    s <- MAST::summary(z, doLRT = "groupB")$datatable
    
    # extract hurdle p-values + coefficient
    p <- s[s$component == "H" & s$contrast == "groupB", c("primerid","Pr(>Chisq)")]
    b <- s[s$component == "logFC" & s$contrast == "groupB", c("primerid","coef")]
    out <- merge(p, b, by = "primerid", all = TRUE)
    colnames(out) <- c("gene","pval","logFC")
    out$stat <- -log10(out$pval)
    return(out)
  }
  
  # paired interaction: ~ group * time + cngeneson
  sca$group <- relevel(factor(sca$group), ref = "A")
  sca$time  <- relevel(factor(sca$time),  ref = "t1")
  
  z <- MAST::zlm(~ group * time + cngeneson, sca, method = "glm", ebayes = FALSE)
  s <- MAST::summary(z, doLRT = "groupB:timet2")$datatable
  
  p <- s[s$component == "H" & s$contrast == "groupB:timet2", c("primerid","Pr(>Chisq)")]
  b <- s[s$component == "logFC" & s$contrast == "groupB:timet2", c("primerid","coef")]
  out <- merge(p, b, by = "primerid", all = TRUE)
  colnames(out) <- c("gene","pval","logFC")
  out$stat <- -log10(out$pval)
  out
}

# Seurat FindMarkers()
run_seurat_findmarkers <- function(sce) {
  counts <- SummarizedExperiment::assay(sce, "counts")
  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  
  if (all(is.na(meta$time))) {
    seu <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
    Seurat::Idents(seu) <- seu$group
    res <- Seurat::FindMarkers(seu, ident.1 = "B", ident.2 = "A", test.use = "wilcox")
    res$gene <- rownames(res)
    out <- data.frame(gene = res$gene, pval = res$p_val, logFC = res$avg_log2FC,
                      stat = if ("statistic" %in% names(res)) res$statistic else -log10(res$p_val),
                      stringsAsFactors = FALSE)
    return(out)
  }
  
  # Paired interaction using per-subject delta pseudo-cells:
  #   For each subject: pseudobulk at t1 and t2 (sum counts), take log1p and difference.
  #   Then FindMarkers compares delta between groups -> tests interaction.
  key <- paste(meta$subject, meta$time, sep = "||")
  pb <- rowsum(t(counts), group = key, reorder = FALSE)  # (subj,time) x genes
  pb <- t(pb)                                            # genes x samples
  
  parts <- do.call(rbind, strsplit(colnames(pb), "\\|\\|"))
  info <- data.frame(subject = parts[,1], time = parts[,2], stringsAsFactors = FALSE)
  
  subj_group <- tapply(as.character(meta$group), meta$subject, `[`, 1)
  
  # build delta matrix (genes x subjects)
  subjects <- sort(unique(info$subject))
  delta <- matrix(NA_real_, nrow = nrow(pb), ncol = length(subjects),
                  dimnames = list(rownames(pb), subjects))
  
  for (sj in subjects) {
    c1 <- paste0(sj, "||t1")
    c2 <- paste0(sj, "||t2")
    if (c1 %in% colnames(pb) && c2 %in% colnames(pb)) {
      delta[, sj] <- log1p(pb[, c2]) - log1p(pb[, c1])
    }
  }
  
  keep <- colSums(is.na(delta)) == 0
  delta <- delta[, keep, drop = FALSE]
  subjects <- colnames(delta)
  
  meta2 <- data.frame(
    subject = subjects,
    group = factor(subj_group[subjects], levels = c("A","B")),
    stringsAsFactors = FALSE
  )
  # "counts" must be nonnegative integers; we store delta in the data slot instead:
  seu <- Seurat::CreateSeuratObject(counts = matrix(0L, nrow(delta), ncol(delta),
                                                    dimnames = dimnames(delta)),
                                    meta.data = meta2)
  seu[["delta"]] <- Seurat::CreateAssayObject(data = delta)  # put delta into assay data
  Seurat::DefaultAssay(seu) <- "delta"
  Seurat::Idents(seu) <- seu$group
  
  res <- Seurat::FindMarkers(seu, ident.1 = "B", ident.2 = "A", test.use = "wilcox",
                             slot = "data")  # use delta assay "data"
  res$gene <- rownames(res)
  out <- data.frame(gene = res$gene, pval = res$p_val,
                    logFC = if ("avg_log2FC" %in% names(res)) res$avg_log2FC else NA_real_,
                    stat = if ("statistic" %in% names(res)) res$statistic else -log10(res$p_val),
                    stringsAsFactors = FALSE)
  out
}

# ---- user settings ----
alpha <- 0.05
save_gene_level <- FALSE   # set TRUE only if you really want per-gene tables saved
seed_checkpoint_every <- 25 # save seed_log every N completed runs

bucket <- "simulation"
prefix <- "scRNA simulation analysis/sim_results"
region <- ""

# RNG for reproducibility + uniqueness
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

total_runs <- nrow(grid) * n_simulations
t0 <- Sys.time()

# We'll keep seed_log small and checkpoint it.
seed_log <- data.frame(
  combo_id = character(total_runs),
  g        = integer(total_runs),
  sim      = integer(total_runs),
  seed1    = integer(total_runs),
  status   = character(total_runs),
  stringsAsFactors = FALSE
)

with_progress({
  p <- progressor(steps = total_runs)
  idx <- 0L
  completed <- 0L
  
  for (g in seq_len(nrow(grid))) {
    row <- as.list(grid[g, , drop = FALSE])
    combo_start <- Sys.time()
    
    for (s in seq_len(n_simulations)) {
      idx <- idx + 1L
      run_id <- sprintf("%s_rep%04d", row$combo_id, s)
      summary_obj_path <- file.path(prefix, "summary", row$combo_id, paste0(run_id, "_summary.rds"))
      
      # If summary already exists, skip (restart-safe)
      if (safe_s3_exists(bucket, summary_obj_path)) {
        seed_log$combo_id[idx] <- row$combo_id
        seed_log$g[idx]        <- g
        seed_log$sim[idx]      <- s
        seed_log$seed1[idx]    <- NA_integer_
        seed_log$status[idx]   <- "skipped_exists"
        p(sprintf("SKIP %s | Elapsed %s", run_id, fmt_time(t0)))
        next
      }
      
      # record seed BEFORE advancing
      seed_log$combo_id[idx] <- row$combo_id
      seed_log$g[idx]        <- g
      seed_log$sim[idx]      <- s
      seed_log$seed1[idx]    <- .Random.seed[1]
      seed_log$status[idx]   <- "started"
      
      # new independent RNG stream
      .Random.seed <- parallel::nextRNGStream(.Random.seed)
      
      # simulate
      sim_start <- Sys.time()
      mem_before_all <- mem_mb()
      
      sce <- NULL
      method_results <- list()
      method_metrics <- list()
      
      # wrapper to make methods error-proof
      run_one_method <- function(name, fun) {
        mem0 <- mem_mb()
        t0m <- Sys.time()
        out <- tryCatch(fun(), error = function(e) e)
        t1m <- Sys.time()
        mem1 <- mem_mb()
        
        if (inherits(out, "error")) {
          list(
            res = NULL,
            metrics = list(
              error = conditionMessage(out),
              runtime_s = as.numeric(difftime(t1m, t0m, units = "secs")),
              mem_delta_mb = max(0, mem1 - mem0),
              mem_total_mb = mem1,
              power = NA_real_, type1 = NA_real_, bias = NA_real_, mse = NA_real_, auc = NA_real_,
              null_ks_p = NA_real_, null_hist = rep(NA_real_, 10)
            )
          )
        } else {
          list(
            res = out,
            metrics = NULL,  # computed after we know truth
            runtime_s = as.numeric(difftime(t1m, t0m, units = "secs")),
            mem_delta_mb = max(0, mem1 - mem0),
            mem_total_mb = mem1
          )
        }
      }
      
      # ---- simulate_one() with tryCatch ----
      sce <- tryCatch(simulate_one(row, base_params, n_genes),
                      error = function(e) e)
      
      if (inherits(sce, "error")) {
        # save failure summary and continue
        fail_summary <- list(
          combo = row, run_id = run_id, seed1 = seed_log$seed1[idx], alpha = alpha,
          error = paste0("simulate_one failed: ", conditionMessage(sce)),
          elapsed_total_s = as.numeric(difftime(Sys.time(), t0, units = "secs"))
        )
        safe_s3saveRDS(fail_summary, bucket, summary_obj_path, region)
        
        seed_log$status[idx] <- "failed_sim"
        p(sprintf("FAIL(sim) %s | Elapsed %s", run_id, fmt_time(t0)))
        next
      }
      
      truth_is_de <- SummarizedExperiment::rowData(sce)$is_de
      true_lfc <- SummarizedExperiment::rowData(sce)$true_lfc
      names(truth_is_de) <- rownames(sce)
      names(true_lfc) <- rownames(sce)
      
      # run methods
      tmp <- run_one_method("nebula", function() run_nebula(sce))
      method_results$nebula <- tmp$res
      method_metrics$nebula <- tmp$metrics %||% list()
      method_metrics$nebula$runtime_s <- tmp$runtime_s
      method_metrics$nebula$mem_delta_mb <- tmp$mem_delta_mb
      method_metrics$nebula$mem_total_mb <- tmp$mem_total_mb
      
      tmp <- run_one_method("pseudobulk_edger", function() run_pseudobulk_edger(sce))
      method_results$pseudobulk_edger <- tmp$res
      method_metrics$pseudobulk_edger <- tmp$metrics %||% list()
      method_metrics$pseudobulk_edger$runtime_s <- tmp$runtime_s
      method_metrics$pseudobulk_edger$mem_delta_mb <- tmp$mem_delta_mb
      method_metrics$pseudobulk_edger$mem_total_mb <- tmp$mem_total_mb
      
      tmp <- run_one_method("mast_zlm", function() run_mast_zlm(sce))
      method_results$mast_zlm <- tmp$res
      method_metrics$mast_zlm <- tmp$metrics %||% list()
      method_metrics$mast_zlm$runtime_s <- tmp$runtime_s
      method_metrics$mast_zlm$mem_delta_mb <- tmp$mem_delta_mb
      method_metrics$mast_zlm$mem_total_mb <- tmp$mem_total_mb
      
      tmp <- run_one_method("seurat_wilcox", function() run_seurat_findmarkers(sce))
      method_results$seurat_wilcox <- tmp$res
      method_metrics$seurat_wilcox <- tmp$metrics %||% list()
      method_metrics$seurat_wilcox$runtime_s <- tmp$runtime_s
      method_metrics$seurat_wilcox$mem_delta_mb <- tmp$mem_delta_mb
      method_metrics$seurat_wilcox$mem_total_mb <- tmp$mem_total_mb
      
      # compute metrics where method returned results (no error)
      for (m in names(method_results)) {
        if (is.null(method_results[[m]])) next
        met <- calc_metrics(method_results[[m]], truth_is_de, true_lfc, alpha = alpha)
        # keep runtime/mem already recorded
        method_metrics[[m]] <- c(method_metrics[[m]], met)
      }
      
      mem_after_all <- mem_mb()
      sim_end <- Sys.time()
      
      summary_obj <- list(
        combo = row,
        run_id = run_id,
        seed1 = seed_log$seed1[idx],
        alpha = alpha,
        elapsed_total_s = as.numeric(difftime(Sys.time(), t0, units = "secs")),
        runtime_total_s = as.numeric(difftime(sim_end, sim_start, units = "secs")),
        mem_start_mb = mem_before_all,
        mem_end_mb   = mem_after_all,
        metrics = method_metrics
      )
      
      # save summary first (checkpoint)
      ok <- safe_s3saveRDS(summary_obj, bucket, summary_obj_path, region)
      
      # optionally save gene-level results
      if (isTRUE(save_gene_level)) {
        for (m in names(method_results)) {
          if (is.null(method_results[[m]])) next
          obj_path <- file.path(prefix, "results", row$combo_id, paste0(run_id, "_", m, "_res.rds"))
          safe_s3saveRDS(method_results[[m]], bucket, obj_path, region)
        }
      }
      
      seed_log$status[idx] <- if (ok) "done" else "done_but_save_failed"
      completed <- completed + 1L
      
      # periodic seed_log checkpoint
      if (completed %% seed_checkpoint_every == 0L) {
        safe_s3saveRDS(seed_log[seq_len(idx), ], bucket,
                       file.path(prefix, "seed_log_checkpoint.rds"), region)
      }
      
      # cleanup
      rm(sce, method_results, method_metrics, summary_obj, tmp)
      gc(FALSE)
      
      p(sprintf("DONE %s | Elapsed %s", run_id, fmt_time(t0)))
    }
    
    cat(sprintf(
      "\n[%s] Finished %s (%d/%d) in %s mins | Total elapsed %s\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      row$combo_id, g, nrow(grid),
      format(as.numeric(difftime(Sys.time(), combo_start, units = "mins")), digits = 3),
      fmt_time(t0)
    ))
  }
})

cat(sprintf("\n==== ALL DONE in %s ====\n", fmt_time(t0)))

# final seed log
safe_s3saveRDS(seed_log, bucket, file.path(prefix, "seed_log_final.rds"), region)

# 
# 
# params <- setParams(params,
#                     nGroups = ifelse(study_design == "paired", 4, 2),
#                     nGenes = n_genes,
#                     batchCells = rep(cells_per_subject[i], n_subjects[i]),  # each batch = one subject
#                     batch.facScale = btw_subject_var[i],
#                     bcv.common = gene_wise_disp[i],
#                     group.prob = c(0.5, 0.5),        # within each subject, half pre half post (approx)
#                     de.prob = effect_size[i]
# )
# 
# # Simulate counts
# sce <- splatSimulate(params,
#                      method = "groups",
#                      verbose = T)
# 
# sce <- logNormCounts(sce)
# sce <- runPCA(sce)
# plotPCA(sce, colour_by = "Group")
# 
# # "Batch" in splatter colData is a factor; we use it as subject id
# meta <- as.data.frame(colData(sce))
# meta$subject <- as.integer(factor(meta$Batch))
# meta$time <- factor(meta$Group, levels = levels(meta$Group), labels = c("pre", "post"))
# 
# table(meta$time)
# table(meta$subject)[1:5]
# 
# de_prop <- 0.10
# lfc <- log(1.25)
# 
# n_de <- ceiling(de_prop * n_genes)
# de_genes <- sample(seq_len(n_genes), n_de)
# 
# post_cells <- which(meta$time == "post")
# 
# counts <- assay(sce, "counts")
# # Multiply counts for DE genes in post cells
# counts[de_genes, post_cells] <-
#   round(counts[de_genes, post_cells] * exp(lfc))
# 
# # run nebula
# pred <- model.matrix(~ time, data = meta)
# 
# lib <- Matrix::colSums(counts)
# lib[lib <= 0] <- 1
# offset <- log(lib)
# 
# fit <- nebula(
#   count  = counts,
#   id     = meta$subject,
#   pred   = pred,
#   offset = offset,
#   method = "LN",
#   ncore  = 2,
#   verbose = FALSE
# )
# 
# sm <- as.data.frame(fit$summary)
# 
# coef_name <- colnames(pred)[grepl("^time", colnames(pred))]
# 
# logfc_col <- paste0("logFC_", coef_name)
# p_col     <- paste0("p_", coef_name)
# 
# res <- data.frame(
#   gene = rownames(sm),
#   est_lfc = sm[[logfc_col]],
#   pval = sm[[p_col]]
# )
# 
# summary(res$pval)
# mean(is.na(res$pval))
# sum(res$pval < 0.05, na.rm = TRUE)
# 
# # Add groups A/B and inject a group × time effect
# 
# # Assign half subjects to A, half to B
# subjects <- sort(unique(meta$subject))
# n_subjects <- length(subjects)
# 
# group_by_subject <- rep(c("A", "B"), length.out = n_subjects)
# group_by_subject <- sample(group_by_subject)  # randomize
# names(group_by_subject) <- subjects
# 
# meta$group <- factor(group_by_subject[as.character(meta$subject)], levels = c("A","B"))
# 
# table(meta$group)
# table(meta$group, meta$time)[, ]
# 
# # Re-simulate a clean dataset with same params
# sce2 <- splatSimulate(params, verbose = FALSE)
# counts2 <- assay(sce2, "counts")
# meta2 <- as.data.frame(colData(sce2))
# 
# meta2$subject <- as.integer(factor(meta2$Batch))
# 
# # Same pre/post assignment scheme
# meta2$time <- unlist(
#   lapply(split(seq_len(nrow(meta2)), meta2$subject), function(idx) {
#     n <- length(idx)
#     sample(rep(c("pre", "post"), length.out = n))
#   })
# )
# meta2$time <- factor(meta2$time, levels = c("pre","post"))
# 
# # Same group assignment scheme (reuse mapping so only counts changed)
# meta2$group <- factor(group_by_subject[as.character(meta2$subject)], levels = c("A","B"))
# 
# table(meta2$group, meta2$time)
# 
# de_prop <- 0.10
# lfc_int <- log(1.25)
# 
# n_genes <- nrow(counts2)
# n_de <- ceiling(de_prop * n_genes)
# de_genes <- sample(seq_len(n_genes), n_de)
# 
# B_post_cells <- which(meta2$group == "B" & meta2$time == "post")
# 
# counts2[de_genes, B_post_cells] <-
#   round(counts2[de_genes, B_post_cells] * exp(lfc_int))
# 
# # sanity check: the effect should appear in B post vs B pre, but not in A
# B_pre_cells <- which(meta2$group == "B" & meta2$time == "pre")
# A_post_cells <- which(meta2$group == "A" & meta2$time == "post")
# A_pre_cells  <- which(meta2$group == "A" & meta2$time == "pre")
# 
# ratio_B <- rowMeans(counts2[de_genes, B_post_cells]) / rowMeans(counts2[de_genes, B_pre_cells])
# ratio_A <- rowMeans(counts2[de_genes, A_post_cells]) / rowMeans(counts2[de_genes, A_pre_cells])
# 
# summary(ratio_B)
# summary(ratio_A)
# 
# pred2 <- model.matrix(~ group * time, data = meta2)
# 
# lib2 <- Matrix::colSums(counts2)
# lib2[lib2 <= 0] <- 1
# offset2 <- log(lib2)
# 
# fit2 <- nebula(
#   count  = counts2,
#   id     = meta2$subject,
#   pred   = pred2,
#   offset = offset2,
#   method = "LN",
#   ncore  = 2,
#   verbose = FALSE
# )
# 
# sm2 <- as.data.frame(fit2$summary)
# 
# # Interaction coefficient name should be "groupB:timepost"
# coef_int <- "groupB:timepost"
# 
# res2 <- data.frame(
#   gene = rownames(sm2),
#   est_lfc = sm2[[paste0("logFC_", coef_int)]],
#   pval = sm2[[paste0("p_", coef_int)]]
# )
# 
# summary(res2$pval)
# mean(is.na(res2$pval))
# sum(res2$pval < 0.05, na.rm = TRUE)
# 
# truth <- rep(FALSE, n_genes)
# truth[de_genes] <- TRUE
# 
# power <- mean(res2$pval[truth] < 0.05)
# type1 <- mean(res2$pval[!truth] < 0.05)
# 
# c(power = power, type1 = type1)