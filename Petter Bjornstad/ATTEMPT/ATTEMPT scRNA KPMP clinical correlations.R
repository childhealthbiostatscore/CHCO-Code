## ---- setup ----
options(stringsAsFactors = FALSE)
.libPaths(c("/mmfs1/gscratch/togo/R", .libPaths()))

suppressPackageStartupMessages({
  library(arsenal); library(Biobase); library(BiocGenerics); library(BiocParallel)
  library(broom.mixed); library(colorspace); library(cowplot); library(data.table)
  library(DirichletReg); library(dplyr); library(edgeR); library(emmeans); library(enrichR)
  library(foreach); library(future); library(future.apply); library(GSEABase)
  library(ggdendro); library(ggpubr); library(glmmTMB); library(harmony); library(jsonlite)
  library(kableExtra); library(limma); library(MAST); library(Matrix); library(msigdbr)
  library(muscat); library(NMF); library(nebula); library(patchwork); library(pheatmap)
  library(readxl); library(REDCapR); library(reshape2); library(rstatix); library(SAVER)
  library(scater); library(scran); library(Seurat); library(SingleCellExperiment)
  library(slingshot); library(tidyverse); library(UpSetR); library(WriteXLS)
  library(parallel); library(doParallel)
  library(reticulate)
})

ts <- function(...){ message(format(Sys.time(), "%F %T"), " — ", sprintf(...)) }
flush_now <- function(){ try(flush.console(), silent = TRUE) }

ts("START: ATTEMPT scRNA KPMP clinical correlations driver")

## ---- cores / workers ----
get_allocated_cores <- function(){
  x <- Sys.getenv("SLURM_CPUS_PER_TASK", "")
  if(nzchar(x)) return(as.integer(x))
  parallel::detectCores(logical = TRUE)
}
alloc_cores <- get_allocated_cores()
ts("Cores: SLURM_CPUS_PER_TASK=%s | detectCores=%s", Sys.getenv("SLURM_CPUS_PER_TASK",""), parallel::detectCores(TRUE))
flush_now()

## You asked for 100; cap at allocation to avoid oversubscription
requested_workers <- 100L
workers_effective <- max(1L, min(requested_workers, alloc_cores))
ts("Workers requested=%d | effective=%d", requested_workers, workers_effective)

## BLAS threading should be 1 (you also set in sbatch env)
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

## ---- functions & config ----
source("/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/ATTEMPT/attempt_functions.R")
ts("Sourced attempt_functions.R")

## ---- reticulate / python env ----
reticulate::use_python("/mmfs1/gscratch/togo/yejichoi/mypy/bin/python", required = TRUE)
ts("Python config:\n%s", capture.output(reticulate::py_config()) %>% paste(collapse="\n"))

if(!reticulate::py_module_available("boto3")) stop("Python module 'boto3' not available in /mypy env")
if(!reticulate::py_module_available("pandas")) ts("WARNING: pandas not found; continuing without it")

boto3 <- reticulate::import("boto3")
# pd <- if (reticulate::py_module_available("pandas")) reticulate::import("pandas") else NULL

## ---- S3 session ----
keys_path <- "/mmfs1/home/yejichoi/keys.json"
stopifnot(file.exists(keys_path))
keys <- jsonlite::fromJSON(keys_path)
stopifnot(all(c("MY_ACCESS_KEY","MY_SECRET_KEY") %in% names(keys)))

session <- boto3$session$Session(
  aws_access_key_id = keys$MY_ACCESS_KEY,
  aws_secret_access_key = keys$MY_SECRET_KEY
)
s3 <- session$client("s3", endpoint_url = "https://s3.kopah.uw.edu")
ts("S3 client initialized (kopah endpoint)")

## ---- read seurat ----
bucket <- "attempt"
key <- "cleaned_data/attempt_clean_so.rds"
tmp <- tempfile(fileext = ".rds")

ts("Downloading Seurat object: s3://%s/%s", bucket, key); flush_now()
tryCatch({
  s3$download_file(bucket, key, tmp)
}, error = function(e){
  stop(sprintf("S3 download failed for s3://%s/%s: %s", bucket, key, conditionMessage(e)))
})
ts("Download complete → %s", tmp); flush_now()

attempt_so <- readRDS(tmp)
unlink(tmp)
ts("Seurat readRDS complete: cells=%d, genes=%d",
   ncol(attempt_so), nrow(attempt_so))
flush_now()

## ---- groups / clinical vars ----
celltype_groups <- list(
  PT           = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL          = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC           = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  EC           = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC"),
  IC           = c("IC-A", "IC-B", "aIC"),
  Immune       = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK"),
  Immune_myeloid  = c("MAC", "MON", "cDC", "pDC"),
  Immune_lymphoid = c("CD4+ T", "CD8+ T", "B", "NK"),
  VSMC_P_FIB   = c("VSMC/P", "FIB"),
  POD          = "POD"
)

## sanity: show which groups have at least one matching cell in metadata
ct_var <- "KPMP_celltype"
present <- lapply(celltype_groups, function(v) intersect(v, unique(attempt_so@meta.data[[ct_var]])))
missing <- lapply(celltype_groups, function(v) setdiff(v, unique(attempt_so@meta.data[[ct_var]])))
ts("Celltype var='%s' — present per group:\n%s",
   ct_var,
   paste(sprintf("  - %s: %s", names(present), vapply(present, function(x) paste(x, collapse=", "), "")),
         collapse = "\n"))
if (any(lengths(missing) > 0L)) {
  ts("NOTE: Missing labels per group (ignored):\n%s",
     paste(sprintf("  - %s: %s", names(missing), vapply(missing, function(x) paste(x, collapse=", "), "")),
           collapse = "\n"))
}
flush_now()

clin_vars <- c(
  "mgfr_jodal_delta", "mgfr_jodal_bsa_delta", "tir_delta",
  "hba1c_delta", "weight_delta", "avg_c_r2_delta", "avg_k_r2_delta",
  "avg_ketones_delta"
)
ts("Clinical variables: %s", paste(clin_vars, collapse=", ")); flush_now()

## ---- run ----
ts("Launching run_nebula_by_groups with workers=%d", workers_effective); flush_now()

results_grid <- run_nebula_by_groups(
  seurat_obj     = attempt_so,
  celltype_groups= celltype_groups,
  celltype_var   = ct_var,
  clin_vars      = clin_vars,
  suffix         = "kpmp",
  n_hvgs         = 2000,
  workers        = workers_effective,
  s3             = s3,
  s3_bucket      = "attempt",
  s3_key_prefix  = "associations/nebula",
  assay_layer    = "counts",
  id_var         = "subject_id",
  offset_var     = "pooled_offset",
  treatment_var  = "treatment",
  extra_covars   = NULL
)

ts("ALL DONE. results_grid keys recorded for %d groups.", length(results_grid))
flush_now()