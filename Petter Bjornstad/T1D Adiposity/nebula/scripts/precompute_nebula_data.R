#!/usr/bin/env Rscript
# =============================================================================
# Precompute per-cell-type NEBULA-ready data objects
# =============================================================================
# Run this ONCE after save_celltype_subsets.R. For each cell-type subset on
# S3, this loads the Seurat subset, derives all binary obesity variables that
# any analysis_type in run_nebula_single.R might need, extracts the count
# matrix + the analysis-relevant metadata columns, computes the raw
# library-size offset (nebula takes the log of this internally — per nebula
# docs), and writes a lean list object to:
#
#   s3://t1d.adiposity/data_clean/nebula_ready/
#       t1d_adiposity_nebula_<celltype>.rds
#
# The precomputed object structure (mirrors what run_nebula_single.R expects
# to feed into nebula()):
#   list(
#     celltype       = <character, original celltype label>,
#     celltype_var   = <"KPMP_celltype" | "KPMP_celltype_general">,
#     count          = dgCMatrix, genes x cells (raw counts, no gene filter),
#     meta           = data.frame, one row per cell (lean column subset),
#     id             = character, record_id per cell,
#     offset         = numeric, colSums(count) — raw library size (NOT log),
#     n_cells        = ncol(count),
#     n_subjects     = length(unique(id)),
#     generated_at   = Sys.time()
#   )
#
# Design notes:
#   - No gene-level filtering applied here. nebula's internal QC
#     (cutoff_count, cpc) does the filtering at run time, matching the
#     current run_nebula_single.R behavior.
#   - Offset is raw library size. nebula logs it internally
#     (https://cran.r-project.org/web/packages/nebula/nebula.pdf).
#   - Binary obesity variables are derived once here so every analysis
#     type can use them without re-deriving per task.
#   - Metadata is trimmed to the columns actually used by any analysis in
#     run_nebula_single.R, dropping the Seurat machinery bulk (assays,
#     reductions, graphs) to shrink the S3 payload.
# =============================================================================

library(aws.s3)
library(jsonlite)
library(Seurat)
library(Matrix)
library(dplyr)

# -----------------------------------------------------------------------------
# AWS / S3
# -----------------------------------------------------------------------------
keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
Sys.setenv(
  "AWS_ACCESS_KEY_ID"     = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION"    = "",
  "AWS_REGION"            = "",
  "AWS_S3_ENDPOINT"       = "s3.kopah.uw.edu"
)
s3_bucket <- "t1d.adiposity"

# -----------------------------------------------------------------------------
# Cell types to precompute
# -----------------------------------------------------------------------------
# Focus on the main lineages of interest for this pass; unused lineages are
# commented out so they're easy to re-enable.

kpmp_celltypes <- c(
  # PT subtypes
  "PT-S1/S2", "PT-S3", "aPT",
  # TAL subtypes
  "C-TAL-1", "C-TAL-2", "aTAL", "dTAL",
  # PC subtypes
  "CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC",
  # IC subtypes
  "IC-A", "IC-B", "aIC",
  # EC subtypes
  "EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A",
  # Immune subtypes (myeloid + lymphoid)
  "MAC", "MON", "cDC", "pDC",
  "CD4+ T", "CD8+ T", "B", "NK", "cycT"

  # ---- Not precomputed this pass (uncomment to re-enable) ----
  # , "DTL", "aDTL", "ATL"
  # , "DCT", "dDCT", "CNT"
  # , "VSMC/P", "FIB"
  # , "POD", "MC", "PEC"
  # , "SchwannCells", "non-specific"
)

kpmp_general <- c(
  "PT", "TAL", "PC", "IC", "EC", "Immune"

  # ---- Not precomputed this pass ----
  # , "DTL_ATL", "DCT_CNT"
  # , "VSMC_P_FIB"
  # , "POD", "MC", "PEC"
  # , "Schwann", "Other"
)

# -----------------------------------------------------------------------------
# Metadata columns any analysis_type may need. Missing columns are skipped
# with a warning (printed but not fatal) so this script tolerates schema
# drift across Seurat updates.
# -----------------------------------------------------------------------------
keep_meta_cols <- c(
  # Identifiers
  "record_id", "barcode", "source", "kit_id", "visit",
  "KPMP_celltype", "KPMP_celltype_general",
  # Grouping variables used by analysis_config
  "group", "study", "age", "sex",
  # BMI-based
  "bmi", "bmi_obesity",
  # DXA-based
  "dxa_obesity",
  # Continuous DXA variables (used by cont_dexa_* analyses)
  "dexa_body_fat", "dexa_bone_mineral_density", "dexa_fat_kg",
  "dexa_lean_mass", "dexa_lean_kg", "dexa_ag_ratio",
  "dexa_est_vat", "dexa_trunk_kg", "dexa_trunk_mass",
  # Step 6 matched-subset flag (unused this pass but cheap to keep)
  "has_both_bmi_dxa"
)

# -----------------------------------------------------------------------------
# Helper: precompute one cell-type subset
# -----------------------------------------------------------------------------
precompute_one <- function(celltype, celltype_var) {
  safe_name <- gsub("/", "_", celltype)
  in_key  <- sprintf("data_clean/subset/t1d_adiposity_subset_%s.rds", safe_name)
  out_key <- sprintf("data_clean/nebula_ready/t1d_adiposity_nebula_%s.rds",
                     safe_name)

  cat(sprintf("  Loading %s ... ", in_key))
  t0 <- Sys.time()
  obj <- tryCatch(
    s3readRDS(object = in_key, bucket = s3_bucket, region = ""),
    error = function(e) { cat(sprintf("SKIP (%s)\n", e$message)); NULL }
  )
  if (is.null(obj)) return(invisible(NULL))
  cat(sprintf("%d cells in %.0fs\n", ncol(obj),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  meta <- obj@meta.data

  # Derive all binary obesity variables upfront
  meta$bmi_obese_binary <- dplyr::case_when(
    meta$bmi_obesity == "Obese"                     ~ "Obese",
    meta$bmi_obesity %in% c("Normal", "Overweight") ~ "Non_Obese",
    TRUE                                            ~ NA_character_
  )
  meta$bmi_ow_obese_binary <- dplyr::case_when(
    meta$bmi_obesity == "Normal"                    ~ "Normal",
    meta$bmi_obesity %in% c("Overweight", "Obese")  ~ "Overweight_Obese",
    TRUE                                            ~ NA_character_
  )
  meta$dxa_obese_binary <- dplyr::case_when(
    meta$dxa_obesity == "Obese"                     ~ "Obese",
    meta$dxa_obesity %in% c("Normal", "Overweight") ~ "Non_Obese",
    TRUE                                            ~ NA_character_
  )
  meta$dxa_ow_obese_binary <- dplyr::case_when(
    meta$dxa_obesity == "Normal"                    ~ "Normal",
    meta$dxa_obesity %in% c("Overweight", "Obese")  ~ "Overweight_Obese",
    TRUE                                            ~ NA_character_
  )

  # Trim metadata to the columns downstream analyses might need
  wanted <- c(keep_meta_cols,
              "bmi_obese_binary", "bmi_ow_obese_binary",
              "dxa_obese_binary", "dxa_ow_obese_binary")
  present <- intersect(wanted, colnames(meta))
  missing <- setdiff(wanted, present)
  if (length(missing) > 0) {
    cat(sprintf("  WARN missing cols for %s: %s\n", celltype,
                paste(missing, collapse = ", ")))
  }
  meta_lean <- meta[, present, drop = FALSE]
  rownames(meta_lean) <- rownames(meta)

  # Counts + offset. offset = colSums(counts) — nebula logs internally.
  counts  <- GetAssayData(obj, assay = "RNA", layer = "counts")
  libsize <- Matrix::colSums(counts)

  out <- list(
    celltype     = celltype,
    celltype_var = celltype_var,
    count        = counts,              # dgCMatrix, genes x cells
    meta         = meta_lean,           # data.frame, one row per cell
    id           = as.character(meta$record_id),
    offset       = as.numeric(libsize), # raw library size; nebula logs it
    n_cells      = ncol(counts),
    n_subjects   = length(unique(meta$record_id)),
    generated_at = Sys.time()
  )

  # Sanity
  stopifnot(ncol(out$count) == nrow(out$meta))
  stopifnot(length(out$id) == ncol(out$count))
  stopifnot(length(out$offset) == ncol(out$count))

  cat(sprintf("  Saving %s ... ", out_key))
  s3saveRDS(out, object = out_key, bucket = s3_bucket, region = "",
            multipart = TRUE)
  cat(sprintf("done (%d cells, %d subjects).\n",
              out$n_cells, out$n_subjects))

  rm(obj, counts, meta, meta_lean, out); gc(verbose = FALSE)
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
cat("\n=== Precomputing KPMP_celltype subsets ===\n")
for (ct in kpmp_celltypes) {
  cat(sprintf("[%s]\n", ct))
  precompute_one(ct, "KPMP_celltype")
}

cat("\n=== Precomputing KPMP_celltype_general subsets ===\n")
for (ct in kpmp_general) {
  cat(sprintf("[%s]\n", ct))
  precompute_one(ct, "KPMP_celltype_general")
}

cat("\nDone.\n")
