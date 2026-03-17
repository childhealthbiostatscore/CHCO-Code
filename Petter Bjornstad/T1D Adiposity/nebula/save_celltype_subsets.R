#!/usr/bin/env Rscript
# =============================================================================
# Save cell type-specific Seurat subsets to S3
# =============================================================================
# Run this ONCE before submitting SLURM array jobs.
# This pre-splits the Seurat object by cell type so each array job
# only needs to load a small subset instead of the full object.
#
# Saves to: s3://t1d.adiposity/data_clean/subset/t1d_adiposity/
# =============================================================================

library(aws.s3)
library(jsonlite)
library(Seurat)
library(dplyr)

# Set up AWS
keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

s3_bucket <- "t1d.adiposity"

# Load the full Seurat object
cat("Loading full Seurat object...\n")
pb90 <- s3readRDS(object = "data_clean/t1d_hc_scrna_w_clinical.rds",
                  bucket = s3_bucket, region = "")
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(pb90), nrow(pb90)))

# Cell type groupings for KPMP_celltype_general
celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),
  DCT_CNT = c("DCT", "dDCT", "CNT"),
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A"),
  Immune_Myeloid = c("MAC", "MON", "cDC", "pDC"),
  Immune_Lymphoid = c("CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",
  PEC = "PEC",
  Schwann = "SchwannCells",
  Other = c("non-specific")
)

# Add KPMP_celltype_general if not present
map_celltype_to_general <- function(celltype, celltype_groups) {
  for (group_name in names(celltype_groups)) {
    if (celltype %in% celltype_groups[[group_name]]) {
      return(group_name)
    }
  }
  return("Other")
}

if (!"KPMP_celltype_general" %in% colnames(pb90@meta.data)) {
  cat("Adding KPMP_celltype_general...\n")
  pb90$KPMP_celltype_general <- sapply(pb90$KPMP_celltype,
                                        map_celltype_to_general,
                                        celltype_groups = celltype_groups)
}

# Save KPMP_celltype subsets
cat("\n--- Saving KPMP_celltype subsets ---\n")
kpmp_types <- unique(pb90$KPMP_celltype)
cat(sprintf("Found %d unique KPMP_celltype values\n", length(kpmp_types)))

for (ct in kpmp_types) {
  cat(sprintf("  Subsetting %s... ", ct))
  sub <- subset(pb90, KPMP_celltype == ct)
  n <- ncol(sub)
  cat(sprintf("%d cells. ", n))

  # Use gsub to make filename safe (replace / with _)
  safe_name <- gsub("/", "_", ct)
  s3_key <- paste0("data_clean/subset/t1d_adiposity/t1d_adiposity_subset_", safe_name, ".rds")
  s3saveRDS(sub, object = s3_key, bucket = s3_bucket, region = "", multipart = TRUE)
  cat(sprintf("Saved to %s\n", s3_key))

  rm(sub); gc()
}

# Save KPMP_celltype_general subsets
cat("\n--- Saving KPMP_celltype_general subsets ---\n")
general_types <- unique(pb90$KPMP_celltype_general)
cat(sprintf("Found %d unique KPMP_celltype_general values\n", length(general_types)))

for (ct in general_types) {
  cat(sprintf("  Subsetting %s... ", ct))
  sub <- subset(pb90, KPMP_celltype_general == ct)
  n <- ncol(sub)
  cat(sprintf("%d cells. ", n))

  safe_name <- gsub("/", "_", ct)
  s3_key <- paste0("data_clean/subset/t1d_adiposity/t1d_adiposity_subset_", safe_name, ".rds")
  s3saveRDS(sub, object = s3_key, bucket = s3_bucket, region = "", multipart = TRUE)
  cat(sprintf("Saved to %s\n", s3_key))

  rm(sub); gc()
}

cat(sprintf("\nDone! Saved %d + %d = %d cell type subsets\n",
            length(kpmp_types), length(general_types),
            length(kpmp_types) + length(general_types)))
cat("These will be loaded by run_nebula_single.R for each array job.\n")
