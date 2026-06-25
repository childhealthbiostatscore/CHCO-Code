#!/usr/bin/env Rscript
# =============================================================================
# 00_split_celltypes.R
# For every dataset referenced by a contrast:
#   1. load the full Seurat object from S3
#   2. derive needed metadata columns (DKD groups, KPMP_celltype_general, ...)
#   3. split into one object per cell type, at BOTH resolutions
#      (general / low-res and specific / high-res), save each to S3
# Then write the job manifest (jobs.tsv): one line per
#   contrast x subset x resolution x cell type  -> one NEBULA run.
# Splitting once here means each NEBULA job loads only its (small) cell-type
# object instead of the whole dataset.
# =============================================================================

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "R"
source(file.path(here, "setup.R"))
source(file.path(dirname(here), "config.R"))
source(file.path(here, "nebula_core.R"))

g   <- CONFIG$global
res <- CONFIG$resolutions
dir.create(g$work_dir, showWarnings = FALSE, recursive = TRUE)

# datasets actually needed (contrast default + any per-subset override)
needed_datasets <- unique(unlist(lapply(CONFIG$contrasts, function(c) {
  c(c$dataset, vapply(c$subsets, function(s) s$dataset %||% c$dataset, character(1)))
})))
needed_datasets <- needed_datasets[!is.na(needed_datasets)]

# celltypes[[dataset]][[resolution_name]] <- character vector of cell types
celltypes <- list()

for (ds in needed_datasets) {
  dcfg <- CONFIG$datasets[[ds]]
  message("\n=== Loading dataset: ", ds, " (", dcfg$s3_object, ") ===")
  obj <- s3readRDS(bucket = dcfg$s3_bucket, object = dcfg$s3_object, region = "")
  obj <- derive_meta(obj, ds, res, g$nebula)

  celltypes[[ds]] <- list()
  for (rname in names(res)) {
    ct_col <- res[[rname]]
    levs   <- sort(unique(as.character(obj@meta.data[[ct_col]])))
    levs   <- levs[!is.na(levs) & levs != "" & levs != "NA"]
    celltypes[[ds]][[rname]] <- levs
    message("  resolution '", rname, "' (", ct_col, "): ", length(levs), " cell types")

    for (ct in levs) {
      ct_safe <- safe_token(ct)
      keep    <- which(obj@meta.data[[ct_col]] %in% ct)
      sub     <- subset(obj, cells = colnames(obj)[keep])
      key     <- paste0(g$splits_prefix, ds, "/", rname, "/", ct_safe, ".rds")
      s3saveRDS(sub, object = key, bucket = g$s3_bucket, region = "", multipart = TRUE)
      message("    saved [", ncol(sub), " cells] -> ", key)
    }
  }
  rm(obj); gc()
}

# ---- build the manifest -----------------------------------------------------
rows <- list()
jid  <- 0L
for (cname in names(CONFIG$contrasts)) {
  ct <- CONFIG$contrasts[[cname]]
  for (skey in names(ct$subsets)) {
    ds <- ct$subsets[[skey]]$dataset %||% ct$dataset   # subset may override dataset
    for (rname in names(res)) {
      ct_col <- res[[rname]]
      for (celltype in celltypes[[ds]][[rname]]) {
        jid <- jid + 1L
        rows[[length(rows) + 1L]] <- data.frame(
          job_id        = jid,
          dataset       = ds,
          contrast      = cname,
          subset        = skey,
          resolution    = rname,
          celltype_col  = ct_col,
          celltype      = celltype,
          celltype_safe = safe_token(celltype),
          split_bucket  = g$s3_bucket,
          split_object  = paste0(g$splits_prefix, ds, "/", rname, "/",
                                  safe_token(celltype), ".rds"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
manifest <- do.call(rbind, rows)

manifest_path <- file.path(g$work_dir, "jobs.tsv")
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
# also keep a copy on S3 for the record (best-effort: the local jobs.tsv above
# is what the array uses, so a failed upload here must not abort the pipeline)
tryCatch(
  s3_put_tsv(manifest, paste0(g$results_prefix, "jobs_manifest.tsv"), g$s3_bucket),
  error = function(e) message("  (manifest S3 upload skipped: ", conditionMessage(e), ")"))

message("\n=== Manifest written: ", nrow(manifest), " runs -> ", manifest_path, " ===")
