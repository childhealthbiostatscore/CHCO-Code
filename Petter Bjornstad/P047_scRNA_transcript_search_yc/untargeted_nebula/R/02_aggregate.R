#!/usr/bin/env Rscript
# =============================================================================
# 02_aggregate.R  —  combine partial CSVs into ONE CSV per contrast.
# Reads all partials under <partials_prefix>/<contrast>/ and row-binds them so
# every subset + cell type + resolution for a contrast lands in one file:
#   <results_prefix>/<contrast>_neb.csv
# Optional arg: a single contrast name (default = all contrasts in config).
# =============================================================================

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "R"
source(file.path(here, "setup.R"))
source(file.path(dirname(here), "config.R"))
source(file.path(here, "nebula_core.R"))   # for s3_get_csv / s3_put_csv

g <- CONFIG$global

args      <- commandArgs(trailingOnly = TRUE)
contrasts <- if (length(args) >= 1) args else names(CONFIG$contrasts)

list_partials <- function(contrast) {
  prefix <- paste0(g$partials_prefix, contrast, "/")
  b <- tryCatch(get_bucket(bucket = g$s3_bucket, prefix = prefix,
                           region = "", max = Inf),
                error = function(e) NULL)
  if (is.null(b)) return(character(0))
  keys <- vapply(b, function(x) x$Key %||% "", character(1))
  keys[grepl("\\.csv$", keys)]
}

for (contrast in contrasts) {
  keys <- list_partials(contrast)
  if (length(keys) == 0) { message("[", contrast, "] no partials found — skipping"); next }

  message("[", contrast, "] combining ", length(keys), " partials")
  dfs <- lapply(keys, function(k)
    s3_get_csv(k, g$s3_bucket, check.names = FALSE, stringsAsFactors = FALSE))
  combined <- dplyr::bind_rows(dfs)   # tolerates differing per-group count cols

  # tidy column order: descriptors first, then everything else
  lead <- c("contrast", "var", "type", "dataset", "subset", "subset_label",
            "resolution", "celltype_col", "celltype", "reference_group",
            "model_formula", "n_genes_tested", "gene")
  lead <- lead[lead %in% names(combined)]
  combined <- combined[, c(lead, setdiff(names(combined), lead))]

  out_key <- paste0(g$results_prefix, contrast, "_neb.csv")
  s3_put_csv(combined, out_key, g$s3_bucket)
  message("[", contrast, "] wrote ", nrow(combined), " rows -> ", out_key)
}

message("=== aggregation complete ===")
