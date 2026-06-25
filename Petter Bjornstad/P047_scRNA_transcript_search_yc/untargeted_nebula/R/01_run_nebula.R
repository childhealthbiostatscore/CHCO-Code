#!/usr/bin/env Rscript
# =============================================================================
# 01_run_nebula.R  —  ONE NEBULA run (one SLURM array task)
# Usage: Rscript 01_run_nebula.R <row_index>
#   <row_index> = SLURM_ARRAY_TASK_ID (1-based line into jobs.tsv).
# Loads the pre-split cell-type object, applies the contrast's subset, fits
# untargeted NEBULA on all sufficiently-expressed genes, writes ONE partial CSV
# (contrast x subset x resolution x cell type) to S3. 02_aggregate.R then
# stitches all partials for a contrast into a single CSV.
# =============================================================================

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "R"
source(file.path(here, "setup.R"))
source(file.path(dirname(here), "config.R"))
source(file.path(here, "nebula_core.R"))

g <- CONFIG$global

args    <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1] %||% Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(task_id)) stop("Provide a row index (or set SLURM_ARRAY_TASK_ID)")
# chunked submission passes a row offset so we can exceed SLURM MaxArraySize
offset <- as.integer(Sys.getenv("NEB_ROW_OFFSET", "0"))
if (is.na(offset)) offset <- 0L
idx    <- offset + task_id

ncore <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", g$nebula$ncore))
if (is.na(ncore) || ncore < 1) ncore <- g$nebula$ncore

manifest <- read.delim(file.path(g$work_dir, "jobs.tsv"), stringsAsFactors = FALSE)
if (idx < 1 || idx > nrow(manifest)) stop("row index ", idx, " out of range (1..", nrow(manifest), ")")
job <- manifest[idx, ]

message("=== run ", idx, "/", nrow(manifest), " : ", job$contrast, " | ",
        job$subset, " | ", job$resolution, " | ", job$celltype, " ===")

contrast <- CONFIG$contrasts[[job$contrast]]
contrast$name <- job$contrast
dcfg     <- CONFIG$datasets[[job$dataset]]
sspec    <- contrast$subsets[[job$subset]]

# load the pre-split cell-type object
obj <- s3readRDS(bucket = job$split_bucket, object = job$split_object, region = "")

out <- run_one(
  obj             = obj,
  contrast        = contrast,
  subset_key      = job$subset,
  subset_spec     = sspec,
  resolution_name = job$resolution,
  ct_col          = job$celltype_col,
  celltype        = job$celltype,
  dataset_cfg     = dcfg,
  nebula_cfg      = g$nebula,
  ncore           = ncore,
  dataset_name    = job$dataset
)

if (is.null(out)) {
  message("=== run ", idx, " produced no output (skipped) ===")
  quit(save = "no", status = 0)
}

partial_key <- paste0(g$partials_prefix, job$contrast, "/",
                      job$subset, "__", job$resolution, "__",
                      job$celltype_safe, ".csv")
s3_put_csv(out, partial_key, g$s3_bucket)
message("=== saved partial -> ", partial_key, " (", nrow(out), " genes) ===")
