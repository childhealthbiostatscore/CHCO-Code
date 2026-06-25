# =============================================================================
# setup.R  —  libraries, user-aware paths, S3 credentials + helpers
# Sourced by every script. Does NOT load Seurat objects.
# =============================================================================

# The container's site-library is read-only. The SLURM scripts pass a writable
# library path via the NEB_RLIB env var (custom name so the container's
# Renviron can't clobber it, unlike R_LIBS_USER). Put it FIRST on .libPaths so
# both installs and library() use it.
.userlib <- Sys.getenv("NEB_RLIB")
if (nzchar(.userlib)) {
  dir.create(.userlib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(.userlib, .libPaths()))
}

.ensure_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", repos = "https://cloud.r-project.org")
    message("Installing ", repo, " into ", .libPaths()[1], " ...")
    remotes::install_github(repo, lib = .libPaths()[1], upgrade = "never")
  }
}
.ensure_github("togolab",       "uwmdi-togo/togolab")
.ensure_github("togolab.scrna", "uwmdi-togo/togolab.scrna")

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(aws.s3)
  library(jsonlite)
  library(nebula)
  library(Matrix)
  library(stringr)
  library(tidyr)
  library(SingleCellExperiment)
  library(togolab)
  library(togolab.scrna)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

togo_paths()

options(future.globals.maxSize = 3000 * 1024^3)

# sanitize a label into a filesystem/S3-safe token
safe_token <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
