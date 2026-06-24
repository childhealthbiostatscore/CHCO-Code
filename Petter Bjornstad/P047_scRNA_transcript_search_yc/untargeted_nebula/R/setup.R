# =============================================================================
# setup.R  —  libraries, user-aware paths, S3 credentials + helpers
# Sourced by every script. Does NOT load Seurat objects.
# =============================================================================

.ensure_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", repos = "https://cloud.r-project.org")
    message("Installing ", repo, " ...")
    remotes::install_github(repo, upgrade = "never")
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
