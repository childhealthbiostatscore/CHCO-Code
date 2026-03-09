suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDesign3)
  library(dplyr)
  library(Matrix)
  library(BiocParallel)
  library(optparse)
  library(tibble)
  library(aws.s3)
  library(DuoClustering2018)
})

# ── S3 / Multi-user setup ─────────────────────────────────────────────────────
# Detect user and set up AWS credentials + endpoint
setup_s3 <- function() {
  user <- Sys.info()[["user"]]
  
  if (user == "choiyej") {
    # --- local (choiyej Mac) ---
    keys_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json"
  } else if (user %in% c("rameshsh", "yejichoi", "pylell")) {
    # --- Hyak HPC ---
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

# pull sample dataset from scDesign3 tutorial (https://songdongyuan1994.github.io/scDesign3/docs/articles/scDesign3.html)
# data from (url("https://figshare.com/ndownloader/files/40581992"))
example_sce <- s3readRDS(bucket = "scrna", object = "Projects/Paired scRNA simulation analysis/example/PANCREAS_sce.rds",
                         region = "")
print(example_sce)

# for computational time, just using 100 genes
example_sce <- example_sce[1:100, ]

# simulation
set.seed(123)
example_simu <- scdesign3(
  sce = example_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = "pseudotime",
  spatial = NULL,
  other_covariates = NULL,
  mu_formula = "s(pseudotime, k = 10, bs = 'cr')",
  sigma_formula = "1", # If you want your dispersion also varies along pseudotime, use "s(pseudotime, k = 5, bs = 'cr')"
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE
)

dim(example_simu$new_count)
# 100 x 2087

# construct sce
logcounts(example_sce) <- log1p(counts(example_sce))
simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count), colData = example_simu$new_covariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))

# add in pseudo group and visit vars
# Create balanced groupings
group <- rep(c("groupA", "groupA", "groupB", "groupB"), length.out = n_cells)
visit <- rep(c("visit1", "visit2", "visit1", "visit2"), length.out = n_cells)

# Shuffle so cells are randomly assigned (not in order)
set.seed(42)
n_cells <- ncol(simu_sce)
idx <- sample(n_cells)
group <- group[order(idx)]
visit <- visit[order(idx)]

# Add to colData
simu_sce$group <- group
simu_sce$visit <- visit

# Verify
table(simu_sce$group, simu_sce$visit)

# run nebula
library(nebula)

# 1. Extract count matrix and metadata
counts <- counts(simu_sce)  # or assay(simu_sce, "counts")
meta <- as.data.frame(colData(simu_sce))

# 2. Create subject IDs (required by NEBULA — needs clustering unit)
# Since this is simulated, create fake subject IDs
# e.g., 10 subjects, cells distributed across them
set.seed(42)
# Create subject IDs where each subject belongs to one group, with both visits
n_subjects_per_group <- 5  # 5 subjects per group, 10 total

subject_group_map <- data.frame(
  subject_id = paste0("subj", 1:10),
  group = rep(c("groupA", "groupB"), each = n_subjects_per_group)
)

# Assign subject IDs consistent with their group
meta$subject_id <- NA

for (subj in subject_group_map$subject_id) {
  grp <- subject_group_map$group[subject_group_map$subject_id == subj]
  # Find cells belonging to this group not yet assigned
  available <- which(meta$group == grp & is.na(meta$subject_id))
  # Assign roughly equal cells per subject
  n_assign <- round(sum(meta$group == grp) / n_subjects_per_group)
  chosen <- sample(available, min(n_assign, length(available)))
  meta$subject_id[chosen] <- subj
}

# Fill any unassigned cells (due to rounding)
unassigned <- which(is.na(meta$subject_id))
for (i in unassigned) {
  grp <- meta$group[i]
  possible_subjs <- subject_group_map$subject_id[subject_group_map$group == grp]
  meta$subject_id[i] <- sample(possible_subjs, 1)
}

# Verify — each subject should only appear in one group
table(meta$subject_id, meta$group)
table(meta$subject_id, meta$visit)  # each subject should have both visits

# 3. Aggregate data by subject (NEBULA's required input format)
data_g <- group_cell(count = counts,
                     id    = meta$subject_id,
                     pred  = meta)

# 4. Build design matrix with group * visit interaction
pred <- model.matrix(~ group * visit, data = meta)

# 5. Run NEBULA
res <- nebula(count = data_g$count,
              id    = data_g$id,
              pred  = pred,
              model = "NBLMM")  # or "NBGMM"

# 6. View results
head(res$summary)
