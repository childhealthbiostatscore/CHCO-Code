library(tidyverse)
library(Seurat)
library(glmmTMB)
library(parallel)
setwd("/gscratch/togo/tvigers")
# Import annotated objects
biopsy = readRDS(
  "./Data/UWMDI/Hailey Hampson/Kidney scRNAseq Project/Data_Clean/cleaned_biopsy.rds"
)
organoid = readRDS(
  "./Data/UWMDI/Hailey Hampson/Kidney scRNAseq Project/Data_Clean/cleaned_organoid.rds"
)
# Create model matrices
biopsy_meta = biopsy@meta.data |>
  select(record_id, group) |>
  mutate(sample_type = "Biopsy")
organoid_meta = organoid@meta.data |>
  select(record_id, group) |>
  mutate(sample_type = "Organoid")
# Combine
mod_mat = rbind(biopsy_meta, organoid_meta)
mod_mat$group = factor(
  mod_mat$group,
  levels = c(
    "Lean_Control",
    "Type_2_Diabetes",
    "Type 2 Diabetes",
    "Lean Control"
  ),
  labels = c(
    "Lean Control",
    "Type 2 Diabetes",
    "Type 2 Diabetes",
    "Lean Control"
  )
)
# Create outcome matrices
biopsy_counts = t(as.matrix(GetAssayData(object = biopsy, layer = "data")))
organoid_counts = t(as.matrix(GetAssayData(object = organoid, layer = "data")))
# Combine
gene_overlap = intersect(colnames(organoid_counts), colnames(biopsy_counts))
counts = rbind(biopsy_counts[, gene_overlap], organoid_counts[, gene_overlap])
# Out everything together and clear memory a bit
df = cbind(mod_mat, counts)
rm(biopsy_meta, organoid_meta, biopsy_counts, organoid_counts, counts, mod_mat)
# Set up parallel processing
# Register cluster
cl <- makeCluster(24, type = "FORK")
# For each gene measured in both biopsies and organoids, fit a mixed model.
# Use the same ZI formula for all models (intercept only for now)
zf = ~1
model_results = parLapply(cl, gene_overlap, function(gene) {
  f1 = as.formula(paste0(
    gene,
    " ~ group * sample_type + (1 | record_id/sample_type)"
  ))
  model = tryCatch(
    glmmTMB(f1, data = df, family = gaussian, ziformula = zf),
    error = function(e) return(NA),
    warning = function(w) return(NA)
  )
  return(model)
})
# Stop the cluster
stopCluster(cl)
# Save results
save(
  model_results,
  "./Data/UWMDI/Hailey Hampson/Kidney scRNAseq Project/Results/model_results.RData"
)
