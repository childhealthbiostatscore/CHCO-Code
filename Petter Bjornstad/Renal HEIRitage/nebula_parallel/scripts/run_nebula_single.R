#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript run_nebula_single.R <analysis_type> <cell_type>")
}

analysis_type <- args[1]
cell_type_group <- args[2]

# Load libraries
library(aws.s3)
library(jsonlite)
library(biomaRt)
library(Seurat)
library(dplyr)

# Log start time
cat(sprintf("Starting analysis: %s for cell type group: %s at %s\n", 
            analysis_type, cell_type_group, Sys.time()))

# Set up user environment
user <- "yejichoi"  # HPC username
root_path <- "/mmfs1/gscratch/togo/yejichoi/"
git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"

# Set up AWS environment for Kopah
keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# Load data and functions
cat("Loading data...\n")
pb90_subset_clean <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset.rds", 
                               bucket = "scrna", region = "")
rh_rh2_croc_improve_unique <- s3readRDS(object = "data_clean/rh_rh2_croc_improve_unique.RDS", 
                                        bucket = "harmonized.dataset", region = "")

source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_scRNA_functions.R"))
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_functions.R"))

# Initialize biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to run specific analysis
run_analysis <- function(analysis_type, cell_type_group) {
  gc()
  
  # Get the actual cell types for this group
  cell_types <- celltype_groups[[cell_type_group]]
  
  if (analysis_type == "DKD_vs_nonDKD_100") {
    # DKD vs non-DKD (ACR >= 100)
    pb90_celltype <- subset(pb90_subset_clean, 
                            KPMP_celltype_general %in% cell_types & 
                              group != "Lean_Control" & 
                              !is.na(pb90_subset_clean@meta.data$dkd_group_100))
    
    cat(sprintf("Sample counts for %s:\n", cell_type_group))
    print(table(pb90_celltype@meta.data %>%
                  distinct(record_id, .keep_all = T) %>%
                  pull(dkd_group_100)))
    
    nebula_res <- run_nebula_parallel(pb90_celltype,
                                      n_cores = 15,
                                      subject_id_col = "record_id",
                                      formula = ~ dkd_group_100,
                                      s3_bucket = "scrna",
                                      s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_vs_nonDKD_100/", 
                                                      cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_dkd100.rds"),
                                      group = F)
    
    processed <- process_nebula_results(nebula_res$results,
                                        pval_col = "p_dkd_group_100DKD",
                                        logfc_col = "logFC_dkd_group_100DKD")
    
  } else if (analysis_type == "DKD_100_vs_HC") {
    # DKD (ACR >= 100) vs HC
    pb90_celltype <- subset(pb90_subset_clean, 
                            KPMP_celltype_general %in% cell_types & 
                              dkd_group_100_hc != "non_DKD" &
                              !is.na(pb90_subset_clean@meta.data$dkd_group_100_hc))
    
    pb90_celltype$dkd_group_100_hc <- factor(pb90_celltype$dkd_group_100_hc)
    pb90_celltype$dkd_group_100_hc <- relevel(pb90_celltype$dkd_group_100_hc, ref = "HC")
    
    cat(sprintf("Sample counts for %s:\n", cell_type_group))
    print(table(pb90_celltype@meta.data %>%
                  distinct(record_id, .keep_all = T) %>%
                  pull(dkd_group_100_hc)))
    
    nebula_res <- run_nebula_parallel(pb90_celltype,
                                      n_cores = 15,
                                      subject_id_col = "record_id",
                                      formula = ~ dkd_group_100_hc,
                                      s3_bucket = "scrna",
                                      s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_100_vs_hc/", 
                                                      cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_dkd100_hc.rds"),
                                      group = F)
    
    processed <- process_nebula_results(nebula_res$results,
                                        pval_col = "p_dkd_group_100_hcDKD",
                                        logfc_col = "logFC_dkd_group_100_hcDKD")
    
  } else if (analysis_type == "nonDKD_100_vs_HC") {
    # non-DKD (ACR >= 100) vs HC
    pb90_celltype <- subset(pb90_subset_clean, 
                            KPMP_celltype_general %in% cell_types & 
                              dkd_group_100_hc != "DKD" &
                              !is.na(pb90_subset_clean@meta.data$dkd_group_100_hc))
    
    pb90_celltype$dkd_group_100_hc <- factor(pb90_celltype$dkd_group_100_hc)
    pb90_celltype$dkd_group_100_hc <- relevel(pb90_celltype$dkd_group_100_hc, ref = "HC")
    
    cat(sprintf("Sample counts for %s:\n", cell_type_group))
    print(table(pb90_celltype@meta.data %>%
                  distinct(record_id, .keep_all = T) %>%
                  pull(dkd_group_100_hc)))
    
    nebula_res <- run_nebula_parallel(pb90_celltype,
                                      n_cores = 15,
                                      subject_id_col = "record_id",
                                      formula = ~ dkd_group_100_hc,
                                      s3_bucket = "scrna",
                                      s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/nonDKD_100_vs_hc/", 
                                                      cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_nondkd100_hc.rds"),
                                      group = F)
    
    processed <- process_nebula_results(nebula_res$results,
                                        pval_col = "p_dkd_group_100_hcnon_DKD",
                                        logfc_col = "logFC_dkd_group_100_hcnon_DKD")
    
  } else if (analysis_type == "DKD_vs_nonDKD_30") {
    # DKD vs non-DKD (ACR >= 30)
    pb90_celltype <- subset(pb90_subset_clean, 
                            KPMP_celltype_general %in% cell_types & 
                              group != "Lean_Control" & 
                              !is.na(pb90_subset_clean@meta.data$dkd_group_30))
    
    cat(sprintf("Sample counts for %s:\n", cell_type_group))
    print(table(pb90_celltype@meta.data %>%
                  distinct(record_id, .keep_all = T) %>%
                  pull(dkd_group_30)))
    
    nebula_res <- run_nebula_parallel(pb90_celltype,
                                      n_cores = 15,
                                      subject_id_col = "record_id",
                                      formula = ~ dkd_group_30,
                                      s3_bucket = "scrna",
                                      s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_vs_nonDKD_30/", 
                                                      cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_dkd30.rds"),
                                      group = F)
    
    processed <- process_nebula_results(nebula_res$results,
                                        pval_col = "p_dkd_group_30DKD",
                                        logfc_col = "logFC_dkd_group_30DKD")
    
  } else if (analysis_type == "DKD_30_vs_HC") {
    # DKD (ACR >= 30) vs HC
    pb90_celltype <- subset(pb90_subset_clean, 
                            KPMP_celltype_general %in% cell_types & 
                              dkd_group_30_hc != "non_DKD" &
                              !is.na(pb90_subset_clean@meta.data$dkd_group_30_hc))
    
    pb90_celltype$dkd_group_30_hc <- factor(pb90_celltype$dkd_group_30_hc)
    pb90_celltype$dkd_group_30_hc <- relevel(pb90_celltype$dkd_group_30_hc, ref = "HC")
    
    cat(sprintf("Sample counts for %s:\n", cell_type_group))
    print(table(pb90_celltype@meta.data %>%
                  distinct(record_id, .keep_all = T) %>%
                  pull(dkd_group_30_hc)))
    
    nebula_res <- run_nebula_parallel(pb90_celltype,
                                      n_cores = 15,
                                      subject_id_col = "record_id",
                                      formula = ~ dkd_group_30_hc,
                                      s3_bucket = "scrna",
                                      s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_30_vs_hc/", 
                                                      cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_dkd30_hc.rds"),
                                      group = F)
    
    processed <- process_nebula_results(nebula_res$results,
                                        pval_col = "p_dkd_group_30_hcDKD",
                                        logfc_col = "logFC_dkd_group_30_hcDKD")
    
  } else if (analysis_type == "nonDKD_30_vs_HC") {
    # non-DKD (ACR >= 30) vs HC
    pb90_celltype <- subset(pb90_subset_clean, 
                            KPMP_celltype_general %in% cell_types & 
                              dkd_group_30_hc != "DKD" &
                              !is.na(pb90_subset_clean@meta.data$dkd_group_30_hc))
    
    pb90_celltype$dkd_group_30_hc <- factor(pb90_celltype$dkd_group_30_hc)
    pb90_celltype$dkd_group_30_hc <- relevel(pb90_celltype$dkd_group_30_hc, ref = "HC")
    
    cat(sprintf("Sample counts for %s:\n", cell_type_group))
    print(table(pb90_celltype@meta.data %>%
                  distinct(record_id, .keep_all = T) %>%
                  pull(dkd_group_30_hc)))
    
    nebula_res <- run_nebula_parallel(pb90_celltype,
                                      n_cores = 15,
                                      subject_id_col = "record_id",
                                      formula = ~ dkd_group_30_hc,
                                      s3_bucket = "scrna",
                                      s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/nonDKD_30_vs_hc/", 
                                                      cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_nondkd30_hc.rds"),
                                      group = F)
    
    processed <- process_nebula_results(nebula_res$results,
                                        pval_col = "p_dkd_group_30_hcnon_DKD",
                                        logfc_col = "logFC_dkd_group_30_hcnon_DKD")
    
  } else if (analysis_type == "GLP_N_vs_HC") {
    # GLP- vs HC
    pb90_celltype <- subset(pb90_subset_clean, 
                            KPMP_celltype_general %in% cell_types & 
                              glp_t2dob != "GLP_Y" & 
                              !is.na(pb90_subset_clean@meta.data$glp_t2dob))
    
    pb90_celltype$glp_t2dob <- factor(pb90_celltype$glp_t2dob)
    pb90_celltype$glp_t2dob <- relevel(pb90_celltype$glp_t2dob, ref = "HC")
    
    cat(sprintf("Sample counts for %s:\n", cell_type_group))
    print(table(pb90_celltype@meta.data %>%
                  distinct(record_id, .keep_all = T) %>%
                  pull(glp_t2dob)))
    
    nebula_res <- run_nebula_parallel(pb90_celltype,
                                      n_cores = 15,
                                      subject_id_col = "record_id",
                                      formula = ~ glp_t2dob,
                                      s3_bucket = "scrna",
                                      s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/GLP_N_vs_HC/", 
                                                      cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_glpn_hc.rds"),
                                      group = F)
    
    processed <- process_nebula_results(nebula_res$results,
                                        pval_col = "p_glp_t2dobGLP_N",
                                        logfc_col = "logFC_glp_t2dobGLP_N")
    
  } else if (analysis_type == "GLP_Y_vs_GLP_N") {
    # GLP+ vs GLP-
    pb90_celltype <- subset(pb90_subset_clean, 
                            KPMP_celltype_general %in% cell_types & 
                              group != "Lean_Control" & 
                              !is.na(pb90_subset_clean@meta.data$glp_t2dob))
    
    pb90_celltype$glp_t2dob <- factor(pb90_celltype$glp_t2dob)
    pb90_celltype$glp_t2dob <- relevel(pb90_celltype$glp_t2dob, ref = "GLP_N")
    
    cat(sprintf("Sample counts for %s:\n", cell_type_group))
    print(table(pb90_celltype@meta.data %>%
                  distinct(record_id, .keep_all = T) %>%
                  pull(glp_t2dob)))
    
    nebula_res <- run_nebula_parallel(pb90_celltype,
                                      n_cores = 15,
                                      subject_id_col = "record_id",
                                      formula = ~ glp_t2dob,
                                      s3_bucket = "scrna",
                                      s3_key = paste0("Projects/CKD/RH_RH2/Results/nebula/GLP_Y_vs_GLP_N/", 
                                                      cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_glpy_glpn.rds"),
                                      group = F)
    
    processed <- process_nebula_results(nebula_res$results,
                                        pval_col = "p_glp_t2dobGLP_Y",
                                        logfc_col = "logFC_glp_t2dobGLP_Y")
  } else {
    stop(sprintf("Unknown analysis type: %s", analysis_type))
  }
  
  return(processed)
}

# Run the analysis
cat(sprintf("Running %s analysis for %s...\n", analysis_type, cell_type_group))
processed <- run_analysis(analysis_type, cell_type_group)

# Annotate with gene information
cat("Adding gene annotations...\n")
gene_symbols <- processed$results$Gene

gene_info <- getBM(
  attributes = c("hgnc_symbol", "description", "gene_biotype"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = mart
) %>%
  dplyr::rename(Gene = hgnc_symbol)

annotated_df <- processed$results %>%
  left_join(gene_info, by = "Gene")

# Save processed results
output_key <- switch(analysis_type,
                     "DKD_vs_nonDKD_100" = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_vs_nonDKD_100/", 
                                                  cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_dkd100_processed.rds"),
                     "DKD_100_vs_HC" = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_100_vs_hc/", 
                                              cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_dkd100_hc_processed.rds"),
                     "nonDKD_100_vs_HC" = paste0("Projects/CKD/RH_RH2/Results/nebula/nonDKD_100_vs_hc/", 
                                                 cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_nondkd100_hc_processed.rds"),
                     "DKD_vs_nonDKD_30" = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_vs_nonDKD_30/", 
                                                 cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_dkd30_processed.rds"),
                     "DKD_30_vs_HC" = paste0("Projects/CKD/RH_RH2/Results/nebula/DKD_30_vs_hc/", 
                                             cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_dkd30_hc_processed.rds"),
                     "nonDKD_30_vs_HC" = paste0("Projects/CKD/RH_RH2/Results/nebula/nonDKD_30_vs_hc/", 
                                                cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_nondkd30_hc_processed.rds"),
                     "GLP_N_vs_HC" = paste0("Projects/CKD/RH_RH2/Results/nebula/GLP_N_vs_HC/", 
                                            cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_glpn_hc_processed.rds"),
                     "GLP_Y_vs_GLP_N" = paste0("Projects/CKD/RH_RH2/Results/nebula/GLP_Y_vs_GLP_N/", 
                                               cell_type_group, "/", cell_type_group, "_rh_rh2_imp_nebula_kpmp_glpy_glpn_processed.rds")
)

cat(sprintf("Saving results to S3: %s\n", output_key))
s3saveRDS(annotated_df, object = output_key, bucket = "scrna", region = "")

# Also save the variable locally with a lowercase name (as in original code)
var_name <- paste0(tolower(cell_type_group), "_kpmp_", 
                   gsub("_", "", tolower(analysis_type)))
assign(var_name, annotated_df, envir = .GlobalEnv)

cat(sprintf("Analysis complete for %s - %s at %s\n", analysis_type, cell_type_group, Sys.time()))