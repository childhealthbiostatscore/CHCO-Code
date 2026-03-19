#!/usr/bin/env Rscript
# =============================================================================
# T1D Adiposity - NEBULA Analysis (Single Job Runner)
# =============================================================================
# Usage: Rscript run_nebula_single.R <analysis_type> <celltype> [celltype_var]
#
# Handles categorical group comparisons, continuous variable associations,
# and interaction models (continuous * group).
# Reads cell type-specific subsets from S3 for efficiency.
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript run_nebula_single.R <analysis_type> <celltype> [celltype_var]")
}

analysis_type <- args[1]
celltype_group_input <- args[2]
celltype_var <- if (length(args) == 3) args[3] else "KPMP_celltype_general"

cat(sprintf("=============================================================\n"))
cat(sprintf("Starting analysis: %s\n", analysis_type))
cat(sprintf("Cell type: %s (using %s)\n", celltype_group_input, celltype_var))
cat(sprintf("Time: %s\n", Sys.time()))
cat(sprintf("=============================================================\n"))

# Load libraries
library(aws.s3)
library(jsonlite)
library(biomaRt)
library(Seurat)
library(dplyr)
library(nebula)

# Set up paths
user <- "yejichoi"
root_path <- "/mmfs1/gscratch/togo/yejichoi/"
git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"

# Set up AWS environment
keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# S3 bucket and base path for this project
s3_bucket <- "t1d.adiposity"
s3_base <- "results/nebula"

# =============================================================================
# ANALYSIS CONFIGURATION
# =============================================================================
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

map_celltype_to_general <- function(celltype, celltype_groups) {
  for (group_name in names(celltype_groups)) {
    if (celltype %in% celltype_groups[[group_name]]) {
      return(group_name)
    }
  }
  return("Other")
}

# =============================================================================
# DEFINE ALL ANALYSES
# =============================================================================
# analysis_mode: "categorical" for group comparisons, "continuous" for variable
#   associations, "interaction" for continuous * group interaction models
# For continuous: continuous_var = the variable, adjust_group = whether to adjust for disease group
# For categorical/continuous: adjust_covariates = optional vector of covariates to include in model
# For interaction: continuous_var * group (tests whether association differs by disease group)

analysis_config <- list(
  
  # =========================================================================
  # CATEGORICAL COMPARISONS (BMI-defined)
  # =========================================================================
  
  # 1. T1D Normal vs T1D Overweight (BMI)
  T1D_normal_vs_overweight_bmi = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity %in% c('Normal', 'Overweight')",
    group_var = "bmi_obesity", ref_level = "Normal",
    pval_col = "p_bmi_obesityOverweight", logfc_col = "logFC_bmi_obesityOverweight",
    s3_subdir = "T1D_normal_vs_overweight_bmi", file_suffix = "t1d_norm_ow_bmi"
  ),
  
  # 2. T1D Non-obese vs T1D Obese (BMI)
  T1D_nonobese_vs_obese_bmi = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity)",
    group_var = "bmi_obese_binary", ref_level = "Non_Obese",
    pval_col = "p_bmi_obese_binaryObese", logfc_col = "logFC_bmi_obese_binaryObese",
    s3_subdir = "T1D_nonobese_vs_obese_bmi", file_suffix = "t1d_nonob_ob_bmi",
    needs_binary_bmi = TRUE
  ),
  
  # 3. T1D Normal vs T1D Obese excluding Overweight (BMI)
  T1D_normal_vs_obese_bmi = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity %in% c('Normal', 'Obese')",
    group_var = "bmi_obesity", ref_level = "Normal",
    pval_col = "p_bmi_obesityObese", logfc_col = "logFC_bmi_obesityObese",
    s3_subdir = "T1D_normal_vs_obese_bmi", file_suffix = "t1d_norm_ob_bmi"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS (DXA-defined)
  # =========================================================================
  
  # 4. T1D Normal vs T1D Overweight (DXA)
  T1D_normal_vs_overweight_dxa = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity %in% c('Normal', 'Overweight')",
    group_var = "dxa_obesity", ref_level = "Normal",
    pval_col = "p_dxa_obesityOverweight", logfc_col = "logFC_dxa_obesityOverweight",
    s3_subdir = "T1D_normal_vs_overweight_dxa", file_suffix = "t1d_norm_ow_dxa"
  ),
  
  # 5. T1D Non-obese vs T1D Obese (DXA)
  T1D_nonobese_vs_obese_dxa = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity)",
    group_var = "dxa_obese_binary", ref_level = "Non_Obese",
    pval_col = "p_dxa_obese_binaryObese", logfc_col = "logFC_dxa_obese_binaryObese",
    s3_subdir = "T1D_nonobese_vs_obese_dxa", file_suffix = "t1d_nonob_ob_dxa",
    needs_binary_dxa = TRUE
  ),
  
  # 6. T1D Normal vs T1D Obese excluding Overweight (DXA)
  T1D_normal_vs_obese_dxa = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity %in% c('Normal', 'Obese')",
    group_var = "dxa_obesity", ref_level = "Normal",
    pval_col = "p_dxa_obesityObese", logfc_col = "logFC_dxa_obesityObese",
    s3_subdir = "T1D_normal_vs_obese_dxa", file_suffix = "t1d_norm_ob_dxa"
  ),
  
  # =========================================================================
  # HC vs T1D COMPARISONS (BMI-defined)
  # =========================================================================
  
  # 7. HC vs T1D Normal weight (BMI)
  HC_vs_T1D_normal_bmi = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Normal'))",
    group_var = "group", ref_level = "Lean Control",
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_normal_bmi", file_suffix = "hc_t1d_norm_bmi"
  ),
  
  # 8. HC vs T1D Overweight (BMI)
  HC_vs_T1D_overweight_bmi = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Overweight'))",
    group_var = "group", ref_level = "Lean Control",
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_overweight_bmi", file_suffix = "hc_t1d_ow_bmi"
  ),
  
  # 9. HC vs T1D Obese (BMI)
  HC_vs_T1D_obese_bmi = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Obese'))",
    group_var = "group", ref_level = "Lean Control",
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_obese_bmi", file_suffix = "hc_t1d_ob_bmi"
  ),
  
  # =========================================================================
  # HC vs T1D COMPARISONS (DXA-defined)
  # =========================================================================
  
  # 7d. HC vs T1D Normal weight (DXA)
  HC_vs_T1D_normal_dxa = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Normal'))",
    group_var = "group", ref_level = "Lean Control",
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_normal_dxa", file_suffix = "hc_t1d_norm_dxa"
  ),
  
  # 8d. HC vs T1D Overweight (DXA)
  HC_vs_T1D_overweight_dxa = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Overweight'))",
    group_var = "group", ref_level = "Lean Control",
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_overweight_dxa", file_suffix = "hc_t1d_ow_dxa"
  ),
  
  # 9d. HC vs T1D Obese (DXA)
  HC_vs_T1D_obese_dxa = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Obese'))",
    group_var = "group", ref_level = "Lean Control",
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_obese_dxa", file_suffix = "hc_t1d_ob_dxa"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS - AGE-ADJUSTED (BMI-defined, T1D within-group)
  # =========================================================================
  
  T1D_normal_vs_overweight_bmi_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity %in% c('Normal', 'Overweight')",
    group_var = "bmi_obesity", ref_level = "Normal",
    adjust_covariates = c("age"),
    pval_col = "p_bmi_obesityOverweight", logfc_col = "logFC_bmi_obesityOverweight",
    s3_subdir = "T1D_normal_vs_overweight_bmi_adj_age", file_suffix = "t1d_norm_ow_bmi_adj_age"
  ),
  
  T1D_nonobese_vs_obese_bmi_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity)",
    group_var = "bmi_obese_binary", ref_level = "Non_Obese",
    adjust_covariates = c("age"),
    pval_col = "p_bmi_obese_binaryObese", logfc_col = "logFC_bmi_obese_binaryObese",
    s3_subdir = "T1D_nonobese_vs_obese_bmi_adj_age", file_suffix = "t1d_nonob_ob_bmi_adj_age",
    needs_binary_bmi = TRUE
  ),
  
  T1D_normal_vs_obese_bmi_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity %in% c('Normal', 'Obese')",
    group_var = "bmi_obesity", ref_level = "Normal",
    adjust_covariates = c("age"),
    pval_col = "p_bmi_obesityObese", logfc_col = "logFC_bmi_obesityObese",
    s3_subdir = "T1D_normal_vs_obese_bmi_adj_age", file_suffix = "t1d_norm_ob_bmi_adj_age"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS - AGE-ADJUSTED (DXA-defined, T1D within-group)
  # =========================================================================
  
  T1D_normal_vs_overweight_dxa_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity %in% c('Normal', 'Overweight')",
    group_var = "dxa_obesity", ref_level = "Normal",
    adjust_covariates = c("age"),
    pval_col = "p_dxa_obesityOverweight", logfc_col = "logFC_dxa_obesityOverweight",
    s3_subdir = "T1D_normal_vs_overweight_dxa_adj_age", file_suffix = "t1d_norm_ow_dxa_adj_age"
  ),
  
  T1D_nonobese_vs_obese_dxa_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity)",
    group_var = "dxa_obese_binary", ref_level = "Non_Obese",
    adjust_covariates = c("age"),
    pval_col = "p_dxa_obese_binaryObese", logfc_col = "logFC_dxa_obese_binaryObese",
    s3_subdir = "T1D_nonobese_vs_obese_dxa_adj_age", file_suffix = "t1d_nonob_ob_dxa_adj_age",
    needs_binary_dxa = TRUE
  ),
  
  T1D_normal_vs_obese_dxa_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity %in% c('Normal', 'Obese')",
    group_var = "dxa_obesity", ref_level = "Normal",
    adjust_covariates = c("age"),
    pval_col = "p_dxa_obesityObese", logfc_col = "logFC_dxa_obesityObese",
    s3_subdir = "T1D_normal_vs_obese_dxa_adj_age", file_suffix = "t1d_norm_ob_dxa_adj_age"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS - AGE-ADJUSTED (HC vs T1D, BMI-defined)
  # =========================================================================
  
  HC_vs_T1D_normal_bmi_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Normal'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_normal_bmi_adj_age", file_suffix = "hc_t1d_norm_bmi_adj_age"
  ),
  
  HC_vs_T1D_overweight_bmi_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Overweight'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_overweight_bmi_adj_age", file_suffix = "hc_t1d_ow_bmi_adj_age"
  ),
  
  HC_vs_T1D_obese_bmi_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Obese'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_obese_bmi_adj_age", file_suffix = "hc_t1d_ob_bmi_adj_age"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS - AGE-ADJUSTED (HC vs T1D, DXA-defined)
  # =========================================================================
  
  HC_vs_T1D_normal_dxa_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Normal'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_normal_dxa_adj_age", file_suffix = "hc_t1d_norm_dxa_adj_age"
  ),
  
  HC_vs_T1D_overweight_dxa_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Overweight'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_overweight_dxa_adj_age", file_suffix = "hc_t1d_ow_dxa_adj_age"
  ),
  
  HC_vs_T1D_obese_dxa_adj_age = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Obese'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_obese_dxa_adj_age", file_suffix = "hc_t1d_ob_dxa_adj_age"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS - AGE+SEX ADJUSTED (BMI-defined, T1D within-group)
  # =========================================================================
  
  T1D_normal_vs_overweight_bmi_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity %in% c('Normal', 'Overweight')",
    group_var = "bmi_obesity", ref_level = "Normal",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_bmi_obesityOverweight", logfc_col = "logFC_bmi_obesityOverweight",
    s3_subdir = "T1D_normal_vs_overweight_bmi_adj_age_sex", file_suffix = "t1d_norm_ow_bmi_adj_age_sex"
  ),
  
  T1D_nonobese_vs_obese_bmi_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity)",
    group_var = "bmi_obese_binary", ref_level = "Non_Obese",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_bmi_obese_binaryObese", logfc_col = "logFC_bmi_obese_binaryObese",
    s3_subdir = "T1D_nonobese_vs_obese_bmi_adj_age_sex", file_suffix = "t1d_nonob_ob_bmi_adj_age_sex",
    needs_binary_bmi = TRUE
  ),
  
  T1D_normal_vs_obese_bmi_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity %in% c('Normal', 'Obese')",
    group_var = "bmi_obesity", ref_level = "Normal",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_bmi_obesityObese", logfc_col = "logFC_bmi_obesityObese",
    s3_subdir = "T1D_normal_vs_obese_bmi_adj_age_sex", file_suffix = "t1d_norm_ob_bmi_adj_age_sex"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS - AGE+SEX ADJUSTED (DXA-defined, T1D within-group)
  # =========================================================================
  
  T1D_normal_vs_overweight_dxa_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity %in% c('Normal', 'Overweight')",
    group_var = "dxa_obesity", ref_level = "Normal",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_dxa_obesityOverweight", logfc_col = "logFC_dxa_obesityOverweight",
    s3_subdir = "T1D_normal_vs_overweight_dxa_adj_age_sex", file_suffix = "t1d_norm_ow_dxa_adj_age_sex"
  ),
  
  T1D_nonobese_vs_obese_dxa_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity)",
    group_var = "dxa_obese_binary", ref_level = "Non_Obese",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_dxa_obese_binaryObese", logfc_col = "logFC_dxa_obese_binaryObese",
    s3_subdir = "T1D_nonobese_vs_obese_dxa_adj_age_sex", file_suffix = "t1d_nonob_ob_dxa_adj_age_sex",
    needs_binary_dxa = TRUE
  ),
  
  T1D_normal_vs_obese_dxa_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity %in% c('Normal', 'Obese')",
    group_var = "dxa_obesity", ref_level = "Normal",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_dxa_obesityObese", logfc_col = "logFC_dxa_obesityObese",
    s3_subdir = "T1D_normal_vs_obese_dxa_adj_age_sex", file_suffix = "t1d_norm_ob_dxa_adj_age_sex"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS - AGE+SEX ADJUSTED (HC vs T1D, BMI-defined)
  # =========================================================================
  
  HC_vs_T1D_normal_bmi_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Normal'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_normal_bmi_adj_age_sex", file_suffix = "hc_t1d_norm_bmi_adj_age_sex"
  ),
  
  HC_vs_T1D_overweight_bmi_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Overweight'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_overweight_bmi_adj_age_sex", file_suffix = "hc_t1d_ow_bmi_adj_age_sex"
  ),
  
  HC_vs_T1D_obese_bmi_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(bmi_obesity) & bmi_obesity == 'Obese'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_obese_bmi_adj_age_sex", file_suffix = "hc_t1d_ob_bmi_adj_age_sex"
  ),
  
  # =========================================================================
  # CATEGORICAL COMPARISONS - AGE+SEX ADJUSTED (HC vs T1D, DXA-defined)
  # =========================================================================
  
  HC_vs_T1D_normal_dxa_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Normal'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_normal_dxa_adj_age_sex", file_suffix = "hc_t1d_norm_dxa_adj_age_sex"
  ),
  
  HC_vs_T1D_overweight_dxa_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Overweight'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_overweight_dxa_adj_age_sex", file_suffix = "hc_t1d_ow_dxa_adj_age_sex"
  ),
  
  HC_vs_T1D_obese_dxa_adj_age_sex = list(
    analysis_mode = "categorical",
    subset_cond = "(group == 'Lean Control' | (group == 'Type 1 Diabetes' & !is.na(dxa_obesity) & dxa_obesity == 'Obese'))",
    group_var = "group", ref_level = "Lean Control",
    adjust_covariates = c("age", "sex"),
    pval_col = "p_groupType 1 Diabetes", logfc_col = "logFC_groupType 1 Diabetes",
    s3_subdir = "HC_vs_T1D_obese_dxa_adj_age_sex", file_suffix = "hc_t1d_ob_dxa_adj_age_sex"
  ),
  
  # =========================================================================
  # CONTINUOUS ASSOCIATIONS - T1D ONLY (unadjusted)
  # =========================================================================
  
  cont_bmi_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi)",
    continuous_var = "bmi", adjust_group = FALSE,
    pval_col = "p_bmi", logfc_col = "logFC_bmi",
    s3_subdir = "continuous/bmi_t1d", file_suffix = "cont_bmi_t1d"
  ),
  cont_dexa_body_fat_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_body_fat)",
    continuous_var = "dexa_body_fat", adjust_group = FALSE,
    pval_col = "p_dexa_body_fat", logfc_col = "logFC_dexa_body_fat",
    s3_subdir = "continuous/dexa_body_fat_t1d", file_suffix = "cont_dexa_body_fat_t1d"
  ),
  cont_dexa_bone_mineral_density_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_bone_mineral_density)",
    continuous_var = "dexa_bone_mineral_density", adjust_group = FALSE,
    pval_col = "p_dexa_bone_mineral_density", logfc_col = "logFC_dexa_bone_mineral_density",
    s3_subdir = "continuous/dexa_bone_mineral_density_t1d", file_suffix = "cont_dexa_bmd_t1d"
  ),
  cont_dexa_fat_kg_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_fat_kg)",
    continuous_var = "dexa_fat_kg", adjust_group = FALSE,
    pval_col = "p_dexa_fat_kg", logfc_col = "logFC_dexa_fat_kg",
    s3_subdir = "continuous/dexa_fat_kg_t1d", file_suffix = "cont_dexa_fat_kg_t1d"
  ),
  cont_dexa_lean_mass_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_lean_mass)",
    continuous_var = "dexa_lean_mass", adjust_group = FALSE,
    pval_col = "p_dexa_lean_mass", logfc_col = "logFC_dexa_lean_mass",
    s3_subdir = "continuous/dexa_lean_mass_t1d", file_suffix = "cont_dexa_lean_mass_t1d"
  ),
  cont_dexa_lean_kg_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_lean_kg)",
    continuous_var = "dexa_lean_kg", adjust_group = FALSE,
    pval_col = "p_dexa_lean_kg", logfc_col = "logFC_dexa_lean_kg",
    s3_subdir = "continuous/dexa_lean_kg_t1d", file_suffix = "cont_dexa_lean_kg_t1d"
  ),
  cont_dexa_ag_ratio_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_ag_ratio)",
    continuous_var = "dexa_ag_ratio", adjust_group = FALSE,
    pval_col = "p_dexa_ag_ratio", logfc_col = "logFC_dexa_ag_ratio",
    s3_subdir = "continuous/dexa_ag_ratio_t1d", file_suffix = "cont_dexa_ag_ratio_t1d"
  ),
  cont_dexa_est_vat_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_est_vat)",
    continuous_var = "dexa_est_vat", adjust_group = FALSE,
    pval_col = "p_dexa_est_vat", logfc_col = "logFC_dexa_est_vat",
    s3_subdir = "continuous/dexa_est_vat_t1d", file_suffix = "cont_dexa_est_vat_t1d"
  ),
  cont_dexa_trunk_kg_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_trunk_kg)",
    continuous_var = "dexa_trunk_kg", adjust_group = FALSE,
    pval_col = "p_dexa_trunk_kg", logfc_col = "logFC_dexa_trunk_kg",
    s3_subdir = "continuous/dexa_trunk_kg_t1d", file_suffix = "cont_dexa_trunk_kg_t1d"
  ),
  cont_dexa_trunk_mass_t1d = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_trunk_mass)",
    continuous_var = "dexa_trunk_mass", adjust_group = FALSE,
    pval_col = "p_dexa_trunk_mass", logfc_col = "logFC_dexa_trunk_mass",
    s3_subdir = "continuous/dexa_trunk_mass_t1d", file_suffix = "cont_dexa_trunk_mass_t1d"
  ),
  
  # =========================================================================
  # CONTINUOUS ASSOCIATIONS - ALL SUBJECTS, UNADJUSTED
  # =========================================================================
  
  cont_bmi_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(bmi)",
    continuous_var = "bmi", adjust_group = FALSE,
    pval_col = "p_bmi", logfc_col = "logFC_bmi",
    s3_subdir = "continuous/bmi_all", file_suffix = "cont_bmi_all"
  ),
  cont_dexa_body_fat_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_body_fat)",
    continuous_var = "dexa_body_fat", adjust_group = FALSE,
    pval_col = "p_dexa_body_fat", logfc_col = "logFC_dexa_body_fat",
    s3_subdir = "continuous/dexa_body_fat_all", file_suffix = "cont_dexa_body_fat_all"
  ),
  cont_dexa_bone_mineral_density_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_bone_mineral_density)",
    continuous_var = "dexa_bone_mineral_density", adjust_group = FALSE,
    pval_col = "p_dexa_bone_mineral_density", logfc_col = "logFC_dexa_bone_mineral_density",
    s3_subdir = "continuous/dexa_bone_mineral_density_all", file_suffix = "cont_dexa_bmd_all"
  ),
  cont_dexa_fat_kg_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_fat_kg)",
    continuous_var = "dexa_fat_kg", adjust_group = FALSE,
    pval_col = "p_dexa_fat_kg", logfc_col = "logFC_dexa_fat_kg",
    s3_subdir = "continuous/dexa_fat_kg_all", file_suffix = "cont_dexa_fat_kg_all"
  ),
  cont_dexa_lean_mass_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_lean_mass)",
    continuous_var = "dexa_lean_mass", adjust_group = FALSE,
    pval_col = "p_dexa_lean_mass", logfc_col = "logFC_dexa_lean_mass",
    s3_subdir = "continuous/dexa_lean_mass_all", file_suffix = "cont_dexa_lean_mass_all"
  ),
  cont_dexa_lean_kg_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_lean_kg)",
    continuous_var = "dexa_lean_kg", adjust_group = FALSE,
    pval_col = "p_dexa_lean_kg", logfc_col = "logFC_dexa_lean_kg",
    s3_subdir = "continuous/dexa_lean_kg_all", file_suffix = "cont_dexa_lean_kg_all"
  ),
  cont_dexa_ag_ratio_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_ag_ratio)",
    continuous_var = "dexa_ag_ratio", adjust_group = FALSE,
    pval_col = "p_dexa_ag_ratio", logfc_col = "logFC_dexa_ag_ratio",
    s3_subdir = "continuous/dexa_ag_ratio_all", file_suffix = "cont_dexa_ag_ratio_all"
  ),
  cont_dexa_est_vat_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_est_vat)",
    continuous_var = "dexa_est_vat", adjust_group = FALSE,
    pval_col = "p_dexa_est_vat", logfc_col = "logFC_dexa_est_vat",
    s3_subdir = "continuous/dexa_est_vat_all", file_suffix = "cont_dexa_est_vat_all"
  ),
  cont_dexa_trunk_kg_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_trunk_kg)",
    continuous_var = "dexa_trunk_kg", adjust_group = FALSE,
    pval_col = "p_dexa_trunk_kg", logfc_col = "logFC_dexa_trunk_kg",
    s3_subdir = "continuous/dexa_trunk_kg_all", file_suffix = "cont_dexa_trunk_kg_all"
  ),
  cont_dexa_trunk_mass_all = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_trunk_mass)",
    continuous_var = "dexa_trunk_mass", adjust_group = FALSE,
    pval_col = "p_dexa_trunk_mass", logfc_col = "logFC_dexa_trunk_mass",
    s3_subdir = "continuous/dexa_trunk_mass_all", file_suffix = "cont_dexa_trunk_mass_all"
  ),
  
  # =========================================================================
  # CONTINUOUS ASSOCIATIONS - ALL SUBJECTS, ADJUSTED FOR GROUP
  # =========================================================================
  
  cont_bmi_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(bmi)",
    continuous_var = "bmi", adjust_group = TRUE,
    pval_col = "p_bmi", logfc_col = "logFC_bmi",
    s3_subdir = "continuous/bmi_all_adj_group", file_suffix = "cont_bmi_all_adj"
  ),
  cont_dexa_body_fat_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_body_fat)",
    continuous_var = "dexa_body_fat", adjust_group = TRUE,
    pval_col = "p_dexa_body_fat", logfc_col = "logFC_dexa_body_fat",
    s3_subdir = "continuous/dexa_body_fat_all_adj_group", file_suffix = "cont_dexa_body_fat_all_adj"
  ),
  cont_dexa_bone_mineral_density_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_bone_mineral_density)",
    continuous_var = "dexa_bone_mineral_density", adjust_group = TRUE,
    pval_col = "p_dexa_bone_mineral_density", logfc_col = "logFC_dexa_bone_mineral_density",
    s3_subdir = "continuous/dexa_bone_mineral_density_all_adj_group", file_suffix = "cont_dexa_bmd_all_adj"
  ),
  cont_dexa_fat_kg_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_fat_kg)",
    continuous_var = "dexa_fat_kg", adjust_group = TRUE,
    pval_col = "p_dexa_fat_kg", logfc_col = "logFC_dexa_fat_kg",
    s3_subdir = "continuous/dexa_fat_kg_all_adj_group", file_suffix = "cont_dexa_fat_kg_all_adj"
  ),
  cont_dexa_lean_mass_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_lean_mass)",
    continuous_var = "dexa_lean_mass", adjust_group = TRUE,
    pval_col = "p_dexa_lean_mass", logfc_col = "logFC_dexa_lean_mass",
    s3_subdir = "continuous/dexa_lean_mass_all_adj_group", file_suffix = "cont_dexa_lean_mass_all_adj"
  ),
  cont_dexa_lean_kg_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_lean_kg)",
    continuous_var = "dexa_lean_kg", adjust_group = TRUE,
    pval_col = "p_dexa_lean_kg", logfc_col = "logFC_dexa_lean_kg",
    s3_subdir = "continuous/dexa_lean_kg_all_adj_group", file_suffix = "cont_dexa_lean_kg_all_adj"
  ),
  cont_dexa_ag_ratio_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_ag_ratio)",
    continuous_var = "dexa_ag_ratio", adjust_group = TRUE,
    pval_col = "p_dexa_ag_ratio", logfc_col = "logFC_dexa_ag_ratio",
    s3_subdir = "continuous/dexa_ag_ratio_all_adj_group", file_suffix = "cont_dexa_ag_ratio_all_adj"
  ),
  cont_dexa_est_vat_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_est_vat)",
    continuous_var = "dexa_est_vat", adjust_group = TRUE,
    pval_col = "p_dexa_est_vat", logfc_col = "logFC_dexa_est_vat",
    s3_subdir = "continuous/dexa_est_vat_all_adj_group", file_suffix = "cont_dexa_est_vat_all_adj"
  ),
  cont_dexa_trunk_kg_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_trunk_kg)",
    continuous_var = "dexa_trunk_kg", adjust_group = TRUE,
    pval_col = "p_dexa_trunk_kg", logfc_col = "logFC_dexa_trunk_kg",
    s3_subdir = "continuous/dexa_trunk_kg_all_adj_group", file_suffix = "cont_dexa_trunk_kg_all_adj"
  ),
  cont_dexa_trunk_mass_all_adj = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_trunk_mass)",
    continuous_var = "dexa_trunk_mass", adjust_group = TRUE,
    pval_col = "p_dexa_trunk_mass", logfc_col = "logFC_dexa_trunk_mass",
    s3_subdir = "continuous/dexa_trunk_mass_all_adj_group", file_suffix = "cont_dexa_trunk_mass_all_adj"
  ),
  
  # =========================================================================
  # CONTINUOUS ASSOCIATIONS - T1D ONLY, ADJUSTED FOR AGE
  # =========================================================================
  
  cont_bmi_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(bmi)",
    continuous_var = "bmi", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_bmi", logfc_col = "logFC_bmi",
    s3_subdir = "continuous/bmi_t1d_adj_age", file_suffix = "cont_bmi_t1d_adj_age"
  ),
  cont_dexa_body_fat_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_body_fat)",
    continuous_var = "dexa_body_fat", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_body_fat", logfc_col = "logFC_dexa_body_fat",
    s3_subdir = "continuous/dexa_body_fat_t1d_adj_age", file_suffix = "cont_dexa_body_fat_t1d_adj_age"
  ),
  cont_dexa_bone_mineral_density_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_bone_mineral_density)",
    continuous_var = "dexa_bone_mineral_density", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_bone_mineral_density", logfc_col = "logFC_dexa_bone_mineral_density",
    s3_subdir = "continuous/dexa_bone_mineral_density_t1d_adj_age", file_suffix = "cont_dexa_bmd_t1d_adj_age"
  ),
  cont_dexa_fat_kg_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_fat_kg)",
    continuous_var = "dexa_fat_kg", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_fat_kg", logfc_col = "logFC_dexa_fat_kg",
    s3_subdir = "continuous/dexa_fat_kg_t1d_adj_age", file_suffix = "cont_dexa_fat_kg_t1d_adj_age"
  ),
  cont_dexa_lean_mass_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_lean_mass)",
    continuous_var = "dexa_lean_mass", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_lean_mass", logfc_col = "logFC_dexa_lean_mass",
    s3_subdir = "continuous/dexa_lean_mass_t1d_adj_age", file_suffix = "cont_dexa_lean_mass_t1d_adj_age"
  ),
  cont_dexa_lean_kg_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_lean_kg)",
    continuous_var = "dexa_lean_kg", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_lean_kg", logfc_col = "logFC_dexa_lean_kg",
    s3_subdir = "continuous/dexa_lean_kg_t1d_adj_age", file_suffix = "cont_dexa_lean_kg_t1d_adj_age"
  ),
  cont_dexa_ag_ratio_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_ag_ratio)",
    continuous_var = "dexa_ag_ratio", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_ag_ratio", logfc_col = "logFC_dexa_ag_ratio",
    s3_subdir = "continuous/dexa_ag_ratio_t1d_adj_age", file_suffix = "cont_dexa_ag_ratio_t1d_adj_age"
  ),
  cont_dexa_est_vat_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_est_vat)",
    continuous_var = "dexa_est_vat", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_est_vat", logfc_col = "logFC_dexa_est_vat",
    s3_subdir = "continuous/dexa_est_vat_t1d_adj_age", file_suffix = "cont_dexa_est_vat_t1d_adj_age"
  ),
  cont_dexa_trunk_kg_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_trunk_kg)",
    continuous_var = "dexa_trunk_kg", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_trunk_kg", logfc_col = "logFC_dexa_trunk_kg",
    s3_subdir = "continuous/dexa_trunk_kg_t1d_adj_age", file_suffix = "cont_dexa_trunk_kg_t1d_adj_age"
  ),
  cont_dexa_trunk_mass_t1d_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "group == 'Type 1 Diabetes' & !is.na(dexa_trunk_mass)",
    continuous_var = "dexa_trunk_mass", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_trunk_mass", logfc_col = "logFC_dexa_trunk_mass",
    s3_subdir = "continuous/dexa_trunk_mass_t1d_adj_age", file_suffix = "cont_dexa_trunk_mass_t1d_adj_age"
  ),
  
  # =========================================================================
  # CONTINUOUS ASSOCIATIONS - ALL SUBJECTS, ADJUSTED FOR AGE (no group)
  # =========================================================================
  
  cont_bmi_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(bmi)",
    continuous_var = "bmi", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_bmi", logfc_col = "logFC_bmi",
    s3_subdir = "continuous/bmi_all_adj_age", file_suffix = "cont_bmi_all_adj_age"
  ),
  cont_dexa_body_fat_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_body_fat)",
    continuous_var = "dexa_body_fat", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_body_fat", logfc_col = "logFC_dexa_body_fat",
    s3_subdir = "continuous/dexa_body_fat_all_adj_age", file_suffix = "cont_dexa_body_fat_all_adj_age"
  ),
  cont_dexa_bone_mineral_density_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_bone_mineral_density)",
    continuous_var = "dexa_bone_mineral_density", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_bone_mineral_density", logfc_col = "logFC_dexa_bone_mineral_density",
    s3_subdir = "continuous/dexa_bone_mineral_density_all_adj_age", file_suffix = "cont_dexa_bmd_all_adj_age"
  ),
  cont_dexa_fat_kg_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_fat_kg)",
    continuous_var = "dexa_fat_kg", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_fat_kg", logfc_col = "logFC_dexa_fat_kg",
    s3_subdir = "continuous/dexa_fat_kg_all_adj_age", file_suffix = "cont_dexa_fat_kg_all_adj_age"
  ),
  cont_dexa_lean_mass_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_lean_mass)",
    continuous_var = "dexa_lean_mass", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_lean_mass", logfc_col = "logFC_dexa_lean_mass",
    s3_subdir = "continuous/dexa_lean_mass_all_adj_age", file_suffix = "cont_dexa_lean_mass_all_adj_age"
  ),
  cont_dexa_lean_kg_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_lean_kg)",
    continuous_var = "dexa_lean_kg", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_lean_kg", logfc_col = "logFC_dexa_lean_kg",
    s3_subdir = "continuous/dexa_lean_kg_all_adj_age", file_suffix = "cont_dexa_lean_kg_all_adj_age"
  ),
  cont_dexa_ag_ratio_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_ag_ratio)",
    continuous_var = "dexa_ag_ratio", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_ag_ratio", logfc_col = "logFC_dexa_ag_ratio",
    s3_subdir = "continuous/dexa_ag_ratio_all_adj_age", file_suffix = "cont_dexa_ag_ratio_all_adj_age"
  ),
  cont_dexa_est_vat_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_est_vat)",
    continuous_var = "dexa_est_vat", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_est_vat", logfc_col = "logFC_dexa_est_vat",
    s3_subdir = "continuous/dexa_est_vat_all_adj_age", file_suffix = "cont_dexa_est_vat_all_adj_age"
  ),
  cont_dexa_trunk_kg_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_trunk_kg)",
    continuous_var = "dexa_trunk_kg", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_trunk_kg", logfc_col = "logFC_dexa_trunk_kg",
    s3_subdir = "continuous/dexa_trunk_kg_all_adj_age", file_suffix = "cont_dexa_trunk_kg_all_adj_age"
  ),
  cont_dexa_trunk_mass_all_adj_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_trunk_mass)",
    continuous_var = "dexa_trunk_mass", adjust_group = FALSE, adjust_covariates = c("age"),
    pval_col = "p_dexa_trunk_mass", logfc_col = "logFC_dexa_trunk_mass",
    s3_subdir = "continuous/dexa_trunk_mass_all_adj_age", file_suffix = "cont_dexa_trunk_mass_all_adj_age"
  ),
  
  # =========================================================================
  # CONTINUOUS ASSOCIATIONS - ALL SUBJECTS, ADJUSTED FOR GROUP + AGE
  # =========================================================================
  
  cont_bmi_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(bmi)",
    continuous_var = "bmi", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_bmi", logfc_col = "logFC_bmi",
    s3_subdir = "continuous/bmi_all_adj_group_age", file_suffix = "cont_bmi_all_adj_group_age"
  ),
  cont_dexa_body_fat_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_body_fat)",
    continuous_var = "dexa_body_fat", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_body_fat", logfc_col = "logFC_dexa_body_fat",
    s3_subdir = "continuous/dexa_body_fat_all_adj_group_age", file_suffix = "cont_dexa_body_fat_all_adj_group_age"
  ),
  cont_dexa_bone_mineral_density_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_bone_mineral_density)",
    continuous_var = "dexa_bone_mineral_density", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_bone_mineral_density", logfc_col = "logFC_dexa_bone_mineral_density",
    s3_subdir = "continuous/dexa_bone_mineral_density_all_adj_group_age", file_suffix = "cont_dexa_bmd_all_adj_group_age"
  ),
  cont_dexa_fat_kg_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_fat_kg)",
    continuous_var = "dexa_fat_kg", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_fat_kg", logfc_col = "logFC_dexa_fat_kg",
    s3_subdir = "continuous/dexa_fat_kg_all_adj_group_age", file_suffix = "cont_dexa_fat_kg_all_adj_group_age"
  ),
  cont_dexa_lean_mass_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_lean_mass)",
    continuous_var = "dexa_lean_mass", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_lean_mass", logfc_col = "logFC_dexa_lean_mass",
    s3_subdir = "continuous/dexa_lean_mass_all_adj_group_age", file_suffix = "cont_dexa_lean_mass_all_adj_group_age"
  ),
  cont_dexa_lean_kg_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_lean_kg)",
    continuous_var = "dexa_lean_kg", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_lean_kg", logfc_col = "logFC_dexa_lean_kg",
    s3_subdir = "continuous/dexa_lean_kg_all_adj_group_age", file_suffix = "cont_dexa_lean_kg_all_adj_group_age"
  ),
  cont_dexa_ag_ratio_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_ag_ratio)",
    continuous_var = "dexa_ag_ratio", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_ag_ratio", logfc_col = "logFC_dexa_ag_ratio",
    s3_subdir = "continuous/dexa_ag_ratio_all_adj_group_age", file_suffix = "cont_dexa_ag_ratio_all_adj_group_age"
  ),
  cont_dexa_est_vat_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_est_vat)",
    continuous_var = "dexa_est_vat", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_est_vat", logfc_col = "logFC_dexa_est_vat",
    s3_subdir = "continuous/dexa_est_vat_all_adj_group_age", file_suffix = "cont_dexa_est_vat_all_adj_group_age"
  ),
  cont_dexa_trunk_kg_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_trunk_kg)",
    continuous_var = "dexa_trunk_kg", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_trunk_kg", logfc_col = "logFC_dexa_trunk_kg",
    s3_subdir = "continuous/dexa_trunk_kg_all_adj_group_age", file_suffix = "cont_dexa_trunk_kg_all_adj_group_age"
  ),
  cont_dexa_trunk_mass_all_adj_group_age = list(
    analysis_mode = "continuous",
    subset_cond = "!is.na(dexa_trunk_mass)",
    continuous_var = "dexa_trunk_mass", adjust_group = TRUE, adjust_covariates = c("age"),
    pval_col = "p_dexa_trunk_mass", logfc_col = "logFC_dexa_trunk_mass",
    s3_subdir = "continuous/dexa_trunk_mass_all_adj_group_age", file_suffix = "cont_dexa_trunk_mass_all_adj_group_age"
  ),
  
  # =========================================================================
  # INTERACTION MODELS - continuous_var * group (HC vs T1D)
  # Tests whether the association of each adiposity variable with gene
  # expression differs between HC and T1D
  # =========================================================================
  
  interaction_bmi = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(bmi) & !is.na(group)",
    continuous_var = "bmi",
    pval_col = "p_bmi:groupType 1 Diabetes",
    logfc_col = "logFC_bmi:groupType 1 Diabetes",
    s3_subdir = "interaction/bmi", file_suffix = "int_bmi"
  ),
  
  interaction_dexa_body_fat = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_body_fat) & !is.na(group)",
    continuous_var = "dexa_body_fat",
    pval_col = "p_dexa_body_fat:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_body_fat:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_body_fat", file_suffix = "int_dexa_body_fat"
  ),
  
  interaction_dexa_bone_mineral_density = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_bone_mineral_density) & !is.na(group)",
    continuous_var = "dexa_bone_mineral_density",
    pval_col = "p_dexa_bone_mineral_density:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_bone_mineral_density:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_bone_mineral_density", file_suffix = "int_dexa_bmd"
  ),
  
  interaction_dexa_fat_kg = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_fat_kg) & !is.na(group)",
    continuous_var = "dexa_fat_kg",
    pval_col = "p_dexa_fat_kg:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_fat_kg:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_fat_kg", file_suffix = "int_dexa_fat_kg"
  ),
  
  interaction_dexa_lean_mass = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_lean_mass) & !is.na(group)",
    continuous_var = "dexa_lean_mass",
    pval_col = "p_dexa_lean_mass:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_lean_mass:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_lean_mass", file_suffix = "int_dexa_lean_mass"
  ),
  
  interaction_dexa_lean_kg = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_lean_kg) & !is.na(group)",
    continuous_var = "dexa_lean_kg",
    pval_col = "p_dexa_lean_kg:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_lean_kg:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_lean_kg", file_suffix = "int_dexa_lean_kg"
  ),
  
  interaction_dexa_ag_ratio = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_ag_ratio) & !is.na(group)",
    continuous_var = "dexa_ag_ratio",
    pval_col = "p_dexa_ag_ratio:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_ag_ratio:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_ag_ratio", file_suffix = "int_dexa_ag_ratio"
  ),
  
  interaction_dexa_est_vat = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_est_vat) & !is.na(group)",
    continuous_var = "dexa_est_vat",
    pval_col = "p_dexa_est_vat:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_est_vat:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_est_vat", file_suffix = "int_dexa_est_vat"
  ),
  
  interaction_dexa_trunk_kg = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_trunk_kg) & !is.na(group)",
    continuous_var = "dexa_trunk_kg",
    pval_col = "p_dexa_trunk_kg:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_trunk_kg:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_trunk_kg", file_suffix = "int_dexa_trunk_kg"
  ),
  
  interaction_dexa_trunk_mass = list(
    analysis_mode = "interaction",
    subset_cond = "!is.na(dexa_trunk_mass) & !is.na(group)",
    continuous_var = "dexa_trunk_mass",
    pval_col = "p_dexa_trunk_mass:groupType 1 Diabetes",
    logfc_col = "logFC_dexa_trunk_mass:groupType 1 Diabetes",
    s3_subdir = "interaction/dexa_trunk_mass", file_suffix = "int_dexa_trunk_mass"
  )
)

# Validate analysis type
if (!analysis_type %in% names(analysis_config)) {
  stop(sprintf("Unknown analysis type: %s\nValid types: %s",
               analysis_type, paste(names(analysis_config), collapse = ", ")))
}

config <- analysis_config[[analysis_type]]

# =============================================================================
# LOAD DATA
# =============================================================================
cat(sprintf("\nLoading cell type-specific data for %s...\n", celltype_group_input))

# Try to load pre-saved cell type subset first
subset_file <- paste0("data_clean/subset/t1d_adiposity_subset_",
                      gsub("/", "_", celltype_group_input), ".rds")

pb90_subset_clean <- tryCatch({
  cat(sprintf("Looking for cell type subset: %s\n", subset_file))
  obj <- s3readRDS(object = subset_file, bucket = s3_bucket, region = "")
  cat(sprintf("Loaded cell type-specific subset with %d cells\n", ncol(obj)))
  obj
}, error = function(e) {
  cat(sprintf("Cell type-specific file not found: %s\n", subset_file))
  cat("Falling back to loading full dataset and subsetting...\n")
  obj <- s3readRDS(object = "data_clean/t1d_hc_scrna_w_clinical.rds",
                   bucket = s3_bucket, region = "")
  cat(sprintf("Loaded full dataset with %d cells\n", ncol(obj)))
  
  # Add KPMP_celltype_general if needed
  if (celltype_var == "KPMP_celltype_general" & !"KPMP_celltype_general" %in% colnames(obj@meta.data)) {
    obj$KPMP_celltype_general <- sapply(obj$KPMP_celltype,
                                        map_celltype_to_general,
                                        celltype_groups = celltype_groups)
  }
  
  obj <- subset(obj, !!sym(celltype_var) == celltype_group_input)
  cat(sprintf("Subset to %s: %d cells\n", celltype_group_input, ncol(obj)))
  obj
})

# Ensure KPMP_celltype_general exists
if (!"KPMP_celltype_general" %in% colnames(pb90_subset_clean@meta.data)) {
  pb90_subset_clean$KPMP_celltype_general <- sapply(
    pb90_subset_clean$KPMP_celltype,
    map_celltype_to_general,
    celltype_groups = celltype_groups
  )
}

# Source helper functions
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_scRNA_functions.R"))
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_functions.R"))

# Initialize biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# =============================================================================
# PREPARE DATA FOR ANALYSIS
# =============================================================================
gc()

# Create derived binary obesity variables if needed
if (isTRUE(config$needs_binary_bmi)) {
  pb90_subset_clean@meta.data <- pb90_subset_clean@meta.data %>%
    mutate(bmi_obese_binary = case_when(
      bmi_obesity %in% c("Normal", "Overweight") ~ "Non_Obese",
      bmi_obesity == "Obese" ~ "Obese",
      TRUE ~ NA_character_
    ))
}

if (isTRUE(config$needs_binary_dxa)) {
  pb90_subset_clean@meta.data <- pb90_subset_clean@meta.data %>%
    mutate(dxa_obese_binary = case_when(
      dxa_obesity %in% c("Normal", "Overweight") ~ "Non_Obese",
      dxa_obesity == "Obese" ~ "Obese",
      TRUE ~ NA_character_
    ))
}

# Build the full subset expression including celltype filter
full_subset_cond <- sprintf("%s == '%s' & %s", celltype_var, celltype_group_input, config$subset_cond)
cat(sprintf("\nApplying subset condition: %s\n", full_subset_cond))

# Apply subset
meta <- pb90_subset_clean@meta.data
keep_cells <- tryCatch({
  with(meta, expr = eval(parse(text = full_subset_cond)))
}, error = function(e) {
  cat(sprintf("ERROR in subset condition: %s\n", e$message))
  stop(e)
})

# Handle NAs in keep_cells
keep_cells[is.na(keep_cells)] <- FALSE

pb90_celltype <- subset(pb90_subset_clean, cells = rownames(meta)[keep_cells])

n_cells_after_subset <- ncol(pb90_celltype)
cat(sprintf("Cells after subsetting: %d\n", n_cells_after_subset))

# Check minimum cells
if (n_cells_after_subset < 30) {
  cat(sprintf("SKIPPING: Too few cells (%d < 30) for %s - %s\n",
              n_cells_after_subset, analysis_type, celltype_group_input))
  quit(save = "no", status = 0)
}

# Remove cells with NA in adjust_covariates (for age-adjusted categorical/continuous models)
if (!is.null(config$adjust_covariates)) {
  cov_meta <- pb90_celltype@meta.data
  cov_complete <- complete.cases(cov_meta[, config$adjust_covariates, drop = FALSE])
  n_before <- ncol(pb90_celltype)
  pb90_celltype <- subset(pb90_celltype, cells = rownames(cov_meta)[cov_complete])
  n_after <- ncol(pb90_celltype)
  cat(sprintf("Removed %d cells with NA in covariates (%s): %d -> %d cells\n",
              n_before - n_after, paste(config$adjust_covariates, collapse = ", "),
              n_before, n_after))
  
  # Re-check minimum cells after covariate NA removal
  if (n_after < 30) {
    cat(sprintf("SKIPPING: Too few cells (%d < 30) after covariate NA removal for %s - %s\n",
                n_after, analysis_type, celltype_group_input))
    quit(save = "no", status = 0)
  }
}

# For continuous analyses, also remove NAs in the continuous variable
if (config$analysis_mode == "continuous") {
  cont_var <- config$continuous_var
  cov_meta <- pb90_celltype@meta.data
  has_value <- !is.na(cov_meta[[cont_var]])
  
  # If adjusting for group, also ensure group is not NA
  if (isTRUE(config$adjust_group)) {
    has_value <- has_value & !is.na(cov_meta[["group"]])
  }
  
  n_before <- ncol(pb90_celltype)
  pb90_celltype <- subset(pb90_celltype, cells = rownames(cov_meta)[has_value])
  n_after <- ncol(pb90_celltype)
  cat(sprintf("Removed %d cells with NA in %s%s: %d -> %d cells\n",
              n_before - n_after, cont_var,
              ifelse(isTRUE(config$adjust_group), " or group", ""),
              n_before, n_after))
  
  # Check for sufficient variance in continuous variable
  unique_vals <- length(unique(pb90_celltype@meta.data[[cont_var]]))
  if (unique_vals < 3) {
    cat(sprintf("SKIPPING: Insufficient variance in %s (%d unique values) for %s\n",
                cont_var, unique_vals, celltype_group_input))
    quit(save = "no", status = 0)
  }
  
  # Scale continuous variable for numerical stability
  pb90_celltype@meta.data[[cont_var]] <- scale(pb90_celltype@meta.data[[cont_var]])[,1]
  cat(sprintf("Scaled %s (mean-centered, SD=1)\n", cont_var))
}

# For interaction analyses, remove NAs in continuous variable and group, then scale
if (config$analysis_mode == "interaction") {
  cont_var <- config$continuous_var
  cov_meta <- pb90_celltype@meta.data
  has_value <- !is.na(cov_meta[[cont_var]]) & !is.na(cov_meta[["group"]])
  
  n_before <- ncol(pb90_celltype)
  pb90_celltype <- subset(pb90_celltype, cells = rownames(cov_meta)[has_value])
  n_after <- ncol(pb90_celltype)
  cat(sprintf("Removed %d cells with NA in %s or group: %d -> %d cells\n",
              n_before - n_after, cont_var, n_before, n_after))
  
  # Check for sufficient variance in continuous variable
  unique_vals <- length(unique(pb90_celltype@meta.data[[cont_var]]))
  if (unique_vals < 3) {
    cat(sprintf("SKIPPING: Insufficient variance in %s (%d unique values) for %s\n",
                cont_var, unique_vals, celltype_group_input))
    quit(save = "no", status = 0)
  }
  
  # Scale continuous variable for numerical stability
  pb90_celltype@meta.data[[cont_var]] <- scale(pb90_celltype@meta.data[[cont_var]])[,1]
  cat(sprintf("Scaled %s (mean-centered, SD=1)\n", cont_var))
}

# Check minimum subjects
n_unique_subjects <- length(unique(pb90_celltype@meta.data$record_id))
if (n_unique_subjects < 3) {
  cat(sprintf("SKIPPING: Too few subjects (%d < 3) for %s - %s\n",
              n_unique_subjects, analysis_type, celltype_group_input))
  quit(save = "no", status = 0)
}

# =============================================================================
# LOG DIAGNOSTIC INFORMATION
# =============================================================================
cat(sprintf("\n=============================================================\n"))
cat(sprintf("DIAGNOSTIC SUMMARY\n"))
cat(sprintf("=============================================================\n"))
cat(sprintf("Analysis: %s\n", analysis_type))
cat(sprintf("Cell type: %s (%s)\n", celltype_group_input, celltype_var))
cat(sprintf("Analysis mode: %s\n", config$analysis_mode))
cat(sprintf("Total cells: %d\n", ncol(pb90_celltype)))
cat(sprintf("Unique record_ids: %d\n", n_unique_subjects))

# Log group distribution
subject_meta <- pb90_celltype@meta.data %>%
  distinct(record_id, .keep_all = TRUE)

cat(sprintf("\nDisease group distribution (subjects):\n"))
print(table(subject_meta$group, useNA = "ifany"))

cat(sprintf("\nBMI obesity distribution (subjects):\n"))
if ("bmi_obesity" %in% colnames(subject_meta)) {
  print(table(subject_meta$bmi_obesity, useNA = "ifany"))
} else {
  cat("  bmi_obesity column not found\n")
}

cat(sprintf("\nDXA obesity distribution (subjects):\n"))
if ("dxa_obesity" %in% colnames(subject_meta)) {
  print(table(subject_meta$dxa_obesity, useNA = "ifany"))
} else {
  cat("  dxa_obesity column not found\n")
}

cat(sprintf("\nSex distribution (subjects):\n"))
print(table(subject_meta$sex, useNA = "ifany"))

# For categorical analyses, log group var distribution
if (config$analysis_mode == "categorical") {
  group_var <- config$group_var
  cat(sprintf("\nGroup variable '%s' distribution (subjects):\n", group_var))
  print(table(subject_meta[[group_var]], useNA = "ifany"))
  
  cat(sprintf("\nGroup variable '%s' distribution (cells):\n", group_var))
  print(table(pb90_celltype@meta.data[[group_var]], useNA = "ifany"))
  
  if (!is.null(config$adjust_covariates)) {
    cat(sprintf("\nAdjusting for covariates: %s\n", paste(config$adjust_covariates, collapse = ", ")))
    for (cov in config$adjust_covariates) {
      cat(sprintf("  Covariate '%s' summary (subjects):\n", cov))
      print(summary(subject_meta[[cov]]))
    }
  }
}

# For continuous, log summary stats
if (config$analysis_mode == "continuous") {
  cont_var <- config$continuous_var
  cat(sprintf("\nContinuous variable '%s' summary (subjects, after scaling):\n", cont_var))
  print(summary(subject_meta[[cont_var]]))
  if (isTRUE(config$adjust_group)) {
    cat("  Adjusting for: group\n")
  }
}

# For interaction, log continuous var summary and group distribution
if (config$analysis_mode == "interaction") {
  cont_var <- config$continuous_var
  cat(sprintf("\nContinuous variable '%s' summary (subjects, after scaling):\n", cont_var))
  print(summary(subject_meta[[cont_var]]))
  cat(sprintf("\nGroup distribution for interaction (subjects):\n"))
  print(table(subject_meta$group, useNA = "ifany"))
  cat(sprintf("\nGroup distribution for interaction (cells):\n"))
  print(table(pb90_celltype@meta.data$group, useNA = "ifany"))
}

# Cells per subject
cells_per_subject <- pb90_celltype@meta.data %>%
  group_by(record_id) %>%
  summarise(n_cells = n(), .groups = "drop")
cat(sprintf("\nCells per subject summary:\n"))
print(summary(cells_per_subject$n_cells))

cat(sprintf("=============================================================\n\n"))

# =============================================================================
# SET UP FORMULA AND FACTOR LEVELS
# =============================================================================
if (config$analysis_mode == "categorical") {
  group_var <- config$group_var
  pb90_celltype@meta.data[[group_var]] <- factor(pb90_celltype@meta.data[[group_var]])
  
  if (!is.null(config$ref_level)) {
    pb90_celltype@meta.data[[group_var]] <- relevel(
      pb90_celltype@meta.data[[group_var]],
      ref = config$ref_level
    )
  }
  
  # Check minimum group sizes
  group_counts <- pb90_celltype@meta.data %>%
    distinct(record_id, .keep_all = TRUE) %>%
    dplyr::count(!!sym(group_var))
  
  min_group <- min(group_counts$n)
  if (min_group < 2) {
    cat(sprintf("SKIPPING: A group has fewer than 2 subjects (min=%d) for %s - %s\n",
                min_group, analysis_type, celltype_group_input))
    print(group_counts)
    quit(save = "no", status = 0)
  }
  
  # Build formula: ~ group_var [+ covariate1 + covariate2 ...]
  if (!is.null(config$adjust_covariates)) {
    formula_obj <- as.formula(paste("~", group_var, "+",
                                    paste(config$adjust_covariates, collapse = " + ")))
  } else {
    formula_obj <- as.formula(paste("~", group_var))
  }
  cat(sprintf("Formula: %s\n", deparse(formula_obj)))
  
} else if (config$analysis_mode == "interaction") {
  # Interaction model: ~ continuous_var * group
  cont_var <- config$continuous_var
  pb90_celltype@meta.data[["group"]] <- factor(pb90_celltype@meta.data[["group"]])
  
  # Check minimum group sizes for interaction
  group_counts <- pb90_celltype@meta.data %>%
    distinct(record_id, .keep_all = TRUE) %>%
    dplyr::count(group)
  
  min_group <- min(group_counts$n)
  if (min_group < 2) {
    cat(sprintf("SKIPPING: A group has fewer than 2 subjects (min=%d) for %s - %s\n",
                min_group, analysis_type, celltype_group_input))
    print(group_counts)
    quit(save = "no", status = 0)
  }
  
  formula_obj <- as.formula(paste("~", cont_var, "* group"))
  cat(sprintf("Formula: %s\n", deparse(formula_obj)))
  
} else {
  # Continuous analysis
  cont_var <- config$continuous_var
  formula_parts <- cont_var
  if (isTRUE(config$adjust_group)) {
    pb90_celltype@meta.data[["group"]] <- factor(pb90_celltype@meta.data[["group"]])
    formula_parts <- c(formula_parts, "group")
  }
  if (!is.null(config$adjust_covariates)) {
    formula_parts <- c(formula_parts, config$adjust_covariates)
  }
  formula_obj <- as.formula(paste("~", paste(formula_parts, collapse = " + ")))
  cat(sprintf("Formula: %s\n", deparse(formula_obj)))
}

# =============================================================================
# CALCULATE SAMPLE/CELL COUNTS FOR METADATA
# =============================================================================
if (config$analysis_mode == "categorical") {
  group_var_sym <- sym(group_var)
  
  sample_counts <- pb90_celltype@meta.data %>%
    distinct(record_id, .keep_all = TRUE) %>%
    dplyr::count(!!group_var_sym, name = "n_samples") %>%
    arrange(!!group_var_sym)
  
  cell_counts <- pb90_celltype@meta.data %>%
    dplyr::count(!!group_var_sym, name = "n_cells") %>%
    arrange(!!group_var_sym)
  
  counts_summary <- sample_counts %>%
    left_join(cell_counts, by = group_var)
} else if (config$analysis_mode == "interaction") {
  # For interaction, report counts by group
  sample_counts <- pb90_celltype@meta.data %>%
    distinct(record_id, .keep_all = TRUE) %>%
    dplyr::count(group, name = "n_samples") %>%
    arrange(group)
  
  cell_counts <- pb90_celltype@meta.data %>%
    dplyr::count(group, name = "n_cells") %>%
    arrange(group)
  
  counts_summary <- sample_counts %>%
    left_join(cell_counts, by = "group")
} else {
  # For continuous, just total counts
  counts_summary <- data.frame(
    group = "all",
    n_samples = n_unique_subjects,
    n_cells = ncol(pb90_celltype)
  )
}

cat("Counts summary:\n")
print(counts_summary)

# =============================================================================
# RUN NEBULA ANALYSIS
# =============================================================================
cat(sprintf("\nRunning NEBULA: %s for %s...\n", analysis_type, celltype_group_input))

# Build S3 paths
s3_key_raw <- sprintf("%s/%s/%s/%s_nebula_%s.rds",
                      s3_base, config$s3_subdir, celltype_group_input,
                      celltype_group_input, config$file_suffix)
s3_key_processed <- sprintf("%s/%s/%s/%s_nebula_%s_processed.rds",
                            s3_base, config$s3_subdir, celltype_group_input,
                            celltype_group_input, config$file_suffix)

# Run nebula
nebula_res <- tryCatch({
  run_nebula_parallel(
    pb90_celltype,
    n_cores = 15,
    subject_id_col = "record_id",
    formula = formula_obj,
    s3_bucket = s3_bucket,
    s3_key = s3_key_raw,
    group = FALSE
  )
}, error = function(e) {
  cat(sprintf("ERROR running NEBULA: %s\n", e$message))
  cat("Attempting to continue with error handling...\n")
  NULL
})

if (is.null(nebula_res)) {
  cat("NEBULA returned NULL - analysis failed. Exiting.\n")
  quit(save = "no", status = 1)
}

# Process results
processed <- process_nebula_results(
  nebula_res$results,
  pval_col = config$pval_col,
  logfc_col = config$logfc_col
)

# =============================================================================
# CREATE AND SAVE RUN METADATA
# =============================================================================
cat("Creating run metadata...\n")

if (config$analysis_mode == "categorical") {
  samples_per_group <- counts_summary %>%
    mutate(label = paste0(!!group_var_sym, ": ", n_samples)) %>%
    pull(label) %>%
    paste(collapse = "; ")
  
  cells_per_group <- counts_summary %>%
    mutate(label = paste0(!!group_var_sym, ": ", n_cells)) %>%
    pull(label) %>%
    paste(collapse = "; ")
} else if (config$analysis_mode == "interaction") {
  samples_per_group <- counts_summary %>%
    mutate(label = paste0(group, ": ", n_samples)) %>%
    pull(label) %>%
    paste(collapse = "; ")
  
  cells_per_group <- counts_summary %>%
    mutate(label = paste0(group, ": ", n_cells)) %>%
    pull(label) %>%
    paste(collapse = "; ")
} else {
  samples_per_group <- paste0("total: ", n_unique_subjects)
  cells_per_group <- paste0("total: ", ncol(pb90_celltype))
}

# Determine the primary variable name for metadata
primary_variable <- if (config$analysis_mode == "categorical") {
  config$group_var
} else {
  config$continuous_var
}

run_metadata <- data.frame(
  analysis_type = analysis_type,
  analysis_mode = config$analysis_mode,
  celltype = celltype_group_input,
  celltype_variable = celltype_var,
  group_variable = primary_variable,
  reference_level = ifelse(is.null(config$ref_level), NA, config$ref_level),
  adjust_group = ifelse(is.null(config$adjust_group), FALSE, config$adjust_group),
  adjust_covariates = ifelse(is.null(config$adjust_covariates), NA,
                             paste(config$adjust_covariates, collapse = ", ")),
  formula_used = deparse(formula_obj),
  pval_column = config$pval_col,
  logfc_column = config$logfc_col,
  n_samples_per_group = samples_per_group,
  n_cells_per_group = cells_per_group,
  total_samples = if (config$analysis_mode == "categorical") {
    sum(counts_summary$n_samples)
  } else if (config$analysis_mode == "interaction") {
    sum(counts_summary$n_samples)
  } else {
    n_unique_subjects
  },
  total_cells = ncol(pb90_celltype),
  subset_condition = full_subset_cond,
  s3_results_key = s3_key_processed,
  run_timestamp = as.character(Sys.time()),
  stringsAsFactors = FALSE
)

# Save metadata
s3_key_metadata <- sprintf("%s/%s/%s/%s_nebula_%s_metadata.csv",
                           s3_base, config$s3_subdir, celltype_group_input,
                           celltype_group_input, config$file_suffix)

cat(sprintf("Saving run metadata to S3: %s\n", s3_key_metadata))

temp_csv <- tempfile(fileext = ".csv")
write.csv(run_metadata, temp_csv, row.names = FALSE)
put_object(file = temp_csv, object = s3_key_metadata, bucket = s3_bucket, region = "")
unlink(temp_csv)

# =============================================================================
# ANNOTATE AND SAVE RESULTS
# =============================================================================
cat("Adding gene annotations...\n")

gene_info <- tryCatch({
  getBM(
    attributes = c("hgnc_symbol", "description", "gene_biotype"),
    filters = "hgnc_symbol",
    values = processed$results$Gene,
    mart = mart
  ) %>%
    dplyr::rename(Gene = hgnc_symbol)
}, error = function(e) {
  cat(sprintf("WARNING: biomaRt annotation failed: %s\n", e$message))
  cat("Saving results without annotations.\n")
  data.frame(Gene = character(0), description = character(0), gene_biotype = character(0))
})

annotated_df <- processed$results %>%
  left_join(gene_info, by = "Gene")

cat(sprintf("Saving processed results to S3: %s\n", s3_key_processed))
s3saveRDS(annotated_df, object = s3_key_processed, bucket = s3_bucket, region = "")

cat(sprintf("\n=============================================================\n"))
cat(sprintf("COMPLETE: %s - %s at %s\n", analysis_type, celltype_group_input, Sys.time()))
cat(sprintf("Total significant genes (FDR < 0.05): %d\n",
            sum(annotated_df$p.adjust < 0.05, na.rm = TRUE)))
cat(sprintf("Results saved to: %s/%s\n", s3_bucket, s3_key_processed))
cat(sprintf("=============================================================\n"))