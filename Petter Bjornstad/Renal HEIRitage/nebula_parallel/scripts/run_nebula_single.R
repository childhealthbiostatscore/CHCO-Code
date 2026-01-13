#!/usr/bin/env Rscript
# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript run_nebula_single.R <analysis_type> <celltype> [celltype_var]")
}

analysis_type <- args[1]
celltype_group_input <- args[2]
celltype_var <- if (length(args) == 3) args[3] else "KPMP_celltype_general"

# Load libraries
library(aws.s3)
library(jsonlite)
library(biomaRt)
library(Seurat)
library(dplyr)

# analysis_type <- "DKD_vs_nonDKD_100"
# celltype_group_input <- "PT-S1/S2"
# celltype_var <- "KPMP_celltype"

cat(sprintf("Starting analysis: %s for cell type: %s (using %s) at %s\n", 
            analysis_type, celltype_group_input, celltype_var, Sys.time()))

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

# =============================================================================
# ANALYSIS CONFIGURATION
# =============================================================================
# Define all analysis configurations in one place for consistency
# subset_cond: filtering conditions EXCLUDING celltype (that's added dynamically)

analysis_config <- list(
  # T2D GLP- vs. HC
  T2D_GLP_N_vs_HC = list(
    subset_cond = "group %in% c('Type_2_Diabetes', 'Lean_Control') & !is.na(glp_t2dob) & glp_t2dob != 'GLP_Y'",
    group_var = "glp_t2dob", ref_level = "HC",
    pval_col = "p_glp_t2dobGLP_N", logfc_col = "logFC_glp_t2dobGLP_N",
    s3_subdir = "T2D_GLP_N_vs_HC", file_suffix = "t2d_glpn_hc"
  ),
  # T2D GLP+ vs T2D GLP-
  T2D_GLP_Y_vs_T2D_GLP_N = list(
    subset_cond = "group == 'Type_2_Diabetes' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "T2D_GLP_Y_vs_T2D_GLP_N", file_suffix = "t2d_glpyn"
  ),
  # DKD vs nonDKD comparisons (ACR >= 100)
  DKD_vs_nonDKD_100 = list(
    subset_cond = "group != 'Lean_Control' & !is.na(dkd_group_100)",
    group_var = "dkd_group_100", ref_level = "non_DKD",
    pval_col = "p_dkd_group_100DKD", logfc_col = "logFC_dkd_group_100DKD",
    s3_subdir = "DKD_vs_nonDKD_100", file_suffix = "dkd100"
  ),
  DKD_vs_nonDKD_100_t2d = list(
    subset_cond = "group == 'Type_2_Diabetes' & !is.na(dkd_group_100)",
    group_var = "dkd_group_100", ref_level = "non_DKD",
    pval_col = "p_dkd_group_100DKD", logfc_col = "logFC_dkd_group_100DKD",
    s3_subdir = "DKD_vs_nonDKD_100_t2d", file_suffix = "dkd100_t2d"
  ),
  DKD_vs_nonDKD_100_nosglt2i = list(
    subset_cond = "group == 'Type_2_Diabetes' & epic_sglti2_1 == 'No' & !is.na(dkd_group_100)",
    group_var = "dkd_group_100", ref_level = "non_DKD",
    pval_col = "p_dkd_group_100DKD", logfc_col = "logFC_dkd_group_100DKD",
    s3_subdir = "DKD_vs_nonDKD_100_nosglt2i", file_suffix = "dkd100_nosglt2i"
  ),
  
  # DKD vs nonDKD comparisons (ACR >= 30)
  DKD_vs_nonDKD_30 = list(
    subset_cond = "group != 'Lean_Control' & !is.na(dkd_group_30)",
    group_var = "dkd_group_30", ref_level = "non_DKD",
    pval_col = "p_dkd_group_30DKD", logfc_col = "logFC_dkd_group_30DKD",
    s3_subdir = "DKD_vs_nonDKD_30", file_suffix = "dkd30"
  ),
  DKD_vs_nonDKD_30_t2d = list(
    subset_cond = "group == 'Type_2_Diabetes' & !is.na(dkd_group_30)",
    group_var = "dkd_group_30", ref_level = "non_DKD",
    pval_col = "p_dkd_group_30DKD", logfc_col = "logFC_dkd_group_30DKD",
    s3_subdir = "DKD_vs_nonDKD_30_t2d", file_suffix = "dkd30_t2d"
  ),
  DKD_vs_nonDKD_30_nosglt2i = list(
    subset_cond = "group == 'Type_2_Diabetes' & epic_sglti2_1 == 'No' & !is.na(dkd_group_30)",
    group_var = "dkd_group_30", ref_level = "non_DKD",
    pval_col = "p_dkd_group_30DKD", logfc_col = "logFC_dkd_group_30DKD",
    s3_subdir = "DKD_vs_nonDKD_30_nosglt2i", file_suffix = "dkd30_nosglt2i"
  ),
  
  # DKD vs HC comparisons (ACR >= 100)
  DKD_100_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcDKD", logfc_col = "logFC_dkd_group_100_hcDKD",
    s3_subdir = "DKD_100_vs_HC", file_suffix = "dkd100_hc"
  ),
  DKD_100_t2d_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & group != 'Obese_Control' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcDKD", logfc_col = "logFC_dkd_group_100_hcDKD",
    s3_subdir = "DKD_100_t2d_vs_HC", file_suffix = "dkd100_t2d_hc"
  ),
  nonDKD_100_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_vs_HC", file_suffix = "nondkd100_hc"
  ),
  nonDKD_100_t2d_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & group != 'Obese_Control' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_t2d_vs_HC", file_suffix = "nondkd100_t2d_hc"
  ),
  
  # DKD vs HC comparisons (ACR >= 30)
  DKD_30_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcDKD", logfc_col = "logFC_dkd_group_30_hcDKD",
    s3_subdir = "DKD_30_vs_HC", file_suffix = "dkd30_hc"
  ),
  DKD_30_t2d_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & group != 'Obese_Control' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcDKD", logfc_col = "logFC_dkd_group_30_hcDKD",
    s3_subdir = "DKD_30_t2d_vs_HC", file_suffix = "dkd30_t2d_hc"
  ),
  nonDKD_30_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_vs_HC", file_suffix = "nondkd30_hc"
  ),
  nonDKD_30_t2d_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & group != 'Obese_Control' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_t2d_vs_HC", file_suffix = "nondkd30_t2d_hc"
  ),
  
  # GLP-1RA comparisons
  GLP_N_vs_HC = list(
    subset_cond = "glp_t2dob != 'GLP_Y' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "HC",
    pval_col = "p_glp_t2dobGLP_N", logfc_col = "logFC_glp_t2dobGLP_N",
    s3_subdir = "GLP_N_vs_HC", file_suffix = "glpn_hc"
  ),
  GLP_Y_vs_GLP_N = list(
    subset_cond = "group != 'Lean_Control' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "GLP_Y_vs_GLP_N", file_suffix = "glpy_glpn"
  ),
  
  # GLP within DKD/nonDKD subgroups (ACR >= 100)
  DKD_100_GLP_Y_vs_DKD_100_GLP_N = list(
    subset_cond = "dkd_group_100_hc == 'DKD' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "DKD_100_GLP_Y_vs_DKD_100_GLP_N", file_suffix = "dkd_100_glpy_glpn"
  ),
  nonDKD_100_GLP_Y_vs_nonDKD_100_GLP_N = list(
    subset_cond = "dkd_group_100_hc == 'non_DKD' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "nonDKD_100_GLP_Y_vs_nonDKD_100_GLP_N", file_suffix = "nondkd_100_glpy_glpn"
  ),
  nonDKD_100_GLP_N_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & glp_t2dob != 'GLP_Y' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_GLP_N_vs_HC", file_suffix = "nondkd100_glpn_hc"
  ),
  nonDKD_100_GLP_Y_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & glp_t2dob != 'GLP_N' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_GLP_Y_vs_HC", file_suffix = "nondkd100_glpy_hc"
  ),
  DKD_100_GLP_Y_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & glp_t2dob != 'GLP_N' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcDKD", logfc_col = "logFC_dkd_group_100_hcDKD",
    s3_subdir = "DKD_100_GLP_Y_vs_HC", file_suffix = "dkd100_glpy_hc"
  ),
  DKD_100_GLP_N_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & glp_t2dob != 'GLP_Y' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcDKD", logfc_col = "logFC_dkd_group_100_hcDKD",
    s3_subdir = "DKD_100_GLP_N_vs_HC", file_suffix = "dkd100_glpn_hc"
  ),
  
  # GLP within DKD/nonDKD subgroups (ACR >= 30)
  DKD_30_GLP_Y_vs_DKD_30_GLP_N = list(
    subset_cond = "dkd_group_30_hc == 'DKD' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "DKD_30_GLP_Y_vs_DKD_30_GLP_N", file_suffix = "dkd_30_glpy_glpn"
  ),
  nonDKD_30_GLP_Y_vs_nonDKD_30_GLP_N = list(
    subset_cond = "dkd_group_30_hc == 'non_DKD' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "nonDKD_30_GLP_Y_vs_nonDKD_30_GLP_N", file_suffix = "nondkd_30_glpy_glpn"
  ),
  nonDKD_30_GLP_N_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & glp_t2dob != 'GLP_Y' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_GLP_N_vs_HC", file_suffix = "nondkd30_glpn_hc"
  ),
  nonDKD_30_GLP_Y_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & glp_t2dob != 'GLP_N' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_GLP_Y_vs_HC", file_suffix = "nondkd30_glpy_hc"
  ),
  DKD_30_GLP_Y_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & glp_t2dob != 'GLP_N' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcDKD", logfc_col = "logFC_dkd_group_30_hcDKD",
    s3_subdir = "DKD_30_GLP_Y_vs_HC", file_suffix = "dkd30_glpy_hc"
  ),
  DKD_30_GLP_N_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & glp_t2dob != 'GLP_Y' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcDKD", logfc_col = "logFC_dkd_group_30_hcDKD",
    s3_subdir = "DKD_30_GLP_N_vs_HC", file_suffix = "dkd30_glpn_hc"
  ),
  
  # SGLT2i comparisons
  SGLT2i_N_vs_HC = list(
    subset_cond = "sglt2i_t2dob != 'SGLT2i_Y' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_N", logfc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    s3_subdir = "SGLT2i_N_vs_HC", file_suffix = "sglt2in_hc"
  ),
  SGLT2i_Y_vs_SGLT2i_N = list(
    subset_cond = "group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "SGLT2i_Y_vs_SGLT2i_N", file_suffix = "sglt2iy_sglt2in"
  ),
  DKD_30_SGLT2i_N_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & sglt2i_t2dob != 'SGLT2i_Y' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_N", logfc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    s3_subdir = "DKD_30_SGLT2i_N_vs_HC", file_suffix = "dkd30_sglt2in_hc"
  ),
  DKD_30_SGLT2i_Y_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & sglt2i_t2dob != 'SGLT2i_N' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "DKD_30_SGLT2i_Y_vs_HC", file_suffix = "dkd30_sglt2iy_hc"
  ),
  DKD_30_SGLT2i_Y_vs_DKD_30_SGLT2i_N = list(
    subset_cond = "dkd_group_30_hc == 'DKD' & group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "DKD_30_SGLT2i_Y_vs_DKD_30_SGLT2i_N", file_suffix = "dkd30_sglt2iy_sglt2in"
  ),
  DKD_100_SGLT2i_Y_vs_DKD_100_SGLT2i_N = list(
    subset_cond = "dkd_group_100_hc == 'DKD' & group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "DKD_100_SGLT2i_Y_vs_DKD_100_SGLT2i_N", file_suffix = "dkd100_sglt2iy_sglt2in"
  ),
  nonDKD_100_SGLT2i_N_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & sglt2i_t2dob != 'SGLT2i_Y' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_SGLT2i_N_vs_HC", file_suffix = "nondkd100_sglt2in_hc"
  ),
  nonDKD_100_SGLT2i_Y_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & sglt2i_t2dob != 'SGLT2i_N' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_SGLT2i_Y_vs_HC", file_suffix = "nondkd100_sglt2iy_hc"
  ),
  nonDKD_30_SGLT2i_N_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & sglt2i_t2dob != 'SGLT2i_Y' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_SGLT2i_N_vs_HC", file_suffix = "nondkd30_sglt2in_hc"
  ),
  nonDKD_30_SGLT2i_Y_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & sglt2i_t2dob != 'SGLT2i_N' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_SGLT2i_Y_vs_HC", file_suffix = "nondkd30_sglt2iy_hc"
  ),
  nonDKD_30_SGLT2i_Y_vs_nonDKD_30_SGLT2i_N = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "nonDKD_30_SGLT2i_Y_vs_nonDKD_30_SGLT2i_N", file_suffix = "nondkd30_sglt2iy_sglt2in"
  ),
  nonDKD_100_SGLT2i_Y_vs_nonDKD_100_SGLT2i_N = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "nonDKD_100_SGLT2i_Y_vs_nonDKD_100_SGLT2i_N", file_suffix = "nondkd100_sglt2iy_sglt2in"
  ),
  DKD_100_SGLT2i_N_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & sglt2i_t2dob != 'SGLT2i_Y' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_N", logfc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    s3_subdir = "DKD_100_SGLT2i_N_vs_HC", file_suffix = "dkd100_sglt2in_hc"
  ),
  DKD_100_SGLT2i_Y_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & sglt2i_t2dob != 'SGLT2i_N' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "DKD_100_SGLT2i_Y_vs_HC", file_suffix = "dkd100_sglt2iy_hc"
  ),
  
  # Combo comparisons - each treatment group vs HC
  Neither_vs_HC = list(
    subset_cond = "combo_t2dob %in% c('Neither', 'HC') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_combo_t2dobNeither", logfc_col = "logFC_combo_t2dobNeither",
    s3_subdir = "Neither_vs_HC", file_suffix = "neither_hc"
  ),
  SGLT2i_only_vs_HC = list(
    subset_cond = "combo_t2dob %in% c('SGLT2i_only', 'HC') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_combo_t2dobSGLT2i_only", logfc_col = "logFC_combo_t2dobSGLT2i_only",
    s3_subdir = "SGLT2i_only_vs_HC", file_suffix = "sglt2ionly_hc"
  ),
  GLP_only_vs_HC = list(
    subset_cond = "combo_t2dob %in% c('GLP_only', 'HC') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_combo_t2dobGLP_only", logfc_col = "logFC_combo_t2dobGLP_only",
    s3_subdir = "GLP_only_vs_HC", file_suffix = "glponly_hc"
  ),
  Both_vs_HC = list(
    subset_cond = "combo_t2dob %in% c('Both', 'HC') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_combo_t2dobBoth", logfc_col = "logFC_combo_t2dobBoth",
    s3_subdir = "Both_vs_HC", file_suffix = "both_hc"
  ),
  
  # Combo comparisons - treatment groups vs Neither (no treatment)
  SGLT2i_only_vs_Neither = list(
    subset_cond = "combo_t2dob %in% c('SGLT2i_only', 'Neither') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "Neither", needs_combo = TRUE,
    pval_col = "p_combo_t2dobSGLT2i_only", logfc_col = "logFC_combo_t2dobSGLT2i_only",
    s3_subdir = "SGLT2i_only_vs_Neither", file_suffix = "sglt2ionly_neither"
  ),
  GLP_only_vs_Neither = list(
    subset_cond = "combo_t2dob %in% c('GLP_only', 'Neither') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "Neither", needs_combo = TRUE,
    pval_col = "p_combo_t2dobGLP_only", logfc_col = "logFC_combo_t2dobGLP_only",
    s3_subdir = "GLP_only_vs_Neither", file_suffix = "glponly_neither"
  ),
  Both_vs_Neither = list(
    subset_cond = "combo_t2dob %in% c('Both', 'Neither') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "Neither", needs_combo = TRUE,
    pval_col = "p_combo_t2dobBoth", logfc_col = "logFC_combo_t2dobBoth",
    s3_subdir = "Both_vs_Neither", file_suffix = "both_neither"
  ),
  
  # Combo comparisons - Both vs single treatments
  Both_vs_SGLT2i_only = list(
    subset_cond = "combo_t2dob %in% c('Both', 'SGLT2i_only') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "SGLT2i_only", needs_combo = TRUE,
    pval_col = "p_combo_t2dobBoth", logfc_col = "logFC_combo_t2dobBoth",
    s3_subdir = "Both_vs_SGLT2i_only", file_suffix = "both_sglt2ionly"
  ),
  Both_vs_GLP_only = list(
    subset_cond = "combo_t2dob %in% c('Both', 'GLP_only') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "GLP_only", needs_combo = TRUE,
    pval_col = "p_combo_t2dobBoth", logfc_col = "logFC_combo_t2dobBoth",
    s3_subdir = "Both_vs_GLP_only", file_suffix = "both_glponly"
  ),
  
  # Single treatments comparison
  SGLT2i_only_vs_GLP_only = list(
    subset_cond = "combo_t2dob %in% c('SGLT2i_only', 'GLP_only') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "GLP_only", needs_combo = TRUE,
    pval_col = "p_combo_t2dobSGLT2i_only", logfc_col = "logFC_combo_t2dobSGLT2i_only",
    s3_subdir = "SGLT2i_only_vs_GLP_only", file_suffix = "sglt2ionly_glponly"
  ),
  
  # Combo comparisons within nonDKD (ACR >= 100) - vs HC
  nonDKD_100_Neither_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & combo_t2dob %in% c('Neither', 'HC') & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_Neither_vs_HC", file_suffix = "nondkd100_neither_hc"
  ),
  nonDKD_100_SGLT2i_only_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & combo_t2dob %in% c('SGLT2i_only', 'HC') & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_SGLT2i_only_vs_HC", file_suffix = "nondkd100_sglt2ionly_hc"
  ),
  nonDKD_100_GLP_only_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & combo_t2dob %in% c('GLP_only', 'HC') & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_GLP_only_vs_HC", file_suffix = "nondkd100_glponly_hc"
  ),
  nonDKD_100_Both_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & combo_t2dob %in% c('Both', 'HC') & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_Both_vs_HC", file_suffix = "nondkd100_both_hc"
  ),
  
  # Combo comparisons within nonDKD (ACR >= 30) - vs HC
  nonDKD_30_Neither_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & combo_t2dob %in% c('Neither', 'HC') & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_Neither_vs_HC", file_suffix = "nondkd30_neither_hc"
  ),
  nonDKD_30_SGLT2i_only_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & combo_t2dob %in% c('SGLT2i_only', 'HC') & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_SGLT2i_only_vs_HC", file_suffix = "nondkd30_sglt2ionly_hc"
  ),
  nonDKD_30_GLP_only_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & combo_t2dob %in% c('GLP_only', 'HC') & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_GLP_only_vs_HC", file_suffix = "nondkd30_glponly_hc"
  ),
  nonDKD_30_Both_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & combo_t2dob %in% c('Both', 'HC') & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_Both_vs_HC", file_suffix = "nondkd30_both_hc"
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
cat(sprintf("Loading cell type-specific data for %s...\n", celltype_group_input))

subset_file <- paste0("data_clean/subset/pb90_ckd_analysis/pb90_ckd_analysis_subset_", 
                      celltype_group_input, ".rds")

pb90_subset_clean <- tryCatch({
  cat(sprintf("Looking for: %s\n", subset_file))
  obj <- s3readRDS(object = subset_file, bucket = "scrna", region = "")
  cat(sprintf("Loaded cell type-specific subset with %d cells\n", ncol(obj)))
  obj
}, error = function(e) {
  cat(sprintf("Cell type-specific file not found: %s\n", subset_file))
  cat("Falling back to loading full dataset and subsetting...\n")
  obj <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset.rds", 
                   bucket = "scrna", region = "")
  cat(sprintf("Loaded full dataset with %d cells\n", ncol(obj)))
  obj <- subset(obj, !!sym(celltype_var) == celltype_group_input)
  cat(sprintf("Subset to %s: %d cells\n", celltype_group_input, ncol(obj)))
  obj
})

celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),   # grouped thin limbs
  DCT_CNT = c("DCT", "dDCT", "CNT"),   # grouped distal tubule/connecting tubule
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A"), 
  # Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  Immune_Myeloid = c("MAC", "MON", "cDC", "pDC"),
  Immune_Lymphoid = c("CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",                         # mesangial cells
  PEC = "PEC",                       # parietal epithelial cells
  Schwann = "SchwannCells",
  Other = c("non-specific")          # catchall
)

map_celltype_to_general <- function(celltype, celltype_groups) {
  for (group_name in names(celltype_groups)) {
    if (celltype %in% celltype_groups[[group_name]]) {
      return(group_name)
    }
  }
  return("Other")  # For any celltype not found in groups
}

pb90_subset_clean$KPMP_celltype_general2 <- sapply(pb90_subset_clean$KPMP_celltype, 
                                            map_celltype_to_general, 
                                            celltype_groups = celltype_groups)

# Load harmonized dataset and source functions
rh_rh2_croc_improve_unique <- s3readRDS(
  object = "data_clean/rh_rh2_croc_improve_unique.RDS", 
  bucket = "harmonized.dataset", region = ""
)
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_scRNA_functions.R"))
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_functions.R"))

# Initialize biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# =============================================================================
# PREPARE DATA FOR ANALYSIS
# =============================================================================
gc()

# Add derived variables if needed
if (isTRUE(config$needs_sglt2i)) {
  pb90_subset_clean@meta.data <- pb90_subset_clean@meta.data %>%
    mutate(sglt2i_t2dob = case_when(
      group == "Lean_Control" ~ "HC",
      epic_sglti2_1 == "Yes" ~ "SGLT2i_Y",
      epic_sglti2_1 == "No" ~ "SGLT2i_N"
    ))
}

if (isTRUE(config$needs_combo)) {
  pb90_subset_clean@meta.data <- pb90_subset_clean@meta.data %>%
    mutate(combo_t2dob = case_when(
      group == "Lean_Control" ~ "HC",
      epic_sglti2_1 == "Yes" & epic_glp1ra_1 == "Yes" ~ "Both",
      epic_sglti2_1 == "Yes" & epic_glp1ra_1 == "No" ~ "SGLT2i_only",
      epic_sglti2_1 == "No" & epic_glp1ra_1 == "Yes" ~ "GLP_only",
      epic_sglti2_1 == "No" & epic_glp1ra_1 == "No" ~ "Neither"
    ))
}

# Build the full subset expression including celltype filter
# The celltype filter uses the user-specified celltype_var column
full_subset_cond <- sprintf("%s == '%s' & %s", celltype_var, celltype_group_input, config$subset_cond)
cat(sprintf("Applying subset condition: %s\n", full_subset_cond))

# Build the subset expression
meta <- pb90_subset_clean@meta.data

# Build logical vector directly
keep_cells <- with(meta, expr = eval(parse(text = full_subset_cond)))

pb90_celltype <- subset(pb90_subset_clean, cells = rownames(meta)[keep_cells])

# Set factor levels
group_var <- config$group_var
pb90_celltype@meta.data[[group_var]] <- factor(pb90_celltype@meta.data[[group_var]])
if (!is.null(config$ref_level)) {
  pb90_celltype@meta.data[[group_var]] <- relevel(
    pb90_celltype@meta.data[[group_var]], 
    ref = config$ref_level
  )
}

# Calculate sample and cell counts per group
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

# Print sample counts
cat(sprintf("Sample counts for %s:\n", celltype_group_input))
print(counts_summary)

# =============================================================================
# RUN NEBULA ANALYSIS
# =============================================================================
cat(sprintf("Running %s analysis for %s...\n", analysis_type, celltype_group_input))

# Build S3 paths
s3_base <- "Projects/CKD/RH_RH2/Results/nebula"
s3_key_raw <- sprintf("%s/%s/%s/%s_rh_rh2_imp_nebula_kpmp_%s.rds",
                      s3_base, config$s3_subdir, celltype_group_input,
                      celltype_group_input, config$file_suffix)
s3_key_processed <- sprintf("%s/%s/%s/%s_rh_rh2_imp_nebula_kpmp_%s_processed.rds",
                            s3_base, config$s3_subdir, celltype_group_input,
                            celltype_group_input, config$file_suffix)

# Build formula
formula_obj <- as.formula(paste("~", group_var))

# Run nebula
nebula_res <- run_nebula_parallel(
  pb90_celltype,
  n_cores = 15,
  subject_id_col = "record_id",
  formula = formula_obj,
  s3_bucket = "scrna",
  s3_key = s3_key_raw,
  group = FALSE
)

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

# Format counts for metadata
samples_per_group <- counts_summary %>%
  mutate(label = paste0(!!group_var_sym, ": ", n_samples)) %>%
  pull(label) %>%
  paste(collapse = "; ")

cells_per_group <- counts_summary %>%
  mutate(label = paste0(!!group_var_sym, ": ", n_cells)) %>%
  pull(label) %>%
  paste(collapse = "; ")

# Create metadata dataframe
run_metadata <- data.frame(
  analysis_type = analysis_type,
  celltype = celltype_group_input,
  celltype_variable = celltype_var,
  group_variable = group_var,
  reference_level = ifelse(is.null(config$ref_level), "default (alphabetical)", config$ref_level),
  pval_column = config$pval_col,
  logfc_column = config$logfc_col,
  n_samples_per_group = samples_per_group,
  n_cells_per_group = cells_per_group,
  total_samples = sum(counts_summary$n_samples),
  total_cells = sum(counts_summary$n_cells),
  subset_condition = full_subset_cond,
  s3_results_key = s3_key_processed,
  run_timestamp = as.character(Sys.time()),
  stringsAsFactors = FALSE
)

# Save metadata CSV to S3
s3_key_metadata <- sprintf("%s/%s/%s/%s_rh_rh2_imp_nebula_kpmp_%s_metadata.csv",
                           s3_base, config$s3_subdir, celltype_group_input,
                           celltype_group_input, config$file_suffix)

cat(sprintf("Saving run metadata to S3: %s\n", s3_key_metadata))

# Write to temp file and upload
temp_csv <- tempfile(fileext = ".csv")
write.csv(run_metadata, temp_csv, row.names = FALSE)
put_object(file = temp_csv, object = s3_key_metadata, bucket = "scrna", region = "")
unlink(temp_csv)

# =============================================================================
# ANNOTATE AND SAVE RESULTS
# =============================================================================
cat("Adding gene annotations...\n")

gene_info <- getBM(
  attributes = c("hgnc_symbol", "description", "gene_biotype"),
  filters = "hgnc_symbol",
  values = processed$results$Gene,
  mart = mart
) %>%
  dplyr::rename(Gene = hgnc_symbol)

annotated_df <- processed$results %>%
  left_join(gene_info, by = "Gene")

cat(sprintf("Saving results to S3: %s\n", s3_key_processed))
s3saveRDS(annotated_df, object = s3_key_processed, bucket = "scrna", region = "")

# Save to global environment
var_name <- paste0(tolower(celltype_group_input), "_kpmp_", 
                   gsub("_", "", tolower(analysis_type)))
assign(var_name, annotated_df, envir = .GlobalEnv)

cat(sprintf("Analysis complete for %s - %s at %s\n", 
            analysis_type, celltype_group_input, Sys.time()))

