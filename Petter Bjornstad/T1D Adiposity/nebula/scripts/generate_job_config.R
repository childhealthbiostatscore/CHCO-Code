#!/usr/bin/env Rscript
# =============================================================================
# Generate job_config.txt for T1D Adiposity NEBULA array jobs
# =============================================================================
# Single unified config replacing:
#   - generate_job_config.R          (this file, now consolidated)
#   - generate_job_config_adj_age_sex.R
#   - generate_job_config_normal_vs_ow_obese.R
#   - generate_job_config_discordance.R
#
# Run locally to generate the config file before submitting SLURM jobs.
# Output: config/job_config.txt (tab-delimited:
#   analysis_type  celltype  celltype_var)
#
# For the initial full-rerun pass, we restrict cell types to the lineages
# of interest: PT, TAL, PC, IC, EC, Immune (and their KPMP subtypes).
# Unused cell types are commented out (not deleted) for easy re-enabling.
#
# Study-adjusted (STEP 3) analyses and matched-subject (STEP 6) analyses
# are intentionally excluded from this pass. If you want to run them, also
# uncomment the corresponding entries in run_nebula_single.R.
# =============================================================================

# -----------------------------------------------------------------------------
# ANALYSIS TYPES
# -----------------------------------------------------------------------------
# All must match names in run_nebula_single.R's analysis_config.
analysis_types <- c(
  # =========================================================================
  # CATEGORICAL (unadjusted) - within-T1D
  # =========================================================================
  "T1D_normal_vs_overweight_bmi",
  "T1D_nonobese_vs_obese_bmi",
  "T1D_normal_vs_obese_bmi",
  "T1D_normal_vs_ow_obese_bmi",
  "T1D_normal_vs_overweight_dxa",
  "T1D_nonobese_vs_obese_dxa",
  "T1D_normal_vs_obese_dxa",
  "T1D_normal_vs_ow_obese_dxa",

  # =========================================================================
  # CATEGORICAL (unadjusted) - HC vs T1D
  # =========================================================================
  "HC_vs_T1D_normal_bmi",
  "HC_vs_T1D_overweight_bmi",
  "HC_vs_T1D_obese_bmi",
  "HC_vs_T1D_normal_dxa",
  "HC_vs_T1D_overweight_dxa",
  "HC_vs_T1D_obese_dxa",

  # =========================================================================
  # CATEGORICAL - age-adjusted (within-T1D and HC vs T1D, BMI + DXA)
  # =========================================================================
  "T1D_normal_vs_overweight_bmi_adj_age",
  "T1D_nonobese_vs_obese_bmi_adj_age",
  "T1D_normal_vs_obese_bmi_adj_age",
  "T1D_normal_vs_ow_obese_bmi_adj_age",
  "T1D_normal_vs_overweight_dxa_adj_age",
  "T1D_nonobese_vs_obese_dxa_adj_age",
  "T1D_normal_vs_obese_dxa_adj_age",
  "T1D_normal_vs_ow_obese_dxa_adj_age",
  "HC_vs_T1D_normal_bmi_adj_age",
  "HC_vs_T1D_overweight_bmi_adj_age",
  "HC_vs_T1D_obese_bmi_adj_age",
  "HC_vs_T1D_normal_dxa_adj_age",
  "HC_vs_T1D_overweight_dxa_adj_age",
  "HC_vs_T1D_obese_dxa_adj_age",

  # =========================================================================
  # CATEGORICAL - age+sex-adjusted (merged from _adj_age_sex config)
  # =========================================================================
  "T1D_normal_vs_overweight_bmi_adj_age_sex",
  "T1D_nonobese_vs_obese_bmi_adj_age_sex",
  "T1D_normal_vs_obese_bmi_adj_age_sex",
  "T1D_normal_vs_ow_obese_bmi_adj_age_sex",
  "T1D_normal_vs_overweight_dxa_adj_age_sex",
  "T1D_nonobese_vs_obese_dxa_adj_age_sex",
  "T1D_normal_vs_obese_dxa_adj_age_sex",
  "T1D_normal_vs_ow_obese_dxa_adj_age_sex",
  "HC_vs_T1D_normal_bmi_adj_age_sex",
  "HC_vs_T1D_overweight_bmi_adj_age_sex",
  "HC_vs_T1D_obese_bmi_adj_age_sex",
  "HC_vs_T1D_normal_dxa_adj_age_sex",
  "HC_vs_T1D_overweight_dxa_adj_age_sex",
  "HC_vs_T1D_obese_dxa_adj_age_sex",

  # =========================================================================
  # STEP 2 DIAGNOSTIC - BMI with ATTEMPT excluded (merged from _discordance)
  # =========================================================================
  "T1D_normal_vs_ow_obese_bmi_noattempt",
  "T1D_normal_vs_ow_obese_bmi_noattempt_adj_age",
  "T1D_normal_vs_ow_obese_bmi_noattempt_adj_age_sex",
  "T1D_nonobese_vs_obese_bmi_noattempt",

  # =========================================================================
  # STEP 3 DIAGNOSTIC - DXA adjusted for study (NOT run this pass)
  # Uncomment these AND the corresponding block in run_nebula_single.R.
  # =========================================================================
  # "T1D_normal_vs_ow_obese_dxa_adj_study",
  # "T1D_normal_vs_ow_obese_dxa_adj_age_study",
  # "T1D_normal_vs_ow_obese_dxa_adj_age_sex_study",
  # "T1D_nonobese_vs_obese_dxa_adj_age_study",
  # "cont_dexa_body_fat_t1d_adj_age_study",
  # "cont_dexa_est_vat_t1d_adj_age_study",
  # "cont_dexa_ag_ratio_t1d_adj_age_study",
  # "cont_dexa_trunk_kg_t1d_adj_age_study",
  # "cont_dexa_lean_mass_t1d_adj_age_study",

  # =========================================================================
  # STEP 6 DIAGNOSTIC - matched-subject (NOT run this pass)
  # Uncomment these AND the corresponding block in run_nebula_single.R.
  # =========================================================================
  # "T1D_normal_vs_ow_obese_bmi_matched",
  # "T1D_normal_vs_ow_obese_dxa_matched",
  # "T1D_normal_vs_ow_obese_bmi_matched_adj_age_study",
  # "T1D_normal_vs_ow_obese_dxa_matched_adj_age_study",

  # =========================================================================
  # CONTINUOUS - T1D only (unadjusted)
  # =========================================================================
  "cont_bmi_t1d",
  "cont_dexa_body_fat_t1d",
  "cont_dexa_bone_mineral_density_t1d",
  "cont_dexa_fat_kg_t1d",
  "cont_dexa_lean_mass_t1d",
  "cont_dexa_lean_kg_t1d",
  "cont_dexa_ag_ratio_t1d",
  "cont_dexa_est_vat_t1d",
  "cont_dexa_trunk_kg_t1d",
  "cont_dexa_trunk_mass_t1d",

  # =========================================================================
  # CONTINUOUS - All subjects, unadjusted
  # =========================================================================
  "cont_bmi_all",
  "cont_dexa_body_fat_all",
  "cont_dexa_bone_mineral_density_all",
  "cont_dexa_fat_kg_all",
  "cont_dexa_lean_mass_all",
  "cont_dexa_lean_kg_all",
  "cont_dexa_ag_ratio_all",
  "cont_dexa_est_vat_all",
  "cont_dexa_trunk_kg_all",
  "cont_dexa_trunk_mass_all",

  # =========================================================================
  # CONTINUOUS - All subjects, adjusted for group
  # =========================================================================
  "cont_bmi_all_adj",
  "cont_dexa_body_fat_all_adj",
  "cont_dexa_bone_mineral_density_all_adj",
  "cont_dexa_fat_kg_all_adj",
  "cont_dexa_lean_mass_all_adj",
  "cont_dexa_lean_kg_all_adj",
  "cont_dexa_ag_ratio_all_adj",
  "cont_dexa_est_vat_all_adj",
  "cont_dexa_trunk_kg_all_adj",
  "cont_dexa_trunk_mass_all_adj",

  # =========================================================================
  # CONTINUOUS - T1D only, adjusted for age
  # =========================================================================
  "cont_bmi_t1d_adj_age",
  "cont_dexa_body_fat_t1d_adj_age",
  "cont_dexa_bone_mineral_density_t1d_adj_age",
  "cont_dexa_fat_kg_t1d_adj_age",
  "cont_dexa_lean_mass_t1d_adj_age",
  "cont_dexa_lean_kg_t1d_adj_age",
  "cont_dexa_ag_ratio_t1d_adj_age",
  "cont_dexa_est_vat_t1d_adj_age",
  "cont_dexa_trunk_kg_t1d_adj_age",
  "cont_dexa_trunk_mass_t1d_adj_age",

  # =========================================================================
  # CONTINUOUS - All subjects, adjusted for age
  # =========================================================================
  "cont_bmi_all_adj_age",
  "cont_dexa_body_fat_all_adj_age",
  "cont_dexa_bone_mineral_density_all_adj_age",
  "cont_dexa_fat_kg_all_adj_age",
  "cont_dexa_lean_mass_all_adj_age",
  "cont_dexa_lean_kg_all_adj_age",
  "cont_dexa_ag_ratio_all_adj_age",
  "cont_dexa_est_vat_all_adj_age",
  "cont_dexa_trunk_kg_all_adj_age",
  "cont_dexa_trunk_mass_all_adj_age",

  # =========================================================================
  # CONTINUOUS - All subjects, adjusted for group + age
  # =========================================================================
  "cont_bmi_all_adj_group_age",
  "cont_dexa_body_fat_all_adj_group_age",
  "cont_dexa_bone_mineral_density_all_adj_group_age",
  "cont_dexa_fat_kg_all_adj_group_age",
  "cont_dexa_lean_mass_all_adj_group_age",
  "cont_dexa_lean_kg_all_adj_group_age",
  "cont_dexa_ag_ratio_all_adj_group_age",
  "cont_dexa_est_vat_all_adj_group_age",
  "cont_dexa_trunk_kg_all_adj_group_age",
  "cont_dexa_trunk_mass_all_adj_group_age",

  # =========================================================================
  # CONTINUOUS - T1D only, adjusted for age + sex (merged from _adj_age_sex)
  # =========================================================================
  "cont_bmi_t1d_adj_age_sex",
  "cont_dexa_body_fat_t1d_adj_age_sex",
  "cont_dexa_bone_mineral_density_t1d_adj_age_sex",
  "cont_dexa_fat_kg_t1d_adj_age_sex",
  "cont_dexa_lean_mass_t1d_adj_age_sex",
  "cont_dexa_lean_kg_t1d_adj_age_sex",
  "cont_dexa_ag_ratio_t1d_adj_age_sex",
  "cont_dexa_est_vat_t1d_adj_age_sex",
  "cont_dexa_trunk_kg_t1d_adj_age_sex",
  "cont_dexa_trunk_mass_t1d_adj_age_sex",

  # =========================================================================
  # CONTINUOUS - All subjects, adjusted for age + sex (merged from _adj_age_sex)
  # =========================================================================
  "cont_bmi_all_adj_age_sex",
  "cont_dexa_body_fat_all_adj_age_sex",
  "cont_dexa_bone_mineral_density_all_adj_age_sex",
  "cont_dexa_fat_kg_all_adj_age_sex",
  "cont_dexa_lean_mass_all_adj_age_sex",
  "cont_dexa_lean_kg_all_adj_age_sex",
  "cont_dexa_ag_ratio_all_adj_age_sex",
  "cont_dexa_est_vat_all_adj_age_sex",
  "cont_dexa_trunk_kg_all_adj_age_sex",
  "cont_dexa_trunk_mass_all_adj_age_sex",

  # =========================================================================
  # CONTINUOUS - All subjects, adjusted for group + age + sex (merged)
  # =========================================================================
  "cont_bmi_all_adj_group_age_sex",
  "cont_dexa_body_fat_all_adj_group_age_sex",
  "cont_dexa_bone_mineral_density_all_adj_group_age_sex",
  "cont_dexa_fat_kg_all_adj_group_age_sex",
  "cont_dexa_lean_mass_all_adj_group_age_sex",
  "cont_dexa_lean_kg_all_adj_group_age_sex",
  "cont_dexa_ag_ratio_all_adj_group_age_sex",
  "cont_dexa_est_vat_all_adj_group_age_sex",
  "cont_dexa_trunk_kg_all_adj_group_age_sex",
  "cont_dexa_trunk_mass_all_adj_group_age_sex",

  # =========================================================================
  # INTERACTION - continuous_var * group (HC vs T1D)
  # =========================================================================
  "interaction_bmi",
  "interaction_dexa_body_fat",
  "interaction_dexa_bone_mineral_density",
  "interaction_dexa_fat_kg",
  "interaction_dexa_lean_mass",
  "interaction_dexa_lean_kg",
  "interaction_dexa_ag_ratio",
  "interaction_dexa_est_vat",
  "interaction_dexa_trunk_kg",
  "interaction_dexa_trunk_mass"
)

# -----------------------------------------------------------------------------
# CELL TYPES
# -----------------------------------------------------------------------------
# Focus on the main lineages of interest for this initial pass:
#   PT, TAL, PC, IC, EC, Immune (+ their KPMP subtypes).
# Unused cell types are commented out (not deleted) — re-enable by
# uncommenting the relevant lines.

# KPMP_celltype (specific subtypes within the lineages of interest)
kpmp_celltypes <- c(
  # PT subtypes
  "PT-S1/S2", "PT-S3", "aPT",
  # TAL subtypes
  "C-TAL-1", "C-TAL-2", "aTAL", "dTAL",
  # PC subtypes
  "CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC",
  # IC subtypes
  "IC-A", "IC-B", "aIC",
  # EC subtypes
  "EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A",
  # Immune subtypes (myeloid + lymphoid)
  "MAC", "MON", "cDC", "pDC",
  "CD4+ T", "CD8+ T", "B", "NK", "cycT"

  # ---- Cell types NOT included in this pass (uncomment to re-enable) ----
  # , "DTL", "aDTL", "ATL"        # thin limbs
  # , "DCT", "dDCT", "CNT"        # distal convoluted / connecting tubule
  # , "VSMC/P", "FIB"             # vascular smooth muscle / fibroblast
  # , "POD", "MC", "PEC"          # glomerular (podocytes, mesangial, PEC)
  # , "SchwannCells", "non-specific"
)

# KPMP_celltype_general (grouped lineages)
kpmp_general <- c(
  "PT", "TAL", "PC", "IC", "EC", "Immune"

  # ---- Grouped lineages NOT included in this pass ----
  # , "DTL_ATL", "DCT_CNT"
  # , "VSMC_P_FIB"
  # , "POD", "MC", "PEC"
  # , "Schwann", "Other"
)

# -----------------------------------------------------------------------------
# BUILD CONFIG
# -----------------------------------------------------------------------------
config_lines <- data.frame(
  analysis_type = character(0),
  celltype = character(0),
  celltype_var = character(0),
  stringsAsFactors = FALSE
)

for (analysis in analysis_types) {
  for (ct in kpmp_celltypes) {
    config_lines <- rbind(config_lines, data.frame(
      analysis_type = analysis,
      celltype = ct,
      celltype_var = "KPMP_celltype",
      stringsAsFactors = FALSE
    ))
  }
  for (ct in kpmp_general) {
    config_lines <- rbind(config_lines, data.frame(
      analysis_type = analysis,
      celltype = ct,
      celltype_var = "KPMP_celltype_general",
      stringsAsFactors = FALSE
    ))
  }
}

# -----------------------------------------------------------------------------
# WRITE
# -----------------------------------------------------------------------------
config_dir <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/T1D Adiposity/nebula/config"
if (!dir.exists(config_dir)) dir.create(config_dir, recursive = TRUE)

config_file <- file.path(config_dir, "job_config.txt")
write.table(config_lines, file = config_file, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat(sprintf("Generated %d jobs in %s\n", nrow(config_lines), config_file))
cat(sprintf("  - %d analysis types\n", length(analysis_types)))
cat(sprintf("  - %d KPMP_celltype subtypes\n", length(kpmp_celltypes)))
cat(sprintf("  - %d KPMP_celltype_general lineages\n", length(kpmp_general)))
cat(sprintf("  - Total jobs = %d analyses x (%d + %d) cell types = %d\n",
            length(analysis_types), length(kpmp_celltypes), length(kpmp_general),
            nrow(config_lines)))

cat("\nAnalysis breakdown:\n")
cat(sprintf("  Categorical, unadjusted:              %d\n",
            sum(!grepl("adj|_noattempt|^cont_|^interaction_", analysis_types))))
cat(sprintf("  Categorical, adj_age:                 %d\n",
            sum(grepl("_adj_age$", analysis_types) & !grepl("^cont_", analysis_types))))
cat(sprintf("  Categorical, adj_age_sex:             %d\n",
            sum(grepl("_adj_age_sex$", analysis_types) & !grepl("^cont_", analysis_types))))
cat(sprintf("  STEP 2 (BMI no-ATTEMPT):              %d\n",
            sum(grepl("_noattempt", analysis_types))))
cat(sprintf("  Continuous associations:              %d\n",
            sum(grepl("^cont_", analysis_types))))
cat(sprintf("  Interaction models:                   %d\n",
            sum(grepl("^interaction_", analysis_types))))

cat(sprintf("\nTo submit: update --array=1-%d in nebula_array.slurm\n",
            nrow(config_lines)))
cat(sprintf("JOB_CONFIG=%s\n", config_file))
