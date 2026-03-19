#!/usr/bin/env Rscript
# =============================================================================
# Generate job_config.txt for T1D Adiposity NEBULA array jobs
# =============================================================================
# Run this locally to generate the config file before submitting SLURM jobs.
# Output: config/job_config.txt (tab-delimited: analysis_type celltype celltype_var)
# =============================================================================

# All analysis types (must match names in run_nebula_single.R)
analysis_types <- c(
  # ---- Categorical: BMI-defined within T1D (unadjusted) ----
  "T1D_normal_vs_overweight_bmi",
  "T1D_nonobese_vs_obese_bmi",
  "T1D_normal_vs_obese_bmi",
  # ---- Categorical: DXA-defined within T1D (unadjusted) ----
  "T1D_normal_vs_overweight_dxa",
  "T1D_nonobese_vs_obese_dxa",
  "T1D_normal_vs_obese_dxa",
  # ---- Categorical: HC vs T1D (BMI-defined, unadjusted) ----
  "HC_vs_T1D_normal_bmi",
  "HC_vs_T1D_overweight_bmi",
  "HC_vs_T1D_obese_bmi",
  # ---- Categorical: HC vs T1D (DXA-defined, unadjusted) ----
  "HC_vs_T1D_normal_dxa",
  "HC_vs_T1D_overweight_dxa",
  "HC_vs_T1D_obese_dxa",

  # ---- Categorical: BMI-defined within T1D (adjusted for age) ----
  "T1D_normal_vs_overweight_bmi_adj_age",
  "T1D_nonobese_vs_obese_bmi_adj_age",
  "T1D_normal_vs_obese_bmi_adj_age",
  # ---- Categorical: DXA-defined within T1D (adjusted for age) ----
  "T1D_normal_vs_overweight_dxa_adj_age",
  "T1D_nonobese_vs_obese_dxa_adj_age",
  "T1D_normal_vs_obese_dxa_adj_age",
  # ---- Categorical: HC vs T1D (BMI-defined, adjusted for age) ----
  "HC_vs_T1D_normal_bmi_adj_age",
  "HC_vs_T1D_overweight_bmi_adj_age",
  "HC_vs_T1D_obese_bmi_adj_age",
  # ---- Categorical: HC vs T1D (DXA-defined, adjusted for age) ----
  "HC_vs_T1D_normal_dxa_adj_age",
  "HC_vs_T1D_overweight_dxa_adj_age",
  "HC_vs_T1D_obese_dxa_adj_age",

  # ---- Continuous: T1D only ----
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
  # ---- Continuous: All subjects, unadjusted ----
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
  # ---- Continuous: All subjects, adjusted for group ----
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

  # ---- Continuous: T1D only, adjusted for age ----
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
  # ---- Continuous: All subjects, adjusted for age ----
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
  # ---- Continuous: All subjects, adjusted for group + age ----
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

  # ---- Interaction: continuous_var * group ----
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

# KPMP_celltype (specific cell types)
kpmp_celltypes <- c(
  "PT-S1/S2", "PT-S3", "aPT",
  "C-TAL-1", "C-TAL-2", "aTAL", "dTAL",
  "CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC",
  "IC-A", "IC-B", "aIC",
  "DTL", "aDTL", "ATL",
  "DCT", "dDCT", "CNT",
  "EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A",
  "MAC", "MON", "cDC", "pDC",
  "CD4+ T", "CD8+ T", "B", "NK", "cycT",
  "VSMC/P", "FIB",
  "POD", "MC", "PEC",
  "SchwannCells", "non-specific"
)

# KPMP_celltype_general (grouped cell types)
kpmp_general <- c(
  "PT", "TAL", "PC", "IC", "DTL_ATL", "DCT_CNT", "EC",
  "Immune_Myeloid", "Immune_Lymphoid",
  "VSMC_P_FIB", "POD", "MC", "PEC", "Schwann", "Other"
)

# Build the config data frame
config_lines <- data.frame(
  analysis_type = character(0),
  celltype = character(0),
  celltype_var = character(0),
  stringsAsFactors = FALSE
)

for (analysis in analysis_types) {
  # Add KPMP_celltype entries
  for (ct in kpmp_celltypes) {
    config_lines <- rbind(config_lines, data.frame(
      analysis_type = analysis,
      celltype = ct,
      celltype_var = "KPMP_celltype",
      stringsAsFactors = FALSE
    ))
  }
  # Add KPMP_celltype_general entries
  for (ct in kpmp_general) {
    config_lines <- rbind(config_lines, data.frame(
      analysis_type = analysis,
      celltype = ct,
      celltype_var = "KPMP_celltype_general",
      stringsAsFactors = FALSE
    ))
  }
}

# Write config file
config_dir <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/T1D Adiposity/nebula/config"
if (!dir.exists(config_dir)) dir.create(config_dir, recursive = TRUE)

config_file <- file.path(config_dir, "job_config.txt")
write.table(config_lines, file = config_file, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat(sprintf("Generated %d jobs in %s\n", nrow(config_lines), config_file))
cat(sprintf("  - %d analysis types\n", length(analysis_types)))
cat(sprintf("  - %d KPMP_celltype levels\n", length(kpmp_celltypes)))
cat(sprintf("  - %d KPMP_celltype_general levels\n", length(kpmp_general)))
cat(sprintf("  - Total: %d analyses x (%d + %d) cell types = %d jobs\n",
            length(analysis_types), length(kpmp_celltypes), length(kpmp_general),
            nrow(config_lines)))

# Print summary
cat("\nAnalysis breakdown:\n")
cat(sprintf("  Categorical comparisons (unadjusted): 12\n"))
cat(sprintf("  Categorical comparisons (adjusted for age): 12\n"))
cat(sprintf("  Continuous T1D only (unadjusted): 10\n"))
cat(sprintf("  Continuous all subjects (unadjusted): 10\n"))
cat(sprintf("  Continuous all subjects (adjusted for group): 10\n"))
cat(sprintf("  Continuous T1D only (adjusted for age): 10\n"))
cat(sprintf("  Continuous all subjects (adjusted for age): 10\n"))
cat(sprintf("  Continuous all subjects (adjusted for group + age): 10\n"))
cat(sprintf("  Interaction (continuous_var * group): 10\n"))
cat(sprintf("  Total unique analyses: %d\n", length(analysis_types)))

cat(sprintf("\nTo submit: update --array=1-%d in nebula_array.slurm\n", nrow(config_lines)))
