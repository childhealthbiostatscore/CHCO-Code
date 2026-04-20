#!/usr/bin/env Rscript
# =============================================================================
# Generate job_config_discordance.txt for T1D Adiposity NEBULA array jobs
# Diagnostic analyses for the BMI-vs-DXA discordance investigation:
#   Step 2: BMI comparisons with ATTEMPT excluded (_noattempt)
#   Step 3: DXA comparisons adjusted for study as a fixed effect (_adj_*_study)
# =============================================================================
# Run this locally to generate the config file before submitting SLURM jobs.
# Output: config/job_config_discordance.txt (tab-delimited:
#   analysis_type  celltype  celltype_var)
# =============================================================================

# Analysis types added in run_nebula_single.R (see "STEP 2 DIAGNOSTIC" and
# "STEP 3 DIAGNOSTIC" sections in analysis_config).
analysis_types <- c(
  # ---- STEP 2: BMI with ATTEMPT excluded ----
  "T1D_normal_vs_ow_obese_bmi_noattempt",
  "T1D_normal_vs_ow_obese_bmi_noattempt_adj_age",
  "T1D_normal_vs_ow_obese_bmi_noattempt_adj_age_sex",
  "T1D_nonobese_vs_obese_bmi_noattempt",

  # ---- STEP 3: DXA adjusted for study ----
  "T1D_normal_vs_ow_obese_dxa_adj_study",
  "T1D_normal_vs_ow_obese_dxa_adj_age_study",
  "T1D_normal_vs_ow_obese_dxa_adj_age_sex_study",
  "T1D_nonobese_vs_obese_dxa_adj_age_study",

  # ---- STEP 3: continuous DXA adjusted for age + study ----
  "cont_dexa_body_fat_t1d_adj_age_study",
  "cont_dexa_est_vat_t1d_adj_age_study",
  "cont_dexa_ag_ratio_t1d_adj_age_study",
  "cont_dexa_trunk_kg_t1d_adj_age_study",
  "cont_dexa_lean_mass_t1d_adj_age_study",

  # ---- STEP 6: matched-subject (has_both_bmi_dxa == TRUE) ----
  "T1D_normal_vs_ow_obese_bmi_matched",
  "T1D_normal_vs_ow_obese_dxa_matched",
  "T1D_normal_vs_ow_obese_bmi_matched_adj_age_study",
  "T1D_normal_vs_ow_obese_dxa_matched_adj_age_study"
)

# Focus on the cell-type groupings most relevant to the kidney question.
# The full KPMP_celltype list is kept for thoroughness; trim if you want a
# smaller job batch.
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
  "POD", "MC", "PEC"
)

kpmp_general <- c(
  "PT", "TAL", "PC", "IC", "DTL_ATL", "DCT_CNT", "EC",
  "Immune_Myeloid", "Immune_Lymphoid",
  "VSMC_P_FIB", "POD", "MC", "PEC"
)

# Build the config data frame
config_lines <- data.frame(
  analysis_type = character(0),
  celltype = character(0),
  celltype_var = character(0),
  stringsAsFactors = FALSE
)

for (analysis in analysis_types) {
  for (ct in kpmp_celltypes) {
    config_lines <- rbind(config_lines, data.frame(
      analysis_type = analysis, celltype = ct,
      celltype_var = "KPMP_celltype", stringsAsFactors = FALSE))
  }
  for (ct in kpmp_general) {
    config_lines <- rbind(config_lines, data.frame(
      analysis_type = analysis, celltype = ct,
      celltype_var = "KPMP_celltype_general", stringsAsFactors = FALSE))
  }
}

config_dir <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad/T1D Adiposity/nebula/config"
if (!dir.exists(config_dir)) dir.create(config_dir, recursive = TRUE)

config_file <- file.path(config_dir, "job_config_discordance.txt")
write.table(config_lines, file = config_file, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat(sprintf("Generated %d jobs in %s\n", nrow(config_lines), config_file))
cat(sprintf("  - %d analysis types\n", length(analysis_types)))
cat(sprintf("  - %d KPMP_celltype levels\n", length(kpmp_celltypes)))
cat(sprintf("  - %d KPMP_celltype_general levels\n", length(kpmp_general)))
cat(sprintf("\nTo submit: copy nebula_array.slurm to nebula_array_discordance.slurm\n"))
cat(sprintf("and set JOB_CONFIG=%s and --array=1-%d\n",
            config_file, nrow(config_lines)))
