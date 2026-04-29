#!/usr/bin/env Rscript
# =============================================================================
# Generate job_config_new_runs.txt — DELTA JOBS ONLY
# =============================================================================
# This config is a TARGETED submission for analyses added after the original
# sweep was already run. It contains ONLY:
#
#   1. All WHtR-based contrasts:
#        - within-T1D categorical (4 contrasts × 3 adjustment families)
#        - HC vs T1D per category (3 categories × 3 adjustments)
#        - HC vs T1D Overweight+Obese combined (× 3 adjustments)
#        - continuous WHtR (T1D-only and all-subjects, ± group, ± age, ± sex)
#
#   2. NEW HC vs T1D Overweight+Obese contrasts for BMI and DXA (these
#      previously existed only as separate Normal/Overweight/Obese vs HC,
#      not as a combined OW+Ob vs HC contrast).
#
# Use generate_job_config.R for the FULL master config (185 analyses × 37
# cell types). Use THIS script when you only need to submit the delta jobs.
#
# Output: config/job_config_new_runs.txt
#         (tab-delimited: analysis_type  celltype  celltype_var)
# =============================================================================

# -----------------------------------------------------------------------------
# ANALYSIS TYPES — new additions only
# -----------------------------------------------------------------------------
analysis_types <- c(
  # =========================================================================
  # WHtR categorical (within T1D) — unadjusted, age, age+sex
  # =========================================================================
  "T1D_normal_vs_overweight_whtr",
  "T1D_nonobese_vs_obese_whtr",
  "T1D_normal_vs_obese_whtr",
  "T1D_normal_vs_ow_obese_whtr",

  "T1D_normal_vs_overweight_whtr_adj_age",
  "T1D_nonobese_vs_obese_whtr_adj_age",
  "T1D_normal_vs_obese_whtr_adj_age",
  "T1D_normal_vs_ow_obese_whtr_adj_age",

  "T1D_normal_vs_overweight_whtr_adj_age_sex",
  "T1D_nonobese_vs_obese_whtr_adj_age_sex",
  "T1D_normal_vs_obese_whtr_adj_age_sex",
  "T1D_normal_vs_ow_obese_whtr_adj_age_sex",

  # =========================================================================
  # WHtR HC vs T1D (per-category) — unadjusted, age, age+sex
  # =========================================================================
  "HC_vs_T1D_normal_whtr",
  "HC_vs_T1D_overweight_whtr",
  "HC_vs_T1D_obese_whtr",

  "HC_vs_T1D_normal_whtr_adj_age",
  "HC_vs_T1D_overweight_whtr_adj_age",
  "HC_vs_T1D_obese_whtr_adj_age",

  "HC_vs_T1D_normal_whtr_adj_age_sex",
  "HC_vs_T1D_overweight_whtr_adj_age_sex",
  "HC_vs_T1D_obese_whtr_adj_age_sex",

  # =========================================================================
  # HC vs T1D Overweight+Obese (combined) — BMI, DXA, WHtR
  # All three are NEW; previously only Normal/Overweight/Obese individually
  # existed against HC.
  # =========================================================================
  "HC_vs_T1D_ow_obese_bmi",
  "HC_vs_T1D_ow_obese_dxa",
  "HC_vs_T1D_ow_obese_whtr",

  "HC_vs_T1D_ow_obese_bmi_adj_age",
  "HC_vs_T1D_ow_obese_dxa_adj_age",
  "HC_vs_T1D_ow_obese_whtr_adj_age",

  "HC_vs_T1D_ow_obese_bmi_adj_age_sex",
  "HC_vs_T1D_ow_obese_dxa_adj_age_sex",
  "HC_vs_T1D_ow_obese_whtr_adj_age_sex",

  # =========================================================================
  # WHtR continuous — T1D-only and all-subjects, ± group, ± age, ± sex
  # =========================================================================
  "cont_whtr_t1d",
  "cont_whtr_all",
  "cont_whtr_all_adj",                # adjusted for group only
  "cont_whtr_t1d_adj_age",
  "cont_whtr_all_adj_age",
  "cont_whtr_all_adj_group_age",
  "cont_whtr_t1d_adj_age_sex",
  "cont_whtr_all_adj_age_sex",
  "cont_whtr_all_adj_group_age_sex"
)

# -----------------------------------------------------------------------------
# CELL TYPES — same focus list as the master config
# -----------------------------------------------------------------------------
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
  # Immune subtypes
  "MAC", "MON", "cDC", "pDC",
  "CD4+ T", "CD8+ T", "B", "NK", "cycT"

  # ---- Cell types NOT included in this pass (uncomment to re-enable) ----
  # , "DTL", "aDTL", "ATL"
  # , "DCT", "dDCT", "CNT"
  # , "VSMC/P", "FIB"
  # , "POD", "MC", "PEC"
  # , "SchwannCells", "non-specific"
)

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

config_file <- file.path(config_dir, "job_config_new_runs.txt")
write.table(config_lines, file = config_file, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# -----------------------------------------------------------------------------
# REPORT
# -----------------------------------------------------------------------------
cat(sprintf("Generated %d delta jobs in %s\n", nrow(config_lines), config_file))
cat(sprintf("  - %d analysis types (NEW additions only)\n", length(analysis_types)))
cat(sprintf("  - %d KPMP_celltype subtypes\n", length(kpmp_celltypes)))
cat(sprintf("  - %d KPMP_celltype_general lineages\n", length(kpmp_general)))
cat(sprintf("  - Total jobs = %d analyses x (%d + %d) cell types = %d\n",
            length(analysis_types), length(kpmp_celltypes), length(kpmp_general),
            nrow(config_lines)))

cat("\nBreakdown by analysis family:\n")
cat(sprintf("  WHtR categorical within-T1D:           %d\n",
            sum(grepl("^T1D.*whtr", analysis_types))))
cat(sprintf("  WHtR HC vs T1D (per-category):         %d\n",
            sum(grepl("^HC_vs_T1D_(normal|overweight|obese)_whtr", analysis_types))))
cat(sprintf("  HC vs T1D OW+Ob (BMI / DXA / WHtR):    %d\n",
            sum(grepl("^HC_vs_T1D_ow_obese_", analysis_types))))
cat(sprintf("  WHtR continuous:                       %d\n",
            sum(grepl("^cont_whtr_", analysis_types))))

cat(sprintf("\nTo submit:\n"))
cat(sprintf("  1. Update --array=1-%d in your sbatch script\n", nrow(config_lines)))
cat(sprintf("  2. Set JOB_CONFIG=%s\n", config_file))
cat(sprintf("  3. TASK_SCRIPT remains nebula/scripts/run_nebula_single.R\n"))
cat(sprintf("  (run_nebula_single.R already contains all listed analysis_config entries)\n"))
