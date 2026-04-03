# =============================================================================
# Compile all NEBULA run metadata CSVs from S3 into one dataframe
# =============================================================================
# Outputs: results/nebula/csv/run_metadata_compiled.csv (on S3)
# =============================================================================
library(aws.s3)
library(jsonlite)
library(dplyr)
# --- S3 setup ---
user <- Sys.info()[["user"]]
if (user == "choiyej") {
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "yejichoi") {
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else {
  stop("Unknown user: please specify keys path.")
}
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)
bucket <- "t1d.adiposity"
s3_base <- "results/nebula"
# --- All analysis config s3_subdirs (must match run_nebula_single.R) ---
analysis_subdirs <- c(
  # =========================================================================
  # CATEGORICAL (unadjusted)
  # =========================================================================
  # BMI within T1D
  "T1D_normal_vs_overweight_bmi",
  "T1D_nonobese_vs_obese_bmi",
  "T1D_normal_vs_obese_bmi",
  "T1D_normal_vs_ow_obese_bmi",
  # DXA within T1D
  "T1D_normal_vs_overweight_dxa",
  "T1D_nonobese_vs_obese_dxa",
  "T1D_normal_vs_obese_dxa",
  "T1D_normal_vs_ow_obese_dxa",
  # HC vs T1D (BMI)
  "HC_vs_T1D_normal_bmi",
  "HC_vs_T1D_overweight_bmi",
  "HC_vs_T1D_obese_bmi",
  # HC vs T1D (DXA)
  "HC_vs_T1D_normal_dxa",
  "HC_vs_T1D_overweight_dxa",
  "HC_vs_T1D_obese_dxa",

  # =========================================================================
  # CATEGORICAL (adj age)
  # =========================================================================
  # BMI within T1D
  "T1D_normal_vs_overweight_bmi_adj_age",
  "T1D_nonobese_vs_obese_bmi_adj_age",
  "T1D_normal_vs_obese_bmi_adj_age",
  "T1D_normal_vs_ow_obese_bmi_adj_age",
  # DXA within T1D
  "T1D_normal_vs_overweight_dxa_adj_age",
  "T1D_nonobese_vs_obese_dxa_adj_age",
  "T1D_normal_vs_obese_dxa_adj_age",
  "T1D_normal_vs_ow_obese_dxa_adj_age",
  # HC vs T1D (BMI)
  "HC_vs_T1D_normal_bmi_adj_age",
  "HC_vs_T1D_overweight_bmi_adj_age",
  "HC_vs_T1D_obese_bmi_adj_age",
  # HC vs T1D (DXA)
  "HC_vs_T1D_normal_dxa_adj_age",
  "HC_vs_T1D_overweight_dxa_adj_age",
  "HC_vs_T1D_obese_dxa_adj_age",

  # =========================================================================
  # CATEGORICAL (adj age + sex)
  # =========================================================================
  # BMI within T1D
  "T1D_normal_vs_overweight_bmi_adj_age_sex",
  "T1D_nonobese_vs_obese_bmi_adj_age_sex",
  "T1D_normal_vs_obese_bmi_adj_age_sex",
  "T1D_normal_vs_ow_obese_bmi_adj_age_sex",
  # DXA within T1D
  "T1D_normal_vs_overweight_dxa_adj_age_sex",
  "T1D_nonobese_vs_obese_dxa_adj_age_sex",
  "T1D_normal_vs_obese_dxa_adj_age_sex",
  "T1D_normal_vs_ow_obese_dxa_adj_age_sex",
  # HC vs T1D (BMI)
  "HC_vs_T1D_normal_bmi_adj_age_sex",
  "HC_vs_T1D_overweight_bmi_adj_age_sex",
  "HC_vs_T1D_obese_bmi_adj_age_sex",
  # HC vs T1D (DXA)
  "HC_vs_T1D_normal_dxa_adj_age_sex",
  "HC_vs_T1D_overweight_dxa_adj_age_sex",
  "HC_vs_T1D_obese_dxa_adj_age_sex",
  
  # =========================================================================
  # CONTINUOUS - T1D only (unadjusted)
  # =========================================================================
  "continuous/bmi_t1d",
  "continuous/dexa_body_fat_t1d",
  "continuous/dexa_bone_mineral_density_t1d",
  "continuous/dexa_fat_kg_t1d",
  "continuous/dexa_lean_mass_t1d",
  "continuous/dexa_lean_kg_t1d",
  "continuous/dexa_ag_ratio_t1d",
  "continuous/dexa_est_vat_t1d",
  "continuous/dexa_trunk_kg_t1d",
  "continuous/dexa_trunk_mass_t1d",
  
  # =========================================================================
  # CONTINUOUS - All subjects (unadjusted)
  # =========================================================================
  "continuous/bmi_all",
  "continuous/dexa_body_fat_all",
  "continuous/dexa_bone_mineral_density_all",
  "continuous/dexa_fat_kg_all",
  "continuous/dexa_lean_mass_all",
  "continuous/dexa_lean_kg_all",
  "continuous/dexa_ag_ratio_all",
  "continuous/dexa_est_vat_all",
  "continuous/dexa_trunk_kg_all",
  "continuous/dexa_trunk_mass_all",
  
  # =========================================================================
  # CONTINUOUS - All subjects (adj group)
  # =========================================================================
  "continuous/bmi_all_adj_group",
  "continuous/dexa_body_fat_all_adj_group",
  "continuous/dexa_bone_mineral_density_all_adj_group",
  "continuous/dexa_fat_kg_all_adj_group",
  "continuous/dexa_lean_mass_all_adj_group",
  "continuous/dexa_lean_kg_all_adj_group",
  "continuous/dexa_ag_ratio_all_adj_group",
  "continuous/dexa_est_vat_all_adj_group",
  "continuous/dexa_trunk_kg_all_adj_group",
  "continuous/dexa_trunk_mass_all_adj_group",
  
  # =========================================================================
  # CONTINUOUS - T1D only (adj age)
  # =========================================================================
  "continuous/bmi_t1d_adj_age",
  "continuous/dexa_body_fat_t1d_adj_age",
  "continuous/dexa_bone_mineral_density_t1d_adj_age",
  "continuous/dexa_fat_kg_t1d_adj_age",
  "continuous/dexa_lean_mass_t1d_adj_age",
  "continuous/dexa_lean_kg_t1d_adj_age",
  "continuous/dexa_ag_ratio_t1d_adj_age",
  "continuous/dexa_est_vat_t1d_adj_age",
  "continuous/dexa_trunk_kg_t1d_adj_age",
  "continuous/dexa_trunk_mass_t1d_adj_age",
  
  # =========================================================================
  # CONTINUOUS - All subjects (adj age)
  # =========================================================================
  "continuous/bmi_all_adj_age",
  "continuous/dexa_body_fat_all_adj_age",
  "continuous/dexa_bone_mineral_density_all_adj_age",
  "continuous/dexa_fat_kg_all_adj_age",
  "continuous/dexa_lean_mass_all_adj_age",
  "continuous/dexa_lean_kg_all_adj_age",
  "continuous/dexa_ag_ratio_all_adj_age",
  "continuous/dexa_est_vat_all_adj_age",
  "continuous/dexa_trunk_kg_all_adj_age",
  "continuous/dexa_trunk_mass_all_adj_age",
  
  # =========================================================================
  # CONTINUOUS - All subjects (adj group + age)
  # =========================================================================
  "continuous/bmi_all_adj_group_age",
  "continuous/dexa_body_fat_all_adj_group_age",
  "continuous/dexa_bone_mineral_density_all_adj_group_age",
  "continuous/dexa_fat_kg_all_adj_group_age",
  "continuous/dexa_lean_mass_all_adj_group_age",
  "continuous/dexa_lean_kg_all_adj_group_age",
  "continuous/dexa_ag_ratio_all_adj_group_age",
  "continuous/dexa_est_vat_all_adj_group_age",
  "continuous/dexa_trunk_kg_all_adj_group_age",
  "continuous/dexa_trunk_mass_all_adj_group_age",
  
  # =========================================================================
  # CONTINUOUS - T1D only (adj age + sex)
  # =========================================================================
  "continuous/bmi_t1d_adj_age_sex",
  "continuous/dexa_body_fat_t1d_adj_age_sex",
  "continuous/dexa_bone_mineral_density_t1d_adj_age_sex",
  "continuous/dexa_fat_kg_t1d_adj_age_sex",
  "continuous/dexa_lean_mass_t1d_adj_age_sex",
  "continuous/dexa_lean_kg_t1d_adj_age_sex",
  "continuous/dexa_ag_ratio_t1d_adj_age_sex",
  "continuous/dexa_est_vat_t1d_adj_age_sex",
  "continuous/dexa_trunk_kg_t1d_adj_age_sex",
  "continuous/dexa_trunk_mass_t1d_adj_age_sex",
  
  # =========================================================================
  # CONTINUOUS - All subjects (adj age + sex)
  # =========================================================================
  "continuous/bmi_all_adj_age_sex",
  "continuous/dexa_body_fat_all_adj_age_sex",
  "continuous/dexa_bone_mineral_density_all_adj_age_sex",
  "continuous/dexa_fat_kg_all_adj_age_sex",
  "continuous/dexa_lean_mass_all_adj_age_sex",
  "continuous/dexa_lean_kg_all_adj_age_sex",
  "continuous/dexa_ag_ratio_all_adj_age_sex",
  "continuous/dexa_est_vat_all_adj_age_sex",
  "continuous/dexa_trunk_kg_all_adj_age_sex",
  "continuous/dexa_trunk_mass_all_adj_age_sex",
  
  # =========================================================================
  # CONTINUOUS - All subjects (adj group + age + sex)
  # =========================================================================
  "continuous/bmi_all_adj_group_age_sex",
  "continuous/dexa_body_fat_all_adj_group_age_sex",
  "continuous/dexa_bone_mineral_density_all_adj_group_age_sex",
  "continuous/dexa_fat_kg_all_adj_group_age_sex",
  "continuous/dexa_lean_mass_all_adj_group_age_sex",
  "continuous/dexa_lean_kg_all_adj_group_age_sex",
  "continuous/dexa_ag_ratio_all_adj_group_age_sex",
  "continuous/dexa_est_vat_all_adj_group_age_sex",
  "continuous/dexa_trunk_kg_all_adj_group_age_sex",
  "continuous/dexa_trunk_mass_all_adj_group_age_sex",
  
  # =========================================================================
  # INTERACTION (continuous_var * group)
  # =========================================================================
  "interaction/bmi",
  "interaction/dexa_body_fat",
  "interaction/dexa_bone_mineral_density",
  "interaction/dexa_fat_kg",
  "interaction/dexa_lean_mass",
  "interaction/dexa_lean_kg",
  "interaction/dexa_ag_ratio",
  "interaction/dexa_est_vat",
  "interaction/dexa_trunk_kg",
  "interaction/dexa_trunk_mass"
)
# --- Find and compile all metadata CSVs ---
cat("Searching for metadata CSVs in S3...\n")
all_meta <- list()
subdirs_found <- c()
subdirs_missing <- c()
for (subdir in analysis_subdirs) {
  prefix <- paste0(s3_base, "/", subdir, "/")
  objs <- tryCatch(
    get_bucket(bucket = bucket, prefix = prefix, region = "", max = 10000),
    error = function(e) { message(paste("Error listing", prefix, ":", e$message)); list() }
  )
  meta_keys <- sapply(objs, function(x) x$Key)
  meta_keys <- meta_keys[grepl("_metadata\\.csv$", meta_keys)]
  if (length(meta_keys) > 0) {
    subdirs_found <- c(subdirs_found, subdir)
  } else {
    subdirs_missing <- c(subdirs_missing, subdir)
  }
  for (key in meta_keys) {
    tryCatch({
      tmp <- tempfile(fileext = ".csv")
      save_object(object = key, bucket = bucket, region = "", file = tmp)
      df <- read.csv(tmp, stringsAsFactors = FALSE)
      all_meta[[length(all_meta) + 1]] <- df
      unlink(tmp)
    }, error = function(e) {
      message(paste("Error reading", key, ":", e$message))
    })
  }
}
cat(sprintf("\n=== S3 scan summary ===\n"))
cat(sprintf("  Analysis types with results: %d / %d\n", length(subdirs_found), length(analysis_subdirs)))
if (length(subdirs_missing) > 0) {
  cat(sprintf("  Analysis types still pending: %d\n", length(subdirs_missing)))
  cat(paste("    -", subdirs_missing, collapse = "\n"))
  cat("\n")
}
cat(sprintf("Found %d metadata files\n", length(all_meta)))
if (length(all_meta) > 0) {
  meta_df_compiled <- bind_rows(all_meta)
  # Save compiled metadata to S3
  tmp_out <- tempfile(fileext = ".csv")
  write.csv(meta_df_compiled, tmp_out, row.names = FALSE)
  put_object(file = tmp_out, object = paste0(s3_base, "/csv/run_metadata_compiled.csv"),
             bucket = bucket, region = "")
  unlink(tmp_out)
  cat(sprintf("Compiled metadata saved: %d rows x %d columns\n",
              nrow(meta_df_compiled), ncol(meta_df_compiled)))
  cat(sprintf("  Unique analysis types: %d\n", n_distinct(meta_df_compiled$analysis_type)))
  cat(sprintf("  Unique cell types: %d\n", n_distinct(meta_df_compiled$celltype)))
} else {
  cat("No metadata files found. Ensure NEBULA jobs have completed.\n")
}