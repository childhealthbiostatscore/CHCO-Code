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
  "interaction/dexa_trunk_mass",

  # =========================================================================
  # DISCORDANCE DIAGNOSTICS (BMI vs DXA investigation)
  # =========================================================================
  # Step 2: BMI with ATTEMPT excluded (ACTIVE)
  "T1D_normal_vs_ow_obese_bmi_noattempt",
  "T1D_normal_vs_ow_obese_bmi_noattempt_adj_age",
  "T1D_normal_vs_ow_obese_bmi_noattempt_adj_age_sex",
  "T1D_nonobese_vs_obese_bmi_noattempt",

  # =========================================================================
  # CATEGORICAL (WHtR-defined, T1D within-group)
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
  # CATEGORICAL (HC vs T1D, WHtR-defined)
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
  # CATEGORICAL (HC vs T1D Overweight+Obese — combined; BMI/DXA/WHtR)
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
  # CONTINUOUS (WHtR)
  # =========================================================================
  "continuous/whtr_t1d",
  "continuous/whtr_all",
  "continuous/whtr_all_adj_group",
  "continuous/whtr_t1d_adj_age",
  "continuous/whtr_all_adj_age",
  "continuous/whtr_all_adj_group_age",
  "continuous/whtr_t1d_adj_age_sex",
  "continuous/whtr_all_adj_age_sex",
  "continuous/whtr_all_adj_group_age_sex"
)
# =============================================================================
# Pass 1: Compile run-level metadata CSVs (one row per analysis x cell type)
# Pass 2: Convert each *_processed.rds (full per-gene results) to a per-run
#         CSV containing logFC / p-value / FDR / Gene / annotations, and
#         upload alongside the rds at:
#             results/nebula/<subdir>/<celltype>/<celltype>_nebula_<suffix>.csv
# =============================================================================

# Toggle the per-run CSV export on/off (off => only the metadata compile runs).
EXPORT_PER_RUN_CSV <- TRUE

# Per-run CSV export filter:
#   - exclude all interaction models
#   - exclude all adjusted runs (categorical or continuous *_adj*)
#   - for continuous, keep only BMI, DXA body fat %, and WHtR (drop other DXA
#     variables like lean_mass, fat_kg, ag_ratio, est_vat, trunk_*, etc.)
#   - keep all unadjusted categorical contrasts
should_export_csv <- function(subdir) {
  if (grepl("^interaction/", subdir)) return(FALSE)
  if (grepl("_adj",          subdir)) return(FALSE)
  if (grepl("_noattempt",          subdir)) return(FALSE)
  if (grepl("^continuous/",  subdir)) {
    return(grepl("^continuous/(bmi|dexa_body_fat|whtr)_(t1d|all)$", subdir))
  }
  TRUE  # everything else (unadjusted categorical) gets exported
}

cat("Searching for NEBULA outputs in S3...\n")
all_meta        <- list()
subdirs_found   <- c()
subdirs_missing <- c()
n_csv_written       <- 0
n_csv_skipped       <- 0     # already exists on S3
n_csv_skipped_filter <- 0    # excluded by should_export_csv()
n_csv_failed        <- 0

for (subdir in analysis_subdirs) {
  prefix <- paste0(s3_base, "/", subdir, "/")
  objs <- tryCatch(
    get_bucket(bucket = bucket, prefix = prefix, region = "", max = 10000),
    error = function(e) { message(paste("Error listing", prefix, ":", e$message)); list() }
  )
  all_keys <- sapply(objs, function(x) x$Key)

  # ---- Pass 1: metadata CSVs ----
  meta_keys <- all_keys[grepl("_metadata\\.csv$", all_keys)]
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

  # ---- Pass 2: per-run CSV export (logFC, p, FDR, Gene per row) ----
  if (isTRUE(EXPORT_PER_RUN_CSV) && should_export_csv(subdir)) {
    rds_keys <- all_keys[grepl("_processed\\.rds$", all_keys)]
    for (rds_key in rds_keys) {
      csv_key <- sub("_processed\\.rds$", ".csv", rds_key)

      # Skip if a fresh CSV already exists and is newer than the rds (cheap rerun guard)
      already_have <- isTRUE(any(all_keys == csv_key))
      if (already_have) {
        n_csv_skipped <- n_csv_skipped + 1
        next
      }

      tryCatch({
        # Download rds, read, write csv, upload
        tmp_rds <- tempfile(fileext = ".rds")
        save_object(object = rds_key, bucket = bucket, region = "", file = tmp_rds)
        df_run <- readRDS(tmp_rds)
        unlink(tmp_rds)

        # Be permissive about classes (data.frame vs tibble vs list-with-results)
        if (is.list(df_run) && !is.data.frame(df_run) && "results" %in% names(df_run)) {
          df_run <- df_run$results
        }
        if (!is.data.frame(df_run)) {
          stop(sprintf("unexpected class for %s: %s",
                       rds_key, paste(class(df_run), collapse = "/")))
        }

        tmp_csv <- tempfile(fileext = ".csv")
        write.csv(df_run, tmp_csv, row.names = FALSE)
        put_object(file = tmp_csv, object = csv_key,
                   bucket = bucket, region = "")
        unlink(tmp_csv)
        n_csv_written <- n_csv_written + 1
      }, error = function(e) {
        n_csv_failed <<- n_csv_failed + 1
        message(sprintf("CSV export failed for %s: %s", rds_key, e$message))
      })
    }
  } else if (isTRUE(EXPORT_PER_RUN_CSV)) {
    # Subdir excluded by should_export_csv() filter (interaction, _adj, or
    # continuous variable not in {bmi, dexa_body_fat, whtr})
    n_csv_skipped_filter <- n_csv_skipped_filter +
      sum(grepl("_processed\\.rds$", all_keys))
  }
}

cat(sprintf("\n=== S3 scan summary ===\n"))
cat(sprintf("  Analysis types with results: %d / %d\n",
            length(subdirs_found), length(analysis_subdirs)))
if (length(subdirs_missing) > 0) {
  cat(sprintf("  Analysis types still pending: %d\n", length(subdirs_missing)))
  cat(paste("    -", subdirs_missing, collapse = "\n"))
  cat("\n")
}
cat(sprintf("Found %d metadata files\n", length(all_meta)))
if (isTRUE(EXPORT_PER_RUN_CSV)) {
  cat(sprintf("Per-run CSV export: %d written, %d already existed, %d filtered out, %d failed\n",
              n_csv_written, n_csv_skipped, n_csv_skipped_filter, n_csv_failed))
}

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