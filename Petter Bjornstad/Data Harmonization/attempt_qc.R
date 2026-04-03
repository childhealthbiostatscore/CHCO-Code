## QC script for attempt_clean.R
## Run after attempt_clean.R to verify all harmonized variables

library(dplyr)
library(tidyr)

egfr_vars    <- c("age", "height", "creatinine_s", "cystatin_c_s")
anthro_vars  <- c("waistcm", "hipcm", "waist_hip_ratio")
blood_vars   <- c("hba1c", "cholesterol", "hdl", "ldl", "triglycerides",
                  "ca_base", "c00009", "hct", "pltct", "sodium_s", "bun")
bold_vars    <- c("bold_l_bl_cortex", "bold_r_bl_cortex",
                  "bold_l_bl_kidney", "bold_r_bl_kidney",
                  "avg_c_r2", "avg_k_r2")
emu_vars     <- c("acr_u", "microalbumin_u", "creatinine_u",
                  "emu_1_acr_u", "emu_1_microalbumin_u", "emu_1_creatinine_u",
                  "emu_2_acr_u", "emu_2_microalbumin_u", "emu_2_creatinine_u",
                  "emu_3_acr_u", "emu_3_microalbumin_u", "emu_3_creatinine_u")
u24_vars     <- c("u24_na", "u24_vl", "u24_mab")
egfr_eq_vars <- c("eGFR_Schwartz", "eGFR_bedside_Schwartz", "eGFR_Zap",
                  "eGFR_fas_cr", "eGFR_fas_cr_cysc", "eGFR_CKD_epi",
                  "eGFR_CKiD_U25_Creat", "eGFR_CKiD_U25_CystatinC", "eGFR_CKiD_U25_avg")
cat_vars     <- c("sex", "ethnicity", "race", "group", "cgm_type",
                  "famhx_t1d", "fam_t1d", "famhx_htn", "fam_htn",
                  "famhx_hyperlipid", "fam_hyperlipid")

all_numeric_vars <- c(egfr_vars, anthro_vars, blood_vars, bold_vars,
                      emu_vars, u24_vars, egfr_eq_vars)
all_expected     <- c(all_numeric_vars, cat_vars, "diabetes_dx_date")

cat("=== VARIABLE PRESENCE CHECK ===\n")
missing_vars <- all_expected[!all_expected %in% names(merged_data)]
if (length(missing_vars) == 0) {
  cat("All expected variables present.\n\n")
} else {
  cat("Missing variables:\n")
  print(missing_vars)
  cat("\n")
}

cat("=== DUPLICATES (record_id x visit) ===\n")
dups <- merged_data %>% count(record_id, visit) %>% filter(n > 1)
if (nrow(dups) == 0) {
  cat("No duplicates found.\n\n")
} else {
  cat(nrow(dups), "duplicate record_id x visit combinations:\n")
  print(dups)
  cat("\n")
}

cat("=== VISITS PER SUBJECT ===\n")
print(table(count(merged_data, record_id)$n))
missing_baseline <- merged_data %>%
  group_by(record_id) %>%
  summarise(has_baseline = any(visit == "baseline", na.rm = TRUE)) %>%
  filter(!has_baseline)
cat(nrow(missing_baseline), "subject(s) missing a baseline visit\n\n")

cat("=== group ===\n")
print(table(merged_data$group, useNA = "ifany"))
cat("\n=== ethnicity ===\n")
print(table(merged_data$ethnicity, useNA = "ifany"))
cat("\n=== race ===\n")
print(table(merged_data$race, useNA = "ifany"))
cat("\n=== sex ===\n")
print(table(merged_data$sex, useNA = "ifany"))
unexpected_sex <- merged_data %>% filter(!sex %in% c("Male", "Female") | is.na(sex))
cat(nrow(unexpected_sex), "row(s) with unexpected/missing sex values\n\n")

cat("=== MISSINGNESS SUMMARY ===\n")
present_numeric <- all_numeric_vars[all_numeric_vars %in% names(merged_data)]
miss_summary <- merged_data %>%
  summarise(across(all_of(present_numeric), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(merged_data) * 100, 1)) %>%
  arrange(desc(pct_missing))
print(miss_summary, n = Inf)
cat("\n")

cat("=== RANGE SUMMARY ===\n")
range_summary <- merged_data %>%
  summarise(across(all_of(present_numeric), list(
    min  = ~ min(as.numeric(.x), na.rm = TRUE),
    max  = ~ max(as.numeric(.x), na.rm = TRUE),
    mean = ~ round(mean(as.numeric(.x), na.rm = TRUE), 2),
    sd   = ~ round(sd(as.numeric(.x), na.rm = TRUE), 2)
  ))) %>%
  pivot_longer(everything(),
               names_to  = c("variable", "stat"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value)
print(range_summary, n = Inf)
cat("\n")

cat("=== UNIT CONVERSION SPOT CHECKS ===\n")
spot_checks <- list(
  "creatinine_s (expect ~0.5-1.5 mg/dL)"   = merged_data$creatinine_s,
  "bun (expect ~5-30 mg/dL)"                = merged_data$bun,
  "hdl (expect ~20-100 mg/dL)"              = merged_data$hdl,
  "ldl (expect ~50-200 mg/dL)"              = merged_data$ldl,
  "triglycerides (expect ~50-500 mg/dL)"    = merged_data$triglycerides,
  "height (expect 100-250 cm)"              = merged_data$height,
  "waistcm (expect 50-200 cm)"              = merged_data$waistcm,
  "acr_u (expect ~0-1000 mg/g)"             = merged_data$acr_u
)
for (label in names(spot_checks)) {
  cat(label, ":\n")
  print(summary(as.numeric(spot_checks[[label]])))
  cat("\n")
}


# Checks that mrn, dob, screen_date, and study IDs are filled for all Denver ATTEMPT rows
if (exists("clean")) {
  cat("=== HARMONIZED DATASET FILL QC (Denver ATTEMPT rows) ===\n")

  id_cols <- c("mrn", "dob", "screen_date",
               "attempt_id", "casper_id", "coffee_id", "croc_id", "improve_id",
               "penguin_id", "rh_id", "rh2_id", "panther_id", "panda_id",
               "rpc2_id", "swht_id", "ultra_id", "co_enroll_id")

  denver_attempt <- clean %>%
    filter(startsWith(record_id, "3"), study == "ATTEMPT")

  cat("Total Denver ATTEMPT rows:", nrow(denver_attempt), "\n\n")

  # Missingness for each id column
  miss <- denver_attempt %>%
    summarise(across(any_of(id_cols), ~ sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
    mutate(pct_missing = round(n_missing / nrow(denver_attempt) * 100, 1)) %>%
    arrange(desc(n_missing))
  print(miss)
  cat("\n")

  # Show subjects where mrn is still missing
  missing_mrn <- denver_attempt %>%
    filter(is.na(mrn)) %>%
    distinct(record_id)
  cat(nrow(missing_mrn), "Denver ATTEMPT record_id(s) with no mrn at all:\n")
  if (nrow(missing_mrn) > 0) print(missing_mrn)

  # Show subjects where mrn is inconsistent within record_id
  mrn_inconsistent <- clean %>%
    filter(startsWith(record_id, "3")) %>%
    group_by(record_id) %>%
    summarise(n_mrn_values = n_distinct(mrn, na.rm = TRUE)) %>%
    filter(n_mrn_values > 1)
  cat("\n", nrow(mrn_inconsistent), "record_id(s) with >1 distinct mrn value:\n")
  if (nrow(mrn_inconsistent) > 0) print(mrn_inconsistent)
  cat("\n")
} else {
  cat("NOTE: 'clean' not found — run data_harmonization_with_soma.R first to check fill QC.\n\n")
}


if (exists("clean")) {
  cat("=== RAW vs HARMONIZED DATA QC ===\n\n")
  attempt_clean <- clean %>% filter(study == "ATTEMPT")

  cat("--- Record ID coverage ---\n")
  raw_ids  <- unique(merged_data$record_id)
  harm_ids <- unique(attempt_clean$record_id)
  missing_from_harm <- setdiff(raw_ids, harm_ids)
  extra_in_harm     <- setdiff(harm_ids, raw_ids)
  cat(length(raw_ids), "unique record_ids in raw merged_data\n")
  cat(length(harm_ids), "unique record_ids in harmonized clean (ATTEMPT)\n")
  cat(length(missing_from_harm), "record_id(s) in raw but missing from harmonized:\n")
  if (length(missing_from_harm) > 0) print(missing_from_harm)
  cat(length(extra_in_harm), "record_id(s) in harmonized but not in raw:\n")
  if (length(extra_in_harm) > 0) print(extra_in_harm)
  cat("\n")

  cat("--- Procedure distribution in harmonized ATTEMPT rows ---\n")
  print(table(attempt_clean$procedure, attempt_clean$visit, useNA = "ifany"))
  cat("\n")

  cat("--- MRI data check (bold_mri rows) ---\n")
  bold_rows <- attempt_clean %>% filter(procedure == "bold_mri")
  cat(nrow(bold_rows), "bold_mri rows\n")
  cat("Missing MRI values in bold_mri rows:\n")
  bold_rows %>%
    summarise(across(any_of(c("bold_l_bl_cortex", "bold_r_bl_cortex",
                               "bold_l_bl_kidney", "bold_r_bl_kidney")),
                     ~ sum(is.na(.)))) %>% print()
  cat("Non-missing MRI values in labs/screening rows (should be 0):\n")
  attempt_clean %>%
    filter(procedure %in% c("labs", "screening")) %>%
    summarise(across(any_of(c("bold_l_bl_cortex", "bold_r_bl_cortex",
                               "bold_l_bl_kidney", "bold_r_bl_kidney")),
                     ~ sum(!is.na(.)))) %>% print()
  cat("\n")

  cat("--- Key lab variables in labs/screening rows ---\n")
  key_labs <- c("creatinine_s", "hba1c", "acr_u", "age", "height", "weight")
  attempt_clean %>%
    filter(procedure %in% c("labs", "screening")) %>%
    summarise(across(any_of(key_labs), ~ sum(!is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_non_missing") %>%
    print()
  cat("\n")

  cat("--- Row counts per visit: raw vs harmonized ---\n")
  raw_counts  <- merged_data %>% count(visit, name = "n_raw")
  harm_counts <- attempt_clean %>% count(visit, procedure, name = "n_harmonized")
  print(left_join(harm_counts, raw_counts, by = "visit"))
  cat("\n")

  cat("--- Non-missing count per variable: raw vs harmonized ---\n")
  shared_vars <- intersect(names(merged_data), names(attempt_clean))
  shared_vars <- shared_vars[!shared_vars %in% c("record_id", "visit", "procedure", "study")]

  raw_counts_var <- merged_data %>%
    summarise(across(all_of(shared_vars), ~ sum(!is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_raw")

  harm_counts_var <- attempt_clean %>%
    summarise(across(all_of(shared_vars), ~ sum(!is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_harmonized")

  var_comparison <- left_join(raw_counts_var, harm_counts_var, by = "variable") %>%
    mutate(diff = n_harmonized - n_raw) %>%
    filter(n_raw > 0 | n_harmonized > 0) %>%
    arrange(diff)

  cat("Variables with fewer non-missing values in harmonized vs raw (possible data loss):\n")
  print(filter(var_comparison, diff < 0), n = Inf)
  cat("\nFull variable comparison (raw > 0):\n")
  print(filter(var_comparison, n_raw > 0), n = Inf)
  cat("\n")

  cat("--- Pass-through variable check (non-renamed/converted variables) ---\n")
  passthrough_vars <- list(
    compliance        = c("compliance_percent", "ondrug_days", "drugstoppage_days",
                          "pillcount_consumed", "pillcount_dispensed", "pillcount_returned",
                          "treatmentperiod_days"),
    diabetes_mgmt     = c("tdid_u", "tdid_basal_u", "tdid_bolus_u", "tdid_u_kg",
                          "tdid_basal_u_kg", "tdid_bolus_u_kg", "icr_average",
                          "isf_average_mmoll", "insulin_regimen", "aid_yn", "controliq_yn",
                          "inuslin_therapy", "csii_device", "aid_algorithm",
                          "blood_glucose_check", "insulin_basal_type", "insulin_bolus_type",
                          "target_glucose_low_mmoll", "target_glucose_high_mmoll"),
    dipstick_urine    = c("glucose_urine_dip", "bilirubin_urine_dip", "ketones_urine_dip",
                          "spgravity_urine_dip", "blood_urine_dip", "ph_urine_dip",
                          "urobilinogen_urine_dip", "protein_urine_dip", "nitrite_urine_dip",
                          "leukocytes_urine_dip"),
    central_blood_lab = c("glucose_serum_mmoll", "albumin_serum_gl", "uricacid_serum_umoll",
                          "magnesium_serum_mmoll", "bhb1_serum_mmoll", "tsh_serum_miul",
                          "pth_serum_pmoll"),
    local_blood_lab   = c("wbc_local", "rbc_local", "hgb_local", "hct_local",
                          "mcv_local", "mch_local", "mchc_local",
                          "potassium_blood_local", "chloride_blood_local",
                          "bicarbonate_local", "alt_local", "alp_local",
                          "bilirubin_local", "lipase_local"),
    mgfr              = c("gfr_raw_plasma", "mgfr_si_adjusted", "mgfr_bm_adult",
                          "mgfr_bm_adult_adjusted", "mgfr_jodal_bsa")
  )

  for (group in names(passthrough_vars)) {
    vars <- passthrough_vars[[group]]
    in_raw  <- vars[vars %in% names(merged_data)]
    in_harm <- vars[vars %in% names(attempt_clean)]
    missing_from_harm <- setdiff(in_raw, names(attempt_clean))
    missing_from_raw  <- setdiff(vars, names(merged_data))

    cat(sprintf("\n[%s]\n", group))
    if (length(missing_from_raw) > 0)
      cat("  Not in raw merged_data:     ", paste(missing_from_raw, collapse = ", "), "\n")
    if (length(missing_from_harm) > 0)
      cat("  Missing from harmonized:    ", paste(missing_from_harm, collapse = ", "), "\n")
    if (length(missing_from_raw) == 0 && length(missing_from_harm) == 0)
      cat("  All present in both raw and harmonized.\n")

    # Non-missing counts for variables present in both
    both <- intersect(in_raw, names(attempt_clean))
    if (length(both) > 0) {
      counts <- tibble(
        variable     = both,
        n_raw        = sapply(both, function(v) sum(!is.na(merged_data[[v]]))),
        n_harmonized = sapply(both, function(v) sum(!is.na(attempt_clean[[v]])))
      ) %>% filter(n_raw > 0 | n_harmonized > 0)
      if (nrow(counts) > 0) print(counts, n = Inf)
    }
  }
  cat("\n")

} else {
  cat("NOTE: 'clean' not found — run data_harmonization_with_soma.R first to check raw vs harmonized QC.\n\n")
}

cat("=== QC COMPLETE ===\n")


