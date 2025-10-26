# Script to create/clean dataset for analysis related to Y-T2D (RH, RH2, IMPROVE, PANTHER, etc)
# Author: Ye Ji Choi

# Import libraries
library(dplyr)
library(growthcleanr)
library(purrr)
library(tidyr)

# specify user for paths
user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
} else {
  stop("Unknown user: please specify root path for this user.")
}

# Pull in harmonized dataset and collapse
harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na.strings = "")

harm_dat_collapsed <- harm_dat %>%
  filter(age < 40) %>%
  mutate(case_when(visit == "screening" ~ "baseline", T ~ visit)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other")) %>%
  arrange(record_id) %>%
  mutate(agem = age*12)

# Add BMI percentile
bmi_percentile = ext_bmiz(data = subset(harm_dat_collapsed, 
                                        select = c("record_id", "visit", "sex", "agem", "weight", "height", "bmi")),
                          wt = "weight", ht = "height") %>%
  dplyr::select(record_id, visit, bmip, bmiz)
harm_dat_collapsed <- left_join(harm_dat_collapsed, bmi_percentile)

# Filter to Y-T2D
yt2dobese <- harm_dat_collapsed %>%
  filter(group != "Type 1 Diabetes") %>%
  filter(study %in% c("IMPROVE", "PANTHER", "RENAL-HEIR", "RENAL-HEIRitage")) %>%
  filter(visit == "baseline")

# Add visits (RH/IMPROVE -> RH2)
yt2dobese <- yt2dobese %>%
  arrange(mrn, date) %>%
  group_by(mrn) %>%
  dplyr::mutate(
    date = as.Date(date),
    # First calculate raw date differences
    raw_date_diff = as.numeric(difftime(date, lag(date), units = "days")),
    # Check if this MRN has any follow-up visits (date diff > 7 days)
    has_followup = any(raw_date_diff > 7, na.rm = TRUE),
    # Set date_diff_days to 0 for baseline only if there's a follow-up
    date_diff_days = case_when(
      is.na(raw_date_diff) & has_followup ~ 0,  # First visit when follow-up exists
      is.na(raw_date_diff) & !has_followup ~ NA_real_,  # First visit when no follow-up
      TRUE ~ raw_date_diff  # All other cases use the actual difference
    ),
    date_diff_years = date_diff_days/365,
    new_cluster = if_else(is.na(lag(date)) | date_diff_days > 7,
                          1, 0),
    cluster_id = cumsum(replace_na(new_cluster, 1)),
    visit_followup = case_when(
      date_diff_days > 7 ~ "follow_up",  # Explicitly a follow-up
      date_diff_days <= 7 & date_diff_days >= 0 ~ "baseline",  # Within 7 days
      is.na(date_diff_days) & has_followup ~ "baseline",  # First visit with follow-up
    )
  ) %>% ungroup() %>%
  group_by(mrn, cluster_id) %>%
  slice_min(order_by = date, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(uacr_100 = case_when(acr_u >= 100 ~ "UACR >=100 mg/g",
                              acr_u < 100 ~ "UACR <100 mg/g"),
         uacr_group = case_when(eGFR_fas_cr < 75 | acr_u >= 100 ~ "eGFR < 75 and/or UACR >=100 mg/g",
                                eGFR_fas_cr >= 75 | acr_u < 100 ~ "eGFR >= 75 and UACR <100 mg/g"),
         uacr_group = factor(uacr_group, levels = c("eGFR >= 75 and UACR <100 mg/g",
                                                    "eGFR < 75 and/or UACR >=100 mg/g")),
         ckd_group = factor(uacr_group, levels = c("eGFR >= 75 and UACR <100 mg/g",
                                                   "eGFR < 75 and/or UACR >=100 mg/g"),
                            labels = c("Non-CKD", "CKD")),
         # eGFR categories
         egfr_cat = case_when(
           eGFR_fas_cr < 60                      ~ "<60",
           eGFR_fas_cr >= 60 & eGFR_fas_cr < 90 ~ "60–<90",
           eGFR_fas_cr >= 90                     ~ "≥90",
           TRUE ~ NA_character_
         ),
         egfr_cat = factor(egfr_cat, levels = c("<60", "60–<90", "≥90")),
         # ACR categories
         acr_u_cat3 = case_when(
           acr_u < 30                  ~ "<30",
           acr_u >= 30 & acr_u < 300   ~ "30–<300",
           acr_u >= 300                ~ "≥300",
           TRUE ~ NA_character_
         ),
         acr_u_cat3 = factor(acr_u_cat3, levels = c("<30", "30–<300", "≥300")),
         bmi_cat = case_when(
           # Adults (≥18 years)
           age >= 18 & bmi < 25                ~ "Lean",
           age >= 18 & bmi >= 25 & bmi < 30    ~ "Overweight",
           age >= 18 & bmi >= 30               ~ "Obese",
           
           # Pediatrics (<18 years)
           age < 18 & bmip < 85      ~ "Lean",
           age < 18 & bmip >= 85 & bmip < 95 ~ "Overweight",
           age < 18 & bmip >= 95     ~ "Obese",
           
           TRUE ~ NA_character_
         ),
         bmi_cat = factor(
           bmi_cat,
           levels = c("Lean", "Overweight", "Obese")
         ),
         group_risk = case_when(bmip >= 95 | hba1c >=6 | group == "Type 2 Diabetes" ~ "High",
                                bmip < 85 & hba1c <= 5.6 ~ "Low"),
         tkv_combined = coalesce(total_kidney_volume_ml, total_kidney_volume_ml_manual),
         htadjtkv_combined = coalesce(ht_adj_tkv, ht_adj_tkv_manual),
         study_pairs = case_when(!is.na(rh_id) & !is.na(rh2_id) & !is.na(visit_followup) ~ "IMPROVE/RH -> RH2",
                           !is.na(panther_id) & !is.na(improve_id) & !is.na(visit_followup) ~ "PANTHER -> IMPROVE",
                           !is.na(panther_id) & !is.na(rh2_id) & !is.na(visit_followup) ~ "PANTHER -> RH2")
  ) %>% ungroup() %>%
  dplyr::summarise(
    # Handle factors separately
    across(where(is.factor), ~ {
      if(all(is.na(.x))) {
        factor(NA, levels = levels(.x))
      } else {
        last(na.omit(.x))
      }
    }),
    # Handle other non-numeric columns (characters)
    across(where(~ !is.numeric(.x) & !is.factor(.x)), ~ {
      ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))
    }),
    # Handle numeric columns
    across(where(is.numeric), ~ {
      ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))
    }),
    .by = c(mrn, visit_followup)
  )

# Determine DKD classification by N (using eGFR and UACR)
summary(yt2dobese$age) # use FAS eGFR
hist(yt2dobese$eGFR_fas_cr)
table(yt2dobese$eGFR_fas_cr < 75 | yt2dobese$acr_u >= 100)

# Save# SaveeGFR_CKD_epi
write.csv(yt2dobese, file.path(root_path, "Y-T2D/Data Clean/y_t2d_obese_dataset.csv"), row.names = F)
