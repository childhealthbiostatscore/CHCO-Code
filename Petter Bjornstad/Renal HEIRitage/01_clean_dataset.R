# clean and create dataset for RH/RH2/IMPROVE analysis (DKD vs. non-DKD)

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "yejichoi") {
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "hhampson") {
  root_path <- "/Volumes/Peds Endo"
} else {
  stop("Unknown user: please specify root path for this user.")
}

library(dplyr)
library(tidyverse)
library(purrr)
library(Hmisc)
library(growthcleanr)
library(aws.s3)
library(jsonlite)

## Create an S3 client

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na.strings = "")

harm_dat_collapsed <- harm_dat %>%
  mutate(visit = case_when(visit == "screening" ~ "baseline", T ~ visit)) %>%
  group_by(record_id, visit) %>%
  fill(date, .direction = "updown") %>% ungroup() %>%
  group_by(mrn, kit_id) %>% 
  fill(c(lv_cardiac_index, lv_cardiac_output, 
         lv_hr, lv_myo_thickness_dias,
         lv_myo_thickness_total, lv_stroke_volume,
         lved_mass, lvedv, lvef, lves_mass, 
         lvesv, rv_cardiac_index, rv_cardiac_output, 
         rv_hr, rv_stroke_volume, rvedv, rvef, rvesv), .direction = "updown") %>% ungroup() %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other"),
                uacr_100 = case_when(acr_u >= 100 ~ "UACR >=100 mg/g",
                                     acr_u < 100 ~ "UACR <100 mg/g"),
                uacr_group = case_when(
                  eGFR_CKD_epi < 90 | acr_u >= 100 ~ "eGFR < 90 and/or UACR >=100 mg/g",
                  eGFR_CKD_epi >= 90 & acr_u < 100 ~ "eGFR >= 90 and UACR <100 mg/g",
                  TRUE ~ NA_character_
                ),
                uacr_group = factor(
                  uacr_group,
                  levels = c("eGFR >= 90 and UACR <100 mg/g",
                             "eGFR < 90 and/or UACR >=100 mg/g")
                ),
                dkd_group_100 = case_when(
                  eGFR_CKD_epi < 90 | acr_u >= 100 ~ "DKD",
                  eGFR_CKD_epi >= 90 & acr_u < 100 ~ "non_DKD",
                  TRUE ~ NA_character_
                ),
                dkd_group_100 = factor(dkd_group_100, levels = c("non_DKD", "DKD")),
                
                dkd_group_30 = case_when(
                  eGFR_CKD_epi < 90 | acr_u >= 30 ~ "DKD",
                  eGFR_CKD_epi >= 90 & acr_u < 30 ~ "non_DKD",
                  TRUE ~ NA_character_
                ),
                dkd_group_30 = factor(dkd_group_30, levels = c("non_DKD", "DKD")),
                # eGFR categories
                egfr_cat = case_when(
                  eGFR_CKD_epi < 60                      ~ "<60",
                  eGFR_CKD_epi >= 60 & eGFR_CKD_epi < 90 ~ "60–<90",
                  eGFR_CKD_epi >= 90                     ~ "≥90",
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
                                       bmip < 85 & hba1c <=5.6 ~ "Low"),
                group_risk_lowyn = case_when(group_risk == "Low" ~ "Low", 
                                             T ~ "Not Low"), 
                tkv_combined = coalesce(total_kidney_volume_ml, total_kidney_volume_ml_manual),
                htadjtkv_combined = coalesce(ht_adj_tkv, ht_adj_tkv_manual),
                glp_t2dob = case_when(group == "Lean Control" ~ "HC",
                                      epic_glp1ra_1 == "Yes" ~ "GLP_Y",
                                      epic_glp1ra_1 == "No" ~ "GLP_N"),
                sglt2i_t2dob = case_when(group == "Lean Control" ~ "HC",
                                         epic_sglti2_1 == "Yes" ~ "SGLT2i_Y",
                                         epic_sglti2_1 == "No" ~ "SGLT2i_N"),
                dkd_group_30_hc = case_when(group == "Lean Control" ~ "HC",
                                            T ~ dkd_group_30),
                dkd_group_100_hc = case_when(group == "Lean Control" ~ "HC",
                                      T ~ dkd_group_100)
  )

# RH/RH2 subset
rh_rh2_croc_panther_unique <- harm_dat_collapsed %>% 
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "CROCODILE", "PANTHER")) %>%
  filter(group %in% c("Type 2 Diabetes", "Obese Control", "Lean Control")) %>%
  filter(visit == "baseline") %>%
  arrange(mrn, date) %>%
  distinct(mrn, date, .keep_all = TRUE) %>%
  group_by(mrn) %>%
  mutate(
    date = as.Date(date),
    new_cluster = if_else(is.na(lag(date)) | difftime(date, lag(date), units = "days") > 7,
                          1, 0),
    cluster_id = cumsum(replace_na(new_cluster, 1))
  ) %>%
  group_by(mrn, cluster_id) %>%
  slice_min(order_by = date, with_ties = FALSE) %>%
  ungroup()

rh_rh2_croc_improve_unique <- harm_dat_collapsed %>% 
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "CROCODILE", "IMPROVE")) %>%
  filter(group %in% c("Type 2 Diabetes", "Obese Control", "Lean Control")) %>%
  filter(visit == "baseline") %>%
  arrange(mrn, date) %>%
  distinct(mrn, date, .keep_all = TRUE) %>%
  group_by(mrn) %>%
  mutate(
    date = as.Date(date),
    new_cluster = if_else(is.na(lag(date)) | difftime(date, lag(date), units = "days") > 7,
                          1, 0),
    cluster_id = cumsum(replace_na(new_cluster, 1))
  ) %>%
  group_by(mrn, cluster_id) %>%
  slice_min(order_by = date, with_ties = FALSE) %>%
  fill(kit_id, .direction = "downup") %>%
  ungroup()

# Save CSV
write.csv(rh_rh2_croc_panther_unique, 
          file = file.path(root_path, "Renal HERITAGE/Data_Cleaned/rh_rh2_croc_panther_unique.csv"),
                           row.names = F)

s3saveRDS(rh_rh2_croc_panther_unique, 
          object = "data_clean/rh_rh2_croc_panther_unique.RDS", 
          bucket = "harmonized.dataset",
          region = "")

write.csv(rh_rh2_croc_improve_unique, 
          file = file.path(root_path, "Renal HERITAGE/Data_Cleaned/rh_rh2_croc_improve_unique.csv"),
          row.names = F)

s3saveRDS(rh_rh2_croc_improve_unique, 
          object = "data_clean/rh_rh2_croc_improve_unique.RDS", 
          bucket = "harmonized.dataset",
          region = "")

rh_rh2_unique <- rh_rh2_croc_panther_unique %>%
  filter(group %nin% c("Lean Control")) %>%
  filter(study != "PANTHER")

write.csv(rh_rh2_unique, 
          file = file.path(root_path, "Renal HERITAGE/Data_Cleaned/rh_rh2_unique.csv"),
          row.names = F)

s3saveRDS(rh_rh2_unique, 
          object = "data_clean/rh_rh2_unique.RDS", 
          bucket = "harmonized.dataset",
          region = "")
