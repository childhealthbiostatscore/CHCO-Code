# clean and create dataset for RH/RH2/IMPROVE analysis (DKD vs. non-DKD)

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "hhampson") {
  root_path <- "/Volumes/Peds Endo/"
} else {
  stop("Unknown user: please specify root path for this user.")
}


library(dplyr)
library(tidyverse)
library(purrr)
library(growthcleanr)

harm_dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

harm_dat_collapsed <- harm_dat %>%
  group_by(record_id, visit) %>%
  fill(date, .direction = "updown") %>% ungroup() %>%
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

# RH/RH2 subset
rh_rh2_croc_panther <- harm_dat_collapsed %>% 
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "CROCODILE", "PANTHER")) %>%
  filter(group %in% c("Type 2 Diabetes", "Obese Control", "Lean Control")) %>%
  filter(visit == "baseline")

bmi_percentile = ext_bmiz(data = subset(rh_rh2_croc_panther, 
                                        select = c("record_id", "sex", "agem", 
                                                   "weight", "height", "bmi")),
                          wt = "weight", ht = "height") %>%
  dplyr::select(record_id, bmip)

rh_rh2_croc_panther_unique <- rh_rh2_croc_panther %>%
  left_join(bmi_percentile) %>%
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
  ungroup() %>%
  #stratify T2D and OB in RH/RH2 by UACR >= 100 mg/g and look at the kidney MRI parameters.. I want to include some of this preliminary data in my slides for Spain
  # eGFR <60 for eGFR_CKD_epi AND/OR UACR>= 100 mg/g
  # in other words, to be in the high risk group EITHER eGFR <60 AND/OR UACR>= 100mg/g
  # so if you have UACR >= 100 mg/g but eGFR >60 you qualify, or if you have UACR <100 but eGFR <60 you qualify, and if you have both, i.e., UACR >= 100 and eGFR <60 you qualify
  mutate(uacr_100 = case_when(acr_u >= 100 ~ "UACR >=100 mg/g",
                              acr_u < 100 ~ "UACR <100 mg/g"),
         uacr_group = case_when(eGFR_CKD_epi < 75 | acr_u >= 100 ~ "eGFR < 75 and/or UACR >=100 mg/g",
                                eGFR_CKD_epi >=75 | acr_u < 100 ~ "eGFR >= 75 and UACR <100 mg/g"),
         uacr_group = factor(uacr_group, levels = c("eGFR >= 75 and UACR <100 mg/g",
                                                    "eGFR < 75 and/or UACR >=100 mg/g")),
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
         tkv_combined = coalesce(total_kidney_volume_ml, total_kidney_volume_ml_manual),
         htadjtkv_combined = coalesce(ht_adj_tkv, ht_adj_tkv_manual),
  )

write.csv(rh_rh2_croc_panther_unique, file.path(root_path, ""))
