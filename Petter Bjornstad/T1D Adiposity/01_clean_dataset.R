# T1D adiposity analysis

# create clinical dataset
library(tidyverse)
library(arsenal)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(Hmisc)
library(readxl)
library(purrr)
library(aws.s3)
library(jsonlite)

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

## Create an S3 client

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/soma_harmonized_dataset.csv"), na.strings = "")

t1d_hc_dat_soma <- harm_dat %>%
  filter(group %in%  c("Type 1 Diabetes","Lean Control") | study == "ATTEMPT") %>%
  filter(study %in% c("CROCODILE", "ATTEMPT", "CASPER", "RENAL-HEIR", "PANDA")) %>%
  dplyr::mutate(visit = case_when(study == "ATTEMPT" & visit == "screening" ~ "baseline",
                                  T ~ visit),
                keep = case_when(study == "RENAL-HEIR" & age > 12 & age < 21 ~ T,
                                 study != "RENAL-HEIR" ~ T,
                                 T ~ F)) %>%
  filter(visit == "baseline") %>%
  filter(keep) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  dplyr::mutate(bmi_obesity = case_when(age > 19 & bmi < 25 ~ "Normal",
                                        age > 19 & bmi >= 25 & bmi < 30 ~ "Overweight",
                                        age > 19 & bmi >= 30 ~ "Obese",
                                        age <= 19 & bmip < 85 ~ "Normal",
                                        age <= 19 & bmip >= 85 & bmip < 95 ~ "Overweight",
                                        age <= 19 & bmip >= 95 ~ "Obese"),
                dxa_obesity = case_when(dexa_body_fat < 25 & sex == "Male" ~ "Normal",
                                        dexa_body_fat < 32 & sex == "Female" ~ "Normal",
                                        dexa_body_fat >= 25 & dexa_body_fat < 30 & sex == "Male" ~ "Overweight",
                                        dexa_body_fat >= 32 & dexa_body_fat < 35 & sex == "Female" ~ "Overweight",
                                        dexa_body_fat >= 30 & sex == "Male" ~ "Obese",
                                        dexa_body_fat >= 35 & sex == "Female" ~ "Obese"),
                bmi_obesity = factor(bmi_obesity, levels = c("Normal", "Overweight", "Obese")),
                dxa_obesity = factor(dxa_obesity, levels = c("Normal", "Overweight", "Obese")),
                # mgfr_raw_combined = coalesce(gfr_raw_plasma, mgfr_jodal),
                # mgfr_bsa_combined = coalesce(gfr_bsa_plasma, mgfr_jodal_bsa)
                ) %>% # need attempt mgfr
  dplyr::select(-c(mrn, dob)) %>%
  filter(record_id %nin% c("PNDA-105", "PNDA-127")) # screen failed

t1d_hc_dat <- t1d_hc_dat_soma %>%
  dplyr::select(-starts_with("seq"))

# test <- t1d_hc_dat %>% select(record_id, age, bmi, bmip, sex, weight, height, dexa_body_fat, bmi_obesity, dxa_obesity)
# current issue (3/10: missing ATTEMPT not from Denver site, checking with Dawn on PNDA 208 who may not have anthro measures)

s3saveRDS(t1d_hc_dat_soma, 
          object = "data_clean/t1d_hc_clinical_data_soma.csv", 
          bucket = "t1d.adiposity",
          region = "",
          multipart = T)

s3saveRDS(t1d_hc_dat, 
          object = "data_clean/t1d_hc_clinical_data.csv", 
          bucket = "t1d.adiposity",
          region = "")
