# ATTEMPT file from Antoine cleaning
library(readxl)
library(dplyr)
library(tidyverse)

attempt_file = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx"
attempt_mri_file = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_MRI_SK_LHSC_TorontoLondon.csv"
attempt_mri_labs_file = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_MRI_HCT_SBP_TorontoLondon.csv"
# attempt data from Antoine on 3/14/25 after SOMAScan results
attempt_031425_raw <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_DenverDataRequest_20250314.csv")
# attempd TIR data from Antoine on 5/28/25
attempt_052825_raw <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_DenverDataRequest_20250528.csv")
attempt_061725_raw <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_DenverDataRequest_20250617.csv")

# non-Denver data formatting to match
attempt_mri <- read.csv(attempt_mri_file) %>%
  dplyr::rename(mri_r2_kidney_r = rk_r2s_mean,
                mri_r2_cortex_r = rc_r2s_mean,
                mri_r2_kidney_l = lk_r2s_mean,
                mri_r2_cortex_l = lc_r2s_mean) %>%
  dplyr::mutate(mri_date = as.Date(mri_date, format = "%m/%d/%y"),
                site = case_when(subject_id < 20000 ~ "Toronto",
                                 subject_id < 30000 ~ "London"))
attempt_mri_labs <- read.csv(attempt_mri_labs_file) %>%
  dplyr::mutate(mri_date = as.Date(mri_date, format = "%m/%d/%y"),
                date_visit = as.Date(date_visit, format = "%m/%d/%y"))
attempt_031425 <- attempt_031425_raw %>%
  dplyr::mutate(date_visit = as.Date(date_visit, format = "%m/%d/%y"),
                urine24h_start_date = as.Date(urine24h_start_date, format = "%m/%d/%y"),
                urine24h_stop_date = as.Date(urine24h_stop_date, format = "%m/%d/%y"),
                target_glucose_low_mmoll = case_when(grepl("<", target_glucose_low_mmoll) ~ "5",
                                                     T ~ target_glucose_low_mmoll),
                target_glucose_low_mmoll = as.numeric(target_glucose_low_mmoll),
                target_glucose_high_mmoll = case_when(grepl("<", target_glucose_high_mmoll) ~ "5",
                                                     T ~ target_glucose_high_mmoll),
                target_glucose_high_mmoll = as.numeric(target_glucose_high_mmoll),
                emu_1_albumin_mgl = as.character(emu_1_albumin_mgl),
                emu_2_albumin_mgl = as.character(emu_2_albumin_mgl),
                emu_3_albumin_mgl = as.character(emu_3_albumin_mgl),
                emu_1_acr_mgmmol = as.character(emu_1_acr_mgmmol),
                emu_2_acr_mgmmol = as.character(emu_2_acr_mgmmol),
                emu_3_acr_mgmmol = as.character(emu_3_acr_mgmmol),
                emu_urine_acr_mean_pooled = as.character(emu_urine_acr_mean_pooled),
                emu_urine_acr_mean = as.character(emu_urine_acr_mean))

dict <- read_excel(attempt_file, sheet = "Data Dictionary", na = c("NA", "")) %>%
  dplyr::select(Variable_Name, Label) %>%
  as.data.frame() 

# pull data from each sheet
demo <- read_excel(attempt_file, sheet = "ATTEMPT_Demographics", na = c("NA", ""))
anthro <- read_excel(attempt_file, sheet = "ATTEMPT_Anthropometrics", na = c("NA", ""))
medfamsmoking <- read_excel(attempt_file, sheet = "ATTEMPT_MedFamSmokingHx", na = c("NA", ""))
diabetesman <- read_excel(attempt_file, sheet = "ATTEMPT_DiabetesManagement", na = c("NA", ""))
glucosemon <- read_excel(attempt_file, sheet = "ATTEMPT_GlucoseMonitoring", na = c("NA", ""))
urinelab <- read_excel(attempt_file, sheet = "ATTEMPT_LocalUrineLabs", na = c("NA", ""))
urine_24h <- read_excel(attempt_file, sheet = "ATTEMPT_Urine24h", na = c("NA", "")) %>%
  mutate(date_visit = urine24h_start_date,
         urine24h_start_time = format(urine24h_start_time, "%H:%M"),
         urine24h_stop_time = format(urine24h_stop_time, "%H:%M"),
         urine24h_volume_litres = case_when(urine24h_volume_litres == "UNK" ~ "",
                                            T ~ urine24h_volume_litres),
         urine24h_volume_litres = as.numeric(urine24h_volume_litres),
         urine24h_volume_24h_litres = case_when(urine24h_volume_24h_litres == "UNK" ~ "",
                                                T ~ urine24h_volume_24h_litres),
         urine24h_volume_24h_litres = as.numeric(urine24h_volume_24h_litres))
urineemu <- read_excel(attempt_file, sheet = "ATTEMPT_UrineEMU", na = c("NA", ""))
bloodlab_local <- read_excel(attempt_file, sheet = "ATTEMPT_LocalBloodLabs", na = c("NA", ""))
bloodlab_central <- read_excel(attempt_file, sheet = "ATTEMPT_CentralBloodLabs", na = c("NA", "")) %>%
  mutate(hba1c_percent = hba1c_percent*100)
compliance <- read_excel(attempt_file, sheet = "ATTEMPT_Compliance", na = c("NA", ""))
egfr <- read_excel(attempt_file, sheet = "ATTEMPT_eGFR", na = c("NA", ""))
# mgfr <- read_excel(attempt_file, sheet = "ATTEMPT_mGFR", na = c("NA", ""))
mgfr <- attempt_061725_raw # newer mgfr file
mri <- read_excel(attempt_file, sheet = "ATTEMPT_BoldMRI", na = c("NA", "")) %>%
  bind_rows(attempt_mri)

randomization <- read_excel("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_Randomization_Denver.xlsx",
                            sheet = "ATTEMPT_Randomization", na = c("NA", ""))

data_frames <- list(attempt_031425 = attempt_031425,
  demo = demo, randomization = randomization, anthro = anthro, 
  medfamsmoking = medfamsmoking, diabetesman = diabetesman, 
  glucosemon = glucosemon, urinelab = urinelab, 
  urineemu = urineemu, urine_24h = urine_24h, 
  bloodlab_local = bloodlab_local, 
  bloodlab_central = bloodlab_central, compliance = compliance, 
  egfr = egfr, mgfr = mgfr, mri = mri, attempt_mri_labs = attempt_mri_labs,
  tir = attempt_052825_raw)

data_frames <- lapply(data_frames, function(df) {
  df %>%
    dplyr::mutate(visit = case_when(str_detect(visit, "V1|R1") ~ "screening", # week -4
                                    str_detect(visit, "V2") ~ "baseline", # week 0
                                    str_detect(visit, "V3") ~ "4_weeks_post", # week 4
                                    str_detect(visit, "V4|R3") ~ "4_months_post", # week 16
                                    str_detect(visit, "V5") ~ "follow_up", # week 18
                                    T ~ visit))
})

# Define categories
# Race/ethnicity
eth_names <- c("White",
               "Black",
               "Latin / South American",
               "Arab / West Asian",
               "Japanese / Korean / Filipino",
               "South East Asian",
               "Mixed")

merged_data <- purrr::reduce(data_frames, ~ full_join(.x, .y)) %>%
  group_by(subject_id) %>%
  fill(treatment_arm, .direction = "downup") %>% ungroup() %>%
  group_by(subject_id, visit) %>%
  fill(., .direction = "downup") %>% ungroup() %>%
  ungroup() %>% rowwise() %>%
  dplyr::mutate(height = height_m * 100,
                randomization_time= format(randomization_time, "%H:%M"), 
                across(where(is.logical), ~ ifelse(.x, "Yes", "No")),
                ethnicity = case_when(ethnicity_categorized == 1 ~ eth_names[1], 
                                      ethnicity_categorized == 2 ~ eth_names[2],
                                      ethnicity_categorized == 3 ~ eth_names[3], 
                                      ethnicity_categorized == 4 ~ eth_names[4], 
                                      ethnicity_categorized == 5 ~ eth_names[5], 
                                      ethnicity_categorized == 6 ~ eth_names[6], 
                                      ethnicity_categorized == 7 ~ eth_names[7]),
                sex = case_when(sex == 1  ~ "Male", sex == 2 ~ "Female"),
                treatment_arm = case_when(treatment_arm == "A" ~ "Placebo",
                                          treatment_arm == "B" ~ "Dapagliflozin 5mg"),
                microalbumin_urine24h = case_when(microalbumin_urine24h == "< 5"  ~ "2.5",
                                                  T ~ microalbumin_urine24h),
                microalbumin_urine24h = as.numeric(microalbumin_urine24h),
                microalbumin_urine24h_mgL = microalbumin_urine24h * 1000,
                creatinine_urine24h_gL = creatinine_urine24h * 0.11312,
                # uaer_24 = microalbumin_urine24h_mgL / creatinine_urine24h_gL,
                mri_r2_cortex_l = case_when(!is.na(mri_r2_cortex_l) ~ as.numeric(mri_r2_cortex_l)),
                mri_r2_cortex_r = case_when(!is.na(mri_r2_cortex_l) ~ as.numeric(mri_r2_cortex_r)),
                mri_r2_kidney_l = case_when(!is.na(mri_r2_cortex_l) ~ as.numeric(mri_r2_kidney_l)),
                mri_r2_kidney_r = case_when(!is.na(mri_r2_cortex_l) ~ as.numeric(mri_r2_kidney_r)),
                avg_c_r2 = mean(c(mri_r2_cortex_l, mri_r2_cortex_r)),
                avg_k_r2  = mean(c(mri_r2_kidney_l, mri_r2_kidney_r)),
                emu_1_albumin_mgl = case_when(grepl("<", emu_1_albumin_mgl) ~ "2.5",
                                              T ~ emu_1_albumin_mgl),
                emu_1_albumin_mgl = as.numeric(emu_1_albumin_mgl), 
                emu_1_acr_mgmmol = emu_1_albumin_mgl/emu_1_creatinine_umoll*1000,
                emu_2_albumin_mgl = case_when(grepl("<", emu_2_albumin_mgl) ~ "2.5",
                                              T ~ emu_2_albumin_mgl),
                emu_2_albumin_mgl = as.numeric(emu_2_albumin_mgl), 
                emu_2_acr_mgmmol = emu_2_albumin_mgl/emu_2_creatinine_umoll*1000,
                emu_3_albumin_mgl = case_when(grepl("<", emu_3_albumin_mgl) ~ "2.5",
                                              T ~ emu_3_albumin_mgl),
                emu_3_albumin_mgl = as.numeric(emu_3_albumin_mgl), 
                emu_3_acr_mgmmol = emu_3_albumin_mgl/emu_3_creatinine_umoll*1000,
                emu_urine_acr_mean = mean(c(emu_1_acr_mgmmol, emu_2_acr_mgmmol, emu_3_acr_mgmmol)),
                subject_id = as.character(subject_id),
                age = coalesce(age, age_baseline),
                sbp = coalesce(sbp_mmhg, mri_sbp),
                dbp = coalesce(dbp_mmhg),
                hba1c = hba1c_percent,
                hct = coalesce(hct_ll, hct_percent_us/100)) %>% ungroup() %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA, last(na.omit(.x)))),
                   across(where(is.numeric),  ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(subject_id, visit)) %>%
  # filter(visit %in% c("baseline", "4_months_post")) %>%
  dplyr::select(-c(age_baseline, sbp_mmhg, mri_sbp, dbp_mmhg, hct_ll, hct_percent_us)) %>%
  dplyr::rename(record_id = subject_id,
                date = date_visit,
                weight = weight_kg, 
                bmi = bmi_kgm2,
                temp = body_temp_celsius, 
                pulse = heart_rate_bpm, 
                diabetes_dx_duration = t1d_duration, 
                creatinine_s = creatinine_blood_local,
                cystatin_c_s = cystatin_c_serum_mgl)

date_cols <- grep("(^date$|_date$|^date_|_date_)", names(merged_data), value = TRUE)
time_cols <- grep("(^time$|time$|^time_|_time_)", names(merged_data), value = TRUE)

# Convert numeric Unix timestamps to formatted date strings
merged_data[date_cols] <- lapply(merged_data[date_cols], function(x) {
  format(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"), "%m/%d/%y")
})

merged_data[time_cols] <- lapply(merged_data[time_cols], function(x) {
  # Only convert if numeric
  if (is.numeric(x)) {
    format(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"), "%H:%M")
  } else {
    x  # leave as-is if not numeric
  }
})

# remove columns that are all NAs
merged_data <- merged_data[, colSums(!is.na(merged_data)) > 0]

merged_data <- merged_data[, !grepl("^ethnicity_ca___|^ethnicity_us___", names(merged_data))]

# merge and save
save(merged_data, file = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Clean/ATTEMPT_AC.RData")

