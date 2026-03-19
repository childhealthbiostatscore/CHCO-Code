# pulling DKA study for clinical trials data reporting
# Diabetic Kidney Alarm (16-1403)
# arm = T1D with DKA defined by ISPAD criteria

library(reticulate)
library(dplyr)
library(tidyr)
library(purrr)
library(REDCapR)

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "hhampson") {
  root_path <- "/Volumes/Peds Endo/"
} else {
  stop("Unknown user: please specify root path for this user.")
}

tokens <- read.csv(file.path(root_path, "/Data Harmonization/api_tokens.csv"))

redcap_uri <- "https://redcap.ucdenver.edu/api/"

# Define a named vector of tokens
study_tokens <- setNames(tokens$Token, tokens$Study)

token <- study_tokens[["DKA"]]

# Read data
dka_data <- redcap_read(
  redcap_uri = redcap_uri,
  token = token
)$data

dka_data_clean <- dka_data %>%
  filter(!is.na(mrn)) %>%
  dplyr::mutate(age_yrs = as.numeric((dt_enroll - dob)/365),
                sex = case_when(sex == 1 ~ "Male", sex == 2 ~ "Female"),
                ethnicity_hisp = case_when(race == 1 ~ "Non-Hispanic or Latino",
                                           race == 0 ~ "Hispanic or Latino",
                                           T ~ "Unknown/Not Reported"),
                qcr = case_when(floor(age) == 8 ~ 0.46,
                                floor(age) == 9 ~ 0.49,
                                floor(age) == 10 ~ 0.51,
                                floor(age) == 11 ~ 0.53,
                                floor(age) == 12 ~ 0.57,
                                floor(age) == 13 ~ 0.59,
                                floor(age) == 14 ~ 0.61, 
                                # Females
                                floor(age) == 15 & sex == "Female" ~ 0.64,
                                floor(age) == 16 & sex == "Female" ~ 0.67,
                                floor(age) == 17 & sex == "Female" ~ 0.69,
                                floor(age) == 18 & sex == "Female" ~ 0.69,
                                floor(age) >= 19 & sex == "Female" ~ 0.70,
                                # Males
                                floor(age) == 15 & sex == "Male" ~ 0.72,
                                floor(age) == 16 & sex == "Male" ~ 0.78,
                                floor(age) == 17 & sex == "Male" ~ 0.82,
                                floor(age) == 18 & sex == "Male" ~ 0.85,
                                floor(age) == 19 & sex == "Male" ~ 0.88,
                                floor(age) > 19 & sex == "Male" ~ 0.90,
                                T ~ floor(age)),
                eGFR_fas_cr_0_8hr = 107.3 / (screatinine_0_8hr / qcr),
                f1_0_8hr = screatinine_0_8hr / qcr,
                f2 = 1 - 0.5,
                f3_0_8hr = cystatinc_0_8hr / 0.82,
                eGFR_fas_cr_cysc_0_8hr = 107.3 / ((0.5 * f1_0_8hr + (f2 * f3_0_8hr))),
                
                eGFR_fas_cr_12_24hr = 107.3 / (screatinine_12_24hr / qcr),
                f1_12_24hr = screatinine_12_24hr / qcr,
                f3_12_24hr = cystatinc_12_24hr / 0.82,
                eGFR_fas_cr_cysc_12_24hr = 107.3 / ((0.5 * f1_12_24hr + (f2 * f3_12_24hr))),
                
                eGFR_fas_cr_3mo = 107.3 / (screatinine_3mo / qcr),
                f1_3mo = screatinine_3mo / qcr,
                f3_3mo = cystatinc_3mo / 0.82,
                eGFR_fas_cr_cysc_3mo = 107.3 / ((0.5 * f1_3mo + (f2 * f3_3mo))),
                aki = case_when()
  )

summary(dka_data_clean$sngal_0_8hr/dka_data_clean$eGFR_fas_cr_cysc_0_8hr)
summary(dka_data_clean$sngal_12_24hr/dka_data_clean$eGFR_fas_cr_cysc_12_24hr)
summary(dka_data_clean$sngal_3mo/dka_data_clean$eGFR_fas_cr_cysc_3mo)

summary(dka_data_clean$skim1_0_8hr)
summary(dka_data_clean$skim1_12_24hr)
summary(dka_data_clean$skim1_3mo)
