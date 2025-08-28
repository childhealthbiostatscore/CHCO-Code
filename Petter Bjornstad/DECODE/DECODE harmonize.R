######################################
######################################

# DECODE KIDNEY Harmonizing Code

# List of studies:
# PANDA, RH, RH2, CASPER, CROC, IMPROVE
# PENGUIN, ATTEMPT, TODAY/TOODAY2,
# TEEN-LABS, T1-DISCO, REMODEL-T1D, RPC2

######################################
######################################

library(Hmisc)
library(dplyr)
library(labelled)

make_data_dictionary <- function(df, source_name = NULL, max_examples = 3) {
  # Coerce to plain data.frame (handles tibble/data.table)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  # Use the passed object name if source_name not provided
  if (is.null(source_name)) source_name <- deparse(substitute(df))
  
  # Helpers
  format_examples <- function(x, k = max_examples) {
    x <- x[!is.na(x)]
    if (inherits(x, "POSIXt")) x <- format(x, "%Y-%m-%d %H:%M:%S")
    if (inherits(x, "Date"))   x <- format(x, "%Y-%m-%d")
    if (is.factor(x))          x <- as.character(x)
    ux <- unique(x)
    if (!length(ux)) return(NA_character_)
    paste(utils::head(ux, k), collapse = ", ")
  }
  
  classes   <- vapply(df, function(x) paste(class(x), collapse = ","), character(1))
  n_missing <- vapply(df, function(x) sum(is.na(x)), integer(1))
  n_unique  <- vapply(df, function(x) length(unique(x)), integer(1))
  examples  <- vapply(df, format_examples, character(1))
  
  out <- data.frame(
    Source         = rep(source_name, ncol(df)),
    Variable       = names(df),
    Class          = classes,
    N_missing      = n_missing,
    N_unique       = n_unique,
    Example_values = examples,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Optional: include a dataset label if stored as an attribute
  ds_attr <- attr(df, "dataset_name")
  if (!is.null(ds_attr)) out$Dataset <- rep(ds_attr, nrow(out))
  
  out
}


# eGFR function
# Floor age and assign qcr based on age & sex
# qcr <- floor(data[[age]])
# 
# # Assign qcr values
# qcr[qcr == 8] <- 0.46
# qcr[qcr == 9] <- 0.49
# qcr[qcr == 10] <- 0.51
# qcr[qcr == 11] <- 0.53
# qcr[qcr == 12] <- 0.57
# qcr[qcr == 13] <- 0.59
# qcr[qcr == 14] <- 0.61
# 
# qcr[data[[sex]] == "F" & qcr == 15] <- 0.64
# qcr[data[[sex]] == "F" & qcr == 16] <- 0.67
# qcr[data[[sex]] == "F" & qcr == 17] <- 0.69
# qcr[data[[sex]] == "F" & qcr == 18] <- 0.69
# qcr[data[[sex]] == "F" & qcr >= 19] <- 0.70
# 
# qcr[data[[sex]] == "M" & qcr == 15] <- 0.72
# qcr[data[[sex]] == "M" & qcr == 16] <- 0.78
# qcr[data[[sex]] == "M" & qcr == 17] <- 0.82
# qcr[data[[sex]] == "M" & qcr == 18] <- 0.85
# qcr[data[[sex]] == "M" & qcr == 19] <- 0.88
# qcr[data[[sex]] == "M" & qcr > 19] <- 0.90

egfr_fas_cr <- function(serum_creatinine, qcr) {
  return(107.3 / (serum_creatinine / qcr))
}

egfr_fas_cr_cysc <- function(serum_creatinine, qcr, cystatin_c, alpha = 0.5) {
  f1 <- serum_creatinine / qcr
  f2 <- 1 - alpha
  f3 <- cystatin_c / 0.82
  return(107.3 / ((0.5 * f1) + (f2 * f3)))
}

egfr_zappitelli <- function(height, cystatin_c, serum_creatinine) {
  return((507.76 * exp(0.003 * height)) /
           ((cystatin_c ^ 0.635) * ((serum_creatinine * 88.4) ^ 0.547)))
}

egfr_schwartz <- function(height, serum_creatinine, cystatin_c, bun, sex) {
  m <- ifelse(sex == "M", 1, 0)
  return(39.1 * ((height / serum_creatinine) ^ 0.516) * 
           ((1.8 / cystatin_c) ^ 0.294) * 
           ((30 / bun) ^ 0.169) * 
           (1.099 ^ m) * ((height / 1.4) ^ 0.188))
}

egfr_bedside_schwartz <- function(height, serum_creatinine) {
  return((41.3 * (height / 100)) / serum_creatinine)
}

egfr_ckd_epi <- function(serum_creatinine, age, sex) {
  f <- ifelse(sex == "F", 1, 0)
  a <- ifelse(sex == "M", -0.302, -0.241)
  k <- ifelse(sex == "M", 0.9, 0.7)
  
  return(142 * 
           (pmin(serum_creatinine / k, 1) ^ a) * 
           (pmax(serum_creatinine / k, 1) ^ -1.200) * 
           (0.9938 ^ age) * (1.012 * f + (1 - f)))
}

# omit ongoing studies for now: REMODEL-T1D, T1D-DISCO and RPC2

# pull existing harmonized dataset
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}


harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na.strings = "")

harm_dat <- harm_dat %>%
  filter(study %in% c("PANDA", "RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE",
                      "CASPER", "CROCODILE", "PENGUIN")) # omit ATTEMPT from this for now since pulling from Antoine's data
harm_dat_dict <- make_data_dictionary(harm_dat, source_name = "harmonized", max_examples = 5)

### To pull from Antoine's dataset (merged_data)
load(file.path(root_path, "ATTEMPT/Data Clean/ATTEMPT_AC.RData"))
attempt <- merged_data %>%
  dplyr::select(-contains("date")) %>%
  mutate(eGFR_CKD_epi = egfr_ckd_epi((creatinine_s/100), age, sex),
         study = "ATTEMPT",
         group = "Type 1 Diabetes",
         waistcm = waist_m * 100,
         hdl = hdl_serum_mmoll * 38.67,
         ldl = ldl_serum_mmoll * 38.67,
         cholesterol = cholesterol_serum_mmoll * 38.67,
         creatinine_s = creatinine_serum_umoll / 10) %>%
  dplyr::rename(diabetes_duration = diabetes_dx_duration,
                waist_hip_ratio = waist_to_hip_ratio,
                eGFR_CKiD_U25_Creat = egfr_ckidu25_cr,
                eGFR_CKiD_U25_CystatinC = egfr_ckidu25_cysc)

attempt_dat_dict <- make_data_dictionary(attempt, source_name = "attempt", max_examples = 5)

## TODAY/TODAY2
# load comorbidity data (comorb)
# load(file.path(root_path, "TODAY subaward/Clinical data/comorb.Rdata"))

# load baseline risk factors (baserisk)
load(file.path(root_path, "TODAY subaward/Clinical data/TODAY/baserisk.Rdata"))

today <- baserisk %>%
  dplyr::select(-race, -sex) %>%
  dplyr::rename(record_id = releaseid, 
                age = AGEBASE,
                race = racedesc,
                cystatin_c_s = serumcystc,
                creatinine_s = SerumCreat,
                hba1c = HbA1c,
                sex = sex_char,
                hdl = HDL,
                eGFR_CKD_epi = EstCreatClear) %>%
  mutate(sex = case_when(sex == "M" ~ "Male", 
                         sex == "F" ~ "Female"),
         eGFR_CKD_epi = egfr_ckd_epi(creatinine_s, age, sex),
         study = "TODAY/TODAY2",
         visit = "baseline",
         group = "Type 2 Diabetes")

today_dat_dict <- make_data_dictionary(today, source_name = "today", max_examples = 5)

## TEEN-LABS
# TEENLABS
load(file.path(root_path, "Teen Labs/Data_Cleaned/analysis_dataset.RData"))
df <- remove_labels(df)

teenlabs <- df %>%
  dplyr::rename(record_id = ID,
                sex = SEX,
                ethnicity = ETHN,
                race = RACE,
                eGFR_fas_cr = eGFR.fas_cr,
                eGFR_fas_cr_cysc = eGFR.fas_cr_cysc,
                acr_u = UACRATIO,
                cystatin_c_s = CYSC,
                creatinine_s = CREAS,
                hba1c = HBA1C,
                hdl = HDL,
                ldl = LDL,
                pulse = hrate) %>%
  dplyr::mutate(group = case_when(diab_baseline == "Yes" ~ "Type 2 Diabetes",
                                  diab_baseline == "No" ~ "Obese Control"),
                acr_u = acr_u*100,
                eGFR_CKD_epi = egfr_ckd_epi(creatinine_s, age, sex),
                race = as.character(race),
                sex = as.character(sex),
                ethnicity = as.character(ethnicity),
                study = "TEEN-LABS") %>%
  dplyr::select(-(starts_with("seq")))

teenlabs_dat_dict <- make_data_dictionary(teenlabs, source_name = "teenlabs", max_examples = 5)

## T1-DISCO
# omit for now (ongoing)

## REMODEL-T1D
# omit for now (ongoing)

## RPC2
# omit for now (ongoing)

######################################


# Dictionary merge for better merging of dataset
decode_dict <- rbind(harm_dat_dict, attempt_dat_dict, today_dat_dict, teenlabs_dat_dict)

write.csv(decode_dict, file.path(root_path, "DECODE/Data Clean/decode_merged_dictionary.csv"), row.names = F)
# Merge

######################################

decode <- bind_rows(today, teenlabs, harm_dat, attempt)

decode <- decode %>%
  dplyr::mutate(race = 
                  case_when(
                    (is.na(race) & ethnicity %in% c("Black")) |
                    (race %in% c("Non-Hispanic Black", "Black or African-American",
                                "Black or African American")) ~ "Black or African American",
                    (is.na(race) & ethnicity %in% c("White")) | 
                      (race %in% c("Non-Hispanic White", "White or Caucasian", "White")) ~ "White",
                    race %in% c("More than one race", "More Than One", 
                                "Asian & Hawaiian/Pacific Islander", "Black/African American & White") ~ "More Than One",
                    is.na(race) | race == "Unknown" ~ "Unknown/Not Reported",
                    T ~ race),
                ethnicity = 
                  case_when(
                    race == "Hispanic" ~ "Hispanic or Latino",
                    ethnicity %in% c("Hispanic or Latino", 
                                     "Latin / South American", "Hispanic") ~ "Hispanic or Latino",
                    ethnicity %in% c("Not Hispanic or Latino", "White", "Black", "Mixed", 
                                     "Japanese / Korean / Filipino", 
                                     "Latin / South American", 
                                     "South East Asian", 
                                     "Arab / West Asian",
                                     "Non-Hispanic") ~ "Not Hispanic or Latino",
                    ethnicity %in% c("Unknown/Not Reported", "Unknown") ~ "Unknown/Not Reported",
                    is.na(ethnicity) ~ "Unknown/Not Reported"),
                race = case_when(race == "Hispanic" ~ "Other", T ~ race)
  )

save(decode, file = file.path(root_path, "DECODE/Data Clean/decode_harmonized.RData"))

decode_collapsed <- decode %>%
  dplyr::summarise(
    across(
      where(~ !any(sapply(.x, is.numeric))), 
      ~ ifelse(all(is.na(.x)), NA_character_, as.character(last(na.omit(.x))))
    ),
    across(
      where(~ all(sapply(.x, is.numeric) | is.na(.x))), 
      ~ ifelse(all(is.na(.x)), NA_real_, mean(as.numeric(.x), na.rm = TRUE))
    ),
    .by = c(record_id, visit)
  )

save(decode_collapsed, file = file.path(root_path, "DECODE/Data Clean/decode_harmonized_collapsed.RData"))

