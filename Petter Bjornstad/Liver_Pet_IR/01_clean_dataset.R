# Could you please make a dataset for Renal heir/heritage participants that includes the following variables 
# "group","age","sex","bmi","hba1c","liver_k1","liver_k2","liver_k3","m_i_gir_190"  
# *Note, shivani made these varialbes: m_i_p2_raw_lean (CROC) --> caluculated using p2_raw_leanm and p2_steady_state_insulin
# m_i_gir_190 (RH) --> calculated using gir_190 and steady_state_insulin. 
# im not sure how she came up with these equations or how to do that but thats what i think im supposed to use. 
# We want to only include participants that had a clamp and liver pet at the same time period, not coenrolled and done at different time periods.

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "hhampson") {
  # root_path <- "/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
  root_path <- "/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive"
} else {
  stop("Unknown user: please specify root path for this user.")
}

library(dplyr)
library(purrr)

harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na.strings = "")

harm_dat_collapsed <- harm_dat %>%
  group_by(record_id, visit) %>%
  tidyr::fill(date, .direction = "updown") %>% ungroup() %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other"))

rh_rh2_df <- harm_dat_collapsed %>%
  filter(!is.na(rh_id) | !is.na(rh2_id)) %>%
  dplyr::select(record_id, visit, mrn, group, age, sex, bmi, hba1c, starts_with("liver"), m_i, gir_190, m_i_gir_190) %>%
  mutate(
    age_liver = if_else(!is.na(liver_k1) & !is.na(liver_k2) & !is.na(liver_k3), age, NA_real_),
    age_ins   = if_else(!is.na(m_i) | !is.na(gir_190), age, NA_real_)
  ) %>%
  group_by(mrn) %>%
  tidyr::fill(starts_with("liver"), m_i, gir_190, m_i_gir_190, age_liver, age_ins, .direction = "updown") %>%
  arrange(mrn) %>%
  filter(!is.na(age_liver) | !is.na(age_ins)) %>%
  distinct(mrn, age_liver, age_ins, .keep_all = T) %>%
  mutate(age_diff = age_liver - age_ins)

rm(harm_dat,harm_dat_collapsed)
