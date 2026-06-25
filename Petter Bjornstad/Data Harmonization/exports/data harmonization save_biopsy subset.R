library(dplyr)
library(tidyr)
library(table1)
library(arsenal)
library(Hmisc)
library(purrr)

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

biopsy_dat <- dat %>%
  filter(procedure == "kidney_biopsy") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(participation_status != "Removed"|is.na(participation_status)) %>%
  mutate(date = as.character(date),
         race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other"),
         obtained_lm_em = if_else(rowSums(across(ends_with("core_spec___1"))) > 0, "Yes", "No"),
         obtained_cryo_hypotherm = if_else(rowSums(across(ends_with("core_spec___3"))) > 0, "Yes", "No"),
         obtained_oct = if_else(rowSums(across(ends_with("core_spec___2"))) > 0, "Yes", "No"),
         obtained_rnalater = if_else(rowSums(across(ends_with("core_spec___4"))) > 0, "Yes", "No"),
         biopsy_yn = if_else(rowSums(!is.na(dplyr::select(., date, kit_id, glut_id, form_id, rnalater_id, cryomold_id, cryostor_id, ln2_id))) > 0, "Yes", "No")) %>%
  dplyr::select(record_id, co_enroll_id, biopsy_yn, date, group, visit, sglt2i_ever, sglti_timepoint, 
                kit_id,glut_id, form_id, rnalater_id, cryomold_id, cryostor_id, ln2_id) %>%
  arrange(desc(biopsy_yn), record_id)

tokens <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
rpc2_token <- subset(tokens, Study == "RPC2")$Token

rpc2 <- REDCapR::redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                     token = rpc2_token)$data

rpc2_biopsy <- rpc2 %>%
  mutate(bx_date = as.character(bx_date),
         group = "Type 2 Diabetes",
         obtained_lm_em = if_else(rowSums(across(ends_with("core_spec___1"))) > 0, "Yes", "No"),
         obtained_cryo_hypotherm = if_else(rowSums(across(ends_with("core_spec___3"))) > 0, "Yes", "No"),
         obtained_oct = if_else(rowSums(across(ends_with("core_spec___2"))) > 0, "Yes", "No"),
         obtained_rnalater = if_else(rowSums(across(ends_with("core_spec___4"))) > 0, "Yes", "No"),
         biopsy_yn = if_else(rowSums(!is.na(dplyr::select(., bx_kit_id, bx_glut_id, bx_form_id, bx_rnalater_id, bx_cryomold_id, bx_cryostor_id))) > 0, "Yes", "No")) %>%
  dplyr::select(subject_id, group, biopsy_yn, bx_date, redcap_event_name, 
                bx_kit_id,bx_glut_id, bx_form_id, bx_rnalater_id, bx_cryomold_id, bx_cryostor_id) %>%
  rename("record_id" = subject_id,
         "date" = bx_date,
         "visit" = redcap_event_name,
         "kit_id" = bx_kit_id,
         "glut_id" = bx_glut_id,
         "form_id" = bx_form_id,
         "rnalater_id" = bx_rnalater_id,
         "cryomold_id" = bx_cryomold_id,
         "cryostor_id" = bx_cryostor_id) %>%
  filter(!is.na(date))
  
biopsy_dat_combined <- full_join(biopsy_dat, rpc2_biopsy) %>%
  arrange(desc(biopsy_yn), record_id)

write.csv(biopsy_dat_combined, 
          "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/harmonized_dataset_biopsy_subset.csv", 
          row.names = F,
          na = "")
