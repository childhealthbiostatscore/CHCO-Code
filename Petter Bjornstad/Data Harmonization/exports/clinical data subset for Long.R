library(dplyr)
library(tidyr)
library(Hmisc)
library(purrr)

dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))
st_ids <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Spatial transcriptomics/spatial_transcriptomics_ids.csv")
pb90_ids <- readRDS("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/scRNA/data_clean/pb90_meta.rds")

st_dat_sub <- dat %>%
  filter(paste0(record_id, visit) %in% paste0(st_ids$record_id, st_ids$visit)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other"),
         sglt2i_ever = case_when(group == "Lean Control" ~ "No", T~ sglt2i_ever)) %>% 
  dplyr::select(record_id, visit, kit_id, cryostor_id, sex, group, diabetes_duration, 
                age, race, ethnicity, weight, height, bmi, acr_u, hba1c, sglt2i_ever, sglti_timepoint,
                eGFR_CKD_epi, sbp, dbp, map)

# write.csv(st_dat_sub, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Spatial transcriptomics/spatial_transcriptomics_clin.csv", row.names = F)


pb90_dat_sub <- dat %>%
  filter(paste0(record_id, visit) %in% paste0(pb90_ids$record_id, pb90_ids$visit)) %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE")) %>%
  # filter(record_id == "IT_07") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" &
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" &
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other"),
         sglt2i_ever = case_when(group == "Lean Control" ~ "No", T~ sglt2i_ever)) %>%
  dplyr::select(-mrn, -procedure, -dob) %>%
  arrange(record_id)
# %>%
#   rename_with(~ paste0(.x, "_x"), .cols = -c(record_id, visit))


# pb90_ids_unique <- pb90_ids %>%
#   distinct(record_id, visit, .keep_all = T) %>%
#   filter(cohort %in% c("RENAL HEIR", "RENAL HEIRITAGE", "IMPROVE")) %>%
#   arrange(record_id)
# 
test_join <- left_join(pb90_dat_sub, pb90_ids_unique) %>%
  mutate(age_x = floor(age_x)) %>%
  select(record_id, visit, eGFR_CKD_epi_x, eGFR_CKD_epi)

write.csv(pb90_dat_sub, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/scRNA/data_clean/pb90_rhrh2improve_clinical_subset.csv", row.names = F, na =)
