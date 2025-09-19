library(dplyr)
library(ggplot2)
library(stringr)


missing_data <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/IDs_forHuntingNewData.txt')



harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))

dat <- dat %>% filter(mrn %in% missing_data$mrn) %>% 
  dplyr::select(record_id, mrn, visit, date, study, gbm_thick_artmean, gbm_thick_harmmean, m_i, acprg, airg, avg_c_k2)


dups_vec <- which(duplicated(dat$mrn))











