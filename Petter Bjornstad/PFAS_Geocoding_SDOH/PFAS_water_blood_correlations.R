###PFAS Correlations



library(dplyr)
library(stringr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(tidygeocoder)



imputed_pfas <- readRDS("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/PANTHER/Data_Cleaned/PFAS_Data_Imputed_11_03.rds")

base.dir <- 'C:/Users/netio/Documents/UofW/Projects/PFAS_Water/'

water_data <- data.table::fread(paste0(base.dir, 'participants_with_study_classification.csv'))

harmonized_path <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"
harmonized_data <- read.csv(harmonized_path, na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))





























