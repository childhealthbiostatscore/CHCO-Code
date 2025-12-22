# load libraries
library(tidyverse)
library(dplyr)
library(REDCapR)
library(readxl)
library(jsonlite)
library(Hmisc)
# Set up environment for Kopah
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else {
  stop("Unknown user: please specify root path for this user.")
}

dat <- read.csv(file.path(root_path, "/Data Harmonization/Data Clean/harmonized_dataset.csv"), na.strings = c(" ", "", "-9999",-9999))

# read in biopsy master spreadsheet (updated real time)
biopsy_master <- read_excel("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/Biopsy master tracker.xlsx", skip = 2) 
colnames(biopsy_master) <- gsub(" ", "_", colnames(biopsy_master))  # Replace spaces with underscores
colnames(biopsy_master) <- gsub("/", "_", colnames(biopsy_master))  # Replace slashes with underscores
colnames(biopsy_master) <- gsub("\\(", "_", colnames(biopsy_master))  # Replace opening parentheses with underscores
colnames(biopsy_master) <- gsub("\\)", "", colnames(biopsy_master))  # Remove closing parentheses
biopsy_master <- biopsy_master %>%
  dplyr::rename(record_id =Study_ID ,
                visit = Visit_ID) %>%
  mutate(visit = case_when(Study == "ATTEMPT" ~ visit,
                           Study == "REMODEL" ~ visit, 
                           Study == "IMPROVE" & visit == "4M" ~ "3_months_post_surgery",
                           visit == "12M" ~ "12_months_post_surgery",
                           visit == "Follow-up" ~ "follow-up",
                           T ~ "baseline"),
         record_id = case_when(startsWith(record_id, "RPC") & !is.na(Coenroll_ID__Same_visit_only) ~ Coenroll_ID__Same_visit_only,
                               T ~ record_id)) # RH2-60-T and RH2-48-T coenrolled into RPC2 as of 02/03/25 data pull

# subset biopsy master to PANDA
biopsy_master_panda <- biopsy_master %>%
  # filter(Study == "PANDA") %>%
  filter(!is.na(record_id))

# subset harm dat to PANDA
dat_panda <- dat %>%
  filter(procedure == "kidney_biopsy") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(mrn, date)) %>%
  filter(!is.na(panda_id)) %>%
  filter(study != "ATTEMPT")
  
# record ID in master spreadsheet but not in REDCap
biopsy_master_panda$record_id[biopsy_master_panda$record_id %nin% dat_panda$record_id]

# record ID in REDCap but not master spreadsheet
dat_panda$record_id[dat_panda$record_id %nin% biopsy_master_panda$record_id & dat_panda$record_id %nin% biopsy_master_panda$Coenroll_ID__Same_visit_only]
