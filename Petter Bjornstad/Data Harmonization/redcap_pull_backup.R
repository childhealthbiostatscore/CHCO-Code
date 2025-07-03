library(reticulate)
library(dplyr)
library(tidyr)
library(purrr)
library(REDCapR)

identifiers <- c("last_name", "first_name", "mr_number", "dob",
                 "name_last", "name_first", "mrn",
                 "phone_number", "phone")

tokens <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive//Data Harmonization/api_tokens.csv")

redcap_uri <- "https://redcap.ucdenver.edu/api/"

# Define a named vector of tokens
study_tokens <- setNames(tokens$Token, tokens$Study)

# Read all REDCap elements (data, metadata, users) for each study
REDCap_backup <- map(names(study_tokens), function(study) {
  token <- study_tokens[[study]]
  
  message("Processing study: ", study)
  
  # Read data
  data <- redcap_read(
    redcap_uri = redcap_uri,
    token = token
  )$data %>%
    dplyr::select(-any_of(identifiers))
  
  # Read metadata
  metadata <- redcap_metadata_read(
    redcap_uri = redcap_uri,
    token = token,
    verbose = FALSE
  )$data
  
  # Read user rights
  users <- redcap_users_export(
    redcap_uri = redcap_uri,
    token = token
  )$data_user
  
  # Read recent log entries (last 30 days)
  log <- redcap_log_read(
    redcap_uri = redcap_uri,
    token = token,
    log_begin_date = Sys.Date() - 30L,
    log_end_date = Sys.Date(),
    verbose = FALSE
  )$data
  
  list(
    data = data,
    metadata = metadata,
    users = users,
    log = log
  )
})

# Name the list elements by study
names(REDCap_backup) <- names(study_tokens)

# Save the back-up as RData
backup_date <- format(Sys.Date(), "%m%d%y")

save(REDCap_backup,
     file = paste0("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/Back Up/REDCap_backup_",
                   backup_date, ".RData"))
save(REDCap_backup,
     file = paste0("/Volumes/Expansion/REDCap/REDCap_backup_",
                   backup_date, ".RData"))
