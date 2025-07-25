# Hi Ye Ji! I'm pulling some samples from CROC and Renal Heir/Heirtage for analysis for Petter 
# and we were wondering if you could send some REDCap lists for the participants who had PI 
# samples collected. Petter was telling me we had 50 people enrolled in CROC but I have 123 
# samples so it seems like I might have some duplicates. Additionally, I heard there were 4 CROC 
# participants who also enrolled in IMPROVE, would you be able to pull their CROC/IT2D participant 
# numbers so I can ensure they are in the analysis? Thank you!

library(reticulate)
library(dplyr)
library(tidyr)
library(purrr)
library(REDCapR)

tokens <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive//Data Harmonization/api_tokens.csv")
tokens <- tokens %>%
  filter(Study %in% c("CROCODILE", "Renal-HEIR", "Renal-HEIRitage"))
redcap_uri <- "https://redcap.ucdenver.edu/api/"

# Define a named vector of tokens
study_tokens <- setNames(tokens$Token, tokens$Study) 

# Read all REDCap elements (data, metadata, users) for each study
data <- map(names(study_tokens), function(study) {
  token <- study_tokens[[study]]
  
  message("Processing study: ", study)
  
  # Read data
  data <- redcap_read(
    redcap_uri = redcap_uri,
    token = token
  )$data

})

rh <- data[[1]]
croc <- data[[2]]
rh2 <- data[[3]]

rh2_pi <- rh2 %>%
  filter(!is.na(pi_copeptin)) %>%
  dplyr::select(record_id, pi_copeptin)

croc_pi <- croc %>%
  filter(!is.na(pi_copeptin)) %>%
  mutate(record_id = case_when(record_id < 10 ~ paste0("CRC-0", record_id),
                               T ~ paste0("CRC-", record_id))) %>%
  dplyr::select(record_id, pi_copeptin)

pi_comp <- rbind(rh2_pi, croc_pi)

harm_dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

rh_rh2_crc <- harm_dat %>% filter(study %in% c("CROCODILE", "RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE")) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>%
  select(record_id, mrn, ends_with("_id")) 

pi_comp_id <- pi_comp %>% left_join(rh_rh2_crc) %>%
  select(record_id, casper_id, coffee_id, croc_id, improve_id, penguin_id, rh_id, rh2_id, panther_id, panda_id)
  
write.csv(pi_comp_id, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Exports/pi_ids_croc_rh_rh2.csv", row.names = F, na = "")
