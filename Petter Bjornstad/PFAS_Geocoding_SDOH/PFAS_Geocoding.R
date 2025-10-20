### Geocoding PFAS, SDOH Analysis 


##libraries
library(dplyr)
library(stringr)
library(tidyverse)
library(purrr)
library(ggplot2)


library(tidygeocoder)
library(dplyr)

# Single address
#address <- "1600 Amphitheatre Parkway, Mountain View, CA"
#geo_result <- geo(address, method = "osm")  # OpenStreetMap/Nominatim - completely free

# Multiple addresses in a dataframe
#addresses_df <- data.frame(
#  id = 1:3,
#  address = c("University of Washington, Seattle, WA",
#              "Space Needle, Seattle, WA", 
#              "Pike Place Market, Seattle, WA")
#)

# Geocode all addresses
#geocoded <- addresses_df %>%
#  geocode(address, method = "osm", lat = latitude, long = longitude)

# maybe use ArcGIS? method = "arcgis" (ArcGIS - free, no API key needed)



#geocoded_zips <- zip_df %>%
#  geocode(zip, method = "osm")

# Or format as full address for better results
#zip_df_formatted <- zip_df %>%
#  mutate(address = paste(zip, "USA")) %>%
#  geocode(address, method = "osm")






### PFAS Data

# Download from: https://www.epa.gov/ground-water-and-drinking-water/safe-drinking-water-information-system-sdwis-federal

# Or access via API
#library(httr)
#library(jsonlite)

# Example API call for a specific PWS
#pws_id <- "WA1234567"
#url <- paste0("https://data.epa.gov/efservice/GEOGRAPHIC_AREA/PWSID/", pws_id, "/JSON")
#response <- GET(url)
#data <- fromJSON(content(response, "text"))



### SDOH Data
setwd('/Users/netio/Documents/UofW/Projects/PFAS_Water/')

all_data <- data.table::fread("ucmr5-occurrence-data/UCMR5_All.txt") %>%
  filter(Contaminant != 'lithium')

#Classifications 
PFCA_long <- c('PFOA', 'PFNA', 'PFDA','PFUnA', 'PFDoA', 'PFTrDA', 'PFTA', 'ADONA')
PFCA_short <- c('PFBA', 'PFPeA', 'PFHxA', 'PFHpA')
PFCA_ultashort <- c('PFMPA', 'PFMBA')
PFCA_full <- c(PFCA_long, PFCA_short, PFCA_ultashort)

PFSA_long <- c('PFHxS', 'PFHpS', 'PFOS')
PFSA_short <- c('PFBS', 'PFPeS')
PFSA_full <- c(PFSA_long, PFSA_short)

FTS <- c("4:2 FTS", "6:2 FTS", "8:2 FTS")

PFAS_deriv <- c('NEtFOSAA', 'NMeFOSAA', 'HFPO-DA', 'ADONA', 'PFEESA', 'NFDHA')

chlorinated <- c("9Cl-PF3ONS", "11Cl-PF3OUdS")


PFAS_all <- c(PFCA_full, PFSA_full, FTS, PFAS_deriv, chlorinated)


zip_codes <- data.table::fread("ucmr5-occurrence-data/UCMR5_ZIPCodes.txt")



#participant data

IMPROVE <- data.table::fread("Participant_Zips/IMPROVET2D-ZipCodes_DATA_2025-10-20_1056.csv") %>%
  filter(!is.na(mr_number)) %>% 
  dplyr::select(record_id = subject_id, mrn = mr_number, zip_code)

RH <- data.table::fread('Participant_Zips/RENALHEIR-ZipCodes_DATA_2025-10-20_1055.csv') %>% 
  filter(!is.na(mr_number)) %>% 
  dplyr::select(record_id = subject_id, mrn = mr_number, zip_code)

full_zip <- bind_rows(IMPROVE, RH)
full_zip$zip_code <- as.character(full_zip$zip_code)


## Zip Code Plotting 

library(tigris)
library(ggplot2)
library(dplyr)
library(sf)
library(data.table)

# Count participants per ZIP code
zip_counts <- full_zip %>%
  group_by(zip_code) %>%
  summarise(n_participants = n()) %>%
  ungroup()

# Get ALL US ZIP codes (no state filter for recent years)
all_zips <- zctas(year = 2020, cb = TRUE)  # cb = TRUE for simplified boundaries

# Join your data
map_data <- all_zips %>%
  left_join(zip_counts, by = c("ZCTA5CE20" = "zip_code"))

# Filter to only show ZIPs with participants (makes map cleaner and faster)
map_data_with_participants <- map_data %>%
  filter(!is.na(n_participants))

# Get state boundaries
states_sf <- states(year = 2020, cb = TRUE) %>%
  filter(STUSPS %in% c("CO", "NM", 'WY', 'NE'))  # Filter after downloading

# Map with state outlines
ggplot() +
  geom_sf(data = map_data_with_participants, 
          aes(fill = n_participants), 
          color = "white", 
          size = 0.1) +
  geom_sf(data = states_sf, 
          fill = NA, 
          color = "black", 
          size = 0.8) +  # State boundaries
  scale_fill_gradient(low = "lightblue", high = "darkblue", 
                      name = "Participants") +
  theme_void() +
  labs(title = "Study Participants by ZIP Code")



## harmonized data 
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

#dat <- harmonized_data %>% dplyr::select(-dob) %>% 
#  arrange(date_of_screen) %>% 
#  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
#                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
#                   .by = c(record_id, visit))



dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))



















