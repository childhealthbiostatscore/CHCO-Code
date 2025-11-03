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

#all_data <- data.table::fread("ucmr5-occurrence-data/UCMR5_All.txt") %>%
#  filter(Contaminant != 'lithium')

#Classifications 
#PFCA_long <- c('PFOA', 'PFNA', 'PFDA','PFUnA', 'PFDoA', 'PFTrDA', 'PFTA', 'ADONA')
#PFCA_short <- c('PFBA', 'PFPeA', 'PFHxA', 'PFHpA')
#PFCA_ultashort <- c('PFMPA', 'PFMBA')
#PFCA_full <- c(PFCA_long, PFCA_short, PFCA_ultashort)

#PFSA_long <- c('PFHxS', 'PFHpS', 'PFOS')
#PFSA_short <- c('PFBS', 'PFPeS')
#PFSA_full <- c(PFSA_long, PFSA_short)

#FTS <- c("4:2 FTS", "6:2 FTS", "8:2 FTS")

#PFAS_deriv <- c('NEtFOSAA', 'NMeFOSAA', 'HFPO-DA', 'ADONA', 'PFEESA', 'NFDHA')

#chlorinated <- c("9Cl-PF3ONS", "11Cl-PF3OUdS")


#PFAS_all <- c(PFCA_full, PFSA_full, FTS, PFAS_deriv, chlorinated)


#zip_codes <- data.table::fread("ucmr5-occurrence-data/UCMR5_ZIPCodes.txt")


#all_data <- all_data %>% left_join(zip_codes)


# 1. Create wide format table with zip codes and contaminant values
#contaminant_wide <- all_data %>%
 # select(ZIPCODE, Contaminant, AnalyticalResultValue) %>%
  # Handle non-detects (NA values) - you can replace with 0 or keep as NA
 # pivot_wider(
 #   names_from = Contaminant,
 #   values_from = AnalyticalResultValue,
 #   values_fn = mean  # If multiple samples per zip/contaminant, take mean
 # ) %>%
 # arrange(ZIPCODE)

# View the result
#head(contaminant_wide)

# 2. Create MRL reference table (one row per contaminant)
#mrl_table <- all_data %>%
#  select(Contaminant, MRL, Units) %>%
#  distinct() %>%
 # arrange(Contaminant)

# View the MRL table
#print(mrl_table)

#remove(all_data)

#participant data

IMPROVE <- data.table::fread("Participant_Zips/IMPROVET2D-ZipCodes_DATA_2025-10-20_1056.csv") %>%
  filter(!is.na(mr_number)) %>% 
  dplyr::select(record_id = subject_id, mrn = mr_number, zip_code)

RH <- data.table::fread('Participant_Zips/RENALHEIR-ZipCodes_DATA_2025-10-20_1055.csv') %>% 
  filter(!is.na(mr_number)) %>% 
  dplyr::select(record_id = subject_id, mrn = mr_number, zip_code)

PANTHER <- data.table::fread('Participant_Zips/PANTHER-ZipCodes_DATA_2025-10-20_1232.csv') %>%
  filter(!is.na(mrn)) %>% 
  dplyr::select(record_id, mrn, zip_code)

RH2 <- readxl::read_excel('Participant_Zips/Final_Renal_Croc_Zip.xlsx') %>%
  dplyr::select(record_id = `Subject ID`, zip_code = `ZIP Code`) %>%
  mutate(mrn = NA)


full_zip <- bind_rows(IMPROVE, RH, PANTHER, RH2)
full_zip$zip_code <- as.character(full_zip$zip_code)



#contaminant_wide_small <- contaminant_wide %>% filter(ZIPCODE %in% full_zip$zip_code)








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
  filter(STUSPS %in% c("CO", "NM", 'WY', 'NE', 'WA'))  # Filter after downloading

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








### Check Colorado data 
library(sf)
co_data <- readxl::read_excel("PFAS2020SamplingProject_DrinkingWaterResults.xlsx")
co_data <- co_data %>% 
  dplyr::select(LOC_NAME, CHEMICAL_NAME, DETECT_FLAG, RESULT_NUMERIC, REPORT_RESULT_LIMIT, METHOD_DETECTION_LIMIT, 
                REPORTING_DETECTION_LIMIT)


pwd_areas <- st_read('CWS_Boundaries_Latest/co-municipal-water-provider-boundaries.geojson')


library(sf)
library(tigris) # Assuming you use tigris for ZCTAs
library(dplyr)

# 1. Load your PWS data (assuming it's loaded as 'pwd_areas')
# pwd_areas <- st_read("path/to/your/PWS_file.geojson") 

# 2. Load ZCTA data using tigris
options(tigris_use_cache = TRUE) # Cache ZCTAs to avoid re-downloading
all_zips <- zctas(year = 2020, cb = TRUE) 

# 3. Make both datasets valid to resolve Topology Exceptions
pwd_areas_valid <- st_make_valid(pwd_areas)
all_zips_valid <- st_make_valid(all_zips)

# 4. Turn off s2 spherical geometry processing (as you did)
sf_use_s2(FALSE)
# "although coordinates are longitude/latitude, st_intersection assumes that they are planar" 
# This message is normal when S2 is FALSE, indicating planar math will be used.

# 5. Ensure both valid datasets have the same Coordinate Reference System (CRS)
pwd_areas_valid <- st_transform(pwd_areas_valid, st_crs(all_zips_valid))

# 6. Perform the intersection with the valid data
# Use the *valid* object names here
pws_zcta_overlap <- st_intersection(pwd_areas_valid, all_zips_valid)



library(dplyr)
library(sf)

# The 'pws_zcta_overlap' object created in the previous steps is what we use here.

# Extract all PWS variables and the corresponding zip code
full_crosswalk_list <- pws_zcta_overlap %>%
  as.data.frame() %>% 
  # Select all columns *except* the geometry columns that we don't need in a flat list
  select(
    MAIL_ZIP, # This is the 5-digit zip code
    everything() # This keeps all other columns from both datasets
  ) %>%
  # Drop the spatial geometry column
  st_drop_geometry() %>%
  # Keep only unique rows (important to de-duplicate based on the combination of 
  # PWS attributes and the zip code it overlaps with)
  distinct() %>%
  # Arrange for better readability
  arrange(MAIL_ZIP, Name) # Replace with the actual PWS ID column name

# View the first few rows and all columns
head(full_crosswalk_list)
names(full_crosswalk_list) # Check all the variable names included

# Save the complete list to a CSV file
write.table(full_crosswalk_list, "PWS_Variables_with_ZipCodes_List.txt", row.names = FALSE, quote=F, sep ='\t')





