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

all_data <- data.table::fread("ucmr5-occurrence-data/UCMR5_All.txt")

#Classifications 
PFCA_long <- c('PFOA', 'PFNA', 'PFDA','PFUnA', 'PFDoA', 'PFTrDA', 'PFTA', 'ADONA')
PFCA_short <- c('PFBA', 'PFPeA', 'PFHxA', 'PFHpA')
PFCA_ultashort <- c('PFMPA', 'PFMBA')
PFCA_full <- c(PFCA_long, PFCA_short, PFCA_ultashort)

PFSA_long <- c('PFHxS', 'PFHpS', 'PFOS')
PFSA_short <- c('PFBS', 'PFPeS')
PFSA_full <- c(PFSA_long, PFSA_short)

FTS <- c("4:2 FTS", "6:2 FTS", "8:2 FTS")

PFAS_deriv <- 

chrlorinated


PFAS_all 


#Analysis
zip_codes <- data.table::fread("ucmr5-occurrence-data/UCMR5_ZIPCodes.txt")






















