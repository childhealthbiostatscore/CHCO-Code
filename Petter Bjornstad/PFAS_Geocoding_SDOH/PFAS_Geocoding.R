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
address <- "1600 Amphitheatre Parkway, Mountain View, CA"
geo_result <- geo(address, method = "osm")  # OpenStreetMap/Nominatim - completely free

# Multiple addresses in a dataframe
addresses_df <- data.frame(
  id = 1:3,
  address = c("University of Washington, Seattle, WA",
              "Space Needle, Seattle, WA", 
              "Pike Place Market, Seattle, WA")
)

# Geocode all addresses
geocoded <- addresses_df %>%
  geocode(address, method = "osm", lat = latitude, long = longitude)

# maybe use ArcGIS? method = "arcgis" (ArcGIS - free, no API key needed)


































