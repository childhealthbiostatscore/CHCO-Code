library(readxl)
library(tidygeocoder)
library(dplyr)
library(purrr)
library(writexl)

# 1. Define the path to your Excel file
excel_file <- "C:/Users/netio/Downloads/AddressData_CHCO-2.xlsx"

# 2. Get the names of all four tabs
tab_names <- excel_sheets(excel_file)

# 3. Create a function to read and geocode a single tab
geocode_tab <- function(sheet_name) {
  cat("Processing sheet:", sheet_name, "\n")
  
  # Read the specific sheet
  df <- read_excel(excel_file, sheet = sheet_name, col_names=F)
  
  df <- df %>% 
    rename(
      study_id = ...1,
      address  = ...2
    )
  
  
  # Geocode the addresses
  df_geocoded <- df %>%
    geocode(
      address = address, 
      method = 'osm',      # 'osm' uses OpenStreetMap (free, no API key needed)
      lat = latitude,      
      long = longitude,    
      unique_only = TRUE   
    )
  
  return(df_geocoded)
}

geocoded_list <- map(tab_names, geocode_tab)
names(geocoded_list) <- tab_names

# 5. Export the results back to a new Excel file with the original tabs
write_xlsx(geocoded_list, "C:/Users/netio/Documents/UofW/Projects/PFAS_Water/geocoded_studies_output.xlsx")

cat("Geocoding complete! File saved to 'geocoded_studies_output.xlsx'\n")


