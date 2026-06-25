library(readxl)
library(tidygeocoder)
library(dplyr)
library(purrr)
library(writexl)

# 1. Define the path to your Excel file
excel_file <- "path/to/your/addresses.xlsx"

# 2. Get the names of all four tabs
tab_names <- excel_sheets(excel_file)

# 3. Create a function to read and geocode a single tab
geocode_tab <- function(sheet_name) {
  cat("Processing sheet:", sheet_name, "\n")
  
  # Read the specific sheet
  df <- read_excel(excel_file, sheet = sheet_name)
  
  # Geocode the addresses
  # Note: Replace 'address' with the actual name of your address column
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

# 4. Loop through all tabs and geocode them
# This stores them in a named list of data frames
geocoded_list <- map(tab_names, geocode_tab)
names(geocoded_list) <- tab_names

# 5. Export the results back to a new Excel file with the original tabs
write_xlsx(geocoded_list, "geocoded_studies_output.xlsx")

cat("Geocoding complete! File saved to 'geocoded_studies_output.xlsx'\n")