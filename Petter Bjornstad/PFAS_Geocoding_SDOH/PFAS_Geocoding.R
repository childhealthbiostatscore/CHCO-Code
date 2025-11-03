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

co_data <- readxl::read_excel("PFAS2020SamplingProject_DrinkingWaterResults.xlsx")
co_data <- co_data %>% 
  dplyr::select(LOC_NAME, CHEMICAL_NAME, DETECT_FLAG, RESULT_NUMERIC, REPORT_RESULT_LIMIT, METHOD_DETECTION_LIMIT, 
                REPORTING_DETECTION_LIMIT)

library(httr)
library(jsonlite)
library(tidyverse)

# ============================================================
# METHOD 1: EPA SDWIS API (Primary - Most Reliable)
# ============================================================

get_pws_from_sdwis <- function(pwsid) {
  base_url <- "https://data.epa.gov/efservice/WATER_SYSTEM"
  
  url <- paste0(base_url, "/PWSID/", pwsid, "/JSON")
  
  cat("Trying SDWIS API for:", pwsid, "\n")
  
  tryCatch({
    response <- GET(url)
    
    if (status_code(response) == 200) {
      data <- fromJSON(content(response, "text", encoding = "UTF-8"))
      
      if (length(data) > 0) {
        return(tibble(
          PWSID = pwsid,
          PWS_Name = data$PWS_NAME[1] %||% NA,
          City = data$CITY_NAME[1] %||% NA,
          State = data$STATE_CODE[1] %||% NA,
          Zip = data$ZIP_CODE[1] %||% NA,
          County = data$COUNTY_SERVED[1] %||% NA,
          Population = data$POPULATION_SERVED_COUNT[1] %||% NA,
          Method = "SDWIS_API"
        ))
      }
    }
    
    cat("  SDWIS API: No results\n")
    return(NULL)
    
  }, error = function(e) {
    cat("  SDWIS API error:", e$message, "\n")
    return(NULL)
  })
}

# ============================================================
# METHOD 2: EPA ECHO API (Alternative endpoint)
# ============================================================

get_pws_from_echo <- function(pwsid) {
  # Try the search endpoint instead
  base_url <- "https://echodata.epa.gov/echo/sdw_rest_services.get_qid"
  
  url <- paste0(
    base_url,
    "?output=JSON",
    "&qcolumns=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15",
    "&p_pwsid=", pwsid
  )
  
  cat("Trying ECHO API for:", pwsid, "\n")
  
  tryCatch({
    response <- GET(url)
    
    if (status_code(response) == 200) {
      data <- fromJSON(content(response, "text", encoding = "UTF-8"))
      
      # Get the query ID
      if (!is.null(data$Results$QueryID)) {
        qid <- data$Results$QueryID
        
        # Now fetch the actual facilities
        facilities_url <- paste0(
          "https://echodata.epa.gov/echo/sdw_rest_services.get_facilities",
          "?output=JSON&qid=", qid
        )
        
        Sys.sleep(0.5)
        fac_response <- GET(facilities_url)
        
        if (status_code(fac_response) == 200) {
          fac_data <- fromJSON(content(fac_response, "text", encoding = "UTF-8"))
          
          if (!is.null(fac_data$Results$Facilities) && nrow(fac_data$Results$Facilities) > 0) {
            facility <- fac_data$Results$Facilities[1, ]
            
            return(tibble(
              PWSID = pwsid,
              PWS_Name = facility$PWSName %||% NA,
              City = facility$CityName %||% NA,
              State = facility$StateCode %||% NA,
              Zip = facility$ZipCode %||% NA,
              County = facility$CountyName %||% NA,
              Population = facility$PopulationServedCount %||% NA,
              Method = "ECHO_API"
            ))
          }
        }
      }
    }
    
    cat("  ECHO API: No results\n")
    return(NULL)
    
  }, error = function(e) {
    cat("  ECHO API error:", e$message, "\n")
    return(NULL)
  })
}

# ============================================================
# METHOD 3: Web Scraping EPA Drinking Water Watch
# ============================================================

get_pws_from_web <- function(pwsid) {
  # Clean PWS ID (remove dashes if present)
  clean_id <- gsub("-", "", pwsid)
  
  url <- paste0(
    "https://enviro.epa.gov/enviro/sdw_query_v3.get_list",
    "?pws_id=", clean_id,
    "&state=CO"
  )
  
  cat("Trying web scraping for:", pwsid, "\n")
  
  tryCatch({
    response <- GET(url)
    
    if (status_code(response) == 200) {
      content_text <- content(response, "text")
      
      # This returns HTML - would need rvest to parse
      # For now, return NULL and we'll use manual lookup
      cat("  Web method needs rvest package for parsing\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("  Web scraping error:", e$message, "\n")
    return(NULL)
  })
}

# ============================================================
# MAIN FUNCTION: Try all methods
# ============================================================

get_pws_details <- function(pwsid) {
  # Try Method 1: SDWIS API
  result <- get_pws_from_sdwis(pwsid)
  if (!is.null(result)) return(result)
  
  Sys.sleep(1)
  
  # Try Method 2: ECHO API
  result <- get_pws_from_echo(pwsid)
  if (!is.null(result)) return(result)
  
  # If all fail, return empty row
  cat("  All methods failed for:", pwsid, "\n\n")
  return(tibble(
    PWSID = pwsid,
    PWS_Name = NA,
    City = NA,
    State = NA,
    Zip = NA,
    County = NA,
    Population = NA,
    Method = "FAILED"
  ))
}

# ============================================================
# MANUAL LOOKUP TABLE (Colorado Water Systems)
# ============================================================

# Based on public records, here are common Colorado systems
manual_lookup <- tribble(
  ~PWSID, ~PWS_Name, ~City, ~State, ~Zip, ~County,
  "CO0101015", "BELLEVIEW CHRISTIAN COLLEGE", "LAKEWOOD", "CO", "80226", "JEFFERSON",
  "CO0101025", "BRIGHTON CITY OF", "BRIGHTON", "CO", "80601", "ADAMS"
)

# Function to use manual lookup first
get_pws_with_fallback <- function(pwsid) {
  # First check manual lookup
  manual_match <- manual_lookup %>% filter(PWSID == pwsid)
  
  if (nrow(manual_match) > 0) {
    cat("Found in manual lookup:", pwsid, "\n")
    return(manual_match %>% mutate(Population = NA, Method = "MANUAL"))
  }
  
  # Otherwise try APIs
  return(get_pws_details(pwsid))
}

# ============================================================
# EXAMPLE USAGE
# ============================================================


cat("Fetching PWS details with fallback methods...\n\n")
pws_details <- map_dfr(example_systems, get_pws_with_fallback)

print(pws_details)




