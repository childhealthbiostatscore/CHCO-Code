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








### Check UCMR-3






### SDOH Data
setwd('/Users/netio/Documents/UofW/Projects/PFAS_Water/')

all_data <- data.table::fread("ucmr3-occurrence-data/UCMR3_All.txt") %>%
  filter(Contaminant != 'lithium')

zip_codes <- data.table::fread("ucmr3-occurrence-data/UCMR3_ZIPCodes.txt")


all_data <- all_data %>% left_join(zip_codes)


# 1. Create wide format table with zip codes and contaminant values
contaminant_wide <- all_data %>%
  select(ZIPCODE, Contaminant, AnalyticalResultValue) %>%
  # Handle non-detects (NA values) - you can replace with 0 or keep as NA
  pivot_wider(
    names_from = Contaminant,
    values_from = AnalyticalResultValue,
    values_fn = mean  # If multiple samples per zip/contaminant, take mean
  ) %>%
  arrange(ZIPCODE)

# View the result
head(contaminant_wide)


IMPROVE <- data.table::fread("Participant_Zips/IMPROVET2D-ZipCodes_DATA_2025-10-20_1056.csv") %>%
  filter(!is.na(mr_number)) %>% 
  dplyr::select(record_id = subject_id, mrn = mr_number, zip_code)

RH <- data.table::fread('Participant_Zips/RENALHEIR-ZipCodes_DATA_2025-10-20_1055.csv') %>% 
  filter(!is.na(mr_number)) %>% 
  dplyr::select(record_id = subject_id, mrn = mr_number, zip_code)

PANTHER <- data.table::fread('Participant_Zips/PANTHER-ZipCodes_DATA_2025-10-20_1232.csv') %>%
  filter(!is.na(mrn)) %>% 
  dplyr::select(record_id, mrn, zip_code)

full_zip <- bind_rows(IMPROVE, RH, PANTHER)
full_zip$zip_code <- as.character(full_zip$zip_code)



contaminant_wide_small <- contaminant_wide %>% filter(ZIPCODE %in% full_zip$zip_code)




########## State-level data 

states <- readxl::read_xlsx("/Users/netio/Downloads/PFAS2020SamplingProject_DrinkingWaterResults.xlsx")

library(zipcodeR)
get_tracts('80301')
all_co_zips <- search_state('CO')


rh <- readxl::read_excel(fs::path(dir_home, 
                                  "drinking_water_pfas", 
                                  "rh_zips.xlsx")) |>
  mutate(zipcode = as.character(zipcode)) |> 
  group_by(zipcode) |> summarize(nparticipants = length(zipcode)) |> ungroup()

# Get city name for each zip code
rh_zips <- tidylog::right_join(rh, all_co_zips)




# Match my zips with list of CO zips
co_zips2 <- co_zips |> 
  tidylog::left_join(rh_zips,
                     by = c("ZCTA5CE10" = "zipcode")) |> 
  mutate(is_in_list = ifelse(!is.na(nparticipants), 
                             "RH participant",
                             "Other Zip Codes"), 
         nparticipants_complete = if_else(is.na(nparticipants), 
                                          0, nparticipants), 
         major_city = toupper(major_city))


# Front Range Counties
frontrange_counties <-  c("Denver County",  
                          "Arapahoe County", 
                          "Jefferson County", 
                          "Adams County", 
                          "Douglas County", 
                          "Broomfield County", 
                          "Elbert County", 
                          # "Park County", 
                          "Gilpin County",
                          "Clear Creek County", 
                          # Colorado Springs Metropolitan Area
                          "El Paso County", 
                          "Teller County", 
                          # Other
                          "Boulder County",
                          # "Eagle County",
                          "Morgan County",
                          "Otero County",
                          "Pueblo County",
                          "Washington County",
                          "Weld County")




participant_counties <- co_zips2 |>
  filter(!is.na(nparticipants))

# Filter only front range counties 
denver_metropolitin_area <- co_zips2 |> 
  filter(county %in% frontrange_counties, 
         zipcode_type != "PO Box") 

# Plot RH Participant Locations
ggplot(data = denver_metropolitin_area |> filter(county == "El Paso County")) + 
  geom_sf(aes(fill = nparticipants_complete), linewidth = .001) + 
  theme_minimal() +
  labs(title = "Participant Locations",
       fill = "Number of Participants")

# Renal-Heir Home prices vs. Colorado Home prices
ggplot(denver_metropolitin_area |> 
         filter(!is.na(median_home_value), 
                zipcode_type != "PO Box"), 
       aes(x = is_in_list, y = median_home_value)) + 
  geom_boxplot()

# Read in Colorado drinking water PFAS results --------------------
co_pfas <- readxl::read_excel(fs::path(dir_home, 
                                       "drinking_water_pfas", 
                                       "PFAS2020SamplingProject_DrinkingWaterResults.xlsx"),
                              col_names = TRUE) |>
  janitor::clean_names() |>
  janitor::remove_constant() |>
  mutate(city_town = if_else(str_detect(loc_desc, " CITY OF") | 
                               str_detect(loc_desc, " TOWN OF"), 
                             str_remove(loc_desc, " CITY OF") |> 
                               str_remove(" TOWN OF"), 
                             "not city water"))


# Filter Detected water systems
co_pfas_ep <- co_pfas |> 
  tidylog::filter(#detect_flag == "Y", 
    loc_type_2 == "ENTRY POINT",  chemical_name == "Total PFAS" 
  ) |> 
  mutate(result_numeric2 = if_else(is.na(result_numeric), 0, result_numeric))


# Summarize by loc_desc 
co_pfas_summary <- co_pfas_ep |> 
  group_by(loc_desc, chemical_name) |> 
  summarise(mean_pfas = mean(result_numeric, na.rm = TRUE) |>
              replace_na(0), 
            n_samples = length(loc_desc),  
            city_town = city_town[1]) |>
  ungroup()

# # Filter not city water
# not_city_water <- co_pfas_ep |> 
#   tidylog::filter(city_town == "not city water")


# Merge participant locations with water districts
participant_pfas <- tidylog::left_join(participant_counties, co_pfas_summary, 
                                       by = c("major_city" = "city_town"))



# Make plot of PK estimated PFAS
# using steady-state ratios implied by the serum PFAS calculator PK models (Lu and Bartell, 2020) 
pfas = c("PFOA", "PFOS", "PFHxS", "PFNA")
S = c(118.35, 129.31, 201.57, 200.67)   
B = c(1.67, 5.20, 1.30, 0.60)   # in ng/mL
W = seq(0,40,10)   # in ng/L


Cinf = S %o% W / 1000 + B           # from Bartell 2017

# Convert to data frame and bind the W vector
df <- as.data.frame(Cinf)
rownames(df) <- pfas
colnames(df) <- W

df <- df |> as_tibble(rownames = "pfas")

# Reshape from wide to long format
df_long <- df %>%
  pivot_longer(cols = -"pfas", 
               names_to = "W", 
               values_to = "concentration") |>
  mutate(W = as.numeric(W))

# 1a: pk modeled PFAS
(p1 <- ggplot(df_long, aes(x=W, y=concentration, color=pfas, group = pfas)) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = 4, linetype = 2, color = "red", linewidth = 1.2) +
    scale_color_brewer(type = "qual", palette = 6) +
    labs(x = "Water PFAS concentration (ng/L)", 
         y = "Expected serum\nPFAS (ng/mL)") +
    xlim(c(0, 40)) +
    theme(legend.position=c(.65, .25), 
          legend.title = element_blank(),
          axis.text.x = element_blank(), 
          axis.title.x = element_blank()))

# 1b: drinking water PFAS
(p2 <- ggplot(participant_pfas |> filter(!is.na(mean_pfas)), 
              aes(x=mean_pfas)) +
    geom_vline(xintercept = 4, linetype = 2, color = "red", linewidth = 1.2) +
    geom_histogram(aes(y = after_stat(count / sum(count))),
                   binwidth = 10, 
                   color = "black", 
                   fill = "grey80", 
                   breaks=seq(0,40,10)) + 
    stat_bin(breaks=seq(0,40,10),
             binwidth = 10, geom = "text", color = "black",
             aes(y = after_stat(count / sum(count)), 
                 label = scales::percent(after_stat(count / sum(count)))),
             position = position_stack(vjust = 0.5)) + 
    xlim(c(0, 40)) +
    scale_y_continuous(labels = scales::percent) + 
    labs(x = "Water PFAS levels (ng/L)", 
         y = "Renal-HEIR\nParticipants (%)"))


# Combine plots 
combined_plot <- cowplot::plot_grid(NULL, p1, 
                                    NULL, p2, 
                                    ncol = 1, 
                                    align = "v", 
                                    rel_heights = c(.05, .5, .05, .5))

ggsave(filename = fs::path(dir_figure, "Figure PK modeled PFAS.jpeg"), 
       height = 5, width = 4)






















