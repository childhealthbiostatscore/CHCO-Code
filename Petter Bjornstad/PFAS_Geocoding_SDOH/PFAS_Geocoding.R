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

save.image('PWS_Variables_withZipCodes_Line295.RData')


full_crosswalk_list_overlap <- full_crosswalk_list %>% filter(MAIL_ZIP %in% full_zip$zip_code) %>% 
  filter(!duplicated(MAIL_ZIP, Name))



pfas_drinking <- readxl::read_excel('PFAS2020SamplingProject_DrinkingWaterResults.xlsx')






library(tidyverse)

# Get all unique PWS IDs from your PFAS data
all_pws <- unique(pfas_drinking$LOC_NAME)

cat("Looking up", length(all_pws), "unique PWS systems...\n\n")

# Match all systems
all_matched_systems <- water_systems %>%
  filter(pwsid %in% all_pws) %>%
  select(pwsid, pws_name, city_name, zip_code, address_line1, state_code, 
         county_code = state_code, population_served_count)  # Include useful columns

cat("=== MATCH RESULTS FOR ALL SYSTEMS ===\n")
cat("Matched:", nrow(all_matched_systems), "out of", length(all_pws), "systems\n")
cat("Match rate:", round(100 * nrow(all_matched_systems) / length(all_pws), 1), "%\n\n")

# Check zip code availability in matched systems
cat("Systems with zip codes:", sum(!is.na(all_matched_systems$zip_code)), 
    "out of", nrow(all_matched_systems), "\n")
cat("Zip code availability:", 
    round(100 * sum(!is.na(all_matched_systems$zip_code)) / nrow(all_matched_systems), 1), "%\n\n")

# Save the complete lookup table
write_csv(all_matched_systems, "PWS_complete_lookup.csv")
cat("✓ Complete lookup table saved to PWS_complete_lookup.csv\n\n")

# Now join with your PFAS data
pfas_with_zips <- pfas_drinking %>%
  left_join(
    all_matched_systems %>% select(pwsid, pws_name, city_name, zip_code),
    by = c("LOC_NAME" = "pwsid")
  )

# Final statistics
final_stats <- pfas_with_zips %>%
  summarise(
    total_rows = n(),
    rows_with_zip = sum(!is.na(zip_code)),
    match_rate = round(100 * rows_with_zip / total_rows, 1)
  )

cat("=== FINAL PFAS DATA WITH ZIP CODES ===\n")
print(final_stats)

# Check which systems didn't match
unmatched_systems <- pfas_drinking %>%
  filter(!LOC_NAME %in% all_matched_systems$pwsid) %>%
  distinct(LOC_NAME, LOC_DESC) %>%
  arrange(LOC_NAME)

if (nrow(unmatched_systems) > 0) {
  cat("\n⚠ Unmatched systems (", nrow(unmatched_systems), "):\n", sep = "")
  print(unmatched_systems)
  write_csv(unmatched_systems, "PWS_unmatched_systems.csv")
  cat("\nUnmatched systems saved to PWS_unmatched_systems.csv\n")
} else {
  cat("\n✓ All systems matched!\n")
}

# Save the final PFAS data with zip codes
write_csv(pfas_with_zips, "pfas_drinking_with_zipcodes_COMPLETE.csv")
cat("\n✓ Final dataset saved to pfas_drinking_with_zipcodes_COMPLETE.csv\n")

# Show a sample of the final data
cat("\nSample of final data:\n")
print(pfas_with_zips %>% 
        select(LOC_NAME, LOC_DESC, CHEMICAL_NAME, RESULT_NUMERIC, 
               city_name, zip_code) %>%
        head(20))


library(tidyverse)

# Pivot wider to get one column per chemical
pfas_wide <- pfas_with_zips %>%
  select(LOC_NAME, LOC_DESC, CHEMICAL_NAME, RESULT_NUMERIC, 
         pws_name, city_name, zip_code) %>%
  pivot_wider(
    names_from = CHEMICAL_NAME,
    values_from = RESULT_NUMERIC,
    values_fn = mean  # In case there are duplicates, take the mean
  )

# Check the result
cat("Original data dimensions:", nrow(pfas_with_zips), "rows x", 
    ncol(pfas_with_zips), "columns\n")
cat("Wide data dimensions:", nrow(pfas_wide), "rows x", 
    ncol(pfas_wide), "columns\n\n")

# Show the chemical columns
chemical_cols <- names(pfas_wide)[!names(pfas_wide) %in% 
                                    c("LOC_NAME", "LOC_DESC", "pws_name", 
                                      "city_name", "zip_code")]
cat("Chemical columns (", length(chemical_cols), "):\n", sep = "")
print(chemical_cols)

# Preview the data
cat("\nPreview of wide format:\n")
print(pfas_wide %>% head(10))

# Save the wide format
write_csv(pfas_wide, "pfas_drinking_wide_format.csv")
cat("\n✓ Wide format saved to pfas_drinking_wide_format.csv\n")

# Optional: Clean up column names (remove spaces, special characters)
pfas_wide_clean <- pfas_wide %>%
  rename_with(
    ~str_replace_all(., "[^[:alnum:]_]", "_") %>%  # Replace special chars with _
      str_replace_all(., "_+", "_") %>%             # Remove multiple underscores
      str_to_lower(),                               # Make lowercase
    .cols = all_of(chemical_cols)
  )

cat("\nCleaned column names:\n")
print(names(pfas_wide_clean))

pfas_wide_clean$zip_code <- str_replace(pfas_wide_clean$zip_code, '-\\d+', '')

write_csv(pfas_wide_clean, "pfas_drinking_wide_format_clean.csv")
cat("\n✓ Clean wide format saved to pfas_drinking_wide_format_clean.csv\n")


## Analysis 


full_zip_analysis <- full_zip %>% left_join(pfas_wide_clean, by='zip_code')


library(tidyverse)
library(data.table)

# Assuming your data is called 'pfas_data' (adjust if different)
# If it's a data.table, convert for easier tidyverse use
pfas_df <- as_tibble(full_zip_analysis)

# 1. Get the PFAS chemical columns (exclude ID and location columns)
pfas_cols <- names(pfas_df)[!names(pfas_df) %in% 
                              c("record_id", "mrn", "zip_code", "LOC_NAME", 
                                "LOC_DESC", "pws_name", "city_name")]

cat("=== PFAS CHEMICALS ANALYZED ===\n")
cat("Total chemicals:", length(pfas_cols), "\n\n")

# 2. Count non-NA (detected) values for each chemical
detection_summary <- pfas_df %>%
  summarise(across(
    all_of(pfas_cols),
    list(
      n_detected = ~sum(!is.na(.)),
      n_participants = ~n(),
      detection_rate = ~round(100 * sum(!is.na(.)) / n(), 1),
      mean_detected = ~mean(., na.rm = TRUE),
      median_detected = ~median(., na.rm = TRUE),
      max_detected = ~max(., na.rm = TRUE)
    )
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("chemical", ".value"),
    names_pattern = "(.+)_(n_detected|n_participants|detection_rate|mean_detected|median_detected|max_detected)"
  ) %>%
  arrange(desc(detection_rate))

cat("=== DETECTION SUMMARY BY CHEMICAL ===\n")
print(detection_summary, n = 25)

write_csv(detection_summary, "pfas_detection_summary.csv")
cat("\n✓ Saved to pfas_detection_summary.csv\n\n")

# 3. Count unique zip codes with data
zip_summary <- pfas_df %>%
  summarise(
    total_participants = n(),
    participants_with_zip = sum(!is.na(zip_code)),
    unique_zips = n_distinct(zip_code, na.rm = TRUE)
  )

cat("=== ZIP CODE SUMMARY ===\n")
print(zip_summary)
cat("\n")

# List of zip codes represented
zip_list <- pfas_df %>%
  filter(!is.na(zip_code)) %>%
  count(zip_code, sort = TRUE)

cat("Zip codes represented:\n")
print(zip_list, n = 20)
write_csv(zip_list, "zip_codes_represented.csv")

# 4. Calculate average per participant (in case there are duplicates)
pfas_participant_avg <- pfas_df %>%
  group_by(record_id, mrn, zip_code) %>%
  summarise(
    across(all_of(pfas_cols), ~mean(., na.rm = TRUE)),
    .groups = "drop"
  )

cat("\n=== PARTICIPANT AVERAGES ===\n")
cat("Unique participants:", nrow(pfas_participant_avg), "\n\n")

write_csv(pfas_participant_avg, "pfas_participant_averages.csv")
cat("✓ Saved to pfas_participant_averages.csv\n\n")

# 5. Create distribution plots for each PFAS
library(ggplot2)
library(patchwork)

# Prepare data for plotting (long format, only detected values)
pfas_long <- pfas_participant_avg %>%
  select(record_id, all_of(pfas_cols)) %>%
  pivot_longer(
    cols = all_of(pfas_cols),
    names_to = "chemical",
    values_to = "concentration"
  ) %>%
  filter(!is.na(concentration)) %>%
  # Clean up chemical names for plotting
  mutate(chemical_clean = str_replace_all(chemical, "_", " ") %>%
           str_to_title() %>%
           str_trunc(40))

# Get top 10 most detected chemicals for plotting
top_chemicals <- detection_summary %>%
  slice_max(n_detected, n = 10) %>%
  pull(chemical)

# Plot 1: Histogram of top chemicals
p1 <- pfas_long %>%
  filter(chemical %in% top_chemicals) %>%
  ggplot(aes(x = concentration)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~chemical_clean, scales = "free", ncol = 3) +
  labs(
    title = "Distribution of PFAS Concentrations",
    subtitle = "Top 10 Most Detected Chemicals",
    x = "Concentration",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    plot.title = element_text(face = "bold")
  )

ggsave("pfas_distributions_histogram.png", p1, width = 12, height = 10, dpi = 300)
cat("✓ Saved pfas_distributions_histogram.png\n")

# Plot 2: Box plots of top chemicals
p2 <- pfas_long %>%
  filter(chemical %in% top_chemicals) %>%
  ggplot(aes(x = reorder(chemical_clean, concentration, median), 
             y = concentration)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(alpha = 0.2, width = 0.2, size = 0.5) +
  coord_flip() +
  labs(
    title = "PFAS Concentration Distributions",
    subtitle = "Box plots showing median and quartiles",
    x = "Chemical",
    y = "Concentration"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("pfas_distributions_boxplot.png", p2, width = 10, height = 8, dpi = 300)
cat("✓ Saved pfas_distributions_boxplot.png\n")

# Plot 3: Detection rates
p3 <- detection_summary %>%
  slice_max(detection_rate, n = 15) %>%
  mutate(chemical_clean = str_replace_all(chemical, "_", " ") %>%
           str_to_title() %>%
           str_trunc(40)) %>%
  ggplot(aes(x = reorder(chemical_clean, detection_rate), 
             y = detection_rate)) +
  geom_col(fill = "coral", alpha = 0.8) +
  geom_text(aes(label = paste0(detection_rate, "%")), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  ylim(0, 110) +
  labs(
    title = "PFAS Detection Rates",
    subtitle = "Percentage of participants with detected levels",
    x = "Chemical",
    y = "Detection Rate (%)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("pfas_detection_rates.png", p3, width = 10, height = 8, dpi = 300)
cat("✓ Saved pfas_detection_rates.png\n")

# Plot 4: Log-scale distributions for chemicals with wide ranges
p4 <- pfas_long %>%
  filter(chemical %in% top_chemicals) %>%
  ggplot(aes(x = concentration)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  facet_wrap(~chemical_clean, scales = "free", ncol = 3) +
  scale_x_log10() +
  labs(
    title = "Distribution of PFAS Concentrations (Log Scale)",
    subtitle = "Top 10 Most Detected Chemicals",
    x = "Concentration (log scale)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    plot.title = element_text(face = "bold")
  )

ggsave("pfas_distributions_log_scale.png", p4, width = 12, height = 10, dpi = 300)
cat("✓ Saved pfas_distributions_log_scale.png\n\n")

cat("=== SUMMARY ===\n")
cat("✓ Generated 4 plots\n")
cat("✓ Created detection summary table\n")
cat("✓ Calculated participant averages\n")
cat("✓ Identified unique zip codes\n")



# Continue from previous code...

# 6. Zip code representation and detection analysis

# Count participants per zip code
zip_participants <- pfas_participant_avg %>%
  filter(!is.na(zip_code)) %>%
  count(zip_code, name = "n_participants") %>%
  arrange(desc(n_participants))

cat("=== ZIP CODE PARTICIPATION ===\n")
print(zip_participants, n = 20)

# Calculate detection rates by zip code for top chemicals
zip_detection <- pfas_participant_avg %>%
  filter(!is.na(zip_code)) %>%
  group_by(zip_code) %>%
  summarise(
    n_participants = n(),
    across(
      all_of(top_chemicals),
      list(
        n_detected = ~sum(!is.na(.)),
        detection_rate = ~round(100 * sum(!is.na(.)) / n(), 1),
        mean_conc = ~mean(., na.rm = TRUE)
      )
    ),
    .groups = "drop"
  )

write_csv(zip_detection, "pfas_detection_by_zipcode.csv")
cat("✓ Saved pfas_detection_by_zipcode.csv\n\n")
# Modified Plot 5: Number of participants per zip code (MORE ZIP CODES)
p5 <- zip_participants %>%
  slice_max(n_participants, n = 50) %>%  # ← Increased from 20 to 50
  ggplot(aes(x = reorder(zip_code, n_participants), y = n_participants)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = n_participants), hjust = -0.2, size = 2.5) +
  coord_flip() +
  labs(
    title = "Participants by Zip Code",
    subtitle = "Top 50 zip codes with most participants",
    x = "Zip Code",
    y = "Number of Participants"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7)  # Smaller text for more zip codes
  )

ggsave("participants_by_zipcode.png", p5, width = 10, height = 14, dpi = 300)
cat("✓ Saved participants_by_zipcode.png\n")

# Modified Plot 6: Detection Y/N by zip code (MORE ZIP CODES, BINARY DETECTION)
top_zips_expanded <- zip_participants %>%
  slice_max(n_participants, n = 50) %>%  # ← Increased from 15 to 50
  pull(zip_code)

zip_detection_binary <- zip_detection %>%
  filter(zip_code %in% top_zips_expanded) %>%
  select(zip_code, n_participants, ends_with("_detection_rate")) %>%
  pivot_longer(
    cols = ends_with("_detection_rate"),
    names_to = "chemical",
    values_to = "detection_rate"
  ) %>%
  mutate(
    chemical = str_remove(chemical, "_detection_rate") %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_trunc(30),
    detected = ifelse(detection_rate > 0, "Detected", "Not Detected")  # ← Binary Y/N
  )

p6 <- zip_detection_binary %>%
  ggplot(aes(x = chemical, y = zip_code, fill = detected)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("Detected" = "red", "Not Detected" = "lightblue"),
    name = "Detection Status"
  ) +
  labs(
    title = "PFAS Detection by Zip Code",
    subtitle = "Top 50 zip codes with most participants (Y/N detection)",
    x = "Chemical",
    y = "Zip Code"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 6),  # Smaller text for more zip codes
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("detection_rates_by_zipcode_heatmap.png", p6, width = 12, height = 16, dpi = 300)
cat("✓ Saved detection_rates_by_zipcode_heatmap.png\n")




#### Participants 


# Fixed study classification based on record_id
study_participants <- full_zip_analysis %>%
  filter(!is.na(city_name)) %>%
  filter(!duplicated(mrn)) %>% 
  mutate(
    study = case_when(
      str_detect(record_id, "^IT[-_]") ~ "IMPROVE",  # ← Fixed: IT- OR IT_
      str_detect(record_id, "^RH2[-_]") ~ "Renal-HEIRitage",  # ← Fixed: RH2- OR RH2_
      str_detect(record_id, "^RH[-_]") ~ "Renal-HEIR",  # ← Fixed: RH- OR RH_
      str_detect(record_id, "^PAN[-_]") ~ "PANTHER",  # ← Fixed: PAN- OR PAN_
      str_detect(record_id, "^\\d+$") ~ "CROCODILE",  # Just numbers
      TRUE ~ "Other"  # Catch any unexpected formats
    )
  )

# Count participants by study
study_counts <- study_participants %>%
  count(study, name = "n_participants") %>%
  arrange(desc(n_participants)) %>%
  mutate(
    percentage = round(100 * n_participants / sum(n_participants), 1),
    label = paste0(study, "\n", n_participants, " (", percentage, "%)")
  )

cat("=== PARTICIPANTS BY STUDY (with city_name) ===\n")
print(study_counts)
cat("\nTotal participants with city_name:", sum(study_counts$n_participants), "\n\n")

# Verify IT_ records are now classified correctly
cat("=== SAMPLE IT RECORDS ===\n")
study_participants %>%
  filter(str_detect(record_id, "^IT")) %>%
  select(record_id, study) %>%
  head(10) %>%
  print()

# Create pie chart
p_pie <- ggplot(study_counts, aes(x = "", y = n_participants, fill = study)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5),
            size = 4,
            fontface = "bold") +
  scale_fill_brewer(palette = "Set2", name = "Study") +
  labs(
    title = "Participant Distribution by Study",
    subtitle = "Participants with city name data"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 10)
  )

ggsave("participants_by_study_pie.png", p_pie, width = 10, height = 8, dpi = 300)
cat("✓ Saved participants_by_study_pie.png\n")

# Bar chart version
p_bar <- ggplot(study_counts, aes(x = reorder(study, n_participants), 
                                  y = n_participants, fill = study)) +
  geom_col(alpha = 0.8, show.legend = FALSE) +
  geom_text(aes(label = paste0(n_participants, "\n(", percentage, "%)")), 
            hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Participant Distribution by Study",
    subtitle = "Participants with city name data",
    x = "Study",
    y = "Number of Participants"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11)
  )

ggsave("participants_by_study_bar.png", p_bar, width = 10, height = 6, dpi = 300)
cat("✓ Saved participants_by_study_bar.png\n")

# Save the study classification data
write_csv(study_participants, "participants_with_study_classification.csv")
cat("\n✓ Saved participants_with_study_classification.csv\n")


