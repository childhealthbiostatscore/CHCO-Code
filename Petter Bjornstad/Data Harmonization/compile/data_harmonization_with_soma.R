# =====================================================================
# DATA HARMONIZATION - SINGLE PIPELINE (clinical + proteomics)
# This is the one script to run to rebuild the harmonized datasets.
# It merges the Python clinical harmonizer and the R proteomics steps
# into one end to end pipeline.
#
# WHAT IT DOES (in order):
#   1. Sources the Python clinical engine (compile/data_harmonization.py)
#      via reticulate and runs harmonize_data(), which pulls every REDCap
#      study, merges them, and returns a temp CSV of the clinical data.
#   2. Adds ATTEMPT data from Antoine (ATTEMPT_AC.RData).
#   3. Computes BMI percentiles and z scores with growthcleanr.
#   4. Merges SomaScan proteomics (soma_combined_anml_2.RData).
#   5. Merges Olink plasma and urine proteomics.
#
# OUTPUTS (to Data Harmonization/Data Clean/):
#   harmonized_dataset.csv               clinical only, with BMI percentiles
#   soma_harmonized_dataset.csv          clinical + SomaScan
#   olink_plasma_harmonized_dataset.csv  clinical + Olink plasma
#   olink_urine_harmonized_dataset.csv   clinical + Olink urine
#   soma_olink_harmonized_dataset.csv    clinical + Olink + SomaScan
#
# REQUIRES: reticulate, dplyr, tidyr, purrr, growthcleanr, lubridate, and
#   the Python clinical engine plus its inputs (REDCap tokens, etc.).
#   The Python engine stays in its own file because it does the heavy
#   REDCap pulls; this R script is the single entry point that drives it.
#
# NOTE: this script writes harmonized_dataset.csv directly. The Python
#   script compile/data harmonization save.py also writes that file and
#   uploads to S3. Decide which one is the source of truth so they do not
#   overwrite each other (see INDEX.md attention notes).
# =====================================================================

library(reticulate)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(togolab)
library(uuid)
library(growthcleanr)
library(readr)
library(Hmisc)
library(aws.s3)
library(powerjoin)
togolab::togo_paths()

# Import python harmonization function & run
source_python(file.path(git_path, 'Data Harmonization/compile/data_harmonization.py'))
temp_path <- harmonize_data()
clean <- read.csv(temp_path, na.strings = c("", "NaN"), check.names = FALSE)
clean <- data.frame(lapply(clean, as.character), check.names = FALSE)
clean[clean == "NaN"] <- NA # Replace NaN from Python to NA
clean[clean == ""] <- NA

# ATTEMPT from Antoine
s3load(object = "attempt/ATTEMPT_AC.RData", bucket = "raw.data", region = "")
names(merged_data)[names(merged_data) %in% names(clean)]
attempt <- merged_data %>%
  mutate(group = "Type 1 Diabetes",
         study = "ATTEMPT",
         hba1c = coalesce(hba1c, hba1c_percent),
         diabetes_duration = coalesce(diabetes_duration, diabetes_dx_duration),
         creatinine_s = creatinine_serum_umoll / 88.42) %>%
  dplyr::select(any_of(names(clean)), record_id, visit, study, af_pwv = PWV, 
                -diabetes_dx_duration,
                treatment_arm, -hba1c_percent, -creatinine_serum_umoll,
                -bmi_z, -bmi_percentile, 
                everything()) %>%
  filter(visit %nin% c("Unscheduled 1", "Unscheduled 2"))
attempt[] <- lapply(attempt, function(x) {
  if (is.numeric(x)) as.character(x) else x
})
attempt$date <- mdy(attempt$date)

# Create screen_date (screening date for each participant or earliest date available)
keys <- c("record_id", "visit", "procedure")

clean <- clean %>%
  dplyr::mutate(date = as.Date(date)) %>%
  power_full_join(
    attempt %>% select(all_of(keys), !any_of(setdiff(names(clean), keys))),
    by = keys,
    conflict = coalesce
  ) %>%
  dplyr::group_by(record_id, visit) %>%
  dplyr::mutate(screen_date = case_when(procedure == "screening" | visit == "screening" ~ date),
                hba1c = coalesce(hba1c, prescreen_a1c),
                hct = coalesce(hct, hematocrit)) %>%
  ungroup() %>% dplyr::group_by(record_id) %>%
  fill(screen_date, .direction = "updown") %>%
  dplyr::mutate(screen_date = coalesce(screen_date, 
                                       if (any(!is.na(date))) min(date, na.rm = TRUE) else as.Date(NA))) %>%
  fill(screen_date, .direction = "updown") %>% ungroup() %>%
  group_by(mrn) %>%
  fill(dob, screen_date,
       attempt_id, casper_id, coffee_id, croc_id, improve_id, penguin_id,
       rh_id, rh2_id, panther_id, panda_id, rpc2_id, swht_id, ultra_id, co_enroll_id,
       .direction = "downup") %>%
  ungroup() %>%
  dplyr::select(record_id, attempt_id, casper_id, coffee_id, croc_id, improve_id, penguin_id,
                rh_id, rh2_id, panther_id, panda_id, rpc2_id, swht_id, ultra_id, co_enroll_id,
                mrn, date, screen_date, everything(),
                -prescreen_a1c, -sphyg_sex)

# Calculate BMI Percentiles using growthcleanr
# Convert date columns safely
clean <- clean %>%
  mutate(
    dob = suppressWarnings(parse_date_time(dob, orders = c("ymd", "mdy", "dmy"))),
    date = suppressWarnings(parse_date_time(date, orders = c("ymd", "mdy", "dmy"))),
    age_mo = coalesce(as.numeric(difftime(date, dob, units = "days")) / 30.44, as.numeric(age)*12),
    age_mo = as.numeric(age_mo),
    weight = as.numeric(weight),
    height = as.numeric(height),
    bmi = as.numeric(bmi)
  )

bmi_input <- clean %>%
  filter(!is.na(sex) & !is.na(age_mo) & !is.na(weight) & !is.na(height) & !is.na(bmi))

# Run growthcleanr BMI z-score and percentile calculation
bmi_percentile <- growthcleanr::ext_bmiz(
  data = subset(bmi_input, select = c(record_id, date, procedure, visit, sex, age_mo, weight, height, bmi)),
  age = "age_mo",
  wt = "weight",
  ht = "height",
  bmi = "bmi",
  adjust.integer.age = FALSE
) %>%
  dplyr::select(record_id, date, procedure, visit, bmip, bmiz) %>%
  filter(!is.na(bmip))

# Merge results into harmonized dataset
clean <- clean %>%
  left_join(bmi_percentile, by = c("record_id", "visit", "date", "procedure"))

# Map to UUID
generate_f_uuid <- function() {
  u <- gsub("-", "", UUIDgenerate())  # 32-char hex, no hyphens
  paste0("f", substring(u, 2))        # swap first char to 'f' for fake
}

# Load existing map if present, otherwise start empty
map_path <- file.path(root_path, "PHI_data/harmonized dataset/mrn_uuid_map.csv")
if (file.exists(map_path)) {
  map_df <- read_csv(map_path, col_types = cols(mrn = col_character(),
                                                uuid = col_character()))
  map_df$mrn <- sub("\\.0+$", "", map_df$mrn)
  mrn_uuid_map <- setNames(map_df$uuid, map_df$mrn)
} else {
  mrn_uuid_map <- character(0)
}

# Unique non-missing MRNs from the harmonized data
mrns <- unique(clean$mrn[!is.na(clean$mrn)])

# Only mint UUIDs for MRNs not already mapped
new_mrns <- setdiff(mrns, names(mrn_uuid_map))

for (m in new_mrns) {
  mrn_uuid_map[m] <- generate_f_uuid()
}

map_df <- tibble(mrn = names(mrn_uuid_map), uuid = unname(mrn_uuid_map))
cat(sprintf("Total: %d mappings (%d newly added)\n",
            nrow(map_df), length(new_mrns)))

# write_csv(map_df, map_path) # check the number of new mappings before saving a new map

clean <- clean %>%
  left_join(map_df) %>%
  dplyr::select(-mrn) # remove mrn from our harmonized datasets

# Save clinical harmonized dataset without mrn
s3write_using_region(x = clean, 
                     FUN = write.csv,
                     object = "harmonized dataset/harmonized_dataset.csv",
                     bucket = "raw.data",
                     row.names = F,
                     na = "",
                     region = "")

cat("Added BMI percentiles and z-scores to harmonized dataset and saved to CSV.\n")

##########################################################################################

# Combine with SOMA
s3load(object = "SomaLogic/soma_combined_anml_2.RData", bucket = "raw.data", region = "")

soma <- soma_combined %>%
  as.matrix() %>% as.data.frame() %>%
  mutate_at(vars(starts_with("seq")), as.numeric) %>%
  dplyr::rename("record_id" = SampleDescription) %>%
  dplyr::mutate(record_id = gsub("IT2D-", "IT_", record_id),
                record_id = gsub("PAN_", "PAN-", record_id), # PAN-44-O entered as PAN_44-O in soma
                record_id = gsub("PAN-14-O", "PAN-14-C", record_id), # PAN-14-O to PAN-14-C (group changed, name should have been changed but changed later)
                record_id = gsub("PAN-103-0", "PAN-103-O", record_id), # PAN-103-O entered as PAN-103-0 in soma
                visit = case_when(grepl("PAN-", record_id) ~ SampleGroup,
                                  T ~ TimePoint),
                visit = case_when(visit == "Baseline" ~ "baseline",
                                  visit == "BL" ~ "baseline",
                                  visit == "Year 1" ~ "year_1",
                                  visit == "Year 2" ~ "year_2",
                                  visit == "Year 3" ~ "year_3",
                                  visit == "3M" ~ "3_months_post_surgery",
                                  visit == "12M" ~ "12_months_post_surgery",
                                  visit == "V1" ~ "baseline",
                                  visit == "V4" ~ "4_months_post"),
                visit = case_when(grepl("PAN-", record_id) & visit == "year_1" ~ "baseline", # PANTHER dates skips baseline in SOMA (shifted)
                                  grepl("PAN-", record_id) & visit == "year_2" ~ "year_1",
                                  grepl("PAN-", record_id) & visit == "year_3" ~ "year_2",
                                  T ~ visit)) %>%
  dplyr::select(record_id, visit, starts_with("seq"), urine_creat_proteomics) %>% 
  mutate_at(vars(starts_with("seq")), as.numeric)

# Check for mismatches between SOMA harmonized and data harmonized (have proteomics but no clinical data)
mismatches <- soma$record_id[!soma$record_id %in% clean$record_id]

# Combine SOMA & data harmonized
soma_clean <- clean %>% 
  full_join(soma, by = c("record_id", "visit")) %>%
  group_by(uuid, screen_date) %>%
  {
    has_uuid_screen <- filter(., !is.na(uuid) & !is.na(screen_date)) %>%
      fill(starts_with("seq"), .direction = "downup")
    bind_rows(
      has_uuid_screen,
      filter(., is.na(uuid) | is.na(screen_date))
    )
  } %>%
  ungroup()

# Save harmonized data with SOMA
s3write_using_region(x = soma_clean,
                     FUN = write.csv,
                     object = "harmonized dataset/soma_harmonized_dataset.csv",
                     bucket = "raw.data",
                     row.names = F,
                     na = "",
                     region = "")

##########################################################################################

# Combine with Olink
# Olink data
## Plasma
olink_plasma <- s3read_using_region(FUN = read.csv,
                                    object = "Olink/plasma_cleaned.csv",
                                    bucket = "raw.data",
                                    region = "")

olink_plasma <- olink_plasma %>%
  mutate(visit = case_when(endsWith(record_id, "BL") ~ "baseline",
                           endsWith(record_id, "12M") ~ "12_months_post_surgery",
                           T ~ "baseline"),
         record_id = gsub("_BL|_12M", "", record_id)) %>%
  rename_with(~ paste0(.x, "_p"), starts_with("OID"))

olink_plasma_clean <- clean %>% left_join(olink_plasma, by = c("record_id", "visit")) %>%
  group_by(uuid, screen_date) %>%
  {
    has_uuid_screen <- filter(., !is.na(uuid) & !is.na(screen_date)) %>%
      fill(starts_with("OID"), .direction = "downup")
    bind_rows(
      has_uuid_screen,
      filter(., is.na(uuid) | is.na(screen_date))
    )
  } %>%
  ungroup()

# Save harmonized data with Olink Plasma
s3write_using_region(x = olink_plasma_clean,
                     FUN = write.csv,
                     object = "harmonized dataset/olink_plasma_harmonized_dataset.csv",
                     bucket = "raw.data",
                     row.names = F,
                     na = "", 
                     region = "")

## Urine
olink_urine <- s3read_using_region(FUN = read.csv,
                                   object = "Olink/urine_cleaned.csv",
                                   bucket = "raw.data",
                                   region = "")
olink_urine <- olink_urine %>%
  mutate(visit = case_when(endsWith(record_id, "BL") ~ "baseline",
                           endsWith(record_id, "12M") ~ "12_months_post_surgery",
                           T ~ "baseline"),
         record_id = gsub("_BL|_12M", "", record_id)) %>%
  rename_with(~ paste0(.x, "_u"), starts_with("OID"))

olink_urine_clean <- clean %>% left_join(olink_urine, by = c("record_id", "visit")) %>%
  group_by(uuid, screen_date) %>%
  fill(starts_with("OID"), .direction = "downup") %>%
  ungroup() 

# Save harmonized data with Olink Urine
s3write_using_region(x = olink_urine_clean,
                     FUN = write.csv,
                     object = "harmonized dataset/olink_urine_harmonized_dataset.csv",
                     bucket = "raw.data",
                     row.names = F,
                     na = "", 
                     region = "")

# Save harmonized data with Olink (both plasma and urine) and SOMA
clean_proteomics_comb <- clean %>% 
  left_join(olink_plasma, by = c("record_id", "visit")) %>%
  left_join(olink_urine, by = c("record_id", "visit")) %>%
  left_join(soma, by = c("record_id", "visit")) %>%
  group_by(uuid, screen_date) %>%
  {
    has_uuid_screen <- filter(., !is.na(uuid) & !is.na(screen_date)) %>%
      fill(starts_with("OID"), .direction = "downup") %>%
      fill(starts_with("seq"), .direction = "downup")
    bind_rows(
      has_uuid_screen,
      filter(., is.na(uuid) | is.na(screen_date))
    )
  } %>%
  ungroup()

s3write_using_region(x = clean_proteomics_comb,
                     FUN = write.csv,
                     object = "harmonized dataset/soma_olink_harmonized_dataset.csv",
                     bucket = "raw.data",
                     row.names = F,
                     na = "", 
                     region = "")
