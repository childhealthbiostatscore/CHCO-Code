###Brain Imaging by CKD Staging 


# Load libraries
library(oro.nifti)  
library(neurobase)  
library(R.matlab)    
library(ggplot2)     
library(dplyr)       
library(corrplot)   
library(psych)        
library(pls)  
library(caret) 
library(randomForest) 
library(igraph)
library(brainGraph)
library(R.matlab)
library(tidyverse)

bucket <- 'brain.mri'


#Identifying groups for analysis 

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))



mri_ids <- c('CRC-10', 'CRC-11', 'CRC-12', 'CRC-13', 'CRC-26', 'CRC-39', 'CRC-46', 'CRC-51', 'CRC-53', 
             'CRC-55', 'CRC-58', 'CRC-60', 
             'RH2-01-O', 'RH2-03-O', 'RH2-08-T', 'RH2-10-L', 'RH2-11-O', 'RH2-13-O', 'RH2-16-O', 'RH2-17-L', 
             'RH2-18-O', 'RH2-19-T', 'RH2-22-T', 'RH2-24-L', 'RH2-27-L', 'RH2-28-L', 'RH2-29-L', 'RH2-33-L', 
             'RH2-34-O', 'RH2-35-T', 'RH2-38-T', 'RH2-39-O', 'RH2-41-T', 'RH2-42-T', 'RH2-43-T', 
             'RH2-44-T', 'RH2-45-T', 'RH2-48-T', 'RH2-49-T', 'RH2-50-L', 'RH2-52-T', 'RH2-53-T', 
             'RH2-55-T')

mri_ids_df <- data.frame(ID = mri_ids, file_id = c(
  'crc_10', 'crc_11', 'crc_12', 'crc_13', 'crc_26', 'crc_39', 'crc_46', 'crc_51', 'crc_53', 
  'crc_55', 'crc_58', 'crc_60', 
  'RH2_01_O', 'RH2_03_O', 'RH2_08_T', 'RH2_10_L', 'RH2_11_O', 'RH2_13_O', 'RH2_16_O', 'RH2_17_L', 
  'RH2_18_O', 'RH2_19_T', 'RH2_22_T', 'RH2_24_L', 'RH2_27_L', 'RH2_28_L', 'RH2_29_L', 'RH2_33_L', 
  'RH2_34_O', 'RH2_35_T', 'RH2_38_T', 'RH2_39_O', 'RH2_41_T', 'RH2_42_T', 'RH2_43_T', 
  'RH2_44_T', 'RH2_45_T', 'RH2_48_T', 'RH2_49_T', 'RH2_50_L', 'RH2_52_T', 'RH2_53_T', 
  'RH2_55_T'
))


small_dat <- dat %>% 
  filter(record_id %in% mri_ids)

#find missing
#mri_ids[which(!mri_ids %in% small_dat$record_id)]

#missing_dat <- dat %>% filter(rh2_id == 'RH2-38-O')


small_dat$group[which(small_dat$record_id == 'RH2-38-T')] <- 'Obese Control'
small_dat$record_id[which(small_dat$record_id == 'RH2-38-T')] <- 'RH2-38-O'





qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")











classify_ckd <- function(gfr = NULL, albumin_mg_g = NULL, albumin_mg_mmol = NULL, 
                         gfr_days_from_baseline = NULL, albumin_days_from_baseline = NULL, 
                         max_days_apart = 365) {
  
  # Convert albumin units if needed
  if (!is.null(albumin_mg_mmol) && is.null(albumin_mg_g)) {
    # Convert mg/mmol to mg/g (multiply by 8.84)
    albumin_mg_g <- albumin_mg_mmol * 8.84
  }
  
  # Handle missing values - allow partial classification
  gfr_available <- !is.null(gfr) && !is.na(gfr)
  albumin_available <- !is.null(albumin_mg_g) && !is.na(albumin_mg_g)
  
  # Check if days from baseline are provided and within acceptable range
  days_available <- !is.null(gfr_days_from_baseline) && !is.null(albumin_days_from_baseline) && 
    !is.na(gfr_days_from_baseline) && !is.na(albumin_days_from_baseline)
  
  within_time_window <- TRUE
  days_apart <- NA
  criteria_met_days_from_baseline <- NA
  
  if (days_available && gfr_available && albumin_available) {
    # Calculate days between tests
    days_apart <- abs(gfr_days_from_baseline - albumin_days_from_baseline)
    within_time_window <- days_apart <= max_days_apart
    
    # Determine when criteria was met (earlier of the two time points)
    if (within_time_window) {
      criteria_met_days_from_baseline <- min(gfr_days_from_baseline, albumin_days_from_baseline)
    }
  } else if ((days_available && (gfr_available || albumin_available)) || 
             (!days_available && gfr_available && albumin_available)) {
    # Handle partial days or no days with values
    
    # For partial classification, use the available days from baseline
    if (gfr_available && !albumin_available && !is.null(gfr_days_from_baseline) && !is.na(gfr_days_from_baseline)) {
      criteria_met_days_from_baseline <- gfr_days_from_baseline
    } else if (!gfr_available && albumin_available && !is.null(albumin_days_from_baseline) && !is.na(albumin_days_from_baseline)) {
      criteria_met_days_from_baseline <- albumin_days_from_baseline
    }
  }
  
  # If both are missing, return all NAs
  if (!gfr_available && !albumin_available) {
    return(list(
      gfr_category = NA,
      gfr_description = NA,
      albumin_category = NA,
      albumin_description = NA,
      risk_level = NA,
      recommendation = NA,
      days_between_tests = NA,
      within_time_window = NA,
      criteria_met_days_from_baseline = NA
    ))
  }
  
  # Classify GFR category if available
  if (gfr_available) {
    gfr_category <- case_when(
      gfr >= 90 ~ "G1",
      gfr >= 60 & gfr < 90 ~ "G2",
      gfr >= 45 & gfr < 60 ~ "G3a",
      gfr >= 30 & gfr < 45 ~ "G3b",
      gfr >= 15 & gfr < 30 ~ "G4",
      gfr < 15 ~ "G5",
      TRUE ~ NA_character_
    )
  } else {
    gfr_category <- NA_character_
  }
  
  # Classify Albuminuria category if available
  if (albumin_available) {
    albumin_category <- case_when(
      albumin_mg_g < 30 ~ "A1",
      albumin_mg_g >= 30 & albumin_mg_g < 300 ~ "A2",
      albumin_mg_g >= 300 ~ "A3",
      TRUE ~ NA_character_
    )
  } else {
    albumin_category <- NA_character_
  }
  
  # Determine risk level and recommendation only if both values are available
  # AND within the time window (or if no days provided)
  if (gfr_available && albumin_available && within_time_window) {
    risk_recommendation <- get_risk_recommendation(gfr_category, albumin_category)
  } else {
    risk_recommendation <- list(risk = NA, recommendation = NA)
  }
  
  return(list(
    gfr_category = gfr_category,
    gfr_description = get_gfr_description(gfr_category),
    albumin_category = albumin_category,
    albumin_description = get_albumin_description(albumin_category),
    risk_level = risk_recommendation$risk,
    recommendation = risk_recommendation$recommendation,
    days_between_tests = if(!is.na(days_apart)) round(days_apart) else NA,
    within_time_window = within_time_window,
    criteria_met_days_from_baseline = if(!is.na(criteria_met_days_from_baseline)) round(criteria_met_days_from_baseline) else NA
  ))
}


get_risk_recommendation <- function(gfr_cat, alb_cat) {
  
  # Create the risk matrix based on the KDIGO guidelines
  risk_matrix <- matrix(
    c(
      # A1         A2              A3
      "Low",      "Moderate",     "High",        # G1
      "Low",      "Moderate",     "High",        # G2
      "Moderate", "High",         "Very high",   # G3a
      "High",     "Very high",    "Very high",   # G3b
      "Very high","Very high",    "Very high",   # G4
      "Very high","Very high",    "Very high"    # G5
    ),
    nrow = 6, ncol = 3, byrow = TRUE,
    dimnames = list(
      c("G1", "G2", "G3a", "G3b", "G4", "G5"),
      c("A1", "A2", "A3")
    )
  )
  
  # Create the recommendation matrix
  recommendation_matrix <- matrix(
    c(
      # A1                    A2                      A3
      "Screen",              "Treat",                "Treat and refer",      # G1
      "Screen",              "Treat",                "Treat and refer",      # G2
      "Treat",               "Treat",                "Treat and refer",      # G3a
      "Treat",               "Treat and refer",      "Treat and refer",      # G3b
      "Treat and refer",     "Treat and refer",      "Treat and refer",      # G4
      "Treat and refer 4+",  "Treat and refer 4+",   "Treat and refer 4+"   # G5
    ),
    nrow = 6, ncol = 3, byrow = TRUE,
    dimnames = list(
      c("G1", "G2", "G3a", "G3b", "G4", "G5"),
      c("A1", "A2", "A3")
    )
  )
  
  # Get risk and recommendation
  if (!is.na(gfr_cat) && !is.na(alb_cat)) {
    risk <- risk_matrix[gfr_cat, alb_cat]
    recommendation <- recommendation_matrix[gfr_cat, alb_cat]
  } else {
    risk <- NA
    recommendation <- NA
  }
  
  return(list(risk = risk, recommendation = recommendation))
}

# Helper function to get GFR category description
get_gfr_description <- function(category) {
  descriptions <- c(
    "G1" = "Normal or high",
    "G2" = "Mildly decreased",
    "G3a" = "Mildly to moderately decreased",
    "G3b" = "Moderately to severely decreased",
    "G4" = "Severely decreased",
    "G5" = "Kidney failure"
  )
  return(descriptions[category])
}

# Helper function to get Albumin category description
get_albumin_description <- function(category) {
  descriptions <- c(
    "A1" = "Normal to mildly increased (<30 mg/g)",
    "A2" = "Moderately increased (30-299 mg/g)",
    "A3" = "Severely increased (≥300 mg/g)"
  )
  return(descriptions[category])
}





##Data for analysis 


# Apply the CKD classification to your MRI subset
small_dat <- small_dat %>%
  rowwise() %>%
  mutate(
    ckd_classification = list(classify_ckd(
      gfr = eGFR_CKD_epi,
      albumin_mg_g = acr_u,
      gfr_days_from_baseline = NULL,
      albumin_days_from_baseline = NULL
    ))
  ) %>%
  ungroup() %>%
  # Unnest the results into separate columns
  mutate(
    gfr_category = sapply(ckd_classification, function(x) x$gfr_category),
    gfr_description = sapply(ckd_classification, function(x) x$gfr_description),
    albumin_category = sapply(ckd_classification, function(x) x$albumin_category),
    albumin_description = sapply(ckd_classification, function(x) x$albumin_description),
    risk_level = sapply(ckd_classification, function(x) x$risk_level),
    recommendation = sapply(ckd_classification, function(x) x$recommendation)
  ) %>%
  select(-ckd_classification)

# Check the results
small_dat %>%
  select(record_id, eGFR_CKD_epi, acr_u, gfr_category, albumin_category, risk_level) %>%
  print(n = 20)

# Summary of CKD stages
table(small_dat$gfr_category, small_dat$albumin_category, useNA = "always")








