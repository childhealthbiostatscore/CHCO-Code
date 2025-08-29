# KDIGO CKD Staging System Implementation in R
# Optimized for Long Format Data
# Based on KDIGO 2012 Clinical Practice Guideline

# Core function to classify CKD stage based on GFR and Albuminuria
# Accepts days from baseline instead of actual dates
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

# Helper function to determine risk level and recommendation
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

# Main function to classify CKD from long format data
classify_ckd_long <- function(data,
                              patient_id_col = "releaseid",
                              days_col = "days",
                              gfr_col = "ckd_gfr",
                              albumin_col = "UAlbCreat_mg_g",
                              max_days_apart = 365) {
  
  library(dplyr)
  library(tidyr)
  
  # Convert columns to standard names for processing
  data_std <- data %>%
    rename(
      patient_id = !!sym(patient_id_col),
      days = !!sym(days_col),
      gfr = !!sym(gfr_col),
      albumin_mg_g = !!sym(albumin_col)
    )
  
  # Option 1: Check for same-day measurements first (both values in same row)
  same_day_measurements <- data_std %>%
    filter(!is.na(gfr) & !is.na(albumin_mg_g)) %>%
    mutate(
      gfr_days = days,
      albumin_days = days,
      days_apart = 0
    ) %>%
    select(patient_id, gfr_days, albumin_days, gfr, albumin_mg_g, days_apart)
  
  # Option 2: Get separate measurements for remaining patients
  patients_with_same_day <- unique(same_day_measurements$patient_id)
  
  # Separate GFR and albumin measurements for patients without same-day measurements
  gfr_data <- data_std %>%
    filter(!patient_id %in% patients_with_same_day & !is.na(gfr)) %>%
    select(patient_id, gfr_days = days, gfr) %>%
    group_by(patient_id, gfr_days) %>%
    slice(1) %>%  # Take first if multiple measurements on same day
    ungroup()
  
  albumin_data <- data_std %>%
    filter(!patient_id %in% patients_with_same_day & !is.na(albumin_mg_g)) %>%
    select(patient_id, albumin_days = days, albumin_mg_g) %>%
    group_by(patient_id, albumin_days) %>%
    slice(1) %>%  # Take first if multiple measurements on same day
    ungroup()
  
  # Create all possible pairs within time window (including same-day pairs)
  different_day_pairs <- gfr_data %>%
    inner_join(albumin_data, by = "patient_id", relationship = "many-to-many") %>%
    mutate(
      days_apart = abs(gfr_days - albumin_days)
    ) %>%
    filter(days_apart <= max_days_apart)
  
  # Combine same-day and different-day measurements
  all_pairs <- bind_rows(same_day_measurements, different_day_pairs)
  
  # Find best pair for each patient (closest in time, preferring same-day)
  best_pairs <- all_pairs %>%
    group_by(patient_id) %>%
    slice_min(days_apart, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Classify each patient
  if(nrow(best_pairs) > 0) {
    results <- best_pairs %>%
      rowwise() %>%
      mutate(
        classification = list(
          classify_ckd(
            gfr = gfr,
            albumin_mg_g = albumin_mg_g,
            gfr_days_from_baseline = gfr_days,
            albumin_days_from_baseline = albumin_days,
            max_days_apart = max_days_apart
          )
        )
      ) %>%
      ungroup() %>%
      mutate(
        gfr_category = sapply(classification, function(x) x$gfr_category),
        gfr_description = sapply(classification, function(x) x$gfr_description),
        albumin_category = sapply(classification, function(x) x$albumin_category),
        albumin_description = sapply(classification, function(x) x$albumin_description),
        risk_level = sapply(classification, function(x) x$risk_level),
        recommendation = sapply(classification, function(x) x$recommendation),
        criteria_met_days_from_baseline = sapply(classification, function(x) x$criteria_met_days_from_baseline),
        stage = paste0(gfr_category, "/", albumin_category)
      ) %>%
      select(-classification)
    
    return(results)
  } else {
    return(data.frame())
  }
}

# Function to get first occurrence of each CKD stage per patient
get_first_ckd_stage <- function(data,
                                patient_id_col = "releaseid",
                                days_col = "days",
                                gfr_col = "ckd_gfr",
                                albumin_col = "UAlbCreat_mg_g",
                                max_days_apart = 365,
                                target_stages = NULL) {
  
  library(dplyr)
  
  # Get all classifications
  all_classifications <- classify_ckd_long(
    data,
    patient_id_col = patient_id_col,
    days_col = days_col,
    gfr_col = gfr_col,
    albumin_col = albumin_col,
    max_days_apart = max_days_apart
  )
  
  if(nrow(all_classifications) == 0) {
    return(data.frame())
  }
  
  # Filter to target stages if specified
  if(!is.null(target_stages)) {
    all_classifications <- all_classifications %>%
      filter(stage %in% target_stages)
  }
  
  # Get first occurrence of each stage per patient
  first_stages <- all_classifications %>%
    group_by(patient_id, stage) %>%
    slice_min(criteria_met_days_from_baseline, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(patient_id, criteria_met_days_from_baseline)
  
  return(first_stages)
}

# Function to track CKD progression over time
track_ckd_progression <- function(data,
                                  patient_id_col = "releaseid",
                                  days_col = "days",
                                  gfr_col = "ckd_gfr",
                                  albumin_col = "UAlbCreat_mg_g",
                                  time_windows = c(90, 180, 365),
                                  max_days_apart = 365) {
  
  library(dplyr)
  library(purrr)
  
  # Get classifications for each time window
  progression_results <- map_df(time_windows, function(window) {
    
    # Filter data to time window
    data_window <- data %>%
      filter(get(days_col) <= window)
    
    # Classify within this window
    classifications <- classify_ckd_long(
      data_window,
      patient_id_col = patient_id_col,
      days_col = days_col,
      gfr_col = gfr_col,
      albumin_col = albumin_col,
      max_days_apart = max_days_apart
    )
    
    if(nrow(classifications) > 0) {
      classifications %>%
        mutate(time_window = window) %>%
        group_by(patient_id) %>%
        slice(1) %>%  # Take first valid classification in window
        ungroup()
    } else {
      data.frame()
    }
  })
  
  return(progression_results)
}

# Helper function to reshape long data for visualization
prepare_ckd_timeline <- function(data,
                                 patient_id_col = "releaseid",
                                 days_col = "days",
                                 gfr_col = "ckd_gfr",
                                 albumin_col = "UAlbCreat_mg_g") {
  
  library(dplyr)
  library(tidyr)
  
  # Standardize column names
  timeline_data <- data %>%
    rename(
      patient_id = !!sym(patient_id_col),
      days = !!sym(days_col),
      gfr = !!sym(gfr_col),
      albumin_mg_g = !!sym(albumin_col)
    ) %>%
    mutate(
      test_type = case_when(
        !is.na(gfr) ~ "GFR",
        !is.na(albumin_mg_g) ~ "Albumin",
        TRUE ~ NA_character_
      ),
      value = coalesce(gfr, albumin_mg_g)
    ) %>%
    filter(!is.na(test_type)) %>%
    select(patient_id, days, test_type, value, gfr, albumin_mg_g)
  
  return(timeline_data)
}

# Example usage
library(dplyr)
library(tidyr)

# Example with your data format
your_data <- data.frame(
  releaseid = c("65-10062", "65-10062", "65-10062", "65-10062", "65-10062", "65-10062",
                "65-10111", "65-10111", "65-10111", "65-10111"),
  days = c(0, 56, 360, 448, 556, 766,
           0, 90, 180, 365),
  ckd_gfr = c(142.2855, NA, 146.7624, NA, NA, 142.0207,
              135.1844, NA, 125.5, NA),
  UAlbCreat_mg_g = c(11, NA, 30, NA, 149, 146,
                     218, 250, NA, 300)
)

# Classify CKD stages
results <- classify_ckd_long(your_data)
print("CKD Classification Results:")
print(results %>% select(patient_id, gfr_days, albumin_days, days_apart,
                         gfr, albumin_mg_g, stage, risk_level, 
                         criteria_met_days_from_baseline))

# Get first occurrence of each stage
first_stages <- get_first_ckd_stage(your_data)
print("\nFirst occurrence of each CKD stage:")
print(first_stages %>% select(patient_id, stage, risk_level, 
                              criteria_met_days_from_baseline))

# Track progression over time windows
progression <- track_ckd_progression(
  your_data,
  time_windows = c(90, 180, 365, 730)
)
print("\nCKD Progression over time:")
print(progression %>% select(patient_id, time_window, stage, risk_level))