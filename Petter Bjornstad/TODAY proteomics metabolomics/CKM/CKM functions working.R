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

# Function to get ALL VALID CKD classifications over time
# Only pairs measurements that make clinical sense
classify_ckd_longitudinal <- function(data,
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
  
  # First priority: Get all same-day measurements (both values in same row)
  same_day_measurements <- data_std %>%
    filter(!is.na(gfr) & !is.na(albumin_mg_g)) %>%
    mutate(
      gfr_days = days,
      albumin_days = days,
      days_apart = 0,
      pair_type = "same_day"
    ) %>%
    select(patient_id, gfr_days, albumin_days, gfr, albumin_mg_g, days_apart, pair_type)
  
  # For each patient, identify which days have same-day measurements
  patients_with_same_day <- same_day_measurements %>%
    select(patient_id, days = gfr_days) %>%
    distinct()
  
  # Second priority: For days without same-day measurements, find nearest pairs
  # Get all GFR and albumin measurements
  gfr_data <- data_std %>%
    filter(!is.na(gfr)) %>%
    select(patient_id, gfr_days = days, gfr) %>%
    distinct()
  
  albumin_data <- data_std %>%
    filter(!is.na(albumin_mg_g)) %>%
    select(patient_id, albumin_days = days, albumin_mg_g) %>%
    distinct()
  
  # For each unique time point, find the closest pairing
  # This prevents creating artificial pairs from distant time points
  nearest_pairs <- NULL
  
  for (pid in unique(data_std$patient_id)) {
    # Get this patient's data
    patient_gfr <- gfr_data %>% filter(patient_id == pid)
    patient_albumin <- albumin_data %>% filter(patient_id == pid)
    patient_same_day <- patients_with_same_day %>% filter(patient_id == pid)
    
    # For each GFR measurement without same-day albumin
    for (i in 1:nrow(patient_gfr)) {
      gfr_day <- patient_gfr$gfr_days[i]
      
      # Skip if this day already has same-day measurement
      if (gfr_day %in% patient_same_day$days) next
      
      # Find closest albumin measurement within window
      albumin_distances <- abs(patient_albumin$albumin_days - gfr_day)
      if (min(albumin_distances) <= max_days_apart) {
        closest_idx <- which.min(albumin_distances)
        
        pair <- data.frame(
          patient_id = pid,
          gfr_days = gfr_day,
          albumin_days = patient_albumin$albumin_days[closest_idx],
          gfr = patient_gfr$gfr[i],
          albumin_mg_g = patient_albumin$albumin_mg_g[closest_idx],
          days_apart = albumin_distances[closest_idx],
          pair_type = "nearest"
        )
        
        nearest_pairs <- bind_rows(nearest_pairs, pair)
      }
    }
  }
  
  # Remove duplicate pairs (keep only one pair per GFR-albumin combination)
  if (!is.null(nearest_pairs)) {
    nearest_pairs <- nearest_pairs %>%
      group_by(patient_id, gfr_days, albumin_days) %>%
      slice(1) %>%
      ungroup()
  }
  
  # Combine same-day and nearest pairs
  all_pairs <- bind_rows(same_day_measurements, nearest_pairs) %>%
    arrange(patient_id, pmin(gfr_days, albumin_days))
  
  # Classify each pair
  if(nrow(all_pairs) > 0) {
    results <- all_pairs %>%
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

# Updated function to track CKD progression events
track_ckd_progression_events <- function(data,
                                         patient_id_col = "releaseid",
                                         days_col = "days",
                                         gfr_col = "ckd_gfr",
                                         albumin_col = "UAlbCreat_mg_g",
                                         max_days_apart = 365) {
  
  library(dplyr)
  
  # Get ALL classifications over time (not just best pair)
  all_classifications <- classify_ckd_longitudinal(
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
  
  # Simplify risk levels: combine Moderate and High
  all_classifications <- all_classifications %>%
    mutate(
      simplified_risk = case_when(
        risk_level == "Low" ~ "Low",
        risk_level %in% c("Moderate", "High") ~ "Moderate/High",
        risk_level == "Very high" ~ "Very high",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(simplified_risk))
  
  # For each patient, track progression
  progression_events <- all_classifications %>%
    group_by(patient_id) %>%
    arrange(criteria_met_days_from_baseline) %>%
    mutate(
      # Get the baseline (first) risk level for each patient
      baseline_risk = first(simplified_risk),
      
      # Identify progression events
      progression_event = case_when(
        # Low -> Moderate/High
        baseline_risk == "Low" & simplified_risk == "Moderate/High" ~ "Low to Moderate/High",
        # Low -> Very high
        baseline_risk == "Low" & simplified_risk == "Very high" ~ "Low to Very high",
        # For patients starting at Moderate/High, only Very high is progression
        baseline_risk == "Moderate/High" & simplified_risk == "Very high" ~ "Moderate/High to Very high",
        TRUE ~ NA_character_
      )
    ) %>%
    # Keep only progression events
    filter(!is.na(progression_event)) %>%
    # Get first occurrence of each type of progression
    group_by(patient_id, progression_event) %>%
    slice_min(criteria_met_days_from_baseline, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  return(progression_events)
}

# Function to get first meaningful progression for each patient
get_first_progression <- function(data,
                                  patient_id_col = "releaseid",
                                  days_col = "days",
                                  gfr_col = "ckd_gfr",
                                  albumin_col = "UAlbCreat_mg_g",
                                  max_days_apart = 365) {
  
  library(dplyr)
  
  # Get all progression events
  all_progressions <- track_ckd_progression_events(
    data,
    patient_id_col = patient_id_col,
    days_col = days_col,
    gfr_col = gfr_col,
    albumin_col = albumin_col,
    max_days_apart = max_days_apart
  )
  
  if(nrow(all_progressions) == 0) {
    return(data.frame())
  }
  
  # Get first progression for each patient (regardless of type)
  first_progression <- all_progressions %>%
    group_by(patient_id) %>%
    slice_min(criteria_met_days_from_baseline, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  return(first_progression)
}

# Updated function to create a summary of progression patterns
summarize_progression <- function(data,
                                  patient_id_col = "releaseid",
                                  days_col = "days", 
                                  gfr_col = "ckd_gfr",
                                  albumin_col = "UAlbCreat_mg_g",
                                  max_days_apart = 365) {
  
  library(dplyr)
  
  # Get ALL classifications over time
  all_classifications <- classify_ckd_longitudinal(
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
  
  # Get progression events
  progression_events <- track_ckd_progression_events(
    data,
    patient_id_col = patient_id_col,
    days_col = days_col,
    gfr_col = gfr_col,
    albumin_col = albumin_col,
    max_days_apart = max_days_apart
  )
  
  # Create summary for each patient
  patient_summary <- all_classifications %>%
    mutate(
      simplified_risk = case_when(
        risk_level == "Low" ~ "Low",
        risk_level %in% c("Moderate", "High") ~ "Moderate/High", 
        risk_level == "Very high" ~ "Very high",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(simplified_risk)) %>%
    group_by(patient_id) %>%
    arrange(criteria_met_days_from_baseline) %>%
    summarise(
      baseline_stage = first(stage),
      baseline_risk = first(simplified_risk),
      baseline_day = first(criteria_met_days_from_baseline),
      final_stage = last(stage),
      final_risk = last(simplified_risk),
      final_day = last(criteria_met_days_from_baseline),
      n_assessments = n(),
      .groups = 'drop'
    )
  
  # Add progression information
  if(nrow(progression_events) > 0) {
    first_progression <- progression_events %>%
      group_by(patient_id) %>%
      slice_min(criteria_met_days_from_baseline, n = 1, with_ties = FALSE) %>%
      select(patient_id, 
             first_progression_event = progression_event,
             first_progression_day = criteria_met_days_from_baseline,
             first_progression_stage = stage)
    
    patient_summary <- patient_summary %>%
      left_join(first_progression, by = "patient_id") %>%
      mutate(
        progressed = !is.na(first_progression_event),
        days_to_progression = first_progression_day - baseline_day
      )
  } else {
    patient_summary <- patient_summary %>%
      mutate(
        progressed = FALSE,
        first_progression_event = NA_character_,
        first_progression_day = NA_real_,
        first_progression_stage = NA_character_,
        days_to_progression = NA_real_
      )
  }
  
  return(patient_summary)
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

# Example with progression scenarios - matching YOUR data structure
# (both GFR and UACR in same row when measured on same day)
progression_data <- data.frame(
  releaseid = c(
    # Patient 1: Low -> Moderate/High -> Very high
    rep("65-10001", 3),
    # Patient 2: Starts at Moderate, progresses to Very high
    rep("65-10002", 3),
    # Patient 3: Stays at Low (no progression)
    rep("65-10003", 2),
    # Patient 4: Low -> Very high (skips moderate)
    rep("65-10004", 2)
  ),
  days = c(
    # Patient 1
    0, 180, 365,
    # Patient 2  
    0, 180, 365,
    # Patient 3
    0, 180,
    # Patient 4
    0, 365
  ),
  ckd_gfr = c(
    # Patient 1: declining GFR (same-day measurements)
    95, 55, 35,
    # Patient 2: starting low
    50, 45, 25,
    # Patient 3: stable normal
    100, 98,
    # Patient 4: rapid decline
    90, 20
  ),
  UAlbCreat_mg_g = c(
    # Patient 1: increasing albuminuria (same-day measurements)
    20, 150, 400,
    # Patient 2: high albuminuria throughout
    180, 250, 500,
    # Patient 3: normal albuminuria
    15, 18,
    # Patient 4: severe albuminuria at end
    25, 600
  )
)

print("Sample progression data (your format - same-day tests in same row):")
print(progression_data)

# Track progression events
progression_events <- track_ckd_progression_events(progression_data)
print("\nProgression events identified:")
print(progression_events %>% 
        select(patient_id, baseline_risk, simplified_risk, progression_event,
               stage, criteria_met_days_from_baseline))

# Get first progression for each patient
first_progression <- get_first_progression(progression_data)
print("\nFirst progression event per patient:")
print(first_progression %>%
        select(patient_id, progression_event, stage, 
               criteria_met_days_from_baseline))

# Get comprehensive summary
summary <- summarize_progression(progression_data)
print("\nProgression summary by patient:")
print(summary %>%
        select(patient_id, baseline_risk, baseline_day, 
               progressed, first_progression_event, days_to_progression))

# Example with mixed scenario (corrected for actual progression)
# Some visits have both tests, some have only one
realistic_data <- data.frame(
  releaseid = c(
    rep("65-10062", 6),
    rep("65-10111", 4),
    rep("65-10222", 5)
  ),
  days = c(
    # Patient 1: progresses from Low to Moderate
    0, 56, 180, 360, 556, 766,
    # Patient 2: progresses from High to Very high
    0, 90, 180, 365,
    # Patient 3: progresses from Low to Very high
    0, 30, 180, 365, 730
  ),
  ckd_gfr = c(
    # Patient 1: Normal -> Mildly decreased
    95, NA, 92, 75, NA, 70,
    # Patient 2: Mildly to moderately decreased -> Moderately to severely decreased
    55, 52, 48, 35,
    # Patient 3: Normal -> Mildly to moderately decreased -> Severely decreased
    105, NA, 85, 45, 25
  ),
  UAlbCreat_mg_g = c(
    # Patient 1: Normal -> Moderately increased (A1 -> A2)
    15, 25, 28, 150, 180, 200,
    # Patient 2: Moderately increased -> Severely increased (A2 -> A3)
    150, 180, 250, 350,
    # Patient 3: Normal -> Severely increased (A1 -> A3)
    10, 15, NA, 350, 500
  )
)

print("\n\nRealistic data with corrected progression events:")
print(realistic_data)

# Test each patient to verify risk levels
test_classifications <- classify_ckd_longitudinal(realistic_data)
print("\nAll classifications to verify risk levels:")
test_classifications %>%
  filter(days_apart == 0) %>%  # Only look at same-day measurements
  select(patient_id, gfr_days, gfr, albumin_mg_g, stage, risk_level) %>%
  print()

# Analyze progression
realistic_summary <- summarize_progression(realistic_data)
print("\nProgression analysis (should show non-zero days):")
print(realistic_summary %>%
        select(patient_id, baseline_stage, baseline_risk, baseline_day,
               progressed, first_progression_event, first_progression_day, days_to_progression))

# Track all progression events
all_progressions <- track_ckd_progression_events(realistic_data)
print("\nAll progression events with timing:")
print(all_progressions %>%
        select(patient_id, baseline_risk, progression_event, stage,
               criteria_met_days_from_baseline, days_from_baseline_to_progression))

# For survival analysis - get time to first progression
survival_data <- get_first_progression(realistic_data)
print("\nTime to first progression (for survival analysis):")
print(survival_data %>%
        select(patient_id, progression_event, stage, criteria_met_days_from_baseline, 
               days_from_baseline_to_progression))