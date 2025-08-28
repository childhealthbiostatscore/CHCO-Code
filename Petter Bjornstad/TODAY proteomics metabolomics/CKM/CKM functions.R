#################################
# CKM CLASSIFICATION            #
#################################

# KDIGO CKD Staging System Implementation in R
# Based on KDIGO 2012 Clinical Practice Guideline

# Function to classify CKD stage based on GFR and Albuminuria
classify_ckd <- function(gfr = NULL, albumin_mg_g = NULL, albumin_mg_mmol = NULL) {
  
  # Convert albumin units if needed
  if (!is.null(albumin_mg_mmol) && is.null(albumin_mg_g)) {
    # Convert mg/mmol to mg/g (multiply by 8.84)
    albumin_mg_g <- albumin_mg_mmol * 8.84
  }
  
  # Handle missing values - allow partial classification
  gfr_available <- !is.null(gfr) && !is.na(gfr)
  albumin_available <- !is.null(albumin_mg_g) && !is.na(albumin_mg_g)
  
  # If both are missing, return all NAs
  if (!gfr_available && !albumin_available) {
    return(list(
      gfr_category = NA,
      gfr_description = NA,
      albumin_category = NA,
      albumin_description = NA,
      risk_level = NA,
      recommendation = NA
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
  if (gfr_available && albumin_available) {
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
    recommendation = risk_recommendation$recommendation
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

# Function to classify multiple patients
classify_ckd_batch <- function(data, gfr_col = "gfr", albumin_mg_g_col = "albumin_mg_g", 
                               albumin_mg_mmol_col = NULL) {
  
  # Load required libraries
  library(dplyr)
  
  results <- data %>%
    rowwise() %>%
    mutate(
      classification = list(
        classify_ckd(
          gfr = get(gfr_col),
          albumin_mg_g = if(!is.null(albumin_mg_g_col)) get(albumin_mg_g_col) else NULL,
          albumin_mg_mmol = if(!is.null(albumin_mg_mmol_col)) get(albumin_mg_mmol_col) else NULL
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
      recommendation = sapply(classification, function(x) x$recommendation)
    ) %>%
    select(-classification)
  
  return(results)
}

# Function to create a summary visualization
plot_ckd_stage <- function(gfr, albumin_mg_g) {
  library(ggplot2)
  
  result <- classify_ckd(gfr, albumin_mg_g)
  
  # Create a simple plot showing the patient's position
  stage_plot <- ggplot() +
    # Add rectangles for each stage
    annotate("rect", xmin = 0, xmax = 30, ymin = 90, ymax = 150, 
             fill = "green", alpha = 0.3) +
    annotate("rect", xmin = 30, xmax = 300, ymin = 90, ymax = 150, 
             fill = "yellow", alpha = 0.3) +
    annotate("rect", xmin = 300, xmax = 1000, ymin = 90, ymax = 150, 
             fill = "orange", alpha = 0.3) +
    
    annotate("rect", xmin = 0, xmax = 30, ymin = 60, ymax = 90, 
             fill = "green", alpha = 0.3) +
    annotate("rect", xmin = 30, xmax = 300, ymin = 60, ymax = 90, 
             fill = "yellow", alpha = 0.3) +
    annotate("rect", xmin = 300, xmax = 1000, ymin = 60, ymax = 90, 
             fill = "orange", alpha = 0.3) +
    
    annotate("rect", xmin = 0, xmax = 30, ymin = 45, ymax = 60, 
             fill = "yellow", alpha = 0.3) +
    annotate("rect", xmin = 30, xmax = 300, ymin = 45, ymax = 60, 
             fill = "orange", alpha = 0.3) +
    annotate("rect", xmin = 300, xmax = 1000, ymin = 45, ymax = 60, 
             fill = "red", alpha = 0.3) +
    
    annotate("rect", xmin = 0, xmax = 30, ymin = 30, ymax = 45, 
             fill = "orange", alpha = 0.3) +
    annotate("rect", xmin = 30, xmax = 300, ymin = 30, ymax = 45, 
             fill = "red", alpha = 0.3) +
    annotate("rect", xmin = 300, xmax = 1000, ymin = 30, ymax = 45, 
             fill = "red", alpha = 0.3) +
    
    annotate("rect", xmin = 0, xmax = 30, ymin = 15, ymax = 30, 
             fill = "red", alpha = 0.3) +
    annotate("rect", xmin = 30, xmax = 300, ymin = 15, ymax = 30, 
             fill = "red", alpha = 0.3) +
    annotate("rect", xmin = 300, xmax = 1000, ymin = 15, ymax = 30, 
             fill = "red", alpha = 0.3) +
    
    annotate("rect", xmin = 0, xmax = 30, ymin = 0, ymax = 15, 
             fill = "darkred", alpha = 0.3) +
    annotate("rect", xmin = 30, xmax = 300, ymin = 0, ymax = 15, 
             fill = "darkred", alpha = 0.3) +
    annotate("rect", xmin = 300, xmax = 1000, ymin = 0, ymax = 15, 
             fill = "darkred", alpha = 0.3) +
    
    # Add patient point
    geom_point(aes(x = albumin_mg_g, y = gfr), size = 5, color = "black") +
    
    # Add labels
    scale_x_log10(breaks = c(30, 300), limits = c(10, 1000)) +
    scale_y_continuous(breaks = c(15, 30, 45, 60, 90), limits = c(0, 150)) +
    
    labs(
      title = paste("CKD Stage:", result$gfr_category, "/", result$albumin_category),
      subtitle = paste("Risk Level:", result$risk_level, "| Recommendation:", result$recommendation),
      x = "Albuminuria (mg/g)",
      y = "GFR (mL/min/1.73 m²)"
    ) +
    
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_blank()
    )
  
  return(stage_plot)
}

# Example usage
# # Single patient classification - both values available
# patient_result1 <- classify_ckd(gfr = 55, albumin_mg_g = 150)
# print("Patient with both values:")
# print(patient_result1)
# 
# # Patient with only GFR available
# patient_result2 <- classify_ckd(gfr = 45, albumin_mg_g = NA)
# print("\nPatient with only GFR:")
# print(patient_result2)
# 
# # Patient with only albumin available
# patient_result3 <- classify_ckd(gfr = NA, albumin_mg_g = 350)
# print("\nPatient with only albumin:")
# print(patient_result3)
# 
# # Multiple patients with some missing values
# sample_data <- data.frame(
#   patient_id = 1:5,
#   gfr = c(95, 75, 50, NA, 20),
#   albumin_mg_g = c(20, NA, 45, 350, 500)
# )
# 
# results <- classify_ckd_batch(sample_data)
# print("\nBatch classification with missing values:")
# print(results)

# Visualize a patient's stage (only works with both values)
# plot_ckd_stage(gfr = 55, albumin_mg_g = 150)

#################################
# BIOMARKER EXTRACTION          #
#################################
# Need to modify function to extract data across both TODAY and TODAY2
extract_biomarker_values <- function(df_list, var) {
  # If a single dataframe is passed, convert to list
  if(is.data.frame(df_list)) {
    df_list <- list(df_list)
  }
  # Combine all dataframes that contain the variable
  combined_df <- NULL
  for(df in df_list) {
    if(var %in% names(df)) {
      # Select only the columns we need
      temp_df <- df %>%
        select(releaseid, days, all_of(var), any_of("mvisit"))
      if(is.null(combined_df)) {
        combined_df <- temp_df
      } else {
        combined_df <- bind_rows(combined_df, temp_df)
      }
    }
  }
  # If variable not found in any dataset, return empty dataframe
  if(is.null(combined_df)) {
    warning(paste("Variable", var, "not found in any of the provided datasets"))
    return(data.frame())
  }
  combined_df %>%
    arrange(releaseid, days) %>%
    group_by(releaseid) %>%
    summarise(
      # Baseline: only M00 visit, NA if not present (unchanged)
      baseline = {
        if("mvisit" %in% names(combined_df)) {
          m00_row <- which(mvisit == "M00" & !is.na(.data[[var]]))
          if(length(m00_row) > 0) {
            .data[[var]][m00_row[1]] 
          } else {
            NA_real_
          }
        } else {
          NA_real_
        }
      },
      # Highest now looks across all combined data
      highest = max(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Clean up infinite values
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.), NA_real_, .))) %>%
    rename_with(~ paste0(var, "_", .), -releaseid)
}

#################################
# FILL IN NA IN EVENTS          #
#################################
replace_missing_events <- function(data) {
  # Check if fup_time column exists
  if (!"fup_time" %in% names(data)) {
    stop("Column 'fup_time' not found in the data")
  }
  # Define the event pairs
  event_pairs <- list(
    c("ARRHYTHMIA", "DAYSTOARRHYTHMIA"),
    c("CAD", "DAYSTOCAD"),
    c("CHF", "DAYSTOCHF"),
    c("LVSD", "DAYSTOLVSD"),
    c("MI", "DAYSTOMI"),
    c("PAD", "DAYSTOPAD"),
    c("DVT", "DAYSTODVT"),
    c("STROKE", "DAYSTOSTROKE"),
    c("TIA", "DAYSTOTIA"),
    c("CKD", "DAYSTOCKD"),
    c("ESKD", "DAYSTOESKD"),
    c("DEATH", "DAYSTODEATH"),
    c("ULCER", "DAYSTOULCER")
  )
  # Loop through each event pair
  for (pair in event_pairs) {
    event_var <- pair[1]
    days_var <- pair[2]
    # Step 1: Set missing event values to 0
    data[[event_var]][is.na(data[[event_var]])] <- 0
    # Step 2: Now set any missing DAYSTO values to fup_time
    # (including those that were missing because the event was missing/0)
    missing_days_rows <- is.na(data[[days_var]])
    data[[days_var]][missing_days_rows] <- data$fup_time[missing_days_rows]
  }
  return(data)
}




