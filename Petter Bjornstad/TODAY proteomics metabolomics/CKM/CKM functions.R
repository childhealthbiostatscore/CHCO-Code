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
  
  # fill down the GFR and albuminuria categories
  
  
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
      albumin_description = sapply(classification, function(x) x$albumin_description)
    ) %>% 
    group_by(releaseid) %>% fill(gfr_category, .direction = "down") %>%
    fill(gfr_description, .direction = "down") %>% 
    fill(albumin_category, .direction = "down") %>% 
    fill(albumin_description, .direction = "down") %>%
    ungroup() %>%
    mutate(
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

# function that returns long df without summarization
extract_biomarker_values_long <- function(df_list, var) {
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
        select(releaseid, days, all_of(var))
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
   # Clean up infinite values
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.), NA_real_, .))) 
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

#################################
# FILL IN NA IN EVENTS  - V2    #
#################################

# Using tidyr's replace_na function
replace_na_tidyverse <- function(df, columns, value = 0) {
  df %>%
    mutate(across(all_of(columns), ~replace_na(.x, value)))
}

# Example usage:
#result <- replace_na_tidyverse(sample_df, c("a", "b"))


# ===========================================================================
# Function: run_fgsea_analysis
# ===========================================================================

run_fgsea_analysis <- function(bg_path = file.path(root_path, "GSEA/"),
                               results_annotated,
                               stat_col = "t",
                               gene_col = "EntrezGeneSymbol",
                               minSize_kegg = 3,
                               maxSize_kegg = 500,
                               minSize_reactome = 3,
                               maxSize_reactome = 500,
                               minSize_go = 5,
                               maxSize_go = 500,
                               minSize_full = 5,
                               maxSize_full = 500,
                               minSize_hallmark = 5,
                               maxSize_hallmark = 500,
                               nPermSimple = 10000,
                               nproc = 1,
                               seed = 1234,
                               references = c("kegg_legacy", "reactome", "go", "full", "hallmark")) {
  
  # --- Prepare GMT files ---
  gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
  
  # Initialize pathway variables
  kegg_legacy <- NULL
  reactome <- NULL
  go <- NULL
  full <- NULL
  hallmark <- NULL
  
  # Only prepare GMT files that are in references
  if ("kegg_legacy" %in% references) {
    kegg_legacy <- prepare_gmt(gmt_files[1], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  if ("reactome" %in% references) {
    reactome <- prepare_gmt(gmt_files[3], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  if ("go" %in% references) {
    go <- prepare_gmt(gmt_files[4], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  if ("full" %in% references) {
    full <- prepare_gmt(gmt_files[5], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  if ("hallmark" %in% references) {
    hallmark <- prepare_gmt(gmt_files[6], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  
  # --- Rank genes by specified statistic ---
  if (!stat_col %in% names(results_annotated)) {
    stop(paste("Column", stat_col, "not found in results_annotated"))
  }
  
  rankings <- results_annotated[[stat_col]]
  names(rankings) <- results_annotated[[gene_col]]
  rankings <- sort(rankings, decreasing = TRUE)
  
  # --- Run FGSEA ---
  set.seed(seed)
  
  # Initialize result variables as NULL
  kegg_res <- NULL
  reactome_res <- NULL
  go_res <- NULL
  full_res <- NULL
  hallmark_res <- NULL
  
  if ("kegg_legacy" %in% references) {
    kegg_res <- fgsea(pathways = kegg_legacy,
                      stats = rankings,
                      scoreType = 'std', 
                      minSize = minSize_kegg,
                      maxSize = maxSize_kegg,
                      nproc = nproc,
                      nPermSimple = nPermSimple)
  }
  
  if ("reactome" %in% references) {
    reactome_res <- fgsea(pathways = reactome,
                          stats = rankings,
                          scoreType = 'std', 
                          minSize = minSize_reactome,
                          maxSize = maxSize_reactome,
                          nproc = nproc,
                          nPermSimple = nPermSimple)
  }
  
  if ("go" %in% references) {
    go_res <- fgsea(pathways = go,
                    stats = rankings,
                    scoreType = "std",
                    minSize = minSize_go,
                    maxSize = maxSize_go,
                    nPermSimple = nPermSimple,
                    nproc = nproc)
  }
  
  if ("full" %in% references) {
    full_res <- fgsea(pathways = full,
                      stats = rankings,
                      scoreType = "std",
                      minSize = minSize_full,
                      maxSize = maxSize_full,
                      nPermSimple = nPermSimple,
                      nproc = nproc)
  }
  
  if ("hallmark" %in% references) {
    hallmark_res <- fgsea(pathways = hallmark,
                          stats = rankings,
                          scoreType = "std",
                          minSize = minSize_hallmark,
                          maxSize = maxSize_hallmark,
                          nPermSimple = nPermSimple,
                          nproc = nproc)
  }
  
  # --- Build summary dataframe dynamically ---
  summary_list <- list()
  
  if ("kegg_legacy" %in% references && !is.null(kegg_res)) {
    summary_list[["KEGG Legacy"]] <- c(
      sum(kegg_res$padj < 0.05, na.rm = TRUE),
      sum(kegg_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  if ("reactome" %in% references && !is.null(reactome_res)) {
    summary_list[["REACTOME"]] <- c(
      sum(reactome_res$padj < 0.05, na.rm = TRUE),
      sum(reactome_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  if ("go" %in% references && !is.null(go_res)) {
    summary_list[["GO"]] <- c(
      sum(go_res$padj < 0.05, na.rm = TRUE),
      sum(go_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  if ("full" %in% references && !is.null(full_res)) {
    summary_list[["FULL"]] <- c(
      sum(full_res$padj < 0.05, na.rm = TRUE),
      sum(full_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  if ("hallmark" %in% references && !is.null(hallmark_res)) {
    summary_list[["HALLMARK"]] <- c(
      sum(hallmark_res$padj < 0.05, na.rm = TRUE),
      sum(hallmark_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  # Convert list to dataframe
  if (length(summary_list) > 0) {
    summary_df <- as.data.frame(summary_list)
    rownames(summary_df) <- c("adj.pval", "p.val")
  } else {
    summary_df <- data.frame()
  }
  
  # --- Return results ---
  return(list(
    summary = summary_df,
    kegg = if("kegg_legacy" %in% references) kegg_res else NULL,
    reactome = if("reactome" %in% references) reactome_res else NULL,
    go = if("go" %in% references) go_res else NULL,
    full = if("full" %in% references) full_res else NULL,
    hallmark = if("hallmark" %in% references) hallmark_res else NULL
  ))
}

# ===========================================================================
# Function: process_enrichr_data
# ===========================================================================

# Function to process enrichr results
process_enrichr_data <- function(enrichr_result, response_type, top_n = 20) {
  # Extract the Reactome 2022 results (or whichever database you used)
  data <- enrichr_result[["Reactome_Pathways_2024"]] # Adjust database name as needed
  
  # Clean and process data
  data_clean <- data %>%
    filter(Adjusted.P.value < 0.05) %>%  # Filter significant pathways
    head(top_n) %>%
    dplyr::mutate(
      neg_log_p = -log10(P.value),
      neg_log_adj_p = -log10(Adjusted.P.value),
      gene_count = as.numeric(str_extract(Overlap, "\\d+")),
      response_type = response_type,
      # Create simplified pathway names for plotting
      pathway_short = str_trunc(Term, 40),
      # Assign categories based on pathway names
      category = case_when(
        str_detect(Term, "Immune|Neutrophil|Cytokine|Interferon|Interleukin") ~ "Immune",
        str_detect(Term, "Metabolism|Metabolic|Amino Acid|Fatty Acid|Glucose") ~ "Metabolism", 
        str_detect(Term, "Signal|Signaling|Receptor|Tyrosine|Growth Factor") ~ "Signaling",
        str_detect(Term, "Proteasome|Autophagy|Ubiquitin|Degradation|Quality") ~ "Quality Control",
        str_detect(Term, "Stress|Oxidative|Response|NFE2L2|KEAP1") ~ "Stress Response",
        str_detect(Term, "Mitochondrial|Respiratory|Electron|ATP") ~ "Mitochondrial",
        str_detect(Term, "Apoptosis|Cell Death|Programmed") ~ "Cell Death",
        str_detect(Term, "Transport|Trafficking|Vesicle|Membrane") ~ "Transport",
        TRUE ~ "Other"
      )
    )
  
  return(data_clean)
}

# ===========================================================================
# Function: prepare_pathway_data
# ===========================================================================

# Process enrichr results
prepare_pathway_data <- function(negative_paths_df, positive_paths_df, discordant_paths_df) {
  # Process each dataset
  negative_paths <- process_enrichr_data(negative_paths_df, "Negative")
  positive_paths <- process_enrichr_data(positive_paths_df, "Positive") 
  discordant_paths <- process_enrichr_data(discordant_paths_df, "Discordant")
  
  return(list(
    negative = negative_paths,
    positive = positive_paths,
    discordant = discordant_paths
  ))
}

# ---- Pathway Plot Functions ----

# ===========================================================================
# Function: matrix_to_list
# ===========================================================================
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
# ===========================================================================
# Function: prepare_gmt
# ===========================================================================
library(fgsea)
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# List of acronyms to uppercase
acronyms <- c("rna", "dna", "mtor", "foxo", "ppar", "nmd", "fgfr", "robo", 
              "bhl", "cov", "jak", "stat", "wnt", "hiv", "bcl", "mapk",
              "pt", "tal", "pc", "ic", "ec", "fibvsmcp")
special_mixed <- c("rrna", "mrna", "trna", "gtpase", "atpase", "robos", "slits", "fibvsmcp")
special_replacements <- c("rRNA", "mRNA", "tRNA", "GTPase", "ATPase", "ROBOs", "SLITs", "FIB/VSMC/P")


# ===========================================================================
# Function: replace_mixed_case
# ===========================================================================

replace_mixed_case <- function(text, from, to) {
  for (i in seq_along(from)) {
    pattern <- paste0("\\b", from[i], "\\b")
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), to[i])
  }
  return(text)
}


# ===========================================================================
# Function: capitalize_acronyms
# ===========================================================================

capitalize_acronyms <- function(text, terms) {
  for (term in terms) {
    pattern <- paste0("\\b", term, "\\b")
    replacement <- toupper(term)
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), replacement)
  }
  return(text)
}

# ===========================================================================
# Function: clean_pathway_names
# ===========================================================================

clean_pathway_names <- function(pathways) {
  # Remove REACTOME_ prefix
  cleaned <- gsub("^REACTOME_", "", pathways)
  cleaned <- gsub("^GOBP_", "", cleaned)
  cleaned <- gsub("^KEGG_", "", cleaned)
  cleaned <- gsub("^HALLMARK_", "", cleaned)
  
  # Replace underscores with spaces
  cleaned <- gsub("_", " ", cleaned)
  
  # Convert to title case (capitalize first letter of each word)
  cleaned <- tools::toTitleCase(tolower(cleaned))
  
  # Define patterns that should be in uppercase
  uppercase_words <- c(
    # Roman numerals
    "\\bI\\b", "\\bIi\\b", "\\bIii\\b", "\\bIv\\b", "\\bV\\b", "\\bVi\\b", 
    "\\bVii\\b", "\\bViii\\b", "\\bIx\\b", "\\bX\\b",
    # Common biological abbreviations
    "\\bTca\\b", "\\bAtp\\b", "\\bAdp\\b", "\\bAmp\\b", "\\bGtp\\b", "\\bGdp\\b",
    "\\bNad\\b", "\\bNadh\\b", "\\bFad\\b", "\\bFadh2\\b", "\\bCoa\\b",
    "\\bDna\\b", "\\bRna\\b", "\\bMrna\\b", "\\bTrna\\b", "\\bRrna\\b",
    "\\bEr\\b", "\\bUpr\\b", "\\bNf\\b", "\\bHif\\b", "\\bMhc\\b",
    "\\bTgf\\b", "\\bEgf\\b", "\\bVegf\\b", "\\bPdgf\\b", "\\bFgf\\b",
    "\\bRos\\b", "\\bRns\\b", "\\bNo\\b", "\\bNos\\b", "\\bInos\\b", "\\bEnos\\b", "\\bNnos\\b",
    "\\Mapk\\b"
  )
  
  # Apply uppercase replacements
  for (pattern in uppercase_words) {
    # Extract the word without the word boundaries
    word <- gsub("\\\\b", "", pattern)
    # Convert to uppercase
    replacement <- toupper(word)
    # Replace in the cleaned string
    cleaned <- gsub(pattern, replacement, cleaned, ignore.case = TRUE)
  }
  
  # Fix specific cases that might need custom handling
  cleaned <- gsub("\\bOf\\b", "of", cleaned)  # lowercase 'of'
  cleaned <- gsub("\\bBy\\b", "by", cleaned)  # lowercase 'by'
  cleaned <- gsub("\\bThe\\b", "the", cleaned)  # lowercase 'the' (except at start)
  cleaned <- gsub("^the", "The", cleaned)  # capitalize 'the' at start
  cleaned <- gsub("\\bO Linked\\b", "O-Linked", cleaned, ignore.case = TRUE)
  
  return(cleaned)
}

# ===========================================================================
# Function: filter_redundant_pathways
# ===========================================================================

filter_redundant_pathways <- function(gsea_result, overlap_pct = 0.3) {
  # Check required columns
  if (!all(c("pathway", "leadingEdge", "padj") %in% colnames(gsea_result))) {
    stop("Input data frame must have 'pathway', 'leadingEdge', and 'padj' columns.")
  }
  
  # Extract leading edge and name them
  leading_edges <- gsea_result$leadingEdge
  names(leading_edges) <- gsea_result$pathway
  
  # Compute Jaccard overlap matrix
  overlap_matrix <- sapply(leading_edges, function(x) 
    sapply(leading_edges, function(y) 
      length(intersect(x, y)) / length(union(x, y))
    )
  )
  
  # Identify pairs with high overlap
  overlap_pairs <- which(overlap_matrix > overlap_pct & lower.tri(overlap_matrix), arr.ind = TRUE)
  
  # If no overlaps found, return original
  if (nrow(overlap_pairs) == 0) {
    message("No redundant pathways found with overlap above threshold.")
    return(gsea_result)
  }
  
  # Create overlap graph
  library(igraph)
  edges <- data.frame(
    from = rownames(overlap_matrix)[overlap_pairs[,1]],
    to   = colnames(overlap_matrix)[overlap_pairs[,2]]
  )
  g <- graph_from_data_frame(edges, directed = FALSE)
  components <- components(g)
  
  # Cluster terms and select representative from each
  redundant_clusters <- split(names(components$membership), components$membership)
  representative_terms <- sapply(redundant_clusters, function(cluster_terms) {
    sub <- gsea_result[gsea_result$pathway %in% cluster_terms, ]
    sub[which.min(sub$padj), "pathway"]
  })
  
  # Keep only representative + non-overlapping terms
  all_overlapping <- unique(unlist(redundant_clusters))
  non_overlapping <- setdiff(gsea_result$pathway, all_overlapping)
  final_terms <- c(non_overlapping, representative_terms)
  
  # Return filtered GSEA result
  gsea_result[gsea_result$pathway %in% final_terms, ]
}

# ===========================================================================
# Function: plot_fgsea_transpose
# ===========================================================================

plot_fgsea_transpose <- function(fgsea_res,
                                 top_n = 30,
                                 title = "Top Enriched Pathways",
                                 xmin = 0,
                                 xmax = 3,
                                 xnudge = (xmax - xmin)/100,
                                 text1 = 6.5,
                                 text2 = 18,
                                 text3 = 20) {
  
  fgsea_res <- fgsea_res %>%
    dplyr::arrange(pval) %>%
    head(top_n) %>%
    dplyr::mutate(
      direction = case_when((NES < 0 & pval <= 0.05 ~ "Negative"), 
                            (NES > 0 & pval <= 0.05 ~ "Positive"),
                            (NES < 0 & pval > 0.05 ~ "Negative p > 0.05"), 
                            (NES > 0 & pval > 0.05 ~ "Positive p > 0.05")),
      face = case_when((NES < 0 & pval <= 0.05 ~ "bold"), 
                       (NES > 0 & pval <= 0.05 ~ "bold"),
                       (NES < 0 & pval > 0.05 ~ "plain"), 
                       (NES > 0 & pval > 0.05 ~ "plain")),
      pathway_clean = str_remove(pathway, "^KEGG_"), 
      pathway_clean = str_remove(pathway_clean, "^REACTOME_"), 
      pathway_clean = str_remove(pathway_clean, "^GOBP_"), 
      pathway_clean = str_remove(pathway_clean, "^GOMF_"), 
      pathway_clean = str_replace_all(pathway_clean, "_", " "),
      pathway_clean = str_to_sentence(pathway_clean),
      pathway_clean = str_replace_all(pathway_clean, "\\bi\\b", "I"),
      pathway_clean = str_replace_all(pathway_clean, "\\bii\\b", "II"),
      pathway_clean = str_replace_all(pathway_clean, "\\biii\\b", "III"),
      pathway_clean = str_replace_all(pathway_clean, "\\biv\\b", "IV"),
      pathway_clean = str_replace_all(pathway_clean, "\\bv\\b", "V"),
      pathway_clean = str_replace_all(pathway_clean, regex("\\(immune\\)", ignore_case = TRUE), "(IMMUNE)"),
      pathway_clean = capitalize_acronyms(pathway_clean, acronyms),
      pathway_clean = replace_mixed_case(pathway_clean, special_mixed, special_replacements),
      pathway_clean = paste0(pathway_clean, " (", size, ")")
    ) %>%
    dplyr::arrange(pval)
  
  fgsea_res$pathway_clean <- reorder(fgsea_res$pathway_clean, -abs(fgsea_res$NES))
  
  fgsea_res %>%
    ggplot(aes(x = abs(NES), y = fct_rev(pathway_clean), label = pathway_clean)) +
    geom_point(aes(size = -log10(pval), color = direction, alpha = 0.8)) +
    # geom_vline(xintercept = 2, linetype = "dashed") +
    geom_text(aes(group = pathway_clean, color = direction, fontface = face), 
              hjust = 0, size = text1, nudge_x = xnudge) +
    scale_size_binned() +
    scale_color_manual(values = c("Positive" = "#c75146", "Negative" = "#2c7da0", 
                                  "Positive p > 0.05" = "#e18c80", "Negative p > 0.05" = "#7ab6d1")) +
    scale_x_continuous(limits = c(xmin, xmax), expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(add = 1)) +
    labs(
      x = "NES",
      y = "Pathways",
      color = "Direction",
      size = "-log(p-value)",
      title = title
    ) +
    guides(alpha = "none") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = text3),
      axis.title = element_text(size = text3),
      axis.ticks.y = element_blank(), 
      legend.position = c(0.9, 0.2),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(color = "black"),
      legend.title = element_text(size = text2),
      legend.text = element_text(size = text2),
      title = element_text(size = text3)
    )
}

