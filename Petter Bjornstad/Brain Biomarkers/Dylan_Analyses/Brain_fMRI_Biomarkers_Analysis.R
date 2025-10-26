############ Brain Biomarkers Analysis 


if (!require("oro.nifti")) install.packages("oro.nifti")
if (!require("neurobase")) install.packages("neurobase")
if (!require("R.matlab")) install.packages("R.matlab")

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("corrplot")) install.packages("corrplot")
if (!require("psych")) install.packages("psych")
if (!require("pls")) install.packages("pls")
if (!require("caret")) install.packages("caret")
if (!require("randomForest")) install.packages("randomForest")

 if (!require("igraph")) install.packages("igraph")
 if (!require("brainGraph")) install.packages("brainGraph")

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

harmonized_data <- read.csv("../OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


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



# Now create the table with proper data types

library(gt)
library(gtsummary)

desc_table1_fixed <- small_dat %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group, ab40_avg_conc, ab42_avg_conc, 
         tau_avg_conc, nfl_avg_conc, gfap_avg_conc, ptau_181_avg_conc, ptau_217_avg_conc) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical",
      ab40_avg_conc ~ 'continuous', 
      ab42_avg_conc ~ 'continuous', 
      tau_avg_conc ~ 'continuous', 
      nfl_avg_conc ~ 'continuous', 
      gfap_avg_conc ~ 'continuous',
      ptau_181_avg_conc ~ 'continuous', 
      ptau_217_avg_conc ~ 'continuous'
      
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study",
      ab40_avg_conc ~ 'Ab40 (pg/mL)', 
      ab42_avg_conc ~ 'Ab42 (pg/mL)', 
      tau_avg_conc ~ 'Tau (pg/mL)', 
      nfl_avg_conc ~ 'NFL (pg/mL)', 
      gfap_avg_conc ~ 'GFAP (pg/mL)',
      ptau_181_avg_conc ~ 'pTau 181 (pg/mL)', 
      ptau_217_avg_conc ~ 'pTau 217 (pg/mL)'
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("Projects/Brain_Imaging/results/BrainfMRI_data_proteomics_demographics.png", 
         vwidth = 1200, vheight = 800)




small_dat <- small_dat %>% 
  dplyr::select(record_id, group, ab40_avg_conc, ab42_avg_conc, tau_avg_conc, 
                nfl_avg_conc, gfap_avg_conc, ptau_181_avg_conc, ptau_217_avg_conc)



t2d_ids <- small_dat$record_id[which(small_dat$group == 'Type 2 Diabetes')]
lc_ids <- small_dat$record_id[which(small_dat$group == 'Lean Control')]








qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")

# Set your base directory
base_dir <- 'Projects/Brain_Imaging/data/'


# Function to extract regional brain volumes from T1 structural MRI
# This creates a "morphometric feature matrix" for each participant
compute_structural_features_from_t1 <- function(base_dir, participant_id) {
  # Look for .nii.gz files matching this participant ID
  all_files <- list.files(base_dir, pattern = paste0("^", participant_id), 
                          full.names = TRUE, ignore.case = TRUE)
  
  # Filter for .nii.gz files
  nii_files <- all_files[grepl("\\.nii\\.gz$", all_files, ignore.case = TRUE)]
  
  if (length(nii_files) == 0) {
    warning(paste("No .nii.gz file found for participant:", participant_id))
    return(NULL)
  }
  
  tryCatch({
    # Load the NIfTI file
    cat("Loading", basename(nii_files[1]), "...\n")
    nii_data <- readNIfTI(nii_files[1], reorient = FALSE)
    
    # Get dimensions
    dims <- dim(nii_data)
    cat("  Dimensions:", paste(dims, collapse = " x "), "\n")
    
    # Create brain parcellation into ROIs
    # Using intensity-based segmentation and spatial parcellation
    
    # Threshold to create brain mask
    brain_data <- as.vector(nii_data)
    threshold <- quantile(brain_data[brain_data > 0], 0.25, na.rm = TRUE)
    brain_mask <- nii_data > threshold
    
    # Get brain voxel coordinates and intensities
    brain_coords <- which(brain_mask, arr.ind = TRUE)
    brain_intensities <- nii_data[brain_mask]
    n_voxels <- nrow(brain_coords)
    cat("  Brain voxels:", n_voxels, "\n")
    
    # Create spatial parcellation (divide brain into ROIs)
    set.seed(123)
    n_rois <- 100  # Standard parcellation size
    
    if (n_voxels < n_rois * 20) {
      n_rois <- max(50, floor(n_voxels / 50))
      cat("  Adjusted n_rois to:", n_rois, "\n")
    }
    
    # Combine spatial and intensity features for parcellation
    # Normalize coordinates
    coords_norm <- scale(brain_coords)
    intensity_norm <- scale(brain_intensities)
    
    # Weight spatial more than intensity (4:1 ratio)
    features_for_clustering <- cbind(coords_norm * 2, intensity_norm)
    
    # Perform k-means clustering
    cat("  Creating", n_rois, "ROI parcellation...\n")
    kmeans_result <- kmeans(features_for_clustering, 
                            centers = n_rois, 
                            iter.max = 100,
                            nstart = 3)
    parcellation <- kmeans_result$cluster
    
    # Extract morphometric features for each ROI
    roi_features <- matrix(NA, nrow = n_rois, ncol = 4)
    colnames(roi_features) <- c("volume", "mean_intensity", 
                                "intensity_sd", "spatial_extent")
    
    cat("  Extracting ROI features...\n")
    for (roi in 1:n_rois) {
      roi_voxels <- parcellation == roi
      roi_coords <- brain_coords[roi_voxels, , drop = FALSE]
      roi_values <- brain_intensities[roi_voxels]
      
      if (length(roi_values) > 0) {
        # Volume (number of voxels)
        roi_features[roi, "volume"] <- length(roi_values)
        
        # Mean intensity
        roi_features[roi, "mean_intensity"] <- mean(roi_values, na.rm = TRUE)
        
        # Intensity variability
        roi_features[roi, "intensity_sd"] <- sd(roi_values, na.rm = TRUE)
        
        # Spatial extent (max distance between voxels)
        if (nrow(roi_coords) > 1) {
          dists <- dist(roi_coords)
          roi_features[roi, "spatial_extent"] <- max(dists)
        } else {
          roi_features[roi, "spatial_extent"] <- 0
        }
      }
    }
    
    # Remove ROIs with missing data
    valid_rois <- complete.cases(roi_features) & 
      roi_features[, "volume"] > 0
    roi_features <- roi_features[valid_rois, ]
    n_valid_rois <- sum(valid_rois)
    cat("  Valid ROIs:", n_valid_rois, "\n")
    
    # Compute morphometric similarity network
    # This creates a connectivity-like matrix based on structural similarity
    cat("  Computing morphometric similarity matrix...\n")
    
    # Normalize features
    roi_features_norm <- scale(roi_features)
    
    # Compute correlation between ROI feature profiles
    # This represents how similar regions are in their morphometry
    similarity_matrix <- cor(t(roi_features_norm), 
                             use = "pairwise.complete.obs")
    
    # Set diagonal to 0
    diag(similarity_matrix) <- 0
    
    # Convert to absolute values (similarity is non-directional)
    similarity_matrix <- abs(similarity_matrix)
    
    cat("  Morphometric similarity matrix:", nrow(similarity_matrix), "x", 
        ncol(similarity_matrix), "\n\n")
    
    return(similarity_matrix)
    
  }, error = function(e) {
    warning(paste("Error processing", participant_id, ":", e$message))
    return(NULL)
  })
}

# Load structural data for all participants
cat("\n=== Computing Morphometric Similarity Networks ===\n")
cat("Using T1-weighted structural MRI data\n")
cat("This may take several minutes...\n\n")

# Process T2D group
cat("Processing T2D group...\n")
t2d_conn <- lapply(mri_ids_df$file_id[which(mri_ids_df$ID %in% t2d_ids)], 
                   function(id) compute_structural_features_from_t1(base_dir, id))

# Process LC group
cat("\nProcessing LC group...\n")
lc_conn <- lapply(mri_ids_df$file_id[which(mri_ids_df$ID %in% lc_ids)], 
                  function(id) compute_structural_features_from_t1(base_dir, id))

# Remove NULL entries (participants without data)
t2d_valid_idx <- !sapply(t2d_conn, is.null)
lc_valid_idx <- !sapply(lc_conn, is.null)

t2d_conn <- t2d_conn[t2d_valid_idx]
lc_conn <- lc_conn[lc_valid_idx]
t2d_ids_valid <- t2d_ids[t2d_valid_idx]
lc_ids_valid <- lc_ids[lc_valid_idx]

cat("\n=== Summary ===\n")
cat("Loaded", length(t2d_conn), "morphometric networks for T2D group\n")
cat("Loaded", length(lc_conn), "morphometric networks for LC group\n")

# Get number of ROIs
if (length(t2d_conn) > 0) {
  n_rois <- nrow(t2d_conn[[1]])
  cat("Number of ROIs:", n_rois, "\n")
} else if (length(lc_conn) > 0) {
  n_rois <- nrow(lc_conn[[1]])
  cat("Number of ROIs:", n_rois, "\n")
} else {
  stop("No morphometric networks computed! Check .nii.gz files.")
}

# Save morphometric similarity networks
cat("\nSaving morphometric similarity networks...\n")
saveRDS(list(
  t2d_conn = t2d_conn,
  lc_conn = lc_conn,
  t2d_ids = t2d_ids_valid,
  lc_ids = lc_ids_valid,
  note = "These are morphometric similarity networks from T1 structural MRI"
), "morphometric_similarity_networks.rds")

cat("Networks saved to: morphometric_similarity_networks.rds\n")
cat("\nNote: These are structural similarity networks, not functional connectivity.\n")
cat("They represent similarity in brain morphometry between regions.\n")






















#### Analysis for 3d Data


library(ggplot2)
library(dplyr)
library(corrplot)
library(tidyr)
library(ggsignif)
library(purrr)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
 
results.dir <- 'Projects/Brain_Imaging/results/'


cat("=== Loading Data ===\n")

# Load connectivity matrices
conn_data <- readRDS(paste0(results.dir, "morphometric_similarity_networks.rds"))
t2d_conn <- conn_data$t2d_conn
lc_conn <- conn_data$lc_conn
t2d_ids_valid <- conn_data$t2d_ids
lc_ids_valid <- conn_data$lc_ids

cat("T2D participants:", length(t2d_conn), "\n")
cat("LC participants:", length(lc_conn), "\n")
cat("ROIs per network:", nrow(t2d_conn[[1]]), "\n\n")

# Load proteomics data (adjust path as needed)

harmonized_data <- read.csv("../OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


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


small_dat <- small_dat %>% dplyr::select(record_id, group, all_of(qx_var))


# If you need to reload: small_dat <- readRDS("path_to_your_proteomics_data.rds")

# ==============================================================================
# 2. EXTRACT NETWORK METRICS FOR EACH PARTICIPANT
# ==============================================================================

cat("=== Extracting Network Metrics ===\n")

# Function to calculate comprehensive network metrics
calculate_network_metrics <- function(conn_matrix) {
  # Set diagonal to NA for calculations
  diag(conn_matrix) <- NA
  
  # Global metrics
  mean_connectivity <- mean(conn_matrix, na.rm = TRUE)
  median_connectivity <- median(conn_matrix, na.rm = TRUE)
  sd_connectivity <- sd(conn_matrix, na.rm = TRUE)
  
  # Network strength (mean absolute connectivity)
  mean_strength <- mean(abs(conn_matrix), na.rm = TRUE)
  
  # Network density (proportion of strong connections)
  threshold_strong <- 0.3
  density_strong <- sum(abs(conn_matrix) > threshold_strong, na.rm = TRUE) / 
    sum(!is.na(conn_matrix))
  
  # Network efficiency (simplified global efficiency)
  # Higher values = more efficient network
  pos_conn <- abs(conn_matrix)
  pos_conn[pos_conn == 0] <- NA
  mean_efficiency <- mean(pos_conn, na.rm = TRUE)
  
  # Modularity-related: variance in connectivity
  # Higher variance suggests more modular structure
  connectivity_variance <- var(as.vector(conn_matrix), na.rm = TRUE)
  
  # Average node strength (sum of connections per ROI)
  node_strengths <- rowSums(abs(conn_matrix), na.rm = TRUE)
  mean_node_strength <- mean(node_strengths, na.rm = TRUE)
  sd_node_strength <- sd(node_strengths, na.rm = TRUE)
  
  # Positive vs negative connections
  prop_positive <- sum(conn_matrix > 0, na.rm = TRUE) / sum(!is.na(conn_matrix))
  
  return(data.frame(
    mean_connectivity = mean_connectivity,
    median_connectivity = median_connectivity,
    sd_connectivity = sd_connectivity,
    mean_strength = mean_strength,
    network_density = density_strong,
    mean_efficiency = mean_efficiency,
    connectivity_variance = connectivity_variance,
    mean_node_strength = mean_node_strength,
    sd_node_strength = sd_node_strength,
    prop_positive = prop_positive
  ))
}

# Calculate metrics for all participants
t2d_metrics <- do.call(rbind, lapply(t2d_conn, calculate_network_metrics))
lc_metrics <- do.call(rbind, lapply(lc_conn, calculate_network_metrics))

# Add IDs and group labels
t2d_metrics$record_id <- t2d_ids_valid
t2d_metrics$group <- "T2D"

lc_metrics$record_id <- lc_ids_valid
lc_metrics$group <- "LC"

# Combine
all_metrics <- rbind(t2d_metrics, lc_metrics)

cat("Network metrics calculated for", nrow(all_metrics), "participants\n\n")

# ==============================================================================
# 3. GROUP COMPARISONS: STATISTICAL TESTS
# ==============================================================================
cat("=== Testing Group Differences ===\n")

# Define brain features
brain_features <- c("mean_connectivity", "median_connectivity", "sd_connectivity",
                    "mean_strength", "network_density", "mean_efficiency",
                    "connectivity_variance", "mean_node_strength", 
                    "sd_node_strength", "prop_positive")

# Perform t-tests
comparison_results <- data.frame(
  Metric = brain_features,
  T2D_mean = NA,
  T2D_sd = NA,
  LC_mean = NA,
  LC_sd = NA,
  t_stat = NA,
  p_value = NA,
  cohens_d = NA,
  interpretation = NA
)

for (i in 1:length(brain_features)) {
  metric <- brain_features[i]
  
  t2d_vals <- all_metrics[[metric]][all_metrics$group == "T2D"]
  lc_vals <- all_metrics[[metric]][all_metrics$group == "LC"]
  
  # Remove NA values
  t2d_vals <- t2d_vals[!is.na(t2d_vals)]
  lc_vals <- lc_vals[!is.na(lc_vals)]
  
  comparison_results$T2D_mean[i] <- mean(t2d_vals, na.rm = TRUE)
  comparison_results$T2D_sd[i] <- sd(t2d_vals, na.rm = TRUE)
  comparison_results$LC_mean[i] <- mean(lc_vals, na.rm = TRUE)
  comparison_results$LC_sd[i] <- sd(lc_vals, na.rm = TRUE)
  
  # Check if we have enough variance to perform t-test
  t2d_var <- var(t2d_vals, na.rm = TRUE)
  lc_var <- var(lc_vals, na.rm = TRUE)
  
  if (length(t2d_vals) < 2 || length(lc_vals) < 2) {
    cat("  Warning: Insufficient sample size for", metric, "\n")
    comparison_results$t_stat[i] <- NA
    comparison_results$p_value[i] <- NA
    comparison_results$cohens_d[i] <- NA
    comparison_results$interpretation[i] <- "Insufficient data"
    next
  }
  
  if (is.na(t2d_var) || is.na(lc_var) || t2d_var < 1e-10 || lc_var < 1e-10) {
    cat("  Warning: No variance in", metric, "(T2D var:", t2d_var, ", LC var:", lc_var, ")\n")
    comparison_results$t_stat[i] <- NA
    comparison_results$p_value[i] <- NA
    comparison_results$cohens_d[i] <- NA
    comparison_results$interpretation[i] <- "No variance"
    next
  }
  
  # T-test with error handling
  tryCatch({
    test <- t.test(t2d_vals, lc_vals)
    comparison_results$t_stat[i] <- test$statistic
    comparison_results$p_value[i] <- test$p.value
    
    # Cohen's d effect size
    pooled_sd <- sqrt(((length(t2d_vals)-1)*var(t2d_vals, na.rm=TRUE) + 
                         (length(lc_vals)-1)*var(lc_vals, na.rm=TRUE)) / 
                        (length(t2d_vals) + length(lc_vals) - 2))
    
    if (pooled_sd > 0) {
      cohens_d <- (comparison_results$T2D_mean[i] - comparison_results$LC_mean[i]) / pooled_sd
      comparison_results$cohens_d[i] <- cohens_d
      
      # Effect size interpretation
      abs_d <- abs(cohens_d)
      if (abs_d < 0.2) {
        comparison_results$interpretation[i] <- "Negligible"
      } else if (abs_d < 0.5) {
        comparison_results$interpretation[i] <- "Small"
      } else if (abs_d < 0.8) {
        comparison_results$interpretation[i] <- "Medium"
      } else {
        comparison_results$interpretation[i] <- "Large"
      }
    } else {
      comparison_results$cohens_d[i] <- NA
      comparison_results$interpretation[i] <- "Cannot compute"
    }
  }, error = function(e) {
    cat("  Error in t-test for", metric, ":", e$message, "\n")
    comparison_results$t_stat[i] <- NA
    comparison_results$p_value[i] <- NA
    comparison_results$cohens_d[i] <- NA
    comparison_results$interpretation[i] <- "Error"
  })
}

# FDR correction
comparison_results$p_fdr <- p.adjust(comparison_results$p_value, method = "fdr")

# Print results
cat("\nGroup Comparison Results:\n")
print(comparison_results %>% 
        select(Metric, T2D_mean, LC_mean, t_stat, p_value, p_fdr, cohens_d, interpretation) %>%
        mutate(across(where(is.numeric), ~round(., 4))))

# Save results
write.csv(comparison_results, paste0(results.dir, "brain_metrics_group_comparison.csv"), row.names = FALSE)
cat("\nResults saved to: brain_metrics_group_comparison.csv\n")

# ==============================================================================
# 4. VISUALIZE GROUP DIFFERENCES
# ==============================================================================


# Prepare data for plotting
plot_data <- all_metrics %>%
  pivot_longer(cols = all_of(brain_features),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = brain_features))

# Create comparison plot
p1 <- ggplot(plot_data, aes(x = group, y = Value, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("T2D" = "#E74C3C", "LC" = "#3498DB")) +
  labs(title = "Brain Network Metrics: T2D vs Lean Control",
       subtitle = "Morphometric Similarity Networks",
       x = "Group", 
       y = "Value",
       fill = "Group") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0(results.dir, "brain_metrics_comparison.pdf"), p1, width = 14, height = 10)
ggsave(paste0(results.dir, "brain_metrics_comparison.png"), p1, width = 14, height = 10, dpi = 300)
cat("Saved: brain_metrics_comparison.pdf/.png\n")

# Create forest plot for effect sizes
effect_plot_data <- comparison_results %>%
  mutate(
    Metric = factor(Metric, levels = rev(brain_features)),
    Significant = ifelse(p_fdr < 0.05, "FDR < 0.05", "Not Sig"),
    ci_lower = cohens_d - 1.96 * sqrt((nrow(t2d_metrics) + nrow(lc_metrics)) / 
                                        (nrow(t2d_metrics) * nrow(lc_metrics))),
    ci_upper = cohens_d + 1.96 * sqrt((nrow(t2d_metrics) + nrow(lc_metrics)) / 
                                        (nrow(t2d_metrics) * nrow(lc_metrics)))
  )

p2 <- ggplot(effect_plot_data, aes(x = cohens_d, y = Metric, color = Significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray70") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = c("FDR < 0.05" = "#E74C3C", "Not Sig" = "gray60")) +
  labs(title = "Effect Sizes: T2D vs LC Brain Network Metrics",
       subtitle = "Cohen's d with 95% CI",
       x = "Cohen's d (T2D - LC)",
       y = NULL,
       color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

ggsave(paste0(results.dir, "effect_sizes_forest_plot.pdf"), p2, width = 10, height = 8)
ggsave(paste0(results.dir, "effect_sizes_forest_plot.png"), p2, width = 10, height = 8, dpi = 300)
cat("Saved: effect_sizes_forest_plot.pdf/.png\n")

# ==============================================================================
# 5. INTEGRATE WITH PROTEOMICS DATA
# ==============================================================================

cat("\n=== Integrating with Proteomics ===\n")

# Merge brain metrics with proteomics
integrated_data <- small_dat %>%
  inner_join(all_metrics, by = c("record_id"))

names(integrated_data)[20] <- 'group'
cat("Integrated dataset: n =", nrow(integrated_data), "\n")
cat("  T2D:", sum(integrated_data$group == "T2D"), "\n")
cat("  LC:", sum(integrated_data$group == "LC"), "\n\n")

# Define protein markers
protein_markers <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc", 
                     "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", 
                     "ptau_217_avg_conc")

# Check which proteins have sufficient data
proteins_available <- protein_markers[protein_markers %in% names(integrated_data)]
cat("Protein markers available:", length(proteins_available), "\n")
print(proteins_available)

# ==============================================================================
# 6. CORRELATION ANALYSIS: BRAIN METRICS vs PROTEOMICS
# ==============================================================================

cat("\n=== Brain-Proteomics Correlations ===\n")

# Create correlation matrix
cor_data <- integrated_data %>%
  select(all_of(c(brain_features, proteins_available))) %>%
  na.omit()

cat("Complete cases for correlation:", nrow(cor_data), "\n")

if (nrow(cor_data) > 5) {
  # Compute correlation matrix
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs", method = "pearson")
  
  # Extract brain-protein correlations
  brain_protein_cors <- cor_matrix[brain_features, proteins_available]
  
  # Visualize
  pdf("brain_protein_correlations.pdf", width = 10, height = 8)
  corrplot(brain_protein_cors, 
           method = "color",
           type = "full",
           tl.col = "black", 
           tl.srt = 45,
           tl.cex = 0.9,
           addCoef.col = "black", 
           number.cex = 0.7,
           col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200),
           title = "Brain Network Metrics vs Proteomic Markers",
           mar = c(0,0,2,0))
  dev.off()
  cat("Saved: brain_protein_correlations.pdf\n")
  
  # Statistical testing for each correlation
  correlation_results <- expand.grid(
    Brain_Feature = brain_features,
    Protein = proteins_available,
    stringsAsFactors = FALSE
  )
  
  correlation_results$Correlation <- NA
  correlation_results$P_value <- NA
  correlation_results$N <- NA
  correlation_results$CI_lower <- NA
  correlation_results$CI_upper <- NA
  
  for (i in 1:nrow(correlation_results)) {
    bf <- correlation_results$Brain_Feature[i]
    pm <- correlation_results$Protein[i]
    
    complete_cases <- complete.cases(integrated_data[[bf]], integrated_data[[pm]])
    
    if (sum(complete_cases) > 3) {
      x <- integrated_data[[bf]][complete_cases]
      y <- integrated_data[[pm]][complete_cases]
      
      cor_test <- cor.test(x, y, method = "pearson")
      
      correlation_results$Correlation[i] <- cor_test$estimate
      correlation_results$P_value[i] <- cor_test$p.value
      correlation_results$N[i] <- sum(complete_cases)
      correlation_results$CI_lower[i] <- cor_test$conf.int[1]
      correlation_results$CI_upper[i] <- cor_test$conf.int[2]
    }
  }
  
  # FDR correction
  correlation_results$P_fdr <- p.adjust(correlation_results$P_value, method = "fdr")
  
  # Add significance flag
  correlation_results$Significant <- ifelse(correlation_results$P_fdr < 0.05, 
                                            "Yes", "No")
  
  # Show significant correlations
  significant_cors <- correlation_results %>%
    filter(P_fdr < 0.05) %>%
    arrange(P_fdr) %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  
  cat("\nSignificant Brain-Protein Correlations (FDR < 0.05):\n")
  if (nrow(significant_cors) > 0) {
    print(significant_cors)
  } else {
    cat("No significant correlations after FDR correction.\n")
  }
  
  # Show top correlations regardless of significance
  top_cors <- correlation_results %>%
    arrange(P_value) %>%
    head(10) %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  
  cat("\nTop 10 Correlations (by p-value):\n")
  print(top_cors)
  
  # Save results
  write.csv(correlation_results, paste0(results.dir, "brain_protein_correlations.csv"), row.names = FALSE)
  cat("\nCorrelation results saved to: brain_protein_correlations.csv\n")
  
  # Create heatmap of p-values
  p_matrix <- matrix(correlation_results$P_value, 
                     nrow = length(brain_features),
                     ncol = length(proteins_available))
  rownames(p_matrix) <- brain_features
  colnames(p_matrix) <- proteins_available
  
  pdf(paste0(results.dir, "brain_protein_pvalues.pdf"), width = 10, height = 8)
  corrplot(-log10(p_matrix), 
           method = "color",
           is.corr = FALSE,
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.9,
           addCoef.col = "black",
           number.cex = 0.7,
           col = colorRampPalette(c("white", "#E74C3C"))(200),
           title = "Brain-Protein Associations (-log10 p-value)",
           mar = c(0,0,2,0))
  dev.off()
  cat("Saved: brain_protein_pvalues.pdf\n")
}

# ==============================================================================
# 7. GROUP-SPECIFIC CORRELATIONS
# ==============================================================================

cat("\n=== Group-Specific Correlations ===\n")

# Function to compute correlations within a group
group_correlations <- function(data, group_name, brain_features, protein_markers) {
  group_data <- data %>% filter(group == group_name)
  
  results <- expand.grid(
    Brain_Feature = brain_features,
    Protein = protein_markers,
    stringsAsFactors = FALSE
  )
  
  results$Correlation <- NA
  results$P_value <- NA
  results$N <- NA
  results$Group <- group_name
  
  for (i in 1:nrow(results)) {
    bf <- results$Brain_Feature[i]
    pm <- results$Protein[i]
    
    if (bf %in% names(group_data) && pm %in% names(group_data)) {
      complete_cases <- complete.cases(group_data[[bf]], group_data[[pm]])
      
      if (sum(complete_cases) > 3) {
        cor_test <- cor.test(group_data[[bf]][complete_cases], 
                             group_data[[pm]][complete_cases], 
                             method = "pearson")
        results$Correlation[i] <- cor_test$estimate
        results$P_value[i] <- cor_test$p.value
        results$N[i] <- sum(complete_cases)
      }
    }
  }
  
  return(results)
}

# Compute group-specific correlations
t2d_cors <- group_correlations(integrated_data, "T2D", brain_features, proteins_available)
lc_cors <- group_correlations(integrated_data, "LC", brain_features, proteins_available)

# Combine
all_group_cors <- rbind(t2d_cors, lc_cors)

# FDR correction within each group
all_group_cors <- all_group_cors %>%
  group_by(Group) %>%
  mutate(P_fdr = p.adjust(P_value, method = "fdr")) %>%
  ungroup()

# Find significant group-specific correlations
sig_group_cors <- all_group_cors %>%
  filter(P_fdr < 0.05) %>%
  arrange(Group, P_fdr)

cat("\nSignificant group-specific correlations (FDR < 0.05):\n")
if (nrow(sig_group_cors) > 0) {
  print(sig_group_cors %>% mutate(across(where(is.numeric), ~round(., 4))))
} else {
  cat("No significant group-specific correlations after FDR correction.\n")
}

# Save group-specific results
write.csv(all_group_cors, paste0(results.dir, "group_specific_correlations.csv"), row.names = FALSE)
cat("\nGroup-specific correlations saved to: group_specific_correlations.csv\n")

# ==============================================================================
# 8. CREATE SCATTER PLOTS FOR TOP ASSOCIATIONS
# ==============================================================================

cat("\n=== Creating Scatter Plots ===\n")

# Get top 6 correlations for plotting
top_assoc <- correlation_results %>%
  filter(!is.na(Correlation)) %>%
  arrange(P_value) %>%
  head(6)

if (nrow(top_assoc) > 0) {
  plot_list <- list()
  
  for (i in 1:nrow(top_assoc)) {
    bf <- top_assoc$Brain_Feature[i]
    pm <- top_assoc$Protein[i]
    r <- round(top_assoc$Correlation[i], 3)
    p <- top_assoc$P_value[i]
    p_fdr <- top_assoc$P_fdr[i]
    
    # Create nice labels
    bf_label <- gsub("_", " ", bf)
    bf_label <- tools::toTitleCase(bf_label)
    pm_label <- gsub("_avg_conc", "", pm)
    pm_label <- toupper(gsub("_", "", pm_label))
    
    subtitle <- sprintf("r = %.3f, p = %.4f, FDR = %.4f", r, p, p_fdr)
    
    p_scatter <- ggplot(integrated_data, 
                        aes_string(x = bf, y = pm, color = "group")) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, aes(group = 1), 
                  color = "black", linetype = "dashed") +
      scale_color_manual(values = c("T2D" = "#E74C3C", "LC" = "#3498DB")) +
      labs(title = paste(bf_label, "vs", pm_label),
           subtitle = subtitle,
           x = bf_label,
           y = paste(pm_label, "(pg/mL)"),
           color = "Group") +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      )
    
    plot_list[[i]] <- p_scatter
  }
  
  # Combine plots
  library(gridExtra)
  combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))
  
  ggsave(paste0(results.dir, "top_brain_protein_associations.pdf"), combined_plot, 
         width = 15, height = 10)
  ggsave(paste0(results.dir, "top_brain_protein_associations.png"), combined_plot, 
         width = 15, height = 10, dpi = 300)
  cat("Saved: top_brain_protein_associations.pdf/.png\n")
}

# ==============================================================================
# 9. SUMMARY REPORT
# ==============================================================================

cat("\n" , "="*80, "\n")
cat("ANALYSIS COMPLETE!\n")
cat("="*80, "\n\n")

cat("Files created:\n")
cat("1. brain_metrics_group_comparison.csv - Statistical comparison of groups\n")
cat("2. brain_metrics_comparison.pdf/.png - Boxplots of all metrics\n")
cat("3. effect_sizes_forest_plot.pdf/.png - Cohen's d effect sizes\n")
cat("4. brain_protein_correlations.csv - All correlation results\n")
cat("5. brain_protein_correlations.pdf - Correlation heatmap\n")
cat("6. brain_protein_pvalues.pdf - P-value heatmap\n")
cat("7. group_specific_correlations.csv - Correlations within T2D and LC\n")
cat("8. top_brain_protein_associations.pdf/.png - Scatter plots\n\n")

cat("Key Findings:\n")
cat("- Significant group differences:", 
    sum(comparison_results$p_fdr < 0.05, na.rm = TRUE), "metrics\n")
cat("- Significant brain-protein correlations:", 
    sum(correlation_results$P_fdr < 0.05, na.rm = TRUE), "associations\n")
cat("- Effect sizes: Large =", sum(abs(comparison_results$cohens_d) > 0.8, na.rm = TRUE),
    ", Medium =", sum(abs(comparison_results$cohens_d) >= 0.5 & 
                        abs(comparison_results$cohens_d) < 0.8, na.rm = TRUE), "\n")








#### FreeSurfer Analysis


# Load FreeSurfer output
thickness_data <- read.table("lh_thickness.txt", header=TRUE)
volume_data <- read.table("volumes.txt", header=TRUE)

# Merge with clinical data
merged <- clinical_data %>%
  left_join(thickness_data, by="subject_id") %>%
  left_join(volume_data, by="subject_id")

# Standard regions to examine in T2D studies:
regions_of_interest <- c(
  "Left.Hippocampus", "Right.Hippocampus",
  "Left.Amygdala", "Right.Amygdala",
  "lh_entorhinal_thickness", "rh_entorhinal_thickness",
  "lh_precuneus_thickness", "rh_precuneus_thickness",
  "TotalGrayVol", "TotalCorticalVol"
)


















