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