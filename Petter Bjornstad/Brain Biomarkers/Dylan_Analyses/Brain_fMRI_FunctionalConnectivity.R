# ==============================================================================
# FUNCTIONAL CONNECTIVITY ANALYSIS FROM RESTING-STATE fMRI
# ==============================================================================

library(oro.nifti)
library(neurobase)
library(ggplot2)
library(dplyr)
library(corrplot)
library(tidyr)

# ==============================================================================
# 1. SETUP AND DATA LOADING
# ==============================================================================



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


small_dat$group[which(small_dat$record_id == 'RH2-38-T')] <- 'Obese Control'
small_dat$record_id[which(small_dat$record_id == 'RH2-38-T')] <- 'RH2-38-O'





qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")





small_dat <- small_dat %>% 
  dplyr::select(record_id, group, ab40_avg_conc, ab42_avg_conc, tau_avg_conc, 
                nfl_avg_conc, gfap_avg_conc, ptau_181_avg_conc, ptau_217_avg_conc)



t2d_ids <- small_dat$record_id[which(small_dat$group == 'Type 2 Diabetes')]
lc_ids <- small_dat$record_id[which(small_dat$group == 'Lean Control')]








qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")


cat("=== Functional Connectivity Analysis ===\n")
cat("Using resting-state fMRI data (rest_brain.nii.gz)\n\n")

base_dir <- '~/Documents/UofW/Projects/Brain_Imaging/fMRI/data/'

# Load your participant IDs and groups (from your earlier script)
# Adjust these based on your actual data structure

# Function to find rest_brain.nii.gz for a participant
find_fmri_file <- function(base_dir, participant_id) {
  # Look for participant folder
  participant_folder <- file.path(base_dir, participant_id)
  
  all_files <- list.files(base_dir, pattern = paste0("^", participant_id), 
                          full.names = TRUE, ignore.case = TRUE)
  
  # Filter for .nii.gz files
  fmri_file <- all_files[grepl("rest_brain\\.nii\\.gz$", all_files, ignore.case = TRUE)]
  
  return(fmri_file)
}

# ==============================================================================
# 2. COMPUTE FUNCTIONAL CONNECTIVITY FROM fMRI
# ==============================================================================

compute_functional_connectivity <- function(fmri_file, participant_id) {
  
  tryCatch({
    cat("\n=== Processing:", participant_id, "===\n")
    cat("Loading fMRI data...\n")
    
    # Load 4D fMRI data (x, y, z, time)
    fmri_data <- readNIfTI(fmri_file, reorient = FALSE)
    dims <- dim(fmri_data)
    
    cat("  Dimensions:", paste(dims, collapse = " x "), "\n")
    
    # Verify it's 4D
    if (length(dims) != 4) {
      warning(paste("Not 4D data for:", participant_id))
      return(NULL)
    }
    
    n_timepoints <- dims[4]
    cat("  Time points:", n_timepoints, "\n")
    
    # Create brain mask (voxels with consistent signal)
    # Calculate temporal mean and standard deviation
    mean_img <- apply(fmri_data, c(1,2,3), mean, na.rm = TRUE)
    sd_img <- apply(fmri_data, c(1,2,3), sd, na.rm = TRUE)
    
    # Mask: voxels with signal and temporal variation
    threshold_mean <- quantile(mean_img[mean_img > 0], 0.25, na.rm = TRUE)
    threshold_sd <- quantile(sd_img[sd_img > 0], 0.1, na.rm = TRUE)
    
    brain_mask <- (mean_img > threshold_mean) & (sd_img > threshold_sd)
    
    # Get brain voxel coordinates
    brain_coords <- which(brain_mask, arr.ind = TRUE)
    n_voxels <- nrow(brain_coords)
    cat("  Brain voxels:", n_voxels, "\n")
    
    if (n_voxels < 1000) {
      warning(paste("Too few brain voxels for:", participant_id))
      return(NULL)
    }
    
    # ==============================================================================
    # 3. CREATE BRAIN PARCELLATION (ROIS)
    # ==============================================================================
    
    cat("  Creating brain parcellation...\n")
    
    # Define number of ROIs (standard: 100-200 for whole brain)
    n_rois <- 100
    
    # Adjust if needed
    if (n_voxels < n_rois * 20) {
      n_rois <- max(50, floor(n_voxels / 50))
      cat("  Adjusted to", n_rois, "ROIs\n")
    }
    
    # Cluster voxels into ROIs using k-means
    # Use spatial coordinates weighted by intensity
    set.seed(123)  # For reproducibility
    
    # Normalize coordinates for clustering
    coords_norm <- scale(brain_coords)
    
    # Perform k-means clustering
    kmeans_result <- kmeans(coords_norm, 
                            centers = n_rois, 
                            iter.max = 100,
                            nstart = 3)
    
    parcellation <- kmeans_result$cluster
    
    # ==============================================================================
    # 4. EXTRACT ROI TIME SERIES
    # ==============================================================================
    
    cat("  Extracting ROI time series...\n")
    
    roi_timeseries <- matrix(NA, nrow = n_timepoints, ncol = n_rois)
    
    for (roi in 1:n_rois) {
      # Get voxels in this ROI
      roi_idx <- which(parcellation == roi)
      
      if (length(roi_idx) > 0) {
        roi_voxels <- brain_coords[roi_idx, , drop = FALSE]
        
        # Extract time series for all voxels in this ROI
        roi_signals <- matrix(NA, nrow = n_timepoints, ncol = nrow(roi_voxels))
        
        for (v in 1:nrow(roi_voxels)) {
          x <- roi_voxels[v, 1]
          y <- roi_voxels[v, 2]
          z <- roi_voxels[v, 3]
          
          # Extract time series for this voxel
          roi_signals[, v] <- fmri_data[x, y, z, ]
        }
        
        # Average time series across all voxels in the ROI
        roi_timeseries[, roi] <- rowMeans(roi_signals, na.rm = TRUE)
      }
    }
    
    # ==============================================================================
    # 5. PREPROCESSING OF TIME SERIES
    # ==============================================================================
    
    cat("  Preprocessing time series...\n")
    
    # Remove linear trend (detrending)
    for (roi in 1:n_rois) {
      if (!all(is.na(roi_timeseries[, roi]))) {
        time_points <- 1:n_timepoints
        lm_fit <- lm(roi_timeseries[, roi] ~ time_points)
        roi_timeseries[, roi] <- residuals(lm_fit)
      }
    }
    
    # Temporal filtering (optional but recommended)
    # Remove very low frequency drift (< 0.01 Hz) and high frequency noise (> 0.1 Hz)
    # This requires knowing your TR (repetition time)
    # For now, we'll skip this step - add if you know your TR
    
    # Remove ROIs with no variance or too many NAs
    valid_rois <- apply(roi_timeseries, 2, function(x) {
      !all(is.na(x)) && 
        sum(!is.na(x)) > (n_timepoints * 0.8) &&
        var(x, na.rm = TRUE) > 1e-10
    })
    
    roi_timeseries <- roi_timeseries[, valid_rois]
    n_valid_rois <- sum(valid_rois)
    cat("  Valid ROIs:", n_valid_rois, "\n")
    
    if (n_valid_rois < 50) {
      warning(paste("Too few valid ROIs for:", participant_id))
      return(NULL)
    }
    
    # ==============================================================================
    # 6. COMPUTE FUNCTIONAL CONNECTIVITY MATRIX
    # ==============================================================================
    
    cat("  Computing functional connectivity matrix...\n")
    
    # Pearson correlation between all pairs of ROI time series
    fc_matrix <- cor(roi_timeseries, use = "pairwise.complete.obs")
    
    # Replace any NaN or Inf with 0
    fc_matrix[is.na(fc_matrix)] <- 0
    fc_matrix[is.infinite(fc_matrix)] <- 0
    
    # Set diagonal to 0 (no self-connections)
    diag(fc_matrix) <- 0
    
    cat("  Final FC matrix:", nrow(fc_matrix), "x", ncol(fc_matrix), "\n")
    cat("  Mean connectivity:", round(mean(fc_matrix, na.rm = TRUE), 3), "\n")
    cat("  SD connectivity:", round(sd(fc_matrix, na.rm = TRUE), 3), "\n")
    
    return(fc_matrix)
    
  }, error = function(e) {
    cat("\n  ERROR processing", participant_id, ":", e$message, "\n")
    return(NULL)
  })
}

# ==============================================================================
# 7. PROCESS ALL PARTICIPANTS
# ==============================================================================

# Find all participant directories with rest_brain.nii.gz
all_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

# Filter for directories that have rest_brain.nii.gz
participants_with_fmri <-  list.files(base_dir, full.names = TRUE, ignore.case = TRUE)

cat("Found", length(participants_with_fmri), "participants with fMRI data\n")
cat("Participants:", paste(head(participants_with_fmri, 10), collapse = ", "), "...\n\n")


all_fc_list <- list()
for (i in c(1:length(participants_with_fmri))) {

  fmri_file <- participants_with_fmri[i]
  participant_id <- basename(fmri_file)  %>%
    str_remove("_rest_brain\\.nii\\.gz$")
  
  fc_matrix <- compute_functional_connectivity(fmri_file, participant_id)
  if (!is.null(fc_matrix)) {
    all_fc_list[[participant_id]] <- fc_matrix
  }
}


# ==============================================================================
# 8. SUMMARY AND SAVE RESULTS
# ==============================================================================

# Save connectivity matrices
output_file <- "~/Documents/UofW/Projects/Brain_Imaging/fMRI/results/functional_connectivity_matrices.rds"

saveRDS(list(
  all_fc = all_fc_list,
  all_ids = names(all_fc_list),
  note = "Functional connectivity from resting-state fMRI (rest_brain.nii.gz)",
  date_created = Sys.time()
), output_file)

cat("\nResults saved to:", output_file, "\n")

# ==============================================================================
# 9. QUICK VISUALIZATION
# ==============================================================================

if (length(t2d_fc_list) > 0) {
  cat("\nCreating example visualization...\n")
  
  # Average connectivity matrix for T2D group
  t2d_fc_mean <- Reduce("+", t2d_fc_list) / length(t2d_fc_list)
  
  # Average connectivity matrix for LC group
  if (length(lc_fc_list) > 0) {
    lc_fc_mean <- Reduce("+", lc_fc_list) / length(lc_fc_list)
    
    # Plot group averages
    pdf("~/Documents/UofW/Projects/Brain_Imaging/fc_group_averages.pdf", 
        width = 12, height = 5)
    par(mfrow = c(1, 2))
    
    corrplot(t2d_fc_mean, method = "color", 
             title = "T2D Group Average FC",
             tl.pos = "n", mar = c(0,0,2,0))
    
    corrplot(lc_fc_mean, method = "color", 
             title = "LC Group Average FC",
             tl.pos = "n", mar = c(0,0,2,0))
    
    dev.off()
    
    cat("Saved visualization to: fc_group_averages.pdf\n")
  }
}

cat("\n=== Next Steps ===\n")
cat("1. Run the connectivity_proteomics_analysis.R script\n")
cat("2. This will use these FC matrices to:\n")
cat("   - Extract network metrics\n")
cat("   - Compare T2D vs LC groups\n")
cat("   - Correlate with proteomics\n")
cat("\nAnalysis complete!\n")