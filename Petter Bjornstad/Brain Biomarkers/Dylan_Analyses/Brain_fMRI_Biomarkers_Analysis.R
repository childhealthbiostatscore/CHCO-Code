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



#Identifying groups for analysis 

harmonized_data <- read.csv("/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))



mri_ids <- c('CRC-10', 'CRC-11', 'CRC-12', 'CRC-13', 'CRC-26', 'CRC-39', 'CRC-46', 'CRC-51', 'CRC-53', 
             'CRC-55', 'CRC-58', 'CRC-60', 
             'RH2-01-O', 'RH2-03-O', 'RH2-08-T', 'RH2-10-L', 'RH2-11-O', 'RH2-13-O', 'RH2-16-O', 'RH2-17-L', 
             'RH2-18-O', 'RH2-19-T', 'RH2-22-T', 'RH2-24-L', 'RH2-27-L', 'RH2-28-L', 'RH2-29-L', 'RH2-33-L', 
             'RH2-34-O', 'RH2-35-T', 'RH2-38-O', 'RH2-39-O', 'RH2-41-T', 'RH2-42-T', 'RH2-43-T', 
             'RH2-44-T', 'RH2-45-T', 'RH2-48-T', 'RH2-49-T', 'RH2-50-L', 'RH2-52-T', 'RH2-53-T', 
             'RH2-55-T')








qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")


# Set your base directory
base_dir <- "/brain.mir/Brain_Functional_Connectivity_Analysis/MRI_preprocess/Analysis"

# Function to get all participant folders
get_participant_folders <- function(base_dir) {
  all_folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  # Filter for folders that contain connectivity data
  participant_folders <- all_folders[grepl("conn_project|RH2_", basename(all_folders))]
  return(participant_folders)
}

# Function to load connectivity matrix from .mat file
load_connectivity_from_mat <- function(participant_folder) {
  # Look for .mat files in the folder
  mat_files <- list.files(participant_folder, pattern = "\\.mat$", full.names = TRUE)
  
  if (length(mat_files) > 0) {
    # Load the first .mat file found
    mat_data <- readMat(mat_files[1])
    # Extract connectivity matrix (adjust field name based on your .mat structure)
    # Common field names: 'conn', 'connectivity', 'corr', 'Z', 'data'
    conn_matrix <- mat_data[[1]]  # Adjust this based on your .mat file structure
    return(conn_matrix)
  } else {
    warning(paste("No .mat file found in", participant_folder))
    return(NULL)
  }
}

# Function to load preprocessed fMRI time series from .nii.gz
load_fmri_timeseries <- function(participant_folder, roi_atlas = NULL) {
  # Look for .nii.gz files
  nii_files <- list.files(participant_folder, pattern = "\\.nii\\.gz$", 
                          full.names = TRUE, recursive = TRUE)
  
  if (length(nii_files) == 0) {
    warning(paste("No .nii.gz files found in", participant_folder))
    return(NULL)
  }
  
  # Load the first 4D fMRI file
  fmri_img <- readNIfTI(nii_files[1])
  
  if (is.null(roi_atlas)) {
    # If no atlas provided, return the full 4D data
    return(fmri_img)
  }
  
  # Extract time series from ROIs
  n_timepoints <- dim(fmri_img)[4]
  n_rois <- max(roi_atlas, na.rm = TRUE)
  
  timeseries_matrix <- matrix(NA, nrow = n_timepoints, ncol = n_rois)
  
  for (roi in 1:n_rois) {
    roi_mask <- roi_atlas == roi
    roi_voxels <- fmri_img[roi_mask, ]
    timeseries_matrix[, roi] <- apply(roi_voxels, 2, mean, na.rm = TRUE)
  }
  
  return(timeseries_matrix)
}

# Get all participant folders
all_participants <- get_participant_folders(base_dir)

# Define your groups (you'll need to specify which participants belong to which group)
# Option 1: Based on folder naming convention
group1_pattern <- "conn_project01|conn_project02"  # Adjust pattern
group2_pattern <- "conn_project03|RH2_"            # Adjust pattern

group1_folders <- all_participants[grepl(group1_pattern, all_participants)]
group2_folders <- all_participants[grepl(group2_pattern, all_participants)]

# Option 2: Load from a CSV file with participant IDs and group labels
# participant_info <- read.csv("participant_groups.csv")
# group1_folders <- all_participants[basename(all_participants) %in% 
#                                    participant_info$ID[participant_info$Group == 1]]
# group2_folders <- all_participants[basename(all_participants) %in% 
#                                    participant_info$ID[participant_info$Group == 2]]

cat("Found", length(group1_folders), "participants in Group 1\n")
cat("Found", length(group2_folders), "participants in Group 2\n")

# Load connectivity matrices if already computed
group1_conn <- lapply(group1_folders, load_connectivity_from_mat)
group2_conn <- lapply(group2_folders, load_connectivity_from_mat)

# Remove NULL entries (participants without data)
group1_conn <- group1_conn[!sapply(group1_conn, is.null)]
group2_conn <- group2_conn[!sapply(group2_conn, is.null)]

# If you need to compute connectivity from raw time series:
# Load ROI atlas if available
# roi_atlas <- readNIfTI("path/to/your/atlas.nii.gz")
# group1_ts <- lapply(group1_folders, load_fmri_timeseries, roi_atlas = roi_atlas)
# group2_ts <- lapply(group2_folders, load_fmri_timeseries, roi_atlas = roi_atlas)

# For now, let's work with the connectivity matrices
# If matrices are correlation values, we'll work with them directly
# Check if we need Fisher Z-transformation
cat("\nLoaded", length(group1_conn), "connectivity matrices for Group 1\n")
cat("Loaded", length(group2_conn), "connectivity matrices for Group 2\n")

# Determine number of ROIs from connectivity matrix
if (length(group1_conn) > 0) {
  n_rois <- nrow(group1_conn[[1]])
  cat("Number of ROIs:", n_rois, "\n")
}

# ============================================================================
# 2. NETWORK ANALYSIS: FUNCTIONAL CONNECTIVITY
# ============================================================================

# Compute correlation matrices (functional connectivity) if needed
compute_fc_matrix <- function(timeseries) {
  cor(timeseries, use = "pairwise.complete.obs")
}

# Fisher Z-transformation for statistical testing
fisher_z <- function(r) {
  # Set diagonal to 0 to avoid Inf values
  diag(r) <- 0
  # Apply Fisher Z transformation
  z <- 0.5 * log((1 + r) / (1 - r))
  # Handle edge cases
  z[is.infinite(z)] <- NA
  z[is.nan(z)] <- NA
  return(z)
}

# If connectivity matrices are already computed, use them directly
# Otherwise compute from time series
if (exists("group1_ts") && !is.null(group1_ts[[1]])) {
  group1_fc <- lapply(group1_ts, compute_fc_matrix)
  group2_fc <- lapply(group2_ts, compute_fc_matrix)
} else {
  # Assume connectivity matrices are already loaded
  group1_fc <- group1_conn
  group2_fc <- group2_conn
}

# Apply Fisher Z-transformation
group1_fc_z <- lapply(group1_fc, fisher_z)
group2_fc_z <- lapply(group2_fc, fisher_z)

# Average connectivity matrices per group
group1_fc_mean <- Reduce("+", group1_fc_z) / length(group1_fc_z)
group2_fc_mean <- Reduce("+", group2_fc_z) / length(group2_fc_z)

# Statistical comparison: permutation test or t-test
compare_networks <- function(group1_fc_z, group2_fc_z, alpha = 0.05) {
  n1 <- length(group1_fc_z)
  n2 <- length(group2_fc_z)
  n_rois <- nrow(group1_fc_z[[1]])
  
  t_matrix <- matrix(NA, n_rois, n_rois)
  p_matrix <- matrix(NA, n_rois, n_rois)
  
  for (i in 1:(n_rois-1)) {
    for (j in (i+1):n_rois) {
      g1_values <- sapply(group1_fc_z, function(x) x[i, j])
      g2_values <- sapply(group2_fc_z, function(x) x[i, j])
      
      test_result <- t.test(g1_values, g2_values)
      t_matrix[i, j] <- test_result$statistic
      p_matrix[i, j] <- test_result$p.value
      
      # Make symmetric
      t_matrix[j, i] <- t_matrix[i, j]
      p_matrix[j, i] <- p_matrix[i, j]
    }
  }
  
  # FDR correction
  p_vector <- p_matrix[upper.tri(p_matrix)]
  p_fdr <- p.adjust(p_vector, method = "fdr")
  p_matrix_fdr <- matrix(NA, n_rois, n_rois)
  p_matrix_fdr[upper.tri(p_matrix_fdr)] <- p_fdr
  p_matrix_fdr[lower.tri(p_matrix_fdr)] <- t(p_matrix_fdr)[lower.tri(p_matrix_fdr)]
  
  return(list(t_stats = t_matrix, p_values = p_matrix_fdr, 
              diff_matrix = group2_fc_mean - group1_fc_mean))
}

network_comparison <- compare_networks(group1_fc_z, group2_fc_z)

# Visualize network differences
corrplot(network_comparison$diff_matrix, 
         method = "color", 
         title = "Group Differences in Functional Connectivity",
         tl.cex = 0.5)

# ============================================================================
# 3. LAG ANALYSIS: TEMPORAL DYNAMICS
# ============================================================================

# Cross-correlation with lags
compute_lag_matrix <- function(timeseries, max_lag = 5) {
  n_rois <- ncol(timeseries)
  lag_matrix <- matrix(NA, n_rois, n_rois)
  strength_matrix <- matrix(NA, n_rois, n_rois)
  
  for (i in 1:n_rois) {
    for (j in 1:n_rois) {
      if (i != j) {
        ccf_result <- ccf(timeseries[, i], timeseries[, j], 
                          lag.max = max_lag, plot = FALSE)
        max_idx <- which.max(abs(ccf_result$acf))
        lag_matrix[i, j] <- ccf_result$lag[max_idx]
        strength_matrix[i, j] <- ccf_result$acf[max_idx]
      }
    }
  }
  
  return(list(lags = lag_matrix, strengths = strength_matrix))
}

# Only compute lags if time series data is available
if (exists("group1_ts") && !is.null(group1_ts[[1]])) {
  cat("\nComputing lag analysis...\n")
  # Compute lags for all subjects
  group1_lags <- lapply(group1_ts, compute_lag_matrix)
  group2_lags <- lapply(group2_ts, compute_lag_matrix)
  
  # Compare lag distributions
  compare_lags <- function(group1_lags, group2_lags) {
    n_rois <- nrow(group1_lags[[1]]$lags)
    t_matrix <- matrix(NA, n_rois, n_rois)
    p_matrix <- matrix(NA, n_rois, n_rois)
    
    for (i in 1:n_rois) {
      for (j in 1:n_rois) {
        if (i != j) {
          g1_lags <- sapply(group1_lags, function(x) x$lags[i, j])
          g2_lags <- sapply(group2_lags, function(x) x$lags[i, j])
          
          test_result <- t.test(g1_lags, g2_lags)
          t_matrix[i, j] <- test_result$statistic
          p_matrix[i, j] <- test_result$p.value
        }
      }
    }
    
    return(list(t_stats = t_matrix, p_values = p_matrix))
  }
  
  lag_comparison <- compare_lags(group1_lags, group2_lags)
  cat("Lag analysis complete.\n")
} else {
  cat("\nSkipping lag analysis - time series data not available.\n")
  cat("Lag analysis requires raw time series data (.nii.gz files).\n")
}

# ============================================================================
# 4. VOLUMETRIC ANALYSIS
# ============================================================================

# Load structural MRI data (T1-weighted)
load_volumetric_data <- function(t1_files, roi_atlas) {
  volumes <- matrix(NA, nrow = length(t1_files), ncol = max(roi_atlas))
  
  for (subj in 1:length(t1_files)) {
    t1_img <- readNIfTI(t1_files[subj])
    
    for (roi in 1:max(roi_atlas)) {
      roi_mask <- roi_atlas == roi
      # Calculate volume (number of voxels * voxel dimensions)
      voxel_dims <- pixdim(t1_img)[2:4]
      voxel_volume <- prod(voxel_dims)
      volumes[subj, roi] <- sum(roi_mask) * voxel_volume
    }
  }
  
  return(volumes)
}

# Load T1 files for both groups
group1_t1_files <- list.files("path/to/group1/structural", 
                              pattern = "*.nii.gz", full.names = TRUE)
group2_t1_files <- list.files("path/to/group2/structural", 
                              pattern = "*.nii.gz", full.names = TRUE)

group1_volumes <- load_volumetric_data(group1_t1_files, roi_atlas)
group2_volumes <- load_volumetric_data(group2_t1_files, roi_atlas)

# Statistical comparison
volume_comparison <- data.frame(
  ROI = 1:ncol(group1_volumes),
  t_stat = NA,
  p_value = NA,
  mean_diff = NA,
  cohens_d = NA
)

for (roi in 1:ncol(group1_volumes)) {
  test_result <- t.test(group1_volumes[, roi], group2_volumes[, roi])
  volume_comparison$t_stat[roi] <- test_result$statistic
  volume_comparison$p_value[roi] <- test_result$p.value
  volume_comparison$mean_diff[roi] <- mean(group2_volumes[, roi]) - 
    mean(group1_volumes[, roi])
  
  # Cohen's d
  pooled_sd <- sqrt(((nrow(group1_volumes)-1) * var(group1_volumes[, roi]) + 
                       (nrow(group2_volumes)-1) * var(group2_volumes[, roi])) / 
                      (nrow(group1_volumes) + nrow(group2_volumes) - 2))
  volume_comparison$cohens_d[roi] <- volume_comparison$mean_diff[roi] / pooled_sd
}

# FDR correction
volume_comparison$p_fdr <- p.adjust(volume_comparison$p_value, method = "fdr")

# Visualize volumetric differences
significant_rois <- volume_comparison$ROI[volume_comparison$p_fdr < 0.05]
ggplot(volume_comparison[significant_rois, ], 
       aes(x = factor(ROI), y = cohens_d)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Volumetric Differences (Significant ROIs)",
       x = "ROI", y = "Cohen's d") +
  theme_minimal()

# ============================================================================
# 5. INTEGRATION WITH PROTEOMIC DATA
# ============================================================================

# Load proteomic data (example format: subjects x proteins)
proteomics <- read.csv("path/to/proteomic_data.csv", row.names = 1)
# Rows: subjects, Columns: protein markers

# Create a combined data frame
create_integrated_dataset <- function(group1_fc_z, group2_fc_z, 
                                      group1_lags, group2_lags,
                                      group1_volumes, group2_volumes,
                                      proteomics) {
  
  # Extract mean FC per subject (average across all connections)
  group1_fc_features <- sapply(group1_fc_z, function(x) {
    mean(x[upper.tri(x)], na.rm = TRUE)
  })
  group2_fc_features <- sapply(group2_fc_z, function(x) {
    mean(x[upper.tri(x)], na.rm = TRUE)
  })
  
  # Extract mean lag per subject
  group1_lag_features <- sapply(group1_lags, function(x) {
    mean(abs(x$lags), na.rm = TRUE)
  })
  group2_lag_features <- sapply(group2_lags, function(x) {
    mean(abs(x$lags), na.rm = TRUE)
  })
  
  # Combine all features
  integrated_data <- data.frame(
    Subject_ID = rownames(proteomics),
    Group = c(rep("Group1", length(group1_fc_features)),
              rep("Group2", length(group2_fc_features))),
    Mean_FC = c(group1_fc_features, group2_fc_features),
    Mean_Lag = c(group1_lag_features, group2_lag_features),
    Total_Volume = c(rowSums(group1_volumes), rowSums(group2_volumes)),
    proteomics
  )
  
  return(integrated_data)
}

integrated_data <- create_integrated_dataset(group1_fc_z, group2_fc_z,
                                             group1_lags, group2_lags,
                                             group1_volumes, group2_volumes,
                                             proteomics)

# Correlation analysis: brain features vs proteomics
brain_features <- c("Mean_FC", "Mean_Lag", "Total_Volume")
protein_markers <- colnames(proteomics)

correlation_results <- data.frame(
  Brain_Feature = rep(brain_features, each = length(protein_markers)),
  Protein = rep(protein_markers, length(brain_features)),
  Correlation = NA,
  P_value = NA
)

idx <- 1
for (bf in brain_features) {
  for (pm in protein_markers) {
    cor_test <- cor.test(integrated_data[[bf]], integrated_data[[pm]], 
                         method = "pearson")
    correlation_results$Correlation[idx] <- cor_test$estimate
    correlation_results$P_value[idx] <- cor_test$p.value
    idx <- idx + 1
  }
}

# FDR correction
correlation_results$P_fdr <- p.adjust(correlation_results$P_value, method = "fdr")

# Visualize significant correlations
significant_cors <- correlation_results[correlation_results$P_fdr < 0.05, ]
print(significant_cors)

# Heatmap of correlations
cor_matrix <- cor(integrated_data[, c(brain_features, protein_markers)], 
                  use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "upper",
         title = "Brain-Proteomics Correlations",
         tl.cex = 0.7, tl.col = "black")

# ============================================================================
# 6. MULTIVARIATE ANALYSIS: PARTIAL LEAST SQUARES (PLS)
# ============================================================================

library(pls)

# PLS regression: predict proteomic markers from brain features
X <- as.matrix(integrated_data[, brain_features])
Y <- as.matrix(integrated_data[, protein_markers])

pls_model <- plsr(Y ~ X, ncomp = 3, validation = "CV")

# Plot explained variance
plot(RMSEP(pls_model), legendpos = "topright")

# Extract loadings
pls_loadings <- loadings(pls_model)
print(pls_loadings)

# ============================================================================
# 7. MACHINE LEARNING: GROUP CLASSIFICATION
# ============================================================================

library(caret)
library(randomForest)

# Prepare data for classification
ml_data <- integrated_data[, c("Group", brain_features, protein_markers)]
ml_data$Group <- as.factor(ml_data$Group)

# Split data
set.seed(123)
train_idx <- createDataPartition(ml_data$Group, p = 0.7, list = FALSE)
train_data <- ml_data[train_idx, ]
test_data <- ml_data[-train_idx, ]

# Train Random Forest model
rf_model <- randomForest(Group ~ ., data = train_data, 
                         importance = TRUE, ntree = 500)

# Predictions
predictions <- predict(rf_model, test_data)
confusionMatrix(predictions, test_data$Group)

# Feature importance
importance_df <- data.frame(
  Feature = rownames(importance(rf_model)),
  Importance = importance(rf_model)[, "MeanDecreaseGini"]
)
importance_df <- importance_df[order(-importance_df$Importance), ]

ggplot(importance_df[1:10, ], aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "coral") +
  coord_flip() +
  labs(title = "Top 10 Feature Importance for Group Classification",
       x = "Feature", y = "Importance") +
  theme_minimal()

# ============================================================================
# 8. EXPORT RESULTS
# ============================================================================

# Save results
write.csv(network_comparison$diff_matrix, "network_differences.csv")
write.csv(volume_comparison, "volume_comparison.csv")
write.csv(correlation_results, "brain_protein_correlations.csv")
write.csv(integrated_data, "integrated_dataset.csv")

# Save plots
ggsave("volumetric_differences.pdf", width = 10, height = 6)

cat("Analysis complete! Results saved.\n")