############ Brain Biomarkers Analysis 



# fMRI Group Comparison Analysis with Proteomic Integration
# Required packages
library(neurobase)      # Neuroimaging analysis
library(oro.nifti)      # NIfTI file handling
library(ANTsR)          # Advanced normalization tools
library(igraph)         # Network analysis
library(brainGraph)     # Brain network analysis
library(ggplot2)        # Visualization
library(dplyr)          # Data manipulation
library(corrplot)       # Correlation plotting
library(lme4)           # Mixed models
library(psych)          # Statistical tools




library(dplyr)
library(stringr)
library(tidyr)






harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))





# ============================================================================
# 1. LOAD AND PREPROCESS fMRI DATA
# ============================================================================

# Function to load preprocessed fMRI time series
load_fmri_timeseries <- function(file_path, roi_atlas) {
  # Load 4D fMRI data
  fmri_img <- readNIfTI(file_path)
  
  # Extract time series from ROIs (regions of interest)
  # Assuming roi_atlas is a 3D image with labeled regions
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

# Load data for both groups
# Example structure: list of subjects per group
group1_files <- list.files("path/to/group1", pattern = "*.nii.gz", full.names = TRUE)
group2_files <- list.files("path/to/group2", pattern = "*.nii.gz", full.names = TRUE)

# Load ROI atlas (e.g., AAL, Schaefer, or custom)
roi_atlas <- readNIfTI("path/to/roi_atlas.nii.gz")
n_rois <- max(roi_atlas, na.rm = TRUE)

# Extract time series for all subjects
group1_ts <- lapply(group1_files, load_fmri_timeseries, roi_atlas = roi_atlas)
group2_ts <- lapply(group2_files, load_fmri_timeseries, roi_atlas = roi_atlas)

# ============================================================================
# 2. NETWORK ANALYSIS: FUNCTIONAL CONNECTIVITY
# ============================================================================

# Compute correlation matrices (functional connectivity)
compute_fc_matrix <- function(timeseries) {
  cor(timeseries, use = "pairwise.complete.obs")
}

# Fisher Z-transformation for statistical testing
fisher_z <- function(r) {
  0.5 * log((1 + r) / (1 - r))
}

# Compute FC for all subjects
group1_fc <- lapply(group1_ts, compute_fc_matrix)
group2_fc <- lapply(group2_ts, compute_fc_matrix)

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













