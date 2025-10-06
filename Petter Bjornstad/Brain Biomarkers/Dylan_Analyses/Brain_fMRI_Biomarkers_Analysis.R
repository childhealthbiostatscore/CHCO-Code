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

# Function to find participant folder by ID
find_participant_folder <- function(base_dir, participant_id) {
  all_folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  
  # Look for folders containing the participant ID
  matching_folders <- all_folders[grepl(participant_id, basename(all_folders), 
                                        ignore.case = TRUE)]
  
  if (length(matching_folders) > 0) {
    return(matching_folders[1])
  } else {
    warning(paste("No folder found for participant:", participant_id))
    return(NULL)
  }
}

# Function to load connectivity matrix from .mat file
load_connectivity_from_mat <- function(participant_folder) {
  if (is.null(participant_folder)) return(NULL)
  
  # Look for .mat files in the folder
  mat_files <- list.files(participant_folder, pattern = "\\.mat$", 
                          full.names = TRUE, recursive = TRUE)
  
  if (length(mat_files) > 0) {
    # Load the first .mat file found
    mat_data <- readMat(mat_files[1])
    # Extract connectivity matrix
    # You may need to adjust this based on your .mat structure
    if (length(mat_data) > 0) {
      conn_matrix <- mat_data[[1]]
      return(conn_matrix)
    }
  }
  
  warning(paste("No connectivity data found in", participant_folder))
  return(NULL)
}

# Load connectivity data for all participants
cat("\nLoading connectivity matrices...\n")

t2d_folders <- lapply(t2d_ids, find_participant_folder, base_dir = base_dir)
lc_folders <- lapply(lc_ids, find_participant_folder, base_dir = base_dir)

t2d_conn <- lapply(t2d_folders, load_connectivity_from_mat)
lc_conn <- lapply(lc_folders, load_connectivity_from_mat)

# Remove NULL entries (participants without data)
t2d_valid_idx <- !sapply(t2d_conn, is.null)
lc_valid_idx <- !sapply(lc_conn, is.null)

t2d_conn <- t2d_conn[t2d_valid_idx]
lc_conn <- lc_conn[lc_valid_idx]
t2d_ids_valid <- t2d_ids[t2d_valid_idx]
lc_ids_valid <- lc_ids[lc_valid_idx]

cat("Loaded", length(t2d_conn), "connectivity matrices for T2D group\n")
cat("Loaded", length(lc_conn), "connectivity matrices for LC group\n")

# Get number of ROIs
if (length(t2d_conn) > 0) {
  n_rois <- nrow(t2d_conn[[1]])
  cat("Number of ROIs:", n_rois, "\n")
}

# ==============================================================================
# 3. NETWORK ANALYSIS: FUNCTIONAL CONNECTIVITY
# ==============================================================================

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

# Apply Fisher Z-transformation
t2d_fc_z <- lapply(t2d_conn, fisher_z)
lc_fc_z <- lapply(lc_conn, fisher_z)

# Average connectivity matrices per group
t2d_fc_mean <- Reduce("+", t2d_fc_z) / length(t2d_fc_z)
lc_fc_mean <- Reduce("+", lc_fc_z) / length(lc_fc_z)

# Statistical comparison: T2D vs LC
compare_networks <- function(t2d_fc_z, lc_fc_z, alpha = 0.05) {
  n1 <- length(t2d_fc_z)
  n2 <- length(lc_fc_z)
  n_rois <- nrow(t2d_fc_z[[1]])
  
  t_matrix <- matrix(NA, n_rois, n_rois)
  p_matrix <- matrix(NA, n_rois, n_rois)
  
  for (i in 1:(n_rois-1)) {
    for (j in (i+1):n_rois) {
      t2d_values <- sapply(t2d_fc_z, function(x) x[i, j])
      lc_values <- sapply(lc_fc_z, function(x) x[i, j])
      
      # Remove NA values
      t2d_values <- t2d_values[!is.na(t2d_values)]
      lc_values <- lc_values[!is.na(lc_values)]
      
      if (length(t2d_values) > 1 && length(lc_values) > 1) {
        test_result <- t.test(t2d_values, lc_values)
        t_matrix[i, j] <- test_result$statistic
        p_matrix[i, j] <- test_result$p.value
        
        # Make symmetric
        t_matrix[j, i] <- t_matrix[i, j]
        p_matrix[j, i] <- p_matrix[i, j]
      }
    }
  }
  
  # FDR correction
  p_vector <- p_matrix[upper.tri(p_matrix)]
  p_vector <- p_vector[!is.na(p_vector)]
  p_fdr <- p.adjust(p_vector, method = "fdr")
  
  p_matrix_fdr <- matrix(NA, n_rois, n_rois)
  p_matrix_fdr[upper.tri(p_matrix_fdr)] <- p_fdr
  p_matrix_fdr[lower.tri(p_matrix_fdr)] <- t(p_matrix_fdr)[lower.tri(p_matrix_fdr)]
  
  # Calculate difference matrix (LC - T2D)
  diff_matrix <- lc_fc_mean - t2d_fc_mean
  
  return(list(
    t_stats = t_matrix, 
    p_values = p_matrix_fdr, 
    diff_matrix = diff_matrix,
    t2d_mean = t2d_fc_mean,
    lc_mean = lc_fc_mean
  ))
}

cat("\nPerforming network comparison...\n")
network_comparison <- compare_networks(t2d_fc_z, lc_fc_z)

# Count significant connections
sig_connections <- sum(network_comparison$p_values < 0.05, na.rm = TRUE) / 2
cat("Significant connections (FDR < 0.05):", sig_connections, "\n")

# Visualize network differences
pdf("network_differences_T2D_vs_LC.pdf", width = 10, height = 10)
corrplot(network_comparison$diff_matrix, 
         method = "color", 
         title = "Group Differences in FC (LC - T2D)",
         tl.cex = 0.5,
         mar = c(0,0,2,0))
dev.off()

# Visualize significant connections only
sig_diff <- network_comparison$diff_matrix
sig_diff[network_comparison$p_values >= 0.05] <- 0

pdf("significant_network_differences.pdf", width = 10, height = 10)
corrplot(sig_diff, 
         method = "color", 
         title = "Significant FC Differences (FDR < 0.05)",
         tl.cex = 0.5,
         mar = c(0,0,2,0))
dev.off()

# ==============================================================================
# 4. EXTRACT NETWORK METRICS FOR EACH PARTICIPANT
# ==============================================================================

# Function to calculate network metrics
calculate_network_metrics <- function(conn_matrix) {
  # Mean connectivity (excluding diagonal)
  diag(conn_matrix) <- NA
  mean_fc <- mean(conn_matrix, na.rm = TRUE)
  
  # Global efficiency (simplified)
  # Convert correlations to distances
  dist_matrix <- 1 - abs(conn_matrix)
  dist_matrix[dist_matrix <= 0] <- NA
  
  # Mean connectivity strength
  mean_strength <- mean(abs(conn_matrix), na.rm = TRUE)
  
  # Network density (proportion of strong connections)
  threshold <- 0.3  # Adjust as needed
  density <- sum(abs(conn_matrix) > threshold, na.rm = TRUE) / 
    sum(!is.na(conn_matrix))
  
  return(data.frame(
    mean_fc = mean_fc,
    mean_strength = mean_strength,
    network_density = density
  ))
}

# Calculate metrics for all participants
t2d_metrics <- do.call(rbind, lapply(t2d_fc_z, calculate_network_metrics))
lc_metrics <- do.call(rbind, lapply(lc_fc_z, calculate_network_metrics))

t2d_metrics$record_id <- t2d_ids_valid
lc_metrics$record_id <- lc_ids_valid

# Combine metrics
all_metrics <- rbind(
  t2d_metrics %>% mutate(group = "T2D"),
  lc_metrics %>% mutate(group = "LC")
)

# ==============================================================================
# 5. INTEGRATE WITH PROTEOMIC DATA
# ==============================================================================

# Merge fMRI metrics with proteomics
integrated_data <- small_dat %>%
  inner_join(all_metrics, by = c("record_id", "group"))

cat("\nIntegrated data: n =", nrow(integrated_data), "\n")
cat("T2D:", sum(integrated_data$group == "T2D"), 
    "| LC:", sum(integrated_data$group == "LC"), "\n")

# Define protein markers
protein_markers <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc", 
                     "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", 
                     "ptau_217_avg_conc")

# Define brain features
brain_features <- c("mean_fc", "mean_strength", "network_density")

# ==============================================================================
# 6. GROUP COMPARISONS: BRAIN METRICS
# ==============================================================================

cat("\n=== Brain Metrics Comparison ===\n")

brain_comparison <- data.frame(
  Metric = brain_features,
  T2D_mean = NA,
  LC_mean = NA,
  t_stat = NA,
  p_value = NA,
  cohens_d = NA
)

for (i in 1:length(brain_features)) {
  metric <- brain_features[i]
  
  t2d_vals <- integrated_data[[metric]][integrated_data$group == "T2D"]
  lc_vals <- integrated_data[[metric]][integrated_data$group == "LC"]
  
  brain_comparison$T2D_mean[i] <- mean(t2d_vals, na.rm = TRUE)
  brain_comparison$LC_mean[i] <- mean(lc_vals, na.rm = TRUE)
  
  test <- t.test(t2d_vals, lc_vals)
  brain_comparison$t_stat[i] <- test$statistic
  brain_comparison$p_value[i] <- test$p.value
  
  # Cohen's d
  pooled_sd <- sqrt(((length(t2d_vals)-1)*var(t2d_vals, na.rm=TRUE) + 
                       (length(lc_vals)-1)*var(lc_vals, na.rm=TRUE)) / 
                      (length(t2d_vals) + length(lc_vals) - 2))
  brain_comparison$cohens_d[i] <- (brain_comparison$LC_mean[i] - 
                                     brain_comparison$T2D_mean[i]) / pooled_sd
}

print(brain_comparison)

# Visualize brain metrics
brain_data_long <- integrated_data %>%
  pivot_longer(cols = all_of(brain_features), 
               names_to = "Metric", 
               values_to = "Value")

ggplot(brain_data_long, aes(x = group, y = Value, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("T2D" = "#E74C3C", "LC" = "#3498DB")) +
  labs(title = "Brain Network Metrics: T2D vs LC",
       x = "Group", y = "Value") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("brain_metrics_comparison.pdf", width = 12, height = 5)

# ==============================================================================
# 7. CORRELATION ANALYSIS: BRAIN FEATURES vs PROTEOMICS
# ==============================================================================

cat("\n=== Brain-Proteomics Correlation Analysis ===\n")

# Create correlation matrix
cor_data <- integrated_data %>%
  select(all_of(c(brain_features, protein_markers))) %>%
  na.omit()

cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

# Focus on brain-protein correlations
brain_protein_cors <- cor_matrix[brain_features, protein_markers]

# Visualize correlations
pdf("brain_protein_correlations.pdf", width = 10, height = 6)
corrplot(brain_protein_cors, method = "color", 
         title = "Brain Metrics vs Proteomic Markers",
         tl.col = "black", tl.srt = 45,
         mar = c(0,0,2,0),
         addCoef.col = "black", number.cex = 0.7)
dev.off()

# Statistical testing for each correlation
correlation_results <- expand.grid(
  Brain_Feature = brain_features,
  Protein = protein_markers,
  stringsAsFactors = FALSE
)

correlation_results$Correlation <- NA
correlation_results$P_value <- NA
correlation_results$N <- NA

for (i in 1:nrow(correlation_results)) {
  bf <- correlation_results$Brain_Feature[i]
  pm <- correlation_results$Protein[i]
  
  complete_cases <- complete.cases(integrated_data[[bf]], integrated_data[[pm]])
  
  if (sum(complete_cases) > 3) {
    cor_test <- cor.test(integrated_data[[bf]][complete_cases], 
                         integrated_data[[pm]][complete_cases], 
                         method = "pearson")
    correlation_results$Correlation[i] <- cor_test$estimate
    correlation_results$P_value[i] <- cor_test$p.value
    correlation_results$N[i] <- sum(complete_cases)
  }
}

# FDR correction
correlation_results$P_fdr <- p.adjust(correlation_results$P_value, method = "fdr")

# Show significant correlations
significant_cors <- correlation_results %>%
  filter(P_fdr < 0.05) %>%
  arrange(P_fdr)

cat("\nSignificant Brain-Protein Correlations (FDR < 0.05):\n")
print(significant_cors)

# ==============================================================================
# 8. GROUP-SPECIFIC CORRELATIONS
# ==============================================================================

cat("\n=== Group-Specific Correlations ===\n")

# Function to compute correlations within a group
group_correlations <- function(data, group_name) {
  group_data <- data %>% filter(group == group_name)
  
  results <- expand.grid(
    Brain_Feature = brain_features,
    Protein = protein_markers,
    stringsAsFactors = FALSE
  )
  
  results$Correlation <- NA
  results$P_value <- NA
  results$Group <- group_name
  
  for (i in 1:nrow(results)) {
    bf <- results$Brain_Feature[i]
    pm <- results$Protein[i]
    
    complete_cases <- complete.cases(group_data[[bf]], group_data[[pm]])
    
    if (sum(complete_cases) > 3) {
      cor_test <- cor.test(group_data[[bf]][complete_cases], 
                           group_data[[pm]][complete_cases], 
                           method = "pearson")
      results$Correlation[i] <- cor_test$estimate
      results$P_value[i] <- cor_test$p.value
    }
  }
  
  return(results)
}

t2d_cors <- group_correlations(integrated_data, "T2D")
lc_cors <- group_correlations(integrated_data, "LC")

# Combine and compare
all_group_cors <- rbind(t2d_cors, lc_cors)
all_group_cors$P_fdr <- ave(all_group_cors$P_value, 
                            all_group_cors$Group, 
                            FUN = function(x) p.adjust(x, method = "fdr"))

# ==============================================================================
# 9. MACHINE LEARNING: GROUP CLASSIFICATION
# ==============================================================================

cat("\n=== Machine Learning Classification ===\n")

# Prepare data for classification
ml_data <- integrated_data %>%
  select(group, all_of(brain_features), all_of(protein_markers)) %>%
  na.omit()

ml_data$group <- as.factor(ml_data$group)

cat("ML dataset: n =", nrow(ml_data), "\n")

if (nrow(ml_data) >= 10) {  # Need sufficient samples
  # Split data
  set.seed(123)
  train_idx <- createDataPartition(ml_data$group, p = 0.7, list = FALSE)
  train_data <- ml_data[train_idx, ]
  test_data <- ml_data[-train_idx, ]
  
  # Train Random Forest
  rf_model <- randomForest(group ~ ., data = train_data, 
                           importance = TRUE, ntree = 500)
  
  # Predictions
  predictions <- predict(rf_model, test_data)
  conf_matrix <- confusionMatrix(predictions, test_data$group)
  
  cat("\nClassification Accuracy:", 
      round(conf_matrix$overall["Accuracy"], 3), "\n")
  print(conf_matrix$table)
  
  # Feature importance
  importance_df <- data.frame(
    Feature = rownames(importance(rf_model)),
    Importance = importance(rf_model)[, "MeanDecreaseGini"]
  ) %>%
    arrange(desc(Importance))
  
  # Plot feature importance
  ggplot(importance_df[1:min(10, nrow(importance_df)), ], 
         aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "coral") +
    coord_flip() +
    labs(title = "Feature Importance for T2D vs LC Classification",
         x = "Feature", y = "Importance (Mean Decrease Gini)") +
    theme_minimal()
  
  ggsave("feature_importance.pdf", width = 8, height = 6)
}

# ==============================================================================
# 10. SAVE RESULTS
# ==============================================================================

cat("\n=== Saving Results ===\n")

# Save all results
write.csv(integrated_data, "integrated_brain_protein_data.csv", row.names = FALSE)
write.csv(brain_comparison, "brain_metrics_comparison.csv", row.names = FALSE)
write.csv(correlation_results, "brain_protein_correlations.csv", row.names = FALSE)
write.csv(all_group_cors, "group_specific_correlations.csv", row.names = FALSE)

# Save network comparison matrices
write.csv(network_comparison$diff_matrix, "fc_difference_matrix.csv")
write.csv(network_comparison$p_values, "fc_pvalues_matrix.csv")

cat("\nAnalysis complete! Results saved.\n")
cat("\nFiles created:\n")
cat("- integrated_brain_protein_data.csv\n")
cat("- brain_metrics_comparison.csv\n")
cat("- brain_protein_correlations.csv\n")
cat("- group_specific_correlations.csv\n")
cat("- fc_difference_matrix.csv\n")
cat("- fc_pvalues_matrix.csv\n")
cat("- network_differences_T2D_vs_LC.pdf\n")
cat("- significant_network_differences.pdf\n")
cat("- brain_metrics_comparison.pdf\n")
cat("- brain_protein_correlations.pdf\n")
cat("- feature_importance.pdf\n")