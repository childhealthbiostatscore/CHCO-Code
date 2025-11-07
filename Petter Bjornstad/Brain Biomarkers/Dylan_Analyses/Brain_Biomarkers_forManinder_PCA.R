


library(scran)
library(future)
library(future.apply)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(REDCapR)
library(data.table)
library(emmeans)
library(NMF)
library(pheatmap)
library(UpSetR)
library(enrichR)
library(WriteXLS)
library(SAVER)
library(readxl)
library(limma)
library(edgeR)
library(BiocGenerics)
library(GSEABase)
library(slingshot)
library(SingleCellExperiment)
library(MAST)
library(muscat)
library(scater)
library(Seurat)
library(jsonlite)
library(dplyr)
library(glmmTMB)
library(reshape2)
library(broom.mixed)
library(nebula)
#library(table1)
library(clusterProfiler)
library('org.Hs.eg.db')
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(readxl)
library(stringr)
library(httr)

qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

names_dat <- names(dat)

dat <- dat %>% filter(study == 'CROCODILE') %>% 
  filter(!is.na(ab40_avg_conc)) %>% filter(visit == 'baseline')



data_dictionary <- readxl::read_xlsx('C:/Users/netio/Downloads/data_dictionary_master.xlsx')
#data_dictionary <- data_dictionary %>% filter(form_name %in% c('clamp', 'UACR', 'fmri', 'eGFR', 'FSOC', 'brain_mri', 'pet_scan', 'dextran_data'))


traits_of_interest <- c('acr_u', 
                        'eGFR_bedside_Schwartz', 'eGFR_CKD_epi', 'eGFR_fas_cr', 'eGFR_fas_cr_cysc','eGFR_Zap','eGFR_Schwartz', 
                        'hba1c')

dat_analysis <- dat %>% 
  dplyr::select(record_id, group, hba1c, #age, sex, bmi, hba1c, study, 
                all_of(qx_var), 
                any_of(data_dictionary$variable_name))


dat_analysis <- dat %>% 
  dplyr::select(record_id, group, all_of(qx_var), acr_u, adipose_ir, cholesterol, hba1c, eGFR_CKD_epi, eGFR_fas_cr, fbg, ldl, hdl, left_kidney_volume_ml, right_kidney_volume_ml, 
                triglycerides, urine_glucose_bl, avg_k_fsoc, avg_c_fsoc, avg_m_fsoc, homa_ir, search_eis)


dat_clean <- dat_analysis[complete.cases(dat_analysis), ]

dat_filtered <- dat %>% filter(record_id %in% dat_clean$record_id)

# Identify columns with NO missing values
complete_cols <- colnames(dat_filtered)[colSums(is.na(dat_filtered)) == 0]

# View the list
print(complete_cols)

# Extract only those columns
data_complete <- dat_filtered[, complete_cols]

# Identify columns with more than 1 unique value (excluding NAs)
varying_cols <- sapply(data_complete, function(x) length(unique(na.omit(x))) > 1)


# Get list of metabolomics variables to remove
metabolomics_vars <- data_dictionary %>%
  filter(form_name %in% c("metabolomics", "metabolomics_blood_raw")) %>%
  pull(variable_name)



# Keep only varying columns
data_varying <- data_complete[, varying_cols]

data_varying <- data_varying %>% 
  dplyr::select(
    -c('croc_id', 'date', 'mrn'),
    -matches("^glucose_\\d+"),
    -matches('^insulin_\\d+'),
    -matches('^ffa_\\d+'),
    -matches("^c0"),             
    -matches("^ac\\d+"), 
    -matches('^c1'), 
    -matches('^hmdb\\d+'), 
    -any_of(metabolomics_vars)
  ) %>% 
  dplyr::select(-hx_met_positive___1, -insulin_injections_timepoint, -insulin_med_timepoint, -insulin_pump_timepoint, 
                -u24_labs, -insulin_minus_10, -insulin_minus_20, -glucose_minus_10, -glucose_minus_20, -ffa_minus_10, -ffa_minus_20, -group_risk, 
                -record_id, -sex, -race, -ethnicity)

data_varying <- data_varying %>% 
  dplyr::select(-length_left, -length_right, 
                -depth_left, -depth_right, 
                -width_left, -width_right)


data_dictionary_small <- data_dictionary %>% filter(variable_name %in% names(data_varying))







##### Differences in Differences 


#### PCA Analysis - Lean Control vs T1D Comparison
set.seed(123)

setwd('C:/Users/netio/Documents/UofW/Projects/Maninder_Data/PCA_Analysis/')
# Load required libraries
library(tidyverse)
library(factoextra)
library(corrplot)
library(pheatmap)
library(RColorBrewer)

# Define your Quanterix biomarkers
qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")

names(data_varying)[which(names(data_varying) == 'weight')] <- 'body_weight'

# Define analysis groups - MODIFIED FOR LEAN CONTROL AND T1D
analysis_groups <- list(
  Lean_Control = data_varying %>% filter(group == "Lean Control"),  # Lean controls only
  T1D = data_varying %>% filter(group == "Type 1 Diabetes")  # T1D only
)

cat("Sample sizes:\n")
cat("Lean Control:", nrow(analysis_groups$Lean_Control), "\n")
cat("Type 1 Diabetes:", nrow(analysis_groups$T1D), "\n\n")

# Initialize storage for all results
all_results <- list()
all_performance <- data.frame()
all_variable_importance <- list()

# ============================================================
# LOOP THROUGH EACH GROUP
# ============================================================

for(group_name in names(analysis_groups)) {
  
  cat("\n\n##########################################################")
  cat("\n### ANALYSIS GROUP:", group_name)
  cat("\n##########################################################\n\n")
  
  # Get data for this group
  current_data <- analysis_groups[[group_name]]
  
  # RENAME 'weight' TO AVOID CONFLICT WITH prcomp()
  if ("weight" %in% names(current_data)) {
    current_data <- current_data %>%
      rename(body_weight = weight)
    cat("Renamed 'weight' to 'body_weight' to avoid prcomp conflict\n\n")
  }
  
  # Step 1: Identify predictor variables (everything except biomarkers and group)
  potential_predictors <- setdiff(names(current_data), c(qx_var, "group"))
  
  # Filter to only numeric columns
  predictor_cols <- potential_predictors[sapply(current_data[, potential_predictors], is.numeric)]
  
  cat("Total columns (excluding biomarkers and group):", length(potential_predictors), "\n")
  cat("Numeric predictor variables:", length(predictor_cols), "\n")
  
  # Show which columns were excluded (non-numeric)
  non_numeric <- setdiff(potential_predictors, predictor_cols)
  if(length(non_numeric) > 0) {
    cat("Non-numeric columns excluded:", paste(head(non_numeric, 10), collapse = ", "), 
        ifelse(length(non_numeric) > 10, "...", ""), "\n")
  }
  
  cat("Sample size:", nrow(current_data), "\n")
  
  # Extract predictors and outcomes
  X <- current_data[, predictor_cols]
  Y <- current_data[, qx_var]
  
  # Handle missing data
  cat("\nMissing values in predictors:", sum(is.na(X)), "\n")
  cat("Missing values in outcomes:", sum(is.na(Y)), "\n")
  
  complete_rows <- complete.cases(X, Y)
  X_complete <- X[complete_rows, ]
  Y_complete <- Y[complete_rows, ]
  
  cat("Sample size after removing missing data:", nrow(X_complete), "\n")
  
  # Check if we have enough data
  if(nrow(X_complete) < 30) {
    cat("\nWARNING: Small sample size for", group_name, "(N =", nrow(X_complete), 
        ") - results may be unreliable\n")
  }
  
  # Remove any columns with zero variance
  col_vars <- sapply(names(X_complete), function(col_name) {
    var(X_complete[[col_name]], na.rm = TRUE)
  })
  names(col_vars) <- names(X_complete)
  
  zero_var_cols <- names(col_vars[col_vars == 0 | is.na(col_vars)])
  
  if(length(zero_var_cols) > 0) {
    cat("\nRemoving", length(zero_var_cols), "zero-variance columns:", 
        paste(head(zero_var_cols, 10), collapse = ", "), 
        ifelse(length(zero_var_cols) > 10, "...", ""), "\n")
    X_complete <- X_complete[, !names(X_complete) %in% zero_var_cols]
  }
  
  cat("Final number of predictor variables:", ncol(X_complete), "\n")
  
  # Step 2: Perform PCA
  pca_result <- prcomp(X_complete, center = TRUE, scale. = TRUE)
  
  # PCA Summary
  cat("\n=== PCA Summary for", group_name, "===\n")
  pca_summary <- summary(pca_result)
  print(pca_summary$importance[, 1:min(10, ncol(pca_summary$importance))])
  
  # Scree plot - save as PNG
  p_scree <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50), 
                      main = paste("Scree Plot:", group_name))
  ggsave(paste0("scree_plot_", group_name, ".png"), p_scree, 
         width = 10, height = 6, dpi = 300, bg = "white")
  print(p_scree)
  
  # Determine number of components
  var_explained <- pca_summary$importance[3, ]
  n_components_80 <- which(var_explained >= 0.80)[1]
  n_components_90 <- which(var_explained >= 0.90)[1]
  
  cat("\nPCs for 80% variance:", n_components_80)
  cat("\nPCs for 90% variance:", n_components_90, "\n")
  
  n_components <- n_components_80  # Choose threshold
  
  # Extract PC scores
  PC_scores <- as.data.frame(pca_result$x[, 1:n_components])
  colnames(PC_scores) <- paste0("PC", 1:n_components)
  
  # Get loadings
  loadings <- pca_result$rotation[, 1:n_components]
  
  # Show top variables for first few PCs
  cat("\n=== TOP VARIABLES FOR EACH PC (", group_name, ") ===\n")
  for(i in 1:min(5, n_components)) {
    cat("\n--- PC", i, "(", round(pca_summary$importance[2, i] * 100, 2), "% variance) ---\n")
    
    sorted_loadings <- sort(abs(loadings[, i]), decreasing = TRUE)
    top_vars <- names(sorted_loadings)[1:min(10, length(sorted_loadings))]
    
    loading_df <- data.frame(
      Variable = top_vars,
      Loading = round(loadings[top_vars, i], 3),
      Abs_Loading = round(abs(loadings[top_vars, i]), 3)
    )
    print(loading_df)
  }
  
  # Save loadings
  write.csv(loadings, paste0("PCA_loadings_", group_name, ".csv"), row.names = TRUE)
  
  # Scale biomarkers for comparability
  cat("\n=== SCALING BIOMARKERS ===\n")
  Y_complete_scaled <- as.data.frame(scale(Y_complete))
  colnames(Y_complete_scaled) <- qx_var
  cat("Biomarkers scaled to mean=0, sd=1 for comparability\n")
  
  # Step 3: Linear Regression for each biomarker
  cat("\n\n=== LINEAR REGRESSION (", group_name, ") ===\n")
  
  group_results <- list()
  group_performance <- data.frame()
  group_var_importance <- list()
  
  for(biomarker in qx_var) {
    cat("\n----------------------------------------")
    cat("\nBiomarker:", biomarker, "(", group_name, ")")
    cat("\n----------------------------------------\n")
    
    # Prepare data - USE SCALED VERSION
    y <- Y_complete_scaled[[biomarker]]
    
    # Create data frame for regression
    reg_data <- cbind(PC_scores, outcome = y)
    
    # Fit linear model
    lm_formula <- as.formula(paste("outcome ~", paste(colnames(PC_scores), collapse = " + ")))
    lm_model <- lm(lm_formula, data = reg_data)
    
    # Get model summary
    model_summary <- summary(lm_model)
    
    # Extract coefficients (excluding intercept)
    coefs <- coef(lm_model)[-1]  # Remove intercept
    
    # Performance metrics
    r2 <- model_summary$r.squared
    adj_r2 <- model_summary$adj.r.squared
    f_stat <- model_summary$fstatistic[1]
    p_value <- pf(f_stat, 
                  model_summary$fstatistic[2], 
                  model_summary$fstatistic[3], 
                  lower.tail = FALSE)
    
    predictions <- predict(lm_model)
    mse <- mean((y - predictions)^2)
    rmse <- sqrt(mse)
    
    # Count significant PCs (p < 0.05)
    pc_pvalues <- model_summary$coefficients[-1, 4]  # Remove intercept row
    n_significant <- sum(pc_pvalues < 0.05)
    
    cat("\nModel Performance:")
    cat("\n  R-squared:", round(r2, 3))
    cat("\n  Adjusted R-squared:", round(adj_r2, 3))
    cat("\n  RMSE:", round(rmse, 3))
    cat("\n  F-statistic:", round(f_stat, 2), "(p =", format.pval(p_value, digits = 3), ")")
    cat("\n  Significant PCs (p < 0.05):", n_significant, "\n")
    
    # Show all PC coefficients with p-values
    cat("\nPC Coefficients:\n")
    coef_table <- model_summary$coefficients[-1, ]  # Remove intercept
    coef_df <- data.frame(
      PC = rownames(coef_table),
      Coefficient = round(coef_table[, 1], 4),
      Std_Error = round(coef_table[, 2], 4),
      t_value = round(coef_table[, 3], 3),
      p_value = format.pval(coef_table[, 4], digits = 3),
      Significant = ifelse(coef_table[, 4] < 0.05, "***", 
                           ifelse(coef_table[, 4] < 0.10, "*", ""))
    )
    print(coef_df)
    
    # Store results
    group_results[[biomarker]] <- list(
      model = lm_model,
      summary = model_summary,
      coefficients = coefs,
      R2 = r2,
      Adj_R2 = adj_r2,
      RMSE = rmse,
      F_stat = f_stat,
      p_value = p_value,
      predictions = predictions
    )
    
    # Performance summary
    perf_row <- data.frame(
      Group = group_name,
      Biomarker = biomarker,
      N = nrow(X_complete),
      N_PCs = n_components,
      R2 = r2,
      Adj_R2 = adj_r2,
      RMSE = rmse,
      F_stat = f_stat,
      Model_p_value = p_value,
      N_Significant_PCs = n_significant
    )
    group_performance <- rbind(group_performance, perf_row)
    
    # Calculate variable importance (loading × coefficient)
    var_importance <- loadings %*% coefs
    
    group_var_importance[[biomarker]] <- var_importance
    
    # Show top variables
    cat("\nTop 20 most important variables:\n")
    sorted_importance <- sort(abs(var_importance), decreasing = TRUE)
    top_n <- min(20, length(sorted_importance))
    top_indices <- order(abs(var_importance), decreasing = TRUE)[1:top_n]
    top_var_names <- rownames(var_importance)[top_indices]
    
    importance_df <- data.frame(
      Variable = top_var_names,
      Importance = round(var_importance[top_var_names, ], 4),
      Abs_Importance = round(abs(var_importance[top_var_names, ]), 4)
    )
    rownames(importance_df) <- NULL
    print(importance_df)
  }
  
  # Store group-level results
  all_results[[group_name]] <- list(
    pca_result = pca_result,
    loadings = loadings,
    PC_scores = PC_scores,
    model_results = group_results,
    variable_importance = group_var_importance,
    predictor_cols_used = names(X_complete)
  )
  
  all_performance <- rbind(all_performance, group_performance)
  all_variable_importance[[group_name]] <- group_var_importance
  
  # Create importance matrix for this group
  importance_matrix <- do.call(cbind, group_var_importance)
  colnames(importance_matrix) <- qx_var
  
  write.csv(importance_matrix, 
            paste0("variable_importance_", group_name, ".csv"), 
            row.names = TRUE)
}

# ============================================================
# SIDE-BY-SIDE COMPARISON HEATMAPS (LEAN CONTROL VS T1D)
# ============================================================

cat("\n\n##########################################################")
cat("\n### CREATING SIDE-BY-SIDE COMPARISON HEATMAPS")
cat("\n##########################################################\n\n")

# Get importance matrices for both groups
lean_importance <- do.call(cbind, all_variable_importance$Lean_Control)
t1d_importance <- do.call(cbind, all_variable_importance$T1D)

colnames(lean_importance) <- qx_var
colnames(t1d_importance) <- qx_var

# Find common variables between both groups
common_vars <- intersect(rownames(lean_importance), rownames(t1d_importance))
cat("Common variables between groups:", length(common_vars), "\n")

if(length(common_vars) == 0) {
  cat("ERROR: No common variables between groups!\n")
} else {
  
  # Create combined matrix for comparison
  lean_subset <- lean_importance[common_vars, ]
  t1d_subset <- t1d_importance[common_vars, ]
  
  # Calculate overall importance for ranking
  overall_importance_lean <- rowMeans(abs(lean_subset))
  overall_importance_t1d <- rowMeans(abs(t1d_subset))
  overall_importance_combined <- (overall_importance_lean + overall_importance_t1d) / 2
  
  # ============================================================
  # 1. TOP 20 VARIABLES - SIDE BY SIDE COMPARISON
  # ============================================================
  
  top_20_vars <- names(sort(overall_importance_combined, decreasing = TRUE)[1:min(20, length(overall_importance_combined))])
  
  # Create combined matrix with group labels
  combined_top20 <- cbind(lean_subset[top_20_vars, ], t1d_subset[top_20_vars, ])
  colnames(combined_top20) <- c(paste0("Lean_", qx_var), paste0("T1D_", qx_var))
  
  # Create annotation for columns
  annotation_col <- data.frame(
    Group = c(rep("Lean Control", length(qx_var)), rep("T1D", length(qx_var)))
  )
  rownames(annotation_col) <- colnames(combined_top20)
  
  # Color scheme for groups
  ann_colors <- list(
    Group = c("Lean Control" = "#4DAF4A", "T1D" = "#E41A1C")
  )
  
  pheatmap(combined_top20,
           scale = "none",
           cluster_cols = FALSE,
           main = "Top 20 Variables: Lean Control vs T1D",
           fontsize_row = 9,
           fontsize_col = 8,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           display_numbers = matrix(sprintf("%.2f", combined_top20), nrow = nrow(combined_top20)),
           number_color = "black",
           fontsize_number = 6,
           filename = "comparison_top20_lean_vs_t1d.png",
           width = 14,
           height = 8)
  
  pheatmap(combined_top20,
           scale = "none",
           cluster_cols = FALSE,
           main = "Top 20 Variables: Lean Control vs T1D",
           fontsize_row = 9,
           fontsize_col = 8,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           display_numbers = matrix(sprintf("%.2f", combined_top20), nrow = nrow(combined_top20)),
           number_color = "black",
           fontsize_number = 6,
           filename = "comparison_top20_lean_vs_t1d.pdf",
           width = 14,
           height = 8)
  
  cat("Created: comparison_top20_lean_vs_t1d.png/.pdf\n")
  
  # ============================================================
  # 2. TOP 50 VARIABLES - SIDE BY SIDE
  # ============================================================
  
  if(length(common_vars) >= 50) {
    top_50_vars <- names(sort(overall_importance_combined, decreasing = TRUE)[1:50])
    combined_top50 <- cbind(lean_subset[top_50_vars, ], t1d_subset[top_50_vars, ])
    colnames(combined_top50) <- c(paste0("Lean_", qx_var), paste0("T1D_", qx_var))
    
    pheatmap(combined_top50,
             scale = "none",
             cluster_cols = FALSE,
             main = "Top 50 Variables: Lean Control vs T1D",
             fontsize_row = 6,
             fontsize_col = 7,
             angle_col = 45,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             filename = "comparison_top50_lean_vs_t1d.png",
             width = 14,
             height = 14)
    
    pheatmap(combined_top50,
             scale = "none",
             cluster_cols = FALSE,
             main = "Top 50 Variables: Lean Control vs T1D",
             fontsize_row = 6,
             fontsize_col = 7,
             angle_col = 45,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             filename = "comparison_top50_lean_vs_t1d.pdf",
             width = 14,
             height = 14)
    
    cat("Created: comparison_top50_lean_vs_t1d.png/.pdf\n")
  }
  
  # ============================================================
  # 3. DIFFERENCE HEATMAP (T1D - Lean)
  # ============================================================
  
  difference_matrix <- t1d_subset - lean_subset
  top_20_diff_vars <- names(sort(rowMeans(abs(difference_matrix)), decreasing = TRUE)[1:min(20, nrow(difference_matrix))])
  
  pheatmap(difference_matrix[top_20_diff_vars, ],
           scale = "none",
           main = "Top 20 Largest Differences (T1D - Lean Control)",
           fontsize_row = 9,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           display_numbers = matrix(sprintf("%.2f", difference_matrix[top_20_diff_vars, ]), 
                                    nrow = length(top_20_diff_vars)),
           number_color = "black",
           fontsize_number = 7,
           filename = "comparison_difference_top20.png",
           width = 10,
           height = 8)
  
  pheatmap(difference_matrix[top_20_diff_vars, ],
           scale = "none",
           main = "Top 20 Largest Differences (T1D - Lean Control)",
           fontsize_row = 9,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           display_numbers = matrix(sprintf("%.2f", difference_matrix[top_20_diff_vars, ]), 
                                    nrow = length(top_20_diff_vars)),
           number_color = "black",
           fontsize_number = 7,
           filename = "comparison_difference_top20.pdf",
           width = 10,
           height = 8)
  
  cat("Created: comparison_difference_top20.png/.pdf\n")
  
  # Save difference matrix
  write.csv(difference_matrix, "importance_difference_t1d_minus_lean.csv", row.names = TRUE)
  
  # ============================================================
  # 4. INDIVIDUAL BIOMARKER COMPARISONS (3-PANEL HEATMAPS)
  # ============================================================
  
  cat("\nCreating individual biomarker comparison heatmaps (LC | T1D | Difference)...\n")
  dir.create("biomarker_comparisons", showWarnings = FALSE)
  
  for(biomarker in qx_var) {
    cat("  Processing", biomarker, "...\n")
    
    # Get top 30 variables for this biomarker across both groups
    lean_bio <- lean_subset[, biomarker]
    t1d_bio <- t1d_subset[, biomarker]
    combined_abs <- (abs(lean_bio) + abs(t1d_bio)) / 2
    
    top_30_bio <- names(sort(combined_abs, decreasing = TRUE)[1:min(30, length(combined_abs))])
    
    # Create three-column matrix: Lean Control | T1D | Difference
    comparison_matrix <- cbind(
      Lean_Control = lean_bio[top_30_bio],
      T1D = t1d_bio[top_30_bio],
      Difference = t1d_bio[top_30_bio] - lean_bio[top_30_bio]
    )
    
    # Create annotation for columns
    annotation_col_bio <- data.frame(
      Type = c("Lean Control", "T1D", "T1D - Lean")
    )
    rownames(annotation_col_bio) <- colnames(comparison_matrix)
    
    # Color scheme
    ann_colors_bio <- list(
      Type = c("Lean Control" = "#4DAF4A", "T1D" = "#E41A1C", "T1D - Lean" = "#984EA3")
    )
    
    # Create heatmap - PNG
    pheatmap(comparison_matrix,
             scale = "none",
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             main = paste(biomarker, ": Lean Control | T1D | Difference"),
             fontsize_row = 8,
             fontsize_col = 10,
             angle_col = 45,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             annotation_col = annotation_col_bio,
             annotation_colors = ann_colors_bio,
             display_numbers = matrix(sprintf("%.2f", comparison_matrix), nrow = nrow(comparison_matrix)),
             number_color = "black",
             fontsize_number = 7,
             filename = paste0("biomarker_comparisons/heatmap_", biomarker, ".png"),
             width = 10,
             height = 12)
    
    # Create heatmap - PDF
    pheatmap(comparison_matrix,
             scale = "none",
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             main = paste(biomarker, ": Lean Control | T1D | Difference"),
             fontsize_row = 8,
             fontsize_col = 10,
             angle_col = 45,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             annotation_col = annotation_col_bio,
             annotation_colors = ann_colors_bio,
             display_numbers = matrix(sprintf("%.2f", comparison_matrix), nrow = nrow(comparison_matrix)),
             number_color = "black",
             fontsize_number = 7,
             filename = paste0("biomarker_comparisons/heatmap_", biomarker, ".pdf"),
             width = 10,
             height = 12)
    
    # Also create a grouped bar plot for alternative visualization
    comparison_df <- data.frame(
      Variable = factor(top_30_bio, levels = rev(top_30_bio)),
      Lean_Control = lean_bio[top_30_bio],
      T1D = t1d_bio[top_30_bio]
    )
    
    comparison_long <- comparison_df %>%
      pivot_longer(cols = c(Lean_Control, T1D), 
                   names_to = "Group", 
                   values_to = "Importance")
    
    p <- ggplot(comparison_long, aes(x = Variable, y = Importance, fill = Group)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Lean_Control" = "#4DAF4A", "T1D" = "#E41A1C"),
                        labels = c("Lean Control", "T1D")) +
      coord_flip() +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8),
            legend.position = "top") +
      labs(title = paste("Top 30 Variables for", biomarker, ": Lean Control vs T1D"),
           x = "Variable",
           y = "Variable Importance",
           fill = "Group") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
    
    ggsave(paste0("biomarker_comparisons/barplot_", biomarker, ".png"),
           p, width = 12, height = 10, dpi = 300, bg = "white")
    
    ggsave(paste0("biomarker_comparisons/barplot_", biomarker, ".pdf"),
           p, width = 12, height = 10)
    
    # Create side-by-side horizontal bar plot (like the overall plot)
    comparison_horizontal <- data.frame(
      Variable = factor(top_30_bio, levels = top_30_bio),  # Keep order from top to bottom
      Lean_Control = lean_bio[top_30_bio],
      T1D = t1d_bio[top_30_bio]
    )
    
    comparison_horizontal_long <- comparison_horizontal %>%
      pivot_longer(cols = c(Lean_Control, T1D), 
                   names_to = "Group", 
                   values_to = "Importance") %>%
      mutate(Group = recode(Group, 
                            "Lean_Control" = "Lean Control",
                            "T1D" = "T1D"))
    
    p_horizontal <- ggplot(comparison_horizontal_long, 
                           aes(x = Importance, y = Variable, fill = Group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
      scale_fill_manual(values = c("Lean Control" = "#4DAF4A", "T1D" = "#E41A1C")) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 9),
            legend.position = "top",
            panel.grid.major.y = element_line(color = "gray90", size = 0.3),
            panel.grid.minor.y = element_blank()) +
      labs(title = paste("Top 30 Variables with Largest Difference in Importance for", biomarker),
           subtitle = "Between Lean Control and T1D",
           y = "Variable",
           x = "Importance",
           fill = "Group")
    
    ggsave(paste0("biomarker_comparisons/horizontal_", biomarker, ".png"),
           p_horizontal, width = 12, height = 10, dpi = 300, bg = "white")
    
    ggsave(paste0("biomarker_comparisons/horizontal_", biomarker, ".pdf"),
           p_horizontal, width = 12, height = 10)
  }
  
  cat("Individual biomarker comparison heatmaps saved in: biomarker_comparisons/\n")
  cat("  - heatmap_[biomarker].png/.pdf (3-column: LC | T1D | Difference)\n")
  cat("  - barplot_[biomarker].png/.pdf (side-by-side bars)\n")
  
  # ============================================================
  # 5. OVERALL IMPORTANCE COMPARISON ACROSS ALL BIOMARKERS
  # ============================================================
  
  cat("\n\n=== CREATING OVERALL IMPORTANCE COMPARISONS ===\n")
  
  # Calculate average importance across all biomarkers for each variable
  lean_overall_importance <- rowMeans(abs(lean_subset))
  t1d_overall_importance <- rowMeans(abs(t1d_subset))
  
  # Create overall comparison data frame
  overall_comparison <- data.frame(
    Variable = common_vars,
    Lean_Avg_Importance = lean_overall_importance,
    T1D_Avg_Importance = t1d_overall_importance,
    Difference = t1d_overall_importance - lean_overall_importance,
    Abs_Difference = abs(t1d_overall_importance - lean_overall_importance)
  )
  
  # Sort by combined importance
  overall_comparison <- overall_comparison %>%
    mutate(Combined_Importance = (Lean_Avg_Importance + T1D_Avg_Importance) / 2) %>%
    arrange(desc(Combined_Importance))
  
  # Save full comparison
  write.csv(overall_comparison, "overall_importance_comparison.csv", row.names = FALSE)
  
  # Top 50 overall comparison heatmap
  top_50_overall <- head(overall_comparison, 50)
  
  overall_matrix <- as.matrix(top_50_overall[, c("Lean_Avg_Importance", "T1D_Avg_Importance", "Difference")])
  rownames(overall_matrix) <- top_50_overall$Variable
  colnames(overall_matrix) <- c("Lean Control", "T1D", "T1D - Lean")
  
  annotation_col_overall <- data.frame(
    Type = c("Lean Control", "T1D", "Difference")
  )
  rownames(annotation_col_overall) <- colnames(overall_matrix)
  
  ann_colors_overall <- list(
    Type = c("Lean Control" = "#4DAF4A", "T1D" = "#E41A1C", "Difference" = "#984EA3")
  )
  
  pheatmap(overall_matrix,
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = "Top 50 Variables: Overall Importance (Averaged Across All Biomarkers)",
           fontsize_row = 7,
           fontsize_col = 10,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           annotation_col = annotation_col_overall,
           annotation_colors = ann_colors_overall,
           display_numbers = matrix(sprintf("%.3f", overall_matrix), nrow = nrow(overall_matrix)),
           number_color = "black",
           fontsize_number = 6,
           filename = "overall_importance_top50_heatmap.png",
           width = 10,
           height = 14)
  
  pheatmap(overall_matrix,
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = "Top 50 Variables: Overall Importance (Averaged Across All Biomarkers)",
           fontsize_row = 7,
           fontsize_col = 10,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           annotation_col = annotation_col_overall,
           annotation_colors = ann_colors_overall,
           display_numbers = matrix(sprintf("%.3f", overall_matrix), nrow = nrow(overall_matrix)),
           number_color = "black",
           fontsize_number = 6,
           filename = "overall_importance_top50_heatmap.pdf",
           width = 10,
           height = 14)
  
  cat("Created: overall_importance_top50_heatmap.png/.pdf\n")
  
  # Scatter plot: Lean vs T1D importance
  top_100_scatter <- head(overall_comparison, 100)
  
  p_scatter <- ggplot(top_100_scatter, aes(x = Lean_Avg_Importance, y = T1D_Avg_Importance)) +
    geom_point(aes(color = Abs_Difference, size = Combined_Importance), alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_text(data = top_100_scatter[1:30, ], 
              aes(label = Variable), 
              size = 2.8, hjust = -0.1, vjust = 0.5, check_overlap = T) +
    scale_color_gradient(low = "lightblue", high = "darkred", name = "Abs(Difference)") +
    scale_size_continuous(name = "Combined\nImportance", range = c(2, 8)) +
    theme_bw() +
    labs(title = "Overall Variable Importance: Lean Control vs T1D",
         subtitle = "Top 100 variables (Top 15 labeled). Points above diagonal line = higher in T1D",
         x = "Average Importance in Lean Control",
         y = "Average Importance in T1D") +
    coord_fixed()
  
  ggsave("overall_importance_scatter.png", p_scatter, width = 12, height = 10, dpi = 300, bg = "white")
  ggsave("overall_importance_scatter.pdf", p_scatter, width = 12, height = 10)
  
  cat("Created: overall_importance_scatter.png/.pdf\n")
  
  # Bar plot showing variables most different between groups
  top_30_different <- overall_comparison %>%
    arrange(desc(Abs_Difference)) %>%
    head(30)
  
  top_30_different_long <- top_30_different %>%
    dplyr::select(Variable, Lean_Avg_Importance, T1D_Avg_Importance) %>%
    pivot_longer(cols = c(Lean_Avg_Importance, T1D_Avg_Importance),
                 names_to = "Group",
                 values_to = "Importance") %>%
    mutate(Group = recode(Group, 
                          "Lean_Avg_Importance" = "Lean Control",
                          "T1D_Avg_Importance" = "T1D"))
  
  p_different <- ggplot(top_30_different_long, aes(x = Importance, y = reorder(Variable, Importance), fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
    scale_fill_manual(values = c("Lean Control" = "#4DAF4A", "T1D" = "#E41A1C")) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          legend.position = "top",
          panel.grid.major.y = element_line(color = "gray90", size = 0.3),
          panel.grid.minor.y = element_blank()) +
    labs(title = "Top 30 Variables with Largest Difference in Importance",
         subtitle = "Between Lean Control and T1D (averaged across all biomarkers)",
         y = "Variable",
         x = "Average Importance",
         fill = "Group")
  
  ggsave("overall_most_different_variables.png", p_different, width = 12, height = 10, dpi = 300, bg = "white")
  ggsave("overall_most_different_variables.pdf", p_different, width = 12, height = 10)
  
  cat("Created: overall_most_different_variables.png/.pdf\n")
  
  # Ranking comparison table
  lean_ranks <- data.frame(
    Variable = names(sort(lean_overall_importance, decreasing = TRUE)),
    Lean_Rank = 1:length(lean_overall_importance)
  )
  
  t1d_ranks <- data.frame(
    Variable = names(sort(t1d_overall_importance, decreasing = TRUE)),
    T1D_Rank = 1:length(t1d_overall_importance)
  )
  
  rank_comparison <- merge(lean_ranks, t1d_ranks, by = "Variable")
  rank_comparison$Rank_Change <- rank_comparison$Lean_Rank - rank_comparison$T1D_Rank
  rank_comparison$Abs_Rank_Change <- abs(rank_comparison$Rank_Change)
  
  # Show top 20 biggest rank changes
  cat("\n=== TOP 20 VARIABLES WITH LARGEST RANK CHANGES ===\n")
  cat("Positive Rank_Change = variable more important in T1D (moved up in ranking)\n")
  cat("Negative Rank_Change = variable more important in Lean Control (moved down in ranking)\n\n")
  
  top_rank_changes <- rank_comparison %>%
    arrange(desc(Abs_Rank_Change)) %>%
    head(20)
  
  print(top_rank_changes)
  
  write.csv(rank_comparison, "overall_rank_comparison.csv", row.names = FALSE)
  cat("\nFull ranking comparison saved to: overall_rank_comparison.csv\n")
  
  # Summary statistics
  cat("\n=== OVERALL IMPORTANCE SUMMARY STATISTICS ===\n")
  cat("Lean Control:\n")
  cat("  Mean importance:", round(mean(lean_overall_importance), 4), "\n")
  cat("  Median importance:", round(median(lean_overall_importance), 4), "\n")
  cat("  SD importance:", round(sd(lean_overall_importance), 4), "\n\n")
  
  cat("T1D:\n")
  cat("  Mean importance:", round(mean(t1d_overall_importance), 4), "\n")
  cat("  Median importance:", round(median(t1d_overall_importance), 4), "\n")
  cat("  SD importance:", round(sd(t1d_overall_importance), 4), "\n\n")
  
  cat("Correlation between Lean and T1D importance rankings:\n")
  cat("  Spearman correlation:", round(cor(lean_overall_importance, t1d_overall_importance, method = "spearman"), 3), "\n")
  cat("  Pearson correlation:", round(cor(lean_overall_importance, t1d_overall_importance, method = "pearson"), 3), "\n")
}

# ============================================================
# PERFORMANCE COMPARISON
# ============================================================

cat("\n\n##########################################################")
cat("\n### PERFORMANCE COMPARISON: LEAN CONTROL VS T1D")
cat("\n##########################################################\n\n")

# Performance comparison
cat("=== PERFORMANCE SUMMARY ===\n")
print(all_performance)

write.csv(all_performance, "performance_summary_lean_vs_t1d.csv", row.names = FALSE)

# Performance comparison plots
p_r2 <- ggplot(all_performance, aes(x = Biomarker, y = R2, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance: R² by Biomarker (Lean Control vs T1D)",
       y = "R-squared", x = "Biomarker") +
  scale_fill_manual(values = c("Lean_Control" = "#4DAF4A", "T1D" = "#E41A1C"),
                    labels = c("Lean Control", "T1D"))

print(p_r2)
ggsave("R2_comparison_lean_vs_t1d.png", p_r2, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("R2_comparison_lean_vs_t1d.pdf", p_r2, width = 10, height = 6)

p_adj_r2 <- ggplot(all_performance, aes(x = Biomarker, y = Adj_R2, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance: Adjusted R² by Biomarker (Lean Control vs T1D)",
       y = "Adjusted R-squared", x = "Biomarker") +
  scale_fill_manual(values = c("Lean_Control" = "#4DAF4A", "T1D" = "#E41A1C"),
                    labels = c("Lean Control", "T1D"))

print(p_adj_r2)
ggsave("Adj_R2_comparison_lean_vs_t1d.png", p_adj_r2, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("Adj_R2_comparison_lean_vs_t1d.pdf", p_adj_r2, width = 10, height = 6)

# Save all results
save(all_results, all_performance, all_variable_importance,
     file = "complete_analysis_results_lean_vs_t1d.RData")

cat("\n\n##########################################################")
cat("\n### ANALYSIS COMPLETE!")
cat("\n##########################################################\n")
cat("\nKey outputs generated:\n")
cat("1. Side-by-side comparison heatmaps (top 20 and top 50)\n")
cat("2. Difference heatmap (Lean - T1D)\n")
cat("3. Individual biomarker comparison plots\n")
cat("4. Performance comparison plots (R², Adjusted R²)\n")
cat("5. CSV files with variable importance and differences\n")
cat("6. Complete analysis results saved in RData file\n")
