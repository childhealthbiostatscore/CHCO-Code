


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













#### PCA Analysis instead 
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
# Define analysis groups
analysis_groups <- list(
  All = data_varying,  # All participants
  T1D = data_varying %>% filter(group == "Type 1 Diabetes")  # Only T1D
)

cat("Sample sizes:\n")
cat("All participants:", nrow(analysis_groups$All), "\n")
cat("Type 1 Diabetes only:", nrow(analysis_groups$T1D), "\n\n")

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
  
  # ============================================================
  # OPTION 3: Multiple heatmaps by magnitude
  # ============================================================
  
  overall_importance <- rowMeans(abs(importance_matrix))
  
  # Heatmap 1: Top 15 highest importance (raw scale with values) - PNG
  top_15_vars <- names(sort(overall_importance, decreasing = TRUE)[1:min(15, length(overall_importance))])
  pheatmap(importance_matrix[top_15_vars, ], 
           scale = "none",
           main = paste("Top 15 Highest Importance:", group_name),
           fontsize_row = 10,
           fontsize_number = 8,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           display_numbers = matrix(sprintf("%.2f", importance_matrix[top_15_vars, ]), 
                                    nrow = length(top_15_vars)),
           number_color = "black",
           filename = paste0("heatmap_top15_", group_name, ".png"),
           width = 10,
           height = 7)
  
  # Also save as PDF
  pheatmap(importance_matrix[top_15_vars, ], 
           scale = "none",
           main = paste("Top 15 Highest Importance:", group_name),
           fontsize_row = 10,
           fontsize_number = 8,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           display_numbers = matrix(sprintf("%.2f", importance_matrix[top_15_vars, ]), 
                                    nrow = length(top_15_vars)),
           number_color = "black",
           filename = paste0("heatmap_top15_", group_name, ".pdf"),
           width = 10,
           height = 7)
  
  # Heatmap 2: Variables ranked 16-50 (better visible without extreme values)
  if(length(overall_importance) >= 50) {
    mid_vars <- names(sort(overall_importance, decreasing = TRUE)[16:50])
    pheatmap(importance_matrix[mid_vars, ], 
             scale = "none",
             main = paste("Variables Ranked 16-50:", group_name),
             fontsize_row = 7,
             angle_col = 45,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             filename = paste0("heatmap_mid_", group_name, ".png"),
             width = 10,
             height = 12)
    
    pheatmap(importance_matrix[mid_vars, ], 
             scale = "none",
             main = paste("Variables Ranked 16-50:", group_name),
             fontsize_row = 7,
             angle_col = 45,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             filename = paste0("heatmap_mid_", group_name, ".pdf"),
             width = 10,
             height = 12)
  } else if(length(overall_importance) > 15) {
    # If less than 50 variables total, show 16 to end
    mid_vars <- names(sort(overall_importance, decreasing = TRUE)[16:length(overall_importance)])
    pheatmap(importance_matrix[mid_vars, ], 
             scale = "none",
             main = paste("Variables Ranked 16+:", group_name),
             fontsize_row = 8,
             angle_col = 45,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             filename = paste0("heatmap_mid_", group_name, ".png"),
             width = 10,
             height = 10)
    
    pheatmap(importance_matrix[mid_vars, ], 
             scale = "none",
             main = paste("Variables Ranked 16+:", group_name),
             fontsize_row = 8,
             angle_col = 45,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             filename = paste0("heatmap_mid_", group_name, ".pdf"),
             width = 10,
             height = 10)
  }
  
  # Heatmap 3: All variables with log scale
  all_vars <- rownames(importance_matrix)
  plot_data_log <- sign(importance_matrix) * log10(abs(importance_matrix) + 1)
  pheatmap(plot_data_log[all_vars, ], 
           scale = "none",
           main = paste("All Variables (log10 scale):", group_name),
           fontsize_row = 6,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           filename = paste0("heatmap_all_log_", group_name, ".png"),
           width = 10,
           height = 16)
  
  pheatmap(plot_data_log[all_vars, ], 
           scale = "none",
           main = paste("All Variables (log10 scale):", group_name),
           fontsize_row = 6,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           filename = paste0("heatmap_all_log_", group_name, ".pdf"),
           width = 10,
           height = 16)
  
  cat("\nMultiple heatmaps saved for", group_name, ":\n")
  cat("  - heatmap_top15_", group_name, ".png/.pdf (with values)\n", sep = "")
  if(length(overall_importance) > 15) {
    cat("  - heatmap_mid_", group_name, ".png/.pdf\n", sep = "")
  }
  cat("  - heatmap_all_log_", group_name, ".png/.pdf\n", sep = "")
  
  # ============================================================
  # INDIVIDUAL BIOMARKER HEATMAPS
  # ============================================================
  
  cat("\nCreating individual biomarker heatmaps for", group_name, "...\n")
  
  # Create a directory for individual biomarker heatmaps
  dir.create(paste0("biomarker_heatmaps_", group_name), showWarnings = FALSE)
  
  for(biomarker in qx_var) {
    cat("  Processing", biomarker, "...\n")
    
    # Get importance for this biomarker
    biomarker_importance <- importance_matrix[, biomarker, drop = FALSE]
    
    # Sort by absolute importance
    biomarker_importance_sorted <- biomarker_importance[order(abs(biomarker_importance[, 1]), decreasing = TRUE), , drop = FALSE]
    
    # Top 30 variables for this biomarker
    n_top <- min(30, nrow(biomarker_importance_sorted))
    top_vars <- rownames(biomarker_importance_sorted)[1:n_top]
    
    # Create data frame for plotting
    plot_df <- data.frame(
      Variable = factor(top_vars, levels = rev(top_vars)),  # Reverse for plotting
      Importance = biomarker_importance_sorted[top_vars, 1],
      Abs_Importance = abs(biomarker_importance_sorted[top_vars, 1])
    )
    
    # Bar plot - PNG
    p <- ggplot(plot_df, aes(x = Variable, y = Importance, fill = Importance)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = 0,
                           name = "Importance") +
      coord_flip() +
      theme_bw() +
      theme(axis.text.y = element_text(size = 9)) +
      labs(title = paste("Top 30 Variables for", biomarker, "-", group_name),
           subtitle = paste("Based on PCA-regression variable importance"),
           x = "Variable",
           y = "Variable Importance") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
    
    ggsave(paste0("biomarker_heatmaps_", group_name, "/barplot_", biomarker, "_", group_name, ".png"),
           p, width = 10, height = 8, dpi = 300, bg = "white")
    
    ggsave(paste0("biomarker_heatmaps_", group_name, "/barplot_", biomarker, "_", group_name, ".pdf"),
           p, width = 10, height = 8)
    
    # Also create a mini heatmap for just this biomarker (top 50)
    n_heatmap <- min(50, nrow(biomarker_importance_sorted))
    heatmap_vars <- rownames(biomarker_importance_sorted)[1:n_heatmap]
    
    pheatmap(biomarker_importance_sorted[heatmap_vars, , drop = FALSE],
             cluster_rows = FALSE,  # Keep sorted by importance
             cluster_cols = FALSE,
             scale = "none",
             main = paste("Top", n_heatmap, "Variables for", biomarker, "-", group_name),
             fontsize_row = 7,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             display_numbers = matrix(sprintf("%.3f", biomarker_importance_sorted[heatmap_vars, 1]), 
                                      ncol = 1),
             number_color = "black",
             fontsize_number = 6,
             filename = paste0("biomarker_heatmaps_", group_name, "/heatmap_", biomarker, "_", group_name, ".png"),
             width = 6,
             height = 12)
    
    pheatmap(biomarker_importance_sorted[heatmap_vars, , drop = FALSE],
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             scale = "none",
             main = paste("Top", n_heatmap, "Variables for", biomarker, "-", group_name),
             fontsize_row = 7,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             display_numbers = matrix(sprintf("%.3f", biomarker_importance_sorted[heatmap_vars, 1]), 
                                      ncol = 1),
             number_color = "black",
             fontsize_number = 6,
             filename = paste0("biomarker_heatmaps_", group_name, "/heatmap_", biomarker, "_", group_name, ".pdf"),
             width = 6,
             height = 12)
  }
  
  cat("Individual biomarker heatmaps saved in:", paste0("biomarker_heatmaps_", group_name, "/\n"))
}

# ============================================================
# COMPARE RESULTS ACROSS GROUPS
# ============================================================

cat("\n\n##########################################################")
cat("\n### COMPARISON ACROSS GROUPS")
cat("\n##########################################################\n\n")

# Performance comparison
cat("=== PERFORMANCE SUMMARY ===\n")
print(all_performance)

write.csv(all_performance, "performance_summary_all_groups.csv", row.names = FALSE)

# Performance comparison plots
library(ggplot2)

p_r2 <- ggplot(all_performance, aes(x = Biomarker, y = R2, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance: R² by Biomarker and Group",
       y = "R-squared", x = "Biomarker") +
  scale_fill_brewer(palette = "Set1")

print(p_r2)
ggsave("R2_comparison.png", p_r2, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("R2_comparison.pdf", p_r2, width = 10, height = 6)

p_adj_r2 <- ggplot(all_performance, aes(x = Biomarker, y = Adj_R2, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance: Adjusted R² by Biomarker and Group",
       y = "Adjusted R-squared", x = "Biomarker") +
  scale_fill_brewer(palette = "Set1")

print(p_adj_r2)
ggsave("Adj_R2_comparison.png", p_adj_r2, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("Adj_R2_comparison.pdf", p_adj_r2, width = 10, height = 6)

p_rmse <- ggplot(all_performance, aes(x = Biomarker, y = RMSE, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance: RMSE by Biomarker and Group (Standardized)",
       y = "Root Mean Squared Error (Standardized Units)", x = "Biomarker") +
  scale_fill_brewer(palette = "Set1")

print(p_rmse)
ggsave("RMSE_comparison.png", p_rmse, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("RMSE_comparison.pdf", p_rmse, width = 10, height = 6)

# Compare variable importance between groups
cat("\n=== COMPARING VARIABLE IMPORTANCE BETWEEN GROUPS ===\n")

for(biomarker in qx_var) {
  cat("\n--- ", biomarker, " ---\n")
  
  all_imp <- all_variable_importance$All[[biomarker]]
  t1d_imp <- all_variable_importance$T1D[[biomarker]]
  
  # Find common variables
  common_vars <- intersect(rownames(all_imp), rownames(t1d_imp))
  
  if(length(common_vars) == 0) {
    cat("WARNING: No common variables between groups for", biomarker, "\n")
    next
  }
  
  # Combine and compare
  comparison_df <- data.frame(
    Variable = common_vars,
    All_Group = all_imp[common_vars, ],
    T1D_Group = t1d_imp[common_vars, ],
    Difference = all_imp[common_vars, ] - t1d_imp[common_vars, ],
    Abs_Difference = abs(all_imp[common_vars, ] - t1d_imp[common_vars, ])
  )
  
  # Show top differences
  comparison_df <- comparison_df[order(-comparison_df$Abs_Difference), ]
  
  cat("\nTop 15 variables with largest importance differences:\n")
  print(head(comparison_df, 15))
  
  write.csv(comparison_df, 
            paste0("importance_comparison_", biomarker, ".csv"),
            row.names = FALSE)
}

# Save all results
save(all_results, all_performance, all_variable_importance,
     file = "complete_analysis_results.RData")


