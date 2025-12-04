### Brain Biomarker Analysis - Clinical Variable Associations
### Simple regressions: Clinical variables (HbA1c, M/I, CGM, BP) predicting brain biomarkers
### Models: Unadjusted and adjusted (+ age + sex + BMI)
### Groups: All participants, T1D, T2D, Lean Control (by sample type)

library(tidyverse)
library(broom)
library(patchwork)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(knitr)
library(kableExtra)
library(moments)  # For skewness calculation

# Define Quanterix brain biomarkers
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Define clinical predictors of interest
clinical_predictors <- list(
  "Glycemic Control" = c("hba1c", "fbg"),
  "Insulin Sensitivity" = c("avg_m_fsoc", "homa_ir", "adipose_ir", "search_eis"),
  "CGM Metrics" = c("cgm_mean_glucose", "cgm_sd", "cgm_cv", "time_in_range", 
                    "time_above_range", "time_below_range"),
  "Blood Pressure" = c("sbp", "dbp", "map")  # Add actual BP variable names from your data
)

# Flatten the list for easier use
all_predictors <- unlist(clinical_predictors, use.names = FALSE)

# ============================================================
# DATA LOADING AND PREPARATION
# ============================================================

cat("\n##########################################################")
cat("\n### LOADING AND PREPARING DATA")
cat("\n##########################################################\n\n")

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

# Summarize data by participant and visit
dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

# ============================================================
# DETERMINE SAMPLE TYPE (SERUM/PLASMA) BY STUDY
# ============================================================

cat("Determining sample type by study/record_id prefix...\n")

dat <- dat %>%
  mutate(
    sample_type = case_when(
      # CROCODILE participants (CRC IDs) - plasma
      grepl("^CRC", record_id, ignore.case = TRUE) ~ "plasma",
      # PENGUIN participants (PEN IDs) - plasma
      grepl("^PEN", record_id, ignore.case = TRUE) ~ "plasma",
      # RH2 participants - serum
      study == "RH2" ~ "serum",
      # Default to study-based assignment
      study == "CROCODILE" ~ "plasma",
      study == "PENGUIN" ~ "plasma",
      TRUE ~ NA_character_
    )
  )

# ============================================================
# FILTER DATA - BASELINE WITH QUANTERIX DATA
# ============================================================

dat_baseline <- dat %>% 
  filter(visit == 'baseline') %>%
  filter(!is.na(ab40_avg_conc)) %>%  # Has Quanterix data
  filter(!is.na(sample_type)) %>%     # Known sample type
  filter(!is.na(age) & !is.na(sex) & !is.na(bmi))  # Has covariates

cat("\nParticipants with Quanterix data and covariates at baseline:\n")
cat("  Total:", nrow(dat_baseline), "\n")
cat("  Plasma:", sum(dat_baseline$sample_type == "plasma", na.rm = TRUE), "\n")
cat("  Serum:", sum(dat_baseline$sample_type == "serum", na.rm = TRUE), "\n\n")

# Check group distribution by sample type
cat("Group distribution by sample type:\n")
print(table(dat_baseline$group, dat_baseline$sample_type))
cat("\n")

# ============================================================
# SETUP OUTPUT DIRECTORY
# ============================================================

setwd('C:/Users/netio/Documents/UofW/Projects/Maninder_Data/Clinical_Biomarker_Associations/')
dir.create(getwd(), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# DEFINE ANALYSIS GROUPS
# ============================================================

# Determine which sample type to use (use the one with more participants)
sample_type_counts <- table(dat_baseline$sample_type)
primary_sample_type <- names(which.max(sample_type_counts))

cat("Using sample type:", primary_sample_type, "(N =", max(sample_type_counts), ")\n\n")

# Filter to primary sample type only
dat_analysis <- dat_baseline %>%
  filter(sample_type == primary_sample_type)

# Define analysis groups
analysis_groups <- list(
  "All_Participants" = dat_analysis,
  "Type_1_Diabetes" = dat_analysis %>% filter(group == "Type 1 Diabetes"),
  "Type_2_Diabetes" = dat_analysis %>% filter(group == "Type 2 Diabetes"),
  "Lean_Control" = dat_analysis %>% filter(group == "Lean Control")
)

# Report sample sizes
cat("Sample sizes by group:\n")
for(group_name in names(analysis_groups)) {
  cat(sprintf("  %-20s: N = %d\n", group_name, nrow(analysis_groups[[group_name]])))
}
cat("\n")

# ============================================================
# FUNCTION TO TEST LOG TRANSFORMATION
# ============================================================

test_log_transformation <- function(x) {
  # Test if log transformation improves normality
  # Returns TRUE if log transformation is beneficial
  
  # Skip if any non-positive values
  if(any(x <= 0, na.rm = TRUE)) {
    return(FALSE)
  }
  
  # Skip if not enough data
  if(sum(!is.na(x)) < 10) {
    return(FALSE)
  }
  
  # Check skewness
  skew_raw <- moments::skewness(x, na.rm = TRUE)
  skew_log <- moments::skewness(log(x), na.rm = TRUE)
  
  # If raw data is highly right-skewed (>1) and log reduces skewness
  if(abs(skew_raw) > 1 && abs(skew_log) < abs(skew_raw)) {
    return(TRUE)
  }
  
  return(FALSE)
}

# ============================================================
# FUNCTION TO RUN REGRESSIONS WITH TRANSFORMATION TESTING AND SCALING
# ============================================================

run_regression_analysis <- function(data, predictor, outcome, covariates = NULL, group_name = "") {
  
  # Check if predictor and outcome exist and have sufficient data
  if(!predictor %in% names(data) || !outcome %in% names(data)) {
    return(NULL)
  }
  
  # Create analysis dataset
  if(is.null(covariates)) {
    analysis_data <- data %>% 
      dplyr::select(all_of(c(predictor, outcome))) %>%
      filter(complete.cases(.))
    
    if(nrow(analysis_data) < 10) return(NULL)
    
    model_type <- "Unadjusted"
    covar_list <- NULL
    
  } else {
    analysis_data <- data %>% 
      dplyr::select(all_of(c(predictor, outcome, covariates))) %>%
      filter(complete.cases(.))
    
    if(nrow(analysis_data) < 10) return(NULL)
    
    model_type <- "Adjusted"
    covar_list <- covariates
  }
  
  # Test log transformation for predictor and outcome
  pred_needs_log <- test_log_transformation(analysis_data[[predictor]])
  outcome_needs_log <- test_log_transformation(analysis_data[[outcome]])
  
  # Create working dataset with transformations
  working_data <- analysis_data
  
  pred_var_name <- predictor
  outcome_var_name <- outcome
  
  if(pred_needs_log) {
    working_data[[predictor]] <- log(working_data[[predictor]])
    pred_var_name <- paste0("log(", predictor, ")")
  }
  
  if(outcome_needs_log) {
    working_data[[outcome]] <- log(working_data[[outcome]])
    outcome_var_name <- paste0("log(", outcome, ")")
  }
  
  # Standardize all variables for standardized coefficients
  # Store original values for model fitting
  original_working_data <- working_data
  
  # Standardize predictor and outcome
  working_data[[predictor]] <- scale(working_data[[predictor]])[,1]
  working_data[[outcome]] <- scale(working_data[[outcome]])[,1]
  
  # Standardize covariates if present (except sex which is binary)
  if(!is.null(covar_list)) {
    for(cov in covar_list) {
      if(cov != "sex" && is.numeric(working_data[[cov]])) {
        working_data[[cov]] <- scale(working_data[[cov]])[,1]
      }
    }
  }
  
  # Build formula
  if(is.null(covar_list)) {
    formula_str <- paste0(outcome, " ~ ", predictor)
  } else {
    formula_str <- paste0(outcome, " ~ ", predictor, " + ", paste(covar_list, collapse = " + "))
  }
  
  # Fit standardized model
  model_std <- lm(as.formula(formula_str), data = working_data)
  model_summary_std <- summary(model_std)
  coef_table_std <- coef(model_summary_std)
  
  # Get predictor coefficient (standardized)
  predictor_row <- which(rownames(coef_table_std) == predictor)
  
  if(length(predictor_row) == 0) return(NULL)
  
  # Also fit unstandardized model for raw beta
  model_raw <- lm(as.formula(formula_str), data = original_working_data)
  model_summary_raw <- summary(model_raw)
  coef_table_raw <- coef(model_summary_raw)
  
  results <- data.frame(
    Group = group_name,
    Sample_Type = unique(data$sample_type)[1],
    Predictor = predictor,
    Outcome = outcome,
    Predictor_Used = pred_var_name,
    Outcome_Used = outcome_var_name,
    Predictor_Log_Transformed = pred_needs_log,
    Outcome_Log_Transformed = outcome_needs_log,
    Model_Type = model_type,
    N = nrow(analysis_data),
    Beta_Raw = coef_table_raw[predictor_row, "Estimate"],
    Beta_Std = coef_table_std[predictor_row, "Estimate"],  # Standardized beta
    SE_Raw = coef_table_raw[predictor_row, "Std. Error"],
    SE_Std = coef_table_std[predictor_row, "Std. Error"],
    t_value = coef_table_std[predictor_row, "t value"],
    p_value = coef_table_std[predictor_row, "Pr(>|t|)"],
    R_squared = model_summary_std$r.squared,
    Adj_R_squared = model_summary_std$adj.r.squared,
    Model_p_value = pf(model_summary_std$fstatistic[1], 
                       model_summary_std$fstatistic[2], 
                       model_summary_std$fstatistic[3], 
                       lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

# ============================================================
# RUN ALL REGRESSIONS
# ============================================================

cat("##########################################################")
cat("\n### RUNNING REGRESSION ANALYSES")
cat("\n##########################################################\n\n")

all_results <- data.frame()

for(group_name in names(analysis_groups)) {
  
  cat("\n=== Analyzing group:", group_name, "===\n")
  current_data <- analysis_groups[[group_name]]
  
  if(nrow(current_data) < 10) {
    cat("  Skipping - insufficient sample size (N =", nrow(current_data), ")\n")
    next
  }
  
  # Check which predictors are available in this dataset
  available_predictors <- intersect(all_predictors, names(current_data))
  
  if(length(available_predictors) == 0) {
    cat("  No clinical predictors available in dataset\n")
    next
  }
  
  cat("  Available predictors:", length(available_predictors), "\n")
  cat("  Sample size:", nrow(current_data), "\n")
  
  # Run regressions for each predictor-outcome combination
  for(predictor in available_predictors) {
    
    # Check if predictor has sufficient non-missing data
    n_available <- sum(!is.na(current_data[[predictor]]))
    if(n_available < 10) next
    
    for(outcome in qx_var) {
      
      # Unadjusted model
      result_unadj <- run_regression_analysis(
        data = current_data,
        predictor = predictor,
        outcome = outcome,
        covariates = NULL,
        group_name = group_name
      )
      
      if(!is.null(result_unadj)) {
        all_results <- rbind(all_results, result_unadj)
      }
      
      # Adjusted model (+ age + sex + BMI)
      result_adj <- run_regression_analysis(
        data = current_data,
        predictor = predictor,
        outcome = outcome,
        covariates = c("age", "sex", "bmi"),
        group_name = group_name
      )
      
      if(!is.null(result_adj)) {
        all_results <- rbind(all_results, result_adj)
      }
    }
  }
  
  cat("  Completed", nrow(all_results), "total regressions so far\n")
}

# Add significance flags
all_results <- all_results %>%
  mutate(
    Significant = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.10 ~ ".",
      TRUE ~ ""
    ),
    Sig_Flag = p_value < 0.05,
    # Create interpretable transformation label
    Transformation = case_when(
      Predictor_Log_Transformed & Outcome_Log_Transformed ~ "Both log",
      Predictor_Log_Transformed ~ "Predictor log",
      Outcome_Log_Transformed ~ "Outcome log",
      TRUE ~ "None"
    )
  )

cat("\n\nTotal regressions completed:", nrow(all_results), "\n")
cat("Significant associations (p < 0.05):", sum(all_results$Sig_Flag), "\n")

# Summary of transformations used
cat("\nTransformation summary:\n")
transformation_summary <- table(all_results$Transformation, all_results$Model_Type)
print(transformation_summary)
cat("\n")

# Save all results
write.csv(all_results, "all_regression_results.csv", row.names = FALSE)

# ============================================================
# CREATE SUMMARY TABLES
# ============================================================

cat("##########################################################")
cat("\n### CREATING SUMMARY TABLES AND VISUALIZATIONS")
cat("\n##########################################################\n\n")

# Summary of significant findings
significant_results <- all_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

write.csv(significant_results, "significant_associations.csv", row.names = FALSE)

cat("Significant associations saved to: significant_associations.csv\n")
cat("Number of significant associations:", nrow(significant_results), "\n\n")

# Print top 20 most significant findings
if(nrow(significant_results) > 0) {
  cat("\n=== TOP 20 MOST SIGNIFICANT ASSOCIATIONS ===\n")
  top_20 <- significant_results %>%
    head(20) %>%
    dplyr::select(Group, Model_Type, Predictor, Outcome, Transformation, N, Beta_Std, p_value, R_squared)
  
  print(kable(top_20, digits = c(0, 0, 0, 0, 0, 0, 3, 4, 3), format = "simple"))
}

# ============================================================
# HEATMAPS OF P-VALUES
# ============================================================

cat("\n\n=== Creating heatmaps ===\n")

# Function to create p-value heatmap
create_pvalue_heatmap <- function(results_subset, title, filename) {
  
  # Pivot to matrix format
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  # Transform p-values for visualization (-log10)
  pval_transformed <- -log10(pval_matrix)
  
  # Cap at -log10(0.001) for visualization
  pval_transformed[pval_transformed > 3] <- 3
  
  # Only create heatmap if we have data
  if(nrow(pval_transformed) > 0 && ncol(pval_transformed) > 0) {
    
    pheatmap(pval_transformed,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
             breaks = seq(0, 3, length.out = 101),
             main = title,
             fontsize_row = 8,
             fontsize_col = 9,
             angle_col = 45,
             display_numbers = matrix(sprintf("%.3f", pval_matrix), nrow = nrow(pval_matrix)),
             number_color = "black",
             fontsize_number = 6,
             legend_breaks = c(0, 0.5, 1, 1.3, 2, 3),
             legend_labels = c("1.0", "0.32", "0.10", "0.05", "0.01", "0.001"),
             filename = filename,
             width = 10,
             height = max(6, nrow(pval_transformed) * 0.3))
    
    return(TRUE)
  }
  return(FALSE)
}

# Create heatmaps for each group and model type
for(group_name in unique(all_results$Group)) {
  for(model_type in c("Unadjusted", "Adjusted")) {
    
    results_subset <- all_results %>%
      filter(Group == group_name, Model_Type == model_type)
    
    if(nrow(results_subset) > 0) {
      title <- paste0(group_name, " - ", model_type, " Models\n(-log10 p-values)")
      filename <- paste0("heatmap_pvalues_", gsub(" ", "_", group_name), "_", model_type, ".png")
      
      success <- create_pvalue_heatmap(results_subset, title, filename)
      if(success) {
        cat("Created:", filename, "\n")
      }
    }
  }
}

# ============================================================
# EFFECT SIZE HEATMAPS (BETA COEFFICIENTS)
# ============================================================

cat("\n=== Creating beta coefficient heatmaps ===\n")

create_beta_heatmap <- function(results_subset, title, filename) {
  
  # Pivot to matrix format - use standardized betas
  beta_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, Beta_Std) %>%
    pivot_wider(names_from = Outcome, values_from = Beta_Std) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  # Get corresponding p-values for significance markers
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  # Create significance markers
  sig_matrix <- matrix("", nrow = nrow(beta_matrix), ncol = ncol(beta_matrix))
  sig_matrix[pval_matrix < 0.001] <- "***"
  sig_matrix[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
  sig_matrix[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"
  
  # Combine beta and significance
  display_matrix <- matrix(
    paste0(sprintf("%.3f", beta_matrix), "\n", sig_matrix),
    nrow = nrow(beta_matrix)
  )
  
  if(nrow(beta_matrix) > 0 && ncol(beta_matrix) > 0) {
    
    pheatmap(beta_matrix,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             main = paste0(title, "\n(Standardized Beta Coefficients)"),
             fontsize_row = 8,
             fontsize_col = 9,
             angle_col = 45,
             display_numbers = display_matrix,
             number_color = "black",
             fontsize_number = 5,
             filename = filename,
             width = 10,
             height = max(6, nrow(beta_matrix) * 0.3))
    
    return(TRUE)
  }
  return(FALSE)
}

for(group_name in unique(all_results$Group)) {
  for(model_type in c("Unadjusted", "Adjusted")) {
    
    results_subset <- all_results %>%
      filter(Group == group_name, Model_Type == model_type)
    
    if(nrow(results_subset) > 0) {
      title <- paste0(group_name, " - ", model_type, " Models\n(Beta Coefficients)")
      filename <- paste0("heatmap_beta_", gsub(" ", "_", group_name), "_", model_type, ".png")
      
      success <- create_beta_heatmap(results_subset, title, filename)
      if(success) {
        cat("Created:", filename, "\n")
      }
    }
  }
}

# ============================================================
# SCATTER PLOTS FOR SIGNIFICANT ASSOCIATIONS
# ============================================================

cat("\n=== Creating scatter plots for top significant associations ===\n")

# Get top 12 most significant associations for plotting
top_for_plots <- significant_results %>%
  arrange(p_value) %>%
  head(12)

if(nrow(top_for_plots) > 0) {
  
  plot_list <- list()
  
  for(i in 1:nrow(top_for_plots)) {
    
    row <- top_for_plots[i, ]
    
    # Get the data for this association
    plot_data <- analysis_groups[[row$Group]] %>%
      dplyr::select(x = all_of(row$Predictor), y = all_of(row$Outcome)) %>%
      filter(complete.cases(.))
    
    if(nrow(plot_data) < 10) next
    
    # Apply transformations if needed
    x_label <- row$Predictor
    y_label <- row$Outcome
    
    if(row$Predictor_Log_Transformed) {
      plot_data$x <- log(plot_data$x)
      x_label <- paste0("log(", x_label, ")")
    }
    
    if(row$Outcome_Log_Transformed) {
      plot_data$y <- log(plot_data$y)
      y_label <- paste0("log(", y_label, ")")
    }
    
    # Create plot
    p <- ggplot(plot_data, aes(x = x, y = y)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "blue") +
      labs(
        title = paste0(row$Predictor, " → ", row$Outcome),
        subtitle = sprintf("%s | %s | β_std=%.3f, p=%s, R²=%.3f",
                           row$Group, row$Model_Type, row$Beta_Std, 
                           format.pval(row$p_value, digits = 2), row$R_squared),
        x = x_label,
        y = y_label
      ) +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    plot_list[[i]] <- p
  }
  
  # Combine plots
  if(length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 3)
    ggsave("top_significant_associations_scatter.png", combined_plot,
           width = 15, height = 12, dpi = 300, bg = "white")
    cat("Created: top_significant_associations_scatter.png\n")
  }
}

# ============================================================
# SUMMARY BY PREDICTOR TYPE
# ============================================================

cat("\n\n=== Summary of associations by predictor category ===\n")

# Add predictor category
all_results <- all_results %>%
  mutate(
    Predictor_Category = case_when(
      Predictor %in% clinical_predictors$`Glycemic Control` ~ "Glycemic Control",
      Predictor %in% clinical_predictors$`Insulin Sensitivity` ~ "Insulin Sensitivity",
      Predictor %in% clinical_predictors$`CGM Metrics` ~ "CGM Metrics",
      Predictor %in% clinical_predictors$`Blood Pressure` ~ "Blood Pressure",
      TRUE ~ "Other"
    )
  )

category_summary <- all_results %>%
  filter(p_value < 0.05) %>%
  group_by(Predictor_Category, Group, Model_Type) %>%
  summarise(
    N_Significant = n(),
    Mean_Abs_Beta = mean(abs(Beta)),
    Mean_R_squared = mean(R_squared),
    .groups = "drop"
  ) %>%
  arrange(desc(N_Significant))

write.csv(category_summary, "summary_by_predictor_category.csv", row.names = FALSE)

cat("\n")
print(kable(category_summary, digits = 3, format = "simple"))

# ============================================================
# FOREST PLOTS - COMPARE EFFECT SIZES ACROSS GROUPS
# ============================================================

cat("\n\n##########################################################")
cat("\n### CREATING FOREST PLOTS FOR CROSS-GROUP COMPARISONS")
cat("\n##########################################################\n\n")

# Focus on adjusted models for cleaner comparison
adjusted_results <- all_results %>%
  filter(Model_Type == "Adjusted")

# Find predictor-outcome combinations that appear in multiple groups
combo_counts <- adjusted_results %>%
  group_by(Predictor, Outcome) %>%
  summarise(N_Groups = n(), .groups = "drop") %>%
  filter(N_Groups >= 2)  # At least in 2 groups

cat("Found", nrow(combo_counts), "predictor-outcome pairs tested in multiple groups\n\n")

# Calculate which associations are significant in T1D or T2D but not in controls
disease_specific <- adjusted_results %>%
  filter(Group %in% c("Type_1_Diabetes", "Type_2_Diabetes", "Lean_Control")) %>%
  dplyr::select(Group, Predictor, Outcome, Beta_Std, p_value, SE_Std) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(Beta_Std, p_value, SE_Std),
    names_sep = "_"
  ) %>%
  mutate(
    # Significant in T1D but not LC
    T1D_specific = (p_value_Type_1_Diabetes < 0.05) & 
      (is.na(p_value_Lean_Control) | p_value_Lean_Control >= 0.05),
    # Significant in T2D but not LC
    T2D_specific = (p_value_Type_2_Diabetes < 0.05) & 
      (is.na(p_value_Lean_Control) | p_value_Lean_Control >= 0.05),
    # Significant in both T1D and T2D
    Both_diabetes = (p_value_Type_1_Diabetes < 0.05) & 
      (p_value_Type_2_Diabetes < 0.05),
    # Different direction of effect (using standardized betas)
    Opposite_direction = sign(Beta_Std_Type_1_Diabetes) != sign(Beta_Std_Lean_Control) |
      sign(Beta_Std_Type_2_Diabetes) != sign(Beta_Std_Lean_Control),
    # Large difference in effect size (>50% difference in standardized betas)
    Large_diff_T1D_vs_LC = abs(Beta_Std_Type_1_Diabetes - 
                                 coalesce(Beta_Std_Lean_Control, 0)) > 0.5 * abs(Beta_Std_Type_1_Diabetes),
    Large_diff_T2D_vs_LC = abs(Beta_Std_Type_2_Diabetes - 
                                 coalesce(Beta_Std_Lean_Control, 0)) > 0.5 * abs(Beta_Std_Type_2_Diabetes)
  )

# Identify disease-specific associations
t1d_specific_assoc <- disease_specific %>%
  filter(T1D_specific) %>%
  arrange(p_value_Type_1_Diabetes)

t2d_specific_assoc <- disease_specific %>%
  filter(T2D_specific) %>%
  arrange(p_value_Type_2_Diabetes)

both_diabetes_assoc <- disease_specific %>%
  filter(Both_diabetes) %>%
  arrange(pmin(p_value_Type_1_Diabetes, p_value_Type_2_Diabetes))

cat("Disease-specific associations found:\n")
cat("  T1D-specific (sig in T1D, not LC):", nrow(t1d_specific_assoc), "\n")
cat("  T2D-specific (sig in T2D, not LC):", nrow(t2d_specific_assoc), "\n")
cat("  Both diabetes types:", nrow(both_diabetes_assoc), "\n\n")

write.csv(t1d_specific_assoc, "T1D_specific_associations.csv", row.names = FALSE)
write.csv(t2d_specific_assoc, "T2D_specific_associations.csv", row.names = FALSE)
write.csv(both_diabetes_assoc, "both_diabetes_associations.csv", row.names = FALSE)

# ============================================================
# CREATE FOREST PLOTS FOR TOP ASSOCIATIONS
# ============================================================

cat("=== Creating forest plots ===\n")

# Function to create forest plot
create_forest_plot <- function(predictor, outcome, results_data, title_suffix = "") {
  
  plot_data <- results_data %>%
    filter(Predictor == predictor, Outcome == outcome, Model_Type == "Adjusted") %>%
    mutate(
      CI_lower = Beta_Std - 1.96 * SE_Std,
      CI_upper = Beta_Std + 1.96 * SE_Std,
      Significant = p_value < 0.05,
      Group_Label = case_when(
        Group == "All_Participants" ~ "All Participants",
        Group == "Type_1_Diabetes" ~ "Type 1 Diabetes",
        Group == "Type_2_Diabetes" ~ "Type 2 Diabetes",
        Group == "Lean_Control" ~ "Lean Control",
        TRUE ~ Group
      )
    ) %>%
    arrange(Beta_Std)
  
  if(nrow(plot_data) == 0) return(NULL)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Beta_Std, y = reorder(Group_Label, Beta_Std))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Significant),
                   height = 0.2, size = 1) +
    geom_point(aes(color = Significant, size = Significant), shape = 18) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3), guide = "none") +
    labs(
      title = paste0(predictor, " → ", outcome, title_suffix),
      subtitle = "Standardized beta coefficients with 95% CI (adjusted for age, sex, BMI)",
      x = "Standardized Beta Coefficient",
      y = NULL,
      color = "Significance"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10)
    )
  
  # Add N and p-value annotations
  p <- p + geom_text(aes(label = sprintf("N=%d\np=%.3f", N, p_value)),
                     x = max(plot_data$CI_upper) * 1.1,
                     size = 3, hjust = 0)
  
  return(p)
}

# Create forest plots for top significant associations
top_associations <- significant_results %>%
  filter(Model_Type == "Adjusted") %>%
  group_by(Predictor, Outcome) %>%
  summarise(
    Min_p = min(p_value),
    N_Groups = n(),
    .groups = "drop"
  ) %>%
  filter(N_Groups >= 2) %>%  # Must appear in at least 2 groups
  arrange(Min_p) %>%
  head(15)

if(nrow(top_associations) > 0) {
  
  forest_plots <- list()
  
  for(i in 1:nrow(top_associations)) {
    predictor <- top_associations$Predictor[i]
    outcome <- top_associations$Outcome[i]
    
    p <- create_forest_plot(predictor, outcome, all_results)
    
    if(!is.null(p)) {
      forest_plots[[i]] <- p
      
      # Save individual plot
      filename <- paste0("forest_plot_", 
                         gsub("[^A-Za-z0-9]", "_", predictor), "_",
                         gsub("[^A-Za-z0-9]", "_", outcome), ".png")
      ggsave(filename, p, width = 10, height = 6, dpi = 300, bg = "white")
    }
  }
  
  cat("Created", length(forest_plots), "forest plots\n")
}

# ============================================================
# CREATE COMBINED FOREST PLOT FOR DISEASE-SPECIFIC ASSOCIATIONS
# ============================================================

cat("\n=== Creating disease-specific comparison plots ===\n")

# T1D-specific top 6
if(nrow(t1d_specific_assoc) > 0) {
  
  t1d_plots <- list()
  
  for(i in 1:min(6, nrow(t1d_specific_assoc))) {
    predictor <- t1d_specific_assoc$Predictor[i]
    outcome <- t1d_specific_assoc$Outcome[i]
    
    p <- create_forest_plot(predictor, outcome, all_results, 
                            title_suffix = " (T1D-specific)")
    
    if(!is.null(p)) {
      t1d_plots[[i]] <- p
    }
  }
  
  if(length(t1d_plots) > 0) {
    combined_t1d <- wrap_plots(t1d_plots, ncol = 2)
    ggsave("forest_plots_T1D_specific.png", combined_t1d,
           width = 16, height = 12, dpi = 300, bg = "white")
    cat("Created: forest_plots_T1D_specific.png\n")
  }
}

# T2D-specific top 6
if(nrow(t2d_specific_assoc) > 0) {
  
  t2d_plots <- list()
  
  for(i in 1:min(6, nrow(t2d_specific_assoc))) {
    predictor <- t2d_specific_assoc$Predictor[i]
    outcome <- t2d_specific_assoc$Outcome[i]
    
    p <- create_forest_plot(predictor, outcome, all_results, 
                            title_suffix = " (T2D-specific)")
    
    if(!is.null(p)) {
      t2d_plots[[i]] <- p
    }
  }
  
  if(length(t2d_plots) > 0) {
    combined_t2d <- wrap_plots(t2d_plots, ncol = 2)
    ggsave("forest_plots_T2D_specific.png", combined_t2d,
           width = 16, height = 12, dpi = 300, bg = "white")
    cat("Created: forest_plots_T2D_specific.png\n")
  }
}

# ============================================================
# EFFECT SIZE COMPARISON HEATMAP ACROSS GROUPS
# ============================================================

cat("\n=== Creating effect size comparison heatmap ===\n")

# Get associations that are significant in at least one group
sig_in_any <- adjusted_results %>%
  group_by(Predictor, Outcome) %>%
  filter(any(p_value < 0.05)) %>%
  ungroup()

# Create wide format for heatmap
beta_wide <- sig_in_any %>%
  mutate(Group_Short = case_when(
    Group == "All_Participants" ~ "All",
    Group == "Type_1_Diabetes" ~ "T1D",
    Group == "Type_2_Diabetes" ~ "T2D",
    Group == "Lean_Control" ~ "LC",
    TRUE ~ Group
  )) %>%
  mutate(Combo = paste(Predictor, Outcome, sep = " → ")) %>%
  dplyr::select(Combo, Group_Short, Beta_Std, p_value) %>%
  pivot_wider(names_from = Group_Short, values_from = c(Beta_Std, p_value))

# Create matrix of betas
beta_cols <- grep("^Beta_Std_", names(beta_wide), value = TRUE)
beta_matrix <- as.matrix(beta_wide[, beta_cols])
rownames(beta_matrix) <- beta_wide$Combo
colnames(beta_matrix) <- gsub("Beta_Std_", "", colnames(beta_matrix))

# Create significance markers
pval_cols <- grep("^p_value_", names(beta_wide), value = TRUE)
pval_matrix <- as.matrix(beta_wide[, pval_cols])
colnames(pval_matrix) <- gsub("p_value_", "", colnames(pval_matrix))

sig_markers <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
sig_markers[pval_matrix < 0.001] <- "***"
sig_markers[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
sig_markers[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"
sig_markers[is.na(pval_matrix)] <- ""

# Limit to top 30 associations by variance in beta across groups
beta_variance <- apply(beta_matrix, 1, var, na.rm = TRUE)
top_30_idx <- order(beta_variance, decreasing = TRUE)[1:min(30, length(beta_variance))]

beta_matrix_top30 <- beta_matrix[top_30_idx, ]
sig_markers_top30 <- sig_markers[top_30_idx, ]

# Create combined display matrix
display_matrix <- matrix(
  paste0(sprintf("%.2f", beta_matrix_top30), sig_markers_top30),
  nrow = nrow(beta_matrix_top30)
)

pheatmap(beta_matrix_top30,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Effect Size Comparison Across Groups\n(Top 30 by variance in standardized effect)",
         fontsize_row = 7,
         fontsize_col = 10,
         angle_col = 0,
         display_numbers = display_matrix,
         number_color = "black",
         fontsize_number = 6,
         breaks = seq(-max(abs(beta_matrix_top30), na.rm = TRUE), 
                      max(abs(beta_matrix_top30), na.rm = TRUE), 
                      length.out = 101),
         filename = "effect_size_comparison_heatmap.png",
         width = 8,
         height = 14)

pheatmap(beta_matrix_top30,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Effect Size Comparison Across Groups\n(Top 30 by variance in standardized effect)",
         fontsize_row = 7,
         fontsize_col = 10,
         angle_col = 0,
         display_numbers = display_matrix,
         number_color = "black",
         fontsize_number = 6,
         breaks = seq(-max(abs(beta_matrix_top30), na.rm = TRUE), 
                      max(abs(beta_matrix_top30), na.rm = TRUE), 
                      length.out = 101),
         filename = "effect_size_comparison_heatmap.pdf",
         width = 8,
         height = 14)

cat("Created: effect_size_comparison_heatmap.png/pdf\n")

# ============================================================
# DOT PLOT SHOWING EFFECT SIZES ACROSS GROUPS
# ============================================================

cat("\n=== Creating dot plot comparison ===\n")

# Top 20 predictor-outcome pairs by overall significance
top_20_pairs <- adjusted_results %>%
  group_by(Predictor, Outcome) %>%
  summarise(
    Min_p = min(p_value, na.rm = TRUE),
    Max_abs_beta = max(abs(Beta), na.rm = TRUE),
    N_Groups = n(),
    .groups = "drop"
  ) %>%
  filter(N_Groups >= 3) %>%  # Must be in at least 3 groups
  arrange(Min_p) %>%
  head(20)

if(nrow(top_20_pairs) > 0) {
  
  plot_data_dots <- adjusted_results %>%
    inner_join(top_20_pairs %>% dplyr::select(Predictor, Outcome), 
               by = c("Predictor", "Outcome")) %>%
    mutate(
      Combo = paste(Predictor, "→", Outcome),
      Group_Label = case_when(
        Group == "All_Participants" ~ "All",
        Group == "Type_1_Diabetes" ~ "T1D",
        Group == "Type_2_Diabetes" ~ "T2D",
        Group == "Lean_Control" ~ "LC",
        TRUE ~ Group
      ),
      Significant = p_value < 0.05
    )
  
  p_dots <- ggplot(plot_data_dots, 
                   aes(x = Beta_Std, y = reorder(Combo, Beta_Std), color = Group_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = ifelse(Significant, 3, 2), 
                   shape = ifelse(Significant, 16, 1)),
               position = position_dodge(width = 0.5), alpha = 0.8) +
    scale_color_manual(values = c("All" = "black", "T1D" = "#E41A1C", 
                                  "T2D" = "#377EB8", "LC" = "#4DAF4A")) +
    scale_size_identity() +
    scale_shape_identity() +
    labs(
      title = "Effect Size Comparison Across Groups (Top 20 Associations)",
      subtitle = "Filled circles = p < 0.05, Open circles = p ≥ 0.05 | Standardized beta coefficients",
      x = "Standardized Beta Coefficient (adjusted for age, sex, BMI)",
      y = NULL,
      color = "Group"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold")
    )
  
  ggsave("effect_size_dot_plot_comparison.png", p_dots,
         width = 12, height = 10, dpi = 300, bg = "white")
  ggsave("effect_size_dot_plot_comparison.pdf", p_dots,
         width = 12, height = 10)
  
  cat("Created: effect_size_dot_plot_comparison.png/pdf\n")
}

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n\n##########################################################")
cat("\n### ANALYSIS COMPLETE")
cat("\n##########################################################\n\n")

cat("Results saved:\n")
cat("  - all_regression_results.csv (all models)\n")
cat("  - significant_associations.csv (p < 0.05 only)\n")
cat("  - summary_by_predictor_category.csv\n")
cat("  - T1D_specific_associations.csv\n")
cat("  - T2D_specific_associations.csv\n")
cat("  - both_diabetes_associations.csv\n")
cat("  - Heatmaps: heatmap_pvalues_*.png and heatmap_beta_*.png\n")
cat("  - Forest plots: forest_plot_*.png and forest_plots_*_specific.png\n")
cat("  - Comparison plots: effect_size_comparison_heatmap.png and effect_size_dot_plot_comparison.png\n")
cat("  - Scatter plots: top_significant_associations_scatter.png\n\n")

cat("Summary statistics:\n")
cat("  Total regressions:", nrow(all_results), "\n")
cat("  Significant (p < 0.05):", sum(all_results$p_value < 0.05), "\n")
cat("  Highly significant (p < 0.001):", sum(all_results$p_value < 0.001), "\n")
cat("  T1D-specific associations:", nrow(t1d_specific_assoc), "\n")
cat("  T2D-specific associations:", nrow(t2d_specific_assoc), "\n")
cat("  Both diabetes types:", nrow(both_diabetes_assoc), "\n")
cat("  Sample type used:", primary_sample_type, "\n\n")

# Save workspace
save.image("clinical_biomarker_analysis.RData")
cat("Workspace saved: clinical_biomarker_analysis.RData\n")











