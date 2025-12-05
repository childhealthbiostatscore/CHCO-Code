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
# Function to create p-value heatmap with consistent ordering
create_pvalue_heatmap <- function(results_subset, title, filename, predictor_order = NULL) {
  
  # Pivot to matrix format
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  # Apply consistent row ordering if provided
  if(!is.null(predictor_order)) {
    # Keep only predictors that exist in this subset
    predictor_order <- predictor_order[predictor_order %in% rownames(pval_matrix)]
    if(length(predictor_order) > 0) {
      pval_matrix <- pval_matrix[predictor_order, , drop = FALSE]
    }
  }
  
  # Transform p-values for visualization (-log10)
  pval_transformed <- -log10(pval_matrix)
  
  # Cap at -log10(0.001) for visualization
  pval_transformed[pval_transformed > 3] <- 3
  
  # Only create heatmap if we have data
  if(nrow(pval_transformed) > 0 && ncol(pval_transformed) > 0) {
    
    pheatmap(pval_transformed,
             cluster_rows = FALSE,  # Don't cluster - keep manual order
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

# ============================================================
# DEFINE CONSISTENT PREDICTOR ORDER FOR ALL HEATMAPS
# ============================================================

# Create a consistent ordering based on predictor categories
predictor_order_all <- all_results %>%
  mutate(
    Predictor_Category = case_when(
      Predictor %in% clinical_predictors$`Glycemic Control` ~ "1_Glycemic",
      Predictor %in% clinical_predictors$`Insulin Sensitivity` ~ "2_Insulin",
      Predictor %in% clinical_predictors$`CGM Metrics` ~ "3_CGM",
      Predictor %in% clinical_predictors$`Blood Pressure` ~ "4_BP",
      TRUE ~ "5_Other"
    )
  ) %>%
  arrange(Predictor_Category, Predictor_Used) %>%
  pull(Predictor_Used) %>%
  unique()

cat("Predictor order established for consistent heatmaps:\n")
print(predictor_order_all)
cat("\n")


# Create heatmaps for each group and model type
for(group_name in unique(all_results$Group)) {
  for(model_type in c("Unadjusted", "Adjusted")) {
    
    results_subset <- all_results %>%
      filter(Group == group_name, Model_Type == model_type)
    
    if(nrow(results_subset) > 0) {
      title <- paste0(group_name, " - ", model_type, " Models\n(-log10 p-values)")
      filename <- paste0("heatmap_pvalues_", gsub(" ", "_", group_name), "_", model_type, ".png")
      
      success <- create_pvalue_heatmap(results_subset, title, filename, predictor_order = predictor_order_all)
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
create_beta_heatmap <- function(results_subset, title, filename, predictor_order = NULL) {
  
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
  
  # Apply consistent row ordering if provided
  if(!is.null(predictor_order)) {
    # Keep only predictors that exist in this subset
    predictor_order <- predictor_order[predictor_order %in% rownames(beta_matrix)]
    if(length(predictor_order) > 0) {
      beta_matrix <- beta_matrix[predictor_order, , drop = FALSE]
      pval_matrix <- pval_matrix[predictor_order, , drop = FALSE]
    }
  }
  
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
             cluster_rows = FALSE,  # Don't cluster - keep manual order
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
      
      success <- create_beta_heatmap(results_subset, title, filename, predictor_order = predictor_order_all)
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
    Mean_Abs_Beta_Std = mean(abs(Beta_Std)),
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
disease_specific_prep <- adjusted_results %>%
  filter(Group %in% c("Type_1_Diabetes", "Type_2_Diabetes", "Lean_Control")) %>%
  dplyr::select(Group, Predictor, Outcome, Beta_Std, p_value, SE_Std)

# Check if we have data
if(nrow(disease_specific_prep) == 0) {
  cat("No data available for disease-specific comparison\n")
  t1d_specific_assoc <- data.frame()
  t2d_specific_assoc <- data.frame()
  both_diabetes_assoc <- data.frame()
} else {
  
  # Check which groups are actually present
  groups_present <- unique(disease_specific_prep$Group)
  cat("\nGroups present in adjusted results:\n")
  print(groups_present)
  cat("\n")
  
  disease_specific <- disease_specific_prep %>%
    pivot_wider(
      names_from = Group,
      values_from = c(Beta_Std, p_value, SE_Std)
    )
  
  # Print column names to debug
  cat("Column names after pivot_wider:\n")
  print(names(disease_specific))
  cat("\n")
  
  # Get actual column names dynamically
  beta_t1d_col <- grep("Beta_Std.*Type.*1.*Diabetes", names(disease_specific), value = TRUE)[1]
  beta_t2d_col <- grep("Beta_Std.*Type.*2.*Diabetes", names(disease_specific), value = TRUE)[1]
  beta_lc_col <- grep("Beta_Std.*Lean.*Control", names(disease_specific), value = TRUE)[1]
  
  pval_t1d_col <- grep("p_value.*Type.*1.*Diabetes", names(disease_specific), value = TRUE)[1]
  pval_t2d_col <- grep("p_value.*Type.*2.*Diabetes", names(disease_specific), value = TRUE)[1]
  pval_lc_col <- grep("p_value.*Lean.*Control", names(disease_specific), value = TRUE)[1]
  
  # Check which columns exist
  has_t1d <- !is.na(beta_t1d_col) && !is.na(pval_t1d_col)
  has_t2d <- !is.na(beta_t2d_col) && !is.na(pval_t2d_col)
  has_lc <- !is.na(beta_lc_col) && !is.na(pval_lc_col)
  
  cat("Groups available:\n")
  cat("  T1D:", has_t1d, "\n")
  cat("  T2D:", has_t2d, "\n")
  cat("  Lean Control:", has_lc, "\n\n")
  
  # Only proceed if we have at least one diabetes group and lean control
  if(has_lc && (has_t1d || has_t2d)) {
    
    # Build mutate expression based on available columns
    if(has_t1d && has_lc) {
      disease_specific <- disease_specific %>%
        mutate(
          T1D_specific = (.data[[pval_t1d_col]] < 0.05) & 
            (is.na(.data[[pval_lc_col]]) | .data[[pval_lc_col]] >= 0.05),
          Large_diff_T1D_vs_LC = abs(.data[[beta_t1d_col]] - 
                                       coalesce(.data[[beta_lc_col]], 0)) > 0.5 * abs(.data[[beta_t1d_col]])
        )
    }
    
    if(has_t2d && has_lc) {
      disease_specific <- disease_specific %>%
        mutate(
          T2D_specific = (.data[[pval_t2d_col]] < 0.05) & 
            (is.na(.data[[pval_lc_col]]) | .data[[pval_lc_col]] >= 0.05),
          Large_diff_T2D_vs_LC = abs(.data[[beta_t2d_col]] - 
                                       coalesce(.data[[beta_lc_col]], 0)) > 0.5 * abs(.data[[beta_t2d_col]])
        )
    }
    
    if(has_t1d && has_t2d) {
      disease_specific <- disease_specific %>%
        mutate(
          Both_diabetes = (.data[[pval_t1d_col]] < 0.05) & 
            (.data[[pval_t2d_col]] < 0.05)
        )
    }
    
    if(has_t1d && has_lc) {
      disease_specific <- disease_specific %>%
        mutate(
          Opposite_direction_T1D = sign(.data[[beta_t1d_col]]) != sign(.data[[beta_lc_col]])
        )
    }
    
    if(has_t2d && has_lc) {
      disease_specific <- disease_specific %>%
        mutate(
          Opposite_direction_T2D = sign(.data[[beta_t2d_col]]) != sign(.data[[beta_lc_col]])
        )
    }
    
    # Identify disease-specific associations
    if(has_t1d && "T1D_specific" %in% names(disease_specific)) {
      t1d_specific_assoc <- disease_specific %>%
        filter(T1D_specific) %>%
        arrange(.data[[pval_t1d_col]])
    } else {
      t1d_specific_assoc <- data.frame()
    }
    
    if(has_t2d && "T2D_specific" %in% names(disease_specific)) {
      t2d_specific_assoc <- disease_specific %>%
        filter(T2D_specific) %>%
        arrange(.data[[pval_t2d_col]])
    } else {
      t2d_specific_assoc <- data.frame()
    }
    
    if(has_t1d && has_t2d && "Both_diabetes" %in% names(disease_specific)) {
      both_diabetes_assoc <- disease_specific %>%
        filter(Both_diabetes) %>%
        arrange(pmin(.data[[pval_t1d_col]], .data[[pval_t2d_col]]))
    } else {
      both_diabetes_assoc <- data.frame()
    }
    
  } else {
    cat("Insufficient groups for comparison (need at least one diabetes group + Lean Control)\n")
    t1d_specific_assoc <- data.frame()
    t2d_specific_assoc <- data.frame()
    both_diabetes_assoc <- data.frame()
  }
}

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
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Significant),
                   height = 0.2, linewidth = 1) +
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
    Max_abs_beta = max(abs(Beta_Std), na.rm = TRUE),
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
  
  # ============================================================
  # CREATE SAME PLOT FOR UNADJUSTED MODELS
  # ============================================================
  
  cat("\n=== Creating unadjusted effect size comparison ===\n")
  
  # Get unadjusted results
  unadjusted_results <- all_results %>%
    filter(Model_Type == "Unadjusted")
  
  # Top 20 predictor-outcome pairs by overall significance (unadjusted)
  top_20_pairs_unadj <- unadjusted_results %>%
    group_by(Predictor, Outcome) %>%
    summarise(
      Min_p = min(p_value, na.rm = TRUE),
      Max_abs_beta = max(abs(Beta_Std), na.rm = TRUE),
      N_Groups = n(),
      .groups = "drop"
    ) %>%
    filter(N_Groups >= 3) %>%  # Must be in at least 3 groups
    arrange(Min_p) %>%
    head(20)
  
  if(nrow(top_20_pairs_unadj) > 0) {
    
    plot_data_dots_unadj <- unadjusted_results %>%
      inner_join(top_20_pairs_unadj %>% dplyr::select(Predictor, Outcome), 
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
    
    p_dots_unadj <- ggplot(plot_data_dots_unadj, 
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
        title = "Effect Size Comparison Across Groups (Top 20 Associations - UNADJUSTED)",
        subtitle = "Filled circles = p < 0.05, Open circles = p ≥ 0.05 | Standardized beta coefficients (no covariates)",
        x = "Standardized Beta Coefficient (unadjusted)",
        y = NULL,
        color = "Group"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold")
      )
    
    ggsave("effect_size_dot_plot_comparison_unadjusted.png", p_dots_unadj,
           width = 12, height = 10, dpi = 300, bg = "white")
    ggsave("effect_size_dot_plot_comparison_unadjusted.pdf", p_dots_unadj,
           width = 12, height = 10)
    
    cat("Created: effect_size_dot_plot_comparison_unadjusted.png/pdf\n")
  } else {
    cat("Not enough unadjusted results to create comparison plot\n")
  }
}

cat("\n\nDot plot comparison complete.\n")
# ============================================================
# INTERACTION ANALYSES - TEST IF EFFECTS DIFFER BY GROUP
# ============================================================

cat("\n\n##########################################################")
cat("\n### INTERACTION ANALYSES")
cat("\n### Testing if predictor effects differ between groups")
cat("\n##########################################################\n\n")

# Function to test interactions between two specific groups
test_interaction_two_groups <- function(data, predictor, outcome, group1, group2, covariates = c("age", "sex", "bmi")) {
  
  # Filter to only the two groups of interest
  analysis_data <- data %>%
    filter(group %in% c(group1, group2)) %>%
    dplyr::select(all_of(c(predictor, outcome, "group", covariates))) %>%
    filter(complete.cases(.))
  
  if(nrow(analysis_data) < 20) return(NULL)
  
  # Make group a factor with group1 as reference
  analysis_data$group <- factor(analysis_data$group, levels = c(group1, group2))
  
  # Check if we need log transformation
  pred_needs_log <- test_log_transformation(analysis_data[[predictor]])
  outcome_needs_log <- test_log_transformation(analysis_data[[outcome]])
  
  pred_var_name <- predictor
  outcome_var_name <- outcome
  
  if(pred_needs_log) {
    analysis_data[[predictor]] <- log(analysis_data[[predictor]])
    pred_var_name <- paste0("log(", predictor, ")")
  }
  
  if(outcome_needs_log) {
    analysis_data[[outcome]] <- log(analysis_data[[outcome]])
    outcome_var_name <- paste0("log(", outcome, ")")
  }
  
  # Standardize continuous variables
  analysis_data[[predictor]] <- scale(analysis_data[[predictor]])[,1]
  analysis_data[[outcome]] <- scale(analysis_data[[outcome]])[,1]
  
  for(cov in covariates) {
    if(cov != "sex" && is.numeric(analysis_data[[cov]])) {
      analysis_data[[cov]] <- scale(analysis_data[[cov]])[,1]
    }
  }
  
  # Main effects model (no interaction)
  formula_main <- as.formula(paste0(outcome, " ~ ", predictor, " + group + ", 
                                    paste(covariates, collapse = " + ")))
  model_main <- lm(formula_main, data = analysis_data)
  
  # Interaction model
  formula_int <- as.formula(paste0(outcome, " ~ ", predictor, " * group + ", 
                                   paste(covariates, collapse = " + ")))
  model_int <- lm(formula_int, data = analysis_data)
  
  # Test if interaction improves model fit
  anova_result <- anova(model_main, model_int)
  interaction_p <- anova_result$`Pr(>F)`[2]
  
  # Extract coefficients
  coef_summary <- summary(model_int)$coefficients
  
  # Get main predictor effect (effect in reference group = group1)
  main_effect_row <- which(rownames(coef_summary) == predictor)
  if(length(main_effect_row) == 0) return(NULL)
  
  main_effect <- coef_summary[main_effect_row, "Estimate"]
  main_effect_p <- coef_summary[main_effect_row, "Pr(>|t|)"]
  
  # Get interaction term (difference in slope between groups)
  interaction_term <- paste0(predictor, ":group", group2)
  int_row <- which(rownames(coef_summary) == interaction_term)
  
  if(length(int_row) == 0) return(NULL)
  
  interaction_beta <- coef_summary[int_row, "Estimate"]
  interaction_se <- coef_summary[int_row, "Std. Error"]
  interaction_t <- coef_summary[int_row, "t value"]
  interaction_p_coef <- coef_summary[int_row, "Pr(>|t|)"]
  
  # Effect in group2 is main_effect + interaction_beta
  group2_effect <- main_effect + interaction_beta
  
  results <- data.frame(
    Predictor = predictor,
    Outcome = outcome,
    Predictor_Used = pred_var_name,
    Outcome_Used = outcome_var_name,
    Predictor_Log_Transformed = pred_needs_log,
    Outcome_Log_Transformed = outcome_needs_log,
    Reference_Group = group1,
    Comparison_Group = group2,
    N_Total = nrow(analysis_data),
    N_Group1 = sum(analysis_data$group == group1),
    N_Group2 = sum(analysis_data$group == group2),
    Effect_in_Group1 = main_effect,
    Effect_in_Group1_p = main_effect_p,
    Effect_in_Group2 = group2_effect,
    Interaction_Beta_Std = interaction_beta,
    Interaction_SE = interaction_se,
    Interaction_t = interaction_t,
    Interaction_p = interaction_p_coef,
    Model_Interaction_p = interaction_p,
    R_squared_main = summary(model_main)$r.squared,
    R_squared_interaction = summary(model_int)$r.squared,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

# Check sample types for each group
cat("Checking sample type distribution by group:\n")
sample_type_check <- dat_analysis %>%
  group_by(group, sample_type) %>%
  summarise(N = n(), .groups = "drop")
print(sample_type_check)
cat("\n")

# Identify which comparisons are valid (same sample type)
group_sample_types <- dat_analysis %>%
  group_by(group) %>%
  summarise(sample_type = unique(sample_type)[1], .groups = "drop")

cat("Group sample types:\n")
print(group_sample_types)
cat("\n")

# Define valid comparisons (only compare groups with same sample type)
valid_comparisons <- list()

# Check LC sample type
lc_sample_type <- group_sample_types %>% 
  filter(group == "Lean Control") %>% 
  pull(sample_type)

# Check T2D sample type
t2d_exists <- "Type 2 Diabetes" %in% group_sample_types$group
if(t2d_exists) {
  t2d_sample_type <- group_sample_types %>% 
    filter(group == "Type 2 Diabetes") %>% 
    pull(sample_type)
  
  if(length(lc_sample_type) > 0 && length(t2d_sample_type) > 0) {
    if(lc_sample_type == t2d_sample_type) {
      valid_comparisons[["LC_vs_T2D"]] <- list(
        group1 = "Lean Control",
        group2 = "Type 2 Diabetes",
        sample_type = lc_sample_type
      )
      cat("✓ LC vs T2D comparison valid (both ", lc_sample_type, ")\n", sep = "")
    } else {
      cat("✗ LC vs T2D comparison INVALID (LC=", lc_sample_type, ", T2D=", t2d_sample_type, ")\n", sep = "")
    }
  }
}

# Check T1D sample type
t1d_exists <- "Type 1 Diabetes" %in% group_sample_types$group
if(t1d_exists) {
  t1d_sample_type <- group_sample_types %>% 
    filter(group == "Type 1 Diabetes") %>% 
    pull(sample_type)
  
  if(length(lc_sample_type) > 0 && length(t1d_sample_type) > 0) {
    if(lc_sample_type == t1d_sample_type) {
      valid_comparisons[["LC_vs_T1D"]] <- list(
        group1 = "Lean Control",
        group2 = "Type 1 Diabetes",
        sample_type = lc_sample_type
      )
      cat("✓ LC vs T1D comparison valid (both ", lc_sample_type, ")\n", sep = "")
    } else {
      cat("✗ LC vs T1D comparison INVALID (LC=", lc_sample_type, ", T1D=", t1d_sample_type, ")\n", sep = "")
    }
  }
}

cat("\n")

if(length(valid_comparisons) == 0) {
  cat("WARNING: No valid group comparisons found (sample type mismatch)\n\n")
}

# Get ALL predictor-outcome combinations that were tested
all_interaction_results <- list()

for(comparison_name in names(valid_comparisons)) {
  
  comparison <- valid_comparisons[[comparison_name]]
  
  cat("\n=================================================\n")
  cat("Testing interactions for:", comparison_name, "\n")
  cat("=================================================\n\n")
  
  # Get predictor-outcome combinations tested in both groups
  predictors_to_test <- all_results %>%
    filter(Model_Type == "Adjusted") %>%
    filter(Group %in% c(comparison$group1, comparison$group2)) %>%
    group_by(Predictor, Outcome) %>%
    filter(n() == 2) %>%  # Must be tested in both groups
    ungroup() %>%
    distinct(Predictor, Outcome)
  
  cat("Testing", nrow(predictors_to_test), "predictor-outcome pairs\n\n")
  
  interaction_results <- data.frame()
  
  if(nrow(predictors_to_test) > 0) {
    
    for(i in 1:nrow(predictors_to_test)) {
      
      predictor <- predictors_to_test$Predictor[i]
      outcome <- predictors_to_test$Outcome[i]
      
      # Test interaction
      int_result <- test_interaction_two_groups(
        data = dat_analysis,
        predictor = predictor,
        outcome = outcome,
        group1 = comparison$group1,
        group2 = comparison$group2,
        covariates = c("age", "sex", "bmi")
      )
      
      if(!is.null(int_result)) {
        # Add comparison name and sample type
        int_result$Comparison <- comparison_name
        int_result$Sample_Type <- comparison$sample_type
        interaction_results <- rbind(interaction_results, int_result)
      }
    }
  }
  
  cat("Completed", nrow(interaction_results), "interaction tests for", comparison_name, "\n")
  
  # Save results for this comparison
  all_interaction_results[[comparison_name]] <- interaction_results
}

# Combine all interaction results
interaction_results_all <- bind_rows(all_interaction_results)

cat("\n\nTotal interaction tests completed:", nrow(interaction_results_all), "\n")

# ============================================================
# VISUALIZE SIGNIFICANT INTERACTIONS
# ============================================================

cat("\n=== Creating interaction plots ===\n")

if(nrow(significant_interactions) > 0) {
  
  # Create interaction plots for top 12 significant interactions
  top_12_interactions <- significant_interactions %>%
    head(12)
  
  interaction_plots <- list()
  
  for(i in 1:nrow(top_12_interactions)) {
    
    row <- top_12_interactions[i, ]
    
    # Get the data
    plot_data <- dat_analysis %>%
      dplyr::select(predictor = all_of(row$Predictor), 
                    outcome = all_of(row$Outcome),
                    group) %>%
      filter(complete.cases(.)) %>%
      filter(group %in% c(row$Reference_Group, row$Comparison_Group))
    
    if(nrow(plot_data) < 20) next
    
    # Apply transformations if needed
    x_label <- row$Predictor
    y_label <- row$Outcome
    
    if(row$Predictor_Log_Transformed) {
      plot_data$predictor <- log(plot_data$predictor)
      x_label <- paste0("log(", x_label, ")")
    }
    
    if(row$Outcome_Log_Transformed) {
      plot_data$outcome <- log(plot_data$outcome)
      y_label <- paste0("log(", y_label, ")")
    }
    
    # Create plot with separate regression lines by group
    p <- ggplot(plot_data, aes(x = predictor, y = outcome, color = group)) +
      geom_point(alpha = 0.5, size = 2) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
      scale_color_manual(values = c("Lean Control" = "#4DAF4A", 
                                    "Type 1 Diabetes" = "#E41A1C",
                                    "Type 2 Diabetes" = "#377EB8")) +
      labs(
        title = paste0(row$Predictor, " → ", row$Outcome),
        subtitle = sprintf("Interaction p = %.4f | Groups differ significantly in slope", 
                           row$Interaction_p),
        x = x_label,
        y = y_label,
        color = "Group"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9)
      )
    
    interaction_plots[[i]] <- p
    
    # Save individual plot
    filename <- paste0("interaction_plot_", 
                       gsub("[^A-Za-z0-9]", "_", row$Predictor), "_",
                       gsub("[^A-Za-z0-9]", "_", row$Outcome), ".png")
    ggsave(filename, p, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  # Combine plots
  if(length(interaction_plots) > 0) {
    combined_interactions <- wrap_plots(interaction_plots, ncol = 3)
    ggsave("significant_interactions_combined.png", combined_interactions,
           width = 18, height = 12, dpi = 300, bg = "white")
    cat("Created: significant_interactions_combined.png\n")
  }
}

# ============================================================
# INTERACTION HEATMAP
# ============================================================

cat("\n=== Creating interaction heatmap ===\n")

if(nrow(interaction_results_all) > 0) {
  
  # Create matrix of interaction p-values
  interaction_matrix <- interaction_results_all %>%
    mutate(Combo = paste(Reference_Group, "vs", Comparison_Group)) %>%
    dplyr::select(Predictor, Outcome, Combo, Interaction_p) %>%
    pivot_wider(names_from = Outcome, values_from = Interaction_p)
  
  # Split by comparison to create separate matrices
  unique_combos <- unique(paste(interaction_results_all$Reference_Group, 
                                "vs", 
                                interaction_results_all$Comparison_Group))
  
  for(combo in unique_combos) {
    
    combo_data <- interaction_results_all %>%
      filter(paste(Reference_Group, "vs", Comparison_Group) == combo)
    
    if(nrow(combo_data) == 0) next
    
    # Create p-value matrix
    pval_wide <- combo_data %>%
      dplyr::select(Predictor, Outcome, Interaction_p) %>%
      pivot_wider(names_from = Outcome, values_from = Interaction_p) %>%
      column_to_rownames("Predictor") %>%
      as.matrix()
    
    # Transform for visualization
    pval_transformed <- -log10(pval_wide)
    pval_transformed[pval_transformed > 3] <- 3
    
    if(nrow(pval_transformed) > 0 && ncol(pval_transformed) > 0) {
      
      pheatmap(pval_transformed,
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
               breaks = seq(0, 3, length.out = 101),
               main = paste0("Interaction P-values: ", combo, "\n(-log10 scale)"),
               fontsize_row = 8,
               fontsize_col = 9,
               angle_col = 45,
               display_numbers = matrix(sprintf("%.3f", pval_wide), nrow = nrow(pval_wide)),
               number_color = "black",
               fontsize_number = 6,
               legend_breaks = c(0, 0.5, 1, 1.3, 2, 3),
               legend_labels = c("1.0", "0.32", "0.10", "0.05", "0.01", "0.001"),
               filename = paste0("interaction_heatmap_", gsub(" ", "_", combo), ".png"),
               width = 10,
               height = max(6, nrow(pval_transformed) * 0.3))
      
      cat("Created: interaction_heatmap_", gsub(" ", "_", combo), ".png\n", sep = "")
    }
  }
}

cat("\n=== Interaction analysis interpretation ===\n")
cat("A significant interaction (p < 0.05) means:\n")
cat("  - The effect of the predictor on the outcome differs between groups\n")
cat("  - The slopes of the regression lines are significantly different\n")
cat("  - This explains why pooled 'All' estimates can differ from individual groups\n\n")

# ============================================================
# FINAL SUMMARY WITH INTERACTIONS
# ============================================================

cat("\n\n##########################################################")
cat("\n### ANALYSIS COMPLETE")
cat("\n##########################################################\n\n")

cat("Results saved:\n")
cat("  - all_regression_results.csv (all models)\n")
cat("  - significant_associations.csv (p < 0.05 only)\n")
cat("  - interaction_analysis_results.csv (all interaction tests)\n")
cat("  - significant_interactions.csv (significant interactions only)\n")
cat("  - summary_by_predictor_category.csv\n")
cat("  - T1D_specific_associations.csv\n")
cat("  - T2D_specific_associations.csv\n")
cat("  - both_diabetes_associations.csv\n")
cat("  - Heatmaps: heatmap_pvalues_*.png and heatmap_beta_*.png\n")
cat("  - Forest plots: forest_plot_*.png and forest_plots_*_specific.png\n")
cat("  - Comparison plots: effect_size_comparison_heatmap.png and effect_size_dot_plot_comparison.png\n")
cat("  - Interaction plots: interaction_plot_*.png and significant_interactions_combined.png\n")
cat("  - Interaction heatmaps: interaction_heatmap_*.png\n")
cat("  - Scatter plots: top_significant_associations_scatter.png\n\n")

cat("Summary statistics:\n")
cat("  Total regressions:", nrow(all_results), "\n")
cat("  Significant (p < 0.05):", sum(all_results$p_value < 0.05), "\n")
cat("  Highly significant (p < 0.001):", sum(all_results$p_value < 0.001), "\n")
if(exists("interaction_results_all") && nrow(interaction_results_all) > 0) {
  cat("  Interaction tests performed:", nrow(interaction_results_all), "\n")
  cat("  Significant interactions (p < 0.05):", sum(interaction_results_all$Interaction_p < 0.05), "\n")
}
cat("  T1D-specific associations:", nrow(t1d_specific_assoc), "\n")
cat("  T2D-specific associations:", nrow(t2d_specific_assoc), "\n")
cat("  Both diabetes types:", nrow(both_diabetes_assoc), "\n")
cat("  Sample type used:", primary_sample_type, "\n\n")

# Save workspace
save.image("clinical_biomarker_analysis.RData")
cat("Workspace saved: clinical_biomarker_analysis.RData\n")




