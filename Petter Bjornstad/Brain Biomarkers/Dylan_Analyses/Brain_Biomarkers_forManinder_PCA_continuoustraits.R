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
# FUNCTION TO RUN REGRESSIONS
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
    
    # Unadjusted model
    formula_str <- paste0(outcome, " ~ ", predictor)
    model_type <- "Unadjusted"
    
  } else {
    analysis_data <- data %>% 
      dplyr::select(all_of(c(predictor, outcome, covariates))) %>%
      filter(complete.cases(.))
    
    if(nrow(analysis_data) < 10) return(NULL)
    
    # Adjusted model
    formula_str <- paste0(outcome, " ~ ", predictor, " + ", paste(covariates, collapse = " + "))
    model_type <- "Adjusted"
  }
  
  # Fit model
  model <- lm(as.formula(formula_str), data = analysis_data)
  
  # Extract predictor results (first coefficient after intercept)
  model_summary <- summary(model)
  coef_table <- coef(model_summary)
  
  # Get predictor coefficient (not intercept, not covariates)
  predictor_row <- which(rownames(coef_table) == predictor)
  
  if(length(predictor_row) == 0) return(NULL)
  
  results <- data.frame(
    Group = group_name,
    Sample_Type = unique(data$sample_type)[1],
    Predictor = predictor,
    Outcome = outcome,
    Model_Type = model_type,
    N = nrow(analysis_data),
    Beta = coef_table[predictor_row, "Estimate"],
    SE = coef_table[predictor_row, "Std. Error"],
    t_value = coef_table[predictor_row, "t value"],
    p_value = coef_table[predictor_row, "Pr(>|t|)"],
    R_squared = model_summary$r.squared,
    Adj_R_squared = model_summary$adj.r.squared,
    Model_p_value = pf(model_summary$fstatistic[1], 
                       model_summary$fstatistic[2], 
                       model_summary$fstatistic[3], 
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
    # Standardized beta (approximate - for comparison across predictors)
    # This is just Beta * SD(X) / SD(Y), but we'll flag it for interpretation
    Sig_Flag = p_value < 0.05
  )

cat("\n\nTotal regressions completed:", nrow(all_results), "\n")
cat("Significant associations (p < 0.05):", sum(all_results$Sig_Flag), "\n\n")

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
    dplyr::select(Group, Model_Type, Predictor, Outcome, N, Beta, p_value, R_squared)
  
  print(kable(top_20, digits = c(0, 0, 0, 0, 0, 4, 4, 3), format = "simple"))
}

# ============================================================
# HEATMAPS OF P-VALUES
# ============================================================

cat("\n\n=== Creating heatmaps ===\n")

# Function to create p-value heatmap
create_pvalue_heatmap <- function(results_subset, title, filename) {
  
  # Pivot to matrix format
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor") %>%
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
  
  # Pivot to matrix format
  beta_matrix <- results_subset %>%
    dplyr::select(Predictor, Outcome, Beta) %>%
    pivot_wider(names_from = Outcome, values_from = Beta) %>%
    column_to_rownames("Predictor") %>%
    as.matrix()
  
  # Get corresponding p-values for significance markers
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor") %>%
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
             main = title,
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
    
    # Create plot
    p <- ggplot(plot_data, aes(x = x, y = y)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "blue") +
      labs(
        title = paste0(row$Predictor, " → ", row$Outcome),
        subtitle = sprintf("%s | %s | β=%.3f, p=%s, R²=%.3f",
                           row$Group, row$Model_Type, row$Beta, 
                           format.pval(row$p_value, digits = 2), row$R_squared),
        x = row$Predictor,
        y = row$Outcome
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
# FINAL SUMMARY
# ============================================================

cat("\n\n##########################################################")
cat("\n### ANALYSIS COMPLETE")
cat("\n##########################################################\n\n")

cat("Results saved:\n")
cat("  - all_regression_results.csv (all models)\n")
cat("  - significant_associations.csv (p < 0.05 only)\n")
cat("  - summary_by_predictor_category.csv\n")
cat("  - Heatmaps: heatmap_pvalues_*.png and heatmap_beta_*.png\n")
cat("  - Scatter plots: top_significant_associations_scatter.png\n\n")

cat("Summary statistics:\n")
cat("  Total regressions:", nrow(all_results), "\n")
cat("  Significant (p < 0.05):", sum(all_results$p_value < 0.05), "\n")
cat("  Highly significant (p < 0.001):", sum(all_results$p_value < 0.001), "\n")
cat("  Sample type used:", primary_sample_type, "\n\n")

# Save workspace
save.image("clinical_biomarker_analysis.RData")
cat("Workspace saved: clinical_biomarker_analysis.RData\n")






































