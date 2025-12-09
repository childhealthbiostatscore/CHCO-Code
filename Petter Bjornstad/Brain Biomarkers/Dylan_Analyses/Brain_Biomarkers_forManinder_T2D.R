### Brain Biomarker Analysis - Clinical Variable Associations
### TYPE 2 DIABETES VS. CONTROL ONLY
### Simple regressions: Clinical variables (HbA1c, M/I, CGM, BP) predicting brain biomarkers
### Models: Unadjusted and adjusted (+ age + sex + BMI)
### Groups: All participants, T2D, Lean Control (by sample type)

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
  "Blood Pressure" = c("sbp", "dbp", "map")
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
# *** MODIFIED: FILTER TO T2D AND LEAN CONTROL ONLY ***
# ============================================================

dat_baseline <- dat %>% 
  filter(visit == 'baseline') %>%
  filter(!is.na(ab40_avg_conc)) %>%  # Has Quanterix data
  filter(!is.na(sample_type)) %>%     # Known sample type
  filter(!is.na(age) & !is.na(sex) & !is.na(bmi)) %>%  # Has covariates
  filter(group %in% c("Type 2 Diabetes", "Lean Control"))  # *** T2D AND LC ONLY ***

cat("\nParticipants with Quanterix data and covariates at baseline (T2D & LC only):\n")
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

output_dir <- 'C:/Users/netio/Documents/UofW/Projects/Maninder_Data/Clinical_Biomarker_Associations_T2D/'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ============================================================
# DEFINE ANALYSIS GROUPS - SEPARATE BY SAMPLE TYPE
# *** MODIFIED: T2D AND LC ONLY ***
# ============================================================

# Analyze each sample type separately since they cannot be combined
sample_types <- unique(dat_baseline$sample_type)

cat("Analyzing sample types separately:", paste(sample_types, collapse = ", "), "\n\n")

# For each sample type, create analysis groups
all_analysis_groups <- list()

for(stype in sample_types) {
  
  cat("=== Sample Type:", stype, "===\n")
  
  dat_sample <- dat_baseline %>% filter(sample_type == stype)
  
  # Check which groups are available
  available_groups <- unique(dat_sample$group)
  cat("Available groups:", paste(available_groups, collapse = ", "), "\n")
  
  # Create group list for this sample type
  analysis_groups <- list(
    "All_Participants" = dat_sample
  )
  
  if("Type 2 Diabetes" %in% available_groups) {
    analysis_groups[["Type_2_Diabetes"]] = dat_sample %>% filter(group == "Type 2 Diabetes")
  }
  
  if("Lean Control" %in% available_groups) {
    analysis_groups[["Lean_Control"]] = dat_sample %>% filter(group == "Lean Control")
  }
  
  # Report sample sizes
  cat("\nSample sizes by group (", stype, "):\n", sep = "")
  for(group_name in names(analysis_groups)) {
    cat(sprintf("  %-20s: N = %d\n", group_name, nrow(analysis_groups[[group_name]])))
  }
  cat("\n")
  
  # Store with sample type prefix
  for(group_name in names(analysis_groups)) {
    key <- paste0(stype, "_", group_name)
    all_analysis_groups[[key]] <- analysis_groups[[group_name]]
  }
}

# ============================================================
# FUNCTION TO TEST LOG TRANSFORMATION
# ============================================================

test_log_transformation <- function(x) {
  if(any(x <= 0, na.rm = TRUE)) return(FALSE)
  if(sum(!is.na(x)) < 10) return(FALSE)
  
  skew_raw <- moments::skewness(x, na.rm = TRUE)
  skew_log <- moments::skewness(log(x), na.rm = TRUE)
  
  if(abs(skew_raw) > 1 && abs(skew_log) < abs(skew_raw)) {
    return(TRUE)
  }
  return(FALSE)
}

# ============================================================
# FUNCTION TO RUN REGRESSIONS
# ============================================================

run_regression_analysis <- function(data, predictor, outcome, covariates = NULL, group_name = "") {
  
  if(!predictor %in% names(data) || !outcome %in% names(data)) {
    return(NULL)
  }
  
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
  
  pred_needs_log <- test_log_transformation(analysis_data[[predictor]])
  outcome_needs_log <- test_log_transformation(analysis_data[[outcome]])
  
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
  
  original_working_data <- working_data
  
  working_data[[predictor]] <- scale(working_data[[predictor]])[,1]
  working_data[[outcome]] <- scale(working_data[[outcome]])[,1]
  
  if(!is.null(covar_list)) {
    for(cov in covar_list) {
      if(cov != "sex" && is.numeric(working_data[[cov]])) {
        working_data[[cov]] <- scale(working_data[[cov]])[,1]
      }
    }
  }
  
  if(is.null(covar_list)) {
    formula_str <- paste0(outcome, " ~ ", predictor)
  } else {
    formula_str <- paste0(outcome, " ~ ", predictor, " + ", paste(covar_list, collapse = " + "))
  }
  
  model_std <- lm(as.formula(formula_str), data = working_data)
  model_summary_std <- summary(model_std)
  coef_table_std <- coef(model_summary_std)
  
  predictor_row <- which(rownames(coef_table_std) == predictor)
  if(length(predictor_row) == 0) return(NULL)
  
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
    Beta_Std = coef_table_std[predictor_row, "Estimate"],
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

for(group_name in names(all_analysis_groups)) {
  
  cat("\n=== Analyzing group:", group_name, "===\n")
  current_data <- all_analysis_groups[[group_name]]
  
  if(nrow(current_data) < 10) {
    cat("  Skipping - insufficient sample size (N =", nrow(current_data), ")\n")
    next
  }
  
  available_predictors <- intersect(all_predictors, names(current_data))
  
  if(length(available_predictors) == 0) {
    cat("  No clinical predictors available in dataset\n")
    next
  }
  
  cat("  Available predictors:", length(available_predictors), "\n")
  cat("  Sample size:", nrow(current_data), "\n")
  
  for(predictor in available_predictors) {
    
    n_available <- sum(!is.na(current_data[[predictor]]))
    if(n_available < 10) next
    
    for(outcome in qx_var) {
      
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
    Transformation = case_when(
      Predictor_Log_Transformed & Outcome_Log_Transformed ~ "Both log",
      Predictor_Log_Transformed ~ "Predictor log",
      Outcome_Log_Transformed ~ "Outcome log",
      TRUE ~ "None"
    )
  )

cat("\n\nTotal regressions completed:", nrow(all_results), "\n")
cat("Significant associations (p < 0.05):", sum(all_results$Sig_Flag), "\n")

cat("\nTransformation summary:\n")
transformation_summary <- table(all_results$Transformation, all_results$Model_Type)
print(transformation_summary)
cat("\n")

write.csv(all_results, "all_regression_results_T2D.csv", row.names = FALSE)

# ============================================================
# CREATE SUMMARY TABLES
# ============================================================

cat("##########################################################")
cat("\n### CREATING SUMMARY TABLES AND VISUALIZATIONS")
cat("\n##########################################################\n\n")

significant_results <- all_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

write.csv(significant_results, "significant_associations_T2D.csv", row.names = FALSE)

cat("Significant associations saved\n")
cat("Number of significant associations:", nrow(significant_results), "\n\n")

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

# Define consistent predictor order
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

# Function to create p-value heatmap
create_pvalue_heatmap <- function(results_subset, title, filename, predictor_order = NULL) {
  
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  if(!is.null(predictor_order)) {
    predictor_order <- predictor_order[predictor_order %in% rownames(pval_matrix)]
    if(length(predictor_order) > 0) {
      pval_matrix <- pval_matrix[predictor_order, , drop = FALSE]
    }
  }
  
  pval_transformed <- -log10(pval_matrix)
  pval_transformed[pval_transformed > 3] <- 3
  
  if(nrow(pval_transformed) > 0 && ncol(pval_transformed) > 0) {
    
    pheatmap(pval_transformed,
             cluster_rows = FALSE,
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
  
  beta_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, Beta_Std) %>%
    pivot_wider(names_from = Outcome, values_from = Beta_Std) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  if(!is.null(predictor_order)) {
    predictor_order <- predictor_order[predictor_order %in% rownames(beta_matrix)]
    if(length(predictor_order) > 0) {
      beta_matrix <- beta_matrix[predictor_order, , drop = FALSE]
      pval_matrix <- pval_matrix[predictor_order, , drop = FALSE]
    }
  }
  
  sig_matrix <- matrix("", nrow = nrow(beta_matrix), ncol = ncol(beta_matrix))
  sig_matrix[pval_matrix < 0.001] <- "***"
  sig_matrix[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
  sig_matrix[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"
  
  display_matrix <- matrix(
    paste0(sprintf("%.3f", beta_matrix), "\n", sig_matrix),
    nrow = nrow(beta_matrix)
  )
  
  if(nrow(beta_matrix) > 0 && ncol(beta_matrix) > 0) {
    
    pheatmap(beta_matrix,
             cluster_rows = FALSE,
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
      title <- paste0(group_name, " - ", model_type, " Models")
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

top_for_plots <- significant_results %>%
  arrange(p_value) %>%
  head(12)

if(nrow(top_for_plots) > 0) {
  
  plot_list <- list()
  
  for(i in 1:nrow(top_for_plots)) {
    
    row <- top_for_plots[i, ]
    
    # Get the correct data group
    plot_data <- all_analysis_groups[[row$Group]] %>%
      dplyr::select(x = all_of(row$Predictor), y = all_of(row$Outcome)) %>%
      filter(complete.cases(.))
    
    if(nrow(plot_data) < 10) next
    
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
  
  if(length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 3)
    ggsave("top_significant_associations_scatter_T2D.png", combined_plot,
           width = 15, height = 12, dpi = 300, bg = "white")
    cat("Created: top_significant_associations_scatter_T2D.png\n")
  }
}

# ============================================================
# SUMMARY BY PREDICTOR TYPE
# ============================================================

cat("\n\n=== Summary of associations by predictor category ===\n")

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

write.csv(category_summary, "summary_by_predictor_category_T2D.csv", row.names = FALSE)

cat("\n")
print(kable(category_summary, digits = 3, format = "simple"))

# ============================================================
# FOREST PLOTS
# ============================================================

cat("\n\n##########################################################")
cat("\n### CREATING FOREST PLOTS FOR CROSS-GROUP COMPARISONS")
cat("\n##########################################################\n\n")

adjusted_results <- all_results %>%
  filter(Model_Type == "Adjusted")

# Function to create forest plot
create_forest_plot <- function(predictor, outcome, results_data, title_suffix = "") {
  
  plot_data <- results_data %>%
    filter(Predictor == predictor, Outcome == outcome, Model_Type == "Adjusted") %>%
    mutate(
      CI_lower = Beta_Std - 1.96 * SE_Std,
      CI_upper = Beta_Std + 1.96 * SE_Std,
      Significant = p_value < 0.05,
      Group_Label = gsub("_", " ", Group)
    ) %>%
    arrange(Beta_Std)
  
  if(nrow(plot_data) == 0) return(NULL)
  
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
    ) +
    geom_text(aes(label = sprintf("N=%d\np=%.3f", N, p_value)),
              x = max(plot_data$CI_upper) * 1.1,
              size = 3, hjust = 0)
  
  return(p)
}

# Create forest plots for top associations
top_associations <- significant_results %>%
  filter(Model_Type == "Adjusted") %>%
  group_by(Predictor, Outcome) %>%
  summarise(
    Min_p = min(p_value),
    N_Groups = n(),
    .groups = "drop"
  ) %>%
  filter(N_Groups >= 2) %>%
  arrange(Min_p) %>%
  head(15)

if(nrow(top_associations) > 0) {
  
  for(i in 1:nrow(top_associations)) {
    predictor <- top_associations$Predictor[i]
    outcome <- top_associations$Outcome[i]
    
    p <- create_forest_plot(predictor, outcome, all_results)
    
    if(!is.null(p)) {
      filename <- paste0("forest_plot_", 
                         gsub("[^A-Za-z0-9]", "_", predictor), "_",
                         gsub("[^A-Za-z0-9]", "_", outcome), "_T2D.png")
      ggsave(filename, p, width = 10, height = 6, dpi = 300, bg = "white")
      cat("Created:", filename, "\n")
    }
  }
}

# ============================================================
# T2D-SPECIFIC ASSOCIATIONS
# ============================================================

cat("\n\n=== Identifying T2D-specific associations ===\n")

# Find associations significant in T2D but not in Lean Control
t2d_lc_comparison <- adjusted_results %>%
  filter(grepl("Type_2_Diabetes|Lean_Control", Group)) %>%
  dplyr::select(Group, Predictor, Outcome, Beta_Std, p_value, SE_Std, Sample_Type) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(Beta_Std, p_value, SE_Std),
    names_sep = "_"
  )

# Get column names dynamically
beta_t2d_cols <- grep("Beta_Std.*Type_2", names(t2d_lc_comparison), value = TRUE)
pval_t2d_cols <- grep("p_value.*Type_2", names(t2d_lc_comparison), value = TRUE)
beta_lc_cols <- grep("Beta_Std.*Lean", names(t2d_lc_comparison), value = TRUE)
pval_lc_cols <- grep("p_value.*Lean", names(t2d_lc_comparison), value = TRUE)

if(length(beta_t2d_cols) > 0 && length(beta_lc_cols) > 0) {
  
  for(sample_type in unique(dat_baseline$sample_type)) {
    
    # Filter to this sample type
    t2d_specific <- t2d_lc_comparison %>%
      filter(Sample_Type == sample_type)
    
    if(nrow(t2d_specific) == 0) next
    
    # Find the columns for this sample type
    beta_t2d_col <- beta_t2d_cols[grepl(sample_type, beta_t2d_cols)][1]
    pval_t2d_col <- pval_t2d_cols[grepl(sample_type, pval_t2d_cols)][1]
    beta_lc_col <- beta_lc_cols[grepl(sample_type, beta_lc_cols)][1]
    pval_lc_col <- pval_lc_cols[grepl(sample_type, pval_lc_cols)][1]
    
    if(!is.na(beta_t2d_col) && !is.na(pval_t2d_col) && 
       !is.na(beta_lc_col) && !is.na(pval_lc_col)) {
      
      t2d_specific <- t2d_specific %>%
        mutate(
          T2D_specific = (.data[[pval_t2d_col]] < 0.05) & 
            (is.na(.data[[pval_lc_col]]) | .data[[pval_lc_col]] >= 0.05)
        ) %>%
        filter(T2D_specific) %>%
        arrange(.data[[pval_t2d_col]])
      
      if(nrow(t2d_specific) > 0) {
        filename <- paste0("T2D_specific_associations_", sample_type, ".csv")
        write.csv(t2d_specific, filename, row.names = FALSE)
        cat("Found", nrow(t2d_specific), "T2D-specific associations for", sample_type, "\n")
      }
    }
  }
}

# ============================================================
# EFFECT SIZE COMPARISON PLOTS
# ============================================================

cat("\n=== Creating effect size comparison plots ===\n")

# Get associations significant in at least one group
sig_in_any <- adjusted_results %>%
  group_by(Predictor, Outcome) %>%
  filter(any(p_value < 0.05)) %>%
  ungroup()

if(nrow(sig_in_any) > 20) {
  
  # Top 20 by minimum p-value
  top_20_pairs <- sig_in_any %>%
    group_by(Predictor, Outcome) %>%
    summarise(
      Min_p = min(p_value, na.rm = TRUE),
      N_Groups = n(),
      .groups = "drop"
    ) %>%
    arrange(Min_p) %>%
    head(20)
  
  plot_data_dots <- adjusted_results %>%
    inner_join(top_20_pairs %>% dplyr::select(Predictor, Outcome), 
               by = c("Predictor", "Outcome")) %>%
    mutate(
      Combo = paste(Predictor, "→", Outcome),
      Group_Label = gsub("_", " ", Group),
      Significant = p_value < 0.05
    )
  
  p_dots <- ggplot(plot_data_dots, 
                   aes(x = Beta_Std, y = reorder(Combo, Beta_Std), color = Group_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = ifelse(Significant, 3, 2), 
                   shape = ifelse(Significant, 16, 1)),
               position = position_dodge(width = 0.5), alpha = 0.8) +
    scale_size_identity() +
    scale_shape_identity() +
    labs(
      title = "Effect Size Comparison: T2D vs Control (Top 20 Associations)",
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
  
  ggsave("effect_size_dot_plot_comparison_T2D.png", p_dots,
         width = 12, height = 10, dpi = 300, bg = "white")
  cat("Created: effect_size_dot_plot_comparison_T2D.png\n")
}

# ============================================================
# INTERACTION ANALYSIS
# ============================================================

cat("\n\n##########################################################")
cat("\n### INTERACTION ANALYSES - T2D vs. CONTROL")
cat("\n##########################################################\n\n")

test_interaction_two_groups <- function(data, predictor, outcome, group1, group2, covariates = c("age", "sex", "bmi")) {
  
  analysis_data <- data %>%
    filter(group %in% c(group1, group2)) %>%
    dplyr::select(all_of(c(predictor, outcome, "group", covariates))) %>%
    filter(complete.cases(.))
  
  if(nrow(analysis_data) < 20) return(NULL)
  
  analysis_data$group <- factor(analysis_data$group, levels = c(group1, group2))
  
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
  
  original_working_data <- analysis_data
  
  analysis_data[[predictor]] <- scale(analysis_data[[predictor]])[,1]
  analysis_data[[outcome]] <- scale(analysis_data[[outcome]])[,1]
  
  for(cov in covariates) {
    if(cov != "sex" && is.numeric(analysis_data[[cov]])) {
      analysis_data[[cov]] <- scale(analysis_data[[cov]])[,1]
    }
  }
  
  formula_main <- as.formula(paste0(outcome, " ~ ", predictor, " + group + ", 
                                    paste(covariates, collapse = " + ")))
  model_main <- lm(formula_main, data = analysis_data)
  
  formula_int <- as.formula(paste0(outcome, " ~ ", predictor, " * group + ", 
                                   paste(covariates, collapse = " + ")))
  model_int <- lm(formula_int, data = analysis_data)
  
  anova_result <- anova(model_main, model_int)
  interaction_p <- anova_result$`Pr(>F)`[2]
  
  coef_summary <- summary(model_int)$coefficients
  
  predictor_row <- which(rownames(coef_summary) == predictor)
  if(length(predictor_row) == 0) return(NULL)
  
  main_effect <- coef_summary[predictor_row, "Estimate"]
  main_effect_p <- coef_summary[predictor_row, "Pr(>|t|)"]
  
  interaction_term <- paste0(predictor, ":group", group2)
  int_row <- which(rownames(coef_summary) == interaction_term)
  
  if(length(int_row) == 0) return(NULL)
  
  interaction_beta <- coef_summary[int_row, "Estimate"]
  interaction_se <- coef_summary[int_row, "Std. Error"]
  interaction_t <- coef_summary[int_row, "t value"]
  interaction_p_coef <- coef_summary[int_row, "Pr(>|t|)"]
  
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

# Run interaction tests for each sample type
interaction_results_all <- data.frame()

for(sample_type in unique(dat_baseline$sample_type)) {
  
  cat("\n=== Testing interactions for sample type:", sample_type, "===\n")
  
  # Get data for this sample type
  sample_data <- dat_baseline %>% filter(sample_type == sample_type)
  
  # Check if both groups exist
  groups_present <- unique(sample_data$group)
  
  if(!"Type 2 Diabetes" %in% groups_present || !"Lean Control" %in% groups_present) {
    cat("Skipping - not both groups present\n")
    next
  }
  
  # Get predictor-outcome combinations tested in both groups
  predictors_to_test <- all_results %>%
    filter(Model_Type == "Adjusted") %>%
    filter(Sample_Type == sample_type) %>%
    filter(grepl("Type_2|Lean", Group)) %>%
    group_by(Predictor, Outcome) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    distinct(Predictor, Outcome)
  
  cat("Testing", nrow(predictors_to_test), "predictor-outcome pairs\n")
  
  for(i in 1:nrow(predictors_to_test)) {
    
    predictor <- predictors_to_test$Predictor[i]
    outcome <- predictors_to_test$Outcome[i]
    
    int_result <- test_interaction_two_groups(
      data = sample_data,
      predictor = predictor,
      outcome = outcome,
      group1 = "Lean Control",
      group2 = "Type 2 Diabetes",
      covariates = c("age", "sex", "bmi")
    )
    
    if(!is.null(int_result)) {
      int_result$Sample_Type <- sample_type
      interaction_results_all <- rbind(interaction_results_all, int_result)
    }
  }
  
  cat("Completed", nrow(filter(interaction_results_all, Sample_Type == sample_type)), 
      "interaction tests for", sample_type, "\n")
}

if(nrow(interaction_results_all) > 0) {
  
  interaction_results_all <- interaction_results_all %>%
    mutate(
      Significant = Interaction_p < 0.05,
      Direction_Change = sign(Effect_in_Group1) != sign(Effect_in_Group2)
    )
  
  write.csv(interaction_results_all, "interaction_analysis_results_T2D.csv", row.names = FALSE)
  
  significant_interactions <- interaction_results_all %>%
    filter(Interaction_p < 0.05) %>%
    arrange(Interaction_p)
  
  write.csv(significant_interactions, "significant_interactions_T2D.csv", row.names = FALSE)
  
  cat("\nInteraction results:\n")
  cat("  Total tests:", nrow(interaction_results_all), "\n")
  cat("  Significant (p < 0.05):", nrow(significant_interactions), "\n\n")
  
  # Plot significant interactions
  if(nrow(significant_interactions) > 0) {
    
    cat("=== Creating interaction plots ===\n")
    
    top_interactions <- significant_interactions %>% head(min(12, nrow(significant_interactions)))
    
    interaction_plots <- list()
    
    for(i in 1:nrow(top_interactions)) {
      
      row <- top_interactions[i, ]
      
      plot_data <- dat_baseline %>%
        filter(sample_type == row$Sample_Type) %>%
        dplyr::select(predictor = all_of(row$Predictor), 
                      outcome = all_of(row$Outcome),
                      group) %>%
        filter(complete.cases(.)) %>%
        filter(group %in% c("Type 2 Diabetes", "Lean Control"))
      
      if(nrow(plot_data) < 20) next
      
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
      
      p <- ggplot(plot_data, aes(x = predictor, y = outcome, color = group)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
        scale_color_manual(values = c("Type 2 Diabetes" = "#377EB8", 
                                      "Lean Control" = "#4DAF4A")) +
        labs(
          title = paste0(row$Predictor, " → ", row$Outcome),
          subtitle = sprintf("Interaction p = %.4f | Sample: %s", 
                             row$Interaction_p, row$Sample_Type),
          x = x_label,
          y = y_label,
          color = "Group"
        ) +
        theme_bw() +
        theme(legend.position = "bottom")
      
      interaction_plots[[i]] <- p
      
      filename <- paste0("interaction_", 
                         gsub("[^A-Za-z0-9]", "_", row$Predictor), "_",
                         gsub("[^A-Za-z0-9]", "_", row$Outcome), "_T2D.png")
      ggsave(filename, p, width = 8, height = 6, dpi = 300, bg = "white")
    }
    
    if(length(interaction_plots) > 0) {
      combined <- wrap_plots(interaction_plots, ncol = 3)
      ggsave("significant_interactions_combined_T2D.png", combined,
             width = 18, height = 12, dpi = 300, bg = "white")
      cat("Created: significant_interactions_combined_T2D.png\n")
    }
  }
}

cat("\n\n##########################################################")
cat("\n### ANALYSIS COMPLETE - T2D VS. CONTROL")
cat("\n##########################################################\n\n")

cat("Results saved:\n")
cat("  - all_regression_results_T2D.csv\n")
cat("  - significant_associations_T2D.csv\n")
cat("  - summary_by_predictor_category_T2D.csv\n")
if(nrow(interaction_results_all) > 0) {
  cat("  - interaction_analysis_results_T2D.csv\n")
  cat("  - significant_interactions_T2D.csv\n")
}
cat("  - Heatmaps (p-values and betas)\n")
cat("  - Forest plots\n")
cat("  - Scatter plots\n")
cat("  - Interaction plots\n\n")

cat("Summary statistics:\n")
cat("  Total regressions:", nrow(all_results), "\n")
cat("  Significant (p < 0.05):", sum(all_results$p_value < 0.05), "\n")
cat("  Highly significant (p < 0.001):", sum(all_results$p_value < 0.001), "\n")
if(nrow(interaction_results_all) > 0) {
  cat("  Interaction tests:", nrow(interaction_results_all), "\n")
  cat("  Significant interactions:", sum(interaction_results_all$Interaction_p < 0.05), "\n")
}

cat("\nSample type distribution:\n")
sample_summary <- dat_baseline %>%
  group_by(group, sample_type) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(sample_type, group)
print(sample_summary)

save.image("clinical_biomarker_analysis_T2D.RData")
cat("\nWorkspace saved: clinical_biomarker_analysis_T2D.RData\n")






### Brain Biomarker Analysis - T2D vs Control
### RH2 STUDY ONLY - SERUM SAMPLES
### Simplified analysis focusing on one study with both groups

library(tidyverse)
library(broom)
library(patchwork)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(knitr)
library(kableExtra)
library(moments)

# Define Quanterix brain biomarkers
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Define clinical predictors
clinical_predictors <- list(
  "Glycemic Control" = c("hba1c", "fbg"),
  "Insulin Sensitivity" = c("avg_m_fsoc", "homa_ir", "adipose_ir", "search_eis"),
  "CGM Metrics" = c("cgm_mean_glucose", "cgm_sd", "cgm_cv", "time_in_range", 
                    "time_above_range", "time_below_range"),
  "Blood Pressure" = c("sbp", "dbp", "map")
)

all_predictors <- unlist(clinical_predictors, use.names = FALSE)

# ============================================================
# DATA LOADING AND PREPARATION
# ============================================================

cat("\n##########################################################")
cat("\n### LOADING DATA - RH2 SERUM ONLY")
cat("\n##########################################################\n\n")

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

# Summarize by participant and visit
dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

# Filter to RH2 study only (serum), baseline, with Quanterix and covariates
dat_baseline <- dat %>% 
  filter(study %in% c("RH2", "RENAL-HEIR", "RENAL-HEIRitage")) %>%  # All serum studies
  filter(visit == 'baseline') %>%
  filter(!is.na(ab40_avg_conc)) %>%
  filter(!is.na(age) & !is.na(sex) & !is.na(bmi)) %>%
  filter(group %in% c("Type 2 Diabetes", "Lean Control")) %>%
  mutate(sample_type = "serum")  # All are serum

cat("RH2/RENAL-HEIR Participants (serum) at baseline:\n")
cat("  Total:", nrow(dat_baseline), "\n\n")

cat("Group distribution:\n")
print(table(dat_baseline$group))
cat("\n")

cat("Study distribution:\n")
print(table(dat_baseline$study, dat_baseline$group))
cat("\n")

# ============================================================
# SETUP OUTPUT DIRECTORY
# ============================================================

output_dir <- 'C:/Users/netio/Documents/UofW/Projects/Maninder_Data/Clinical_Biomarker_Associations_T2D_RH2/'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ============================================================
# DEFINE ANALYSIS GROUPS
# ============================================================

analysis_groups <- list(
  "All_Participants" = dat_baseline,
  "Type_2_Diabetes" = dat_baseline %>% filter(group == "Type 2 Diabetes"),
  "Lean_Control" = dat_baseline %>% filter(group == "Lean Control")
)

cat("Sample sizes:\n")
for(group_name in names(analysis_groups)) {
  cat(sprintf("  %-20s: N = %d\n", group_name, nrow(analysis_groups[[group_name]])))
}
cat("\n")

# ============================================================
# FUNCTIONS
# ============================================================

test_log_transformation <- function(x) {
  if(any(x <= 0, na.rm = TRUE)) return(FALSE)
  if(sum(!is.na(x)) < 5) return(FALSE)  # Changed from 10 to 5
  
  skew_raw <- moments::skewness(x, na.rm = TRUE)
  skew_log <- moments::skewness(log(x), na.rm = TRUE)
  
  if(abs(skew_raw) > 1 && abs(skew_log) < abs(skew_raw)) {
    return(TRUE)
  }
  return(FALSE)
}

run_regression_analysis <- function(data, predictor, outcome, covariates = NULL, group_name = "") {
  
  if(!predictor %in% names(data) || !outcome %in% names(data)) {
    return(NULL)
  }
  
  if(is.null(covariates)) {
    analysis_data <- data %>% 
      dplyr::select(all_of(c(predictor, outcome))) %>%
      filter(complete.cases(.))
    
    if(nrow(analysis_data) < 5) return(NULL)  # Changed from 10 to 5
    
    model_type <- "Unadjusted"
    covar_list <- NULL
    
  } else {
    analysis_data <- data %>% 
      dplyr::select(all_of(c(predictor, outcome, covariates))) %>%
      filter(complete.cases(.))
    
    if(nrow(analysis_data) < 5) return(NULL)  # Changed from 10 to 5
    
    model_type <- "Adjusted"
    covar_list <- covariates
  }
  
  pred_needs_log <- test_log_transformation(analysis_data[[predictor]])
  outcome_needs_log <- test_log_transformation(analysis_data[[outcome]])
  
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
  
  original_working_data <- working_data
  
  working_data[[predictor]] <- scale(working_data[[predictor]])[,1]
  working_data[[outcome]] <- scale(working_data[[outcome]])[,1]
  
  if(!is.null(covar_list)) {
    for(cov in covar_list) {
      if(cov != "sex" && is.numeric(working_data[[cov]])) {
        working_data[[cov]] <- scale(working_data[[cov]])[,1]
      }
    }
  }
  
  if(is.null(covar_list)) {
    formula_str <- paste0(outcome, " ~ ", predictor)
  } else {
    formula_str <- paste0(outcome, " ~ ", predictor, " + ", paste(covar_list, collapse = " + "))
  }
  
  model_std <- lm(as.formula(formula_str), data = working_data)
  model_summary_std <- summary(model_std)
  coef_table_std <- coef(model_summary_std)
  
  predictor_row <- which(rownames(coef_table_std) == predictor)
  if(length(predictor_row) == 0) return(NULL)
  
  model_raw <- lm(as.formula(formula_str), data = original_working_data)
  model_summary_raw <- summary(model_raw)
  coef_table_raw <- coef(model_summary_raw)
  
  results <- data.frame(
    Group = group_name,
    Sample_Type = "serum",
    Predictor = predictor,
    Outcome = outcome,
    Predictor_Used = pred_var_name,
    Outcome_Used = outcome_var_name,
    Predictor_Log_Transformed = pred_needs_log,
    Outcome_Log_Transformed = outcome_needs_log,
    Model_Type = model_type,
    N = nrow(analysis_data),
    Beta_Raw = coef_table_raw[predictor_row, "Estimate"],
    Beta_Std = coef_table_std[predictor_row, "Estimate"],
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
  
  if(nrow(current_data) < 5) {  # Changed from 10 to 5
    cat("  Skipping - insufficient sample size\n")
    next
  }
  
  available_predictors <- intersect(all_predictors, names(current_data))
  cat("  Available predictors:", length(available_predictors), "\n")
  cat("  Sample size:", nrow(current_data), "\n")
  
  for(predictor in available_predictors) {
    
    n_available <- sum(!is.na(current_data[[predictor]]))
    if(n_available < 5) next  # Changed from 10 to 5
    
    for(outcome in qx_var) {
      
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
    Transformation = case_when(
      Predictor_Log_Transformed & Outcome_Log_Transformed ~ "Both log",
      Predictor_Log_Transformed ~ "Predictor log",
      Outcome_Log_Transformed ~ "Outcome log",
      TRUE ~ "None"
    )
  )

cat("\n\nTotal regressions completed:", nrow(all_results), "\n")
cat("Significant associations (p < 0.05):", sum(all_results$Sig_Flag), "\n\n")

write.csv(all_results, "all_regression_results_T2D_RH2_serum.csv", row.names = FALSE)

# ============================================================
# CREATE SUMMARY TABLES
# ============================================================

significant_results <- all_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

write.csv(significant_results, "significant_associations_T2D_RH2.csv", row.names = FALSE)

cat("Significant associations:", nrow(significant_results), "\n\n")

if(nrow(significant_results) > 0) {
  cat("=== TOP 20 MOST SIGNIFICANT ASSOCIATIONS ===\n")
  top_20 <- significant_results %>%
    head(20) %>%
    dplyr::select(Group, Model_Type, Predictor, Outcome, N, Beta_Std, p_value, R_squared)
  
  print(kable(top_20, digits = c(0, 0, 0, 0, 0, 3, 4, 3), format = "simple"))
}

cat("\n\n##########################################################")
cat("\n### ANALYSIS COMPLETE - RH2 SERUM ONLY")
cat("\n##########################################################\n\n")

cat("Results saved to:\n")
cat("  - all_regression_results_T2D_RH2_serum.csv\n")
cat("  - significant_associations_T2D_RH2.csv\n\n")

cat("Groups analyzed:\n")
for(group_name in names(analysis_groups)) {
  n <- nrow(analysis_groups[[group_name]])
  n_results <- sum(all_results$Group == group_name)
  cat(sprintf("  %-20s: N = %2d participants, %4d regressions\n", 
              group_name, n, n_results))
}

save.image("clinical_biomarker_analysis_T2D_RH2_serum.RData")
cat("\nWorkspace saved\n")






















































### Create Visualizations from T2D RH2 Results
### Standalone script to generate heatmaps, forest plots, and scatter plots

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(patchwork)

# ============================================================
# SETUP
# ============================================================

setwd('C:/Users/netio/Documents/UofW/Projects/Maninder_Data/Clinical_Biomarker_Associations_T2D_RH2/')

cat("\n##########################################################")
cat("\n### CREATING VISUALIZATIONS FROM RESULTS")
cat("\n##########################################################\n\n")

# Read results
all_results <- read.csv("all_regression_results_T2D_RH2_serum.csv", stringsAsFactors = FALSE)

cat("Loaded results:\n")
cat("  Total rows:", nrow(all_results), "\n")
cat("  Groups:", paste(unique(all_results$Group), collapse = ", "), "\n")
cat("  Significant (p<0.05):", sum(all_results$p_value < 0.05, na.rm = TRUE), "\n\n")

# Define clinical predictors for ordering
clinical_predictors <- list(
  "Glycemic Control" = c("hba1c", "fbg"),
  "Insulin Sensitivity" = c("avg_m_fsoc", "homa_ir", "adipose_ir", "search_eis"),
  "CGM Metrics" = c("cgm_mean_glucose", "cgm_sd", "cgm_cv", "time_in_range", 
                    "time_above_range", "time_below_range"),
  "Blood Pressure" = c("sbp", "dbp", "map")
)

# Create consistent predictor ordering
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

# ============================================================
# HEATMAP FUNCTIONS
# ============================================================

create_pvalue_heatmap <- function(results_subset, title, filename, predictor_order = NULL) {
  
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  if(!is.null(predictor_order)) {
    predictor_order <- predictor_order[predictor_order %in% rownames(pval_matrix)]
    if(length(predictor_order) > 0) {
      pval_matrix <- pval_matrix[predictor_order, , drop = FALSE]
    }
  }
  
  pval_transformed <- -log10(pval_matrix)
  pval_transformed[pval_transformed > 3] <- 3
  
  if(nrow(pval_transformed) > 0 && ncol(pval_transformed) > 0) {
    
    pheatmap(pval_transformed,
             cluster_rows = FALSE,
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

create_beta_heatmap <- function(results_subset, title, filename, predictor_order = NULL) {
  
  beta_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, Beta_Std) %>%
    pivot_wider(names_from = Outcome, values_from = Beta_Std) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  pval_matrix <- results_subset %>%
    dplyr::select(Predictor_Used, Outcome, p_value) %>%
    pivot_wider(names_from = Outcome, values_from = p_value) %>%
    column_to_rownames("Predictor_Used") %>%
    as.matrix()
  
  if(!is.null(predictor_order)) {
    predictor_order <- predictor_order[predictor_order %in% rownames(beta_matrix)]
    if(length(predictor_order) > 0) {
      beta_matrix <- beta_matrix[predictor_order, , drop = FALSE]
      pval_matrix <- pval_matrix[predictor_order, , drop = FALSE]
    }
  }
  
  sig_matrix <- matrix("", nrow = nrow(beta_matrix), ncol = ncol(beta_matrix))
  sig_matrix[pval_matrix < 0.001] <- "***"
  sig_matrix[pval_matrix >= 0.001 & pval_matrix < 0.01] <- "**"
  sig_matrix[pval_matrix >= 0.01 & pval_matrix < 0.05] <- "*"
  
  display_matrix <- matrix(
    paste0(sprintf("%.3f", beta_matrix), "\n", sig_matrix),
    nrow = nrow(beta_matrix)
  )
  
  if(nrow(beta_matrix) > 0 && ncol(beta_matrix) > 0) {
    
    pheatmap(beta_matrix,
             cluster_rows = FALSE,
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

# ============================================================
# CREATE HEATMAPS
# ============================================================

cat("\n##########################################################")
cat("\n### CREATING HEATMAPS")
cat("\n##########################################################\n\n")

for(group_name in unique(all_results$Group)) {
  for(model_type in c("Unadjusted", "Adjusted")) {
    
    results_subset <- all_results %>%
      filter(Group == group_name, Model_Type == model_type)
    
    if(nrow(results_subset) > 0) {
      
      # P-value heatmap
      title_pval <- paste0(group_name, " - ", model_type, "\n(-log10 p-values)")
      filename_pval <- paste0("heatmap_pvalues_", gsub(" ", "_", group_name), "_", model_type, ".png")
      
      success_pval <- create_pvalue_heatmap(results_subset, title_pval, filename_pval, predictor_order = predictor_order_all)
      if(success_pval) {
        cat("Created:", filename_pval, "\n")
      }
      
      # Beta heatmap
      title_beta <- paste0(group_name, " - ", model_type)
      filename_beta <- paste0("heatmap_beta_", gsub(" ", "_", group_name), "_", model_type, ".png")
      
      success_beta <- create_beta_heatmap(results_subset, title_beta, filename_beta, predictor_order = predictor_order_all)
      if(success_beta) {
        cat("Created:", filename_beta, "\n")
      }
    }
  }
}

# ============================================================
# CREATE FOREST PLOTS
# ============================================================

cat("\n\n##########################################################")
cat("\n### CREATING FOREST PLOTS")
cat("\n##########################################################\n\n")

create_forest_plot <- function(predictor, outcome, results_data, model_type = "Adjusted") {
  
  plot_data <- results_data %>%
    filter(Predictor == predictor, Outcome == outcome, Model_Type == model_type) %>%
    mutate(
      CI_lower = Beta_Std - 1.96 * SE_Std,
      CI_upper = Beta_Std + 1.96 * SE_Std,
      Significant = p_value < 0.05,
      Group_Label = case_when(
        Group == "All_Participants" ~ "All Participants",
        Group == "Type_2_Diabetes" ~ "Type 2 Diabetes",
        Group == "Lean_Control" ~ "Lean Control",
        TRUE ~ Group
      )
    ) %>%
    arrange(Beta_Std)
  
  if(nrow(plot_data) == 0) return(NULL)
  
  p <- ggplot(plot_data, aes(x = Beta_Std, y = reorder(Group_Label, Beta_Std))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Significant),
                   height = 0.2, linewidth = 1) +
    geom_point(aes(color = Significant, size = Significant), shape = 18) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3), guide = "none") +
    labs(
      title = paste0(predictor, " → ", outcome),
      subtitle = paste0("Standardized beta coefficients with 95% CI (", model_type, " model)"),
      x = "Standardized Beta Coefficient",
      y = NULL,
      color = "Significance"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 10)
    ) +
    geom_text(aes(label = sprintf("N=%d\np=%.3f", N, p_value)),
              x = max(plot_data$CI_upper) * 1.15,
              size = 3, hjust = 0)
  
  return(p)
}

# Get adjusted results
adjusted_results <- all_results %>%
  filter(Model_Type == "Adjusted")

# Top predictor-outcome pairs (by minimum p-value across groups)
top_pairs <- adjusted_results %>%
  group_by(Predictor, Outcome) %>%
  summarise(
    Min_p = min(p_value, na.rm = TRUE),
    Max_abs_beta = max(abs(Beta_Std), na.rm = TRUE),
    N_Groups = n(),
    .groups = "drop"
  ) %>%
  filter(N_Groups >= 2) %>%  # Must be tested in at least 2 groups
  arrange(Min_p) %>%
  head(15)

cat("Creating forest plots for top", nrow(top_pairs), "associations...\n\n")

forest_plots <- list()

if(nrow(top_pairs) > 0) {
  
  for(i in 1:nrow(top_pairs)) {
    
    predictor <- top_pairs$Predictor[i]
    outcome <- top_pairs$Outcome[i]
    
    p <- create_forest_plot(predictor, outcome, all_results, model_type = "Adjusted")
    
    if(!is.null(p)) {
      forest_plots[[i]] <- p
      
      # Save individual plot
      filename <- paste0("forest_plot_", 
                         gsub("[^A-Za-z0-9]", "_", predictor), "_",
                         gsub("[^A-Za-z0-9]", "_", outcome), ".png")
      ggsave(filename, p, width = 10, height = 5, dpi = 300, bg = "white")
      cat("Created:", filename, "\n")
    }
  }
  
  # Create combined forest plot grid
  if(length(forest_plots) > 0) {
    n_plots <- min(12, length(forest_plots))
    combined_forest <- wrap_plots(forest_plots[1:n_plots], ncol = 2)
    
    ggsave("forest_plots_combined.png", combined_forest,
           width = 16, height = 4 * ceiling(n_plots/2), dpi = 300, bg = "white")
    
    cat("\nCreated: forest_plots_combined.png\n")
  }
} else {
  cat("No predictor-outcome pairs tested in multiple groups\n")
}

# ============================================================
# CREATE COMPARISON SUMMARY TABLE
# ============================================================

cat("\n\n##########################################################")
cat("\n### CREATING COMPARISON SUMMARY")
cat("\n##########################################################\n\n")

comparison_summary <- adjusted_results %>%
  mutate(Group_Type = case_when(
    Group == "Type_2_Diabetes" ~ "T2D",
    Group == "Lean_Control" ~ "LC",
    Group == "All_Participants" ~ "All"
  )) %>%
  dplyr::select(Predictor, Outcome, Group_Type, Beta_Std, p_value, N) %>%
  pivot_wider(
    names_from = Group_Type,
    values_from = c(Beta_Std, p_value, N)
  )

# Add comparison flags
if("Beta_Std_T2D" %in% names(comparison_summary) && "Beta_Std_LC" %in% names(comparison_summary)) {
  comparison_summary <- comparison_summary %>%
    mutate(
      T2D_sig = !is.na(p_value_T2D) & p_value_T2D < 0.05,
      LC_sig = !is.na(p_value_LC) & p_value_LC < 0.05,
      T2D_specific = T2D_sig & !LC_sig,
      Both_sig = T2D_sig & LC_sig,
      Beta_difference = abs(Beta_Std_T2D - coalesce(Beta_Std_LC, 0))
    ) %>%
    arrange(p_value_T2D)
  
  write.csv(comparison_summary, "T2D_vs_LC_comparison.csv", row.names = FALSE)
  
  cat("Comparison summary:\n")
  cat("  Total pairs compared:", nrow(comparison_summary), "\n")
  cat("  Significant in T2D:", sum(comparison_summary$T2D_sig, na.rm = TRUE), "\n")
  cat("  Significant in LC:", sum(comparison_summary$LC_sig, na.rm = TRUE), "\n")
  cat("  Significant in both:", sum(comparison_summary$Both_sig, na.rm = TRUE), "\n")
  cat("  T2D-specific:", sum(comparison_summary$T2D_specific, na.rm = TRUE), "\n\n")
  
  cat("Saved: T2D_vs_LC_comparison.csv\n")
  
  # Top T2D-specific associations
  t2d_specific <- comparison_summary %>%
    filter(T2D_specific) %>%
    head(10)
  
  if(nrow(t2d_specific) > 0) {
    cat("\nTop 10 T2D-specific associations:\n")
    print(t2d_specific %>% 
            dplyr::select(Predictor, Outcome, Beta_Std_T2D, p_value_T2D, N_T2D))
  }
}

# ============================================================
# EFFECT SIZE COMPARISON DOT PLOT - ADJUSTED
# ============================================================

cat("\n\n##########################################################")
cat("\n### CREATING EFFECT SIZE COMPARISON PLOT - ADJUSTED")
cat("\n##########################################################\n\n")

if(nrow(top_pairs) > 0) {
  
  plot_data_dots <- adjusted_results %>%
    inner_join(top_pairs %>% dplyr::select(Predictor, Outcome), 
               by = c("Predictor", "Outcome")) %>%
    mutate(
      Combo = paste(Predictor, "→", Outcome),
      Group_Label = case_when(
        Group == "All_Participants" ~ "All",
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
    scale_color_manual(values = c("All" = "black", "T2D" = "#377EB8", "LC" = "#4DAF4A")) +
    scale_size_identity() +
    scale_shape_identity() +
    labs(
      title = "Effect Size Comparison: T2D vs Control (Adjusted)",
      subtitle = "Filled circles = p < 0.05, Open circles = p ≥ 0.05 | Adjusted for age, sex, BMI",
      x = "Standardized Beta Coefficient",
      y = NULL,
      color = "Group"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold")
    )
  
  ggsave("effect_size_comparison_adjusted.png", p_dots,
         width = 12, height = 10, dpi = 300, bg = "white")
  cat("Created: effect_size_comparison_adjusted.png\n")
}

# ============================================================
# EFFECT SIZE COMPARISON DOT PLOT - UNADJUSTED
# ============================================================

cat("\n##########################################################")
cat("\n### CREATING EFFECT SIZE COMPARISON PLOT - UNADJUSTED")
cat("\n##########################################################\n\n")

# Get unadjusted results
unadjusted_results <- all_results %>%
  filter(Model_Type == "Unadjusted")

# Use same top pairs for consistency
if(nrow(top_pairs) > 0) {
  
  plot_data_unadj <- unadjusted_results %>%
    inner_join(top_pairs %>% dplyr::select(Predictor, Outcome), 
               by = c("Predictor", "Outcome")) %>%
    mutate(
      Combo = paste(Predictor, "→", Outcome),
      Group_Label = case_when(
        Group == "All_Participants" ~ "All",
        Group == "Type_2_Diabetes" ~ "T2D",
        Group == "Lean_Control" ~ "LC",
        TRUE ~ Group
      ),
      Significant = p_value < 0.05
    )
  
  p_dots_unadj <- ggplot(plot_data_unadj, 
                         aes(x = Beta_Std, y = reorder(Combo, Beta_Std), color = Group_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = ifelse(Significant, 3, 2), 
                   shape = ifelse(Significant, 16, 1)),
               position = position_dodge(width = 0.5), alpha = 0.8) +
    scale_color_manual(values = c("All" = "black", "T2D" = "#377EB8", "LC" = "#4DAF4A")) +
    scale_size_identity() +
    scale_shape_identity() +
    labs(
      title = "Effect Size Comparison: T2D vs Control (Unadjusted)",
      subtitle = "Filled circles = p < 0.05, Open circles = p ≥ 0.05 | Unadjusted models",
      x = "Standardized Beta Coefficient",
      y = NULL,
      color = "Group"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold")
    )
  
  ggsave("effect_size_comparison_unadjusted.png", p_dots_unadj,
         width = 12, height = 10, dpi = 300, bg = "white")
  cat("Created: effect_size_comparison_unadjusted.png\n")
}

# ============================================================
# SIDE-BY-SIDE COMPARISON: ADJUSTED vs UNADJUSTED
# ============================================================

cat("\n##########################################################")
cat("\n### CREATING SIDE-BY-SIDE COMPARISON PLOT")
cat("\n##########################################################\n\n")

if(nrow(top_pairs) > 0) {
  
  # Combine adjusted and unadjusted
  plot_data_both <- all_results %>%
    inner_join(top_pairs %>% dplyr::select(Predictor, Outcome), 
               by = c("Predictor", "Outcome")) %>%
    mutate(
      Combo = paste(Predictor, "→", Outcome),
      Group_Label = case_when(
        Group == "All_Participants" ~ "All",
        Group == "Type_2_Diabetes" ~ "T2D",
        Group == "Lean_Control" ~ "LC",
        TRUE ~ Group
      ),
      Significant = p_value < 0.05
    )
  
  p_both <- ggplot(plot_data_both, 
                   aes(x = Beta_Std, y = reorder(Combo, Beta_Std), color = Group_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = ifelse(Significant, 3, 2), 
                   shape = ifelse(Significant, 16, 1)),
               position = position_dodge(width = 0.5), alpha = 0.8) +
    scale_color_manual(values = c("All" = "black", "T2D" = "#377EB8", "LC" = "#4DAF4A")) +
    scale_size_identity() +
    scale_shape_identity() +
    facet_wrap(~Model_Type, ncol = 2) +
    labs(
      title = "Effect Size Comparison: Adjusted vs Unadjusted Models",
      subtitle = "Filled circles = p < 0.05, Open circles = p ≥ 0.05",
      x = "Standardized Beta Coefficient",
      y = NULL,
      color = "Group"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold")
    )
  
  ggsave("effect_size_comparison_both_models.png", p_both,
         width = 16, height = 10, dpi = 300, bg = "white")
  cat("Created: effect_size_comparison_both_models.png\n")
}

cat("\n\n##########################################################")
cat("\n### CREATING ALL SIGNIFICANT ASSOCIATIONS PLOT")
cat("\n##########################################################\n\n")

# Get ALL significant associations (in any group, any model)
all_significant_pairs <- all_results %>%
  filter(p_value < 0.05) %>%
  distinct(Predictor, Outcome) %>%
  arrange(Predictor, Outcome)

cat("Found", nrow(all_significant_pairs), "unique predictor-outcome pairs with p < 0.05 in any group/model\n\n")

if(nrow(all_significant_pairs) > 0) {
  
  # For adjusted models
  plot_data_all_sig_adj <- adjusted_results %>%
    inner_join(all_significant_pairs, by = c("Predictor", "Outcome")) %>%
    mutate(
      Combo = paste(Predictor, "→", Outcome),
      Group_Label = case_when(
        Group == "All_Participants" ~ "All",
        Group == "Type_2_Diabetes" ~ "T2D",
        Group == "Lean_Control" ~ "LC",
        TRUE ~ Group
      ),
      Significant = p_value < 0.05
    )
  
  # Sort by mean beta across groups
  combo_order <- plot_data_all_sig_adj %>%
    group_by(Combo) %>%
    summarise(Mean_Beta = mean(Beta_Std, na.rm = TRUE), .groups = "drop") %>%
    arrange(Mean_Beta) %>%
    pull(Combo)
  
  plot_data_all_sig_adj$Combo <- factor(plot_data_all_sig_adj$Combo, levels = combo_order)
  
  p_all_sig_adj <- ggplot(plot_data_all_sig_adj, 
                          aes(x = Beta_Std, y = Combo, color = Group_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = ifelse(Significant, 3, 2), 
                   shape = ifelse(Significant, 16, 1)),
               position = position_dodge(width = 0.5), alpha = 0.8) +
    scale_color_manual(values = c("All" = "black", "T2D" = "#377EB8", "LC" = "#4DAF4A")) +
    scale_size_identity() +
    scale_shape_identity() +
    labs(
      title = "All Significant Associations (Adjusted Models)",
      subtitle = "Filled circles = p < 0.05, Open circles = p ≥ 0.05 | Sorted by mean effect size",
      x = "Standardized Beta Coefficient (adjusted for age, sex, BMI)",
      y = NULL,
      color = "Group"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 6),
      plot.title = element_text(face = "bold")
    )
  
  ggsave("effect_size_all_significant_adjusted.png", p_all_sig_adj,
         width = 12, height = max(10, nrow(all_significant_pairs) * 0.3), dpi = 300, bg = "white")
  cat("Created: effect_size_all_significant_adjusted.png\n")
  
  # For unadjusted models
  plot_data_all_sig_unadj <- unadjusted_results %>%
    inner_join(all_significant_pairs, by = c("Predictor", "Outcome")) %>%
    mutate(
      Combo = paste(Predictor, "→", Outcome),
      Group_Label = case_when(
        Group == "All_Participants" ~ "All",
        Group == "Type_2_Diabetes" ~ "T2D",
        Group == "Lean_Control" ~ "LC",
        TRUE ~ Group
      ),
      Significant = p_value < 0.05
    )
  
  plot_data_all_sig_unadj$Combo <- factor(plot_data_all_sig_unadj$Combo, levels = combo_order)
  
  p_all_sig_unadj <- ggplot(plot_data_all_sig_unadj, 
                            aes(x = Beta_Std, y = Combo, color = Group_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = ifelse(Significant, 3, 2), 
                   shape = ifelse(Significant, 16, 1)),
               position = position_dodge(width = 0.5), alpha = 0.8) +
    scale_color_manual(values = c("All" = "black", "T2D" = "#377EB8", "LC" = "#4DAF4A")) +
    scale_size_identity() +
    scale_shape_identity() +
    labs(
      title = "All Significant Associations (Unadjusted Models)",
      subtitle = "Filled circles = p < 0.05, Open circles = p ≥ 0.05 | Sorted by mean effect size",
      x = "Standardized Beta Coefficient",
      y = NULL,
      color = "Group"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 6),
      plot.title = element_text(face = "bold")
    )
  
  ggsave("effect_size_all_significant_unadjusted.png", p_all_sig_unadj,
         width = 12, height = max(10, nrow(all_significant_pairs) * 0.3), dpi = 300, bg = "white")
  cat("Created: effect_size_all_significant_unadjusted.png\n")
  
  # Summary of negative vs positive associations
  cat("\nBreakdown of significant associations by direction:\n")
  direction_summary <- plot_data_all_sig_adj %>%
    filter(Significant) %>%
    mutate(Direction = ifelse(Beta_Std < 0, "Negative", "Positive")) %>%
    group_by(Direction, Group_Label) %>%
    summarise(N = n(), .groups = "drop")
  
  print(direction_summary)
}

# ============================================================
# SUMMARY
# ============================================================

cat("\n\n##########################################################")
cat("\n### VISUALIZATION COMPLETE")
cat("\n##########################################################\n\n")

cat("Files created:\n")
cat("  - Heatmaps (p-values and betas) for each group and model type\n")
cat("  - Forest plots comparing effect sizes across groups\n")
cat("  - Combined forest plot grid\n")
cat("  - Effect size comparison dot plot\n")
cat("  - T2D vs LC comparison summary table\n\n")
















