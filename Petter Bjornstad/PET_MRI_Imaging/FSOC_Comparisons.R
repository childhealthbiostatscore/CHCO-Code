# ============================================================================
# FSOC Analysis: Tubular Oxygen Consumption Across Disease Groups
# COMPLETE ANALYSIS SCRIPT
# ============================================================================

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(emmeans)
library(car)
library(gtsummary)
library(corrplot)
library(rstatix)
library(ggpubr)
library(patchwork)
library(broom)

# ============================================================================
# 1. DATA LOADING
# ============================================================================

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

# ============================================================================
# 2. CALCULATE AVERAGED FSOC ENDPOINTS
# ============================================================================

calculate_fsoc_averages <- function(df) {
  df <- df %>%
    mutate(
      # Average whole-kidney FSOC (absolute) from left and right
      whole_kidney_fsoc_abs = case_when(
        !is.na(fsoc_l_kidney) & !is.na(fsoc_r_kidney) ~ 
          (fsoc_l_kidney + fsoc_r_kidney) / 2,
        !is.na(fsoc_l_kidney) ~ fsoc_l_kidney,
        !is.na(fsoc_r_kidney) ~ fsoc_r_kidney,
        TRUE ~ NA_real_
      ),
      
      # Average medullary FSOC (absolute) from left and right
      medullary_fsoc_abs = case_when(
        !is.na(fsoc_l_medulla) & !is.na(fsoc_r_medulla) ~ 
          (fsoc_l_medulla + fsoc_r_medulla) / 2,
        !is.na(fsoc_l_medulla) ~ fsoc_l_medulla,
        !is.na(fsoc_r_medulla) ~ fsoc_r_medulla,
        TRUE ~ NA_real_
      ),
      
      # Average cortical FSOC (absolute) from left and right
      cortical_fsoc_abs = case_when(
        !is.na(fsoc_l_cortex) & !is.na(fsoc_r_cortex) ~ 
          (fsoc_l_cortex + fsoc_r_cortex) / 2,
        !is.na(fsoc_l_cortex) ~ fsoc_l_cortex,
        !is.na(fsoc_r_cortex) ~ fsoc_r_cortex,
        TRUE ~ NA_real_
      ),
      
      # Log-transform UACR
      log_UACR = log(acr_u + 1)
    )
  
  return(df)
}

# ============================================================================
# 3. NEGATIVE FSOC VALUES ANALYSIS
# ============================================================================

analyze_negative_fsoc <- function(dat) {
  
  cat("\n========================================\n")
  cat("NEGATIVE FSOC VALUES ANALYSIS\n")
  cat("========================================\n\n")
  
  # Overall summary
  negative_summary <- dat %>%
    filter(!is.na(group)) %>%
    summarise(
      n_total = n(),
      n_wk_negative = sum(whole_kidney_fsoc_abs < 0, na.rm = TRUE),
      pct_wk_negative = round(100 * n_wk_negative / sum(!is.na(whole_kidney_fsoc_abs)), 1),
      n_med_negative = sum(medullary_fsoc_abs < 0, na.rm = TRUE),
      pct_med_negative = round(100 * n_med_negative / sum(!is.na(medullary_fsoc_abs)), 1)
    )
  
  cat("OVERALL SUMMARY:\n")
  cat("Total subjects:", negative_summary$n_total, "\n")
  cat("Whole Kidney FSOC < 0:", negative_summary$n_wk_negative, 
      "(", negative_summary$pct_wk_negative, "%)\n")
  cat("Medullary FSOC < 0:", negative_summary$n_med_negative, 
      "(", negative_summary$pct_med_negative, "%)\n\n")
  
  # Breakdown by disease group (exclude groups with no FSOC data)
  negative_by_group <- dat %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    filter(sum(!is.na(whole_kidney_fsoc_abs)) > 0 | sum(!is.na(medullary_fsoc_abs)) > 0) %>%
    summarise(
      n = n(),
      n_wk_available = sum(!is.na(whole_kidney_fsoc_abs)),
      n_wk_neg = sum(whole_kidney_fsoc_abs < 0, na.rm = TRUE),
      pct_wk_neg = ifelse(n_wk_available > 0, 
                          round(100 * n_wk_neg / n_wk_available, 1), NA_real_),
      n_med_available = sum(!is.na(medullary_fsoc_abs)),
      n_med_neg = sum(medullary_fsoc_abs < 0, na.rm = TRUE),
      pct_med_neg = ifelse(n_med_available > 0,
                           round(100 * n_med_neg / n_med_available, 1), NA_real_),
      min_wk_fsoc = ifelse(n_wk_available > 0,
                           round(min(whole_kidney_fsoc_abs, na.rm = TRUE), 2), NA_real_),
      min_med_fsoc = ifelse(n_med_available > 0,
                            round(min(medullary_fsoc_abs, na.rm = TRUE), 2), NA_real_),
      .groups = "drop"
    ) %>%
    select(-n_wk_available, -n_med_available)
  
  cat("BY DISEASE GROUP:\n")
  print(as.data.frame(negative_by_group))
  
  # Breakdown by sex
  negative_by_sex <- dat %>%
    filter(!is.na(group), !is.na(sex)) %>%
    group_by(sex) %>%
    summarise(
      n = n(),
      n_wk_neg = sum(whole_kidney_fsoc_abs < 0, na.rm = TRUE),
      pct_wk_neg = round(100 * n_wk_neg / sum(!is.na(whole_kidney_fsoc_abs)), 1),
      n_med_neg = sum(medullary_fsoc_abs < 0, na.rm = TRUE),
      pct_med_neg = round(100 * n_med_neg / sum(!is.na(medullary_fsoc_abs)), 1),
      .groups = "drop"
    )
  
  cat("\nBY SEX:\n")
  print(as.data.frame(negative_by_sex))
  
  # Identify specific subjects with negative values
  negative_subjects <- dat %>%
    filter(!is.na(group)) %>%
    filter(whole_kidney_fsoc_abs < 0 | medullary_fsoc_abs < 0) %>%
    select(record_id, group, sex, age, 
           whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    arrange(medullary_fsoc_abs)
  
  cat("\nSUBJECTS WITH NEGATIVE FSOC VALUES:\n")
  cat("(n =", nrow(negative_subjects), ")\n")
  print(as.data.frame(negative_subjects))
  
  # Chi-square test: Is negative FSOC associated with group?
  cat("\n--- Statistical Tests ---\n")
  
  if (sum(dat$medullary_fsoc_abs < 0, na.rm = TRUE) >= 5) {
    dat_test <- dat %>%
      filter(!is.na(group), !is.na(medullary_fsoc_abs)) %>%
      mutate(med_negative = ifelse(medullary_fsoc_abs < 0, "Negative", "Positive"))
    
    chi_result <- chisq.test(table(dat_test$group, dat_test$med_negative))
    cat("\nChi-square: Medullary FSOC negative vs Disease Group\n")
    print(chi_result)
  }
  
  # Compare characteristics of negative vs positive FSOC subjects
  cat("\n--- Characteristics Comparison ---\n")
  
  # Check which variables exist in the dataset
  available_vars <- names(dat)
  
  # Find eGFR variable (could be named differently)
  egfr_var <- case_when(
    "eGFR_CKD_epi" %in% available_vars ~ "eGFR_CKD_epi",
    "gfr_fas_cr" %in% available_vars ~ "gfr_fas_cr",
    "eGFR" %in% available_vars ~ "eGFR",
    "egfr" %in% available_vars ~ "egfr",
    "gfr_raw_plasma" %in% available_vars ~ "gfr_raw_plasma",
    "gfr_bsa_plasma" %in% available_vars ~ "gfr_bsa_plasma",
    TRUE ~ NA_character_
  )
  
  comparison_table <- dat %>%
    filter(!is.na(group), !is.na(medullary_fsoc_abs)) %>%
    mutate(fsoc_status = ifelse(medullary_fsoc_abs < 0, 
                                "Negative FSOC", "Positive FSOC")) %>%
    group_by(fsoc_status) %>%
    summarise(
      n = n(),
      mean_age = round(mean(age, na.rm = TRUE), 1),
      pct_male = round(100 * sum(sex == "Male", na.rm = TRUE) / n(), 1),
      mean_weight = round(mean(weight, na.rm = TRUE), 1),
      mean_bmi = if ("bmi" %in% available_vars) round(mean(bmi, na.rm = TRUE), 1) else NA_real_,
      mean_hba1c = if ("hba1c" %in% available_vars) round(mean(hba1c, na.rm = TRUE), 1) else NA_real_,
      .groups = "drop"
    )
  
  # Add eGFR if available
  if (!is.na(egfr_var)) {
    egfr_summary <- dat %>%
      filter(!is.na(group), !is.na(medullary_fsoc_abs)) %>%
      mutate(fsoc_status = ifelse(medullary_fsoc_abs < 0, 
                                  "Negative FSOC", "Positive FSOC")) %>%
      group_by(fsoc_status) %>%
      summarise(mean_egfr = round(mean(.data[[egfr_var]], na.rm = TRUE), 1), .groups = "drop")
    
    comparison_table <- left_join(comparison_table, egfr_summary, by = "fsoc_status")
    cat("Note: Using", egfr_var, "for eGFR\n")
  } else {
    cat("Note: No eGFR variable found in dataset\n")
  }
  
  print(as.data.frame(comparison_table))
  
  results <- list(
    overall_summary = negative_summary,
    by_group = negative_by_group,
    by_sex = negative_by_sex,
    negative_subjects = negative_subjects,
    comparison = comparison_table
  )
  
  return(invisible(results))
}

# Visualize negative FSOC distribution
plot_negative_fsoc_distribution <- function(dat) {
  
  fsoc_long <- dat %>%
    filter(!is.na(group)) %>%
    dplyr::select(record_id, group,
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    pivot_longer(cols = contains("fsoc"),
                 names_to = "fsoc_type",
                 values_to = "fsoc_value") %>%
    filter(!is.na(fsoc_value)) %>%
    mutate(
      fsoc_type = case_match(fsoc_type,
                             "whole_kidney_fsoc_abs" ~ "Whole Kidney FSOC",
                             "medullary_fsoc_abs" ~ "Medullary FSOC"),
      is_negative = fsoc_value < 0
    )
  
  # Histogram with negative values highlighted
  p1 <- ggplot(fsoc_long, aes(x = fsoc_value, fill = is_negative)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    facet_wrap(~ fsoc_type, scales = "free", ncol = 2) +
    scale_fill_manual(values = c("FALSE" = "#66C2A5", "TRUE" = "#FC8D62"),
                      labels = c("Positive", "Negative"),
                      name = "FSOC Value") +
    theme_bw(base_size = 12) +
    labs(x = "FSOC Value (s⁻¹)", y = "Count",
         title = "Distribution of FSOC Values",
         subtitle = "Red dashed line = 0; Orange = negative values") +
    theme(legend.position = "bottom")
  
  # Bar plot: proportion negative by group
  prop_negative <- fsoc_long %>%
    group_by(group, fsoc_type) %>%
    summarise(
      n_total = n(),
      n_negative = sum(is_negative),
      pct_negative = 100 * n_negative / n_total,
      .groups = "drop"
    )
  
  p2 <- ggplot(prop_negative, aes(x = group, y = pct_negative, fill = group)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = paste0(n_negative, "/", n_total)), 
              vjust = -0.3, size = 3) +
    facet_wrap(~ fsoc_type, ncol = 2) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = "Disease Group", y = "% with Negative FSOC",
         title = "Proportion of Negative FSOC Values by Group")
  
  combined <- p1 / p2 + plot_layout(heights = c(1, 1))
  
  return(combined)
}

# ============================================================================
# 4. OBJECTIVE 1.1: GROUP COMPARISONS (ADJUSTED)
# ============================================================================

analyze_group_differences <- function(dat, outcome_var) {
  
  formula_str <- paste(outcome_var, "~ group + age + sex + weight")
  model <- lm(as.formula(formula_str), data = dat)
  
  anova_result <- Anova(model, type = "III")
  
  emm <- emmeans(model, specs = "group")
  pairwise <- pairs(emm, adjust = "tukey")
  
  list(
    model = model,
    anova = anova_result,
    emmeans = emm,
    pairwise = pairwise,
    summary = summary(model)
  )
}

run_all_group_comparisons <- function(dat) {
  fsoc_endpoints <- c("whole_kidney_fsoc_abs", "medullary_fsoc_abs")
  
  results <- list()
  for (endpoint in fsoc_endpoints) {
    cat("\n========================================\n")
    cat("Analysis for:", endpoint, "\n")
    cat("========================================\n")
    
    results[[endpoint]] <- analyze_group_differences(dat, endpoint)
    
    print(results[[endpoint]]$anova)
    cat("\nPairwise comparisons:\n")
    print(results[[endpoint]]$pairwise)
  }
  
  return(results)
}

# ============================================================================
# 5. OBJECTIVE 1.2: SEX DIFFERENCES
# ============================================================================

analyze_sex_differences <- function(dat, outcome_var) {
  
  formula_str <- paste(outcome_var, "~ group * sex + age + weight")
  model_interaction <- lm(as.formula(formula_str), data = dat)
  
  anova_interaction <- Anova(model_interaction, type = "III")
  
  stratified_results <- dat %>%
    group_by(group) %>%
    do({
      model <- lm(as.formula(paste(outcome_var, "~ sex + age + weight")), 
                  data = .)
      tidy_result <- broom::tidy(model) %>%
        filter(term == "sexMale")
      tidy_result
    })
  
  list(
    interaction_model = model_interaction,
    interaction_anova = anova_interaction,
    stratified = stratified_results
  )
}

run_sex_analyses <- function(dat) {
  fsoc_endpoints <- c("whole_kidney_fsoc_abs", "medullary_fsoc_abs")
  
  results <- list()
  for (endpoint in fsoc_endpoints) {
    cat("\n========================================\n")
    cat("Sex analysis for:", endpoint, "\n")
    cat("========================================\n")
    
    results[[endpoint]] <- analyze_sex_differences(dat, endpoint)
    
    print(results[[endpoint]]$interaction_anova)
    cat("\nStratified results:\n")
    print(results[[endpoint]]$stratified)
  }
  
  return(results)
}

# ============================================================================
# 6. OBJECTIVE 1.3: DXA ASSOCIATIONS
# ============================================================================

analyze_dxa_associations <- function(dat) {
  
  fsoc_endpoints <- c("whole_kidney_fsoc_abs", "medullary_fsoc_abs")
  dxa_vars <- c("dexa_fat_kg", "dexa_trunk_kg", "dexa_lean_kg", "dexa_body_fat")
  
  cor_matrix <- matrix(NA, nrow = length(fsoc_endpoints), ncol = length(dxa_vars))
  p_matrix <- matrix(NA, nrow = length(fsoc_endpoints), ncol = length(dxa_vars))
  
  rownames(cor_matrix) <- fsoc_endpoints
  colnames(cor_matrix) <- dxa_vars
  rownames(p_matrix) <- fsoc_endpoints
  colnames(p_matrix) <- dxa_vars
  
  detailed_results <- list()
  
  for (i in 1:length(fsoc_endpoints)) {
    for (j in 1:length(dxa_vars)) {
      formula_str <- paste(fsoc_endpoints[i], "~", dxa_vars[j], "+ age + sex + group")
      
      model <- lm(as.formula(formula_str), data = dat)
      coef_summary <- summary(model)$coefficients
      
      dxa_row <- which(rownames(coef_summary) == dxa_vars[j])
      if (length(dxa_row) > 0) {
        cor_matrix[i, j] <- coef_summary[dxa_row, "Estimate"]
        p_matrix[i, j] <- coef_summary[dxa_row, "Pr(>|t|)"]
      }
      
      detailed_results[[paste(fsoc_endpoints[i], dxa_vars[j], sep = "_")]] <- model
    }
  }
  
  list(
    coefficient_matrix = cor_matrix,
    p_value_matrix = p_matrix,
    detailed_models = detailed_results
  )
}

# ============================================================================
# 6B. OBJECTIVE 1.4: CLINICAL PARAMETER ASSOCIATIONS
# ============================================================================

analyze_clinical_associations <- function(dat) {
  
  cat("\n========================================\n")
  cat("CLINICAL PARAMETER ASSOCIATIONS\n")
  cat("========================================\n\n")
  
  fsoc_endpoints <- c("whole_kidney_fsoc_abs", "medullary_fsoc_abs")
  
  # Define clinical variables (check which exist in dataset)
  available_vars <- names(dat)
  
  # Map of desired variables to possible names in dataset
  clinical_var_map <- list(
    "eGFR" = c("eGFR_CKD_epi", "gfr_fas_cr", "eGFR", "egfr"),
    "mGFR" = c("gfr_raw_plasma", "gfr_bsa_plasma", "mGFR", "mgfr"),
    "HbA1c" = c("hba1c", "HbA1c", "a1c"),
    "SBP" = c("sbp", "SBP", "systolic_bp", "map_sys"),
    "DBP" = c("dbp", "DBP", "diastolic_bp", "map_dia"),
    "RPF" = c("erpf_raw_plasma", "erpf_bsa_plasma", "rpf", "RPF"),
    "FF" = c("ff", "FF", "filtration_fraction"),
    "ACR" = c("acr_u", "uacr", "ACR", "log_UACR"),
    "Diabetes_Duration" = c("diabetes_duration", "dm_duration", "duration_dm")
  )
  
  # Find which variables exist
  clinical_vars <- c()
  clinical_labels <- c()
  
  for (var_name in names(clinical_var_map)) {
    possible_names <- clinical_var_map[[var_name]]
    found_var <- possible_names[possible_names %in% available_vars][1]
    if (!is.na(found_var)) {
      clinical_vars <- c(clinical_vars, found_var)
      clinical_labels <- c(clinical_labels, var_name)
    }
  }
  
  cat("Clinical variables found:", paste(clinical_labels, collapse = ", "), "\n")
  cat("Variable names used:", paste(clinical_vars, collapse = ", "), "\n\n")
  
  if (length(clinical_vars) == 0) {
    cat("No clinical variables found in dataset.\n")
    return(NULL)
  }
  
  # Create matrices
  cor_matrix <- matrix(NA, nrow = length(fsoc_endpoints), ncol = length(clinical_vars))
  p_matrix <- matrix(NA, nrow = length(fsoc_endpoints), ncol = length(clinical_vars))
  
  rownames(cor_matrix) <- fsoc_endpoints
  colnames(cor_matrix) <- clinical_vars
  rownames(p_matrix) <- fsoc_endpoints
  colnames(p_matrix) <- clinical_vars
  
  detailed_results <- list()
  
  for (i in 1:length(fsoc_endpoints)) {
    for (j in 1:length(clinical_vars)) {
      
      # Use log-transformed ACR if it's the ACR variable
      if (clinical_labels[j] == "ACR" && clinical_vars[j] != "log_UACR") {
        dat_temp <- dat %>% 
          mutate(temp_var = log(.data[[clinical_vars[j]]] + 1))
        formula_str <- paste(fsoc_endpoints[i], "~ temp_var + age + sex + group")
      } else {
        dat_temp <- dat
        formula_str <- paste(fsoc_endpoints[i], "~", clinical_vars[j], "+ age + sex + group")
      }
      
      # Fit model
      tryCatch({
        model <- lm(as.formula(formula_str), data = dat_temp)
        coef_summary <- summary(model)$coefficients
        
        # Get the clinical variable row
        var_row <- ifelse(clinical_labels[j] == "ACR" && clinical_vars[j] != "log_UACR",
                          "temp_var", clinical_vars[j])
        coef_row <- which(rownames(coef_summary) == var_row)
        
        if (length(coef_row) > 0) {
          cor_matrix[i, j] <- coef_summary[coef_row, "Estimate"]
          p_matrix[i, j] <- coef_summary[coef_row, "Pr(>|t|)"]
        }
        
        detailed_results[[paste(fsoc_endpoints[i], clinical_vars[j], sep = "_")]] <- model
        
      }, error = function(e) {
        cat("Error fitting model for", fsoc_endpoints[i], "~", clinical_vars[j], ":", e$message, "\n")
      })
    }
  }
  
  # Use clean labels for output
  colnames(cor_matrix) <- clinical_labels
  colnames(p_matrix) <- clinical_labels
  
  # Print summary
  cat("\nCoefficients (adjusted for age, sex, disease group):\n")
  print(round(cor_matrix, 4))
  cat("\nP-values:\n")
  print(round(p_matrix, 4))
  
  list(
    coefficient_matrix = cor_matrix,
    p_value_matrix = p_matrix,
    detailed_models = detailed_results,
    variables_used = setNames(clinical_vars, clinical_labels)
  )
}

# Clinical heatmap visualization
plot_clinical_heatmap <- function(clinical_results) {
  
  if (is.null(clinical_results)) {
    cat("No clinical results to plot.\n")
    return(NULL)
  }
  
  p_mat <- clinical_results$p_value_matrix
  coef_mat <- clinical_results$coefficient_matrix
  
  # Create labels with coefficients and significance stars
  sig_labels <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
  for (i in 1:nrow(p_mat)) {
    for (j in 1:ncol(p_mat)) {
      if (!is.na(coef_mat[i, j])) {
        coef_val <- round(coef_mat[i, j], 3)
        stars <- ""
        if (!is.na(p_mat[i, j])) {
          if (p_mat[i, j] < 0.001) stars <- "***"
          else if (p_mat[i, j] < 0.01) stars <- "**"
          else if (p_mat[i, j] < 0.05) stars <- "*"
        }
        sig_labels[i, j] <- paste0(coef_val, stars)
      }
    }
  }
  
  # Softer color palette
  pastel_colors <- colorRampPalette(c("#6699CC", "#99CCDD", "#CCE5EE", 
                                      "#F5F5F5", "#FFEEDD", "#FFCC99", 
                                      "#E07766"))(100)
  
  # Clean row labels
  row_labels <- c("Whole Kidney FSOC", "Medullary FSOC")
  rownames(coef_mat) <- row_labels
  rownames(sig_labels) <- row_labels
  
  pheatmap(coef_mat,
           display_numbers = sig_labels,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = pastel_colors,
           main = "FSOC Associations with Clinical Parameters\nAdjusted for age, sex, disease group",
           fontsize = 11,
           fontsize_number = 9,
           number_color = "black",
           border_color = "grey80",
           cellwidth = 60,
           cellheight = 40,
           angle_col = 45,
           na_col = "grey90")
}

# ============================================================================
# 7. IMPROVED VISUALIZATIONS
# ============================================================================

# DXA Heatmap - Softer colors
plot_dxa_heatmap <- function(dxa_results) {
  
  p_mat <- dxa_results$p_value_matrix
  coef_mat <- dxa_results$coefficient_matrix
  
  sig_labels <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
  for (i in 1:nrow(p_mat)) {
    for (j in 1:ncol(p_mat)) {
      coef_val <- round(coef_mat[i, j], 3)
      stars <- ""
      if (!is.na(p_mat[i, j])) {
        if (p_mat[i, j] < 0.001) stars <- "***"
        else if (p_mat[i, j] < 0.01) stars <- "**"
        else if (p_mat[i, j] < 0.05) stars <- "*"
      }
      sig_labels[i, j] <- paste0(coef_val, stars)
    }
  }
  
  pastel_colors <- colorRampPalette(c("#6699CC", "#99CCDD", "#CCE5EE", 
                                      "#F5F5F5", "#FFEEDD", "#FFCC99", 
                                      "#E07766"))(100)
  
  row_labels <- c("Whole Kidney FSOC", "Medullary FSOC")
  col_labels <- c("Fat Mass (kg)", "Trunk Fat (kg)", "Lean Mass (kg)", "Body Fat %")
  
  rownames(coef_mat) <- row_labels
  colnames(coef_mat) <- col_labels
  rownames(sig_labels) <- row_labels
  colnames(sig_labels) <- col_labels
  
  pheatmap(coef_mat,
           display_numbers = sig_labels,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = pastel_colors,
           main = "FSOC Associations with Body Composition (DXA)\nAdjusted for age, sex, disease group",
           fontsize = 11,
           fontsize_number = 9,
           number_color = "black",
           border_color = "grey80",
           cellwidth = 80,
           cellheight = 40,
           angle_col = 45)
}

# FSOC by Group - Unadjusted p-values (simplified)
plot_fsoc_by_group <- function(dat) {
  
  fsoc_long <- dat %>%
    filter(!is.na(group)) %>%
    dplyr::select(record_id, group, sex,
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    pivot_longer(cols = contains("fsoc"),
                 names_to = "fsoc_type",
                 values_to = "fsoc_value") %>%
    filter(!is.na(fsoc_value)) %>%
    mutate(fsoc_type = case_match(fsoc_type,
                                  "whole_kidney_fsoc_abs" ~ "Whole Kidney FSOC",
                                  "medullary_fsoc_abs" ~ "Medullary FSOC"))
  
  # Unadjusted Kruskal-Wallis - calculate separately for subtitle
  kw_whole <- kruskal.test(fsoc_value ~ group, 
                           data = fsoc_long %>% filter(fsoc_type == "Whole Kidney FSOC"))
  kw_med <- kruskal.test(fsoc_value ~ group, 
                         data = fsoc_long %>% filter(fsoc_type == "Medullary FSOC"))
  
  # Format p-values for subtitle
  format_p <- function(p) {
    if (p < 0.001) return("p < 0.001")
    else return(paste0("p = ", round(p, 3)))
  }
  
  subtitle_text <- paste0("Whole Kidney: ", format_p(kw_whole$p.value), 
                          ";  Medullary: ", format_p(kw_med$p.value),
                          " (Kruskal-Wallis)")
  
  # Print pairwise comparisons to console
  cat("\n--- Pairwise Wilcoxon Tests (Unadjusted) ---\n")
  pairwise_results <- fsoc_long %>%
    group_by(fsoc_type) %>%
    wilcox_test(fsoc_value ~ group)
  print(as.data.frame(pairwise_results %>% filter(p < 0.05)))
  
  # Base plot - no geom_text, use subtitle instead
  p <- ggplot(fsoc_long, aes(x = group, y = fsoc_value, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    facet_wrap(~ fsoc_type, scales = "free_y", ncol = 2) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(x = "Disease Group", 
         y = "FSOC Value (s⁻¹)",
         title = "FSOC Endpoints by Disease Group",
         subtitle = subtitle_text) +
    scale_fill_brewer(palette = "Set2")
  
  return(p)
}

# FSOC by Group and Sex - Unadjusted p-values
plot_fsoc_by_group_sex <- function(dat) {
  
  fsoc_long <- dat %>%
    filter(!is.na(group), !is.na(sex)) %>%
    dplyr::select(record_id, group, sex,
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    pivot_longer(cols = contains("fsoc"),
                 names_to = "fsoc_type",
                 values_to = "fsoc_value") %>%
    filter(!is.na(fsoc_value)) %>%
    mutate(fsoc_type = case_match(fsoc_type,
                                  "whole_kidney_fsoc_abs" ~ "Whole Kidney FSOC",
                                  "medullary_fsoc_abs" ~ "Medullary FSOC"))
  
  # Unadjusted Wilcoxon for sex
  sex_stats <- fsoc_long %>%
    group_by(fsoc_type) %>%
    wilcox_test(fsoc_value ~ sex) %>%
    mutate(p_label = case_when(
      p < 0.001 ~ "Sex: p < 0.001",
      p < 0.05 ~ paste0("Sex: p = ", round(p, 3)),
      TRUE ~ paste0("Sex: p = ", round(p, 2))
    ))
  
  # Interaction test
  interaction_stats <- fsoc_long %>%
    group_by(fsoc_type) %>%
    do({
      model <- aov(fsoc_value ~ group * sex, data = .)
      anova_res <- summary(model)[[1]]
      interaction_p <- anova_res["group:sex", "Pr(>F)"]
      data.frame(interaction_p = interaction_p)
    }) %>%
    mutate(interaction_label = case_when(
      interaction_p < 0.001 ~ "Interaction: p < 0.001",
      interaction_p < 0.05 ~ paste0("Interaction: p = ", round(interaction_p, 3)),
      TRUE ~ paste0("Interaction: p = ", round(interaction_p, 2))
    ))
  
  combined_stats <- left_join(sex_stats, interaction_stats, by = "fsoc_type") %>%
    mutate(full_label = paste0(p_label, "\n", interaction_label))
  
  p <- ggplot(fsoc_long, aes(x = group, y = fsoc_value, fill = sex)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8)) +
    geom_point(aes(color = sex), 
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
               alpha = 0.4, size = 1.5) +
    facet_wrap(~ fsoc_type, scales = "free_y", ncol = 2) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(x = "Disease Group", y = "FSOC Value (s⁻¹)",
         title = "FSOC Endpoints by Disease Group and Sex",
         fill = "Sex", color = "Sex",
         caption = "Unadjusted p-values") +
    scale_fill_manual(values = c("Female" = "#E8A0A0", "Male" = "#8CB4D0")) +
    scale_color_manual(values = c("Female" = "#C06060", "Male" = "#4080A0"))
  
  p <- p + geom_text(data = combined_stats,
                     aes(x = 2.5, y = Inf, label = full_label),
                     inherit.aes = FALSE,
                     vjust = 1.3, hjust = 0.5,
                     size = 3, fontface = "italic",
                     lineheight = 0.9)
  
  return(p)
}

# ============================================================================
# 8. DESCRIPTIVE STATISTICS (TABLE 1)
# ============================================================================

create_table1 <- function(dat) {
  dat %>%
    select(group, age, sex, weight,
           whole_kidney_fsoc_abs, medullary_fsoc_abs,
           dexa_fat_kg, dexa_lean_kg) %>%
    tbl_summary(
      by = group,
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 2
    ) %>%
    add_p() %>%
    add_overall() %>%
    modify_header(label ~ "**Variable**") %>%
    bold_labels()
}

# ============================================================================
# 9. EXPORT RESULTS
# ============================================================================

OUTPUT_DIR <- "C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani"

export_results <- function(group_results, sex_results, dxa_results, 
                           clinical_results, negative_results, output_dir = OUTPUT_DIR) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Export group comparison results
  sink(file.path(output_dir, "group_comparisons.txt"))
  cat("========================================\n")
  cat("GROUP COMPARISONS: FSOC BY DISEASE GROUP\n")
  cat("========================================\n\n")
  print(group_results)
  sink()
  
  # Export sex difference results
  if (!is.null(sex_results)) {
    sink(file.path(output_dir, "sex_differences.txt"))
    cat("========================================\n")
    cat("SEX DIFFERENCES IN FSOC\n")
    cat("========================================\n\n")
    print(sex_results)
    sink()
  }
  
  # Export negative FSOC analysis
  sink(file.path(output_dir, "negative_fsoc_analysis.txt"))
  cat("========================================\n")
  cat("NEGATIVE FSOC VALUES ANALYSIS\n")
  cat("========================================\n\n")
  cat("Overall Summary:\n")
  print(as.data.frame(negative_results$overall_summary))
  cat("\nBy Disease Group:\n")
  print(as.data.frame(negative_results$by_group))
  cat("\nBy Sex:\n")
  print(as.data.frame(negative_results$by_sex))
  cat("\nCharacteristics Comparison:\n")
  print(as.data.frame(negative_results$comparison))
  sink()
  
  # Export negative subjects list
  write.csv(negative_results$negative_subjects, 
            file.path(output_dir, "negative_fsoc_subjects.csv"), 
            row.names = FALSE)
  
  # Export DXA heatmap
  pdf(file.path(output_dir, "dxa_heatmap.pdf"), width = 8, height = 6)
  plot_dxa_heatmap(dxa_results)
  dev.off()
  
  tiff(file.path(output_dir, "dxa_heatmap.tiff"), width = 8, height = 6, 
       units = "in", res = 300, compression = "lzw")
  plot_dxa_heatmap(dxa_results)
  dev.off()
  
  # Export clinical heatmap
  if (!is.null(clinical_results)) {
    pdf(file.path(output_dir, "clinical_heatmap.pdf"), width = 10, height = 6)
    plot_clinical_heatmap(clinical_results)
    dev.off()
    
    tiff(file.path(output_dir, "clinical_heatmap.tiff"), width = 10, height = 6, 
         units = "in", res = 300, compression = "lzw")
    plot_clinical_heatmap(clinical_results)
    dev.off()
    
    # Export clinical coefficient matrices
    write.csv(clinical_results$coefficient_matrix, 
              file.path(output_dir, "clinical_coefficients.csv"), row.names = TRUE)
    write.csv(clinical_results$p_value_matrix, 
              file.path(output_dir, "clinical_pvalues.csv"), row.names = TRUE)
  }
  
  # Export FSOC boxplots with p-values
  p1 <- plot_fsoc_by_group(dat)
  ggsave(file.path(output_dir, "fsoc_by_group.pdf"), p1, width = 10, height = 6)
  ggsave(file.path(output_dir, "fsoc_by_group.tiff"), p1, width = 10, height = 6, dpi = 300)
  
  p2 <- plot_fsoc_by_group_sex(dat)
  ggsave(file.path(output_dir, "fsoc_by_group_sex.pdf"), p2, width = 12, height = 6)
  ggsave(file.path(output_dir, "fsoc_by_group_sex.tiff"), p2, width = 12, height = 6, dpi = 300)
  
  # Export negative FSOC distribution plot
  p3 <- plot_negative_fsoc_distribution(dat)
  ggsave(file.path(output_dir, "negative_fsoc_distribution.pdf"), p3, width = 12, height = 10)
  ggsave(file.path(output_dir, "negative_fsoc_distribution.tiff"), p3, width = 12, height = 10, dpi = 300)
  
  # Export DXA coefficient matrices
  write.csv(dxa_results$coefficient_matrix, 
            file.path(output_dir, "dxa_coefficients.csv"), row.names = TRUE)
  write.csv(dxa_results$p_value_matrix, 
            file.path(output_dir, "dxa_pvalues.csv"), row.names = TRUE)
  
  cat("\nResults exported to:", output_dir, "\n")
  cat("\nFiles created:\n")
  cat("- group_comparisons.txt\n")
  if (!is.null(sex_results)) cat("- sex_differences.txt\n")
  cat("- negative_fsoc_analysis.txt\n")
  cat("- negative_fsoc_subjects.csv\n")
  cat("- negative_fsoc_distribution.pdf/.tiff\n")
  cat("- dxa_heatmap.pdf/.tiff\n")
  cat("- dxa_coefficients.csv and dxa_pvalues.csv\n")
  if (!is.null(clinical_results)) {
    cat("- clinical_heatmap.pdf/.tiff\n")
    cat("- clinical_coefficients.csv and clinical_pvalues.csv\n")
  }
  cat("- fsoc_by_group.pdf/.tiff\n")
  cat("- fsoc_by_group_sex.pdf/.tiff\n")
}

# ============================================================================
# 10. MAIN EXECUTION WORKFLOW
# ============================================================================

cat("========================================\n")
cat("FSOC ANALYSIS - STARTING\n")
cat("========================================\n\n")

# Step 1: Calculate averaged FSOC endpoints
dat <- calculate_fsoc_averages(dat)
cat("Step 1: FSOC averages calculated\n")

# Step 2: Analyze negative FSOC values
negative_results <- analyze_negative_fsoc(dat)
cat("\nStep 2: Negative FSOC analysis complete\n")

# Step 3: Run group comparisons (adjusted)
group_results <- run_all_group_comparisons(dat)
cat("\nStep 3: Group comparisons complete\n")

# Step 4: Run sex analyses
sex_results <- tryCatch({
  run_sex_analyses(dat)
}, error = function(e) {
  cat("Sex analysis failed - insufficient variation in sex across groups\n")
  return(NULL)
})
cat("Step 4: Sex analyses complete\n")

# Step 5: Run DXA analyses
dxa_results <- analyze_dxa_associations(dat)
cat("Step 5: DXA analyses complete\n")

# Step 6: Run clinical parameter analyses
clinical_results <- tryCatch({
  analyze_clinical_associations(dat)
}, error = function(e) {
  cat("Clinical analysis failed:", e$message, "\n")
  return(NULL)
})
cat("Step 6: Clinical parameter analyses complete\n")

# Step 7: Generate descriptive table
table1 <- create_table1(dat)
print(table1)
cat("Step 7: Table 1 generated\n")

# Step 8: Save Table 1
table1_df <- as.data.frame(table1)
write.csv(table1_df, file.path(OUTPUT_DIR, "table1_descriptives.csv"))

# Step 9: Create and display plots
p_group <- plot_fsoc_by_group(dat)
print(p_group)

p_sex <- plot_fsoc_by_group_sex(dat)
print(p_sex)

p_negative <- plot_negative_fsoc_distribution(dat)
print(p_negative)

# Step 10: Export all results
export_results(group_results, sex_results, dxa_results, clinical_results, negative_results)

cat("\n========================================\n")
cat("ALL ANALYSES COMPLETE!\n")
cat("Results saved to:", OUTPUT_DIR, "\n")
cat("========================================\n")









# ============================================================================
# FIND SUBJECTS WITH NEGATIVE FSOC VALUES
# ============================================================================

find_negative_fsoc <- function(dat) {
  
  result <- dat %>%
    filter(!is.na(group)) %>%
    filter(whole_kidney_fsoc_abs < 0 | medullary_fsoc_abs < 0) %>%
    dplyr::select(record_id, group, sex, age, 
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    arrange(whole_kidney_fsoc_abs)
  
  cat("Found", nrow(result), "subjects with negative FSOC values:\n\n")
  print(as.data.frame(result))
  
  return(invisible(result))
}

 negative_subjects <- find_negative_fsoc(dat)