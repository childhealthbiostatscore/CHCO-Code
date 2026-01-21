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

# FSOC by Group - Unadjusted p-values (excludes negative FSOC values)
plot_fsoc_by_group <- function(dat) {
  
  fsoc_long <- dat %>%
    filter(!is.na(group)) %>%
    dplyr::select(record_id, group, sex,
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    pivot_longer(cols = contains("fsoc"),
                 names_to = "fsoc_type",
                 values_to = "fsoc_value") %>%
    filter(!is.na(fsoc_value)) %>%
    # REMOVE NEGATIVE FSOC VALUES
    filter(fsoc_value >= 0) %>%
    mutate(fsoc_type = case_match(fsoc_type,
                                  "whole_kidney_fsoc_abs" ~ "Whole Kidney FSOC",
                                  "medullary_fsoc_abs" ~ "Medullary FSOC"))
  
  # Count excluded
  n_excluded <- dat %>%
    filter(!is.na(group)) %>%
    summarise(
      n_wk_neg = sum(whole_kidney_fsoc_abs < 0, na.rm = TRUE),
      n_med_neg = sum(medullary_fsoc_abs < 0, na.rm = TRUE)
    )
  cat("\nExcluded negative FSOC values: WK =", n_excluded$n_wk_neg, 
      ", Medullary =", n_excluded$n_med_neg, "\n")
  
  # Unadjusted Kruskal-Wallis
  kw_whole <- kruskal.test(fsoc_value ~ group, 
                           data = fsoc_long %>% filter(fsoc_type == "Whole Kidney FSOC"))
  kw_med <- kruskal.test(fsoc_value ~ group, 
                         data = fsoc_long %>% filter(fsoc_type == "Medullary FSOC"))
  
  format_p <- function(p) {
    if (p < 0.001) return("p < 0.001")
    else return(paste0("p = ", round(p, 3)))
  }
  
  subtitle_text <- paste0("Whole Kidney: ", format_p(kw_whole$p.value), 
                          ";  Medullary: ", format_p(kw_med$p.value),
                          " (Kruskal-Wallis)")
  
  # Pairwise Wilcoxon tests - RAW p-values
  cat("\n--- Pairwise Wilcoxon Tests (Raw/Unadjusted p-values) ---\n")
  cat("Negative FSOC values excluded from analysis\n\n")
  
  pairwise_results <- fsoc_long %>%
    group_by(fsoc_type) %>%
    wilcox_test(fsoc_value ~ group, p.adjust.method = "none") %>%
    mutate(p_formatted = case_when(
      p < 0.001 ~ "<0.001",
      p < 0.01 ~ as.character(round(p, 4)),
      TRUE ~ as.character(round(p, 3))
    ))
  
  # Print all pairwise comparisons
  cat("WHOLE KIDNEY FSOC:\n")
  wk_pairs <- pairwise_results %>% 
    filter(fsoc_type == "Whole Kidney FSOC") %>%
    dplyr::select(group1, group2, p, p_formatted)
  print(as.data.frame(wk_pairs))
  
  cat("\nMEDULLARY FSOC:\n")
  med_pairs <- pairwise_results %>% 
    filter(fsoc_type == "Medullary FSOC") %>%
    dplyr::select(group1, group2, p, p_formatted)
  print(as.data.frame(med_pairs))
  
  # Base plot
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
         subtitle = subtitle_text,
         caption = "Negative FSOC values excluded") +
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






# Install ggsignif if needed
# install.packages("ggsignif")
library(ggsignif)
library(patchwork)

# FSOC by Group with ALL pairwise p-values on plot
plot_fsoc_by_group <- function(dat) {
  
  fsoc_long <- dat %>%
    filter(!is.na(group)) %>%
    dplyr::select(record_id, group, sex,
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    filter(whole_kidney_fsoc_abs < 15) %>% 
    pivot_longer(cols = contains("fsoc"),
                 names_to = "fsoc_type",
                 values_to = "fsoc_value") %>%
    filter(!is.na(fsoc_value)) %>%
    filter(fsoc_value >= 0) %>%
    mutate(fsoc_type = case_match(fsoc_type,
                                  "whole_kidney_fsoc_abs" ~ "Whole Kidney FSOC",
                                  "medullary_fsoc_abs" ~ "Medullary FSOC"))
  
  # Get unique groups (excluding PKD if no data)
  groups <- fsoc_long %>% 
    filter(fsoc_type == "Whole Kidney FSOC") %>%
    pull(group) %>% 
    unique() %>%
    as.character()
  
  # Create all pairwise comparisons
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  # Calculate pairwise p-values
  calc_pairwise_p <- function(data, fsoc_name) {
    subset_data <- data %>% filter(fsoc_type == fsoc_name)
    p_values <- sapply(comparisons, function(pair) {
      g1 <- subset_data %>% filter(group == pair[1]) %>% pull(fsoc_value)
      g2 <- subset_data %>% filter(group == pair[2]) %>% pull(fsoc_value)
      if (length(g1) > 1 & length(g2) > 1) {
        wilcox.test(g1, g2)$p.value
      } else {
        NA
      }
    })
    return(p_values)
  }
  
  wk_pvals <- calc_pairwise_p(fsoc_long, "Whole Kidney FSOC")
  med_pvals <- calc_pairwise_p(fsoc_long, "Medullary FSOC")
  
  # Format p-values for plot
  format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("p<0.001")
    if (p < 0.01) return(sprintf("p=%.3f", p))
    return(sprintf("p=%.2f", p))
  }
  
  # Print to console
  cat("\n--- Pairwise Wilcoxon Tests (Raw p-values) ---\n")
  cat("WHOLE KIDNEY:\n")
  for (i in seq_along(comparisons)) {
    cat(sprintf("  %s vs %s: %s\n", comparisons[[i]][1], comparisons[[i]][2], format_p(wk_pvals[i])))
  }
  cat("\nMEDULLARY:\n")
  for (i in seq_along(comparisons)) {
    cat(sprintf("  %s vs %s: %s\n", comparisons[[i]][1], comparisons[[i]][2], format_p(med_pvals[i])))
  }
  
  # Whole Kidney plot
  wk_data <- fsoc_long %>% filter(fsoc_type == "Whole Kidney FSOC")
  
  p_wk <- ggplot(wk_data, aes(x = group, y = fsoc_value, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    geom_signif(
      comparisons = comparisons,
      annotations = sapply(wk_pvals, format_p),
      step_increase = 0.12,
      tip_length = 0.02,
      textsize = 2.8,
      vjust = 0.3,
      map_signif_level = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    ) +
    labs(x = "Disease Group", 
         y = "FSOC Value (s⁻¹)",
         title = "Whole Kidney FSOC") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  # Medullary plot
  med_data <- fsoc_long %>% filter(fsoc_type == "Medullary FSOC")
  
  p_med <- ggplot(med_data, aes(x = group, y = fsoc_value, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    geom_signif(
      comparisons = comparisons,
      annotations = sapply(med_pvals, format_p),
      step_increase = 0.12,
      tip_length = 0.02,
      textsize = 2.8,
      vjust = 0.3,
      map_signif_level = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    ) +
    labs(x = "Disease Group", 
         y = "FSOC Value (s⁻¹)",
         title = "Medullary FSOC") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  # Combine
  combined <- p_wk + p_med +
    plot_annotation(
      title = "FSOC Endpoints by Disease Group",
      subtitle = "Raw pairwise Wilcoxon p-values; negative FSOC values excluded",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, face = "italic", hjust = 0.5)
      )
    )
  
  return(combined)
}

# Run it
p_group <- plot_fsoc_by_group(dat)
print(p_group)

# Save
ggsave(file.path(OUTPUT_DIR, "fsoc_by_group_pairwise.pdf"), p_group, width = 14, height = 8)
ggsave(file.path(OUTPUT_DIR, "fsoc_by_group_pairwise.tiff"), p_group, width = 14, height = 8, dpi = 300)







# FSOC by Group with Sex side-by-side and pairwise p-values
library(ggsignif)
library(patchwork)

plot_fsoc_by_group_sex <- function(dat) {
  
  # Simple filtering
  dat_clean <- dat %>%
    filter(!is.na(group)) %>%
    filter(!is.na(sex)) %>%
    filter(whole_kidney_fsoc_abs >= 0 | is.na(whole_kidney_fsoc_abs)) %>%
    filter(whole_kidney_fsoc_abs <= 15 | is.na(whole_kidney_fsoc_abs)) %>%
    filter(medullary_fsoc_abs >= 0 | is.na(medullary_fsoc_abs))
  
  cat("Excluded", nrow(dat) - nrow(dat_clean), "subjects with invalid FSOC or missing sex\n")
  
  fsoc_long <- dat_clean %>%
    dplyr::select(record_id, group, sex,
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    pivot_longer(cols = contains("fsoc"),
                 names_to = "fsoc_type",
                 values_to = "fsoc_value") %>%
    filter(!is.na(fsoc_value)) %>%
    mutate(fsoc_type = case_match(fsoc_type,
                                  "whole_kidney_fsoc_abs" ~ "Whole Kidney FSOC",
                                  "medullary_fsoc_abs" ~ "Medullary FSOC"))
  
  # Get unique groups
  groups <- fsoc_long %>% 
    filter(fsoc_type == "Whole Kidney FSOC") %>%
    pull(group) %>% 
    unique() %>%
    as.character()
  
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  # Calculate pairwise p-values (comparing groups, ignoring sex)
  calc_pairwise_p <- function(data, fsoc_name) {
    subset_data <- data %>% filter(fsoc_type == fsoc_name)
    p_values <- sapply(comparisons, function(pair) {
      g1 <- subset_data %>% filter(group == pair[1]) %>% pull(fsoc_value)
      g2 <- subset_data %>% filter(group == pair[2]) %>% pull(fsoc_value)
      if (length(g1) > 1 & length(g2) > 1) {
        wilcox.test(g1, g2)$p.value
      } else {
        NA
      }
    })
    return(p_values)
  }
  
  wk_pvals <- calc_pairwise_p(fsoc_long, "Whole Kidney FSOC")
  med_pvals <- calc_pairwise_p(fsoc_long, "Medullary FSOC")
  
  # Calculate sex p-value within each group
  calc_sex_p <- function(data, fsoc_name) {
    subset_data <- data %>% filter(fsoc_type == fsoc_name)
    sex_pvals <- subset_data %>%
      group_by(group) %>%
      summarise(
        p = tryCatch(wilcox.test(fsoc_value ~ sex)$p.value, error = function(e) NA),
        .groups = "drop"
      )
    return(sex_pvals)
  }
  
  wk_sex_p <- calc_sex_p(fsoc_long, "Whole Kidney FSOC")
  med_sex_p <- calc_sex_p(fsoc_long, "Medullary FSOC")
  
  format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("p<0.001")
    if (p < 0.01) return(sprintf("p=%.3f", p))
    return(sprintf("p=%.2f", p))
  }
  
  # Print to console
  cat("\n--- Pairwise Wilcoxon Tests Between Groups (Raw p-values) ---\n")
  cat("WHOLE KIDNEY:\n")
  for (i in seq_along(comparisons)) {
    cat(sprintf("  %s vs %s: %s\n", comparisons[[i]][1], comparisons[[i]][2], format_p(wk_pvals[i])))
  }
  cat("\nMEDULLARY:\n")
  for (i in seq_along(comparisons)) {
    cat(sprintf("  %s vs %s: %s\n", comparisons[[i]][1], comparisons[[i]][2], format_p(med_pvals[i])))
  }
  
  cat("\n--- Sex Differences Within Each Group (Wilcoxon) ---\n")
  cat("WHOLE KIDNEY:\n")
  for (i in 1:nrow(wk_sex_p)) {
    cat(sprintf("  %s: %s\n", wk_sex_p$group[i], format_p(wk_sex_p$p[i])))
  }
  cat("\nMEDULLARY:\n")
  for (i in 1:nrow(med_sex_p)) {
    cat(sprintf("  %s: %s\n", med_sex_p$group[i], format_p(med_sex_p$p[i])))
  }
  
  # Whole Kidney plot
  wk_data <- fsoc_long %>% filter(fsoc_type == "Whole Kidney FSOC")
  
  p_wk <- ggplot(wk_data, aes(x = group, y = fsoc_value, fill = sex)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8)) +
    geom_point(aes(color = sex),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               alpha = 0.4, size = 1.5) +
    geom_signif(
      comparisons = comparisons,
      annotations = sapply(wk_pvals, format_p),
      step_increase = 0.12,
      tip_length = 0.02,
      textsize = 2.8,
      vjust = 0.3,
      map_signif_level = FALSE,
      y_position = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    ) +
    labs(x = "Disease Group", 
         y = "FSOC Value (s⁻¹)",
         title = "Whole Kidney FSOC",
         fill = "Sex", color = "Sex") +
    scale_fill_manual(values = c("Female" = "#E8A0A0", "Male" = "#8CB4D0")) +
    scale_color_manual(values = c("Female" = "#C06060", "Male" = "#4080A0")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  # Medullary plot
  med_data <- fsoc_long %>% filter(fsoc_type == "Medullary FSOC")
  
  p_med <- ggplot(med_data, aes(x = group, y = fsoc_value, fill = sex)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8)) +
    geom_point(aes(color = sex),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               alpha = 0.4, size = 1.5) +
    geom_signif(
      comparisons = comparisons,
      annotations = sapply(med_pvals, format_p),
      step_increase = 0.12,
      tip_length = 0.02,
      textsize = 2.8,
      vjust = 0.3,
      map_signif_level = FALSE,
      y_position = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    ) +
    labs(x = "Disease Group", 
         y = "FSOC Value (s⁻¹)",
         title = "Medullary FSOC",
         fill = "Sex", color = "Sex") +
    scale_fill_manual(values = c("Female" = "#E8A0A0", "Male" = "#8CB4D0")) +
    scale_color_manual(values = c("Female" = "#C06060", "Male" = "#4080A0")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  # Combine with shared legend
  combined <- p_wk + p_med +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "FSOC Endpoints by Disease Group and Sex",
      subtitle = "Raw pairwise Wilcoxon p-values (between groups); negative & WK FSOC >15 excluded",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, face = "italic", hjust = 0.5),
        legend.position = "bottom"
      )
    )
  
  return(combined)
}

# Run it
p_sex <- plot_fsoc_by_group_sex(dat)
print(p_sex)

# Save
ggsave(file.path(OUTPUT_DIR, "fsoc_by_group_sex_pairwise.pdf"), p_sex, width = 14, height = 8)
ggsave(file.path(OUTPUT_DIR, "fsoc_by_group_sex_pairwise.tiff"), p_sex, width = 14, height = 8, dpi = 300)



# ============================================================================
# UNADJUSTED CORRELATION HEATMAPS BY GROUP
# ============================================================================

library(pheatmap)
library(tidyverse)

# ============================================================================
# 1. DXA CORRELATIONS BY GROUP (UNADJUSTED)
# ============================================================================

create_dxa_heatmap_by_group <- function(dat) {
  
  # Filter invalid FSOC
  dat_clean <- dat %>%
    filter(!is.na(group)) %>%
    filter(whole_kidney_fsoc_abs >= 0 | is.na(whole_kidney_fsoc_abs)) %>%
    filter(whole_kidney_fsoc_abs <= 15 | is.na(whole_kidney_fsoc_abs)) %>%
    filter(medullary_fsoc_abs >= 0 | is.na(medullary_fsoc_abs))
  
  # Define variables
  fsoc_vars <- c("whole_kidney_fsoc_abs", "medullary_fsoc_abs")
  dxa_vars <- c("dexa_fat_kg", "dexa_trunk_kg", "dexa_lean_kg", "dexa_body_fat")
  dxa_labels <- c("Fat Mass", "Trunk Fat", "Lean Mass", "Body Fat %")
  
  groups <- unique(dat_clean$group[!is.na(dat_clean$group)])
  
  # Create matrices for correlations and p-values
  # Rows: Group_WK, Group_Med for each group
  # Cols: DXA variables
  
  row_names <- c()
  for (g in groups) {
    row_names <- c(row_names, paste0(g, " - WK"), paste0(g, " - Med"))
  }
  
  cor_matrix <- matrix(NA, nrow = length(row_names), ncol = length(dxa_vars))
  p_matrix <- matrix(NA, nrow = length(row_names), ncol = length(dxa_vars))
  rownames(cor_matrix) <- row_names
  colnames(cor_matrix) <- dxa_labels
  rownames(p_matrix) <- row_names
  colnames(p_matrix) <- dxa_labels
  
  # Calculate correlations for each group
  row_idx <- 1
  for (g in groups) {
    group_data <- dat_clean %>% filter(group == g)
    
    for (j in seq_along(dxa_vars)) {
      # Whole kidney correlation
      wk_data <- group_data %>% 
        filter(!is.na(whole_kidney_fsoc_abs), !is.na(.data[[dxa_vars[j]]]))
      if (nrow(wk_data) >= 5) {
        cor_test <- cor.test(wk_data$whole_kidney_fsoc_abs, wk_data[[dxa_vars[j]]], method = "spearman")
        cor_matrix[row_idx, j] <- cor_test$estimate
        p_matrix[row_idx, j] <- cor_test$p.value
      }
      
      # Medullary correlation
      med_data <- group_data %>% 
        filter(!is.na(medullary_fsoc_abs), !is.na(.data[[dxa_vars[j]]]))
      if (nrow(med_data) >= 5) {
        cor_test <- cor.test(med_data$medullary_fsoc_abs, med_data[[dxa_vars[j]]], method = "spearman")
        cor_matrix[row_idx + 1, j] <- cor_test$estimate
        p_matrix[row_idx + 1, j] <- cor_test$p.value
      }
    }
    row_idx <- row_idx + 2
  }
  
  # Create significance labels
  sig_labels <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  for (i in 1:nrow(p_matrix)) {
    for (j in 1:ncol(p_matrix)) {
      if (!is.na(cor_matrix[i, j])) {
        r_val <- sprintf("%.2f", cor_matrix[i, j])
        stars <- ""
        if (!is.na(p_matrix[i, j])) {
          if (p_matrix[i, j] < 0.001) stars <- "***"
          else if (p_matrix[i, j] < 0.01) stars <- "**"
          else if (p_matrix[i, j] < 0.05) stars <- "*"
        }
        sig_labels[i, j] <- paste0(r_val, stars)
      }
    }
  }
  
  # Row annotation for group
  row_annotation <- data.frame(
    Group = rep(groups, each = 2),
    FSOC = rep(c("Whole Kidney", "Medulla"), length(groups))
  )
  rownames(row_annotation) <- row_names
  
  # Color palette
  pastel_colors <- colorRampPalette(c("#4575B4", "#91BFDB", "#E0F3F8", 
                                      "#FFFFBF", "#FEE090", "#FC8D59", 
                                      "#D73027"))(100)
  
  # Print correlations to console
  cat("\n========================================\n")
  cat("DXA CORRELATIONS BY GROUP (Spearman, Unadjusted)\n")
  cat("========================================\n")
  cat("* p<0.05, ** p<0.01, *** p<0.001\n\n")
  
  for (g in groups) {
    cat(paste0("\n", g, ":\n"))
    idx_wk <- which(row_names == paste0(g, " - WK"))
    idx_med <- which(row_names == paste0(g, " - Med"))
    cat("  Whole Kidney: ")
    cat(paste(dxa_labels, "=", sig_labels[idx_wk, ], collapse = ", "), "\n")
    cat("  Medulla:      ")
    cat(paste(dxa_labels, "=", sig_labels[idx_med, ], collapse = ", "), "\n")
  }
  
  # Create heatmap
  p <- pheatmap(cor_matrix,
                display_numbers = sig_labels,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                color = pastel_colors,
                breaks = seq(-1, 1, length.out = 101),
                main = "FSOC-DXA Correlations by Group\n(Spearman, Unadjusted)",
                fontsize = 10,
                fontsize_number = 8,
                number_color = "black",
                border_color = "grey80",
                cellwidth = 60,
                cellheight = 25,
                angle_col = 45,
                gaps_row = seq(2, length(row_names) - 2, by = 2),
                annotation_row = row_annotation,
                na_col = "grey90")
  
  return(list(cor_matrix = cor_matrix, p_matrix = p_matrix, plot = p))
}

# ============================================================================
# 2. CLINICAL CORRELATIONS BY GROUP (UNADJUSTED)
# ============================================================================

create_clinical_heatmap_by_group <- function(dat) {
  
  # Filter invalid FSOC
  dat_clean <- dat %>%
    filter(!is.na(group)) %>%
    filter(whole_kidney_fsoc_abs >= 0 | is.na(whole_kidney_fsoc_abs)) %>%
    filter(whole_kidney_fsoc_abs <= 15 | is.na(whole_kidney_fsoc_abs)) %>%
    filter(medullary_fsoc_abs >= 0 | is.na(medullary_fsoc_abs))
  
  # Find available clinical variables
  available_vars <- names(dat_clean)
  
  clinical_var_map <- list(
    "eGFR" = c("eGFR_CKD_epi", "gfr_fas_cr", "eGFR", "egfr"),
    "mGFR" = c("gfr_raw_plasma", "gfr_bsa_plasma", "mGFR"),
    "HbA1c" = c("hba1c", "HbA1c", "a1c"),
    "SBP" = c("sbp", "SBP", "systolic_bp"),
    "DBP" = c("dbp", "DBP", "diastolic_bp"),
    "BMI" = c("bmi", "BMI"),
    "log(UACR)" = c("acr_u", "uacr", "ACR")
  )
  
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
  
  cat("\nClinical variables found:", paste(clinical_labels, collapse = ", "), "\n")
  
  if (length(clinical_vars) == 0) {
    cat("No clinical variables found.\n")
    return(NULL)
  }
  
  groups <- unique(dat_clean$group[!is.na(dat_clean$group)])
  
  # Create row names
  row_names <- c()
  for (g in groups) {
    row_names <- c(row_names, paste0(g, " - WK"), paste0(g, " - Med"))
  }
  
  cor_matrix <- matrix(NA, nrow = length(row_names), ncol = length(clinical_vars))
  p_matrix <- matrix(NA, nrow = length(row_names), ncol = length(clinical_vars))
  rownames(cor_matrix) <- row_names
  colnames(cor_matrix) <- clinical_labels
  rownames(p_matrix) <- row_names
  colnames(p_matrix) <- clinical_labels
  
  # Calculate correlations
  row_idx <- 1
  for (g in groups) {
    group_data <- dat_clean %>% filter(group == g)
    
    for (j in seq_along(clinical_vars)) {
      # Log-transform ACR if needed
      if (clinical_labels[j] == "log(ACR)") {
        group_data <- group_data %>%
          mutate(temp_clinical = log(.data[[clinical_vars[j]]] + 1))
      } else {
        group_data <- group_data %>%
          mutate(temp_clinical = .data[[clinical_vars[j]]])
      }
      
      # Whole kidney correlation
      wk_data <- group_data %>% 
        filter(!is.na(whole_kidney_fsoc_abs), !is.na(temp_clinical))
      if (nrow(wk_data) >= 5) {
        cor_test <- cor.test(wk_data$whole_kidney_fsoc_abs, wk_data$temp_clinical, method = "spearman")
        cor_matrix[row_idx, j] <- cor_test$estimate
        p_matrix[row_idx, j] <- cor_test$p.value
      }
      
      # Medullary correlation
      med_data <- group_data %>% 
        filter(!is.na(medullary_fsoc_abs), !is.na(temp_clinical))
      if (nrow(med_data) >= 5) {
        cor_test <- cor.test(med_data$medullary_fsoc_abs, med_data$temp_clinical, method = "spearman")
        cor_matrix[row_idx + 1, j] <- cor_test$estimate
        p_matrix[row_idx + 1, j] <- cor_test$p.value
      }
    }
    row_idx <- row_idx + 2
  }
  
  # Create significance labels
  sig_labels <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  for (i in 1:nrow(p_matrix)) {
    for (j in 1:ncol(p_matrix)) {
      if (!is.na(cor_matrix[i, j])) {
        r_val <- sprintf("%.2f", cor_matrix[i, j])
        stars <- ""
        if (!is.na(p_matrix[i, j])) {
          if (p_matrix[i, j] < 0.001) stars <- "***"
          else if (p_matrix[i, j] < 0.01) stars <- "**"
          else if (p_matrix[i, j] < 0.05) stars <- "*"
        }
        sig_labels[i, j] <- paste0(r_val, stars)
      }
    }
  }
  
  # Row annotation
  row_annotation <- data.frame(
    Group = rep(groups, each = 2),
    FSOC = rep(c("Whole Kidney", "Medulla"), length(groups))
  )
  rownames(row_annotation) <- row_names
  
  # Color palette
  pastel_colors <- colorRampPalette(c("#4575B4", "#91BFDB", "#E0F3F8", 
                                      "#FFFFBF", "#FEE090", "#FC8D59", 
                                      "#D73027"))(100)
  
  # Print to console
  cat("\n========================================\n")
  cat("CLINICAL CORRELATIONS BY GROUP (Spearman, Unadjusted)\n")
  cat("========================================\n")
  cat("* p<0.05, ** p<0.01, *** p<0.001\n")
  cat("Note: ACR is log-transformed\n\n")
  
  for (g in groups) {
    cat(paste0("\n", g, ":\n"))
    idx_wk <- which(row_names == paste0(g, " - WK"))
    idx_med <- which(row_names == paste0(g, " - Med"))
    cat("  Whole Kidney: ")
    cat(paste(clinical_labels, "=", sig_labels[idx_wk, ], collapse = ", "), "\n")
    cat("  Medulla:      ")
    cat(paste(clinical_labels, "=", sig_labels[idx_med, ], collapse = ", "), "\n")
  }
  
  # Create heatmap
  p <- pheatmap(cor_matrix,
                display_numbers = sig_labels,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                color = pastel_colors,
                breaks = seq(-1, 1, length.out = 101),
                main = "FSOC-Clinical Correlations by Group\n(Spearman, Unadjusted)",
                fontsize = 10,
                fontsize_number = 8,
                number_color = "black",
                border_color = "grey80",
                cellwidth = 50,
                cellheight = 25,
                angle_col = 45,
                gaps_row = seq(2, length(row_names) - 2, by = 2),
                annotation_row = row_annotation,
                na_col = "grey90")
  
  return(list(cor_matrix = cor_matrix, p_matrix = p_matrix, plot = p))
}

# ============================================================================
# RUN AND SAVE
# ============================================================================

# DXA heatmap by group
dxa_by_group <- create_dxa_heatmap_by_group(dat)

pdf(file.path(OUTPUT_DIR, "dxa_heatmap_by_group.pdf"), width = 10, height = 10)
create_dxa_heatmap_by_group(dat)
dev.off()

# Clinical heatmap by group
clinical_by_group <- create_clinical_heatmap_by_group(dat)

pdf(file.path(OUTPUT_DIR, "clinical_heatmap_by_group.pdf"), width = 12, height = 10)
create_clinical_heatmap_by_group(dat)
dev.off()

# Export correlation matrices
write.csv(dxa_by_group$cor_matrix, file.path(OUTPUT_DIR, "dxa_correlations_by_group.csv"))
write.csv(dxa_by_group$p_matrix, file.path(OUTPUT_DIR, "dxa_pvalues_by_group.csv"))

if (!is.null(clinical_by_group)) {
  write.csv(clinical_by_group$cor_matrix, file.path(OUTPUT_DIR, "clinical_correlations_by_group.csv"))
  write.csv(clinical_by_group$p_matrix, file.path(OUTPUT_DIR, "clinical_pvalues_by_group.csv"))
}

cat("\n\nHeatmaps saved to:", OUTPUT_DIR, "\n")


























find_high_fsoc <- function(dat, wk_threshold = 15, med_threshold = 30) {
  
  cat("\n========================================\n")
  cat("HIGH FSOC VALUES\n")
  cat("Thresholds: Whole Kidney >", wk_threshold, ", Medullary >", med_threshold, "\n")
  cat("========================================\n\n")
  
  # Find high whole kidney FSOC
  high_wk <- dat %>%
    filter(!is.na(group)) %>%
    filter(whole_kidney_fsoc_abs > wk_threshold) %>%
    dplyr::select(record_id, group, sex, age, 
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    arrange(desc(whole_kidney_fsoc_abs))
  
  cat("WHOLE KIDNEY FSOC >", wk_threshold, ":\n")
  cat("Found", nrow(high_wk), "subjects\n")
  if (nrow(high_wk) > 0) {
    print(as.data.frame(high_wk))
  }
  
  # Find high medullary FSOC
  high_med <- dat %>%
    filter(!is.na(group)) %>%
    filter(medullary_fsoc_abs > med_threshold) %>%
    dplyr::select(record_id, group, sex, age, 
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    arrange(desc(medullary_fsoc_abs))
  
  cat("\nMEDULLARY FSOC >", med_threshold, ":\n")
  cat("Found", nrow(high_med), "subjects\n")
  if (nrow(high_med) > 0) {
    print(as.data.frame(high_med))
  }
  
  # Combined: any extreme value
  high_any <- dat %>%
    filter(!is.na(group)) %>%
    filter(whole_kidney_fsoc_abs > wk_threshold | medullary_fsoc_abs > med_threshold) %>%
    dplyr::select(record_id, group, sex, age, 
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    arrange(desc(whole_kidney_fsoc_abs))
  
  cat("\n--- Summary ---\n")
  cat("Total subjects with any high FSOC:", nrow(high_any), "\n")
  
  result <- list(
    high_whole_kidney = high_wk,
    high_medullary = high_med,
    high_any = high_any
  )
  
  return(invisible(result))
}

find_high_fsoc(dat)


















#### Updated analyses


# ============================================================================
# 11. VARIABLE IDENTIFICATION
# ============================================================================
# FSOC = Furosemide Stress Oxygen Consumption (from BOLD MRI R2*)
#   Variables: fsoc_l_kidney, fsoc_r_kidney, fsoc_l_medulla, fsoc_r_medulla
#   Already calculated as: pre-furosemide R2* - post-furosemide R2*
#
# K2 = TCA cycle metabolism (from C-11 acetate PET) - if available
#
# BOLD R2* baseline: bold_l_bl_medulla, bold_r_bl_medulla, etc.
#   (Higher R2* = lower oxygenation / more deoxyhemoglobin)
# ============================================================================

identify_all_variables <- function(dat) {
  available_vars <- names(dat)
  
  cat("\n========================================\n")
  cat("VARIABLE IDENTIFICATION FOR PHENOTYPING\n")
  cat("========================================\n\n")
  
  # K2 (PET-derived TCA metabolism)
  k2_vars <- available_vars[grepl("^k2|_k2|cortical_k2|medullary_k2", 
                                  available_vars, ignore.case = TRUE)]
  
  cat("--- K2 (PET TCA metabolism) variables ---\n")
  if (length(k2_vars) > 0) {
    for (v in k2_vars) {
      n_nm <- sum(!is.na(dat[[v]]))
      if (n_nm > 0) {
        cat(sprintf("  %s: n=%d, range=[%.2f, %.2f]\n", 
                    v, n_nm, min(dat[[v]], na.rm=T), max(dat[[v]], na.rm=T)))
      }
    }
  } else {
    cat("  No K2 variables found.\n")
  }
  
  # BOLD R2* baseline
  bold_vars <- available_vars[grepl("bold_.*_bl_|bold.*baseline", 
                                    available_vars, ignore.case = TRUE)]
  
  cat("\n--- BOLD R2* baseline variables ---\n")
  if (length(bold_vars) > 0) {
    for (v in bold_vars) {
      n_nm <- sum(!is.na(dat[[v]]))
      if (n_nm > 0) {
        cat(sprintf("  %s: n=%d, range=[%.2f, %.2f]\n", 
                    v, n_nm, min(dat[[v]], na.rm=T), max(dat[[v]], na.rm=T)))
      }
    }
  } else {
    cat("  No BOLD baseline variables found.\n")
  }
  
  # FSOC variables
  fsoc_vars <- available_vars[grepl("fsoc_", available_vars, ignore.case = TRUE)]
  
  cat("\n--- FSOC variables ---\n")
  if (length(fsoc_vars) > 0) {
    for (v in fsoc_vars) {
      n_nm <- sum(!is.na(dat[[v]]))
      if (n_nm > 0) {
        cat(sprintf("  %s: n=%d, range=[%.2f, %.2f]\n", 
                    v, n_nm, min(dat[[v]], na.rm=T), max(dat[[v]], na.rm=T)))
      }
    }
  } else {
    cat("  No FSOC variables found.\n")
  }
  
  return(list(k2 = k2_vars, bold = bold_vars, fsoc = fsoc_vars))
}

# ============================================================================
# Calculate averaged BOLD baseline R2* if needed
# ============================================================================

calculate_bold_averages <- function(dat) {
  dat <- dat %>%
    mutate(
      medullary_r2star_bl = case_when(
        !is.na(bold_l_bl_medulla) & !is.na(bold_r_bl_medulla) ~ 
          (bold_l_bl_medulla + bold_r_bl_medulla) / 2,
        !is.na(bold_l_bl_medulla) ~ bold_l_bl_medulla,
        !is.na(bold_r_bl_medulla) ~ bold_r_bl_medulla,
        TRUE ~ NA_real_
      ),
      cortical_r2star_bl = case_when(
        !is.na(bold_l_bl_cortex) & !is.na(bold_r_bl_cortex) ~ 
          (bold_l_bl_cortex + bold_r_bl_cortex) / 2,
        !is.na(bold_l_bl_cortex) ~ bold_l_bl_cortex,
        !is.na(bold_r_bl_cortex) ~ bold_r_bl_cortex,
        TRUE ~ NA_real_
      )
    )
  return(dat)
}

# ============================================================================
# 12. PHENOTYPE CLASSIFICATION (High K2/R2*, Low FSOC)
# ============================================================================

create_fsoc_phenotypes <- function(dat, 
                                   secondary_var = NULL,
                                   fsoc_var = "medullary_fsoc_abs",
                                   method = "control_reference",
                                   control_group = "Lean Control",
                                   n_sd = 1,
                                   fsoc_cutoff = NULL,
                                   secondary_cutoff = NULL) {
  
  cat("\n========================================\n")
  cat("FSOC PHENOTYPE CLASSIFICATION\n")
  cat("========================================\n\n")
  
  # Auto-detect secondary variable if not specified
  if (is.null(secondary_var)) {
    available_vars <- names(dat)
    # Try K2 first, then BOLD baseline R2*
    k2_candidates <- c("cortical_k2", "medullary_k2", "k2_cortex", "k2_medulla")
    bold_candidates <- c("medullary_r2star_bl", "cortical_r2star_bl")
    
    secondary_var <- k2_candidates[k2_candidates %in% available_vars][1]
    if (is.na(secondary_var)) {
      secondary_var <- bold_candidates[bold_candidates %in% available_vars][1]
    }
    
    if (is.na(secondary_var)) {
      cat("ERROR: No K2 or BOLD baseline variable found.\n")
      cat("Please specify secondary_var parameter.\n")
      return(NULL)
    }
  }
  
  cat("Using FSOC variable:", fsoc_var, "\n")
  cat("Using secondary variable:", secondary_var, "\n")
  cat("Classification method:", method, "\n")
  cat("Control group:", control_group, "\n\n")
  
  # Filter to subjects with both measurements
  dat_pheno <- dat %>%
    filter(!is.na(group)) %>%
    filter(!is.na(.data[[fsoc_var]])) %>%
    filter(!is.na(.data[[secondary_var]])) %>%
    filter(.data[[fsoc_var]] >= 0)  # Exclude negative FSOC
  
  cat("Subjects with both measurements:", nrow(dat_pheno), "\n")
  
  # Check if control group exists
  if (!control_group %in% unique(dat_pheno$group)) {
    # Try alternative names
    possible_controls <- c("Lean Control", "LC", "Lean_Control", "Control", "Healthy Control")
    control_group <- possible_controls[possible_controls %in% unique(dat_pheno$group)][1]
    if (is.na(control_group)) {
      cat("WARNING: No control group found. Falling back to median method.\n")
      method <- "median"
    } else {
      cat("Using control group:", control_group, "\n")
    }
  }
  
  # Determine cutoffs based on method
  if (method == "control_reference") {
    # Get control group data
    control_data <- dat_pheno %>% filter(group == control_group)
    
    cat("\n--- Control Group Reference (", control_group, ") ---\n")
    cat("N in control group:", nrow(control_data), "\n")
    
    # FSOC cutoff: Control mean - n_sd * SD
    # Low FSOC = below this threshold (impaired response)
    control_fsoc_mean <- mean(control_data[[fsoc_var]], na.rm = TRUE)
    control_fsoc_sd <- sd(control_data[[fsoc_var]], na.rm = TRUE)
    
    if (is.null(fsoc_cutoff)) {
      fsoc_cutoff <- control_fsoc_mean - n_sd * control_fsoc_sd
    }
    
    cat(sprintf("Control FSOC: mean = %.3f, SD = %.3f\n", control_fsoc_mean, control_fsoc_sd))
    cat(sprintf("FSOC cutoff (mean - %d SD): %.3f\n", n_sd, fsoc_cutoff))
    
    # Secondary variable cutoff: Control mean + n_sd * SD
    # High K2/R2* = above this threshold (elevated metabolism/hypoxia)
    control_sec_mean <- mean(control_data[[secondary_var]], na.rm = TRUE)
    control_sec_sd <- sd(control_data[[secondary_var]], na.rm = TRUE)
    
    if (is.null(secondary_cutoff)) {
      secondary_cutoff <- control_sec_mean + n_sd * control_sec_sd
    }
    
    cat(sprintf("Control %s: mean = %.3f, SD = %.3f\n", secondary_var, control_sec_mean, control_sec_sd))
    cat(sprintf("Secondary cutoff (mean + %d SD): %.3f\n", n_sd, secondary_cutoff))
    
  } else if (method == "median") {
    if (is.null(fsoc_cutoff)) fsoc_cutoff <- median(dat_pheno[[fsoc_var]], na.rm = TRUE)
    if (is.null(secondary_cutoff)) secondary_cutoff <- median(dat_pheno[[secondary_var]], na.rm = TRUE)
    cat("Using median-based cutoffs\n")
  } else if (method == "tertile") {
    if (is.null(fsoc_cutoff)) fsoc_cutoff <- quantile(dat_pheno[[fsoc_var]], 0.33, na.rm = TRUE)
    if (is.null(secondary_cutoff)) secondary_cutoff <- quantile(dat_pheno[[secondary_var]], 0.67, na.rm = TRUE)
    cat("Using tertile-based cutoffs\n")
  }
  
  cat(sprintf("\n  FINAL FSOC cutoff: %.3f (values BELOW = Low FSOC)\n", fsoc_cutoff))
  cat(sprintf("  FINAL Secondary cutoff: %.3f (values ABOVE = High K2/R2*)\n\n", secondary_cutoff))
  
  # Create phenotype classifications
  dat_pheno <- dat_pheno %>%
    mutate(
      fsoc_level = ifelse(.data[[fsoc_var]] <= fsoc_cutoff, "Low FSOC", "High FSOC"),
      secondary_level = ifelse(.data[[secondary_var]] >= secondary_cutoff, "High", "Low"),
      
      phenotype_4 = case_when(
        fsoc_level == "Low FSOC" & secondary_level == "High" ~ "Low FSOC / High K2-R2*",
        fsoc_level == "Low FSOC" & secondary_level == "Low" ~ "Low FSOC / Low K2-R2*",
        fsoc_level == "High FSOC" & secondary_level == "High" ~ "High FSOC / High K2-R2*",
        fsoc_level == "High FSOC" & secondary_level == "Low" ~ "High FSOC / Low K2-R2*"
      ),
      
      discordant = ifelse(fsoc_level == "Low FSOC" & secondary_level == "High",
                          "Discordant", "Concordant")
    )
  
  # Phenotype distribution
  cat("PHENOTYPE DISTRIBUTION:\n")
  pheno_table <- dat_pheno %>%
    group_by(phenotype_4) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(pct = round(100 * n / sum(n), 1))
  print(as.data.frame(pheno_table))
  
  # By disease group
  cat("\nPHENOTYPE BY DISEASE GROUP:\n")
  pheno_by_group <- dat_pheno %>%
    count(group, phenotype_4) %>%
    group_by(group) %>%
    mutate(pct = round(100 * n / sum(n), 1))
  print(as.data.frame(pheno_by_group))
  
  # Chi-square test
  cat("\nCHI-SQUARE TEST: Phenotype vs Disease Group\n")
  chi_test <- chisq.test(table(dat_pheno$group, dat_pheno$phenotype_4))
  print(chi_test)
  
  # Discordant phenotype by group
  cat("\nDISCORDANT PHENOTYPE BY GROUP:\n")
  discordant_by_group <- dat_pheno %>%
    group_by(group) %>%
    summarise(
      n_total = n(),
      n_discordant = sum(discordant == "Discordant"),
      pct_discordant = round(100 * n_discordant / n_total, 1),
      .groups = "drop"
    )
  print(as.data.frame(discordant_by_group))
  
  return(list(
    data = dat_pheno,
    cutoffs = list(fsoc = fsoc_cutoff, secondary = secondary_cutoff),
    secondary_var = secondary_var,
    phenotype_table = pheno_table,
    phenotype_by_group = pheno_by_group,
    discordant_by_group = discordant_by_group,
    chi_test = chi_test
  ))
}

# ============================================================================
# PHENOTYPE VISUALIZATION
# ============================================================================

plot_phenotype_scatter <- function(pheno_results, fsoc_var = "medullary_fsoc_abs") {
  
  dat_pheno <- pheno_results$data
  fsoc_cut <- pheno_results$cutoffs$fsoc
  sec_cut <- pheno_results$cutoffs$secondary
  sec_var <- pheno_results$secondary_var
  
  # Determine y-axis label
  if (grepl("r2star", sec_var, ignore.case = TRUE)) {
    y_label <- "Baseline R2* (s⁻¹)\n(Higher = more hypoxic)"
  } else {
    y_label <- "K2 - TCA metabolism (s⁻¹)"
  }
  
  # Count subjects in each quadrant
  n_discordant <- sum(dat_pheno$discordant == "Discordant", na.rm = TRUE)
  n_total <- nrow(dat_pheno)
  
  p <- ggplot(dat_pheno, aes(x = .data[[fsoc_var]], y = .data[[sec_var]])) +
    # Quadrant shading - highlight the "at risk" quadrant
    annotate("rect", xmin = -Inf, xmax = fsoc_cut, ymin = sec_cut, ymax = Inf,
             fill = "#FFCCCC", alpha = 0.4) +
    annotate("rect", xmin = fsoc_cut, xmax = Inf, ymin = -Inf, ymax = sec_cut,
             fill = "#CCFFCC", alpha = 0.3) +
    annotate("rect", xmin = -Inf, xmax = fsoc_cut, ymin = -Inf, ymax = sec_cut,
             fill = "#FFFFCC", alpha = 0.2) +
    annotate("rect", xmin = fsoc_cut, xmax = Inf, ymin = sec_cut, ymax = Inf,
             fill = "#CCE5FF", alpha = 0.2) +
    # Points
    geom_point(aes(color = group, shape = discordant), size = 3.5, alpha = 0.8) +
    # Cutoff lines
    geom_vline(xintercept = fsoc_cut, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = sec_cut, linetype = "dashed", color = "red", linewidth = 0.8) +
    # Quadrant labels
    annotate("label", x = min(dat_pheno[[fsoc_var]], na.rm = TRUE), 
             y = max(dat_pheno[[sec_var]], na.rm = TRUE),
             label = "DISCORDANT\n(At Risk)", fontface = "bold", size = 3, 
             color = "#CC0000", fill = "#FFEEEE", hjust = 0, vjust = 1) +
    annotate("label", x = max(dat_pheno[[fsoc_var]], na.rm = TRUE), 
             y = min(dat_pheno[[sec_var]], na.rm = TRUE),
             label = "Normal\nResponse", fontface = "bold", size = 3, 
             color = "#006600", fill = "#EEFFEE", hjust = 1, vjust = 0) +
    # Theme
    theme_bw(base_size = 12) +
    theme(legend.position = "right", panel.grid.minor = element_blank()) +
    labs(
      x = "FSOC - Furosemide O₂ Consumption (s⁻¹)",
      y = y_label,
      title = "FSOC Phenotype Classification (Lean Control Reference)",
      subtitle = paste0(
        "Cutoffs based on Lean Control mean ± 1 SD | ",
        "Discordant: n=", n_discordant, " (", round(100*n_discordant/n_total, 1), "%)\n",
        "FSOC cutoff: ", round(fsoc_cut, 2), " | K2/R2* cutoff: ", round(sec_cut, 2)
      ),
      color = "Disease Group",
      shape = "Phenotype"
    ) +
    scale_color_brewer(palette = "Set2") +
    scale_shape_manual(values = c("Concordant" = 16, "Discordant" = 17),
                       labels = c("Concordant", "Discordant (High K2/R2*, Low FSOC)"))
  
  return(p)
}

plot_discordant_by_group <- function(pheno_results) {
  
  discordant_df <- pheno_results$discordant_by_group
  
  p <- ggplot(discordant_df, aes(x = reorder(group, -pct_discordant), 
                                 y = pct_discordant, fill = group)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = paste0(n_discordant, "/", n_total)), vjust = -0.3, size = 3.5) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Disease Group", y = "% Discordant (High K2/R2*, Low FSOC)",
         title = "Prevalence of Discordant Phenotype by Group") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  return(p)
}

# ============================================================================
# 12B. CLINICAL CHARACTERISTICS BY PHENOTYPE
# ============================================================================

compare_phenotype_characteristics <- function(pheno_results, dat_full = NULL) {
  
  cat("\n========================================\n")
  cat("CLINICAL CHARACTERISTICS BY PHENOTYPE\n")
  cat("========================================\n\n")
  
  dat_pheno <- pheno_results$data
  
  # If full dataset provided, merge additional variables
  if (!is.null(dat_full)) {
    merge_vars <- setdiff(names(dat_full), names(dat_pheno))
    merge_vars <- c("record_id", merge_vars)
    dat_pheno <- left_join(dat_pheno, dat_full[, merge_vars, drop = FALSE], by = "record_id")
  }
  
  available_vars <- names(dat_pheno)
  
  # Define clinical variables to compare
  clinical_var_map <- list(
    "Age" = c("age"),
    "BMI" = c("bmi", "BMI"),
    "Weight" = c("weight"),
    "HbA1c" = c("hba1c", "HbA1c"),
    "Fasting Glucose" = c("fasting_glucose", "glucose_fasting", "fbg"),
    "Diabetes Duration" = c("diabetes_duration", "dm_duration"),
    "eGFR (FAS)" = c("gfr_fas_cr", "eGFR_CKD_epi", "egfr"),
    "mGFR (Iohexol)" = c("gfr_raw_plasma", "gfr_bsa_plasma"),
    "ERPF" = c("erpf_raw_plasma", "erpf_bsa_plasma"),
    "ACR" = c("acr_u", "uacr"),
    "Filtration Fraction" = c("ff", "FF"),
    "SBP" = c("sbp", "SBP", "map_sys"),
    "DBP" = c("dbp", "DBP", "map_dia"),
    "MAP" = c("map", "MAP"),
    "Fat Mass (kg)" = c("dexa_fat_kg", "fat_mass"),
    "Lean Mass (kg)" = c("dexa_lean_kg", "lean_mass"),
    "Trunk Fat (kg)" = c("dexa_trunk_kg", "trunk_fat"),
    "Body Fat %" = c("dexa_body_fat", "body_fat_pct"),
    "Total Cholesterol" = c("cholesterol", "total_cholesterol"),
    "LDL" = c("ldl", "LDL"),
    "HDL" = c("hdl", "HDL"),
    "Triglycerides" = c("triglycerides", "tg", "TG"),
    "Medullary FSOC" = c("medullary_fsoc_abs"),
    "Whole Kidney FSOC" = c("whole_kidney_fsoc_abs"),
    "Cortical FSOC" = c("cortical_fsoc_abs"),
    "Baseline R2* (Medulla)" = c("medullary_r2star_bl", "bold_medulla_bl"),
    "Baseline R2* (Cortex)" = c("cortical_r2star_bl", "bold_cortex_bl")
  )
  
  # Find which variables exist
  vars_to_compare <- list()
  for (var_name in names(clinical_var_map)) {
    possible_names <- clinical_var_map[[var_name]]
    found_var <- possible_names[possible_names %in% available_vars][1]
    if (!is.na(found_var)) {
      vars_to_compare[[var_name]] <- found_var
    }
  }
  
  cat("Variables found for comparison:", length(vars_to_compare), "\n")
  cat(paste(names(vars_to_compare), collapse = ", "), "\n\n")
  
  # Check if we have discordant subjects
  n_discord <- sum(dat_pheno$discordant == "Discordant", na.rm = TRUE)
  n_concord <- sum(dat_pheno$discordant == "Concordant", na.rm = TRUE)
  
  # =========================================
  # COMPARISON 1: Discordant vs Concordant
  # =========================================
  cat("-------------------------------------------\n")
  cat("COMPARISON 1: DISCORDANT vs CONCORDANT\n")
  cat("-------------------------------------------\n\n")
  
  cat(sprintf("Discordant: n = %d\n", n_discord))
  cat(sprintf("Concordant: n = %d\n\n", n_concord))
  
  results_binary <- data.frame()
  
  if (n_discord >= 2 & n_concord >= 2) {
    # Proceed with comparison
    for (var_name in names(vars_to_compare)) {
      var_col <- vars_to_compare[[var_name]]
      
      discord_vals <- dat_pheno[[var_col]][dat_pheno$discordant == "Discordant"]
      concord_vals <- dat_pheno[[var_col]][dat_pheno$discordant == "Concordant"]
      
      discord_vals <- discord_vals[!is.na(discord_vals)]
      concord_vals <- concord_vals[!is.na(concord_vals)]
      
      if (length(discord_vals) >= 2 & length(concord_vals) >= 2) {
        wilcox_p <- tryCatch({
          wilcox.test(discord_vals, concord_vals)$p.value
        }, error = function(e) NA)
        
        results_binary <- rbind(results_binary, data.frame(
          Variable = var_name,
          Discordant_mean = round(mean(discord_vals), 2),
          Discordant_sd = round(sd(discord_vals), 2),
          Concordant_mean = round(mean(concord_vals), 2),
          Concordant_sd = round(sd(concord_vals), 2),
          Difference = round(mean(discord_vals) - mean(concord_vals), 2),
          P_value = round(wilcox_p, 4)
        ))
      }
    }
    
    if (nrow(results_binary) > 0) {
      results_binary <- results_binary %>% 
        arrange(P_value) %>%
        mutate(Sig = case_when(
          P_value < 0.001 ~ "***",
          P_value < 0.01 ~ "**",
          P_value < 0.05 ~ "*",
          TRUE ~ ""
        ))
      print(as.data.frame(results_binary))
    }
  } else {
    cat("** NOT ENOUGH SUBJECTS IN ONE OR BOTH GROUPS FOR COMPARISON **\n")
    cat("Consider using a less strict cutoff (e.g., n_sd = 0.5) or median-based method.\n\n")
  }
  
  # Sex distribution (with error handling)
  cat("\n\nSEX DISTRIBUTION:\n")
  sex_table <- dat_pheno %>% count(discordant, sex)
  print(as.data.frame(sex_table))
  
  if (n_discord >= 1 & n_concord >= 1 & length(unique(dat_pheno$sex)) >= 2) {
    sex_fisher <- tryCatch({
      fisher.test(table(dat_pheno$discordant, dat_pheno$sex))
    }, error = function(e) NULL)
    
    if (!is.null(sex_fisher)) {
      cat(sprintf("\nFisher's exact test p-value: %.4f\n", sex_fisher$p.value))
    }
  }
  
  # Disease group distribution
  cat("\n\nDISEASE GROUP DISTRIBUTION:\n")
  group_table <- dat_pheno %>% count(discordant, group)
  print(as.data.frame(group_table))
  
  if (n_discord >= 1 & n_concord >= 1 & length(unique(dat_pheno$group)) >= 2) {
    group_fisher <- tryCatch({
      fisher.test(table(dat_pheno$discordant, dat_pheno$group))
    }, error = function(e) NULL)
    
    if (!is.null(group_fisher)) {
      cat(sprintf("\nFisher's exact test p-value: %.4f\n", group_fisher$p.value))
    }
  }
  
  # =========================================
  # COMPARISON 2: All 4 Phenotype Groups
  # =========================================
  cat("\n\n-------------------------------------------\n")
  cat("COMPARISON 2: ALL PHENOTYPE GROUPS (Kruskal-Wallis)\n")
  cat("-------------------------------------------\n\n")
  
  pheno_n <- dat_pheno %>% count(phenotype_4)
  print(as.data.frame(pheno_n))
  
  results_4group <- data.frame(Variable = character(), KW_pvalue = numeric())
  
  # Only run if we have at least 2 groups with data
  n_groups_with_data <- sum(pheno_n$n >= 2)
  
  if (n_groups_with_data >= 2) {
    for (var_name in names(vars_to_compare)) {
      var_col <- vars_to_compare[[var_name]]
      
      kw_p <- tryCatch({
        kruskal.test(dat_pheno[[var_col]] ~ dat_pheno$phenotype_4)$p.value
      }, error = function(e) NA)
      
      results_4group <- rbind(results_4group, data.frame(
        Variable = var_name,
        KW_pvalue = round(kw_p, 4)
      ))
    }
    
    results_4group <- results_4group %>% arrange(KW_pvalue)
    
    cat("\nKruskal-Wallis test results (sorted by p-value):\n")
    print(as.data.frame(results_4group %>% filter(!is.na(KW_pvalue)) %>% head(15)))
  } else {
    cat("Not enough groups with sufficient sample size for comparison.\n")
  }
  
  return(list(
    binary_comparison = results_binary,
    fourgroup_comparison = results_4group,
    sex_distribution = sex_table,
    group_distribution = group_table
  ))
}

# ============================================================================
# VISUALIZATION: Clinical Characteristics Forest Plot
# ============================================================================

plot_phenotype_characteristics <- function(comparison_results) {
  
  binary_df <- comparison_results$binary_comparison
  
  # Filter to variables with sufficient data and reasonable effect sizes
  plot_df <- binary_df %>%
    filter(!is.na(P_value)) %>%
    mutate(
      # Standardized difference (Cohen's d approximation)
      pooled_sd = sqrt((Discordant_sd^2 + Concordant_sd^2) / 2),
      effect_size = Difference / pooled_sd,
      Variable = factor(Variable, levels = rev(Variable))  # Reverse for plotting
    ) %>%
    filter(abs(effect_size) < 5)  # Remove extreme outliers
  
  # Create forest plot
  p <- ggplot(plot_df, aes(x = effect_size, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = -log10(P_value), color = P_value < 0.05), alpha = 0.7) +
    geom_segment(aes(x = 0, xend = effect_size, y = Variable, yend = Variable,
                     color = P_value < 0.05), alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#E41A1C"),
                       labels = c("p ≥ 0.05", "p < 0.05"),
                       name = "Significance") +
    scale_size_continuous(range = c(2, 6), name = "-log10(p)") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 9)
    ) +
    labs(
      x = "Standardized Difference (Discordant - Concordant)",
      y = "",
      title = "Clinical Characteristics: Discordant vs Concordant Phenotype",
      subtitle = "Positive values = higher in discordant group"
    )
  
  return(p)
}

# ============================================================================
# VISUALIZATION: Boxplots of Key Variables by Phenotype
# ============================================================================

plot_key_variables_by_phenotype <- function(pheno_results, vars_to_plot = NULL) {
  
  dat_pheno <- pheno_results$data
  available_vars <- names(dat_pheno)
  
  # Default variables to plot
  if (is.null(vars_to_plot)) {
    default_vars <- c("gfr_fas_cr", "hba1c", "acr_u", "bmi", "sbp", 
                      "dexa_fat_kg", "medullary_fsoc_abs")
    vars_to_plot <- default_vars[default_vars %in% available_vars]
  }
  
  plots <- list()
  
  for (var in vars_to_plot) {
    # Skip if too many NAs
    if (sum(!is.na(dat_pheno[[var]])) < 10) next
    
    # Clean variable name for title
    var_label <- gsub("_", " ", var)
    var_label <- tools::toTitleCase(var_label)
    
    # Kruskal-Wallis p-value
    kw_p <- tryCatch({
      kruskal.test(dat_pheno[[var]] ~ dat_pheno$phenotype_4)$p.value
    }, error = function(e) NA)
    
    p <- ggplot(dat_pheno, aes(x = phenotype_4, y = .data[[var]], fill = phenotype_4)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
      theme_bw(base_size = 10) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
      ) +
      labs(
        x = "",
        y = var_label,
        title = var_label,
        subtitle = ifelse(!is.na(kw_p), sprintf("KW p = %.3f", kw_p), "")
      ) +
      scale_fill_manual(values = c(
        "Low FSOC / High K2-R2*" = "#E41A1C",
        "Low FSOC / Low K2-R2*" = "#377EB8",
        "High FSOC / High K2-R2*" = "#FF7F00",
        "High FSOC / Low K2-R2*" = "#4DAF4A"
      ))
    
    plots[[var]] <- p
  }
  
  if (length(plots) >= 4) {
    combined <- wrap_plots(plots, ncol = 2) +
      plot_annotation(title = "Key Variables by FSOC-K2/R2* Phenotype")
  } else if (length(plots) > 0) {
    combined <- wrap_plots(plots, ncol = length(plots))
  } else {
    return(NULL)
  }
  
  return(combined)
}

# ============================================================================
# 12C. SENSITIVITY ANALYSIS: DIFFERENT CUTOFF THRESHOLDS
# ============================================================================

run_cutoff_sensitivity <- function(dat, 
                                   secondary_var,
                                   fsoc_var = "medullary_fsoc_abs",
                                   control_group = "Lean Control",
                                   sd_values = c(0.5, 1, 1.5, 2)) {
  
  cat("\n========================================\n")
  cat("SENSITIVITY ANALYSIS: CUTOFF THRESHOLDS\n")
  cat("========================================\n\n")
  
  results <- data.frame()
  
  for (n_sd in sd_values) {
    
    pheno_temp <- create_fsoc_phenotypes(
      dat, 
      secondary_var = secondary_var,
      fsoc_var = fsoc_var,
      method = "control_reference",
      control_group = control_group,
      n_sd = n_sd
    )
    
    if (!is.null(pheno_temp)) {
      n_total <- nrow(pheno_temp$data)
      n_discordant <- sum(pheno_temp$data$discordant == "Discordant")
      
      # Get counts by phenotype
      pheno_counts <- pheno_temp$data %>%
        count(phenotype_4) %>%
        pivot_wider(names_from = phenotype_4, values_from = n, values_fill = 0)
      
      results <- rbind(results, data.frame(
        SD_threshold = n_sd,
        FSOC_cutoff = round(pheno_temp$cutoffs$fsoc, 3),
        K2_cutoff = round(pheno_temp$cutoffs$secondary, 3),
        N_total = n_total,
        N_discordant = n_discordant,
        Pct_discordant = round(100 * n_discordant / n_total, 1)
      ))
    }
  }
  
  cat("\n\n========================================\n")
  cat("SENSITIVITY ANALYSIS SUMMARY\n")
  cat("========================================\n\n")
  print(as.data.frame(results))
  
  # Recommendation
  cat("\n\nRECOMMENDATION:\n")
  if (any(results$N_discordant >= 5)) {
    best_sd <- results$SD_threshold[which(results$N_discordant >= 5)[1]]
    cat(sprintf("Use n_sd = %.1f to get at least 5 discordant subjects for meaningful comparison.\n", best_sd))
  } else if (any(results$N_discordant >= 3)) {
    best_sd <- results$SD_threshold[which(results$N_discordant >= 3)[1]]
    cat(sprintf("Use n_sd = %.1f to get at least 3 discordant subjects (limited power).\n", best_sd))
  } else {
    cat("No threshold yields sufficient discordant subjects.\n")
    cat("Consider: (1) using median-based cutoffs, or (2) analyzing K2 and FSOC as continuous variables.\n")
  }
  
  return(results)
}

# ============================================================================
# 12D. CONTINUOUS ANALYSIS: K2 vs FSOC CORRELATION
# ============================================================================

analyze_k2_fsoc_continuous <- function(dat, 
                                       k2_var,
                                       fsoc_var = "medullary_fsoc_abs") {
  
  cat("\n========================================\n")
  cat("K2 vs FSOC: CONTINUOUS ANALYSIS\n")
  cat("========================================\n\n")
  
  dat_clean <- dat %>%
    filter(!is.na(group)) %>%
    filter(!is.na(.data[[fsoc_var]])) %>%
    filter(!is.na(.data[[k2_var]])) %>%
    filter(.data[[fsoc_var]] >= 0)
  
  cat("N subjects:", nrow(dat_clean), "\n\n")
  
  # Overall correlation
  cor_overall <- cor.test(dat_clean[[fsoc_var]], dat_clean[[k2_var]], 
                          method = "spearman", exact = FALSE)
  
  cat("OVERALL CORRELATION (Spearman):\n")
  cat(sprintf("  rho = %.3f, p = %.4f\n\n", cor_overall$estimate, cor_overall$p.value))
  
  # By group
  cat("CORRELATION BY GROUP:\n")
  cor_by_group <- dat_clean %>%
    group_by(group) %>%
    summarise(
      n = n(),
      rho = cor(.data[[fsoc_var]], .data[[k2_var]], method = "spearman", use = "complete.obs"),
      .groups = "drop"
    )
  print(as.data.frame(cor_by_group))
  
  # Linear regression: FSOC ~ K2 + covariates
  cat("\n\nLINEAR REGRESSION: FSOC ~ K2 + age + sex + group\n")
  
  formula_str <- paste(fsoc_var, "~", k2_var, "+ age + sex + group")
  model <- lm(as.formula(formula_str), data = dat_clean)
  
  cat("\nModel summary:\n")
  print(summary(model)$coefficients)
  cat(sprintf("\nR-squared: %.3f\n", summary(model)$r.squared))
  
  # Scatter plot
  p <- ggplot(dat_clean, aes(x = .data[[k2_var]], y = .data[[fsoc_var]])) +
    geom_point(aes(color = group), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linetype = "dashed") +
    geom_smooth(aes(color = group), method = "lm", se = FALSE, linewidth = 0.8) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    labs(
      x = paste0(k2_var, " (K2 - TCA metabolism)"),
      y = paste0(fsoc_var, " (Furosemide O2 Consumption)"),
      title = "K2 vs FSOC Relationship",
      subtitle = sprintf("Overall: rho = %.2f, p = %.3f", cor_overall$estimate, cor_overall$p.value),
      color = "Disease Group"
    ) +
    scale_color_brewer(palette = "Set2")
  
  return(list(
    correlation_overall = cor_overall,
    correlation_by_group = cor_by_group,
    model = model,
    plot = p
  ))
}

# ============================================================================
# 13. FSOC CUTOFF ANALYSIS
# ============================================================================

analyze_fsoc_cutoffs <- function(dat, fsoc_var = "medullary_fsoc_abs") {
  
  cat("\n========================================\n")
  cat("FSOC CUTOFF ANALYSIS:", fsoc_var, "\n")
  cat("========================================\n\n")
  
  dat_clean <- dat %>%
    filter(!is.na(group), !is.na(.data[[fsoc_var]]), .data[[fsoc_var]] >= 0)
  
  cat("N subjects:", nrow(dat_clean), "\n\n")
  
  # Distribution statistics
  fsoc_values <- dat_clean[[fsoc_var]]
  
  cat("DISTRIBUTION STATISTICS:\n")
  stats_df <- data.frame(
    Stat = c("Mean", "SD", "Median", "IQR", "Min", "Max", "Q1", "Q3"),
    Value = c(
      round(mean(fsoc_values), 3), round(sd(fsoc_values), 3),
      round(median(fsoc_values), 3), round(IQR(fsoc_values), 3),
      round(min(fsoc_values), 3), round(max(fsoc_values), 3),
      round(quantile(fsoc_values, 0.25), 3), round(quantile(fsoc_values, 0.75), 3)
    )
  )
  print(stats_df)
  
  # By disease group
  cat("\nBY DISEASE GROUP:\n")
  stats_by_group <- dat_clean %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = round(mean(.data[[fsoc_var]]), 3),
      sd = round(sd(.data[[fsoc_var]]), 3),
      median = round(median(.data[[fsoc_var]]), 3),
      Q1 = round(quantile(.data[[fsoc_var]], 0.25), 3),
      Q3 = round(quantile(.data[[fsoc_var]], 0.75), 3),
      .groups = "drop"
    )
  print(as.data.frame(stats_by_group))
  
  # Tertile analysis
  cat("\nTERTILE ANALYSIS:\n")
  t1 <- quantile(fsoc_values, 1/3)
  t2 <- quantile(fsoc_values, 2/3)
  cat(sprintf("Tertile cutoffs: T1 < %.3f, T2 < %.3f, T3 >= %.3f\n", t1, t2, t2))
  
  dat_clean <- dat_clean %>%
    mutate(fsoc_tertile = cut(.data[[fsoc_var]],
                              breaks = c(-Inf, t1, t2, Inf),
                              labels = c("Low (T1)", "Medium (T2)", "High (T3)")))
  
  tertile_table <- dat_clean %>%
    count(group, fsoc_tertile) %>%
    group_by(group) %>%
    mutate(pct = round(100 * n / sum(n), 1))
  print(as.data.frame(tertile_table))
  
  # Chi-square
  cat("\nChi-square: Tertile vs Group\n")
  chi_test <- chisq.test(table(dat_clean$group, dat_clean$fsoc_tertile))
  print(chi_test)
  
  return(list(
    stats = stats_df,
    stats_by_group = stats_by_group,
    tertile_cutoffs = c(t1, t2),
    data_with_tertiles = dat_clean
  ))
}

# ============================================================================
# 14. CONTINUOUS ASSOCIATION ANALYSIS
# ============================================================================

analyze_fsoc_continuous <- function(dat, fsoc_var = "medullary_fsoc_abs") {
  
  cat("\n========================================\n")
  cat("FSOC CONTINUOUS ASSOCIATIONS\n")
  cat("========================================\n\n")
  
  dat_clean <- dat %>%
    filter(!is.na(group), !is.na(.data[[fsoc_var]]), .data[[fsoc_var]] >= 0)
  
  available_vars <- names(dat_clean)
  
  # Define outcomes to test
  outcomes <- list(
    "eGFR" = c("eGFR_CKD_epi", "gfr_fas_cr"),
    "mGFR" = c("gfr_raw_plasma", "gfr_bsa_plasma"),
    "HbA1c" = c("hba1c"),
    "ACR" = c("acr_u")
  )
  
  results <- list()
  
  for (outcome_name in names(outcomes)) {
    var_found <- outcomes[[outcome_name]][outcomes[[outcome_name]] %in% available_vars][1]
    
    if (!is.na(var_found)) {
      cat("\n---", outcome_name, "---\n")
      
      # Prepare outcome
      if (outcome_name == "ACR") {
        dat_model <- dat_clean %>% mutate(outcome = log(.data[[var_found]] + 1))
      } else {
        dat_model <- dat_clean %>% mutate(outcome = .data[[var_found]])
      }
      
      # Spearman correlation
      cor_test <- cor.test(dat_model[[fsoc_var]], dat_model$outcome, 
                           method = "spearman", exact = FALSE)
      cat(sprintf("Spearman rho = %.3f, p = %.4f\n", cor_test$estimate, cor_test$p.value))
      
      # Linear model (adjusted)
      formula_adj <- as.formula(paste("outcome ~", fsoc_var, "+ age + sex + group"))
      model_adj <- lm(formula_adj, data = dat_model)
      coef_adj <- summary(model_adj)$coefficients
      
      if (fsoc_var %in% rownames(coef_adj)) {
        cat(sprintf("Adjusted beta = %.4f, p = %.4f\n",
                    coef_adj[fsoc_var, "Estimate"], coef_adj[fsoc_var, "Pr(>|t|)"]))
      }
      
      results[[outcome_name]] <- list(correlation = cor_test, model = model_adj)
    }
  }
  
  return(results)
}

# ============================================================================
# 15. VISUALIZATION FUNCTIONS
# ============================================================================

plot_fsoc_distribution <- function(dat, fsoc_var = "medullary_fsoc_abs") {
  
  dat_clean <- dat %>%
    filter(!is.na(group), !is.na(.data[[fsoc_var]]), .data[[fsoc_var]] >= 0)
  
  median_val <- median(dat_clean[[fsoc_var]])
  
  p1 <- ggplot(dat_clean, aes(x = .data[[fsoc_var]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, 
                   fill = "#66C2A5", color = "white", alpha = 0.7) +
    geom_density(color = "#1B7837", linewidth = 1) +
    geom_vline(xintercept = median_val, linetype = "dashed", color = "#E7298A", linewidth = 1) +
    theme_bw(base_size = 12) +
    labs(x = "FSOC (s⁻¹)", y = "Density", title = "Overall Distribution",
         subtitle = paste0("Median = ", round(median_val, 2)))
  
  p2 <- ggplot(dat_clean, aes(x = .data[[fsoc_var]], fill = group)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = median_val, linetype = "dashed", color = "grey40") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    labs(x = "FSOC (s⁻¹)", y = "Density", title = "By Disease Group") +
    scale_fill_brewer(palette = "Set2")
  
  p3 <- ggplot(dat_clean, aes(x = group, y = .data[[fsoc_var]], fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    geom_hline(yintercept = median_val, linetype = "dashed", color = "#E7298A") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Disease Group", y = "FSOC (s⁻¹)", title = "Boxplot by Group") +
    scale_fill_brewer(palette = "Set2")
  
  combined <- (p1 | p2) / p3 +
    plot_annotation(title = paste0("FSOC Distribution: ", fsoc_var))
  
  return(combined)
}

# ============================================================================
# 16. RUN ALL ANALYSES
# ============================================================================

cat("\n########################################################\n")
cat("# STARTING PHENOTYPE AND CUTOFF ANALYSES               #\n")
cat("########################################################\n\n")

# Step 1: Identify variables
vars_found <- identify_all_variables(dat)

# Step 2: Calculate BOLD averages if needed
if (any(grepl("bold_.*_bl_", names(dat), ignore.case = TRUE))) {
  dat <- calculate_bold_averages(dat)
  cat("\nBOLD baseline averages calculated.\n")
}

# Step 3: Run phenotype analysis
# Try K2 first, then BOLD baseline R2*
k2_vars <- vars_found$k2
if (length(k2_vars) > 0) {
  secondary_var <- k2_vars[1]
  cat("\nUsing K2 variable for phenotyping:", secondary_var, "\n")
} else if ("medullary_r2star_bl" %in% names(dat)) {
  secondary_var <- "medullary_r2star_bl"
  cat("\nUsing baseline R2* for phenotyping:", secondary_var, "\n")
} else {
  secondary_var <- NULL
  cat("\nNo K2 or R2* baseline variable found for phenotyping.\n")
}

if (!is.null(secondary_var)) {
  # =====================================================
  # PHENOTYPE ANALYSIS WITH LEAN CONTROL REFERENCE
  # =====================================================
  
  pheno_results <- create_fsoc_phenotypes(
    dat, 
    secondary_var = secondary_var,
    method = "control_reference",
    control_group = "Lean Control",
    n_sd = 1
  )
  
  if (!is.null(pheno_results)) {
    # Check if we have any discordant subjects
    n_discordant <- sum(pheno_results$data$discordant == "Discordant", na.rm = TRUE)
    
    if (n_discordant == 0) {
      cat("\n** NO DISCORDANT SUBJECTS WITH 1 SD CUTOFF **\n")
      cat("Running sensitivity analysis with different thresholds...\n\n")
      
      # Run sensitivity analysis
      sensitivity_results <- run_cutoff_sensitivity(
        dat, 
        secondary_var = secondary_var,
        control_group = "Lean Control",
        sd_values = c(0.5, 0.75, 1, 1.5)
      )
      
      # Save sensitivity results
      write.csv(sensitivity_results, 
                file.path(OUTPUT_DIR, "phenotype_sensitivity_analysis.csv"), row.names = FALSE)
      
      # If any threshold gives discordant subjects, re-run with that threshold
      if (any(sensitivity_results$N_discordant >= 3)) {
        best_sd <- sensitivity_results$SD_threshold[which(sensitivity_results$N_discordant >= 3)[1]]
        cat(sprintf("\n** Re-running phenotype analysis with n_sd = %.2f **\n\n", best_sd))
        
        pheno_results <- create_fsoc_phenotypes(
          dat, 
          secondary_var = secondary_var,
          method = "control_reference",
          control_group = "Lean Control",
          n_sd = best_sd
        )
      }
    }
    
    # Scatter plot
    p_scatter <- plot_phenotype_scatter(pheno_results)
    print(p_scatter)
    ggsave(file.path(OUTPUT_DIR, "fsoc_phenotype_scatter.pdf"), p_scatter, width = 10, height = 8)
    
    # Discordant by group
    p_discordant <- plot_discordant_by_group(pheno_results)
    print(p_discordant)
    ggsave(file.path(OUTPUT_DIR, "discordant_phenotype_by_group.pdf"), p_discordant, width = 8, height = 6)
    
    # Clinical characteristics comparison
    cat("\n\n--- Comparing Clinical Characteristics by Phenotype ---\n")
    char_comparison <- compare_phenotype_characteristics(pheno_results, dat_full = dat)
    
    # Export comparison results if we have data
    if (nrow(char_comparison$fourgroup_comparison) > 0) {
      write.csv(char_comparison$fourgroup_comparison,
                file.path(OUTPUT_DIR, "phenotype_group_comparison.csv"), row.names = FALSE)
    }
    
    # Export phenotype data
    write.csv(pheno_results$data %>% 
                dplyr::select(record_id, group, sex, age, medullary_fsoc_abs, 
                              all_of(secondary_var), phenotype_4, discordant),
              file.path(OUTPUT_DIR, "fsoc_phenotypes.csv"), row.names = FALSE)
    
    # =====================================================
    # CONTINUOUS ANALYSIS: K2 vs FSOC
    # =====================================================
    cat("\n\n--- Running K2 vs FSOC Continuous Analysis ---\n")
    k2_fsoc_continuous <- analyze_k2_fsoc_continuous(dat, k2_var = secondary_var)
    
    print(k2_fsoc_continuous$plot)
    ggsave(file.path(OUTPUT_DIR, "k2_vs_fsoc_scatter.pdf"), k2_fsoc_continuous$plot, 
           width = 9, height = 7)
  }
}

if (!is.null(pheno_results)) {
  # Visualizations
  p_scatter <- plot_phenotype_scatter(pheno_results)
  print(p_scatter)
  ggsave(file.path(OUTPUT_DIR, "fsoc_phenotype_scatter.pdf"), p_scatter, width = 10, height = 8)
  
  p_discordant <- plot_discordant_by_group(pheno_results)
  print(p_discordant)
  ggsave(file.path(OUTPUT_DIR, "discordant_phenotype_by_group.pdf"), p_discordant, width = 8, height = 6)
  
  # =====================================================
  # NEW: Clinical characteristics comparison
  # =====================================================
  cat("\n\n--- Comparing Clinical Characteristics by Phenotype ---\n")
  char_comparison <- compare_phenotype_characteristics(pheno_results, dat_full = dat)
  
  # Forest plot of differences
  p_forest <- plot_phenotype_characteristics(char_comparison)
  print(p_forest)
  ggsave(file.path(OUTPUT_DIR, "phenotype_characteristics_forest.pdf"), p_forest, 
         width = 10, height = 8)
  
  # Boxplots of key variables
  p_boxplots <- plot_key_variables_by_phenotype(pheno_results)
  if (!is.null(p_boxplots)) {
    print(p_boxplots)
    ggsave(file.path(OUTPUT_DIR, "phenotype_key_variables_boxplots.pdf"), p_boxplots, 
           width = 12, height = 10)
  }
  
  # Export comparison results
  write.csv(char_comparison$binary_comparison, 
            file.path(OUTPUT_DIR, "phenotype_discordant_vs_concordant.csv"), row.names = FALSE)
  write.csv(char_comparison$fourgroup_comparison,
            file.path(OUTPUT_DIR, "phenotype_4group_comparison.csv"), row.names = FALSE)
  
  # Export phenotype data
  write.csv(pheno_results$data %>% 
              select(record_id, group, sex, age, medullary_fsoc_abs, 
                     all_of(secondary_var), phenotype_4, discordant),
            file.path(OUTPUT_DIR, "fsoc_phenotypes.csv"), row.names = FALSE)
}
}

# Step 4: Cutoff analysis
cat("\n\n--- Running FSOC Cutoff Analysis ---\n")
cutoff_med <- analyze_fsoc_cutoffs(dat, "medullary_fsoc_abs")
cutoff_wk <- analyze_fsoc_cutoffs(dat, "whole_kidney_fsoc_abs")

# Step 5: Continuous associations
cat("\n\n--- Running Continuous Association Analysis ---\n")
continuous_results <- analyze_fsoc_continuous(dat, "medullary_fsoc_abs")

# Step 6: Distribution plots
p_dist <- plot_fsoc_distribution(dat, "medullary_fsoc_abs")
print(p_dist)
ggsave(file.path(OUTPUT_DIR, "fsoc_distribution.pdf"), p_dist, width = 12, height = 10)

# Export results
write.csv(cutoff_med$stats_by_group, file.path(OUTPUT_DIR, "fsoc_stats_by_group.csv"), row.names = FALSE)

cat("\n########################################################\n")
cat("# ANALYSES COMPLETE                                    #\n")
cat("########################################################\n")
cat("\nFiles saved to:", OUTPUT_DIR, "\n")










