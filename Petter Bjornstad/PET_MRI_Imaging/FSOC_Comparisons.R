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
    # Filter to subjects with at least one FSOC measurement
    filter(!is.na(whole_kidney_fsoc_abs) | !is.na(medullary_fsoc_abs)) %>%
    # Also filter out invalid FSOC values (matching other analyses in script)
    filter(is.na(whole_kidney_fsoc_abs) | whole_kidney_fsoc_abs < 15) %>%
    filter(is.na(medullary_fsoc_abs) | medullary_fsoc_abs >= 0) %>%
    select(group, age, sex, bmi, hba1c, diabetes_duration,
           eGFR_CKD_epi, acr_u,
           whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    tbl_summary(
      by = group,
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)",
        # Override UACR to show median (IQR)
        acr_u ~ "{median} ({p25}, {p75})"
      ),
      digits = list(
        all_continuous() ~ 1,
        acr_u ~ 1
      ),
      label = list(
        age ~ "Age, years",
        sex ~ "Sex",
        bmi ~ "BMI, kg/m²",
        hba1c ~ "HbA1c, %",
        diabetes_duration ~ "Diabetes Duration, years",
        eGFR_CKD_epi ~ "eGFR, mL/min/1.73m²",
        acr_u ~ "UACR, mg/g",
        whole_kidney_fsoc_abs ~ "Whole Kidney FSOC, s⁻¹",
        medullary_fsoc_abs ~ "Medullary FSOC, s⁻¹"
      ),
      missing_text = "Missing"
    ) %>%
    add_p(
      test = list(
        acr_u ~ "kruskal.test"
      )
    ) %>%
    add_overall() %>%
    modify_header(label ~ "**Characteristic**") %>%
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
# DISCORDANT ANALYSIS WITH LEAN CONTROL REFERENCE CUTOFFS
# UPDATED: Removed lipids, Added DEXA parameters and mGFR
# UPDATED: Filter out variables with insufficient sample sizes
# ============================================================================

# Load ggrepel if available, otherwise we'll use geom_text
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  cat("Note: ggrepel package not installed. Using geom_text instead.\n")
  cat("Install with: install.packages('ggrepel')\n")
  use_repel <- FALSE
} else {
  library(ggrepel)
  use_repel <- TRUE
}

# Define custom color palette for groups
group_colors <- c(
  "Type 2 Diabetes" = "#E41A1C",
  "Type 1 Diabetes" = "#984EA3", 
  "Lean Control" = "#4DAF4A",
  "Obese Control" = "#377EB8"
)

# Phenotype colors
phenotype_colors <- c(
  "High K2 / Low FSOC" = "#E41A1C",
  "High K2 / Normal FSOC" = "#FF7F00",
  "Normal K2 / Low FSOC" = "#FFFF33",
  "Normal K2 / Normal FSOC" = "#4DAF4A"
)

# ============================================================================
# STEP 1: CREATE LINKED DATASET
# UPDATED: Added DEXA and mGFR variables to possible_demo_vars
# ============================================================================

create_linked_k2_fsoc_dataset <- function(dat) {
  
  cat("\n========================================\n")
  cat("CREATING LINKED K2-FSOC DATASET\n")
  cat("========================================\n\n")
  
  # UPDATED: Added DEXA and mGFR, removed lipids
  possible_demo_vars <- c("sex", "age", "bmi", "hba1c", "eGFR_CKD_epi", "gfr_fas_cr",
                          "gfr_raw_plasma", "gfr_bsa_plasma",  # mGFR variables
                          "acr_u", "uacr", "diabetes_duration", "dm_duration",
                          "sbp", "dbp", "weight", "height", "waist", "hip",
                          "dexa_fat_kg", "dexa_trunk_kg", "dexa_lean_kg", "dexa_body_fat",  # DEXA
                          "creatinine", "cystatin_c")
  available_demo_vars <- possible_demo_vars[possible_demo_vars %in% names(dat)]
  
  # Extract K2 data
  k2_select_vars <- c("record_id", "mrn", "group", available_demo_vars, "avg_c_k2", "avg_m_k2")
  k2_select_vars <- unique(k2_select_vars[k2_select_vars %in% names(dat)])
  
  k2_data <- dat %>%
    filter(!is.na(avg_c_k2) | !is.na(avg_m_k2)) %>%
    filter(!is.na(mrn) & mrn != "") %>%
    select(all_of(k2_select_vars)) %>%
    rename(k2_record_id = record_id, k2_group = group) %>%
    distinct(mrn, .keep_all = TRUE)
  
  cat("Subjects with K2 data and MRN:", nrow(k2_data), "\n")
  
  # Extract FSOC data
  fsoc_select_vars <- c("record_id", "mrn", "group", 
                        "fsoc_l_medulla", "fsoc_r_medulla", 
                        "fsoc_l_kidney", "fsoc_r_kidney")
  fsoc_select_vars <- fsoc_select_vars[fsoc_select_vars %in% names(dat)]
  
  fsoc_data <- dat %>%
    filter(!is.na(fsoc_l_medulla) | !is.na(fsoc_r_medulla)) %>%
    filter(!is.na(mrn) & mrn != "") %>%
    select(all_of(fsoc_select_vars)) %>%
    rename(fsoc_record_id = record_id, fsoc_group = group) %>%
    distinct(mrn, .keep_all = TRUE)
  
  cat("Subjects with FSOC data and MRN:", nrow(fsoc_data), "\n")
  
  # Join by MRN
  linked_data <- inner_join(k2_data, fsoc_data, by = "mrn")
  
  cat("\n*** LINKED SUBJECTS:", nrow(linked_data), "***\n")
  
  if (nrow(linked_data) == 0) return(NULL)
  
  # Calculate averaged FSOC and log-transform ACR
  linked_data <- linked_data %>%
    mutate(
      fsoc_medulla = case_when(
        !is.na(fsoc_l_medulla) & !is.na(fsoc_r_medulla) ~ (fsoc_l_medulla + fsoc_r_medulla) / 2,
        !is.na(fsoc_l_medulla) ~ fsoc_l_medulla,
        !is.na(fsoc_r_medulla) ~ fsoc_r_medulla,
        TRUE ~ NA_real_
      ),
      group = k2_group,
      # Log10 transform ACR
      log_acr = ifelse(!is.na(acr_u) & acr_u > 0, log10(acr_u), NA_real_)
    )
  
  return(linked_data)
}

# ============================================================================
# STEP 2: PHENOTYPE WITH LEAN CONTROL CUTOFFS
# ============================================================================

run_lean_control_phenotype <- function(linked_data, 
                                       k2_var = "avg_c_k2",
                                       fsoc_var = "fsoc_medulla",
                                       n_sd = 1) {
  
  cat("\n========================================\n")
  cat("PHENOTYPE ANALYSIS: LEAN CONTROL REFERENCE\n")
  cat("========================================\n\n")
  
  dat_valid <- linked_data %>%
    filter(!is.na(.data[[k2_var]])) %>%
    filter(!is.na(.data[[fsoc_var]]))
  
  cat("Total subjects with both K2 and FSOC:", nrow(dat_valid), "\n\n")
  
  lean_controls <- dat_valid %>%
    filter(grepl("Lean", group, ignore.case = TRUE))
  
  cat("Lean Control subjects:", nrow(lean_controls), "\n")
  
  if (nrow(lean_controls) < 2) {
    cat("WARNING: Not enough Lean Controls. Using all non-diabetic controls.\n")
    lean_controls <- dat_valid %>%
      filter(grepl("Control", group, ignore.case = TRUE))
    cat("Control subjects (all):", nrow(lean_controls), "\n")
  }
  
  lc_k2_mean <- mean(lean_controls[[k2_var]], na.rm = TRUE)
  lc_k2_sd <- sd(lean_controls[[k2_var]], na.rm = TRUE)
  lc_fsoc_mean <- mean(lean_controls[[fsoc_var]], na.rm = TRUE)
  lc_fsoc_sd <- sd(lean_controls[[fsoc_var]], na.rm = TRUE)
  
  cat(sprintf("\nLean Control K2: mean = %.4f, SD = %.4f\n", lc_k2_mean, lc_k2_sd))
  cat(sprintf("Lean Control FSOC: mean = %.2f, SD = %.2f\n", lc_fsoc_mean, lc_fsoc_sd))
  
  k2_cutoff <- lc_k2_mean + n_sd * lc_k2_sd
  fsoc_cutoff <- lc_fsoc_mean - n_sd * lc_fsoc_sd
  
  cat(sprintf("\nCUTOFFS (Lean Control mean ± %.1f SD):\n", n_sd))
  cat(sprintf("  K2 cutoff: %.4f (above = High K2)\n", k2_cutoff))
  cat(sprintf("  FSOC cutoff: %.2f (below = Low FSOC)\n", fsoc_cutoff))
  
  dat_valid <- dat_valid %>%
    mutate(
      k2_level = ifelse(.data[[k2_var]] > k2_cutoff, "High K2", "Normal K2"),
      fsoc_level = ifelse(.data[[fsoc_var]] < fsoc_cutoff, "Low FSOC", "Normal FSOC"),
      phenotype_4 = case_when(
        k2_level == "High K2" & fsoc_level == "Low FSOC" ~ "High K2 / Low FSOC",
        k2_level == "High K2" & fsoc_level == "Normal FSOC" ~ "High K2 / Normal FSOC",
        k2_level == "Normal K2" & fsoc_level == "Low FSOC" ~ "Normal K2 / Low FSOC",
        k2_level == "Normal K2" & fsoc_level == "Normal FSOC" ~ "Normal K2 / Normal FSOC"
      ),
      discordant = ifelse(k2_level == "High K2" & fsoc_level == "Low FSOC",
                          "Discordant", "Concordant")
    )
  
  cat("\n\n--- PHENOTYPE DISTRIBUTION ---\n")
  print(table(dat_valid$phenotype_4))
  cat("\n--- BY GROUP ---\n")
  print(table(dat_valid$group, dat_valid$phenotype_4))
  cat("\n--- DISCORDANT SUMMARY ---\n")
  disc_summary <- dat_valid %>%
    group_by(group) %>%
    summarise(n = n(), n_discordant = sum(discordant == "Discordant"),
              pct_discordant = round(100 * n_discordant / n, 1), .groups = "drop")
  print(as.data.frame(disc_summary))
  
  n_disc <- sum(dat_valid$discordant == "Discordant")
  cat(sprintf("\nOVERALL: %d / %d (%.1f%%) are DISCORDANT\n", 
              n_disc, nrow(dat_valid), 100 * n_disc / nrow(dat_valid)))
  
  return(list(
    data = dat_valid,
    cutoffs = list(k2 = k2_cutoff, fsoc = fsoc_cutoff),
    lean_control_stats = list(k2_mean = lc_k2_mean, k2_sd = lc_k2_sd,
                              fsoc_mean = lc_fsoc_mean, fsoc_sd = lc_fsoc_sd),
    n_sd = n_sd
  ))
}

# ============================================================================
# STEP 3: CLINICAL COMPARISON - DISCORDANT vs CONCORDANT
# UPDATED: Removed lipids, Added DEXA and mGFR, Filter insufficient sample sizes
# ============================================================================

compare_clinical_characteristics <- function(pheno_results, linked_data, min_n = 2) {
  
  cat("\n========================================\n")
  cat("CLINICAL COMPARISON: DISCORDANT vs CONCORDANT\n")
  cat("========================================\n\n")
  
  dat <- pheno_results$data
  n_disc <- sum(dat$discordant == "Discordant")
  n_conc <- sum(dat$discordant == "Concordant")
  
  cat(sprintf("Discordant (High K2 / Low FSOC): n = %d\n", n_disc))
  cat(sprintf("Concordant: n = %d\n\n", n_conc))
  
  if (n_disc < 1) {
    cat("No discordant subjects found.\n")
    return(NULL)
  }
  
  # UPDATED: Removed lipids, Added DEXA and mGFR, use log10 ACR
  var_labels <- c(
    "age" = "Age (years)",
    "bmi" = "BMI (kg/m²)",
    "hba1c" = "HbA1c (%)",
    "eGFR_CKD_epi" = "eGFR (mL/min/1.73m²)",
    "gfr_raw_plasma" = "mGFR raw (mL/min)",
    "gfr_bsa_plasma" = "mGFR BSA (mL/min/1.73m²)",
    "log_acr" = "log₁₀(ACR) (mg/g)",
    "diabetes_duration" = "Diabetes Duration (years)",
    "sbp" = "Systolic BP (mmHg)",
    "dbp" = "Diastolic BP (mmHg)",
    "weight" = "Weight (kg)",
    "dexa_fat_kg" = "Fat Mass (kg)",
    "dexa_trunk_kg" = "Trunk Fat (kg)",
    "dexa_lean_kg" = "Lean Mass (kg)",
    "dexa_body_fat" = "Body Fat (%)",
    "avg_c_k2" = "Cortical K2 (s⁻¹)",
    "avg_m_k2" = "Medullary K2 (s⁻¹)",
    "fsoc_medulla" = "Medullary FSOC (s⁻¹)"
  )
  
  available_vars <- names(var_labels)[names(var_labels) %in% names(dat)]
  cat("Variables available for comparison:", length(available_vars), "\n\n")
  
  results <- data.frame()
  skipped_vars <- c()
  
  for (var in available_vars) {
    disc_vals <- dat[[var]][dat$discordant == "Discordant"]
    conc_vals <- dat[[var]][dat$discordant == "Concordant"]
    disc_vals <- disc_vals[!is.na(disc_vals)]
    conc_vals <- conc_vals[!is.na(conc_vals)]
    
    # UPDATED: Skip if either group has insufficient sample size
    if (length(disc_vals) < min_n | length(conc_vals) < min_n) {
      skipped_vars <- c(skipped_vars, 
                        sprintf("%s (Discordant n=%d, Concordant n=%d)", 
                                var_labels[var], length(disc_vals), length(conc_vals)))
      next
    }
    
    if (is.numeric(dat[[var]])) {
      p_val <- NA
      if (length(disc_vals) >= 2 & length(conc_vals) >= 2) {
        p_val <- tryCatch({
          wilcox.test(disc_vals, conc_vals, exact = FALSE)$p.value
        }, error = function(e) NA)
      }
      
      results <- rbind(results, data.frame(
        Variable = var_labels[var],
        Discordant_n = length(disc_vals),
        Discordant_mean = mean(disc_vals),
        Discordant_sd = sd(disc_vals),
        Concordant_n = length(conc_vals),
        Concordant_mean = mean(conc_vals),
        Concordant_sd = sd(conc_vals),
        Difference = mean(disc_vals) - mean(conc_vals),
        P_value = p_val,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Report skipped variables
  if (length(skipped_vars) > 0) {
    cat("\n--- VARIABLES SKIPPED (insufficient sample size, min n =", min_n, ") ---\n")
    for (sv in skipped_vars) {
      cat("  -", sv, "\n")
    }
    cat("\n")
  }
  
  if (nrow(results) > 0) {
    results <- results %>%
      mutate(
        Discordant = sprintf("%.2f ± %.2f", Discordant_mean, Discordant_sd),
        Concordant = sprintf("%.2f ± %.2f", Concordant_mean, Concordant_sd),
        P_value_fmt = ifelse(is.na(P_value), "-", 
                             ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))),
        Sig = case_when(is.na(P_value) ~ "", P_value < 0.01 ~ "**",
                        P_value < 0.05 ~ "*", P_value < 0.1 ~ "†", TRUE ~ "")
      ) %>%
      arrange(P_value)
    
    cat("\n--- CLINICAL CHARACTERISTICS TABLE ---\n\n")
    print_table <- results %>% select(Variable, Discordant, Concordant, P_value_fmt, Sig)
    print(as.data.frame(print_table), row.names = FALSE)
    cat("\n* p < 0.05, ** p < 0.01, † p < 0.10\n")
  }
  
  if ("sex" %in% names(dat)) {
    cat("\n\n--- SEX DISTRIBUTION ---\n")
    sex_table <- table(dat$discordant, dat$sex)
    print(sex_table)
    if (min(dim(sex_table)) >= 2 & all(sex_table > 0)) {
      fisher_p <- tryCatch({ fisher.test(sex_table)$p.value }, error = function(e) NA)
      if (!is.na(fisher_p)) cat(sprintf("Fisher's exact test p = %.4f\n", fisher_p))
    }
  }
  
  cat("\n\n--- GROUP DISTRIBUTION ---\n")
  print(table(dat$discordant, dat$group))
  
  return(results)
}

# ============================================================================
# STEP 6D: DISCORDANT vs ALL OTHERS COMBINED
# UPDATED: Removed lipids, Added DEXA and mGFR, Filter insufficient sample sizes
# ============================================================================

compare_discordant_vs_others <- function(pheno_results, min_n = 2) {
  
  cat("\n========================================\n")
  cat("DISCORDANT vs ALL OTHERS COMBINED\n")
  cat("========================================\n\n")
  
  dat <- pheno_results$data
  n_disc <- sum(dat$discordant == "Discordant")
  n_others <- sum(dat$discordant == "Concordant")
  
  cat(sprintf("Discordant (High K2 / Low FSOC): n = %d\n", n_disc))
  cat(sprintf("All Others Combined: n = %d\n", n_others))
  cat(sprintf("  - High K2 / Normal FSOC: n = %d\n", sum(dat$phenotype_4 == "High K2 / Normal FSOC")))
  cat(sprintf("  - Normal K2 / Low FSOC: n = %d\n", sum(dat$phenotype_4 == "Normal K2 / Low FSOC")))
  cat(sprintf("  - Normal K2 / Normal FSOC: n = %d\n\n", sum(dat$phenotype_4 == "Normal K2 / Normal FSOC")))
  
  if (n_disc < 1) {
    cat("No discordant subjects found.\n")
    return(NULL)
  }
  
  # UPDATED: Removed lipids, Added DEXA and mGFR, use log10 ACR
  var_labels <- c(
    "age" = "Age (years)",
    "bmi" = "BMI (kg/m²)",
    "hba1c" = "HbA1c (%)",
    "eGFR_CKD_epi" = "eGFR (mL/min/1.73m²)",
    "gfr_raw_plasma" = "mGFR raw (mL/min)",
    "gfr_bsa_plasma" = "mGFR BSA (mL/min/1.73m²)",
    "log_acr" = "log₁₀(ACR) (mg/g)",
    "diabetes_duration" = "Diabetes Duration (years)",
    "sbp" = "Systolic BP (mmHg)",
    "dbp" = "Diastolic BP (mmHg)",
    "weight" = "Weight (kg)",
    "dexa_fat_kg" = "Fat Mass (kg)",
    "dexa_trunk_kg" = "Trunk Fat (kg)",
    "dexa_lean_kg" = "Lean Mass (kg)",
    "dexa_body_fat" = "Body Fat (%)",
    "avg_c_k2" = "Cortical K2 (s⁻¹)",
    "avg_m_k2" = "Medullary K2 (s⁻¹)",
    "fsoc_medulla" = "Medullary FSOC (s⁻¹)"
  )
  
  available_vars <- names(var_labels)[names(var_labels) %in% names(dat)]
  results <- data.frame()
  skipped_vars <- c()
  
  for (var in available_vars) {
    if (!is.numeric(dat[[var]])) next
    
    disc_vals <- dat[[var]][dat$discordant == "Discordant"]
    other_vals <- dat[[var]][dat$discordant == "Concordant"]
    disc_vals <- disc_vals[!is.na(disc_vals)]
    other_vals <- other_vals[!is.na(other_vals)]
    
    # UPDATED: Skip if either group has insufficient sample size
    if (length(disc_vals) < min_n | length(other_vals) < min_n) {
      skipped_vars <- c(skipped_vars, 
                        sprintf("%s (Discordant n=%d, Others n=%d)", 
                                var_labels[var], length(disc_vals), length(other_vals)))
      next
    }
    
    p_val <- NA
    if (length(disc_vals) >= 2 & length(other_vals) >= 2) {
      p_val <- tryCatch({
        wilcox.test(disc_vals, other_vals, exact = FALSE)$p.value
      }, error = function(e) NA)
    }
    
    pooled_sd <- sqrt(((length(disc_vals)-1)*sd(disc_vals)^2 + 
                         (length(other_vals)-1)*sd(other_vals)^2) / 
                        (length(disc_vals) + length(other_vals) - 2))
    cohens_d <- ifelse(pooled_sd > 0, (mean(disc_vals) - mean(other_vals)) / pooled_sd, NA)
    
    results <- rbind(results, data.frame(
      Variable = var_labels[var],
      Discordant_n = length(disc_vals),
      Discordant_mean = mean(disc_vals),
      Discordant_sd = sd(disc_vals),
      Others_n = length(other_vals),
      Others_mean = mean(other_vals),
      Others_sd = sd(other_vals),
      Difference = mean(disc_vals) - mean(other_vals),
      Cohens_d = cohens_d,
      P_value = p_val,
      stringsAsFactors = FALSE
    ))
  }
  
  # Report skipped variables
  if (length(skipped_vars) > 0) {
    cat("\n--- VARIABLES SKIPPED (insufficient sample size, min n =", min_n, ") ---\n")
    for (sv in skipped_vars) {
      cat("  -", sv, "\n")
    }
    cat("\n")
  }
  
  if (nrow(results) > 0) {
    results <- results %>%
      mutate(
        Discordant = sprintf("%.2f ± %.2f", Discordant_mean, Discordant_sd),
        Others = sprintf("%.2f ± %.2f", Others_mean, Others_sd),
        Effect = sprintf("%.2f", Cohens_d),
        P_value_fmt = ifelse(is.na(P_value), "-", 
                             ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))),
        Sig = case_when(is.na(P_value) ~ "", P_value < 0.01 ~ "**",
                        P_value < 0.05 ~ "*", P_value < 0.1 ~ "†", TRUE ~ ""),
        Effect_size = case_when(is.na(Cohens_d) ~ "", abs(Cohens_d) >= 0.8 ~ "Large",
                                abs(Cohens_d) >= 0.5 ~ "Medium", abs(Cohens_d) >= 0.2 ~ "Small",
                                TRUE ~ "Negligible")
      ) %>%
      arrange(P_value)
    
    cat("\n--- DISCORDANT vs ALL OTHERS ---\n\n")
    print_tbl <- results %>%
      select(Variable, Discordant, Others, Effect, Effect_size, P_value_fmt, Sig)
    names(print_tbl) <- c("Variable", paste0("Discordant (n=", n_disc, ")"), 
                          paste0("All Others (n=", n_others, ")"), 
                          "Cohen's d", "Effect Size", "p-value", "")
    print(as.data.frame(print_tbl), row.names = FALSE)
    cat("\n* p < 0.05, ** p < 0.01, † p < 0.10\n")
    cat("Cohen's d: |d| >= 0.8 Large, >= 0.5 Medium, >= 0.2 Small\n")
  }
  
  return(results)
}

# ============================================================================
# STEP 6E: FOREST PLOT - DISCORDANT VS OTHERS
# UPDATED: Filter out variables with insufficient sample sizes
# ============================================================================

plot_discordant_vs_others_forest <- function(comparison_results, min_n_per_group = 2) {
  
  if (is.null(comparison_results) || nrow(comparison_results) == 0) {
    cat("No comparison results to plot.\n")
    return(NULL)
  }
  
  # Filter out rows with:
  # 1. Missing Cohen's d
  # 2. Insufficient sample size in either group
  # 3. Extreme effect sizes (likely artifacts)
  plot_data <- comparison_results %>%
    filter(!is.na(Cohens_d)) %>%
    filter(Discordant_n >= min_n_per_group & Others_n >= min_n_per_group) %>%
    filter(abs(Cohens_d) < 5) %>%
    mutate(
      Variable = factor(Variable, levels = rev(Variable)),
      significant = !is.na(P_value) & P_value < 0.05,
      direction = ifelse(Cohens_d > 0, "Higher in Discordant", "Lower in Discordant")
    )
  
  if (nrow(plot_data) == 0) {
    cat("No valid data for forest plot after filtering.\n")
    cat("Check that both groups have sufficient sample sizes (n >=", min_n_per_group, ")\n")
    return(NULL)
  }
  
  cat("\nForest plot includes", nrow(plot_data), "variables\n")
  cat("(Variables with n <", min_n_per_group, "in either group were excluded)\n\n")
  
  p <- ggplot(plot_data, aes(x = Cohens_d, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.8, -0.5, 0.5, 0.8), linetype = "dotted", color = "grey70") +
    geom_segment(aes(x = 0, xend = Cohens_d, y = Variable, yend = Variable,
                     color = direction), linewidth = 2, alpha = 0.8) +
    geom_point(aes(fill = direction, size = -log10(P_value + 0.001)), 
               shape = 21, color = "black", stroke = 0.8) +
    geom_text(aes(label = Sig, x = Cohens_d + sign(Cohens_d) * 0.15), 
              hjust = 0.5, vjust = 0.5, size = 5, fontface = "bold") +
    scale_fill_manual(values = c("Higher in Discordant" = "#E41A1C", 
                                 "Lower in Discordant" = "#4DAF4A"), name = "Direction") +
    scale_color_manual(values = c("Higher in Discordant" = "#E41A1C", 
                                  "Lower in Discordant" = "#4DAF4A"), guide = "none") +
    scale_size_continuous(range = c(4, 10), name = "-log10(p)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 11), plot.title = element_text(face = "bold")) +
    labs(x = "Cohen's d (Standardized Effect Size)", y = "",
         title = "Discordant vs All Others: Clinical Characteristics",
         subtitle = paste0("Positive values = higher in Discordant | Dotted lines at d = ±0.5, ±0.8\n",
                           "Variables with n < ", min_n_per_group, " in either group excluded"),
         caption = "* p < 0.05, ** p < 0.01, † p < 0.10") +
    coord_cartesian(xlim = c(-2.5, 2.5))
  
  return(p)
}

# ============================================================================
# STEP 4: VISUALIZATION - SCATTER PLOT
# ============================================================================

plot_phenotype_scatter_lc <- function(pheno_results, k2_var = "avg_c_k2", fsoc_var = "fsoc_medulla") {
  
  dat <- pheno_results$data
  k2_cut <- pheno_results$cutoffs$k2
  fsoc_cut <- pheno_results$cutoffs$fsoc
  lc_stats <- pheno_results$lean_control_stats
  n_sd <- pheno_results$n_sd
  
  n_disc <- sum(dat$discordant == "Discordant")
  n_total <- nrow(dat)
  
  x_range <- range(dat[[fsoc_var]], na.rm = TRUE)
  y_range <- range(dat[[k2_var]], na.rm = TRUE)
  x_pad <- diff(x_range) * 0.1
  y_pad <- diff(y_range) * 0.1
  
  p <- ggplot(dat, aes(x = .data[[fsoc_var]], y = .data[[k2_var]])) +
    annotate("rect", xmin = -Inf, xmax = fsoc_cut, ymin = k2_cut, ymax = Inf,
             fill = "#FF6B6B", alpha = 0.25) +
    annotate("rect", xmin = fsoc_cut, xmax = Inf, ymin = -Inf, ymax = k2_cut,
             fill = "#4ECDC4", alpha = 0.25) +
    annotate("rect", xmin = -Inf, xmax = fsoc_cut, ymin = -Inf, ymax = k2_cut,
             fill = "#FFE66D", alpha = 0.2) +
    annotate("rect", xmin = fsoc_cut, xmax = Inf, ymin = k2_cut, ymax = Inf,
             fill = "#95E1D3", alpha = 0.2) +
    geom_vline(xintercept = fsoc_cut, linetype = "dashed", color = "#E63946", linewidth = 1) +
    geom_hline(yintercept = k2_cut, linetype = "dashed", color = "#E63946", linewidth = 1) +
    geom_vline(xintercept = lc_stats$fsoc_mean, linetype = "dotted", color = "blue", linewidth = 0.7) +
    geom_hline(yintercept = lc_stats$k2_mean, linetype = "dotted", color = "blue", linewidth = 0.7) +
    geom_point(aes(fill = group, shape = discordant), size = 5, alpha = 0.85, color = "black", stroke = 0.5) +
    {if (use_repel) {
      geom_text_repel(aes(label = k2_record_id), size = 2.5, max.overlaps = 20,
                      box.padding = 0.3, point.padding = 0.2)
    } else {
      geom_text(aes(label = k2_record_id), hjust = -0.15, vjust = -0.4, size = 2.5)
    }} +
    annotate("text", x = x_range[1] + x_pad, y = y_range[2] - y_pad,
             label = paste0("DISCORDANT\n(High K2 / Low FSOC)\nn = ", 
                            sum(dat$phenotype_4 == "High K2 / Low FSOC")),
             fontface = "bold", size = 3.5, color = "#C92A2A", hjust = 0, vjust = 1) +
    annotate("text", x = x_range[2] - x_pad, y = y_range[1] + y_pad,
             label = paste0("NORMAL\n(Normal K2 / Normal FSOC)\nn = ", 
                            sum(dat$phenotype_4 == "Normal K2 / Normal FSOC")),
             fontface = "bold", size = 3.5, color = "#2B8A3E", hjust = 1, vjust = 0) +
    scale_fill_manual(values = group_colors, name = "Group") +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_shape_manual(values = c("Concordant" = 21, "Discordant" = 24), name = "Phenotype") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)), color = "none") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 10)) +
    labs(x = expression("Medullary FSOC (s"^-1*")"),
         y = expression("Cortical K2 - TCA Metabolism (s"^-1*")"),
         title = "K2 vs FSOC Phenotype Classification",
         subtitle = paste0("Cutoffs: Lean Control mean ± ", n_sd, " SD\n",
                           "K2 cutoff: ", round(k2_cut, 4), " | FSOC cutoff: ", round(fsoc_cut, 2), "\n",
                           "Discordant: ", n_disc, "/", n_total, " (", round(100*n_disc/n_total, 1), "%)"),
         caption = "Blue dotted lines = Lean Control means; Red dashed lines = cutoffs")
  
  return(p)
}

# ============================================================================
# STEP 5: VISUALIZATION - CLINICAL COMPARISON BAR PLOT
# ============================================================================

plot_clinical_comparison <- function(comparison_results, pheno_results, min_n_per_group = 2) {
  
  if (is.null(comparison_results) || nrow(comparison_results) == 0) {
    cat("No comparison results to plot.\n")
    return(NULL)
  }
  
  plot_data <- comparison_results %>%
    filter(!is.na(Difference)) %>%
    filter(Discordant_n >= min_n_per_group & Concordant_n >= min_n_per_group) %>%
    mutate(
      pooled_sd = sqrt((Discordant_sd^2 + Concordant_sd^2) / 2),
      effect_size = ifelse(pooled_sd > 0, Difference / pooled_sd, 0),
      Variable = factor(Variable, levels = rev(Variable)),
      significant = P_value < 0.05,
      direction = ifelse(Difference > 0, "Higher in Discordant", "Lower in Discordant")
    ) %>%
    filter(abs(effect_size) < 5)
  
  if (nrow(plot_data) == 0) {
    cat("No valid data for comparison plot.\n")
    return(NULL)
  }
  
  p <- ggplot(plot_data, aes(x = effect_size, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "grey70") +
    geom_segment(aes(x = 0, xend = effect_size, y = Variable, yend = Variable,
                     color = direction), linewidth = 1.5, alpha = 0.7) +
    geom_point(aes(fill = direction, size = -log10(P_value + 0.001)), 
               shape = 21, color = "black", stroke = 0.5) +
    geom_text(aes(label = Sig), hjust = -0.5, vjust = 0.5, size = 4) +
    scale_fill_manual(values = c("Higher in Discordant" = "#E63946", 
                                 "Lower in Discordant" = "#457B9D"), name = "Direction") +
    scale_color_manual(values = c("Higher in Discordant" = "#E63946", 
                                  "Lower in Discordant" = "#457B9D"), guide = "none") +
    scale_size_continuous(range = c(3, 8), name = "-log10(p)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right", panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10)) +
    labs(x = "Standardized Difference (Cohen's d)", y = "",
         title = "Clinical Characteristics: Discordant vs Concordant",
         subtitle = paste0("Positive = higher in Discordant group | Dotted lines = |d| = 0.5 (medium effect)\n",
                           "Variables with n < ", min_n_per_group, " in either group excluded"),
         caption = "* p < 0.05, ** p < 0.01, † p < 0.10") +
    coord_cartesian(xlim = c(-2.5, 2.5))
  
  return(p)
}

# ============================================================================
# STEP 6: BOXPLOTS OF KEY VARIABLES
# UPDATED: Added DEXA and mGFR, Filter insufficient sample sizes
# ============================================================================

plot_key_variables_boxplot <- function(pheno_results, min_n_per_group = 2) {
  
  dat <- pheno_results$data
  
  # UPDATED: Added DEXA and mGFR variables
  vars_to_plot <- c("avg_c_k2", "avg_m_k2", "fsoc_medulla", "bmi", "hba1c", 
                    "eGFR_CKD_epi", "gfr_bsa_plasma", "age",
                    "dexa_fat_kg", "dexa_lean_kg", "dexa_body_fat")
  vars_to_plot <- vars_to_plot[vars_to_plot %in% names(dat)]
  
  var_labels <- c(
    "avg_c_k2" = "Cortical K2",
    "avg_m_k2" = "Medullary K2", 
    "fsoc_medulla" = "Medullary FSOC",
    "bmi" = "BMI",
    "hba1c" = "HbA1c",
    "eGFR_CKD_epi" = "eGFR",
    "gfr_bsa_plasma" = "mGFR (BSA)",
    "age" = "Age",
    "dexa_fat_kg" = "Fat Mass (kg)",
    "dexa_lean_kg" = "Lean Mass (kg)",
    "dexa_body_fat" = "Body Fat (%)"
  )
  
  plots <- list()
  skipped_vars <- c()
  
  for (var in vars_to_plot) {
    # Check sample sizes per group
    n_disc <- sum(!is.na(dat[[var]]) & dat$discordant == "Discordant")
    n_conc <- sum(!is.na(dat[[var]]) & dat$discordant == "Concordant")
    
    if (n_disc < min_n_per_group | n_conc < min_n_per_group) {
      skipped_vars <- c(skipped_vars, sprintf("%s (Disc n=%d, Conc n=%d)", 
                                              var_labels[var], n_disc, n_conc))
      next
    }
    
    p_val <- tryCatch({
      wilcox.test(dat[[var]] ~ dat$discordant, exact = FALSE)$p.value
    }, error = function(e) NA)
    
    p_label <- ifelse(is.na(p_val), "", 
                      ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val)))
    
    p <- ggplot(dat, aes(x = discordant, y = .data[[var]], fill = discordant)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
      geom_jitter(aes(color = group), width = 0.15, size = 2.5, alpha = 0.8) +
      scale_fill_manual(values = c("Concordant" = "#4ECDC4", "Discordant" = "#FF6B6B")) +
      scale_color_manual(values = group_colors) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", axis.title.x = element_blank()) +
      labs(y = var_labels[var], subtitle = p_label)
    
    plots[[var]] <- p
  }
  
  if (length(skipped_vars) > 0) {
    cat("\n--- BOXPLOT VARIABLES SKIPPED (insufficient sample size, min n =", min_n_per_group, ") ---\n")
    for (sv in skipped_vars) {
      cat("  -", sv, "\n")
    }
    cat("\n")
  }
  
  if (length(plots) >= 2) {
    combined <- wrap_plots(plots, ncol = 3) +
      plot_annotation(
        title = "Key Variables by Phenotype (Discordant vs Concordant)",
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
    return(combined)
  }
  return(NULL)
}

# ============================================================================
# STEP 6B: COMPARISON ACROSS ALL 4 PHENOTYPE GROUPS
# UPDATED: Removed lipids, Added DEXA and mGFR, Filter insufficient sample sizes
# ============================================================================

compare_all_phenotypes <- function(pheno_results, min_n = 2) {
  
  cat("\n========================================\n")
  cat("COMPARISON ACROSS ALL 4 PHENOTYPE GROUPS\n")
  cat("========================================\n\n")
  
  dat <- pheno_results$data
  cat("PHENOTYPE GROUP SIZES:\n")
  pheno_counts <- dat %>% count(phenotype_4) %>% arrange(desc(n))
  print(as.data.frame(pheno_counts))
  
  # UPDATED: Removed lipids, Added DEXA and mGFR, use log10 ACR
  var_labels <- c(
    "age" = "Age (years)",
    "bmi" = "BMI (kg/m²)",
    "hba1c" = "HbA1c (%)",
    "eGFR_CKD_epi" = "eGFR (mL/min/1.73m²)",
    "gfr_raw_plasma" = "mGFR raw (mL/min)",
    "gfr_bsa_plasma" = "mGFR BSA (mL/min/1.73m²)",
    "log_acr" = "log₁₀(ACR) (mg/g)",
    "diabetes_duration" = "Diabetes Duration (years)",
    "sbp" = "Systolic BP (mmHg)",
    "dbp" = "Diastolic BP (mmHg)",
    "dexa_fat_kg" = "Fat Mass (kg)",
    "dexa_trunk_kg" = "Trunk Fat (kg)",
    "dexa_lean_kg" = "Lean Mass (kg)",
    "dexa_body_fat" = "Body Fat (%)",
    "avg_c_k2" = "Cortical K2 (s⁻¹)",
    "avg_m_k2" = "Medullary K2 (s⁻¹)",
    "fsoc_medulla" = "Medullary FSOC (s⁻¹)"
  )
  
  available_vars <- names(var_labels)[names(var_labels) %in% names(dat)]
  results <- data.frame()
  skipped_vars <- c()
  
  for (var in available_vars) {
    if (!is.numeric(dat[[var]])) next
    
    high_high <- dat[[var]][dat$phenotype_4 == "High K2 / Low FSOC"]
    high_norm <- dat[[var]][dat$phenotype_4 == "High K2 / Normal FSOC"]
    norm_low <- dat[[var]][dat$phenotype_4 == "Normal K2 / Low FSOC"]
    norm_norm <- dat[[var]][dat$phenotype_4 == "Normal K2 / Normal FSOC"]
    
    # Remove NAs
    high_high <- high_high[!is.na(high_high)]
    high_norm <- high_norm[!is.na(high_norm)]
    norm_low <- norm_low[!is.na(norm_low)]
    norm_norm <- norm_norm[!is.na(norm_norm)]
    
    # Check if at least 2 groups have sufficient data for KW test
    group_ns <- c(length(high_high), length(high_norm), length(norm_low), length(norm_norm))
    groups_with_data <- sum(group_ns >= min_n)
    
    if (groups_with_data < 2) {
      skipped_vars <- c(skipped_vars, 
                        sprintf("%s (groups with n>=%d: %d)", var_labels[var], min_n, groups_with_data))
      next
    }
    
    kw_p <- tryCatch({ kruskal.test(dat[[var]] ~ dat$phenotype_4)$p.value }, error = function(e) NA)
    
    disc_vs_norm_p <- tryCatch({
      disc <- dat[[var]][dat$phenotype_4 == "High K2 / Low FSOC"]
      norm <- dat[[var]][dat$phenotype_4 == "Normal K2 / Normal FSOC"]
      disc <- disc[!is.na(disc)]
      norm <- norm[!is.na(norm)]
      if (length(disc) >= min_n & length(norm) >= min_n) wilcox.test(disc, norm, exact = FALSE)$p.value else NA
    }, error = function(e) NA)
    
    results <- rbind(results, data.frame(
      Variable = var_labels[var],
      High_K2_Low_FSOC = sprintf("%.2f ± %.2f", mean(high_high, na.rm = TRUE), sd(high_high, na.rm = TRUE)),
      n1 = length(high_high),
      High_K2_Norm_FSOC = sprintf("%.2f ± %.2f", mean(high_norm, na.rm = TRUE), sd(high_norm, na.rm = TRUE)),
      n2 = length(high_norm),
      Norm_K2_Low_FSOC = sprintf("%.2f ± %.2f", mean(norm_low, na.rm = TRUE), sd(norm_low, na.rm = TRUE)),
      n3 = length(norm_low),
      Norm_K2_Norm_FSOC = sprintf("%.2f ± %.2f", mean(norm_norm, na.rm = TRUE), sd(norm_norm, na.rm = TRUE)),
      n4 = length(norm_norm),
      KW_p = kw_p,
      Disc_vs_Norm_p = disc_vs_norm_p,
      stringsAsFactors = FALSE
    ))
  }
  
  # Report skipped variables
  if (length(skipped_vars) > 0) {
    cat("\n--- VARIABLES SKIPPED (insufficient sample size, min n =", min_n, ") ---\n")
    for (sv in skipped_vars) {
      cat("  -", sv, "\n")
    }
    cat("\n")
  }
  
  if (nrow(results) > 0) {
    results <- results %>%
      mutate(
        KW_p_fmt = ifelse(is.na(KW_p), "-", ifelse(KW_p < 0.001, "<0.001", sprintf("%.3f", KW_p))),
        Disc_vs_Norm_fmt = ifelse(is.na(Disc_vs_Norm_p), "-",
                                  ifelse(Disc_vs_Norm_p < 0.001, "<0.001", sprintf("%.3f", Disc_vs_Norm_p))),
        Sig = case_when(is.na(KW_p) ~ "", KW_p < 0.01 ~ "**", KW_p < 0.05 ~ "*", KW_p < 0.1 ~ "†", TRUE ~ "")
      ) %>%
      arrange(KW_p)
    
    cat("\n--- CLINICAL CHARACTERISTICS BY PHENOTYPE GROUP ---\n\n")
    print_tbl <- results %>%
      select(Variable, High_K2_Low_FSOC, High_K2_Norm_FSOC, Norm_K2_Low_FSOC, Norm_K2_Norm_FSOC, KW_p_fmt, Sig)
    names(print_tbl) <- c("Variable", "High K2/Low FSOC\n(Discordant)", "High K2/Norm FSOC",
                          "Norm K2/Low FSOC", "Norm K2/Norm FSOC\n(Normal)", "KW p", "")
    print(as.data.frame(print_tbl), row.names = FALSE)
    cat("\n* p < 0.05, ** p < 0.01, † p < 0.10 (Kruskal-Wallis test)\n")
    
    cat("\n\n--- DISCORDANT vs NORMAL (pairwise comparison) ---\n\n")
    print_tbl2 <- results %>% select(Variable, High_K2_Low_FSOC, Norm_K2_Norm_FSOC, Disc_vs_Norm_fmt)
    names(print_tbl2) <- c("Variable", "Discordant", "Normal", "p-value")
    print(as.data.frame(print_tbl2), row.names = FALSE)
  }
  
  return(results)
}

# ============================================================================
# STEP 6C: BOXPLOTS BY ALL 4 PHENOTYPE GROUPS
# UPDATED: Added DEXA and mGFR, Filter insufficient sample sizes
# ============================================================================

plot_by_phenotype_group <- function(pheno_results, min_n_per_group = 2) {
  
  dat <- pheno_results$data
  
  # UPDATED: Added DEXA and mGFR
  vars_to_plot <- c("avg_c_k2", "avg_m_k2", "fsoc_medulla", "bmi", "hba1c", 
                    "eGFR_CKD_epi", "gfr_bsa_plasma",
                    "dexa_fat_kg", "dexa_lean_kg", "dexa_body_fat")
  vars_to_plot <- vars_to_plot[vars_to_plot %in% names(dat)]
  
  var_labels <- c(
    "avg_c_k2" = "Cortical K2 (s⁻¹)",
    "avg_m_k2" = "Medullary K2 (s⁻¹)", 
    "fsoc_medulla" = "Medullary FSOC (s⁻¹)",
    "bmi" = "BMI (kg/m²)",
    "hba1c" = "HbA1c (%)",
    "eGFR_CKD_epi" = "eGFR (mL/min/1.73m²)",
    "gfr_bsa_plasma" = "mGFR BSA (mL/min/1.73m²)",
    "dexa_fat_kg" = "Fat Mass (kg)",
    "dexa_lean_kg" = "Lean Mass (kg)",
    "dexa_body_fat" = "Body Fat (%)"
  )
  
  dat <- dat %>%
    mutate(phenotype_short = case_when(
      phenotype_4 == "High K2 / Low FSOC" ~ "High K2\nLow FSOC",
      phenotype_4 == "High K2 / Normal FSOC" ~ "High K2\nNorm FSOC",
      phenotype_4 == "Normal K2 / Low FSOC" ~ "Norm K2\nLow FSOC",
      phenotype_4 == "Normal K2 / Normal FSOC" ~ "Norm K2\nNorm FSOC"
    ))
  
  dat$phenotype_short <- factor(dat$phenotype_short, 
                                levels = c("High K2\nLow FSOC", "High K2\nNorm FSOC",
                                           "Norm K2\nLow FSOC", "Norm K2\nNorm FSOC"))
  
  plots <- list()
  skipped_vars <- c()
  
  for (var in vars_to_plot) {
    # Check sample sizes per phenotype group
    group_ns <- dat %>%
      filter(!is.na(.data[[var]])) %>%
      count(phenotype_4) %>%
      pull(n)
    
    groups_with_data <- sum(group_ns >= min_n_per_group)
    
    if (groups_with_data < 2) {
      skipped_vars <- c(skipped_vars, sprintf("%s (groups with n>=%d: %d)", 
                                              var_labels[var], min_n_per_group, groups_with_data))
      next
    }
    
    kw_p <- tryCatch({ kruskal.test(dat[[var]] ~ dat$phenotype_4)$p.value }, error = function(e) NA)
    p_label <- ifelse(is.na(kw_p), "", ifelse(kw_p < 0.001, "KW p < 0.001", sprintf("KW p = %.3f", kw_p)))
    
    p <- ggplot(dat, aes(x = phenotype_short, y = .data[[var]], fill = phenotype_short)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(aes(color = group), width = 0.2, size = 2, alpha = 0.8) +
      scale_fill_manual(values = c("High K2\nLow FSOC" = "#E41A1C", "High K2\nNorm FSOC" = "#FF7F00",
                                   "Norm K2\nLow FSOC" = "#FFFF33", "Norm K2\nNorm FSOC" = "#4DAF4A")) +
      scale_color_manual(values = group_colors) +
      theme_bw(base_size = 9) +
      theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 7)) +
      labs(y = var_labels[var], subtitle = p_label)
    
    plots[[var]] <- p
  }
  
  if (length(skipped_vars) > 0) {
    cat("\n--- PHENOTYPE BOXPLOT VARIABLES SKIPPED (insufficient sample size, min n =", min_n_per_group, ") ---\n")
    for (sv in skipped_vars) {
      cat("  -", sv, "\n")
    }
    cat("\n")
  }
  
  if (length(plots) >= 2) {
    combined <- wrap_plots(plots, ncol = 3) +
      plot_annotation(
        title = "Clinical Variables by K2-FSOC Phenotype Group",
        subtitle = "Colors: Red=T2D, Purple=T1D, Green=Lean Control, Blue=Obese Control",
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
    return(combined)
  }
  return(NULL)
}

# ============================================================================
# STEP 7: RUN COMPLETE ANALYSIS
# ============================================================================

cat("\n################################################################\n")
cat("# DISCORDANT ANALYSIS WITH LEAN CONTROL CUTOFFS               #\n")
cat("################################################################\n")

linked_data <- create_linked_k2_fsoc_dataset(dat)

if (!is.null(linked_data) && nrow(linked_data) > 0) {
  
  cat("\n\nLINKED DATA:\n")
  print(linked_data %>% 
          select(k2_record_id, fsoc_record_id, group, avg_c_k2, avg_m_k2, fsoc_medulla) %>%
          arrange(group))
  
  # Run with 1 SD cutoff
  cat("\n\n############################################\n")
  cat("# ANALYSIS WITH 1 SD CUTOFF               #\n")
  cat("############################################\n")
  
  pheno_1sd <- run_lean_control_phenotype(linked_data, n_sd = 1)
  
  p_scatter_1sd <- plot_phenotype_scatter_lc(pheno_1sd)
  print(p_scatter_1sd)
  ggsave(file.path(OUTPUT_DIR, "phenotype_scatter_LC_1sd.pdf"), p_scatter_1sd, width = 11, height = 9)
  
  comparison_1sd <- compare_clinical_characteristics(pheno_1sd, linked_data, min_n = 2)
  
  if (!is.null(comparison_1sd) && nrow(comparison_1sd) > 0) {
    p_clinical_1sd <- plot_clinical_comparison(comparison_1sd, pheno_1sd, min_n_per_group = 2)
    if (!is.null(p_clinical_1sd)) {
      print(p_clinical_1sd)
      ggsave(file.path(OUTPUT_DIR, "clinical_comparison_LC_1sd.pdf"), p_clinical_1sd, width = 10, height = 7)
    }
  }
  
  p_box_1sd <- plot_key_variables_boxplot(pheno_1sd, min_n_per_group = 2)
  if (!is.null(p_box_1sd)) {
    print(p_box_1sd)
    ggsave(file.path(OUTPUT_DIR, "boxplots_discordant_LC_1sd.pdf"), p_box_1sd, width = 12, height = 8)
  }
  
  # Compare all 4 phenotype groups
  cat("\n\n############################################\n")
  cat("# COMPARISON ACROSS ALL 4 PHENOTYPES      #\n")
  cat("############################################\n")
  
  comparison_4groups <- compare_all_phenotypes(pheno_1sd, min_n = 2)
  
  p_4groups <- plot_by_phenotype_group(pheno_1sd, min_n_per_group = 2)
  if (!is.null(p_4groups)) {
    print(p_4groups)
    ggsave(file.path(OUTPUT_DIR, "boxplots_4_phenotype_groups.pdf"), p_4groups, width = 12, height = 8)
  }
  
  # Discordant vs All Others Combined
  cat("\n\n############################################\n")
  cat("# DISCORDANT vs ALL OTHERS COMBINED       #\n")
  cat("############################################\n")
  
  comparison_disc_vs_others <- compare_discordant_vs_others(pheno_1sd, min_n = 2)
  
  p_forest <- plot_discordant_vs_others_forest(comparison_disc_vs_others, min_n_per_group = 2)
  if (!is.null(p_forest)) {
    print(p_forest)
    ggsave(file.path(OUTPUT_DIR, "forest_discordant_vs_others.pdf"), p_forest, width = 10, height = 7)
  }
  
  # Run with 0.5 SD cutoff
  cat("\n\n############################################\n")
  cat("# ANALYSIS WITH 0.5 SD CUTOFF             #\n")
  cat("############################################\n")
  
  pheno_05sd <- run_lean_control_phenotype(linked_data, n_sd = 0.5)
  
  p_scatter_05sd <- plot_phenotype_scatter_lc(pheno_05sd)
  print(p_scatter_05sd)
  ggsave(file.path(OUTPUT_DIR, "phenotype_scatter_LC_05sd.pdf"), p_scatter_05sd, width = 11, height = 9)
  
  comparison_05sd <- compare_clinical_characteristics(pheno_05sd, linked_data, min_n = 2)
  
  if (!is.null(comparison_05sd) && nrow(comparison_05sd) > 0) {
    p_clinical_05sd <- plot_clinical_comparison(comparison_05sd, pheno_05sd, min_n_per_group = 2)
    if (!is.null(p_clinical_05sd)) {
      print(p_clinical_05sd)
      ggsave(file.path(OUTPUT_DIR, "clinical_comparison_LC_05sd.pdf"), p_clinical_05sd, width = 10, height = 7)
    }
  }
  
  # Export results
  export_vars <- c("mrn", "k2_record_id", "fsoc_record_id", "group", 
                   "avg_c_k2", "avg_m_k2", "fsoc_medulla",
                   "k2_level", "fsoc_level", "phenotype_4", "discordant")
  export_vars <- export_vars[export_vars %in% names(pheno_1sd$data)]
  
  clinical_vars <- c("age", "sex", "bmi", "hba1c", "eGFR_CKD_epi", "gfr_bsa_plasma",
                     "acr_u", "diabetes_duration", "dexa_fat_kg", "dexa_lean_kg", "dexa_body_fat")
  clinical_vars <- clinical_vars[clinical_vars %in% names(pheno_1sd$data)]
  export_vars <- c(export_vars, clinical_vars)
  
  write.csv(pheno_1sd$data %>% select(all_of(unique(export_vars))),
            file.path(OUTPUT_DIR, "phenotypes_LC_cutoff.csv"), row.names = FALSE)
  
  if (!is.null(comparison_1sd)) {
    write.csv(comparison_1sd %>% select(Variable, Discordant, Concordant, P_value_fmt, Sig),
              file.path(OUTPUT_DIR, "clinical_comparison_table.csv"), row.names = FALSE)
  }
  
  if (!is.null(comparison_4groups)) {
    write.csv(comparison_4groups, file.path(OUTPUT_DIR, "clinical_comparison_4groups.csv"), row.names = FALSE)
  }
  
  if (!is.null(comparison_disc_vs_others)) {
    write.csv(comparison_disc_vs_others %>% 
                select(Variable, Discordant, Others, Effect, Effect_size, P_value_fmt, Sig),
              file.path(OUTPUT_DIR, "clinical_discordant_vs_others.csv"), row.names = FALSE)
  }
  
  cat("\n\n========================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("========================================\n\n")
  cat("Files saved to:", OUTPUT_DIR, "\n")
  cat("  - phenotype_scatter_LC_1sd.pdf\n")
  cat("  - phenotype_scatter_LC_05sd.pdf\n")
  cat("  - clinical_comparison_LC_1sd.pdf\n")
  cat("  - clinical_comparison_LC_05sd.pdf\n")
  cat("  - boxplots_discordant_LC_1sd.pdf\n")
  cat("  - boxplots_4_phenotype_groups.pdf\n")
  cat("  - forest_discordant_vs_others.pdf\n")
  cat("  - phenotypes_LC_cutoff.csv\n")
  cat("  - clinical_comparison_table.csv\n")
  cat("  - clinical_comparison_4groups.csv\n")
  cat("  - clinical_discordant_vs_others.csv\n")
  
} else {
  cat("\nERROR: Could not create linked dataset.\n")
}
