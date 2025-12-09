#### T2D vs. LC for Maninder


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

setwd('C:/Users/netio/Documents/UofW/Projects/Maninder_Data/Clinical_Biomarker_Associations_T2D/')
dir.create(getwd(), showWarnings = FALSE, recursive = TRUE)

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

# Continue with rest of analysis code...
# (Heatmaps, forest plots, etc. using the same structure as original)

cat("\n\n##########################################################")
cat("\n### ANALYSIS COMPLETE - T2D VS. CONTROL")
cat("\n##########################################################\n\n")

cat("Sample type distribution in analysis:\n")
sample_summary <- all_results %>%
  distinct(Group, Sample_Type, .keep_all = TRUE) %>%
  dplyr::select(Group, Sample_Type, N) %>%
  arrange(Sample_Type, Group)
print(sample_summary)

save.image("clinical_biomarker_analysis_T2D.RData")
cat("\nWorkspace saved: clinical_biomarker_analysis_T2D.RData\n")































