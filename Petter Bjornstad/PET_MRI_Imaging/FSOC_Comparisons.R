#####FSOC Analysis 



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




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )









# ============================================================================
# FSOC Analysis: Tubular Oxygen Consumption Across Disease Groups
# ============================================================================

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(emmeans)
library(car)
library(gtsummary)
library(corrplot)

# ============================================================================
# 1. DATA PREPARATION
# ============================================================================

# Your data is already loaded as 'dat'
# Data structure expected:
# - group (LC, OC, T1D, T2D)
# - sex (Male, Female)
# - age, weight
# - FSOC variables: fsoc_l_cortex, fsoc_l_kidney, fsoc_l_medulla,
#                   fsoc_r_cortex, fsoc_r_kidney, fsoc_r_medulla
# - DXA: dexa_fat_kg, dexa_trunk_kg, dexa_lean_kg, dexa_body_fat

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
      log_UACR = log(acr_u + 1)  # Add 1 to handle zeros
    )
  
  # If you have pre/post R2 values and need to calculate % change:
  # Uncomment and modify the following section
  
  # df <- df %>%
  #   mutate(
  #     # Whole-kidney FSOC (% change)
  #     left_wk_fsoc_pct = ((left_kidney_pre_R2 - left_kidney_post_R2) / 
  #                          left_kidney_pre_R2) * 100,
  #     right_wk_fsoc_pct = ((right_kidney_pre_R2 - right_kidney_post_R2) / 
  #                           right_kidney_pre_R2) * 100,
  #     whole_kidney_fsoc_pct = case_when(
  #       !is.na(left_wk_fsoc_pct) & !is.na(right_wk_fsoc_pct) ~ 
  #         (left_wk_fsoc_pct + right_wk_fsoc_pct) / 2,
  #       !is.na(left_wk_fsoc_pct) ~ left_wk_fsoc_pct,
  #       !is.na(right_wk_fsoc_pct) ~ right_wk_fsoc_pct,
  #       TRUE ~ NA_real_
  #     ),
  #     
  #     # Medullary FSOC (% change)
  #     left_med_fsoc_pct = ((left_medulla_pre_R2 - left_medulla_post_R2) / 
  #                           left_medulla_pre_R2) * 100,
  #     right_med_fsoc_pct = ((right_medulla_pre_R2 - right_medulla_post_R2) / 
  #                            right_medulla_pre_R2) * 100,
  #     medullary_fsoc_pct = case_when(
  #       !is.na(left_med_fsoc_pct) & !is.na(right_med_fsoc_pct) ~ 
  #         (left_med_fsoc_pct + right_med_fsoc_pct) / 2,
  #       !is.na(left_med_fsoc_pct) ~ left_med_fsoc_pct,
  #       !is.na(right_med_fsoc_pct) ~ right_med_fsoc_pct,
  #       TRUE ~ NA_real_
  #     )
  #   )
  
  return(df)
}

# Apply FSOC averaging
# df <- calculate_fsoc_averages(df)

# ============================================================================
# 3. OBJECTIVE 1.1: GROUP COMPARISONS
# ============================================================================

# Function for group comparisons
analyze_group_differences <- function(dat, outcome_var) {
  
  # Fit GLM adjusted for age, sex, and body weight
  formula_str <- paste(outcome_var, "~ group + age + sex + weight")
  model <- lm(as.formula(formula_str), data = dat)
  
  # ANOVA for disease group effect
  anova_result <- Anova(model, type = "III")
  
  # Pairwise comparisons with Tukey adjustment
  emm <- emmeans(model, specs = "group")
  pairwise <- pairs(emm, adjust = "tukey")
  
  # Return results
  list(
    model = model,
    anova = anova_result,
    emmeans = emm,
    pairwise = pairwise,
    summary = summary(model)
  )
}

# Run group comparisons for all FSOC endpoints
run_all_group_comparisons <- function(dat) {
  # Update to use only absolute FSOC values
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

# Execute group comparisons
# group_results <- run_all_group_comparisons(dat)

# ============================================================================
# 4. OBJECTIVE 1.2: SEX DIFFERENCES
# ============================================================================

# Fixed sex analysis function
analyze_sex_differences <- function(dat, outcome_var) {
  
  # Model with sex × disease group interaction
  formula_str <- paste(outcome_var, 
                       "~ group * sex + age + weight")
  
  # Check if sex has sufficient variation
  sex_table <- table(dat$sex, dat$group)
  
  # Only run if we have both sexes in multiple groups
  if (min(rowSums(sex_table > 0)) < 2) {
    cat("Warning: Insufficient sex variation for interaction analysis\n")
    return(list(
      interaction_model = NULL,
      interaction_anova = NULL,
      stratified = NULL,
      note = "Insufficient data for sex analysis"
    ))
  }
  
  model_interaction <- lm(as.formula(formula_str), data = dat)
  
  # Test interaction
  anova_interaction <- Anova(model_interaction, type = "III")
  
  # Stratified analyses by disease group (only for groups with both sexes)
  stratified_results <- dat %>%
    group_by(group) %>%
    filter(n_distinct(sex) > 1) %>%  # Only groups with both sexes
    do({
      if(nrow(.) > 5) {  # Need sufficient sample size
        model <- lm(as.formula(paste(outcome_var, "~ sex + age + weight")), 
                    data = .)
        tidy_result <- broom::tidy(model) %>%
          filter(grepl("sex", term, ignore.case = TRUE))
        tidy_result
      } else {
        data.frame()
      }
    })
  
  list(
    interaction_model = model_interaction,
    interaction_anova = anova_interaction,
    stratified = stratified_results
  )
}


# Run sex difference analyses
run_sex_analyses <- function(dat) {
  # Update to use only absolute FSOC values
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

# Execute sex analyses
# sex_results <- run_sex_analyses(dat)

# ============================================================================
# 5. OBJECTIVE 1.3: DXA ASSOCIATIONS
# ============================================================================

analyze_dxa_associations <- function(dat) {
  
  # Update to use only absolute FSOC values
  fsoc_endpoints <- c("whole_kidney_fsoc_abs", "medullary_fsoc_abs")
  
  dxa_vars <- c("dexa_fat_kg", "dexa_trunk_kg", "dexa_lean_kg", "dexa_body_fat")
  
  # Create matrix for correlation coefficients and p-values
  cor_matrix <- matrix(NA, nrow = length(fsoc_endpoints), 
                       ncol = length(dxa_vars))
  p_matrix <- matrix(NA, nrow = length(fsoc_endpoints), 
                     ncol = length(dxa_vars))
  
  rownames(cor_matrix) <- fsoc_endpoints
  colnames(cor_matrix) <- dxa_vars
  rownames(p_matrix) <- fsoc_endpoints
  colnames(p_matrix) <- dxa_vars
  
  # GLM for each combination
  detailed_results <- list()
  
  for (i in 1:length(fsoc_endpoints)) {
    for (j in 1:length(dxa_vars)) {
      formula_str <- paste(fsoc_endpoints[i], "~", dxa_vars[j], 
                           "+ age + sex + group")
      
      model <- lm(as.formula(formula_str), data = dat)
      coef_summary <- summary(model)$coefficients
      
      # Extract coefficient and p-value for DXA variable
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

# Visualize DXA associations
plot_dxa_heatmap <- function(dxa_results) {
  
  # Create annotation for significance
  p_mat <- dxa_results$p_value_matrix
  sig_labels <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
  sig_labels[p_mat < 0.001] <- "***"
  sig_labels[p_mat >= 0.001 & p_mat < 0.01] <- "**"
  sig_labels[p_mat >= 0.01 & p_mat < 0.05] <- "*"
  
  # Create heatmap
  pheatmap(dxa_results$coefficient_matrix,
           display_numbers = sig_labels,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "FSOC Associations with DXA Parameters\n(Adjusted for age, sex, disease group)",
           fontsize = 10,
           fontsize_number = 12)
}

# Execute DXA analyses
# dxa_results <- analyze_dxa_associations(dat)
# plot_dxa_heatmap(dxa_results)

# ============================================================================
# 6. DESCRIPTIVE STATISTICS
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

# Generate Table 1
# table1 <- create_table1(dat)
# print(table1)

# ============================================================================
# 7. VISUALIZATION: FSOC BY DISEASE GROUP
# ============================================================================

plot_fsoc_by_group <- function(dat) {
  
  fsoc_long <- dat %>%
    select(record_id, group, sex,
           whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    pivot_longer(cols = contains("fsoc"),
                 names_to = "fsoc_type",
                 values_to = "fsoc_value")
  
  ggplot(fsoc_long, aes(x = group, y = fsoc_value, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    facet_wrap(~ fsoc_type, scales = "free_y", ncol = 2) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Disease Group", y = "FSOC Value",
         title = "FSOC Endpoints by Disease Group") +
    scale_fill_brewer(palette = "Set2")
}

# Create plot
# plot_fsoc_by_group(dat)

# ============================================================================
# 8. EXPORT RESULTS
# ============================================================================

# Set output directory
OUTPUT_DIR <- "C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani"

# Function to export all results (DXA only, no clinical)
export_results <- function(group_results, sex_results, dxa_results, 
                           output_dir = OUTPUT_DIR) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Export group comparison results
  sink(file.path(output_dir, "group_comparisons.txt"))
  cat("========================================\n")
  cat("GROUP COMPARISONS: FSOC BY DISEASE GROUP\n")
  cat("========================================\n\n")
  print(group_results)
  sink()
  
  # Export sex difference results
  sink(file.path(output_dir, "sex_differences.txt"))
  cat("========================================\n")
  cat("SEX DIFFERENCES IN FSOC\n")
  cat("========================================\n\n")
  print(sex_results)
  sink()
  
  # Export heatmap
  pdf(file.path(output_dir, "dxa_heatmap.pdf"), width = 8, height = 6)
  plot_dxa_heatmap(dxa_results)
  dev.off()
  
  # Export FSOC boxplots (note: dat must be in global environment)
  pdf(file.path(output_dir, "fsoc_by_group.pdf"), width = 10, height = 6)
  print(plot_fsoc_by_group(dat))
  dev.off()
  
  # Export coefficient matrices
  write.csv(dxa_results$coefficient_matrix, 
            file.path(output_dir, "dxa_coefficients.csv"), row.names = TRUE)
  write.csv(dxa_results$p_value_matrix, 
            file.path(output_dir, "dxa_pvalues.csv"), row.names = TRUE)
  
  cat("Results exported to:", output_dir, "\n")
}

# Helper function to save individual plots
save_plot <- function(plot_obj, filename, width = 10, height = 6) {
  pdf(file.path(OUTPUT_DIR, filename), width = width, height = height)
  print(plot_obj)
  dev.off()
  cat("Saved:", filename, "\n")
}

# ============================================================================
# 9. MAIN EXECUTION WORKFLOW
# ============================================================================

# Output directory is set globally
OUTPUT_DIR <- "C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani"

 # Step 2: Calculate averaged FSOC endpoints
 dat <- calculate_fsoc_averages(dat)
 
 # Step 3: Run all analyses
 group_results <- run_all_group_comparisons(dat)
 sex_results <- run_sex_analyses(dat)
 dxa_results <- analyze_dxa_associations(dat)
 
 # Step 4: Create and save visualizations
 save_plot(plot_fsoc_by_group(dat), "fsoc_by_group.pdf")
 
 # Step 5: Generate descriptive table
 table1 <- create_table1(dat)
 print(table1)
 
 # Step 6: Save Table 1 as CSV
 table1_df <- as.data.frame(table1)
 write.csv(table1_df, file.path(OUTPUT_DIR, "table1_descriptives.csv"))
 
 # Step 7: Export all results (heatmaps and matrices)
 export_results(group_results, sex_results, dxa_results)
 
 cat("\n========================================\n")
 cat("All analyses complete!\n")
 cat("Results saved to:", OUTPUT_DIR, "\n")
 cat("========================================\n")

cat("FSOC Analysis Script Loaded Successfully!\n")
cat("Output directory set to:", OUTPUT_DIR, "\n")




























