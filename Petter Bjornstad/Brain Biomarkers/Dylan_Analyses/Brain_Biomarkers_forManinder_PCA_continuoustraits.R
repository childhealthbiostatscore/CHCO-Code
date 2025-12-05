# ============================================================
# INTERACTION ANALYSES - TEST IF EFFECTS DIFFER BY GROUP
# ============================================================

cat("\n\n##########################################################")
cat("\n### INTERACTION ANALYSES")
cat("\n### Testing if predictor effects differ between groups")
cat("\n##########################################################\n\n")

# Function to test interactions
test_interaction <- function(data, predictor, outcome, covariates = c("age", "sex", "bmi")) {
  
  # Need at least 2 groups to test interaction
  if(length(unique(data$group)) < 2) {
    return(NULL)
  }
  
  # Create analysis dataset
  analysis_vars <- c(predictor, outcome, "group", covariates)
  analysis_data <- data %>% 
    dplyr::select(all_of(analysis_vars)) %>%
    filter(complete.cases(.))
  
  if(nrow(analysis_data) < 20) return(NULL)
  
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
  
  # Extract interaction coefficients
  coef_summary <- summary(model_int)$coefficients
  interaction_terms <- grep(":", rownames(coef_summary), value = TRUE)
  
  if(length(interaction_terms) == 0) return(NULL)
  
  # Get main predictor effect (reference group)
  main_effect_row <- which(rownames(coef_summary) == predictor)
  if(length(main_effect_row) == 0) return(NULL)
  
  main_effect <- coef_summary[main_effect_row, "Estimate"]
  main_effect_p <- coef_summary[main_effect_row, "Pr(>|t|)"]
  
  # Extract interaction effects
  interaction_results <- data.frame()
  
  for(int_term in interaction_terms) {
    # Extract group name from interaction term
    group_name <- gsub(paste0(predictor, ":group"), "", int_term)
    
    int_row <- which(rownames(coef_summary) == int_term)
    
    result_row <- data.frame(
      Predictor = predictor,
      Outcome = outcome,
      Predictor_Used = pred_var_name,
      Outcome_Used = outcome_var_name,
      Predictor_Log_Transformed = pred_needs_log,
      Outcome_Log_Transformed = outcome_needs_log,
      Reference_Group = levels(factor(analysis_data$group))[1],
      Comparison_Group = group_name,
      N_Total = nrow(analysis_data),
      Main_Effect_Beta_Std = main_effect,
      Main_Effect_p = main_effect_p,
      Interaction_Beta_Std = coef_summary[int_row, "Estimate"],
      Interaction_SE = coef_summary[int_row, "Std. Error"],
      Interaction_t = coef_summary[int_row, "t value"],
      Interaction_p = coef_summary[int_row, "Pr(>|t|)"],
      Model_Interaction_p = interaction_p,
      R_squared_main = summary(model_main)$r.squared,
      R_squared_interaction = summary(model_int)$r.squared,
      stringsAsFactors = FALSE
    )
    
    interaction_results <- rbind(interaction_results, result_row)
  }
  
  return(interaction_results)
}

# Run interaction tests for predictors that showed group differences
cat("Testing interactions for clinical predictors with brain biomarkers...\n\n")

# Get predictors that were significant in at least one group
predictors_to_test <- significant_results %>%
  filter(Model_Type == "Adjusted") %>%
  group_by(Predictor, Outcome) %>%
  filter(n() >= 2) %>%  # Present in at least 2 groups
  ungroup() %>%
  distinct(Predictor, Outcome)

cat("Testing", nrow(predictors_to_test), "predictor-outcome pairs for interactions\n\n")

interaction_results_all <- data.frame()

if(nrow(predictors_to_test) > 0) {
  
  for(i in 1:nrow(predictors_to_test)) {
    
    predictor <- predictors_to_test$Predictor[i]
    outcome <- predictors_to_test$Outcome[i]
    
    # Test interaction using all participants data with group variable
    int_result <- test_interaction(
      data = dat_analysis,
      predictor = predictor,
      outcome = outcome,
      covariates = c("age", "sex", "bmi")
    )
    
    if(!is.null(int_result)) {
      interaction_results_all <- rbind(interaction_results_all, int_result)
    }
  }
}

cat("Completed", nrow(interaction_results_all), "interaction tests\n\n")

# Add significance flags
if(nrow(interaction_results_all) > 0) {
  interaction_results_all <- interaction_results_all %>%
    mutate(
      Interaction_Significant = Interaction_p < 0.05,
      Model_Significant = Model_Interaction_p < 0.05,
      Sig_Flag = case_when(
        Interaction_p < 0.001 ~ "***",
        Interaction_p < 0.01 ~ "**",
        Interaction_p < 0.05 ~ "*",
        Interaction_p < 0.10 ~ ".",
        TRUE ~ ""
      )
    )
  
  # Save results
  write.csv(interaction_results_all, "interaction_analysis_results.csv", row.names = FALSE)
  
  # Identify significant interactions
  significant_interactions <- interaction_results_all %>%
    filter(Interaction_p < 0.05) %>%
    arrange(Interaction_p)
  
  write.csv(significant_interactions, "significant_interactions.csv", row.names = FALSE)
  
  cat("Significant interactions (p < 0.05):", nrow(significant_interactions), "\n")
  
  if(nrow(significant_interactions) > 0) {
    cat("\n=== TOP 20 SIGNIFICANT INTERACTIONS ===\n")
    cat("(Interaction term tests if the predictor effect differs between groups)\n\n")
    
    top_interactions <- significant_interactions %>%
      head(20) %>%
      dplyr::select(Predictor, Outcome, Reference_Group, Comparison_Group, 
                    Main_Effect_Beta_Std, Interaction_Beta_Std, Interaction_p, 
                    R_squared_interaction)
    
    print(kable(top_interactions, digits = c(0, 0, 0, 0, 3, 3, 4, 3), format = "simple"))
    
    cat("\nInterpretation:\n")
    cat("- Main_Effect_Beta_Std: Effect of predictor in the reference group\n")
    cat("- Interaction_Beta_Std: Additional effect in comparison group\n")
    cat("- Total effect in comparison group = Main_Effect + Interaction\n\n")
  }
  
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
  
} else {
  cat("No predictor-outcome pairs available for interaction testing\n")
}

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




