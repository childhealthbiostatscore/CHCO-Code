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
      study == "RENAL-HEIRitage" ~ "serum",
      study == "RENAL-HEIR" ~ "serum",
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


dup_ids <- dat_baseline$record_id[which(dat_baseline$mrn %in% dat_baseline$mrn[which(duplicated(dat_baseline$mrn))])]


dat_baseline_dup <- dat_baseline %>% filter(record_id %in% dup_ids) %>% 
  dplyr::select(record_id, mrn, group, study, ab40_avg_conc, ab42_avg_conc, adipose_ir, search_eis) %>%
  arrange(by = mrn)



library(ggplot2)
library(tidyr)

# Define your variables
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Reshape data to long format
dat_long <- dat_baseline %>%
  pivot_longer(cols = all_of(qx_var),
               names_to = "biomarker",
               values_to = "concentration")

# Create the plot
p <- ggplot(dat_long, aes(x = group, y = concentration, fill = sample_type)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("plasma" = "#56B4E9", "serum" = "#E69F00")) +
  theme_minimal() +
  labs(x = "Group", 
       y = "Concentration",
       fill = "Sample Type") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10))

# Save as JPEG with high resolution to specified path
ggsave("C:/Users/netio/Downloads/violin_plot_biomarkers.jpeg", 
       plot = p,
       width = 14, 
       height = 10, 
       dpi = 300,
       units = "in")



library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# Define your variables
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Filter for lean controls and reshape data to long format
dat_long <- dat_baseline %>%
  filter(group == "Lean Control") %>%
  pivot_longer(cols = all_of(qx_var),
               names_to = "biomarker",
               values_to = "concentration")

# Create the plot
p <- ggplot(dat_long, aes(x = sample_type, y = concentration, fill = sample_type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("plasma" = "#56B4E9", "serum" = "#E69F00")) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("plasma", "serum")),
                     label = "p.signif",
                     tip.length = 0.02) +
  theme_minimal() +
  labs(x = "Sample Type", 
       y = "Concentration",
       fill = "Sample Type",
       title = "Lean Controls - Biomarker Concentrations by Sample Type") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10))

# Save as JPEG with high resolution to specified path
ggsave("C:/Users/netio/Downloads/violin_plot_lean_controls_biomarkers.jpeg", 
       plot = p,
       width = 14, 
       height = 8, 
       dpi = 300,
       units = "in")










library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# Define your variables
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Filter for serum only and reshape data to long format
dat_long <- dat_baseline %>%
  filter(sample_type == "serum") %>%
  pivot_longer(cols = all_of(qx_var),
               names_to = "biomarker",
               values_to = "concentration") %>%
  filter(!is.na(concentration))

# Simplified comparisons - only compare each group to Lean Control
my_comparisons <- list(
  c("Lean Control", "Obese Control"),
  c("Lean Control", "Type 2 Diabetes")
)

# Create the plot with custom colors
p <- ggplot(dat_long, aes(x = group, y = concentration, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Lean Control" = "#56B4E9",
                               "Obese Control" = "#F0E442", 
                               "Type 2 Diabetes" = "#D55E00")) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     method.args = list(exact = FALSE),
                     tip.length = 0.01) +
  theme_minimal() +
  labs(x = "Group", 
       y = "Concentration",
       fill = "Group",
       title = "Serum Biomarker Concentrations by Group") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))

# Save as JPEG with high resolution to specified path
ggsave("C:/Users/netio/Downloads/violin_plot_serum_groups_biomarkers.jpeg", 
       plot = p,
       width = 16, 
       height = 10, 
       dpi = 300,
       units = "in")








### eGFR adjustment

library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(ggpubr)

# Define your variables
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Filter for serum only
dat_serum <- dat_baseline %>%
  filter(sample_type == "serum") %>%
  filter(!is.na(eGFR_CKD_epi))

# Initialize lists to store results
emm_results <- list()
pvalue_results <- list()
plot_data_list <- list()
adjusted_values_list <- list()

# Loop through each biomarker
for (biomarker in qx_var) {
  
  cat("\n========================================\n")
  cat("Biomarker:", biomarker, "\n")
  cat("========================================\n")
  
  # Create formula
  formula_str <- paste(biomarker, "~ group + eGFR_CKD_epi")
  
  # Fit model
  model <- lm(as.formula(formula_str), data = dat_serum)
  
  # Get adjusted means
  emm <- emmeans(model, ~ group)
  
  # Store EMM results
  emm_results[[biomarker]] <- summary(emm)
  
  # Get pairwise comparisons
  pairs_result <- pairs(emm, adjust = "none")  # uncorrected p-values
  pvalue_results[[biomarker]] <- summary(pairs_result)
  
  # Print results
  cat("\nAdjusted Means (controlling for eGFR):\n")
  print(summary(emm))
  
  cat("\nPairwise Comparisons (uncorrected p-values):\n")
  print(summary(pairs_result))
  
  # Prepare data for plotting adjusted means
  emm_df <- as.data.frame(emm)
  emm_df$biomarker <- biomarker
  plot_data_list[[biomarker]] <- emm_df
  
  # Calculate adjusted values for violin plots
  # These are the observed values adjusted to the mean eGFR
  dat_serum_subset <- dat_serum %>%
    filter(!is.na(.data[[biomarker]]))
  
  # Get the mean eGFR
  mean_egfr <- mean(dat_serum_subset$eGFR_CKD_epi, na.rm = TRUE)
  
  # Fit full model
  model_full <- lm(as.formula(paste(biomarker, "~ group + eGFR_CKD_epi")), 
                   data = dat_serum_subset)
  
  # Calculate adjusted values: observed - eGFR effect + eGFR effect at mean
  egfr_coef <- coef(model_full)["eGFR_CKD_epi"]
  adjusted_values <- dat_serum_subset[[biomarker]] - 
    egfr_coef * (dat_serum_subset$eGFR_CKD_epi - mean_egfr)
  
  adjusted_df <- data.frame(
    group = dat_serum_subset$group,
    adjusted_value = adjusted_values,
    biomarker = biomarker
  )
  
  adjusted_values_list[[biomarker]] <- adjusted_df
}

# Combine all EMM data for plotting
plot_data_all <- bind_rows(plot_data_list)

# Combine all adjusted values data
adjusted_data_all <- bind_rows(adjusted_values_list)

# Create simplified comparisons for annotation
my_comparisons <- list(
  c("Lean Control", "Obese Control"),
  c("Lean Control", "Type 2 Diabetes")
)

# Create violin plot with adjusted values (eGFR-adjusted)
p_violin_adjusted <- ggplot(adjusted_data_all, aes(x = group, y = adjusted_value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Lean Control" = "#56B4E9",
                               "Obese Control" = "#F0E442", 
                               "Type 2 Diabetes" = "#D55E00")) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     method.args = list(exact = FALSE),
                     tip.length = 0.01) +
  theme_minimal() +
  labs(x = "Group", 
       y = "Adjusted Concentration (eGFR-Adjusted)",
       fill = "Group",
       title = "Serum Biomarker Concentrations (Adjusted for eGFR)") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))

# Save violin plot
ggsave("C:/Users/netio/Downloads/violin_plot_adjusted_serum_biomarkers.jpeg", 
       plot = p_violin_adjusted,
       width = 16, 
       height = 10, 
       dpi = 300,
       units = "in")

# Create bar plot with LSMeans/EMMs
p_adjusted_bars <- ggplot(plot_data_all, aes(x = group, y = emmean, fill = group)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Lean Control" = "#56B4E9",
                               "Obese Control" = "#F0E442", 
                               "Type 2 Diabetes" = "#D55E00")) +
  theme_minimal() +
  labs(x = "Group", 
       y = "Adjusted Mean Concentration (LSMeans, eGFR-adjusted)",
       fill = "Group",
       title = "Serum Biomarker Concentrations (Estimated Marginal Means, Adjusted for eGFR)") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))

# Save bar plot
ggsave("C:/Users/netio/Downloads/adjusted_means_bars_serum_biomarkers.jpeg", 
       plot = p_adjusted_bars,
       width = 16, 
       height = 10, 
       dpi = 300,
       units = "in")

# Create summary tables
cat("\n\n========================================")
cat("\n=== SUMMARY TABLE: P-VALUES (SERUM) ===")
cat("\n========================================\n\n")

# Format p-value table
pvalue_table <- bind_rows(lapply(names(pvalue_results), function(bio) {
  df <- as.data.frame(pvalue_results[[bio]])
  df$biomarker <- bio
  df
}))

# Print p-value table
pvalue_table %>%
  select(biomarker, contrast, estimate, SE, p.value) %>%
  mutate(p.value = format.pval(p.value, digits = 3)) %>%
  kable(format = "simple", 
        col.names = c("Biomarker", "Contrast", "Difference", "SE", "P-value"),
        caption = "Pairwise Comparisons - Serum (Adjusted for eGFR)") %>%
  print()

# Create adjusted means table
cat("\n\n========================================")
cat("\n=== SUMMARY TABLE: ADJUSTED MEANS (SERUM) ===")
cat("\n========================================\n\n")

emm_table <- bind_rows(lapply(names(emm_results), function(bio) {
  df <- as.data.frame(emm_results[[bio]])
  df$biomarker <- bio
  df
}))

# Print adjusted means table
emm_table %>%
  select(biomarker, group, emmean, SE, lower.CL, upper.CL) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  kable(format = "simple",
        col.names = c("Biomarker", "Group", "Adjusted Mean", "SE", "Lower CI", "Upper CI"),
        caption = "Estimated Marginal Means - Serum (Adjusted for eGFR)") %>%
  print()

# Save tables as CSV
write.csv(pvalue_table, "C:/Users/netio/Downloads/pvalue_table_adjusted_eGFR_serum.csv", row.names = FALSE)
write.csv(emm_table, "C:/Users/netio/Downloads/adjusted_means_table_eGFR_serum.csv", row.names = FALSE)

cat("\n\nTables and plots saved to Downloads folder!\n")
cat("- violin_plot_adjusted_serum_biomarkers.jpeg (violin plots with adjusted values)\n")
cat("- adjusted_means_bars_serum_biomarkers.jpeg (bar plot with LSMeans/EMMs)\n")
cat("- pvalue_table_adjusted_eGFR_serum.csv\n")
cat("- adjusted_means_table_eGFR_serum.csv\n")



























#####Plasma analysis 




library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# Define your variables
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Filter for plasma only and reshape data to long format
dat_long_plasma <- dat_baseline %>%
  filter(sample_type == "plasma") %>%
  pivot_longer(cols = all_of(qx_var),
               names_to = "biomarker",
               values_to = "concentration") %>%
  filter(!is.na(concentration))

# Simplified comparisons - compare T1D and PKD to Lean Control
my_comparisons_plasma <- list(
  c("Lean Control", "Type 1 Diabetes"),
  c("Lean Control", "PKD")
)

# Create the plot with custom colors
p_plasma <- ggplot(dat_long_plasma, aes(x = group, y = concentration, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Lean Control" = "#56B4E9",
                               "Type 1 Diabetes" = "#E69F00", 
                               "PKD" = "#009E73")) +
  stat_compare_means(comparisons = my_comparisons_plasma,
                     method = "wilcox.test",
                     label = "p.format",
                     method.args = list(exact = FALSE),
                     tip.length = 0.01) +
  theme_minimal() +
  labs(x = "Group", 
       y = "Concentration",
       fill = "Group",
       title = "Plasma Biomarker Concentrations by Group") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))

# Save as JPEG with high resolution to specified path
ggsave("C:/Users/netio/Downloads/violin_plot_plasma_groups_biomarkers.jpeg", 
       plot = p_plasma,
       width = 16, 
       height = 10, 
       dpi = 300,
       units = "in")







#Emmeans corrections

library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(ggpubr)

# Define your variables
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Filter for plasma only
dat_plasma <- dat_baseline %>%
  filter(sample_type == "plasma") %>%
  filter(!is.na(eGFR_CKD_epi))

# Initialize lists to store results
emm_results_plasma <- list()
pvalue_results_plasma <- list()
plot_data_list_plasma <- list()
adjusted_values_list_plasma <- list()

# Loop through each biomarker
for (biomarker in qx_var) {
  
  cat("\n========================================\n")
  cat("Biomarker:", biomarker, "\n")
  cat("========================================\n")
  
  # Create formula
  formula_str <- paste(biomarker, "~ group + eGFR_CKD_epi")
  
  # Fit model
  model <- lm(as.formula(formula_str), data = dat_plasma)
  
  # Get adjusted means
  emm <- emmeans(model, ~ group)
  
  # Store EMM results
  emm_results_plasma[[biomarker]] <- summary(emm)
  
  # Get pairwise comparisons
  pairs_result <- pairs(emm, adjust = "none")  # uncorrected p-values
  pvalue_results_plasma[[biomarker]] <- summary(pairs_result)
  
  # Print results
  cat("\nAdjusted Means (controlling for eGFR):\n")
  print(summary(emm))
  
  cat("\nPairwise Comparisons (uncorrected p-values):\n")
  print(summary(pairs_result))
  
  # Prepare data for plotting adjusted means
  emm_df <- as.data.frame(emm)
  emm_df$biomarker <- biomarker
  plot_data_list_plasma[[biomarker]] <- emm_df
  
  # Calculate adjusted values for violin plots
  # These are the observed values adjusted to the mean eGFR
  dat_plasma_subset <- dat_plasma %>%
    filter(!is.na(.data[[biomarker]]))
  
  # Get the mean eGFR
  mean_egfr <- mean(dat_plasma_subset$eGFR_CKD_epi, na.rm = TRUE)
  
  # Fit full model
  model_full <- lm(as.formula(paste(biomarker, "~ group + eGFR_CKD_epi")), 
                   data = dat_plasma_subset)
  
  # Calculate adjusted values: observed - eGFR effect + eGFR effect at mean
  egfr_coef <- coef(model_full)["eGFR_CKD_epi"]
  adjusted_values <- dat_plasma_subset[[biomarker]] - 
    egfr_coef * (dat_plasma_subset$eGFR_CKD_epi - mean_egfr)
  
  adjusted_df <- data.frame(
    group = dat_plasma_subset$group,
    adjusted_value = adjusted_values,
    biomarker = biomarker
  )
  
  adjusted_values_list_plasma[[biomarker]] <- adjusted_df
}

# Combine all EMM data for plotting
plot_data_all_plasma <- bind_rows(plot_data_list_plasma)

# Combine all adjusted values data
adjusted_data_all_plasma <- bind_rows(adjusted_values_list_plasma)

# Create comparisons including PKD vs T1D
my_comparisons_plasma <- list(
  c("Lean Control", "Type 1 Diabetes"),
  c("Lean Control", "PKD"),
  c("PKD", "Type 1 Diabetes")
)

# Create violin plot with adjusted values (eGFR-adjusted)
p_violin_adjusted_plasma <- ggplot(adjusted_data_all_plasma, aes(x = group, y = adjusted_value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Lean Control" = "#56B4E9",
                               "Type 1 Diabetes" = "#E69F00", 
                               "PKD" = "#009E73")) +
  stat_compare_means(comparisons = my_comparisons_plasma,
                     method = "wilcox.test",
                     label = "p.format",
                     method.args = list(exact = FALSE),
                     tip.length = 0.01,
                     step.increase = 0.08) +
  theme_minimal() +
  labs(x = "Group", 
       y = "Adjusted Concentration (eGFR-Adjusted)",
       fill = "Group",
       title = "Plasma Biomarker Concentrations (Adjusted for eGFR)") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))

# Save violin plot
ggsave("C:/Users/netio/Downloads/violin_plot_adjusted_plasma_biomarkers.jpeg", 
       plot = p_violin_adjusted_plasma,
       width = 16, 
       height = 10, 
       dpi = 300,
       units = "in")

# Create bar plot with LSMeans/EMMs
p_adjusted_bars_plasma <- ggplot(plot_data_all_plasma, aes(x = group, y = emmean, fill = group)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Lean Control" = "#56B4E9",
                               "Type 1 Diabetes" = "#E69F00", 
                               "PKD" = "#009E73")) +
  theme_minimal() +
  labs(x = "Group", 
       y = "Adjusted Mean Concentration (LSMeans, eGFR-adjusted)",
       fill = "Group",
       title = "Plasma Biomarker Concentrations (Estimated Marginal Means, Adjusted for eGFR)") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))

# Save bar plot
ggsave("C:/Users/netio/Downloads/adjusted_means_bars_plasma_biomarkers.jpeg", 
       plot = p_adjusted_bars_plasma,
       width = 16, 
       height = 10, 
       dpi = 300,
       units = "in")

# Create summary tables
cat("\n\n========================================")
cat("\n=== SUMMARY TABLE: P-VALUES (PLASMA) ===")
cat("\n========================================\n\n")

# Format p-value table
pvalue_table_plasma <- bind_rows(lapply(names(pvalue_results_plasma), function(bio) {
  df <- as.data.frame(pvalue_results_plasma[[bio]])
  df$biomarker <- bio
  df
}))

# Print p-value table
pvalue_table_plasma %>%
  select(biomarker, contrast, estimate, SE, p.value) %>%
  mutate(p.value = format.pval(p.value, digits = 3)) %>%
  kable(format = "simple", 
        col.names = c("Biomarker", "Contrast", "Difference", "SE", "P-value"),
        caption = "Pairwise Comparisons - Plasma (Adjusted for eGFR)") %>%
  print()

# Create adjusted means table
cat("\n\n========================================")
cat("\n=== SUMMARY TABLE: ADJUSTED MEANS (PLASMA) ===")
cat("\n========================================\n\n")

emm_table_plasma <- bind_rows(lapply(names(emm_results_plasma), function(bio) {
  df <- as.data.frame(emm_results_plasma[[bio]])
  df$biomarker <- bio
  df
}))

# Print adjusted means table
emm_table_plasma %>%
  select(biomarker, group, emmean, SE, lower.CL, upper.CL) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  kable(format = "simple",
        col.names = c("Biomarker", "Group", "Adjusted Mean", "SE", "Lower CI", "Upper CI"),
        caption = "Estimated Marginal Means - Plasma (Adjusted for eGFR)") %>%
  print()

# Save tables as CSV
write.csv(pvalue_table_plasma, "C:/Users/netio/Downloads/pvalue_table_adjusted_eGFR_plasma.csv", row.names = FALSE)
write.csv(emm_table_plasma, "C:/Users/netio/Downloads/adjusted_means_table_eGFR_plasma.csv", row.names = FALSE)

cat("\n\nTables and plots saved to Downloads folder!\n")
cat("- violin_plot_adjusted_plasma_biomarkers.jpeg (violin plots with adjusted values)\n")
cat("- adjusted_means_bars_plasma_biomarkers.jpeg (bar plot with LSMeans/EMMs)\n")
cat("- pvalue_table_adjusted_eGFR_plasma.csv\n")
cat("- adjusted_means_table_eGFR_plasma.csv\n")










#### demographics 
library(dplyr)
library(knitr)
library(kableExtra)

# ============================================================
# SERUM DEMOGRAPHICS TABLE
# ============================================================

cat("\n========================================")
cat("\n=== SERUM DEMOGRAPHICS TABLE ===")
cat("\n========================================\n\n")

# Filter for serum samples
dat_serum_demo <- dat_baseline %>%
  filter(sample_type == "serum")

# Create demographics table for serum
serum_demo <- dat_serum_demo %>%
  group_by(group) %>%
  summarise(
    N = n(),
    Age_mean = mean(age, na.rm = TRUE),
    Age_sd = sd(age, na.rm = TRUE),
    Sex_Female_n = sum(sex == "Female", na.rm = TRUE),
    Sex_Female_pct = mean(sex == "Female", na.rm = TRUE) * 100,
    BMI_mean = mean(bmi, na.rm = TRUE),
    BMI_sd = sd(bmi, na.rm = TRUE),
    HbA1c_mean = mean(hba1c, na.rm = TRUE),
    HbA1c_sd = sd(hba1c, na.rm = TRUE),
    eGFR_mean = mean(eGFR_CKD_epi, na.rm = TRUE),
    eGFR_sd = sd(eGFR_CKD_epi, na.rm = TRUE)
  ) %>%
  mutate(
    Age = sprintf("%.1f ± %.1f", Age_mean, Age_sd),
    Sex_Female = sprintf("%d (%.1f%%)", Sex_Female_n, Sex_Female_pct),
    BMI = sprintf("%.1f ± %.1f", BMI_mean, BMI_sd),
    HbA1c = sprintf("%.1f ± %.1f", HbA1c_mean, HbA1c_sd),
    eGFR = sprintf("%.1f ± %.1f", eGFR_mean, eGFR_sd)
  ) %>%
  select(group, N, Age, Sex_Female, BMI, HbA1c, eGFR)

# Print serum demographics table
print(serum_demo %>%
        kable(format = "simple",
              col.names = c("Group", "N", "Age (years)", "Female n (%)", "BMI (kg/m²)", "HbA1c (%)", "eGFR (mL/min/1.73m²)"),
              caption = "Serum Sample Demographics") %>%
        kable_styling())

# Save serum demographics table
write.csv(serum_demo, "C:/Users/netio/Downloads/serum_demographics_table.csv", row.names = FALSE)


# ============================================================
# PLASMA DEMOGRAPHICS TABLE
# ============================================================

cat("\n\n========================================")
cat("\n=== PLASMA DEMOGRAPHICS TABLE ===")
cat("\n========================================\n\n")

# Filter for plasma samples
dat_plasma_demo <- dat_baseline %>%
  filter(sample_type == "plasma")

# Create demographics table for plasma
plasma_demo <- dat_plasma_demo %>%
  group_by(group) %>%
  summarise(
    N = n(),
    Age_mean = mean(age, na.rm = TRUE),
    Age_sd = sd(age, na.rm = TRUE),
    Sex_Female_n = sum(sex == "Female", na.rm = TRUE),
    Sex_Female_pct = mean(sex == "Female", na.rm = TRUE) * 100,
    BMI_mean = mean(bmi, na.rm = TRUE),
    BMI_sd = sd(bmi, na.rm = TRUE),
    HbA1c_mean = mean(hba1c, na.rm = TRUE),
    HbA1c_sd = sd(hba1c, na.rm = TRUE),
    eGFR_mean = mean(eGFR_CKD_epi, na.rm = TRUE),
    eGFR_sd = sd(eGFR_CKD_epi, na.rm = TRUE)
  ) %>%
  mutate(
    Age = sprintf("%.1f ± %.1f", Age_mean, Age_sd),
    Sex_Female = sprintf("%d (%.1f%%)", Sex_Female_n, Sex_Female_pct),
    BMI = sprintf("%.1f ± %.1f", BMI_mean, BMI_sd),
    HbA1c = sprintf("%.1f ± %.1f", HbA1c_mean, HbA1c_sd),
    eGFR = sprintf("%.1f ± %.1f", eGFR_mean, eGFR_sd)
  ) %>%
  select(group, N, Age, Sex_Female, BMI, HbA1c, eGFR)

# Print plasma demographics table
print(plasma_demo %>%
        kable(format = "simple",
              col.names = c("Group", "N", "Age (years)", "Female n (%)", "BMI (kg/m²)", "HbA1c (%)", "eGFR (mL/min/1.73m²)"),
              caption = "Plasma Sample Demographics") %>%
        kable_styling())

# Save plasma demographics table
write.csv(plasma_demo, "C:/Users/netio/Downloads/plasma_demographics_table.csv", row.names = FALSE)

cat("\n\nDemographics tables saved to Downloads folder!\n")
cat("- serum_demographics_table.csv\n")
cat("- plasma_demographics_table.csv\n")


