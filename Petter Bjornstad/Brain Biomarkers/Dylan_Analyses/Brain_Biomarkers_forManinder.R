######Brain biomarkers and clinical characteristics 








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
#library(table1)
library(clusterProfiler)
library('org.Hs.eg.db')
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(readxl)
library(stringr)
library(httr)

qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

names_dat <- names(dat)

dat <- dat %>% filter(study == 'CROCODILE') %>% 
  filter(!is.na(ab40_avg_conc)) %>% filter(visit == 'baseline')

### Demographics 


library(gt)
library(gtsummary)


dat_small <- dat %>% filter(group %in% c('Lean Control', 'Type 1 Diabetes')) %>%
  filter(!is.na(ab40_avg_conc))


desc_table1_fixed <- dat_small %>%
  select(age, sex, group, race_ethnicity, bmi, hba1c, study, diabetes_duration) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      sex ~ 'categorical',
      hba1c ~ "continuous",
      diabetes_duration ~ 'continuous',
      race_ethnicity ~ "categorical",
      study ~ "categorical"
      
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      race_ethnicity ~ "Race/Ethnicity",
      sex ~ 'Sex', 
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      diabetes_duration ~ 'Diabetes Duration, years',
      race_ethnicity ~ 'Race/Ethnicity', 
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave(paste0("/Users/netio/Documents/UofW/Projects/Maninder_Data/CROCODILE_demographics.png"), 
         vwidth = 1200, vheight = 800)

















##### Seven Markers comparisons 

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Define biomarkers
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Create nice labels for the biomarkers
biomarker_labels <- c(
  "ab40_avg_conc" = "Aβ40",
  "ab42_avg_conc" = "Aβ42",
  "tau_avg_conc" = "Tau",
  "nfl_avg_conc" = "NfL",
  "gfap_avg_conc" = "GFAP",
  "ptau_181_avg_conc" = "pTau-181",
  "ptau_217_avg_conc" = "pTau-217"
)

# Set color palette
group_colors <- c("Lean Control" = "#4DBBD5", "Type 1 Diabetes" = "#E64B35")

setwd('/Users/netio/Documents/UofW/Projects/Maninder_Data/')

# Function to test normality and recommend transformation
test_normality <- function(data, biomarker_col, group_col) {
  
  results <- data.frame(
    biomarker = character(),
    group = character(),
    n = integer(),
    shapiro_p_raw = numeric(),
    shapiro_p_log = numeric(),
    skewness_raw = numeric(),
    skewness_log = numeric(),
    transform_recommended = character(),
    stringsAsFactors = FALSE
  )
  
  biomarkers <- unique(data[[biomarker_col]])
  groups <- unique(data[[group_col]])
  
  for (marker in biomarkers) {
    for (grp in groups) {
      
      subset_data <- data %>% 
        filter(.data[[biomarker_col]] == marker, .data[[group_col]] == grp) %>%
        pull(concentration)
      
      # Remove NA values
      subset_data <- subset_data[!is.na(subset_data)]
      
      if (length(subset_data) < 3) next
      
      # Test raw data
      shapiro_raw <- shapiro.test(subset_data)
      skew_raw <- (mean(subset_data) - median(subset_data)) / sd(subset_data)
      
      # Test log-transformed data (add small constant if any zeros)
      min_val <- min(subset_data[subset_data > 0], na.rm = TRUE)
      log_data <- log(subset_data + min_val * 0.01)
      shapiro_log <- shapiro.test(log_data)
      skew_log <- (mean(log_data) - median(log_data)) / sd(log_data)
      
      results <- rbind(results, data.frame(
        biomarker = marker,
        group = grp,
        n = length(subset_data),
        shapiro_p_raw = shapiro_raw$p.value,
        shapiro_p_log = shapiro_log$p.value,
        skewness_raw = abs(skew_raw),
        skewness_log = abs(skew_log),
        transform_recommended = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Determine transformation recommendation per biomarker
  # Criteria: if log transformation improves normality (higher p-value) 
  # and reduces skewness, recommend log
  transform_decision <- results %>%
    group_by(biomarker) %>%
    summarize(
      mean_shapiro_raw = mean(shapiro_p_raw),
      mean_shapiro_log = mean(shapiro_p_log),
      mean_skew_raw = mean(skewness_raw),
      mean_skew_log = mean(skewness_log),
      .groups = 'drop'
    ) %>%
    mutate(
      transform = case_when(
        mean_shapiro_log > mean_shapiro_raw & mean_skew_log < mean_skew_raw ~ "log",
        TRUE ~ "none"
      )
    )
  
  list(
    detailed_results = results,
    transformation_decision = transform_decision
  )
}

# Function to apply transformation
apply_transformation <- function(data, transform_decisions) {
  data_transformed <- data
  
  for (i in 1:nrow(transform_decisions)) {
    marker <- transform_decisions$biomarker[i]
    trans_type <- transform_decisions$transform[i]
    
    if (trans_type == "log") {
      # Get minimum positive value for this biomarker
      min_val <- min(data$concentration[data$biomarker == marker & data$concentration > 0], 
                     na.rm = TRUE)
      
      # Apply log transformation
      data_transformed <- data_transformed %>%
        mutate(concentration = ifelse(
          biomarker == marker,
          log(concentration + min_val * 0.01),
          concentration
        ))
    }
  }
  
  data_transformed
}

# Function to perform Mann-Whitney U test and format p-values
perform_stats_tests <- function(data, biomarker_col, group_col) {
  
  results <- data.frame(
    biomarker = character(),
    test = character(),
    p_value = numeric(),
    p_formatted = character(),
    significance = character(),
    stringsAsFactors = FALSE
  )
  
  biomarkers <- unique(data[[biomarker_col]])
  groups <- unique(data[[group_col]])
  
  if (length(groups) != 2) {
    stop("This function requires exactly 2 groups for comparison")
  }
  
  for (marker in biomarkers) {
    subset_data <- data %>% filter(.data[[biomarker_col]] == marker)
    
    group1_data <- subset_data %>% 
      filter(.data[[group_col]] == groups[1]) %>% 
      pull(concentration)
    
    group2_data <- subset_data %>% 
      filter(.data[[group_col]] == groups[2]) %>% 
      pull(concentration)
    
    # Mann-Whitney U test (Wilcoxon rank-sum test)
    wilcox_test <- wilcox.test(group1_data, group2_data, exact = FALSE)
    
    # Format p-value
    p_val <- wilcox_test$p.value
    p_formatted <- case_when(
      p_val < 0.001 ~ "p < 0.001",
      p_val < 0.01 ~ sprintf("p = %.3f", p_val),
      p_val < 0.05 ~ sprintf("p = %.3f", p_val),
      TRUE ~ sprintf("p = %.3f", p_val)
    )
    
    # Significance stars
    sig <- case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    
    results <- rbind(results, data.frame(
      biomarker = marker,
      test = "Mann-Whitney U",
      p_value = p_val,
      p_formatted = p_formatted,
      significance = sig,
      stringsAsFactors = FALSE
    ))
  }
  
  results
}

# =============================================================================
# MAIN ANALYSIS BEGINS HERE
# =============================================================================

# Prepare data in long format
df_long <- dat_small %>%
  select(group, all_of(qx_var)) %>%
  pivot_longer(cols = all_of(qx_var),
               names_to = "biomarker",
               values_to = "concentration") %>%
  mutate(biomarker_label = biomarker_labels[biomarker])

set.seed(123)

# Test normality and get transformation recommendations
cat("\n==========================================\n")
cat("TESTING NORMALITY AND TRANSFORMATION NEEDS\n")
cat("==========================================\n\n")

normality_results <- test_normality(df_long, "biomarker", "group")

# Print summary
cat("TRANSFORMATION RECOMMENDATIONS:\n")
cat("--------------------------------\n")
for (i in 1:nrow(normality_results$transformation_decision)) {
  marker <- normality_results$transformation_decision$biomarker[i]
  trans <- normality_results$transformation_decision$transform[i]
  marker_label <- biomarker_labels[marker]
  
  cat(sprintf("%-15s (%s): %s\n", 
              marker_label,
              marker,
              ifelse(trans == "log", "LOG TRANSFORM", "NO TRANSFORM")))
}

cat("\nDETAILED RESULTS BY GROUP:\n")
cat("---------------------------\n")
print(normality_results$detailed_results, row.names = FALSE)

# Apply transformations
df_long_transformed <- apply_transformation(df_long, 
                                            normality_results$transformation_decision)

# Add transformation info to labels
df_long_transformed <- df_long_transformed %>%
  left_join(
    normality_results$transformation_decision %>% 
      select(biomarker, transform),
    by = "biomarker"
  ) %>%
  mutate(
    biomarker_label = factor(
      paste0(biomarker_labels[biomarker], 
             ifelse(transform == "log", "\n(log-transformed)", "")),
      levels = paste0(biomarker_labels, 
                      ifelse(normality_results$transformation_decision$transform == "log", 
                             "\n(log-transformed)", ""))
    ),
    y_label = ifelse(transform == "log", 
                     "Log Concentration", 
                     "Concentration (pg/mL)")
  )

# Perform statistical tests
cat("\n\n==========================================\n")
cat("STATISTICAL COMPARISONS (Mann-Whitney U)\n")
cat("==========================================\n\n")

stats_results <- perform_stats_tests(df_long_transformed, "biomarker", "group")

# Print statistical results
for (i in 1:nrow(stats_results)) {
  marker_label <- biomarker_labels[stats_results$biomarker[i]]
  trans_info <- normality_results$transformation_decision %>% 
    filter(biomarker == stats_results$biomarker[i])
  
  trans_label <- ifelse(trans_info$transform == "log", " (log-transformed)", "")
  
  cat(sprintf("%-15s%s: %s %s\n",
              marker_label,
              trans_label,
              stats_results$p_formatted[i],
              stats_results$significance[i]))
}

# =============================================================================
# CREATE PLOTS
# =============================================================================

# Create individual plots for each biomarker
plot_list <- list()

# Ensure we iterate through ALL biomarkers in qx_var
for (marker_name in qx_var) {
  
  # Check if this biomarker exists in the data
  if (!marker_name %in% unique(df_long_transformed$biomarker)) {
    cat(sprintf("WARNING: %s not found in transformed data!\n", marker_name))
    next
  }
  
  marker_label <- biomarker_labels[marker_name]
  trans_info <- normality_results$transformation_decision %>% 
    filter(biomarker == marker_name)
  
  # Get statistical test result
  stat_result <- stats_results %>% filter(biomarker == marker_name)
  
  if (nrow(stat_result) == 0) {
    cat(sprintf("WARNING: No stats result for %s\n", marker_name))
    next
  }
  
  # Create label with transformation info
  plot_title <- ifelse(nrow(trans_info) > 0 && trans_info$transform == "log",
                       paste0(marker_label, " (log-transformed)"),
                       marker_label)
  
  y_axis_label <- ifelse(nrow(trans_info) > 0 && trans_info$transform == "log",
                         "Log Concentration",
                         "Concentration (pg/mL)")
  
  df_subset <- df_long_transformed %>% filter(biomarker == marker_name)
  
  if (nrow(df_subset) == 0) {
    cat(sprintf("WARNING: No data for %s after filtering\n", marker_name))
    next
  }
  
  # Calculate y position for p-value annotation
  y_max <- max(df_subset$concentration, na.rm = TRUE)
  y_min <- min(df_subset$concentration, na.rm = TRUE)
  y_range <- y_max - y_min
  y_pos <- y_max + y_range * 0.05
  
  p <- ggplot(df_subset, aes(x = group, y = concentration, fill = group)) +
    # Violin plot
    geom_violin(alpha = 0.6, trim = TRUE, scale = "width") +
    # Box plot overlaid
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA,
                 position = position_dodge(0.9)) +
    # Individual points (jittered)
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.8, color = "black") +
    # Add significance bracket
    annotate("segment", x = 1, xend = 2, 
             y = y_pos, yend = y_pos,
             color = "black", size = 0.5) +
    annotate("segment", x = 1, xend = 1, 
             y = y_pos, yend = y_pos - y_range * 0.02,
             color = "black", size = 0.5) +
    annotate("segment", x = 2, xend = 2, 
             y = y_pos, yend = y_pos - y_range * 0.02,
             color = "black", size = 0.5) +
    # Add p-value text
    annotate("text", x = 1.5, y = y_pos + y_range * 0.03,
             label = stat_result$p_formatted,
             size = 3.5, fontface = "bold") +
    # Color scheme
    scale_fill_manual(values = group_colors) +
    # Labels
    labs(
      title = plot_title,
      x = NULL,
      y = y_axis_label
    ) +
    # Expand y-axis to fit p-value
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    # Theme
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11, face = "bold"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  plot_list[[marker_label]] <- p
}

cat(sprintf("\nCreated %d individual plots\n", length(plot_list)))

# Combine all plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title = "Brain Biomarkers: Type 1 Diabetes vs Lean Control",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

# Display the combined plot
print(combined_plot)

# Save high-resolution figures
cat("\nSaving high-resolution figures...\n")

# PDF version (vector graphics - best for publications)
ggsave("brain_biomarkers_violin_box.pdf", combined_plot, 
       width = 16, height = 10, dpi = 300, device = "pdf")
cat("Saved: brain_biomarkers_violin_box.pdf\n")

# PNG version (high quality raster)
ggsave("brain_biomarkers_violin_box.png", combined_plot, 
       width = 16, height = 10, dpi = 600, bg = "white")
cat("Saved: brain_biomarkers_violin_box.png (600 DPI)\n")

# =============================================================================
# CREATE FACETED PLOT
# =============================================================================

# Create a dataframe with p-value positions for each biomarker
p_value_data <- df_long_transformed %>%
  group_by(biomarker) %>%
  summarize(
    y_max = max(concentration, na.rm = TRUE),
    y_min = min(concentration, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_pos = y_max + y_range * 0.1,
    x_pos = 1.5
  ) %>%
  left_join(stats_results, by = "biomarker") %>%
  left_join(
    normality_results$transformation_decision %>% select(biomarker, transform),
    by = "biomarker"
  ) %>%
  mutate(
    biomarker_label = factor(
      paste0(biomarker_labels[biomarker], 
             ifelse(transform == "log", "\n(log-transformed)", "")),
      levels = paste0(biomarker_labels[qx_var], 
                      ifelse(normality_results$transformation_decision$transform[
                        match(qx_var, normality_results$transformation_decision$biomarker)] == "log", 
                        "\n(log-transformed)", ""))
    )
  )

# Ensure biomarker_label factor levels match in df_long_transformed
df_long_transformed <- df_long_transformed %>%
  mutate(
    biomarker_label = factor(
      paste0(biomarker_labels[biomarker], 
             ifelse(transform == "log", "\n(log-transformed)", "")),
      levels = paste0(biomarker_labels[qx_var], 
                      ifelse(normality_results$transformation_decision$transform[
                        match(qx_var, normality_results$transformation_decision$biomarker)] == "log", 
                        "\n(log-transformed)", ""))
    )
  )

faceted_plot <- ggplot(df_long_transformed, aes(x = group, y = concentration, fill = group)) +
  geom_violin(alpha = 0.6, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.5, color = "black") +
  # Add p-value text
  geom_text(data = p_value_data, 
            aes(x = x_pos, y = y_pos, label = p_formatted),
            inherit.aes = FALSE,
            size = 3, fontface = "bold") +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~ biomarker_label, scales = "free_y", ncol = 4) +
  labs(
    title = "Brain Biomarkers: Type 1 Diabetes vs Lean Control",
    x = NULL,
    y = "Concentration",
    fill = "Group"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey90", size = 0.3)
  )

print(faceted_plot)

# Save faceted version
cat("\nSaving faceted plot...\n")
ggsave("brain_biomarkers_faceted.pdf", faceted_plot, 
       width = 14, height = 8, dpi = 300, device = "pdf")
cat("Saved: brain_biomarkers_faceted.pdf\n")

ggsave("brain_biomarkers_faceted.png", faceted_plot, 
       width = 14, height = 8, dpi = 600, bg = "white")
cat("Saved: brain_biomarkers_faceted.png (600 DPI)\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n\n==========================================\n")
cat("SUMMARY\n")
cat("==========================================\n")
cat(sprintf("Total biomarkers analyzed: %d\n", length(qx_var)))
cat(sprintf("Log-transformed: %d\n", 
            sum(normality_results$transformation_decision$transform == "log")))
cat(sprintf("No transformation: %d\n", 
            sum(normality_results$transformation_decision$transform == "none")))

cat("\nStatistical significance:\n")
cat(sprintf("  p < 0.001 (***): %d\n", sum(stats_results$p_value < 0.001)))
cat(sprintf("  p < 0.01  (**):  %d\n", sum(stats_results$p_value >= 0.001 & stats_results$p_value < 0.01)))
cat(sprintf("  p < 0.05  (*):   %d\n", sum(stats_results$p_value >= 0.01 & stats_results$p_value < 0.05)))
cat(sprintf("  p >= 0.05 (ns):  %d\n", sum(stats_results$p_value >= 0.05)))

cat("\nPlots have been generated with appropriate transformations and statistical comparisons.\n")
cat("Statistical test: Mann-Whitney U (Wilcoxon rank-sum test)\n")




## t-tests
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Define biomarkers
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Create nice labels for the biomarkers
biomarker_labels <- c(
  "ab40_avg_conc" = "Aβ40",
  "ab42_avg_conc" = "Aβ42",
  "tau_avg_conc" = "Tau",
  "nfl_avg_conc" = "NfL",
  "gfap_avg_conc" = "GFAP",
  "ptau_181_avg_conc" = "pTau-181",
  "ptau_217_avg_conc" = "pTau-217"
)

# Set color palette
group_colors <- c("Lean Control" = "#4DBBD5", "Type 1 Diabetes" = "#E64B35")

# Function to test normality and recommend transformation
test_normality <- function(data, biomarker_col, group_col) {
  
  results <- data.frame(
    biomarker = character(),
    group = character(),
    n = integer(),
    shapiro_p_raw = numeric(),
    shapiro_p_log = numeric(),
    skewness_raw = numeric(),
    skewness_log = numeric(),
    transform_recommended = character(),
    stringsAsFactors = FALSE
  )
  
  biomarkers <- unique(data[[biomarker_col]])
  groups <- unique(data[[group_col]])
  
  for (marker in biomarkers) {
    for (grp in groups) {
      
      subset_data <- data %>% 
        filter(.data[[biomarker_col]] == marker, .data[[group_col]] == grp) %>%
        pull(concentration)
      
      # Remove NA values
      subset_data <- subset_data[!is.na(subset_data)]
      
      if (length(subset_data) < 3) next
      
      # Test raw data
      shapiro_raw <- shapiro.test(subset_data)
      skew_raw <- (mean(subset_data) - median(subset_data)) / sd(subset_data)
      
      # Test log-transformed data (add small constant if any zeros)
      min_val <- min(subset_data[subset_data > 0], na.rm = TRUE)
      log_data <- log(subset_data + min_val * 0.01)
      shapiro_log <- shapiro.test(log_data)
      skew_log <- (mean(log_data) - median(log_data)) / sd(log_data)
      
      results <- rbind(results, data.frame(
        biomarker = marker,
        group = grp,
        n = length(subset_data),
        shapiro_p_raw = shapiro_raw$p.value,
        shapiro_p_log = shapiro_log$p.value,
        skewness_raw = abs(skew_raw),
        skewness_log = abs(skew_log),
        transform_recommended = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Determine transformation recommendation per biomarker
  # Criteria: if log transformation improves normality (higher p-value) 
  # and reduces skewness, recommend log
  transform_decision <- results %>%
    group_by(biomarker) %>%
    summarize(
      mean_shapiro_raw = mean(shapiro_p_raw),
      mean_shapiro_log = mean(shapiro_p_log),
      mean_skew_raw = mean(skewness_raw),
      mean_skew_log = mean(skewness_log),
      .groups = 'drop'
    ) %>%
    mutate(
      transform = case_when(
        mean_shapiro_log > mean_shapiro_raw & mean_skew_log < mean_skew_raw ~ "log",
        TRUE ~ "none"
      )
    )
  
  list(
    detailed_results = results,
    transformation_decision = transform_decision
  )
}

# Function to apply transformation
apply_transformation <- function(data, transform_decisions) {
  data_transformed <- data
  
  for (i in 1:nrow(transform_decisions)) {
    marker <- transform_decisions$biomarker[i]
    trans_type <- transform_decisions$transform[i]
    
    if (trans_type == "log") {
      # Get minimum positive value for this biomarker
      min_val <- min(data$concentration[data$biomarker == marker & data$concentration > 0], 
                     na.rm = TRUE)
      
      # Apply log transformation
      data_transformed <- data_transformed %>%
        mutate(concentration = ifelse(
          biomarker == marker,
          log(concentration + min_val * 0.01),
          concentration
        ))
    }
  }
  
  data_transformed
}

# Function to perform t-tests and format p-values
perform_stats_tests <- function(data, biomarker_col, group_col) {
  
  results <- data.frame(
    biomarker = character(),
    test = character(),
    p_value = numeric(),
    p_formatted = character(),
    significance = character(),
    mean_diff = numeric(),
    stringsAsFactors = FALSE
  )
  
  biomarkers <- unique(data[[biomarker_col]])
  groups <- unique(data[[group_col]])
  
  if (length(groups) != 2) {
    stop("This function requires exactly 2 groups for comparison")
  }
  
  for (marker in biomarkers) {
    subset_data <- data %>% filter(.data[[biomarker_col]] == marker)
    
    group1_data <- subset_data %>% 
      filter(.data[[group_col]] == groups[1]) %>% 
      pull(concentration)
    
    group2_data <- subset_data %>% 
      filter(.data[[group_col]] == groups[2]) %>% 
      pull(concentration)
    
    # Two-sample t-test (Welch's t-test by default, which doesn't assume equal variances)
    t_test <- t.test(group1_data, group2_data)
    
    # Calculate mean difference
    mean_diff <- mean(group2_data, na.rm = TRUE) - mean(group1_data, na.rm = TRUE)
    
    # Format p-value
    p_val <- t_test$p.value
    p_formatted <- case_when(
      p_val < 0.001 ~ "p < 0.001",
      p_val < 0.01 ~ sprintf("p = %.3f", p_val),
      p_val < 0.05 ~ sprintf("p = %.3f", p_val),
      TRUE ~ sprintf("p = %.3f", p_val)
    )
    
    # Significance stars
    sig <- case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    
    results <- rbind(results, data.frame(
      biomarker = marker,
      test = "Welch's t-test",
      p_value = p_val,
      p_formatted = p_formatted,
      significance = sig,
      mean_diff = mean_diff,
      stringsAsFactors = FALSE
    ))
  }
  
  results
}

# =============================================================================
# MAIN ANALYSIS BEGINS HERE
# =============================================================================

# Prepare data in long format
# NOTE: Replace 'dat_small' with your actual dataset name!
df_long <- dat_small %>%
  select(group, all_of(qx_var)) %>%
  pivot_longer(cols = all_of(qx_var),
               names_to = "biomarker",
               values_to = "concentration") %>%
  mutate(biomarker_label = biomarker_labels[biomarker])

set.seed(123)

# Test normality and get transformation recommendations
cat("\n==========================================\n")
cat("TESTING NORMALITY AND TRANSFORMATION NEEDS\n")
cat("==========================================\n\n")

normality_results <- test_normality(df_long, "biomarker", "group")

# Print summary
cat("TRANSFORMATION RECOMMENDATIONS:\n")
cat("--------------------------------\n")
for (i in 1:nrow(normality_results$transformation_decision)) {
  marker <- normality_results$transformation_decision$biomarker[i]
  trans <- normality_results$transformation_decision$transform[i]
  marker_label <- biomarker_labels[marker]
  
  cat(sprintf("%-15s (%s): %s\n", 
              marker_label,
              marker,
              ifelse(trans == "log", "LOG TRANSFORM", "NO TRANSFORM")))
}

cat("\nDETAILED RESULTS BY GROUP:\n")
cat("---------------------------\n")
print(normality_results$detailed_results, row.names = FALSE)

# Apply transformations
df_long_transformed <- apply_transformation(df_long, 
                                            normality_results$transformation_decision)

# Add transformation info to labels
df_long_transformed <- df_long_transformed %>%
  left_join(
    normality_results$transformation_decision %>% 
      select(biomarker, transform),
    by = "biomarker"
  ) %>%
  mutate(
    biomarker_label = factor(
      paste0(biomarker_labels[biomarker], 
             ifelse(transform == "log", "\n(log-transformed)", "")),
      levels = paste0(biomarker_labels, 
                      ifelse(normality_results$transformation_decision$transform == "log", 
                             "\n(log-transformed)", ""))
    ),
    y_label = ifelse(transform == "log", 
                     "Log Concentration", 
                     "Concentration (pg/mL)")
  )

# Perform statistical tests
cat("\n\n==========================================\n")
cat("STATISTICAL COMPARISONS (Welch's t-test)\n")
cat("==========================================\n\n")

stats_results <- perform_stats_tests(df_long_transformed, "biomarker", "group")

# Print statistical results
for (i in 1:nrow(stats_results)) {
  marker_label <- biomarker_labels[stats_results$biomarker[i]]
  trans_info <- normality_results$transformation_decision %>% 
    filter(biomarker == stats_results$biomarker[i])
  
  trans_label <- ifelse(trans_info$transform == "log", " (log-transformed)", "")
  
  cat(sprintf("%-15s%s: %s %s (mean diff = %.3f)\n",
              marker_label,
              trans_label,
              stats_results$p_formatted[i],
              stats_results$significance[i],
              stats_results$mean_diff[i]))
}

# =============================================================================
# CREATE PLOTS
# =============================================================================

# Create individual plots for each biomarker
plot_list <- list()

# Ensure we iterate through ALL biomarkers in qx_var
for (marker_name in qx_var) {
  
  # Check if this biomarker exists in the data
  if (!marker_name %in% unique(df_long_transformed$biomarker)) {
    cat(sprintf("WARNING: %s not found in transformed data!\n", marker_name))
    next
  }
  
  marker_label <- biomarker_labels[marker_name]
  trans_info <- normality_results$transformation_decision %>% 
    filter(biomarker == marker_name)
  
  # Get statistical test result
  stat_result <- stats_results %>% filter(biomarker == marker_name)
  
  if (nrow(stat_result) == 0) {
    cat(sprintf("WARNING: No stats result for %s\n", marker_name))
    next
  }
  
  # Create label with transformation info
  plot_title <- ifelse(nrow(trans_info) > 0 && trans_info$transform == "log",
                       paste0(marker_label, " (log-transformed)"),
                       marker_label)
  
  y_axis_label <- ifelse(nrow(trans_info) > 0 && trans_info$transform == "log",
                         "Log Concentration",
                         "Concentration (pg/mL)")
  
  df_subset <- df_long_transformed %>% filter(biomarker == marker_name)
  
  if (nrow(df_subset) == 0) {
    cat(sprintf("WARNING: No data for %s after filtering\n", marker_name))
    next
  }
  
  # Calculate y position for p-value annotation
  y_max <- max(df_subset$concentration, na.rm = TRUE)
  y_min <- min(df_subset$concentration, na.rm = TRUE)
  y_range <- y_max - y_min
  y_pos <- y_max + y_range * 0.05
  
  p <- ggplot(df_subset, aes(x = group, y = concentration, fill = group)) +
    # Violin plot
    geom_violin(alpha = 0.6, trim = TRUE, scale = "width") +
    # Box plot overlaid
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA,
                 position = position_dodge(0.9)) +
    # Individual points (jittered)
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.8, color = "black") +
    # Add significance bracket
    annotate("segment", x = 1, xend = 2, 
             y = y_pos, yend = y_pos,
             color = "black", size = 0.5) +
    annotate("segment", x = 1, xend = 1, 
             y = y_pos, yend = y_pos - y_range * 0.02,
             color = "black", size = 0.5) +
    annotate("segment", x = 2, xend = 2, 
             y = y_pos, yend = y_pos - y_range * 0.02,
             color = "black", size = 0.5) +
    # Add p-value text
    annotate("text", x = 1.5, y = y_pos + y_range * 0.03,
             label = stat_result$p_formatted,
             size = 3.5, fontface = "bold") +
    # Color scheme
    scale_fill_manual(values = group_colors) +
    # Labels
    labs(
      title = plot_title,
      x = NULL,
      y = y_axis_label
    ) +
    # Expand y-axis to fit p-value
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    # Theme
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11, face = "bold"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  plot_list[[marker_label]] <- p
}

cat(sprintf("\nCreated %d individual plots\n", length(plot_list)))

# Combine all plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title = "Brain Biomarkers: Type 1 Diabetes vs Lean Control",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

# Display the combined plot
print(combined_plot)

# Save high-resolution figures
cat("\nSaving high-resolution figures...\n")

# PDF version (vector graphics - best for publications)
ggsave("brain_biomarkers_violin_box_ttest.pdf", combined_plot, 
       width = 16, height = 10, dpi = 300, device = "pdf")
cat("Saved: brain_biomarkers_violin_box_ttest.pdf\n")

# PNG version (high quality raster)
ggsave("brain_biomarkers_violin_box_ttest.png", combined_plot, 
       width = 16, height = 10, dpi = 600, bg = "white")
cat("Saved: brain_biomarkers_violin_box_ttest.png (600 DPI)\n")

# =============================================================================
# CREATE FACETED PLOT
# =============================================================================

# Create a dataframe with p-value positions for each biomarker
p_value_data <- df_long_transformed %>%
  group_by(biomarker) %>%
  summarize(
    y_max = max(concentration, na.rm = TRUE),
    y_min = min(concentration, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_pos = y_max + y_range * 0.1,
    x_pos = 1.5
  ) %>%
  left_join(stats_results, by = "biomarker") %>%
  left_join(
    normality_results$transformation_decision %>% select(biomarker, transform),
    by = "biomarker"
  ) %>%
  mutate(
    biomarker_label = factor(
      paste0(biomarker_labels[biomarker], 
             ifelse(transform == "log", "\n(log-transformed)", "")),
      levels = paste0(biomarker_labels[qx_var], 
                      ifelse(normality_results$transformation_decision$transform[
                        match(qx_var, normality_results$transformation_decision$biomarker)] == "log", 
                        "\n(log-transformed)", ""))
    )
  )

# Ensure biomarker_label factor levels match in df_long_transformed
df_long_transformed <- df_long_transformed %>%
  mutate(
    biomarker_label = factor(
      paste0(biomarker_labels[biomarker], 
             ifelse(transform == "log", "\n(log-transformed)", "")),
      levels = paste0(biomarker_labels[qx_var], 
                      ifelse(normality_results$transformation_decision$transform[
                        match(qx_var, normality_results$transformation_decision$biomarker)] == "log", 
                        "\n(log-transformed)", ""))
    )
  )

faceted_plot <- ggplot(df_long_transformed, aes(x = group, y = concentration, fill = group)) +
  geom_violin(alpha = 0.6, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.5, color = "black") +
  # Add p-value text
  geom_text(data = p_value_data, 
            aes(x = x_pos, y = y_pos, label = p_formatted),
            inherit.aes = FALSE,
            size = 3, fontface = "bold") +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~ biomarker_label, scales = "free_y", ncol = 4) +
  labs(
    title = "Brain Biomarkers: Type 1 Diabetes vs Lean Control",
    x = NULL,
    y = "Concentration",
    fill = "Group"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey90", size = 0.3)
  )

print(faceted_plot)

# Save faceted version
cat("\nSaving faceted plot...\n")
ggsave("brain_biomarkers_faceted_ttest.pdf", faceted_plot, 
       width = 14, height = 8, dpi = 300, device = "pdf")
cat("Saved: brain_biomarkers_faceted_ttest.pdf\n")

ggsave("brain_biomarkers_faceted_ttest.png", faceted_plot, 
       width = 14, height = 8, dpi = 600, bg = "white")
cat("Saved: brain_biomarkers_faceted_ttest.png (600 DPI)\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n\n==========================================\n")
cat("SUMMARY\n")
cat("==========================================\n")
cat(sprintf("Total biomarkers analyzed: %d\n", length(qx_var)))
cat(sprintf("Log-transformed: %d\n", 
            sum(normality_results$transformation_decision$transform == "log")))
cat(sprintf("No transformation: %d\n", 
            sum(normality_results$transformation_decision$transform == "none")))

cat("\nStatistical significance:\n")
cat(sprintf("  p < 0.001 (***): %d\n", sum(stats_results$p_value < 0.001)))
cat(sprintf("  p < 0.01  (**):  %d\n", sum(stats_results$p_value >= 0.001 & stats_results$p_value < 0.01)))
cat(sprintf("  p < 0.05  (*):   %d\n", sum(stats_results$p_value >= 0.01 & stats_results$p_value < 0.05)))
cat(sprintf("  p >= 0.05 (ns):  %d\n", sum(stats_results$p_value >= 0.05)))

cat("\nPlots have been generated with appropriate transformations and statistical comparisons.\n")
cat("Statistical test: Welch's t-test (two-sample t-test)\n")























#Clinical characteristics (eGFR, UACR, HBA1C, Clamp, PET, DEXA, glycemia, insulin sensitivity)


data_dictionary <- readxl::read_xlsx('C:/Users/netio/Downloads/data_dictionary_master.xlsx')
#data_dictionary <- data_dictionary %>% filter(form_name %in% c('clamp', 'UACR', 'fmri', 'eGFR', 'FSOC', 'brain_mri', 'pet_scan', 'dextran_data'))


traits_of_interest <- c('acr_u', 
                        'eGFR_bedside_Schwartz', 'eGFR_CKD_epi', 'eGFR_fas_cr', 'eGFR_fas_cr_cysc','eGFR_Zap','eGFR_Schwartz', 
                        'hba1c')

dat_analysis <- dat %>% 
  dplyr::select(record_id, group, hba1c, #age, sex, bmi, hba1c, study, 
                all_of(qx_var), 
                any_of(data_dictionary$variable_name))


dat_analysis <- dat %>% 
  dplyr::select(record_id, group, all_of(qx_var), acr_u, adipose_ir, cholesterol, hba1c, eGFR_CKD_epi, eGFR_fas_cr, fbg, ldl, hdl, left_kidney_volume_ml, right_kidney_volume_ml, 
                triglycerides, urine_glucose_bl, avg_k_fsoc, avg_c_fsoc, avg_m_fsoc, homa_ir, search_eis)


dat_clean <- dat_analysis[complete.cases(dat_analysis), ]

library(corrplot)

dat_filtered <- dat %>% filter(record_id %in% dat_clean$record_id)

# Identify columns with NO missing values
complete_cols <- colnames(dat_filtered)[colSums(is.na(dat_filtered)) == 0]

# View the list
print(complete_cols)

# Extract only those columns
data_complete <- dat_filtered[, complete_cols]

# Identify columns with more than 1 unique value (excluding NAs)
varying_cols <- sapply(data_complete, function(x) length(unique(na.omit(x))) > 1)


# Get list of metabolomics variables to remove
metabolomics_vars <- data_dictionary %>%
  filter(form_name %in% c("metabolomics", "metabolomics_blood_raw")) %>%
  pull(variable_name)



# Keep only varying columns
data_varying <- data_complete[, varying_cols]

data_varying <- data_varying %>% 
  dplyr::select(
    -c('croc_id', 'date', 'mrn'),
    -matches("^glucose_\\d+"),
    -matches('^insulin_\\d+'),
    -matches('^ffa_\\d+'),
    -matches("^c0"),             
    -matches("^ac\\d+"), 
    -matches('^c1'), 
    -matches('^hmdb\\d+'), 
    -any_of(metabolomics_vars)
  ) %>% 
  dplyr::select(-hx_met_positive___1, -insulin_injections_timepoint, -insulin_med_timepoint, -insulin_pump_timepoint, 
                -u24_labs, -insulin_minus_10, -insulin_minus_20, -glucose_minus_10, -glucose_minus_20, -ffa_minus_10, -ffa_minus_20, -group_risk, 
                -record_id, -sex, -race, -ethnicity)


data_dictionary_small <- data_dictionary %>% filter(variable_name %in% names(data_varying))


#demographics (only all data) 

desc_table1_fixed <- data_complete %>%
  select(age, sex, group, race_ethnicity, bmi, hba1c, study) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      sex ~ 'categorical',
      hba1c ~ "continuous",
 #     diabetes_duration ~ 'continuous',
      race_ethnicity ~ "categorical",
      study ~ "categorical"
      
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      race_ethnicity ~ "Race/Ethnicity",
      sex ~ 'Sex', 
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
 #     diabetes_duration ~ 'Diabetes Duration, years',
      race_ethnicity ~ 'Race/Ethnicity', 
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave(paste0("/Users/netio/Documents/UofW/Projects/Maninder_Data/CROCODILE_demographics_noNA.png"), 
         vwidth = 1200, vheight = 800)























#### Elastic net models 
library(readxl)
library(tidyr)
library(stringr)
library(dplyr)
library(glmnet)
library(ggplot2)
library(patchwork)

set.seed(123)

# Define biomarkers (outcomes)
biomarkers <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
                "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", 
                "ptau_217_avg_conc")

biomarker_labels <- c(
  "ab40_avg_conc" = "Aβ40",
  "ab42_avg_conc" = "Aβ42",
  "tau_avg_conc" = "Tau",
  "nfl_avg_conc" = "NfL",
  "gfap_avg_conc" = "GFAP",
  "ptau_181_avg_conc" = "pTau-181",
  "ptau_217_avg_conc" = "pTau-217"
)


#### Functions 

test_transformation <- function(x, var_name = "") {
  # Remove NAs
  x <- x[!is.na(x)]
  
  if (length(x) < 3) {
    return(list(transform = "none", reason = "insufficient_data"))
  }
  
  # Test raw data
  shapiro_raw <- tryCatch(shapiro.test(x)$p.value, error = function(e) NA)
  skew_raw <- (mean(x) - median(x)) / sd(x)
  
  # Test log-transformed (add small constant if zeros)
  if (min(x) <= 0) {
    x_log <- log(x - min(x) + 0.01)
  } else {
    x_log <- log(x)
  }
  
  shapiro_log <- tryCatch(shapiro.test(x_log)$p.value, error = function(e) NA)
  skew_log <- (mean(x_log) - median(x_log)) / sd(x_log)
  
  # Decision: use log if it improves normality and reduces skewness
  if (!is.na(shapiro_log) && !is.na(shapiro_raw)) {
    if (shapiro_log > shapiro_raw && abs(skew_log) < abs(skew_raw)) {
      return(list(
        transform = "log",
        shapiro_raw = shapiro_raw,
        shapiro_log = shapiro_log,
        skew_raw = abs(skew_raw),
        skew_log = abs(skew_log)
      ))
    }
  }
  
  return(list(
    transform = "none",
    shapiro_raw = shapiro_raw,
    shapiro_log = ifelse(is.na(shapiro_log), NA, shapiro_log),
    skew_raw = abs(skew_raw),
    skew_log = abs(skew_log)
  ))
}


apply_transformations <- function(data, transform_decisions) {
  data_transformed <- data
  
  for (var in names(transform_decisions)) {
    if (transform_decisions[[var]]$transform == "log") {
      x <- data_transformed[[var]]
      if (min(x, na.rm = TRUE) <= 0) {
        data_transformed[[var]] <- log(x - min(x, na.rm = TRUE) + 0.01)
      } else {
        data_transformed[[var]] <- log(x)
      }
    }
  }
  
  data_transformed
}


run_elastic_net <- function(data, outcome_var, predictor_vars, cohort_name) {
  
  cat(sprintf("\n=== Running Elastic Net: %s in %s ===\n", 
              outcome_var, cohort_name))
  
  # Prepare data
  complete_data <- data %>%
    select(all_of(c(outcome_var, predictor_vars))) %>%
    drop_na()
  
  if (nrow(complete_data) < 10) {
    cat("Insufficient data (n < 10). Skipping.\n")
    return(NULL)
  }
  
  cat(sprintf("Sample size: %d\n", nrow(complete_data)))
  
  # Prepare matrices
  X <- as.matrix(complete_data[, predictor_vars])
  y <- complete_data[[outcome_var]]
  
  # Standardize predictors (glmnet does this internally, but we track it)
  X_scaled <- scale(X)
  
  # Set up cross-validation folds
  nFolds <- min(5, nrow(complete_data))
  foldid <- sample(rep(seq(nFolds), length.out = nrow(complete_data)))
  
  # Run elastic net with alpha = 0.5 (equal mix of L1 and L2)
  cv_fit <- cv.glmnet(
    x = X_scaled,
    y = y,
    family = "gaussian",
    alpha = 0.5,
    nfolds = nFolds,
    foldid = foldid,
    keep = TRUE
  )
  
  # Extract coefficients at lambda.min
  coef_min <- coef(cv_fit, s = "lambda.min")
  coef_df <- data.frame(
    Variable = rownames(coef_min),
    Coefficient = as.numeric(coef_min),
    stringsAsFactors = FALSE
  ) %>%
    filter(Variable != "(Intercept)", Coefficient != 0) %>%
    arrange(desc(abs(Coefficient)))
  
  # Extract coefficients at lambda.1se (more regularized)
  coef_1se <- coef(cv_fit, s = "lambda.1se")
  coef_1se_df <- data.frame(
    Variable = rownames(coef_1se),
    Coefficient = as.numeric(coef_1se),
    stringsAsFactors = FALSE
  ) %>%
    filter(Variable != "(Intercept)", Coefficient != 0) %>%
    arrange(desc(abs(Coefficient)))
  
  # Calculate R-squared
  predictions_min <- predict(cv_fit, newx = X_scaled, s = "lambda.min")
  r2_min <- 1 - sum((y - predictions_min)^2) / sum((y - mean(y))^2)
  
  predictions_1se <- predict(cv_fit, newx = X_scaled, s = "lambda.1se")
  r2_1se <- 1 - sum((y - predictions_1se)^2) / sum((y - mean(y))^2)
  
  cat(sprintf("Lambda.min: %.6f (R² = %.4f, %d variables)\n", 
              cv_fit$lambda.min, r2_min, nrow(coef_df)))
  cat(sprintf("Lambda.1se: %.6f (R² = %.4f, %d variables)\n", 
              cv_fit$lambda.1se, r2_1se, nrow(coef_1se_df)))
  
  list(
    cv_fit = cv_fit,
    coef_min = coef_df,
    coef_1se = coef_1se_df,
    r2_min = r2_min,
    r2_1se = r2_1se,
    n = nrow(complete_data),
    outcome = outcome_var,
    cohort = cohort_name
  )
}


plot_cv_results <- function(result, biomarker_label) {
  
  cv_fit <- result$cv_fit
  
  # Lambda plot
  lambda_df <- data.frame(
    lambda = cv_fit$lambda,
    mse = cv_fit$cvm,
    mse_lower = cv_fit$cvlo,
    mse_upper = cv_fit$cvup,
    nzero = cv_fit$nzero
  )
  
  p1 <- ggplot(lambda_df, aes(x = log(lambda), y = mse)) +
    geom_point(color = "#E64B35", size = 2) +
    geom_errorbar(aes(ymin = mse_lower, ymax = mse_upper), 
                  width = 0.1, color = "#E64B35", alpha = 0.5) +
    geom_vline(xintercept = log(cv_fit$lambda.min), 
               linetype = "dashed", color = "blue", size = 1) +
    geom_vline(xintercept = log(cv_fit$lambda.1se), 
               linetype = "dashed", color = "darkgreen", size = 1) +
    labs(
      title = sprintf("%s - %s", biomarker_label, result$cohort),
      x = "Log(Lambda)",
      y = "Mean Squared Error",
      subtitle = sprintf("λ.min (blue): %.4f | λ.1se (green): %.4f", 
                         cv_fit$lambda.min, cv_fit$lambda.1se)
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  
  # Coefficient plot at lambda.min
  if (nrow(result$coef_min) > 0) {
    coef_plot_df <- result$coef_min %>%
      arrange(Coefficient) %>%
      mutate(Variable = factor(Variable, levels = Variable))
    
    p2 <- ggplot(coef_plot_df, aes(x = Coefficient, y = Variable)) +
      geom_segment(aes(x = 0, xend = Coefficient, y = Variable, yend = Variable),
                   color = "grey50", size = 1) +
      geom_point(color = "#4DBBD5", size = 4) +
      geom_vline(xintercept = 0, linetype = "solid", color = "black") +
      labs(
        title = "Selected Features (λ.min)",
        x = "Standardized Coefficient",
        y = NULL,
        subtitle = sprintf("R² = %.3f, n = %d", result$r2_min, result$n)
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9)
      )
  } else {
    p2 <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = "No variables selected", size = 6) +
      theme_void()
  }
  
  p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))
}


dat <- dat_clean
# Define predictors (exclude biomarkers and group)
all_vars <- setdiff(names(dat), c(biomarkers, "group"))
predictor_vars <- all_vars

cat("Biomarkers (outcomes):", paste(biomarker_labels[biomarkers], collapse = ", "), "\n")
cat("Number of predictors:", length(predictor_vars), "\n")
cat("Predictors:", paste(predictor_vars, collapse = ", "), "\n\n")


transform_decisions <- list()

# Test biomarkers
cat("BIOMARKERS:\n")
for (var in biomarkers) {
  result <- test_transformation(dat[[var]], var)
  transform_decisions[[var]] <- result
  cat(sprintf("%-20s: %s", biomarker_labels[var], 
              ifelse(result$transform == "log", "LOG TRANSFORM", "NO TRANSFORM")))
  if (!is.null(result$shapiro_raw)) {
    cat(sprintf(" (Shapiro raw: %.3f, log: %.3f)\n", 
                result$shapiro_raw, result$shapiro_log))
  } else {
    cat("\n")
  }
}

# Test predictors
cat("\nPREDICTORS:\n")
for (var in predictor_vars) {
  result <- test_transformation(dat[[var]], var)
  transform_decisions[[var]] <- result
  cat(sprintf("%-25s: %s\n", var, 
              ifelse(result$transform == "log", "LOG TRANSFORM", "NO TRANSFORM")))
}

# Apply transformations
dat_transformed <- apply_transformations(dat, transform_decisions)

cat("\nTransformations applied.\n")
cat(sprintf("Log-transformed variables: %d\n", 
            sum(sapply(transform_decisions, function(x) x$transform == "log"))))



#Run elastic net models

# Initialize results lists HERE - BEFORE the loop
results_all <- list()
plot_list <- list()

# Define cohorts
cohorts <- list(
  "All" = dat_transformed,
  "T1D_only" = dat_transformed %>% filter(group == "Type 1 Diabetes")
)

counter <- 1

for (cohort_name in names(cohorts)) {
  cohort_data <- cohorts[[cohort_name]]
  
  cat(sprintf("\n--- COHORT: %s (n = %d) ---\n", cohort_name, nrow(cohort_data)))
  
  for (biomarker in biomarkers) {
    
    result <- tryCatch({
      run_elastic_net(
        data = cohort_data,
        outcome_var = biomarker,
        predictor_vars = predictor_vars,
        cohort_name = cohort_name
      )
    }, error = function(e) {
      cat(sprintf("ERROR in %s - %s: %s\n", biomarker, cohort_name, e$message))
      return(NULL)
    })
    
    if (!is.null(result)) {
      results_all[[counter]] <- result
      
      # Create plot
      tryCatch({
        p <- plot_cv_results(result, biomarker_labels[biomarker])
        plot_list[[counter]] <- p
      }, error = function(e) {
        cat(sprintf("ERROR creating plot for %s - %s: %s\n", 
                    biomarker, cohort_name, e$message))
      })
      
      counter <- counter + 1
    }
  }
}

cat(sprintf("\n\nTotal models successfully run: %d\n", length(results_all)))

# Create output directory
dir.create("elastic_net_results", showWarnings = FALSE)

# Save individual plots
for (i in seq_along(plot_list)) {
  result <- results_all[[i]]
  filename <- sprintf("elastic_net_results/%s_%s_plot.png",
                      result$outcome, 
                      gsub(" ", "_", result$cohort))
  
  ggsave(filename, plot_list[[i]], width = 12, height = 6, dpi = 300)
  cat(sprintf("Saved: %s\n", filename))
}

# Save coefficients
for (i in seq_along(results_all)) {
  result <- results_all[[i]]
  
  # Lambda.min coefficients
  if (nrow(result$coef_min) > 0) {
    filename <- sprintf("elastic_net_results/%s_%s_coefficients_lambda_min.csv",
                        result$outcome, 
                        gsub(" ", "_", result$cohort))
    
    write.csv(result$coef_min, filename, row.names = FALSE)
    cat(sprintf("Saved: %s\n", filename))
  }
  
  # Lambda.1se coefficients
  if (nrow(result$coef_1se) > 0) {
    filename <- sprintf("elastic_net_results/%s_%s_coefficients_lambda_1se.csv",
                        result$outcome, 
                        gsub(" ", "_", result$cohort))
    
    write.csv(result$coef_1se, filename, row.names = FALSE)
  }
}

# Create summary table
summary_df <- data.frame(
  Biomarker = sapply(results_all, function(x) biomarker_labels[x$outcome]),
  Cohort = sapply(results_all, function(x) x$cohort),
  N = sapply(results_all, function(x) x$n),
  Lambda_min = sapply(results_all, function(x) x$cv_fit$lambda.min),
  R2_min = sapply(results_all, function(x) x$r2_min),
  N_vars_min = sapply(results_all, function(x) nrow(x$coef_min)),
  Lambda_1se = sapply(results_all, function(x) x$cv_fit$lambda.1se),
  R2_1se = sapply(results_all, function(x) x$r2_1se),
  N_vars_1se = sapply(results_all, function(x) nrow(x$coef_1se))
)

write.csv(summary_df, "elastic_net_results/summary_all_models.csv", row.names = FALSE)
cat("Saved: elastic_net_results/summary_all_models.csv\n")

cat("\n==========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==========================================\n")
cat(sprintf("Total models run: %d\n", length(results_all)))
print(summary_df)











#Follow-up Analysis

# ==========================================
# COMPARISON: ELASTIC NET VS LINEAR REGRESSION
# ==========================================

cat("\n==========================================\n")
cat("ELASTIC NET VS LINEAR REGRESSION COMPARISON\n")
cat("==========================================\n\n")

comparison_table <- data.frame(
  Biomarker = character(),
  Cohort = character(),
  N = integer(),
  ElasticNet_R2 = numeric(),
  LinearReg_R2 = numeric(),
  LinearReg_AdjR2 = numeric(),
  N_predictors = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(results_all)) {
  result <- results_all[[i]]
  
  # Find matching regression result
  matching_idx <- which(sapply(regression_results_all, function(x) 
    x$outcome == result$outcome && x$cohort == result$cohort))
  
  if (length(matching_idx) > 0) {
    reg_result <- regression_results_all[[matching_idx[1]]]
    
    comparison_table <- rbind(comparison_table, data.frame(
      Biomarker = biomarker_labels[result$outcome],
      Cohort = result$cohort,
      N = result$n,
      ElasticNet_R2 = result$r2_min,
      LinearReg_R2 = reg_result$model_summary$r.squared,
      LinearReg_AdjR2 = reg_result$model_summary$adj.r.squared,
      N_predictors = nrow(result$coef_min),
      stringsAsFactors = FALSE
    ))
  } else {
    # No regression was run (no variables selected)
    comparison_table <- rbind(comparison_table, data.frame(
      Biomarker = biomarker_labels[result$outcome],
      Cohort = result$cohort,
      N = result$n,
      ElasticNet_R2 = result$r2_min,
      LinearReg_R2 = NA,
      LinearReg_AdjR2 = NA,
      N_predictors = nrow(result$coef_min),
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(comparison_table, 
          "regression_results/elasticnet_vs_linear_comparison.csv", 
          row.names = FALSE)
cat("Saved: regression_results/elasticnet_vs_linear_comparison.csv\n")

cat("\nComparison Table:\n")
print(comparison_table)

cat("\n==========================================\n")
cat("REGRESSION ANALYSIS COMPLETE\n")
cat("==========================================\n")
cat(sprintf("Total elastic net models: %d\n", length(results_all)))
cat(sprintf("Total regression models: %d\n", length(regression_results_all)))
cat("Results saved in: regression_results/\n")


















#### Using far more variables 

#### Elastic net models 
library(readxl)
library(tidyr)
library(stringr)
library(dplyr)
library(glmnet)
library(ggplot2)
library(patchwork)

set.seed(123)

# Define biomarkers (outcomes)
biomarkers <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
                "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", 
                "ptau_217_avg_conc")

biomarker_labels <- c(
  "ab40_avg_conc" = "Aβ40",
  "ab42_avg_conc" = "Aβ42",
  "tau_avg_conc" = "Tau",
  "nfl_avg_conc" = "NfL",
  "gfap_avg_conc" = "GFAP",
  "ptau_181_avg_conc" = "pTau-181",
  "ptau_217_avg_conc" = "pTau-217"
)


#### Functions 

test_transformation <- function(x, var_name = "") {
  # Remove NAs
  x <- x[!is.na(x)]
  
  # Check if numeric
  if (!is.numeric(x)) {
    return(list(transform = "none", reason = "non_numeric"))
  }
  
  if (length(x) < 3) {
    return(list(transform = "none", reason = "insufficient_data"))
  }
  
  # Test raw data
  shapiro_raw <- tryCatch(shapiro.test(x)$p.value, error = function(e) NA)
  skew_raw <- (mean(x) - median(x)) / sd(x)
  
  # Test log-transformed (add small constant if zeros)
  if (min(x) <= 0) {
    x_log <- log(x - min(x) + 0.01)
  } else {
    x_log <- log(x)
  }
  
  shapiro_log <- tryCatch(shapiro.test(x_log)$p.value, error = function(e) NA)
  skew_log <- (mean(x_log) - median(x_log)) / sd(x_log)
  
  # Decision: use log if it improves normality and reduces skewness
  if (!is.na(shapiro_log) && !is.na(shapiro_raw)) {
    if (shapiro_log > shapiro_raw && abs(skew_log) < abs(skew_raw)) {
      return(list(
        transform = "log",
        shapiro_raw = shapiro_raw,
        shapiro_log = shapiro_log,
        skew_raw = abs(skew_raw),
        skew_log = abs(skew_log)
      ))
    }
  }
  
  return(list(
    transform = "none",
    shapiro_raw = shapiro_raw,
    shapiro_log = ifelse(is.na(shapiro_log), NA, shapiro_log),
    skew_raw = abs(skew_raw),
    skew_log = abs(skew_log)
  ))
}


apply_transformations <- function(data, transform_decisions) {
  data_transformed <- data
  
  for (var in names(transform_decisions)) {
    if (transform_decisions[[var]]$transform == "log") {
      x <- data_transformed[[var]]
      if (min(x, na.rm = TRUE) <= 0) {
        data_transformed[[var]] <- log(x - min(x, na.rm = TRUE) + 0.01)
      } else {
        data_transformed[[var]] <- log(x)
      }
    }
  }
  
  data_transformed
}


run_elastic_net <- function(data, outcome_var, predictor_vars, cohort_name) {
  
  cat(sprintf("\n=== Running Elastic Net: %s in %s ===\n", 
              outcome_var, cohort_name))
  
  # First, remove predictors that have too many NAs or no variance in this cohort
  valid_predictors <- predictor_vars[sapply(predictor_vars, function(var) {
    x <- data[[var]]
    # Keep if: numeric, less than 50% NAs, and has variance
    is.numeric(x) && sum(!is.na(x)) >= 10 && sd(x, na.rm = TRUE) > 0
  })]
  
  if (length(valid_predictors) == 0) {
    cat("No valid predictors for this cohort. Skipping.\n")
    return(NULL)
  }
  
  cat(sprintf("Valid predictors in cohort: %d (from %d total)\n", 
              length(valid_predictors), length(predictor_vars)))
  
  # Prepare data - drop rows with NA in outcome or any valid predictor
  complete_data <- data %>%
    select(all_of(c(outcome_var, valid_predictors))) %>%
    drop_na()
  
  if (nrow(complete_data) < 10) {
    cat(sprintf("Insufficient data after removing NAs (n = %d < 10). Skipping.\n", 
                nrow(complete_data)))
    return(NULL)
  }
  
  cat(sprintf("Sample size: %d\n", nrow(complete_data)))
  
  # Prepare matrices
  X <- as.matrix(complete_data[, valid_predictors])
  y <- complete_data[[outcome_var]]
  
  # Standardize predictors (glmnet does this internally, but we track it)
  X_scaled <- scale(X)
  
  # Set up cross-validation folds
  nFolds <- min(5, nrow(complete_data))
  foldid <- sample(rep(seq(nFolds), length.out = nrow(complete_data)))
  
  # Run elastic net with alpha = 0.5 (equal mix of L1 and L2)
  cv_fit <- cv.glmnet(
    x = X_scaled,
    y = y,
    family = "gaussian",
    alpha = 0.5,
    nfolds = nFolds,
    foldid = foldid,
    keep = TRUE
  )
  
  # Extract coefficients at lambda.min
  coef_min <- coef(cv_fit, s = "lambda.min")
  coef_df <- data.frame(
    Variable = rownames(coef_min),
    Coefficient = as.numeric(coef_min),
    stringsAsFactors = FALSE
  ) %>%
    filter(Variable != "(Intercept)", Coefficient != 0) %>%
    arrange(desc(abs(Coefficient)))
  
  # Extract coefficients at lambda.1se (more regularized)
  coef_1se <- coef(cv_fit, s = "lambda.1se")
  coef_1se_df <- data.frame(
    Variable = rownames(coef_1se),
    Coefficient = as.numeric(coef_1se),
    stringsAsFactors = FALSE
  ) %>%
    filter(Variable != "(Intercept)", Coefficient != 0) %>%
    arrange(desc(abs(Coefficient)))
  
  # Calculate R-squared
  predictions_min <- predict(cv_fit, newx = X_scaled, s = "lambda.min")
  r2_min <- 1 - sum((y - predictions_min)^2) / sum((y - mean(y))^2)
  
  predictions_1se <- predict(cv_fit, newx = X_scaled, s = "lambda.1se")
  r2_1se <- 1 - sum((y - predictions_1se)^2) / sum((y - mean(y))^2)
  
  cat(sprintf("Lambda.min: %.6f (R² = %.4f, %d variables)\n", 
              cv_fit$lambda.min, r2_min, nrow(coef_df)))
  cat(sprintf("Lambda.1se: %.6f (R² = %.4f, %d variables)\n", 
              cv_fit$lambda.1se, r2_1se, nrow(coef_1se_df)))
  
  list(
    cv_fit = cv_fit,
    coef_min = coef_df,
    coef_1se = coef_1se_df,
    r2_min = r2_min,
    r2_1se = r2_1se,
    n = nrow(complete_data),
    outcome = outcome_var,
    cohort = cohort_name
  )
}


plot_cv_results <- function(result, biomarker_label) {
  
  cv_fit <- result$cv_fit
  
  # Lambda plot
  lambda_df <- data.frame(
    lambda = cv_fit$lambda,
    mse = cv_fit$cvm,
    mse_lower = cv_fit$cvlo,
    mse_upper = cv_fit$cvup,
    nzero = cv_fit$nzero
  )
  
  p1 <- ggplot(lambda_df, aes(x = log(lambda), y = mse)) +
    geom_point(color = "#E64B35", size = 2) +
    geom_errorbar(aes(ymin = mse_lower, ymax = mse_upper), 
                  width = 0.1, color = "#E64B35", alpha = 0.5) +
    geom_vline(xintercept = log(cv_fit$lambda.min), 
               linetype = "dashed", color = "blue", size = 1) +
    geom_vline(xintercept = log(cv_fit$lambda.1se), 
               linetype = "dashed", color = "darkgreen", size = 1) +
    labs(
      title = sprintf("%s - %s", biomarker_label, result$cohort),
      x = "Log(Lambda)",
      y = "Mean Squared Error",
      subtitle = sprintf("λ.min (blue): %.4f | λ.1se (green): %.4f", 
                         cv_fit$lambda.min, cv_fit$lambda.1se)
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  
  # Coefficient plot at lambda.min
  if (nrow(result$coef_min) > 0) {
    coef_plot_df <- result$coef_min %>%
      arrange(Coefficient) %>%
      mutate(Variable = factor(Variable, levels = Variable))
    
    p2 <- ggplot(coef_plot_df, aes(x = Coefficient, y = Variable)) +
      geom_segment(aes(x = 0, xend = Coefficient, y = Variable, yend = Variable),
                   color = "grey50", size = 1) +
      geom_point(color = "#4DBBD5", size = 4) +
      geom_vline(xintercept = 0, linetype = "solid", color = "black") +
      labs(
        title = "Selected Features (λ.min)",
        x = "Standardized Coefficient",
        y = NULL,
        subtitle = sprintf("R² = %.3f, n = %d", result$r2_min, result$n)
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9)
      )
  } else {
    p2 <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = "No variables selected", size = 6) +
      theme_void()
  }
  
  p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))
}


# MODIFIED: Use data_varying instead of dat_clean
dat <- data_varying

# Define predictors (exclude biomarkers and group)
all_vars <- setdiff(names(dat), c(biomarkers, "group"))

# FILTER OUT NON-NUMERIC VARIABLES
numeric_vars <- sapply(all_vars, function(var) is.numeric(dat[[var]]))
predictor_vars <- all_vars[numeric_vars]

cat("Total variables (excluding biomarkers and group):", length(all_vars), "\n")
cat("Numeric predictor variables:", length(predictor_vars), "\n")

# Print non-numeric variables if any
non_numeric_vars <- all_vars[!numeric_vars]
if (length(non_numeric_vars) > 0) {
  cat("Non-numeric variables excluded:", paste(non_numeric_vars, collapse = ", "), "\n\n")
}

cat("Biomarkers (outcomes):", paste(biomarker_labels[biomarkers], collapse = ", "), "\n")
cat("Number of predictors:", length(predictor_vars), "\n")
cat("Predictors:", paste(predictor_vars, collapse = ", "), "\n\n")


transform_decisions <- list()

# Test biomarkers
cat("BIOMARKERS:\n")
for (var in biomarkers) {
  result <- test_transformation(dat[[var]], var)
  transform_decisions[[var]] <- result
  cat(sprintf("%-20s: %s", biomarker_labels[var], 
              ifelse(result$transform == "log", "LOG TRANSFORM", "NO TRANSFORM")))
  if (!is.null(result$shapiro_raw)) {
    cat(sprintf(" (Shapiro raw: %.3f, log: %.3f)\n", 
                result$shapiro_raw, result$shapiro_log))
  } else {
    cat("\n")
  }
}

# Test predictors
cat("\nPREDICTORS:\n")
for (var in predictor_vars) {
  result <- test_transformation(dat[[var]], var)
  transform_decisions[[var]] <- result
  cat(sprintf("%-25s: %s\n", var, 
              ifelse(result$transform == "log", "LOG TRANSFORM", "NO TRANSFORM")))
}

# Apply transformations
dat_transformed <- apply_transformations(dat, transform_decisions)

cat("\nTransformations applied.\n")
cat(sprintf("Log-transformed variables: %d\n", 
            sum(sapply(transform_decisions, function(x) x$transform == "log"))))



#Run elastic net models

# Initialize results lists HERE - BEFORE the loop
results_all <- list()
plot_list <- list()

# Define cohorts
cohorts <- list(
  "All" = dat_transformed,
  "T1D_only" = dat_transformed %>% filter(group == "Type 1 Diabetes")
)

# Check what values exist in the group column
cat("\n=== CHECKING GROUP VARIABLE ===\n")
cat("Unique values in 'group' column:\n")
print(unique(dat_transformed$group))
cat(sprintf("Total rows in full dataset: %d\n", nrow(dat_transformed)))
cat(sprintf("Rows where group == 'Type 1 Diabetes': %d\n", 
            sum(dat_transformed$group == "Type 1 Diabetes", na.rm = TRUE)))
cat("\n")

counter <- 1

for (cohort_name in names(cohorts)) {
  cohort_data <- cohorts[[cohort_name]]
  
  cat(sprintf("\n--- COHORT: %s (n = %d) ---\n", cohort_name, nrow(cohort_data)))
  
  for (biomarker in biomarkers) {
    
    result <- tryCatch({
      run_elastic_net(
        data = cohort_data,
        outcome_var = biomarker,
        predictor_vars = predictor_vars,
        cohort_name = cohort_name
      )
    }, error = function(e) {
      cat(sprintf("ERROR in %s - %s: %s\n", biomarker, cohort_name, e$message))
      return(NULL)
    })
    
    if (!is.null(result)) {
      results_all[[counter]] <- result
      
      # Create plot
      tryCatch({
        p <- plot_cv_results(result, biomarker_labels[biomarker])
        plot_list[[counter]] <- p
      }, error = function(e) {
        cat(sprintf("ERROR creating plot for %s - %s: %s\n", 
                    biomarker, cohort_name, e$message))
      })
      
      counter <- counter + 1
    }
  }
}

cat(sprintf("\n\nTotal models successfully run: %d\n", length(results_all)))

# MODIFIED: Create output directory with new name
dir.create("elastic_net_results_largeranalysis", showWarnings = FALSE)

# MODIFIED: Save individual plots with new naming
for (i in seq_along(plot_list)) {
  result <- results_all[[i]]
  filename <- sprintf("elastic_net_results_largeranalysis/%s_%s_plot_largeranalysis.png",
                      result$outcome, 
                      gsub(" ", "_", result$cohort))
  
  ggsave(filename, plot_list[[i]], width = 12, height = 6, dpi = 300)
  cat(sprintf("Saved: %s\n", filename))
}

# MODIFIED: Save coefficients with new naming
for (i in seq_along(results_all)) {
  result <- results_all[[i]]
  
  # Lambda.min coefficients
  if (nrow(result$coef_min) > 0) {
    filename <- sprintf("elastic_net_results_largeranalysis/%s_%s_coefficients_lambda_min_largeranalysis.csv",
                        result$outcome, 
                        gsub(" ", "_", result$cohort))
    
    write.csv(result$coef_min, filename, row.names = FALSE)
    cat(sprintf("Saved: %s\n", filename))
  }
  
  # Lambda.1se coefficients
  if (nrow(result$coef_1se) > 0) {
    filename <- sprintf("elastic_net_results_largeranalysis/%s_%s_coefficients_lambda_1se_largeranalysis.csv",
                        result$outcome, 
                        gsub(" ", "_", result$cohort))
    
    write.csv(result$coef_1se, filename, row.names = FALSE)
  }
}

# MODIFIED: Create summary table with new naming
summary_df <- data.frame(
  Biomarker = sapply(results_all, function(x) biomarker_labels[x$outcome]),
  Cohort = sapply(results_all, function(x) x$cohort),
  N = sapply(results_all, function(x) x$n),
  Lambda_min = sapply(results_all, function(x) x$cv_fit$lambda.min),
  R2_min = sapply(results_all, function(x) x$r2_min),
  N_vars_min = sapply(results_all, function(x) nrow(x$coef_min)),
  Lambda_1se = sapply(results_all, function(x) x$cv_fit$lambda.1se),
  R2_1se = sapply(results_all, function(x) x$r2_1se),
  N_vars_1se = sapply(results_all, function(x) nrow(x$coef_1se))
)

write.csv(summary_df, "elastic_net_results_largeranalysis/summary_all_models_largeranalysis.csv", row.names = FALSE)
cat("Saved: elastic_net_results_largeranalysis/summary_all_models_largeranalysis.csv\n")

cat("\n==========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==========================================\n")
cat(sprintf("Total models run: %d\n", length(results_all)))
print(summary_df)


# ==========================================
# LINEAR REGRESSION ON SELECTED VARIABLES
# ==========================================

cat("\n==========================================\n")
cat("RUNNING LINEAR REGRESSIONS ON SELECTED VARIABLES\n")
cat("==========================================\n\n")

# Create output directory
dir.create("regression_results_largeranalysis", showWarnings = FALSE)

regression_results_all <- list()

for (i in seq_along(results_all)) {
  result <- results_all[[i]]
  
  # Only run regression if variables were selected
  if (nrow(result$coef_min) > 0) {
    
    cat(sprintf("\n=== Linear Regression: %s in %s ===\n", 
                biomarker_labels[result$outcome], result$cohort))
    
    selected_vars <- result$coef_min$Variable
    
    # Get the appropriate cohort data - MORE FLEXIBLE MATCHING
    if (result$cohort == "All") {
      cohort_data <- dat_transformed
    } else if (result$cohort == "T1D_only") {
      cohort_data <- dat_transformed %>% filter(group == "Type 1 Diabetes")
    } else {
      # Fallback for any other cohort name
      cohort_data <- dat_transformed
      cat(sprintf("Warning: Unknown cohort '%s', using all data\n", result$cohort))
    }
    
    # Check if selected variables exist in cohort_data
    missing_vars <- selected_vars[!selected_vars %in% names(cohort_data)]
    if (length(missing_vars) > 0) {
      cat(sprintf("Warning: Variables not found in cohort data: %s\n", 
                  paste(missing_vars, collapse = ", ")))
      selected_vars <- selected_vars[selected_vars %in% names(cohort_data)]
    }
    
    if (length(selected_vars) == 0) {
      cat("No valid variables for regression. Skipping.\n")
      regression_results_all[[i]] <- NULL
      next
    }
    
    # Prepare regression data
    reg_data <- cohort_data %>%
      select(all_of(c(result$outcome, selected_vars))) %>%
      drop_na()
    
    if (nrow(reg_data) < length(selected_vars) + 2) {
      cat(sprintf("Insufficient data for regression (n = %d, need > %d). Skipping.\n",
                  nrow(reg_data), length(selected_vars) + 1))
      regression_results_all[[i]] <- NULL
      next
    }
    
    cat(sprintf("Variables: %s\n", paste(selected_vars, collapse = ", ")))
    cat(sprintf("Sample size: %d\n", nrow(reg_data)))
    
    # Build formula
    formula_str <- sprintf("%s ~ %s", result$outcome, paste(selected_vars, collapse = " + "))
    
    # Run linear regression with error handling
    tryCatch({
      lm_fit <- lm(as.formula(formula_str), data = reg_data)
      lm_summary <- summary(lm_fit)
      
      cat(sprintf("R²: %.4f, Adjusted R²: %.4f\n", 
                  lm_summary$r.squared, lm_summary$adj.r.squared))
      
      # Extract coefficients
      coef_table <- as.data.frame(lm_summary$coefficients)
      coef_table$Variable <- rownames(coef_table)
      rownames(coef_table) <- NULL
      coef_table <- coef_table[, c("Variable", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
      names(coef_table) <- c("Variable", "Estimate", "Std_Error", "t_value", "p_value")
      
      # Save to list
      regression_results_all[[i]] <- list(
        outcome = result$outcome,
        cohort = result$cohort,
        n = nrow(reg_data),
        model = lm_fit,
        model_summary = lm_summary,
        coefficients = coef_table,
        selected_vars = selected_vars
      )
      
      # Save coefficients to CSV
      filename <- sprintf("regression_results_largeranalysis/%s_%s_regression_coefficients_largeranalysis.csv",
                          result$outcome, 
                          gsub(" ", "_", result$cohort))
      write.csv(coef_table, filename, row.names = FALSE)
      cat(sprintf("Saved: %s\n", filename))
      
    }, error = function(e) {
      cat(sprintf("ERROR in regression: %s\n", e$message))
      regression_results_all[[i]] <<- NULL
    })
    
  } else {
    cat(sprintf("\n=== Skipping %s in %s (no variables selected) ===\n", 
                biomarker_labels[result$outcome], result$cohort))
    regression_results_all[[i]] <- NULL
  }
}

# Remove NULL entries
regression_results_all <- regression_results_all[!sapply(regression_results_all, is.null)]

cat(sprintf("\nTotal regression models run: %d\n", length(regression_results_all)))











#Follow-up Analysis

# ==========================================
# COMPARISON: ELASTIC NET VS LINEAR REGRESSION
# ==========================================

# Only run comparison if regression_results_all exists
if (exists("regression_results_all")) {
  
  cat("\n==========================================\n")
  cat("ELASTIC NET VS LINEAR REGRESSION COMPARISON\n")
  cat("==========================================\n\n")
  
  # Create output directory for comparison
  dir.create("regression_results_largeranalysis", showWarnings = FALSE)
  
  comparison_table <- data.frame(
    Biomarker = character(),
    Cohort = character(),
    N = integer(),
    ElasticNet_R2 = numeric(),
    LinearReg_R2 = numeric(),
    LinearReg_AdjR2 = numeric(),
    N_predictors = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(results_all)) {
    result <- results_all[[i]]
    
    # Find matching regression result
    matching_idx <- which(sapply(regression_results_all, function(x) 
      x$outcome == result$outcome && x$cohort == result$cohort))
    
    if (length(matching_idx) > 0) {
      reg_result <- regression_results_all[[matching_idx[1]]]
      
      comparison_table <- rbind(comparison_table, data.frame(
        Biomarker = biomarker_labels[result$outcome],
        Cohort = result$cohort,
        N = result$n,
        ElasticNet_R2 = result$r2_min,
        LinearReg_R2 = reg_result$model_summary$r.squared,
        LinearReg_AdjR2 = reg_result$model_summary$adj.r.squared,
        N_predictors = nrow(result$coef_min),
        stringsAsFactors = FALSE
      ))
    } else {
      # No regression was run (no variables selected)
      comparison_table <- rbind(comparison_table, data.frame(
        Biomarker = biomarker_labels[result$outcome],
        Cohort = result$cohort,
        N = result$n,
        ElasticNet_R2 = result$r2_min,
        LinearReg_R2 = NA,
        LinearReg_AdjR2 = NA,
        N_predictors = nrow(result$coef_min),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # MODIFIED: Save comparison table with new naming
  write.csv(comparison_table, 
            "regression_results_largeranalysis/elasticnet_vs_linear_comparison_largeranalysis.csv", 
            row.names = FALSE)
  cat("Saved: regression_results_largeranalysis/elasticnet_vs_linear_comparison_largeranalysis.csv\n")
  
  cat("\nComparison Table:\n")
  print(comparison_table)
  
  cat("\n==========================================\n")
  cat("REGRESSION ANALYSIS COMPLETE\n")
  cat("==========================================\n")
  cat(sprintf("Total elastic net models: %d\n", length(results_all)))
  cat(sprintf("Total regression models: %d\n", length(regression_results_all)))
  cat("Results saved in: regression_results_largeranalysis/\n")
  
} else {
  cat("\n==========================================\n")
  cat("SKIPPING COMPARISON SECTION\n")
  cat("==========================================\n")
  cat("Note: regression_results_all not found.\n")
  cat("Run linear regression analysis first if you want to compare results.\n")
  cat(sprintf("Elastic net analysis complete with %d models.\n", length(results_all)))
}





cat("\n==========================================\n")
cat("RUNNING LINEAR REGRESSIONS WITH LAMBDA.1SE\n")
cat("For models where lambda.min selected too many variables\n")
cat("==========================================\n\n")

# Find Tau T1D_only and pTau-181 T1D_only results
special_models <- list(
  list(outcome = "tau_avg_conc", cohort = "T1D_only", biomarker = "Tau"),
  list(outcome = "ptau_181_avg_conc", cohort = "T1D_only", biomarker = "pTau-181")
)

for (model_spec in special_models) {
  # Find the result
  result_idx <- which(sapply(results_all, function(x) 
    x$outcome == model_spec$outcome && x$cohort == model_spec$cohort))
  
  if (length(result_idx) == 0) {
    cat(sprintf("No results found for %s in %s\n", model_spec$biomarker, model_spec$cohort))
    next
  }
  
  result <- results_all[[result_idx[1]]]
  
  # Check if lambda.1se selected any variables
  if (nrow(result$coef_1se) == 0) {
    cat(sprintf("\n=== %s in %s ===\n", model_spec$biomarker, model_spec$cohort))
    cat("No variables selected at lambda.1se. Skipping.\n")
    next
  }
  
  cat(sprintf("\n=== Linear Regression (lambda.1se): %s in %s ===\n", 
              model_spec$biomarker, model_spec$cohort))
  
  selected_vars <- result$coef_1se$Variable
  
  # Get T1D cohort data
  cohort_data <- dat_transformed %>% filter(group == "Type 1 Diabetes")
  
  # Prepare regression data
  reg_data <- cohort_data %>%
    select(all_of(c(result$outcome, selected_vars))) %>%
    drop_na()
  
  cat(sprintf("Variables selected (n=%d): %s\n", 
              length(selected_vars), paste(selected_vars, collapse = ", ")))
  cat(sprintf("Sample size: %d\n", nrow(reg_data)))
  
  if (nrow(reg_data) < length(selected_vars) + 2) {
    cat(sprintf("Insufficient data for regression (n = %d, need > %d). Skipping.\n",
                nrow(reg_data), length(selected_vars) + 1))
    next
  }
  
  # Build formula
  formula_str <- sprintf("%s ~ %s", result$outcome, paste(selected_vars, collapse = " + "))
  
  # Run linear regression
  tryCatch({
    lm_fit <- lm(as.formula(formula_str), data = reg_data)
    lm_summary <- summary(lm_fit)
    
    cat(sprintf("R²: %.4f, Adjusted R²: %.4f\n", 
                lm_summary$r.squared, lm_summary$adj.r.squared))
    cat(sprintf("F-statistic: %.2f on %d and %d DF, p-value: %.4f\n",
                lm_summary$fstatistic[1], lm_summary$fstatistic[2], 
                lm_summary$fstatistic[3], 
                pf(lm_summary$fstatistic[1], lm_summary$fstatistic[2], 
                   lm_summary$fstatistic[3], lower.tail = FALSE)))
    
    # Extract and display coefficients
    coef_table <- as.data.frame(lm_summary$coefficients)
    coef_table$Variable <- rownames(coef_table)
    rownames(coef_table) <- NULL
    coef_table <- coef_table[, c("Variable", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
    names(coef_table) <- c("Variable", "Estimate", "Std_Error", "t_value", "p_value")
    
    cat("\nCoefficients:\n")
    print(coef_table, row.names = FALSE, digits = 4)
    
    # Save coefficients to CSV
    filename <- sprintf("regression_results_largeranalysis/%s_%s_regression_lambda1se_largeranalysis.csv",
                        result$outcome, 
                        gsub(" ", "_", result$cohort))
    write.csv(coef_table, filename, row.names = FALSE)
    cat(sprintf("\nSaved: %s\n", filename))
    
  }, error = function(e) {
    cat(sprintf("ERROR in regression: %s\n", e$message))
  })
}










##### Correlation matrix 
library(corrplot)

# Prepare data for correlation - include biomarkers and predictors
cor_data <- dat_transformed %>%
  select(all_of(c(biomarkers, predictor_vars))) %>%
  select(where(is.numeric))

# Calculate correlation matrix with pairwise complete observations
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

# Check for any issues
cat(sprintf("Correlation matrix dimensions: %d x %d\n", nrow(cor_matrix), ncol(cor_matrix)))
cat(sprintf("Variables included: %d\n", ncol(cor_data)))

# Create correlation plot with hierarchical clustering - LARGER SIZE
png("elastic_net_results_largeranalysis/correlation_plot_clustered_largeranalysis.png", 
    width = 4500, height = 4500, res = 300)

corrplot(cor_matrix, 
         method = "color",           # Color-coded correlations
         type = "full",              # Show full matrix
         order = "hclust",           # Hierarchical clustering
         hclust.method = "complete", # Complete linkage clustering
         addrect = 5,                # Add 5 rectangles around clusters
         tl.col = "black",           # Text label color
         tl.cex = 0.8,               # Text label size (increased)
         tl.srt = 45,                # Text label rotation
         col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200), # Color palette
         addCoef.col = NULL,         # Don't add correlation coefficients (too crowded)
         cl.cex = 1.0,               # Color legend text size (increased)
         title = "Correlation Matrix with Hierarchical Clustering\n(Transformed Variables Including Biomarkers)",
         mar = c(0, 0, 3, 0))        # Margins

dev.off()

cat("Saved: elastic_net_results_largeranalysis/correlation_plot_clustered_largeranalysis.png\n")

# Also create a version with only biomarkers for clearer visualization
cor_matrix_biomarkers <- cor(dat_transformed[, biomarkers], 
                             use = "pairwise.complete.obs")

png("elastic_net_results_largeranalysis/correlation_plot_biomarkers_largeranalysis.png", 
    width = 2000, height = 2000, res = 300)

corrplot(cor_matrix_biomarkers, 
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",      # Add correlation coefficients for biomarkers
         number.cex = 0.8,
         tl.col = "black",
         tl.cex = 1,
         tl.srt = 45,
         col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
         cl.cex = 1,
         title = "Biomarker Correlations with Clustering\n(Transformed Variables)",
         mar = c(0, 0, 2, 0))

dev.off()


