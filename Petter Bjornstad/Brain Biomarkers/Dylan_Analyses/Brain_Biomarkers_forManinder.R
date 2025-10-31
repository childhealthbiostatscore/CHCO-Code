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

# Reshape data to long format
# Assuming your data frame is called 'df' with a column 'group' 
# containing "Type 1 Diabetes" and "Lean Control"
# df_long <- df %>%
#   select(group, all_of(qx_var)) %>%
#   pivot_longer(cols = all_of(qx_var),
#                names_to = "biomarker",
#                values_to = "concentration") %>%
#   mutate(biomarker_label = biomarker_labels[biomarker])

# Example with simulated data (replace with your actual data)
set.seed(123)
n_per_group <- 50

df_long <- data.frame(
  group = rep(c("Lean Control", "Type 1 Diabetes"), each = n_per_group * length(qx_var)),
  biomarker = rep(rep(qx_var, each = n_per_group), 2),
  concentration = c(
    rnorm(n_per_group, 150, 30),  # ab40 - Control
    rnorm(n_per_group, 45, 10),   # ab42 - Control
    rexp(n_per_group, 1/12) + 5,  # tau - Control (skewed)
    rexp(n_per_group, 1/20) + 5,  # nfl - Control (skewed)
    rnorm(n_per_group, 180, 40),  # gfap - Control
    rexp(n_per_group, 1/2.5),     # ptau_181 - Control (skewed)
    rexp(n_per_group, 1/1.8),     # ptau_217 - Control (skewed)
    rnorm(n_per_group, 160, 35),  # ab40 - T1D
    rnorm(n_per_group, 42, 9),    # ab42 - T1D
    rexp(n_per_group, 1/14) + 5,  # tau - T1D (skewed)
    rexp(n_per_group, 1/24) + 5,  # nfl - T1D (skewed)
    rnorm(n_per_group, 210, 45),  # gfap - T1D
    rexp(n_per_group, 1/3.0),     # ptau_181 - T1D (skewed)
    rexp(n_per_group, 1/2.1)      # ptau_217 - T1D (skewed)
  )
)

cat("\n==========================================\n")
cat("DATA CHECK\n")
cat("==========================================\n")
cat(sprintf("Total rows in dataset: %d\n", nrow(df_long)))
cat(sprintf("Unique biomarkers: %d\n", length(unique(df_long$biomarker))))
cat("Biomarkers present:\n")
print(unique(df_long$biomarker))
cat("\n")

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

# Set color palette
group_colors <- c("Lean Control" = "#4DBBD5", "Type 1 Diabetes" = "#E64B35")

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

# Alternative: Create a single faceted plot with p-values
# First, create a dataframe with p-value positions for each biomarker
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















#Clinical characteristics (eGFR, UACR, HBA1C, Clamp, PET, DEXA, glycemia, insulin sensitivity)

traits_of_interest <- c('acr_u', 
                        'eGFR_bedside_Schwartz', 'eGFR_CKD_epi', 'eGFR_fas_cr', 'eGFR_fas_cr_cysc','eGFR_Zap','eGFR_Schwartz', 
                        'hba1c')


dat <- dat %>% dplyr::select(record_id, group, age, sex, bmi, hba1c, study, all_of(qx_var))






#### Elastic net models 





















