######FSOC with immune analysis 


library(dplyr)
library(stringr)




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))



# ============================================================================
# PROTEOMICS ANALYSIS BY FSOC STATUS
# Analyzing inflammatory markers in relation to medullary FSOC
# ============================================================================

library(tidyverse)
library(limma)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(readxl)
library(gtsummary)
library(ggrepel)

# ============================================================================
# SETUP
# ============================================================================

OUTPUT_DIR <- "C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani/FSOC_Proteomics"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

# Load harmonized data
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

# Aggregate by subject
dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

# Calculate FSOC averages
dat <- dat %>%
  mutate(
    # Average medullary FSOC from left and right
    medullary_fsoc_abs = case_when(
      !is.na(fsoc_l_medulla) & !is.na(fsoc_r_medulla) ~ 
        (fsoc_l_medulla + fsoc_r_medulla) / 2,
      !is.na(fsoc_l_medulla) ~ fsoc_l_medulla,
      !is.na(fsoc_r_medulla) ~ fsoc_r_medulla,
      TRUE ~ NA_real_
    ),
    # Average whole kidney FSOC
    whole_kidney_fsoc_abs = case_when(
      !is.na(fsoc_l_kidney) & !is.na(fsoc_r_kidney) ~ 
        (fsoc_l_kidney + fsoc_r_kidney) / 2,
      !is.na(fsoc_l_kidney) ~ fsoc_l_kidney,
      !is.na(fsoc_r_kidney) ~ fsoc_r_kidney,
      TRUE ~ NA_real_
    )
  )

# Load data dictionary
data_dictionary <- readxl::read_xlsx('C:/Users/netio/Downloads/data_dictionary_master.xlsx')

# Get proteomics variables
form_names <- unique(data_dictionary$form_name)
proteo <- form_names[str_which(form_names, pattern = 'proteom')]

variables_class <- c('proteomics', 'metabolomics', 'metabolomics_blood_raw',
                     'az_urine_metabolites', 'metabolomics_aq')

data_dictionary_small <- data_dictionary %>% 
  filter(form_name %in% variables_class)

# Get existing proteomics columns
existing_cols <- intersect(data_dictionary_small$variable_name, names(dat))

cat("Found:", length(existing_cols), "proteomics/metabolomics variables\n")

# ============================================================================
# 2. FILTER TO SUBJECTS WITH VALID FSOC AND PROTEOMICS
# ============================================================================

# Filter to baseline visit with valid FSOC and proteomics
dat_fsoc_proteomics <- dat %>%
  filter(visit == 'baseline') %>%
  filter(!is.na(medullary_fsoc_abs)) %>%
  filter(medullary_fsoc_abs >= 0) %>%  # Exclude negative FSOC
  filter(whole_kidney_fsoc_abs < 15 | is.na(whole_kidney_fsoc_abs)) %>%  # Exclude extreme WK values
  # Filter to subjects with at least some proteomics data
  filter(rowSums(!is.na(select(., all_of(existing_cols)))) > 10)

cat("\nSubjects with valid FSOC and proteomics:", nrow(dat_fsoc_proteomics), "\n")

# Classify FSOC status
dat_fsoc_proteomics <- dat_fsoc_proteomics %>%
  mutate(
    # Binary classification at median
    fsoc_binary = ifelse(medullary_fsoc_abs < median(medullary_fsoc_abs, na.rm = TRUE), 
                         "Impaired", "Normal"),
    # Tertile classification
    fsoc_tertile = case_when(
      medullary_fsoc_abs <= quantile(medullary_fsoc_abs, 1/3, na.rm = TRUE) ~ "Low",
      medullary_fsoc_abs <= quantile(medullary_fsoc_abs, 2/3, na.rm = TRUE) ~ "Medium",
      TRUE ~ "High"
    ),
    # Quartile classification for extreme comparison
    fsoc_quartile = case_when(
      medullary_fsoc_abs <= quantile(medullary_fsoc_abs, 0.25, na.rm = TRUE) ~ "Q1_Impaired",
      medullary_fsoc_abs >= quantile(medullary_fsoc_abs, 0.75, na.rm = TRUE) ~ "Q4_Normal",
      TRUE ~ "Q2-Q3_Intermediate"
    ),
    fsoc_binary = factor(fsoc_binary, levels = c("Normal", "Impaired")),
    fsoc_tertile = factor(fsoc_tertile, levels = c("Low", "Medium", "High")),
    fsoc_quartile = factor(fsoc_quartile, levels = c("Q1_Impaired", "Q2-Q3_Intermediate", "Q4_Normal"))
  )

cat("\nFSOC Classification:\n")
cat("Binary:\n")
print(table(dat_fsoc_proteomics$fsoc_binary))
cat("\nTertile:\n")
print(table(dat_fsoc_proteomics$fsoc_tertile))
cat("\nQuartile:\n")
print(table(dat_fsoc_proteomics$fsoc_quartile))

# ============================================================================
# 3. DEMOGRAPHICS TABLE BY FSOC STATUS
# ============================================================================

demographics <- dat_fsoc_proteomics %>%
  select(age, sex, race_ethnicity, bmi, hba1c, group, 
         medullary_fsoc_abs, whole_kidney_fsoc_abs, fsoc_binary) %>%
  tbl_summary(
    by = fsoc_binary,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      medullary_fsoc_abs ~ "continuous",
      whole_kidney_fsoc_abs ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      group ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      medullary_fsoc_abs ~ 3,
      whole_kidney_fsoc_abs ~ 3,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      group ~ "Disease Group",
      medullary_fsoc_abs ~ "Medullary FSOC (s⁻¹)",
      whole_kidney_fsoc_abs ~ "Whole Kidney FSOC (s⁻¹)"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "wilcox.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**FSOC Status**")

print(demographics)

demographics %>%
  as_gt() %>%
  gtsave(file.path(OUTPUT_DIR, "fsoc_proteomics_demographics.png"), 
         vwidth = 1400, vheight = 800)

# ============================================================================
# 4. DEFINE INFLAMMATORY/KIDNEY INJURY MARKERS
# ============================================================================

# Key markers to focus on (if available in proteomics)
inflammatory_markers_priority <- c(
  # Kidney injury
  "HAVCR1", "KIM1", "TIMD1",  # KIM-1
  "LCN2", "NGAL",              # NGAL
  "CST3", "CYSC",              # Cystatin C
  "CLU",                        # Clusterin
  
  # Cytokines
  "IL6", "IL1B", "TNF", "TNFA",
  "IL8", "CXCL8",
  "IL10",
  "IFNG", "IFN",
  
  # Chemokines
  "CCL2", "MCP1",              # MCP-1
  "CCL5", "RANTES",
  
  # Acute phase
  "CRP", 
  "SAA1", "SAA2", "SAA",
  "FGB", "FIBRINOGEN",
  "HP", "HPT", "HAPTOGLOBIN",
  
  # Complement
  "C3", "C5", "CFB",
  
  # Endothelial
  "ICAM1", "CD54",
  "VCAM1", "CD106",
  "SELE", "SELECTIN",
  
  # Other
  "MPO", "MYELOPEROXIDASE",
  "PTX3", "PENTRAXIN",
  "UMOD", "UROMODULIN",
  "SPP1", "OSTEOPONTIN",
  "FGF23"
)

# Find which markers are available (case-insensitive)
available_inflammatory <- existing_cols[
  toupper(existing_cols) %in% toupper(inflammatory_markers_priority)
]

cat("\n=== Inflammatory Markers Available ===\n")
cat("Found:", length(available_inflammatory), "of", 
    length(inflammatory_markers_priority), "priority markers\n")
if (length(available_inflammatory) > 0) {
  print(available_inflammatory)
}

# If specific markers not available, use all proteomics
proteomics_to_test <- if (length(available_inflammatory) > 5) {
  available_inflammatory
} else {
  # Get all proteomics variables
  proteo_vars <- data_dictionary_small %>%
    filter(form_name == 'proteomics') %>%
    pull(variable_name)
  intersect(proteo_vars, existing_cols)
}

cat("\nTotal proteomics variables to analyze:", length(proteomics_to_test), "\n")

# ============================================================================
# 5. ANALYSIS 1: BINARY FSOC COMPARISON (IMPAIRED VS NORMAL)
# ============================================================================

cat("\n========================================\n")
cat("BINARY FSOC ANALYSIS\n")
cat("========================================\n\n")

results_binary <- data.frame(
  variable = character(),
  variable_label = character(),
  test_used = character(),
  p_value = numeric(),
  mean_impaired = numeric(),
  mean_normal = numeric(),
  median_impaired = numeric(),
  median_normal = numeric(),
  stringsAsFactors = FALSE
)

for(var in proteomics_to_test) {
  
  var_label <- data_dictionary$label[data_dictionary$variable_name == var]
  if(length(var_label) == 0) var_label <- var
  
  temp_data <- dat_fsoc_proteomics %>%
    select(fsoc_binary, all_of(var), age, sex, group) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value), !is.na(fsoc_binary))
  
  if(nrow(temp_data) < 6 || length(unique(temp_data$fsoc_binary)) < 2) next
  
  impaired_data <- temp_data$value[temp_data$fsoc_binary == "Impaired"]
  normal_data <- temp_data$value[temp_data$fsoc_binary == "Normal"]
  
  if(length(impaired_data) < 3 || length(normal_data) < 3) next
  
  # Use Mann-Whitney for robustness
  test_result <- wilcox.test(value ~ fsoc_binary, data = temp_data)
  test_name <- "Mann-Whitney U"
  
  results_binary <- rbind(results_binary, data.frame(
    variable = var,
    variable_label = var_label,
    test_used = test_name,
    p_value = test_result$p.value,
    mean_impaired = mean(impaired_data, na.rm = TRUE),
    mean_normal = mean(normal_data, na.rm = TRUE),
    median_impaired = median(impaired_data, na.rm = TRUE),
    median_normal = median(normal_data, na.rm = TRUE)
  ))
}

# Adjust for multiple testing
results_binary <- results_binary %>%
  mutate(
    p_bonferroni = p.adjust(p_value, method = "bonferroni"),
    p_fdr = p.adjust(p_value, method = "fdr"),
    fold_change = mean_impaired / mean_normal,
    log2_fc = log2(fold_change),
    difference = mean_impaired - mean_normal
  ) %>%
  arrange(p_value)

cat("Top 20 proteins by p-value:\n")
print(head(results_binary %>% select(variable_label, p_value, p_fdr, 
                                     mean_impaired, mean_normal, fold_change), 20))

cat("\nSignificant at FDR < 0.05:", sum(results_binary$p_fdr < 0.05), "\n")
cat("Significant at FDR < 0.10:", sum(results_binary$p_fdr < 0.10), "\n")

write.csv(results_binary, file.path(OUTPUT_DIR, 'fsoc_binary_proteomics_results.csv'), 
          row.names = FALSE)

# ============================================================================
# 6. ANALYSIS 2: CORRELATION WITH CONTINUOUS FSOC
# ============================================================================

cat("\n========================================\n")
cat("CONTINUOUS FSOC CORRELATION ANALYSIS\n")
cat("========================================\n\n")

results_correlation <- data.frame(
  variable = character(),
  variable_label = character(),
  spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(var in proteomics_to_test) {
  
  var_label <- data_dictionary$label[data_dictionary$variable_name == var]
  if(length(var_label) == 0) var_label <- var
  
  temp_data <- dat_fsoc_proteomics %>%
    select(medullary_fsoc_abs, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value), !is.na(medullary_fsoc_abs))
  
  if(nrow(temp_data) < 10) next
  
  cor_test <- cor.test(temp_data$medullary_fsoc_abs, temp_data$value, 
                       method = "spearman")
  
  results_correlation <- rbind(results_correlation, data.frame(
    variable = var,
    variable_label = var_label,
    spearman_rho = cor_test$estimate,
    p_value = cor_test$p.value
  ))
}

results_correlation <- results_correlation %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "fdr"),
    abs_rho = abs(spearman_rho)
  ) %>%
  arrange(p_value)

cat("Top 20 correlations:\n")
print(head(results_correlation %>% select(variable_label, spearman_rho, 
                                          p_value, p_fdr), 20))

write.csv(results_correlation, file.path(OUTPUT_DIR, 'fsoc_continuous_correlations.csv'), 
          row.names = FALSE)

# ============================================================================
# 7. VISUALIZATIONS
# ============================================================================

# 7A. Volcano Plot
sig_binary <- results_binary %>% filter(p_fdr < 0.10)

p_volcano <- ggplot(results_binary, 
                    aes(x = log2_fc, y = -log10(p_value), 
                        color = p_fdr < 0.10)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"),
                     name = "FDR < 0.10") +
  theme_bw(base_size = 12) +
  labs(x = "Log2 Fold Change (Impaired / Normal FSOC)",
       y = "-log10(P-value)",
       title = "Proteomics: Impaired vs Normal Medullary FSOC") +
  theme(legend.position = "top")

if (nrow(sig_binary) > 0 && nrow(sig_binary) <= 20) {
  p_volcano <- p_volcano +
    geom_text_repel(data = sig_binary,
                    aes(label = variable_label), 
                    size = 3, max.overlaps = 20)
}

ggsave(file.path(OUTPUT_DIR, "volcano_plot_fsoc.pdf"), 
       p_volcano, width = 10, height = 8)

# 7B. Top Markers Boxplots
if (nrow(sig_binary) > 0) {
  top_n <- min(9, nrow(sig_binary))
  top_vars <- head(results_binary$variable, top_n)
  
  plot_list <- list()
  
  for(i in 1:length(top_vars)) {
    var <- top_vars[i]
    var_info <- results_binary %>% filter(variable == var)
    var_label <- var_info$variable_label
    p_val <- var_info$p_value
    p_fdr <- var_info$p_fdr
    
    plot_data <- dat_fsoc_proteomics %>%
      select(fsoc_binary, all_of(var)) %>%
      rename(value = all_of(var)) %>%
      filter(!is.na(value), !is.na(fsoc_binary))
    
    p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.4f", p_val))
    p_fdr_label <- ifelse(p_fdr < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p_fdr))
    
    p <- ggplot(plot_data, aes(x = fsoc_binary, y = value, fill = fsoc_binary)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
      scale_fill_manual(values = c("Impaired" = "#E6550D", "Normal" = "#3182BD")) +
      theme_bw(base_size = 10) +
      labs(title = var_label,
           subtitle = paste0(p_label, " | ", p_fdr_label),
           x = "FSOC Status",
           y = "Expression") +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold"),
            plot.subtitle = element_text(size = 7))
    
    plot_list[[i]] <- p
  }
  
  png(file.path(OUTPUT_DIR, 'top_markers_by_fsoc.png'), 
      height = 12, width = 12, units = 'in', res = 300)
  grid.arrange(grobs = plot_list, ncol = 3)
  dev.off()
}

# 7C. Correlation Scatter Plots for Top Hits
top_cors <- head(results_correlation, 9)

if (nrow(top_cors) > 0) {
  plot_list_cor <- list()
  
  for(i in 1:nrow(top_cors)) {
    var <- top_cors$variable[i]
    var_label <- top_cors$variable_label[i]
    rho <- top_cors$spearman_rho[i]
    p_val <- top_cors$p_value[i]
    
    plot_data <- dat_fsoc_proteomics %>%
      select(medullary_fsoc_abs, all_of(var), fsoc_binary) %>%
      rename(value = all_of(var)) %>%
      filter(!is.na(value), !is.na(medullary_fsoc_abs))
    
    p_label <- sprintf("ρ = %.3f, p = %.4f", rho, p_val)
    
    p <- ggplot(plot_data, aes(x = medullary_fsoc_abs, y = value)) +
      geom_point(aes(color = fsoc_binary), alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
      scale_color_manual(values = c("Impaired" = "#E6550D", "Normal" = "#3182BD"),
                         name = "FSOC") +
      theme_bw(base_size = 10) +
      labs(title = var_label,
           subtitle = p_label,
           x = "Medullary FSOC (s⁻¹)",
           y = "Expression") +
      theme(plot.title = element_text(size = 9, face = "bold"),
            plot.subtitle = element_text(size = 7),
            legend.position = "right")
    
    plot_list_cor[[i]] <- p
  }
  
  png(file.path(OUTPUT_DIR, 'top_correlations_scatter.png'), 
      height = 12, width = 14, units = 'in', res = 300)
  grid.arrange(grobs = plot_list_cor, ncol = 3)
  dev.off()
}

# ============================================================================
# 8. HEATMAP OF TOP MARKERS
# ============================================================================

if (nrow(sig_binary) > 0) {
  top_heatmap <- min(30, nrow(sig_binary))
  top_vars_heatmap <- head(results_binary$variable, top_heatmap)
  
  # Prepare matrix
  heatmap_data <- dat_fsoc_proteomics %>%
    select(record_id, fsoc_binary, group, sex, all_of(top_vars_heatmap)) %>%
    column_to_rownames("record_id")
  
  # Separate annotations from data
  annotation_col <- heatmap_data %>% select(fsoc_binary, group, sex)
  heatmap_matrix <- heatmap_data %>% 
    select(-fsoc_binary, -group, -sex) %>%
    as.matrix() %>%
    t()
  
  # Scale by row (z-score)
  heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
  
  # Annotation colors
  annot_colors <- list(
    fsoc_binary = c("Impaired" = "#E6550D", "Normal" = "#3182BD"),
    sex = c("Female" = "#E8A0A0", "Male" = "#8CB4D0"),
    group = c("Lean Control" = "#3182BD",
              "Obese Control" = "#9ECAE1",
              "Type 1 Diabetes" = "#FDAE6B",
              "Type 2 Diabetes" = "#E6550D",
              "PKD" = "#9E9AC8")
  )
  
  # Row labels (use variable labels)
  rownames(heatmap_matrix_scaled) <- sapply(rownames(heatmap_matrix_scaled), function(v) {
    label <- data_dictionary$label[data_dictionary$variable_name == v]
    if(length(label) == 0) return(v)
    return(label)
  })
  
  pdf(file.path(OUTPUT_DIR, "heatmap_top_markers.pdf"), width = 12, height = 10)
  pheatmap(heatmap_matrix_scaled,
           annotation_col = annotation_col,
           annotation_colors = annot_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_colnames = FALSE,
           fontsize_row = 8,
           main = "Top Proteomics Markers by FSOC Status\n(Z-score scaled)",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
  dev.off()
}

# ============================================================================
# 9. SUMMARY REPORT
# ============================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat("========================================\n\n")

cat("Total subjects analyzed:", nrow(dat_fsoc_proteomics), "\n")
cat("  - Impaired FSOC:", sum(dat_fsoc_proteomics$fsoc_binary == "Impaired"), "\n")
cat("  - Normal FSOC:", sum(dat_fsoc_proteomics$fsoc_binary == "Normal"), "\n\n")

cat("Proteomics variables tested:", nrow(results_binary), "\n")
cat("Significant (FDR < 0.05):", sum(results_binary$p_fdr < 0.05), "\n")
cat("Significant (FDR < 0.10):", sum(results_binary$p_fdr < 0.10), "\n")
cat("Nominal significance (p < 0.05):", sum(results_binary$p_value < 0.05), "\n\n")

cat("Top 5 upregulated in Impaired FSOC:\n")
up_reg <- results_binary %>% filter(log2_fc > 0) %>% head(5)
print(up_reg %>% select(variable_label, fold_change, p_value, p_fdr))

cat("\nTop 5 downregulated in Impaired FSOC:\n")
down_reg <- results_binary %>% filter(log2_fc < 0) %>% head(5)
print(down_reg %>% select(variable_label, fold_change, p_value, p_fdr))

cat("\n\nFiles saved to:", OUTPUT_DIR, "\n")
cat("  - fsoc_binary_proteomics_results.csv\n")
cat("  - fsoc_continuous_correlations.csv\n")
cat("  - volcano_plot_fsoc.pdf\n")
cat("  - top_markers_by_fsoc.png\n")
cat("  - top_correlations_scatter.png\n")
cat("  - heatmap_top_markers.pdf\n")
cat("  - fsoc_proteomics_demographics.png\n")





































