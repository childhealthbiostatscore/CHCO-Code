### Brain Imaging by CKD Staging - Split by Serum/Plasma
### Updated to show serum and plasma side-by-side

# Load libraries
library(oro.nifti)  
library(neurobase)  
library(R.matlab)    
library(ggplot2)     
library(dplyr)       
library(corrplot)   
library(psych)        
library(pls)  
library(caret) 
library(randomForest) 
library(igraph)
library(brainGraph)
library(R.matlab)
library(tidyverse)
library(patchwork)  # For combining plots side-by-side

bucket <- 'brain.mri'

# ============================================================
# DATA LOADING AND PREPARATION
# ============================================================

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

# ============================================================
# DETERMINE SAMPLE TYPE (SERUM/PLASMA) BY STUDY
# ============================================================

dat <- dat %>%
  mutate(
    sample_type = case_when(
      # CROCODILE participants (CRC IDs) - plasma
      grepl("^CRC", record_id, ignore.case = TRUE) ~ "plasma",
      # PENGUIN participants (PEN IDs) - plasma
      grepl("^PEN", record_id, ignore.case = TRUE) ~ "plasma",
      # RH2 participants - serum
      study == "RH2" ~ "serum",
      grepl("^RH2", record_id, ignore.case = TRUE) ~ "serum",
      # RENAL-HEIR studies - serum
      study == "RENAL-HEIR" ~ "serum",
      study == "RENAL-HEIRitage" ~ "serum",
      # Default to study-based assignment
      study == "CROCODILE" ~ "plasma",
      study == "PENGUIN" ~ "plasma",
      TRUE ~ NA_character_
    )
  )

# Check sample type assignment
cat("\nSample type distribution:\n")
print(table(dat$sample_type, useNA = "always"))

# MRI IDs (keeping for reference)
mri_ids <- c('CRC-10', 'CRC-11', 'CRC-12', 'CRC-13', 'CRC-26', 'CRC-39', 'CRC-46', 'CRC-51', 'CRC-53', 
             'CRC-55', 'CRC-58', 'CRC-60', 
             'RH2-01-O', 'RH2-03-O', 'RH2-08-T', 'RH2-10-L', 'RH2-11-O', 'RH2-13-O', 'RH2-16-O', 'RH2-17-L', 
             'RH2-18-O', 'RH2-19-T', 'RH2-22-T', 'RH2-24-L', 'RH2-27-L', 'RH2-28-L', 'RH2-29-L', 'RH2-33-L', 
             'RH2-34-O', 'RH2-35-T', 'RH2-38-T', 'RH2-39-O', 'RH2-41-T', 'RH2-42-T', 'RH2-43-T', 
             'RH2-44-T', 'RH2-45-T', 'RH2-48-T', 'RH2-49-T', 'RH2-50-L', 'RH2-52-T', 'RH2-53-T', 
             'RH2-55-T')

small_dat <- dat

# Fix the RH2-38-T record
small_dat$group[which(small_dat$record_id == 'RH2-38-T')] <- 'Obese Control'
small_dat$record_id[which(small_dat$record_id == 'RH2-38-T')] <- 'RH2-38-O'

# Define biomarkers
qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")

# ============================================================
# CKD CLASSIFICATION FUNCTIONS
# ============================================================

classify_ckd <- function(gfr = NULL, albumin_mg_g = NULL, albumin_mg_mmol = NULL, 
                         gfr_days_from_baseline = NULL, albumin_days_from_baseline = NULL, 
                         max_days_apart = 365) {
  
  if (!is.null(albumin_mg_mmol) && is.null(albumin_mg_g)) {
    albumin_mg_g <- albumin_mg_mmol * 8.84
  }
  
  gfr_available <- !is.null(gfr) && !is.na(gfr)
  albumin_available <- !is.null(albumin_mg_g) && !is.na(albumin_mg_g)
  
  days_available <- !is.null(gfr_days_from_baseline) && !is.null(albumin_days_from_baseline) && 
    !is.na(gfr_days_from_baseline) && !is.na(albumin_days_from_baseline)
  
  within_time_window <- TRUE
  days_apart <- NA
  criteria_met_days_from_baseline <- NA
  
  if (days_available && gfr_available && albumin_available) {
    days_apart <- abs(gfr_days_from_baseline - albumin_days_from_baseline)
    within_time_window <- days_apart <= max_days_apart
    if (within_time_window) {
      criteria_met_days_from_baseline <- min(gfr_days_from_baseline, albumin_days_from_baseline)
    }
  } else if ((days_available && (gfr_available || albumin_available)) || 
             (!days_available && gfr_available && albumin_available)) {
    if (gfr_available && !albumin_available && !is.null(gfr_days_from_baseline) && !is.na(gfr_days_from_baseline)) {
      criteria_met_days_from_baseline <- gfr_days_from_baseline
    } else if (!gfr_available && albumin_available && !is.null(albumin_days_from_baseline) && !is.na(albumin_days_from_baseline)) {
      criteria_met_days_from_baseline <- albumin_days_from_baseline
    }
  }
  
  if (!gfr_available && !albumin_available) {
    return(list(gfr_category = NA, gfr_description = NA, albumin_category = NA,
                albumin_description = NA, risk_level = NA, recommendation = NA,
                days_between_tests = NA, within_time_window = NA, criteria_met_days_from_baseline = NA))
  }
  
  if (gfr_available) {
    gfr_category <- case_when(
      gfr >= 90 ~ "G1", gfr >= 60 & gfr < 90 ~ "G2", gfr >= 45 & gfr < 60 ~ "G3a",
      gfr >= 30 & gfr < 45 ~ "G3b", gfr >= 15 & gfr < 30 ~ "G4", gfr < 15 ~ "G5",
      TRUE ~ NA_character_
    )
  } else { gfr_category <- NA_character_ }
  
  if (albumin_available) {
    albumin_category <- case_when(
      albumin_mg_g < 30 ~ "A1", albumin_mg_g >= 30 & albumin_mg_g < 300 ~ "A2",
      albumin_mg_g >= 300 ~ "A3", TRUE ~ NA_character_
    )
  } else { albumin_category <- NA_character_ }
  
  if (gfr_available && albumin_available && within_time_window) {
    risk_recommendation <- get_risk_recommendation(gfr_category, albumin_category)
  } else { risk_recommendation <- list(risk = NA, recommendation = NA) }
  
  return(list(
    gfr_category = gfr_category, gfr_description = get_gfr_description(gfr_category),
    albumin_category = albumin_category, albumin_description = get_albumin_description(albumin_category),
    risk_level = risk_recommendation$risk, recommendation = risk_recommendation$recommendation,
    days_between_tests = if(!is.na(days_apart)) round(days_apart) else NA,
    within_time_window = within_time_window,
    criteria_met_days_from_baseline = if(!is.na(criteria_met_days_from_baseline)) round(criteria_met_days_from_baseline) else NA
  ))
}

get_risk_recommendation <- function(gfr_cat, alb_cat) {
  risk_matrix <- matrix(
    c("Low","Moderate","High","Low","Moderate","High","Moderate","High","Very high",
      "High","Very high","Very high","Very high","Very high","Very high","Very high","Very high","Very high"),
    nrow = 6, ncol = 3, byrow = TRUE,
    dimnames = list(c("G1","G2","G3a","G3b","G4","G5"), c("A1","A2","A3"))
  )
  recommendation_matrix <- matrix(
    c("Screen","Treat","Treat and refer","Screen","Treat","Treat and refer",
      "Treat","Treat","Treat and refer","Treat","Treat and refer","Treat and refer",
      "Treat and refer","Treat and refer","Treat and refer",
      "Treat and refer 4+","Treat and refer 4+","Treat and refer 4+"),
    nrow = 6, ncol = 3, byrow = TRUE,
    dimnames = list(c("G1","G2","G3a","G3b","G4","G5"), c("A1","A2","A3"))
  )
  if (!is.na(gfr_cat) && !is.na(alb_cat)) {
    risk <- risk_matrix[gfr_cat, alb_cat]
    recommendation <- recommendation_matrix[gfr_cat, alb_cat]
  } else { risk <- NA; recommendation <- NA }
  return(list(risk = risk, recommendation = recommendation))
}

get_gfr_description <- function(category) {
  descriptions <- c("G1"="Normal or high","G2"="Mildly decreased","G3a"="Mildly to moderately decreased",
                    "G3b"="Moderately to severely decreased","G4"="Severely decreased","G5"="Kidney failure")
  return(descriptions[category])
}

get_albumin_description <- function(category) {
  descriptions <- c("A1"="Normal to mildly increased (<30 mg/g)",
                    "A2"="Moderately increased (30-299 mg/g)","A3"="Severely increased (≥300 mg/g)")
  return(descriptions[category])
}

# ============================================================
# APPLY CKD CLASSIFICATION
# ============================================================

small_dat <- small_dat %>%
  rowwise() %>%
  mutate(
    ckd_classification = list(classify_ckd(gfr = eGFR_CKD_epi, albumin_mg_g = acr_u,
                                           gfr_days_from_baseline = NULL, albumin_days_from_baseline = NULL))
  ) %>%
  ungroup() %>%
  mutate(
    gfr_category = sapply(ckd_classification, function(x) x$gfr_category),
    gfr_description = sapply(ckd_classification, function(x) x$gfr_description),
    albumin_category = sapply(ckd_classification, function(x) x$albumin_category),
    albumin_description = sapply(ckd_classification, function(x) x$albumin_description),
    risk_level = sapply(ckd_classification, function(x) x$risk_level),
    recommendation = sapply(ckd_classification, function(x) x$recommendation)
  ) %>%
  select(-ckd_classification)

# Check results
cat("\nCKD Classification Results:\n")
small_dat %>%
  select(record_id, sample_type, eGFR_CKD_epi, acr_u, gfr_category, albumin_category, risk_level) %>%
  print(n = 20)

# Summary tables
cat("\nCKD stages by sample type:\n")
print(table(small_dat$sample_type, small_dat$gfr_category, useNA = "always"))

# ============================================================
# CREATE SIDE-BY-SIDE PLOTS BY SERUM/PLASMA
# ============================================================

output_path <- "C:/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/CKD"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
  cat("Created directory:", output_path, "\n")
}

# Biomarker labels
biomarker_labels <- c(
  "ab40_avg_conc" = "Aβ40", "ab42_avg_conc" = "Aβ42", "tau_avg_conc" = "Tau",
  "nfl_avg_conc" = "NFL", "gfap_avg_conc" = "GFAP",
  "ptau_181_avg_conc" = "pTau-181", "ptau_217_avg_conc" = "pTau-217"
)

# Create long format data with sample type
plot_data <- small_dat %>%
  filter(!is.na(sample_type)) %>%  # Only include samples with known sample type
  select(record_id, sample_type, gfr_category, albumin_category, all_of(qx_var)) %>%
  pivot_longer(cols = all_of(qx_var), names_to = "biomarker", values_to = "concentration") %>%
  filter(!is.na(concentration)) %>%
  mutate(
    biomarker_label = biomarker_labels[biomarker],
    sample_type = factor(sample_type, levels = c("serum", "plasma"))
  )

# Check data availability
cat("\nData availability by sample type and GFR category:\n")
plot_data %>%
  filter(!is.na(gfr_category)) %>%
  distinct(record_id, sample_type, gfr_category) %>%
  count(sample_type, gfr_category) %>%
  print()

# ============================================================
# PLOT 1: GFR Category - Side by Side Serum/Plasma
# ============================================================

# Function to create single biomarker plot for GFR
create_gfr_plot <- function(data, biomarker_name) {
  bm_data <- data %>% filter(biomarker == biomarker_name, !is.na(gfr_category))
  
  if(nrow(bm_data) == 0) return(NULL)
  
  ggplot(bm_data, aes(x = gfr_category, y = concentration, fill = sample_type)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, position = position_dodge(width = 0.8)) +
    geom_point(aes(group = sample_type), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
               alpha = 0.4, size = 1) +
    scale_fill_manual(values = c("serum" = "#E69F00", "plasma" = "#56B4E9"),
                      labels = c("serum" = "Serum", "plasma" = "Plasma")) +
    labs(title = biomarker_labels[biomarker_name],
         x = "GFR Category", y = "Concentration", fill = "Sample Type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(face = "bold", size = 11),
          legend.position = "bottom")
}

# Create all GFR plots
gfr_plots <- lapply(qx_var, function(bm) create_gfr_plot(plot_data, bm))
gfr_plots <- gfr_plots[!sapply(gfr_plots, is.null)]  # Remove NULL plots

# Combine plots
if(length(gfr_plots) > 0) {
  p1_combined <- wrap_plots(gfr_plots, ncol = 2) +
    plot_annotation(
      title = "Brain Biomarkers by GFR Category: Serum vs Plasma",
      subtitle = "Orange = Serum (RH2 participants), Blue = Plasma (CRC/PEN participants)",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  print(p1_combined)
  ggsave(file.path(output_path, "biomarkers_by_gfr_serum_plasma.png"), 
         p1_combined, width = 12, height = 14, dpi = 300)
  cat("\nSaved: biomarkers_by_gfr_serum_plasma.png\n")
}

# ============================================================
# PLOT 2: Albumin Category - Side by Side Serum/Plasma
# ============================================================

create_albumin_plot <- function(data, biomarker_name) {
  bm_data <- data %>% filter(biomarker == biomarker_name, !is.na(albumin_category))
  
  if(nrow(bm_data) == 0) return(NULL)
  
  ggplot(bm_data, aes(x = albumin_category, y = concentration, fill = sample_type)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, position = position_dodge(width = 0.8)) +
    geom_point(aes(group = sample_type), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
               alpha = 0.4, size = 1) +
    scale_fill_manual(values = c("serum" = "#E69F00", "plasma" = "#56B4E9"),
                      labels = c("serum" = "Serum", "plasma" = "Plasma")) +
    labs(title = biomarker_labels[biomarker_name],
         x = "Albumin Category", y = "Concentration", fill = "Sample Type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(face = "bold", size = 11),
          legend.position = "bottom")
}

# Create all Albumin plots
alb_plots <- lapply(qx_var, function(bm) create_albumin_plot(plot_data, bm))
alb_plots <- alb_plots[!sapply(alb_plots, is.null)]

if(length(alb_plots) > 0) {
  p2_combined <- wrap_plots(alb_plots, ncol = 2) +
    plot_annotation(
      title = "Brain Biomarkers by Albumin Category: Serum vs Plasma",
      subtitle = "Orange = Serum (RH2 participants), Blue = Plasma (CRC/PEN participants)",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  print(p2_combined)
  ggsave(file.path(output_path, "biomarkers_by_albumin_serum_plasma.png"), 
         p2_combined, width = 12, height = 14, dpi = 300)
  cat("\nSaved: biomarkers_by_albumin_serum_plasma.png\n")
}

# ============================================================
# ALTERNATIVE: Faceted by Sample Type (Separate Panels)
# ============================================================

# GFR - Faceted version
p3_faceted <- ggplot(plot_data %>% filter(!is.na(gfr_category)), 
                     aes(x = gfr_category, y = concentration, fill = gfr_category)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_grid(biomarker_label ~ sample_type, scales = "free_y",
             labeller = labeller(sample_type = c("serum" = "Serum", "plasma" = "Plasma"))) +
  labs(title = "Brain Biomarkers by GFR Category",
       subtitle = "Split by Sample Type (Serum vs Plasma)",
       x = "GFR Category", y = "Concentration", fill = "GFR Category") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

print(p3_faceted)
ggsave(file.path(output_path, "biomarkers_by_gfr_faceted.png"), 
       p3_faceted, width = 10, height = 16, dpi = 300)

# Albumin - Faceted version
p4_faceted <- ggplot(plot_data %>% filter(!is.na(albumin_category)), 
                     aes(x = albumin_category, y = concentration, fill = albumin_category)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_grid(biomarker_label ~ sample_type, scales = "free_y",
             labeller = labeller(sample_type = c("serum" = "Serum", "plasma" = "Plasma"))) +
  labs(title = "Brain Biomarkers by Albumin Category",
       subtitle = "Split by Sample Type (Serum vs Plasma)",
       x = "Albumin Category", y = "Concentration", fill = "Albumin Category") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

print(p4_faceted)
ggsave(file.path(output_path, "biomarkers_by_albumin_faceted.png"), 
       p4_faceted, width = 10, height = 16, dpi = 300)

# ============================================================
# SUMMARY STATISTICS BY SAMPLE TYPE AND CKD CATEGORY
# ============================================================

cat("\n\n========== SUMMARY STATISTICS ==========\n")

# Sample sizes
cat("\nSample sizes by GFR category and sample type:\n")
plot_data %>%
  filter(!is.na(gfr_category)) %>%
  distinct(record_id, sample_type, gfr_category) %>%
  count(sample_type, gfr_category) %>%
  pivot_wider(names_from = sample_type, values_from = n, values_fill = 0) %>%
  print()

cat("\nSample sizes by Albumin category and sample type:\n")
plot_data %>%
  filter(!is.na(albumin_category)) %>%
  distinct(record_id, sample_type, albumin_category) %>%
  count(sample_type, albumin_category) %>%
  pivot_wider(names_from = sample_type, values_from = n, values_fill = 0) %>%
  print()

# Mean concentrations by biomarker, sample type, and GFR category
cat("\nMean biomarker concentrations by GFR category and sample type:\n")
summary_stats <- plot_data %>%
  filter(!is.na(gfr_category)) %>%
  group_by(biomarker_label, sample_type, gfr_category) %>%
  summarise(
    n = n(),
    mean = mean(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats, n = 50)

# Save summary stats
write.csv(summary_stats, file.path(output_path, "biomarker_summary_by_ckd_sample_type.csv"), row.names = FALSE)
cat("\nSaved summary statistics to: biomarker_summary_by_ckd_sample_type.csv\n")

cat("\n========== ANALYSIS COMPLETE ==========\n")




























##########################DKD


### Brain Imaging by Healthy vs DKD Classification
### Healthy: eGFR 90-120 AND UACR < 30
### DKD: eGFR < 90 or > 120, AND/OR UACR > 300

# Load libraries
library(oro.nifti)  
library(neurobase)  
library(R.matlab)    
library(ggplot2)     
library(dplyr)       
library(corrplot)   
library(psych)        
library(pls)  
library(caret) 
library(randomForest) 
library(igraph)
library(brainGraph)
library(R.matlab)
library(tidyverse)
library(patchwork)

bucket <- 'brain.mri'

# ============================================================
# DATA LOADING AND PREPARATION
# ============================================================

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

# ============================================================
# DETERMINE SAMPLE TYPE (SERUM/PLASMA) BY STUDY
# ============================================================

dat <- dat %>%
  mutate(
    sample_type = case_when(
      grepl("^CRC", record_id, ignore.case = TRUE) ~ "plasma",
      grepl("^PEN", record_id, ignore.case = TRUE) ~ "plasma",
      study == "RH2" ~ "serum",
      grepl("^RH2", record_id, ignore.case = TRUE) ~ "serum",
      study == "RENAL-HEIR" ~ "serum",
      study == "RENAL-HEIRitage" ~ "serum",
      study == "CROCODILE" ~ "plasma",
      study == "PENGUIN" ~ "plasma",
      TRUE ~ NA_character_
    )
  )

cat("\nSample type distribution:\n")
print(table(dat$sample_type, useNA = "always"))

small_dat <- dat

# Fix the RH2-38-T record
small_dat$group[which(small_dat$record_id == 'RH2-38-T')] <- 'Obese Control'
small_dat$record_id[which(small_dat$record_id == 'RH2-38-T')] <- 'RH2-38-O'

# Define biomarkers
qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")

# ============================================================
# NEW CLASSIFICATION: HEALTHY vs DKD
# ============================================================
# Healthy: eGFR 90-120 AND UACR < 30
# DKD: eGFR < 90 or > 120, AND/OR UACR > 300

small_dat <- small_dat %>%
  mutate(
    kidney_status = case_when(
      # Healthy: eGFR between 90-120 AND UACR < 30
      (eGFR_CKD_epi >= 90 & eGFR_CKD_epi <= 120) & (acr_u < 30) ~ "Healthy",
      # DKD: eGFR outside 90-120 range OR UACR > 300
      (eGFR_CKD_epi < 90 | eGFR_CKD_epi > 120) | (acr_u > 300) ~ "DKD",
      # Everyone else (e.g., normal eGFR but UACR 30-300) - intermediate/unclassified
      TRUE ~ NA_character_
    )
  )

# Check classification results
cat("\n========== KIDNEY STATUS CLASSIFICATION ==========\n")
cat("\nOverall distribution:\n")
print(table(small_dat$kidney_status, useNA = "always"))

cat("\nBy sample type:\n")
print(table(small_dat$sample_type, small_dat$kidney_status, useNA = "always"))

# Summary of eGFR and UACR by group
cat("\nSummary statistics by kidney status:\n")
small_dat %>%
  filter(!is.na(kidney_status)) %>%
  group_by(kidney_status) %>%
  summarise(
    n = n(),
    eGFR_mean = mean(eGFR_CKD_epi, na.rm = TRUE),
    eGFR_sd = sd(eGFR_CKD_epi, na.rm = TRUE),
    eGFR_min = min(eGFR_CKD_epi, na.rm = TRUE),
    eGFR_max = max(eGFR_CKD_epi, na.rm = TRUE),
    UACR_mean = mean(acr_u, na.rm = TRUE),
    UACR_sd = sd(acr_u, na.rm = TRUE),
    UACR_median = median(acr_u, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# Show which participants fall into each category
cat("\nParticipants by kidney status:\n")
participant_summary <- small_dat %>%
  filter(!is.na(kidney_status)) %>%
  select(record_id, sample_type, eGFR_CKD_epi, acr_u, kidney_status) %>%
  arrange(kidney_status, record_id) %>%
  as.data.frame()
print(participant_summary)

# ============================================================
# OUTPUT PATH
# ============================================================

output_path <- "C:/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/CKD"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
  cat("Created directory:", output_path, "\n")
}

# Biomarker labels
biomarker_labels <- c(
  "ab40_avg_conc" = "Aβ40", "ab42_avg_conc" = "Aβ42", "tau_avg_conc" = "Tau",
  "nfl_avg_conc" = "NFL", "gfap_avg_conc" = "GFAP",
  "ptau_181_avg_conc" = "pTau-181", "ptau_217_avg_conc" = "pTau-217"
)

# Create long format data
plot_data <- small_dat %>%
  filter(!is.na(sample_type), !is.na(kidney_status)) %>%
  select(record_id, sample_type, kidney_status, eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
  pivot_longer(cols = all_of(qx_var), names_to = "biomarker", values_to = "concentration") %>%
  filter(!is.na(concentration)) %>%
  mutate(
    biomarker_label = biomarker_labels[biomarker],
    sample_type = factor(sample_type, levels = c("serum", "plasma")),
    kidney_status = factor(kidney_status, levels = c("Healthy", "DKD"))
  )

# Check data availability
cat("\nData availability by sample type and kidney status:\n")
plot_data %>%
  distinct(record_id, sample_type, kidney_status) %>%
  count(sample_type, kidney_status) %>%
  print()

# ============================================================
# PLOT 1: Healthy vs DKD - Side by Side Serum/Plasma
# ============================================================

create_kidney_plot <- function(data, biomarker_name) {
  bm_data <- data %>% filter(biomarker == biomarker_name)
  
  if(nrow(bm_data) == 0) return(NULL)
  
  ggplot(bm_data, aes(x = kidney_status, y = concentration, fill = sample_type)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, position = position_dodge(width = 0.8)) +
    geom_point(aes(group = sample_type), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
               alpha = 0.4, size = 1.5) +
    scale_fill_manual(values = c("serum" = "#E69F00", "plasma" = "#56B4E9"),
                      labels = c("serum" = "Serum", "plasma" = "Plasma")) +
    labs(title = biomarker_labels[biomarker_name],
         x = "Kidney Status", y = "Concentration (pg/mL)", fill = "Sample Type") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 11),
          legend.position = "bottom")
}

# Create all plots
kidney_plots <- lapply(qx_var, function(bm) create_kidney_plot(plot_data, bm))
kidney_plots <- kidney_plots[!sapply(kidney_plots, is.null)]

if(length(kidney_plots) > 0) {
  p1_combined <- wrap_plots(kidney_plots, ncol = 2) +
    plot_annotation(
      title = "Brain Biomarkers: Healthy vs DKD",
      subtitle = "Healthy: eGFR 90-120 & UACR <30 | DKD: eGFR <90 or >120, and/or UACR >300\nOrange = Serum, Blue = Plasma",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  print(p1_combined)
  ggsave(file.path(output_path, "biomarkers_healthy_vs_dkd_serum_plasma.png"), 
         p1_combined, width = 12, height = 14, dpi = 300)
  cat("\nSaved: biomarkers_healthy_vs_dkd_serum_plasma.png\n")
}

# ============================================================
# PLOT 2: Faceted by Sample Type
# ============================================================

p2_faceted <- ggplot(plot_data, 
                     aes(x = kidney_status, y = concentration, fill = kidney_status)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  facet_grid(biomarker_label ~ sample_type, scales = "free_y",
             labeller = labeller(sample_type = c("serum" = "Serum", "plasma" = "Plasma"))) +
  scale_fill_manual(values = c("Healthy" = "#2E86AB", "DKD" = "#E94F37")) +
  labs(title = "Brain Biomarkers: Healthy vs DKD",
       subtitle = "Healthy: eGFR 90-120 & UACR <30 | DKD: eGFR <90 or >120, and/or UACR >300",
       x = "Kidney Status", y = "Concentration (pg/mL)", fill = "Status") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")

print(p2_faceted)
ggsave(file.path(output_path, "biomarkers_healthy_vs_dkd_faceted.png"), 
       p2_faceted, width = 10, height = 16, dpi = 300)

# ============================================================
# PLOT 3: Combined (ignore sample type)
# ============================================================

p3_combined_all <- ggplot(plot_data, 
                          aes(x = kidney_status, y = concentration, fill = kidney_status)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  facet_wrap(~ biomarker_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Healthy" = "#2E86AB", "DKD" = "#E94F37")) +
  labs(title = "Brain Biomarkers: Healthy vs DKD (All Samples Combined)",
       subtitle = "Healthy: eGFR 90-120 & UACR <30 | DKD: eGFR <90 or >120, and/or UACR >300",
       x = "Kidney Status", y = "Concentration (pg/mL)", fill = "Status") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        axis.text.x = element_text(size = 10))

print(p3_combined_all)
ggsave(file.path(output_path, "biomarkers_healthy_vs_dkd_combined.png"), 
       p3_combined_all, width = 10, height = 12, dpi = 300)

# ============================================================
# STATISTICAL TESTS: Healthy vs DKD
# ============================================================

cat("\n\n========== STATISTICAL COMPARISONS ==========\n")

# Wilcoxon tests for each biomarker (overall)
cat("\nWilcoxon rank-sum tests (Healthy vs DKD - all samples):\n")
stat_results_all <- plot_data %>%
  group_by(biomarker_label) %>%
  summarise(
    n_healthy = sum(kidney_status == "Healthy"),
    n_dkd = sum(kidney_status == "DKD"),
    median_healthy = median(concentration[kidney_status == "Healthy"], na.rm = TRUE),
    median_dkd = median(concentration[kidney_status == "DKD"], na.rm = TRUE),
    p_value = wilcox.test(concentration ~ kidney_status)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         significant = ifelse(p_adj < 0.05, "*", ""))

print(stat_results_all)

# By sample type
cat("\nWilcoxon tests by sample type:\n")
stat_results_by_sample <- plot_data %>%
  group_by(sample_type, biomarker_label) %>%
  summarise(
    n_healthy = sum(kidney_status == "Healthy"),
    n_dkd = sum(kidney_status == "DKD"),
    median_healthy = median(concentration[kidney_status == "Healthy"], na.rm = TRUE),
    median_dkd = median(concentration[kidney_status == "DKD"], na.rm = TRUE),
    p_value = tryCatch(wilcox.test(concentration ~ kidney_status)$p.value, error = function(e) NA),
    .groups = "drop"
  )

print(stat_results_by_sample, n = 20)

# ============================================================
# SUMMARY STATISTICS
# ============================================================

cat("\n\n========== SUMMARY STATISTICS ==========\n")

# Sample sizes
cat("\nSample sizes by kidney status and sample type:\n")
plot_data %>%
  distinct(record_id, sample_type, kidney_status) %>%
  count(sample_type, kidney_status) %>%
  pivot_wider(names_from = sample_type, values_from = n, values_fill = 0) %>%
  print()

# Mean concentrations
cat("\nMean biomarker concentrations by kidney status:\n")
summary_stats <- plot_data %>%
  group_by(biomarker_label, kidney_status) %>%
  summarise(
    n = n(),
    mean = mean(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    IQR = IQR(concentration, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# Save results
write.csv(stat_results_all, file.path(output_path, "stats_healthy_vs_dkd.csv"), row.names = FALSE)
write.csv(summary_stats, file.path(output_path, "summary_healthy_vs_dkd.csv"), row.names = FALSE)

cat("\n========== ANALYSIS COMPLETE ==========\n")























# ============================================================
# ADDITIONAL ANALYSES: GROUP-ADJUSTED PLOTS AND CORRELATIONS
# ============================================================

# ============================================================
# PLOT 4: Healthy vs DKD - Adjusted by Group (T1D/T2D/Control)
# ============================================================

cat("\n\n========== GROUP-ADJUSTED ANALYSES ==========\n")

# Check group distribution
cat("\nGroup distribution by kidney status:\n")
small_dat %>%
  filter(!is.na(kidney_status)) %>%
  count(kidney_status, group) %>%
  pivot_wider(names_from = kidney_status, values_from = n, values_fill = 0) %>%
  print()

# Create plot data with group
plot_data_group <- small_dat %>%
  filter(!is.na(sample_type), !is.na(kidney_status), !is.na(group)) %>%
  select(record_id, sample_type, kidney_status, group, eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
  pivot_longer(cols = all_of(qx_var), names_to = "biomarker", values_to = "concentration") %>%
  filter(!is.na(concentration)) %>%
  mutate(
    biomarker_label = biomarker_labels[biomarker],
    sample_type = factor(sample_type, levels = c("serum", "plasma")),
    kidney_status = factor(kidney_status, levels = c("Healthy", "DKD")),
    group = factor(group)
  )

# Plot: Faceted by biomarker and group
p4_group_adjusted <- ggplot(plot_data_group, 
                            aes(x = kidney_status, y = concentration, fill = kidney_status)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  facet_grid(biomarker_label ~ group, scales = "free_y") +
  scale_fill_manual(values = c("Healthy" = "#2E86AB", "DKD" = "#E94F37")) +
  labs(title = "Brain Biomarkers: Healthy vs DKD by Group",
       subtitle = "Stratified by T1D, T2D, and Control groups",
       x = "Kidney Status", y = "Concentration (pg/mL)", fill = "Status") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

print(p4_group_adjusted)
ggsave(file.path(output_path, "biomarkers_healthy_vs_dkd_by_group.png"), 
       p4_group_adjusted, width = 14, height = 16, dpi = 300)

# Alternative: Each biomarker separately with group coloring
create_group_biomarker_plot <- function(data, biomarker_name) {
  bm_data <- data %>% filter(biomarker == biomarker_name)
  
  if(nrow(bm_data) == 0) return(NULL)
  
  ggplot(bm_data, aes(x = kidney_status, y = concentration, fill = group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, position = position_dodge(width = 0.8)) +
    geom_point(aes(group = group), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
               alpha = 0.5, size = 1.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = biomarker_labels[biomarker_name],
         x = "Kidney Status", y = "Concentration (pg/mL)", fill = "Group") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "bottom")
}

group_plots <- lapply(qx_var, function(bm) create_group_biomarker_plot(plot_data_group, bm))
group_plots <- group_plots[!sapply(group_plots, is.null)]

if(length(group_plots) > 0) {
  p5_group_combined <- wrap_plots(group_plots, ncol = 2) +
    plot_annotation(
      title = "Brain Biomarkers by Group: Healthy vs DKD",
      subtitle = "Each panel shows distribution across T1D, T2D, and Control groups",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  print(p5_group_combined)
  ggsave(file.path(output_path, "biomarkers_healthy_vs_dkd_group_colored.png"), 
         p5_group_combined, width = 12, height = 14, dpi = 300)
}

# ============================================================
# PLOT 5: UACR-ADJUSTED PLOTS (Continuous UACR)
# ============================================================

cat("\n\n========== UACR-ADJUSTED ANALYSES ==========\n")

# Create correlation data
uacr_data <- small_dat %>%
  filter(!is.na(sample_type), !is.na(acr_u)) %>%
  select(record_id, sample_type, group, kidney_status, eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
  pivot_longer(cols = all_of(qx_var), names_to = "biomarker", values_to = "concentration") %>%
  filter(!is.na(concentration)) %>%
  mutate(
    biomarker_label = biomarker_labels[biomarker],
    log_uacr = log10(acr_u + 1),  # Log transform UACR
    sample_type = factor(sample_type, levels = c("serum", "plasma"))
  )

# Scatter plots: Biomarker vs log(UACR)
p6_uacr_scatter <- ggplot(uacr_data, 
                          aes(x = log_uacr, y = concentration, color = kidney_status)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  facet_wrap(~ biomarker_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Healthy" = "#2E86AB", "DKD" = "#E94F37"), 
                     na.value = "gray50") +
  labs(title = "Brain Biomarkers vs Log(UACR)",
       subtitle = "Linear relationship with albumin-to-creatinine ratio",
       x = "Log10(UACR + 1) [mg/g]", y = "Concentration (pg/mL)", 
       color = "Kidney Status") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom")

print(p6_uacr_scatter)
ggsave(file.path(output_path, "biomarkers_vs_uacr_scatter.png"), 
       p6_uacr_scatter, width = 12, height = 14, dpi = 300)

# Facet by sample type as well
p7_uacr_by_sample <- ggplot(uacr_data, 
                            aes(x = log_uacr, y = concentration, color = sample_type)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  facet_wrap(~ biomarker_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("serum" = "#E69F00", "plasma" = "#56B4E9")) +
  labs(title = "Brain Biomarkers vs Log(UACR) by Sample Type",
       x = "Log10(UACR + 1) [mg/g]", y = "Concentration (pg/mL)", 
       color = "Sample Type") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom")

print(p7_uacr_by_sample)
ggsave(file.path(output_path, "biomarkers_vs_uacr_by_sample.png"), 
       p7_uacr_by_sample, width = 12, height = 14, dpi = 300)

# ============================================================
# PLOT 6: eGFR-ADJUSTED PLOTS (Continuous eGFR)
# ============================================================

egfr_data <- small_dat %>%
  filter(!is.na(sample_type), !is.na(eGFR_CKD_epi)) %>%
  select(record_id, sample_type, group, kidney_status, eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
  pivot_longer(cols = all_of(qx_var), names_to = "biomarker", values_to = "concentration") %>%
  filter(!is.na(concentration)) %>%
  mutate(
    biomarker_label = biomarker_labels[biomarker],
    sample_type = factor(sample_type, levels = c("serum", "plasma"))
  )

p8_egfr_scatter <- ggplot(egfr_data, 
                          aes(x = eGFR_CKD_epi, y = concentration, color = kidney_status)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  geom_vline(xintercept = c(90, 120), linetype = "dashed", color = "gray40", alpha = 0.6) +
  facet_wrap(~ biomarker_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Healthy" = "#2E86AB", "DKD" = "#E94F37"), 
                     na.value = "gray50") +
  labs(title = "Brain Biomarkers vs eGFR",
       subtitle = "Dashed lines indicate healthy range (90-120 mL/min/1.73m²)",
       x = "eGFR (mL/min/1.73m²)", y = "Concentration (pg/mL)", 
       color = "Kidney Status") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom")

print(p8_egfr_scatter)
ggsave(file.path(output_path, "biomarkers_vs_egfr_scatter.png"), 
       p8_egfr_scatter, width = 12, height = 14, dpi = 300)

# ============================================================
# CORRELATION ANALYSES
# ============================================================

cat("\n\n========== CORRELATION ANALYSES ==========\n")

# 1. Biomarker intercorrelations
cat("\n1. BIOMARKER INTERCORRELATIONS\n")

correlation_data <- small_dat %>%
  filter(!is.na(kidney_status)) %>%
  dplyr::select(all_of(qx_var)) %>%
  drop_na()

if(nrow(correlation_data) > 0 && ncol(correlation_data) > 1) {
  cor_matrix <- cor(correlation_data, use = "pairwise.complete.obs", method = "spearman")
  
  # Rename columns for better display
  cor_matrix_display <- cor_matrix
  rownames(cor_matrix_display) <- biomarker_labels[rownames(cor_matrix_display)]
  colnames(cor_matrix_display) <- biomarker_labels[colnames(cor_matrix_display)]
  
  # Open PNG device
  png(file.path(output_path, "biomarker_intercorrelations.png"), 
      width = 10, height = 10, units = "in", res = 300)
  
  # Create corrplot
  corrplot::corrplot(cor_matrix_display, 
                     method = "color", 
                     type = "upper", 
                     tl.col = "black", 
                     tl.srt = 45,
                     tl.cex = 0.9,
                     addCoef.col = "black", 
                     number.cex = 0.7,
                     col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
                     diag = TRUE,
                     cl.cex = 0.8,
                     mar = c(0, 0, 2, 0),
                     title = "Biomarker Intercorrelations (All Samples)")
  
  # Close device
  dev.off()
  
  cat("\nSpearman correlations between biomarkers:\n")
  print(round(cor_matrix, 3))
  
  # Save correlation matrix as CSV
  write.csv(cor_matrix, 
            file.path(output_path, "biomarker_intercorrelations.csv"))
  
} else {
  cat("Insufficient data for biomarker intercorrelations.\n")
}

# 2. Correlations with kidney function (eGFR and UACR)
cat("\n\n2. CORRELATIONS WITH KIDNEY FUNCTION\n")

kidney_cor_data <- small_dat %>%
  filter(!is.na(kidney_status)) %>%
  dplyr::select(eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
  mutate(log_uacr = log10(acr_u + 1))

# Calculate correlations
kidney_correlations <- data.frame(
  biomarker = qx_var,
  biomarker_label = biomarker_labels[qx_var]
)

for(bm in qx_var) {
  # eGFR correlation
  if(sum(!is.na(kidney_cor_data[[bm]]) & !is.na(kidney_cor_data$eGFR_CKD_epi)) > 3) {
    egfr_test <- cor.test(kidney_cor_data[[bm]], kidney_cor_data$eGFR_CKD_epi, 
                          method = "spearman", exact = FALSE)
    kidney_correlations[kidney_correlations$biomarker == bm, "egfr_rho"] <- egfr_test$estimate
    kidney_correlations[kidney_correlations$biomarker == bm, "egfr_p"] <- egfr_test$p.value
  } else {
    kidney_correlations[kidney_correlations$biomarker == bm, "egfr_rho"] <- NA
    kidney_correlations[kidney_correlations$biomarker == bm, "egfr_p"] <- NA
  }
  
  # UACR correlation
  if(sum(!is.na(kidney_cor_data[[bm]]) & !is.na(kidney_cor_data$log_uacr)) > 3) {
    uacr_test <- cor.test(kidney_cor_data[[bm]], kidney_cor_data$log_uacr, 
                          method = "spearman", exact = FALSE)
    kidney_correlations[kidney_correlations$biomarker == bm, "log_uacr_rho"] <- uacr_test$estimate
    kidney_correlations[kidney_correlations$biomarker == bm, "log_uacr_p"] <- uacr_test$p.value
  } else {
    kidney_correlations[kidney_correlations$biomarker == bm, "log_uacr_rho"] <- NA
    kidney_correlations[kidney_correlations$biomarker == bm, "log_uacr_p"] <- NA
  }
}

kidney_correlations <- kidney_correlations %>%
  mutate(
    egfr_sig = ifelse(!is.na(egfr_p) & egfr_p < 0.05, "*", ""),
    log_uacr_sig = ifelse(!is.na(log_uacr_p) & log_uacr_p < 0.05, "*", "")
  )

cat("\nCorrelations with kidney function:\n")
print(kidney_correlations %>% dplyr::select(-biomarker), row.names = FALSE)

write.csv(kidney_correlations, 
          file.path(output_path, "kidney_function_correlations.csv"), 
          row.names = FALSE)

# 3. Stratified correlations by kidney status
cat("\n\n3. STRATIFIED CORRELATIONS BY KIDNEY STATUS\n")

for(status in c("Healthy", "DKD")) {
  cat(paste0("\n", status, " group:\n"))
  
  status_data <- small_dat %>%
    filter(kidney_status == status) %>%
    dplyr::select(eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
    mutate(log_uacr = log10(acr_u + 1)) %>%
    dplyr::select(-acr_u) %>%
    drop_na()
  
  if(nrow(status_data) > 5 && ncol(status_data) > 2) {
    status_cor <- cor(status_data, use = "pairwise.complete.obs", method = "spearman")
    
    # Rename for display
    status_cor_display <- status_cor
    colnames(status_cor_display)[1] <- "eGFR"
    rownames(status_cor_display)[1] <- "eGFR"
    colnames(status_cor_display)[2] <- "Log10(UACR)"
    rownames(status_cor_display)[2] <- "Log10(UACR)"
    for(i in 3:ncol(status_cor_display)) {
      old_name <- colnames(status_cor_display)[i]
      if(old_name %in% names(biomarker_labels)) {
        colnames(status_cor_display)[i] <- biomarker_labels[old_name]
        rownames(status_cor_display)[i] <- biomarker_labels[old_name]
      }
    }
    
    png(file.path(output_path, paste0("correlations_", tolower(status), ".png")), 
        width = 12, height = 10, units = "in", res = 300)
    
    corrplot::corrplot(status_cor_display, 
                       method = "color", 
                       type = "upper", 
                       tl.col = "black", 
                       tl.srt = 45,
                       tl.cex = 0.9,
                       addCoef.col = "black", 
                       number.cex = 0.7,
                       col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
                       diag = TRUE,
                       cl.cex = 0.8,
                       mar = c(0, 0, 2, 0),
                       title = paste("Correlations -", status, "Group"))
    
    dev.off()
    
    cat("Correlation matrix saved.\n")
    # Print eGFR and log(UACR) correlations with biomarkers
    if(ncol(status_cor) > 2) {
      cat("eGFR and Log10(UACR) correlations with biomarkers:\n")
      print(round(status_cor[1:2, 3:ncol(status_cor)], 3))
    }
  } else {
    cat("Insufficient data for correlation analysis.\n")
  }
}

# 4. Correlations by group (T1D/T2D/Control)
cat("\n\n4. CORRELATIONS BY DIABETES GROUP\n")

for(grp in unique(small_dat$group[!is.na(small_dat$group)])) {
  cat(paste0("\n", grp, ":\n"))
  
  group_data <- small_dat %>%
    filter(group == grp, !is.na(kidney_status)) %>%
    dplyr::select(eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
    mutate(log_uacr = log10(acr_u + 1)) %>%
    dplyr::select(-acr_u) %>%
    drop_na()
  
  if(nrow(group_data) > 5 && ncol(group_data) > 2) {
    group_cor <- cor(group_data, use = "pairwise.complete.obs", method = "spearman")
    
    # Rename for display
    group_cor_display <- group_cor
    colnames(group_cor_display)[1] <- "eGFR"
    rownames(group_cor_display)[1] <- "eGFR"
    colnames(group_cor_display)[2] <- "Log10(UACR)"
    rownames(group_cor_display)[2] <- "Log10(UACR)"
    for(i in 3:ncol(group_cor_display)) {
      old_name <- colnames(group_cor_display)[i]
      if(old_name %in% names(biomarker_labels)) {
        colnames(group_cor_display)[i] <- biomarker_labels[old_name]
        rownames(group_cor_display)[i] <- biomarker_labels[old_name]
      }
    }
    
    safe_grp_name <- gsub("[^A-Za-z0-9]", "_", grp)
    png(file.path(output_path, paste0("correlations_", tolower(safe_grp_name), ".png")), 
        width = 12, height = 10, units = "in", res = 300)
    
    corrplot::corrplot(group_cor_display, 
                       method = "color", 
                       type = "upper", 
                       tl.col = "black", 
                       tl.srt = 45,
                       tl.cex = 0.9,
                       addCoef.col = "black", 
                       number.cex = 0.7,
                       col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
                       diag = TRUE,
                       cl.cex = 0.8,
                       mar = c(0, 0, 2, 0),
                       title = paste("Correlations -", grp))
    
    dev.off()
    
    cat("Correlation matrix saved.\n")
  } else {
    cat("Insufficient data.\n")
  }
}

# 5. Partial correlations (adjusting for age, BMI, etc.)
cat("\n\n5. PARTIAL CORRELATIONS (if demographic data available)\n")

if(all(c("age", "bmi") %in% names(small_dat))) {
  # Check if ppcor is installed
  if(!require(ppcor, quietly = TRUE)) {
    cat("Installing ppcor package...\n")
    install.packages("ppcor")
    library(ppcor)
  }
  
  partial_cor_data <- small_dat %>%
    filter(!is.na(kidney_status)) %>%
    dplyr::select(age, bmi, eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
    mutate(log_uacr = log10(acr_u + 1)) %>%
    dplyr::select(-acr_u) %>%
    drop_na()
  
  if(nrow(partial_cor_data) > 10) {
    # Partial correlation controlling for age and BMI
    pcor_results <- pcor(partial_cor_data, method = "spearman")
    
    cat("\nPartial correlations (controlling for age, BMI):\n")
    cat("eGFR and Log10(UACR) with biomarkers:\n")
    print(round(pcor_results$estimate[3:4, 5:ncol(pcor_results$estimate)], 3))
    
    cat("\nP-values:\n")
    print(round(pcor_results$p.value[3:4, 5:ncol(pcor_results$p.value)], 4))
    
    # Save partial correlation results
    partial_results <- data.frame(
      biomarker = colnames(pcor_results$estimate)[5:ncol(pcor_results$estimate)],
      egfr_partial_r = pcor_results$estimate[3, 5:ncol(pcor_results$estimate)],
      egfr_partial_p = pcor_results$p.value[3, 5:ncol(pcor_results$p.value)],
      log_uacr_partial_r = pcor_results$estimate[4, 5:ncol(pcor_results$estimate)],
      log_uacr_partial_p = pcor_results$p.value[4, 5:ncol(pcor_results$p.value)]
    )
    
    write.csv(partial_results, 
              file.path(output_path, "partial_correlations_age_bmi_adjusted.csv"), 
              row.names = FALSE)
    cat("\nPartial correlation results saved.\n")
  } else {
    cat("Insufficient data for partial correlation analysis.\n")
  }
} else {
  cat("Age or BMI not available for adjustment.\n")
}

# 6. Create a comprehensive correlation heatmap
cat("\n\n6. COMPREHENSIVE CORRELATION HEATMAP\n")

comprehensive_data <- small_dat %>%
  filter(!is.na(kidney_status)) %>%
  dplyr::select(eGFR_CKD_epi, acr_u, all_of(qx_var)) %>%
  mutate(log_uacr = log10(acr_u + 1)) %>%
  dplyr::select(-acr_u) %>%
  drop_na()

if(nrow(comprehensive_data) > 5 && ncol(comprehensive_data) > 2) {
  comp_cor <- cor(comprehensive_data, use = "pairwise.complete.obs", method = "spearman")
  
  # Rename for better labels
  comp_cor_display <- comp_cor
  colnames(comp_cor_display)[1] <- "eGFR"
  rownames(comp_cor_display)[1] <- "eGFR"
  colnames(comp_cor_display)[2] <- "Log10(UACR)"
  rownames(comp_cor_display)[2] <- "Log10(UACR)"
  for(i in 3:ncol(comp_cor_display)) {
    old_name <- colnames(comp_cor_display)[i]
    if(old_name %in% names(biomarker_labels)) {
      colnames(comp_cor_display)[i] <- biomarker_labels[old_name]
      rownames(comp_cor_display)[i] <- biomarker_labels[old_name]
    }
  }
  
  png(file.path(output_path, "comprehensive_correlations_all.png"), 
      width = 12, height = 11, units = "in", res = 300)
  
  corrplot::corrplot(comp_cor_display, 
                     method = "color", 
                     type = "upper", 
                     tl.col = "black", 
                     tl.srt = 45,
                     tl.cex = 1,
                     addCoef.col = "black", 
                     number.cex = 0.7,
                     col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
                     diag = TRUE,
                     cl.cex = 0.8,
                     mar = c(0, 0, 3, 0),
                     title = "Comprehensive Correlations: Kidney Function & Brain Biomarkers")
  
  dev.off()
  
  cat("Comprehensive correlation heatmap saved.\n")
  
  # Print the full correlation matrix
  cat("\nFull correlation matrix:\n")
  print(round(comp_cor, 3))
  
  # Save as CSV
  write.csv(comp_cor, 
            file.path(output_path, "comprehensive_correlations.csv"))
}

# ============================================================
# SUMMARY TABLE: Group-stratified statistics
# ============================================================

cat("\n\n========== GROUP-STRATIFIED SUMMARY ==========\n")

group_summary <- plot_data_group %>%
  group_by(biomarker_label, group, kidney_status) %>%
  summarise(
    n = n(),
    mean = mean(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    IQR = IQR(concentration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(biomarker_label, group, kidney_status)

print(group_summary, n = 50)

write.csv(group_summary, 
          file.path(output_path, "summary_by_group_and_kidney_status.csv"), 
          row.names = FALSE)

cat("\n========== ALL ANALYSES COMPLETE ==========\n")
cat(paste0("All plots and tables saved to: ", output_path, "\n"))



