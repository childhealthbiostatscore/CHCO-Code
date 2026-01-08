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

