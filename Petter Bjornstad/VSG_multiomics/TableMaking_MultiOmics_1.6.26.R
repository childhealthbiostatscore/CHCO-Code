##### Making tables for multiomics VSG paper 


library(readr)
library(dplyr)
library(gridExtra)
library(Seurat)
library(lme4)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(stats) 
#
library(lmerTest)
library(broom.mixed)
library(emmeans)
#
library(UpSetR)
library(tidyr)
#
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggridges)
library(tidytext)
library(msigdbr)
library(fgsea)
library(stringr)
library(purrr)
library(rlang)
#
library(KEGGREST)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)
library(edgeR)
library(limma)
library(ggrepel)
#

base_dir <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/scrna/"

scrna <- readRDS(paste0(base_dir, "data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds"))
scrna <- subset(scrna, subset = percent.mt < 50 & nFeature_RNA < 5000 & nFeature_RNA > 500) 

#
clin <- read_csv(paste0(base_dir, "data_clean/pb90_rhrh2improve_clinical_subset.csv"), show_col_types = FALSE) %>%
  dplyr::select(record_id, visit,acr_u, gfr_bsa_plasma, gfr_raw_plasma, bmi, eGFR_CKiD_U25_avg, eGFR_fas_cr_cysc) %>%
  dplyr::mutate(visit = dplyr::recode(visit,"baseline" = "pre","12_months_post_surgery" = "post"))


#############
## IMPROVE ##
#############
improve <- subset(
  scrna, 
  subset = cohort == "IMPROVE" | 
    record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T"))
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$pre_post <- 'pre'
improve@meta.data <- improve@meta.data %>%
  tibble::rownames_to_column("cell") %>%   
  left_join(clin, by = c("record_id" = "record_id", "pre_post" = "visit")) %>%
  tibble::column_to_rownames("cell") 
improve@meta.data[improve@meta.data$record_id == "RH-59-T", ]$record_id <- 'IT_07'
improve@meta.data[improve@meta.data$record_id == "RH-60-T", ]$record_id <- 'IT_08'
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$record_id <- 'IT_09'
improve@meta.data[improve@meta.data$record_id == "RH-66-T", ]$record_id <- 'IT_10'

########
## RH ##
########
rh <- subset(
  scrna,
  subset = cohort %in% c('RENAL HEIR', 'RENAL HEIRITAGE') &
    !scrna$record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T", "RH2-07-O")
)
rh$pre_post <- ifelse(rh$cohort == "RENAL HEIR", "pre", "post")
rh@meta.data <- rh@meta.data %>%
  tibble::rownames_to_column("cell") %>%   
  left_join(clin, by = c("record_id" = "record_id")) %>%
  tibble::column_to_rownames("cell") 
#
rh@meta.data[rh@meta.data$record_id == "RH2-14-T", ]$record_id <- 'RH-23-T'
rh@meta.data[rh@meta.data$record_id == "RH2-19-T", ]$record_id <- 'RH-67-T'

#
improve@meta.data$treatment <- "VSG"
rh@meta.data$treatment <- "Standard"

healthy <- subset(scrna, 
                  subset = record_id %in% c("CRC-03","CRC-14","CRC-02","CRC-13",
                                            "CRC-10","CRC-39","CRC-40","CRC-11",
                                            "CRC-46","CRC-58","CRC-56","CRC-54"))
healthy@meta.data$treatment <- "Healthy"
healthy@meta.data$pre_post <- "healthy"
#
combined <- merge(improve, rh)








###harmonized 


harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

dat <- dat %>% filter(visit == 'baseline') %>% 
  dplyr::select(record_id, race_ethnicity, diabetes_duration, sbp, dbp, hba1c, epic_insulin_1, 
                epic_sglti2_1, epic_glp1ra_1, epic_statin_1,  epic_fibrate_1,  
                ace_inhibitor, epic_raasi_1, epic_mfm_1, epic_tzd_1)







library(gtsummary)
library(dplyr)

# Select variables from dat
dat_subset <- dat %>%
  select(
    record_id,
    race_ethnicity,
    diabetes_duration,
    sbp,
    dbp,
    hba1c,
    epic_insulin_1,
    epic_sglti2_1,
    epic_glp1ra_1,
    epic_statin_1,
    epic_fibrate_1,
    epic_raasi_1,
    epic_mfm_1,
    epic_tzd_1
  )

# Get baseline data from combined Seurat object
baseline_data <- combined@meta.data %>%
  filter(pre_post == "pre", treatment %in% c("VSG", "Standard")) %>%
  distinct(record_id, .keep_all = TRUE)

# Merge with dat
patient_level <- baseline_data %>%
  left_join(dat_subset, by = "record_id")

# Check sample sizes
table(patient_level$treatment)



library(gtsummary)
library(dplyr)


# Create Table 1 with updated formatting
table1 <- patient_level %>%
  mutate(
    age_numeric = as.numeric(as.character(age)),
    microalbuminuria = ifelse(acr_u > 30, "Yes", "No")
  ) %>%
  select(
    treatment,
    age_numeric,
    sex,
    race_ethnicity,
    diabetes_duration,
    bmi,
    hba1c.y,
    sbp,
    dbp,
    eGFR_fas_cr_cysc,
    gfr_raw_plasma,
    gfr_bsa_plasma,
    acr_u,
    microalbuminuria,
    epic_insulin_1,
    epic_sglti2_1,
    epic_glp1ra_1,
    epic_tzd_1,
    epic_mfm_1,
    epic_raasi_1,
    epic_statin_1,
    epic_fibrate_1
  ) %>%
  tbl_summary(
    by = treatment,
    type = list(
      age_numeric ~ "continuous",
      acr_u ~ "continuous",
      microalbuminuria ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      acr_u ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age_numeric ~ 0,
      hba1c.y ~ 1,
      dbp ~ 0,
      acr_u ~ 1,
      gfr_raw_plasma ~ 0,
      gfr_bsa_plasma ~ 0
    ),
    label = list(
      age_numeric ~ "Age, years",
      sex ~ "Sex",
      race_ethnicity ~ "Race/Ethnicity",
      diabetes_duration ~ "Diabetes duration, years",
      bmi ~ "BMI, kg/m²",
      hba1c.y ~ "HbA1c, %",
      sbp ~ "SBP, mmHg",
      dbp ~ "DBP, mmHg",
      eGFR_fas_cr_cysc ~ "eGFR, mL/min/1.73m²",
      gfr_raw_plasma ~ "Measured GFR, mL/min",
      gfr_bsa_plasma ~ "Measured GFR (BSA-adjusted), mL/min/1.73m²",
      acr_u ~ "UACR, mg/g",
      microalbuminuria ~ "Microalbuminuria (UACR >30 mg/g)",
      epic_insulin_1 ~ "Insulin",
      epic_sglti2_1 ~ "SGLT2i",
      epic_glp1ra_1 ~ "GLP-1RA",
      epic_tzd_1 ~ "PPARγ agonist",
      epic_mfm_1 ~ "Metformin",
      epic_raasi_1 ~ "ACEi/ARB",
      epic_statin_1 ~ "Statin",
      epic_fibrate_1 ~ "Fenofibrate"
    ),
    missing = "no"
  ) %>%
  add_p() %>%
  bold_labels()

table1



output_dir <- "C:/Users/netio/Documents/UofW/Projects/Long_Paper/"
# Export to Word
table1 %>% 
  as_gt() %>% 
  gt::gtsave(paste0(output_dir, "Table1_VSG_vs_SMT.docx"))



















#Updated

# Study Design Overview Figure - VSG vs SMT vs Healthy Controls
# With complete participant verification

library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(dplyr)
library(readr)
library(Seurat)
library(gtsummary)
library(tidyr)
library(purrr)  # For negate() function

# ============================================================================
# SECTION 1: LOAD AND PROCESS DATA
# ============================================================================

base_dir <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/scrna/"
output_dir <- "C:/Users/netio/Documents/UofW/Projects/Long_Paper/"

# Load Seurat object
scrna <- readRDS(paste0(base_dir, "data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds"))
scrna <- subset(scrna, subset = percent.mt < 50 & nFeature_RNA < 5000 & nFeature_RNA > 500) 

# Load clinical data
clin <- read_csv(paste0(base_dir, "data_clean/pb90_rhrh2improve_clinical_subset.csv"), show_col_types = FALSE) %>%
  dplyr::select(record_id, visit, acr_u, gfr_bsa_plasma, gfr_raw_plasma, bmi, eGFR_CKiD_U25_avg, eGFR_fas_cr_cysc) %>%
  dplyr::mutate(visit = dplyr::recode(visit, "baseline" = "pre", "12_months_post_surgery" = "post"))

# Load harmonized data
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

dat <- dat %>% 
  filter(visit == 'baseline') %>% 
  dplyr::select(record_id, race_ethnicity, diabetes_duration, sbp, dbp, hba1c, 
                epic_insulin_1, epic_sglti2_1, epic_glp1ra_1, epic_statin_1, epic_fibrate_1,  
                ace_inhibitor, epic_raasi_1, epic_mfm_1, epic_tzd_1)

# ============================================================================
# SECTION 2: DEFINE COHORTS PROPERLY
# ============================================================================

# --- IMPROVE (VSG) cohort ---
improve <- subset(
  scrna, 
  subset = cohort == "IMPROVE" | 
    record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T")
)
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$pre_post <- 'pre'
improve@meta.data <- improve@meta.data %>%
  tibble::rownames_to_column("cell") %>%   
  left_join(clin, by = c("record_id" = "record_id", "pre_post" = "visit")) %>%
  tibble::column_to_rownames("cell") 

# Rename RH samples to IT
improve@meta.data[improve@meta.data$record_id == "RH-59-T", ]$record_id <- 'IT_07'
improve@meta.data[improve@meta.data$record_id == "RH-60-T", ]$record_id <- 'IT_08'
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$record_id <- 'IT_09'
improve@meta.data[improve@meta.data$record_id == "RH-66-T", ]$record_id <- 'IT_10'
improve@meta.data$treatment <- "VSG"

# --- RENAL HEIR (SMT) cohort ---
rh <- subset(
  scrna,
  subset = cohort %in% c('RENAL HEIR', 'RENAL HEIRITAGE') &
    !scrna$record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T", "RH2-07-O")
)
rh$pre_post <- ifelse(rh$cohort == "RENAL HEIR", "pre", "post")
rh@meta.data <- rh@meta.data %>%
  tibble::rownames_to_column("cell") %>%   
  left_join(clin, by = c("record_id" = "record_id")) %>%
  tibble::column_to_rownames("cell") 

# Rename follow-up samples
rh@meta.data[rh@meta.data$record_id == "RH2-14-T", ]$record_id <- 'RH-23-T'
rh@meta.data[rh@meta.data$record_id == "RH2-19-T", ]$record_id <- 'RH-67-T'
rh@meta.data$treatment <- "Standard"

# --- Healthy Controls ---
healthy_ids <- c("CRC-03", "CRC-14", "CRC-02", "CRC-13",
                 "CRC-10", "CRC-39", "CRC-40", "CRC-11",
                 "CRC-46", "CRC-58", "CRC-56", "CRC-54")

healthy <- subset(scrna, subset = record_id %in% healthy_ids)
healthy@meta.data$treatment <- "Healthy"
healthy@meta.data$pre_post <- "healthy"

# ============================================================================
# SECTION 3: VERIFY PARTICIPANT COUNTS
# ============================================================================

cat("=== PARTICIPANT VERIFICATION ===\n\n")

# VSG participants (baseline only)
vsg_ids <- improve@meta.data %>%
  filter(pre_post == "pre") %>%
  distinct(record_id) %>%
  pull(record_id)
cat("VSG (IMPROVE) baseline participants:", length(vsg_ids), "\n")
cat("IDs:", paste(sort(vsg_ids), collapse = ", "), "\n\n")

# SMT participants (baseline only)
smt_ids <- rh@meta.data %>%
  filter(pre_post == "pre") %>%
  distinct(record_id) %>%
  pull(record_id)
cat("SMT (RENAL HEIR) baseline participants:", length(smt_ids), "\n")
cat("IDs:", paste(sort(smt_ids), collapse = ", "), "\n\n")

# Healthy controls
hc_ids <- healthy@meta.data %>%
  distinct(record_id) %>%
  pull(record_id)
cat("Healthy Controls:", length(hc_ids), "\n")
cat("IDs:", paste(sort(hc_ids), collapse = ", "), "\n\n")

cat("TOTAL UNIQUE PARTICIPANTS:", length(vsg_ids) + length(smt_ids) + length(hc_ids), "\n")

# ============================================================================
# SECTION 4: CREATE PATIENT-LEVEL DATA FOR ALL GROUPS
# ============================================================================

# Combine all three groups
combined_all <- merge(improve, y = list(rh, healthy))

# Get baseline/single timepoint data for each group
patient_level_all <- combined_all@meta.data %>%
  filter(
    (treatment == "VSG" & pre_post == "pre") |
      (treatment == "Standard" & pre_post == "pre") |
      (treatment == "Healthy")
  ) %>%
  distinct(record_id, .keep_all = TRUE)

# Merge with harmonized data (for T2D groups)
dat_subset <- dat %>%
  dplyr::select(
    record_id, race_ethnicity, diabetes_duration, sbp, dbp, hba1c,
    epic_insulin_1, epic_sglti2_1, epic_glp1ra_1, epic_statin_1,
    epic_fibrate_1, epic_raasi_1, epic_mfm_1, epic_tzd_1
  )

patient_level_all <- patient_level_all %>%
  left_join(dat_subset, by = "record_id")

# Verify final counts
cat("\n=== FINAL PATIENT-LEVEL COUNTS ===\n")
print(table(patient_level_all$treatment))

# ============================================================================
# SECTION 5: CREATE DEMOGRAPHICS TABLE
# ============================================================================

# Create gtsummary table
table1_full <- patient_level_all %>%
  mutate(
    treatment = factor(treatment, levels = c("Healthy", "VSG", "Standard")),
    age_numeric = as.numeric(as.character(age)),
    microalbuminuria = ifelse(acr_u > 30, "Yes", "No"),
    sex_male = ifelse(sex == "Male", 1, 0)
  ) %>%
  dplyr::select(
    treatment,
    age_numeric,
    sex,
    race_ethnicity,
    diabetes_duration,
    bmi,
    hba1c.y,
    sbp,
    dbp,
    eGFR_fas_cr_cysc,
    gfr_raw_plasma,
    gfr_bsa_plasma,
    acr_u,
    microalbuminuria,
    epic_insulin_1,
    epic_sglti2_1,
    epic_glp1ra_1,
    epic_tzd_1,
    epic_mfm_1,
    epic_raasi_1,
    epic_statin_1,
    epic_fibrate_1
  ) %>%
  tbl_summary(
    by = treatment,
    type = list(
      age_numeric ~ "continuous",
      acr_u ~ "continuous",
      microalbuminuria ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      acr_u ~ "{median} [{p25}, {p75}]",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age_numeric ~ 1,
      diabetes_duration ~ 1,
      bmi ~ 1,
      hba1c.y ~ 1,
      eGFR_fas_cr_cysc ~ 0,
      gfr_raw_plasma ~ 0,
      gfr_bsa_plasma ~ 0,
      acr_u ~ 1
    ),
    label = list(
      age_numeric ~ "Age (years)",
      sex ~ "Sex",
      race_ethnicity ~ "Race/Ethnicity",
      diabetes_duration ~ "Diabetes duration (years)",
      bmi ~ "BMI (kg/m²)",
      hba1c.y ~ "HbA1c (%)",
      sbp ~ "SBP (mmHg)",
      dbp ~ "DBP (mmHg)",
      eGFR_fas_cr_cysc ~ "eGFR (mL/min/1.73m²)",
      gfr_raw_plasma ~ "Measured GFR (mL/min)",
      gfr_bsa_plasma ~ "Measured GFR BSA-adjusted (mL/min/1.73m²)",
      acr_u ~ "UACR (mg/g)",
      microalbuminuria ~ "Microalbuminuria (UACR >30)",
      epic_insulin_1 ~ "Insulin",
      epic_sglti2_1 ~ "SGLT2i",
      epic_glp1ra_1 ~ "GLP-1RA",
      epic_tzd_1 ~ "PPARγ agonist",
      epic_mfm_1 ~ "Metformin",
      epic_raasi_1 ~ "ACEi/ARB",
      epic_statin_1 ~ "Statin",
      epic_fibrate_1 ~ "Fenofibrate"
    ),
    missing = "no"
  ) %>%
  add_p() %>%
  bold_labels()

# Print table
table1_full

# ============================================================================
# SECTION 6: EXTRACT VALUES FOR FIGURE
# ============================================================================

# Helper functions
mean_sd <- function(x) {
  if (all(is.na(x))) return("—")
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

median_iqr <- function(x) {
  if (all(is.na(x))) return("—")
  sprintf("%.1f [%.1f, %.1f]", 
          median(x, na.rm = TRUE), 
          quantile(x, 0.25, na.rm = TRUE), 
          quantile(x, 0.75, na.rm = TRUE))
}

n_pct <- function(x, condition = TRUE) {
  if (all(is.na(x))) return("—")
  valid <- x[!is.na(x)]
  if (is.logical(condition)) {
    n <- sum(valid == condition | valid == 1, na.rm = TRUE)
  } else {
    n <- sum(valid == condition, na.rm = TRUE)
  }
  pct <- 100 * n / length(valid)
  sprintf("%d (%.0f%%)", n, pct)
}

# Split by group
hc <- patient_level_all %>% filter(treatment == "Healthy")
vsg <- patient_level_all %>% filter(treatment == "VSG")
smt <- patient_level_all %>% filter(treatment == "Standard")

# ============================================================================
# DEBUG: Check medication columns and values
# ============================================================================

cat("\n=== CHECKING MEDICATION DATA ===\n")

# Check what columns exist
cat("\nColumns in patient_level_all containing 'epic':\n")
print(grep("epic", names(patient_level_all), value = TRUE))

# Check unique values in medication columns
cat("\nUnique values in epic_mfm_1 (Metformin):\n")
print(table(patient_level_all$epic_mfm_1, useNA = "always"))

cat("\nUnique values in epic_insulin_1:\n")
print(table(patient_level_all$epic_insulin_1, useNA = "always"))

# Check if medications are coded as "Yes"/"No" instead of 1/0
cat("\nSample of medication data:\n
")
print(patient_level_all %>% 
        filter(treatment != "Healthy") %>% 
        dplyr::select(record_id, treatment, starts_with("epic_")) %>% 
        head(10))

# Check HC BMI from original metadata
cat("\n=== CHECKING HEALTHY CONTROL DATA ===\n")
cat("\nHealthy control BMI values from Seurat metadata:\n")
print(hc %>% dplyr::select(record_id, bmi, age, sex) %>% head(15))

# ============================================================================
# FIXED: Helper function for medications (handles different codings)
# ============================================================================

n_pct_med <- function(x) {
  if (all(is.na(x))) return("—")
  valid <- x[!is.na(x)]
  # Handle both numeric (1/0) and character ("Yes"/"No") coding
  if (is.numeric(valid)) {
    n <- sum(valid == 1, na.rm = TRUE)
  } else {
    n <- sum(valid %in% c("Yes", "yes", "Y", "1", "TRUE", "True"), na.rm = TRUE)
  }
  pct <- 100 * n / length(valid)
  sprintf("%d (%.0f%%)", n, pct)
}

# Build demographics tribble - BASELINE characteristics
demographics <- tribble(
  ~Variable,                              ~HC,                                      ~VSG,                                     ~SMT,
  "N (baseline)",                         as.character(nrow(hc)),                   as.character(nrow(vsg)),                  as.character(nrow(smt)),
  "Age (years)",                          mean_sd(as.numeric(hc$age)),              mean_sd(as.numeric(vsg$age)),             mean_sd(as.numeric(smt$age)),
  "Sex, Male",                            n_pct(hc$sex, "Male"),                    n_pct(vsg$sex, "Male"),                   n_pct(smt$sex, "Male"),
  "Diabetes duration (years)",            "N/A",                                    mean_sd(vsg$diabetes_duration),           mean_sd(smt$diabetes_duration),
  "BMI (kg/m²)",                          mean_sd(hc$bmi),                          mean_sd(vsg$bmi),                         mean_sd(smt$bmi),
  "HbA1c (%)",                            mean_sd(hc$hba1c.y),                      mean_sd(vsg$hba1c.y),                     mean_sd(smt$hba1c.y),
  "eGFR (mL/min/1.73m²)",                 mean_sd(hc$eGFR_fas_cr_cysc),             mean_sd(vsg$eGFR_fas_cr_cysc),            mean_sd(smt$eGFR_fas_cr_cysc),
  "UACR (mg/g)",                          median_iqr(hc$acr_u),                     median_iqr(vsg$acr_u),                    median_iqr(smt$acr_u),
  "Metformin",                            "N/A",                                    n_pct_med(vsg$epic_mfm_1),                n_pct_med(smt$epic_mfm_1),
  "Insulin",                              "N/A",                                    n_pct_med(vsg$epic_insulin_1),            n_pct_med(smt$epic_insulin_1),
  "SGLT2i",                               "N/A",                                    n_pct_med(vsg$epic_sglti2_1),             n_pct_med(smt$epic_sglti2_1),
  "ACEi/ARB",                             "N/A",                                    n_pct_med(vsg$epic_raasi_1),              n_pct_med(smt$epic_raasi_1)
)

# ============================================================================
# SECTION 6B: STUDY DESIGN INFO - for inclusion criteria box
# ============================================================================

# Count paired samples (pre AND post)
n_vsg_paired <- improve@meta.data %>%
  filter(treatment == "VSG") %>%
  group_by(record_id) %>%
  summarise(has_both = n_distinct(pre_post) == 2) %>%
  filter(has_both) %>%
  nrow()

n_smt_paired <- rh@meta.data %>%
  group_by(record_id) %>%
  summarise(has_both = n_distinct(pre_post) == 2) %>%
  filter(has_both) %>%
  nrow()

cat("\n=== PAIRED SAMPLES FOR LONGITUDINAL ANALYSIS ===\n")
cat("VSG paired (pre + post):", n_vsg_paired, "\n")
cat("SMT paired (pre + post):", n_smt_paired, "\n")

print(demographics)

# ============================================================================
# SECTION 7: COLOR DEFINITIONS
# ============================================================================

col_hc <- "#81C784"
col_hc_light <- "#C8E6C9"
col_vsg <- "#7986CB" 
col_vsg_light <- "#C5CAE9"
col_smt <- "#FFB74D"
col_smt_light <- "#FFE0B2"

# ============================================================================
# SECTION 8: CREATE FIGURE COMPONENTS
# ============================================================================

create_person_shape <- function() {
  theta <- seq(0, 2 * pi, length.out = 50)
  head <- data.frame(x = 0.25 * cos(theta), y = 0.25 * sin(theta) + 1.6, part = "head")
  body <- data.frame(
    x = c(-0.35, 0.35, 0.25, -0.25, -0.35),
    y = c(0.2, 0.2, 1.25, 1.25, 0.2),
    part = "body"
  )
  list(head = head, body = body)
}

create_group_icon <- function(fill_color, n_people = 3, label = "", label_color = "black") {
  shapes <- create_person_shape()
  
  all_data <- lapply(1:n_people, function(i) {
    offset <- (i - (n_people + 1) / 2) * 0.7
    head_i <- shapes$head %>% mutate(x = x + offset, person = i)
    body_i <- shapes$body %>% mutate(x = x + offset, person = i)
    bind_rows(head_i, body_i)
  }) %>% bind_rows()
  
  ggplot() +
    geom_polygon(data = all_data %>% filter(part == "head"),
                 aes(x = x, y = y, group = person),
                 fill = fill_color, color = NA) +
    geom_polygon(data = all_data %>% filter(part == "body"),
                 aes(x = x, y = y, group = person),
                 fill = fill_color, color = NA) +
    annotate("text", x = 0, y = -0.3, label = label, 
             size = 5, fontface = "bold", color = label_color) +
    coord_fixed(ratio = 1, xlim = c(-1.5, 1.5), ylim = c(-0.8, 2.2)) +
    theme_void()
}

create_demographics_table <- function(demo_data, col_hc, col_vsg, col_smt) {
  n_rows <- nrow(demo_data)
  n_cols <- ncol(demo_data)
  
  fill_matrix <- matrix("white", nrow = n_rows, ncol = n_cols)
  fill_matrix[, 2] <- col_hc
  fill_matrix[, 3] <- col_vsg
  fill_matrix[, 4] <- col_smt
  
  for (i in seq(2, n_rows, by = 2)) {
    fill_matrix[i, 1] <- "gray95"
  }
  
  tt <- ttheme_minimal(
    core = list(
      bg_params = list(fill = fill_matrix, col = "gray60", lwd = 0.5),
      fg_params = list(fontsize = 10, hjust = 0.5, x = 0.5),
      padding = unit(c(8, 6), "mm")  # More padding for readability
    ),
    colhead = list(
      bg_params = list(fill = c("gray90", col_hc, col_vsg, col_smt), 
                       col = "gray60", lwd = 0.5),
      fg_params = list(fontsize = 11, fontface = "bold", hjust = 0.5, x = 0.5),
      padding = unit(c(8, 6), "mm")
    )
  )
  
  tableGrob(demo_data, rows = NULL, theme = tt)
}

# ============================================================================
# SECTION 9: ASSEMBLE FINAL FIGURE
# ============================================================================

create_study_design_figure <- function(demographics, n_hc, n_vsg, n_smt,
                                       save_path = NULL, width = 14, height = 10) {
  
  # Create icons with actual N
  icon_hc <- create_group_icon(col_hc, n_people = 3, 
                               label = paste0("Healthy Controls\n(n=", n_hc, ")"), 
                               label_color = "#2E7D32")
  
  icon_vsg <- create_group_icon(col_vsg, n_people = 3, 
                                label = paste0("VSG\n(n=", n_vsg, ")"), 
                                label_color = "#303F9F")
  
  icon_smt <- create_group_icon(col_smt, n_people = 3, 
                                label = paste0("Standard Therapy\n(n=", n_smt, ")"), 
                                label_color = "#E65100")
  
  table_grob <- create_demographics_table(demographics, col_hc_light, col_vsg_light, col_smt_light)
  
  icons_panel <- plot_grid(icon_hc, icon_vsg, icon_smt, nrow = 1, rel_widths = c(1, 1, 1))
  
  # Create table as ggplot for better positioning control
  table_plot <- ggdraw() + 
    draw_grob(table_grob, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5)
  
  # Combine icons and table vertically
  top_section <- ggdraw() +
    draw_label("Bariatric Surgery Drives Molecular Kidney Remodeling Compared to\nStandard Medical Therapy in Youth-Onset Type 2 Diabetes",
               x = 0.5, y = 0.85, size = 13, fontface = "bold", hjust = 0.5) +
    draw_label("Multi-omic kidney profiling from IMPROVE-T2D (VSG), RENAL-HEIR/HEIRitage (SMT), and healthy controls",
               x = 0.5, y = 0.70, size = 10, fontface = "italic", hjust = 0.5, color = "gray40") +
    draw_plot(icons_panel, x = 0.05, y = 0.02, width = 0.9, height = 0.62)
  
  # Stack vertically using plot_grid
  final_figure <- plot_grid(
    top_section,
    table_plot,
    ncol = 1,
    rel_heights = c(0.4, 0.6)
  )
  
  if (!is.null(save_path)) {
    ggsave(save_path, final_figure, width = width, height = height, dpi = 300, bg = "white")
    message("Figure saved to: ", save_path)
  }
  
  return(final_figure)
}

# ============================================================================
# SECTION 10: GENERATE FIGURE
# ============================================================================

# Check if paired counts exist, otherwise set to NA
if (!exists("n_vsg_paired")) n_vsg_paired <- NA
if (!exists("n_smt_paired")) n_smt_paired <- NA

fig <- create_study_design_figure(
  demographics = demographics,
  n_hc = nrow(hc),
  n_vsg = nrow(vsg),
  n_smt = nrow(smt),
  save_path = paste0(output_dir, "Figure0_Study_Design.pdf"),
  width = 14,
  height = 12
)

print(fig)

# Save table to Word
table1_full %>% 
  as_gt() %>% 
  gt::gtsave(paste0(output_dir, "Table1_HC_VSG_SMT.docx"))




