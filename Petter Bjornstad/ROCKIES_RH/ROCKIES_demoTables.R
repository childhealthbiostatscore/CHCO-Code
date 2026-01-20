###########ROCKIES Make Tables as .docx

#######################################################################
# ROCKIES Publication - Demographics Tables for Word Documents
# Based on existing analysis script data preparation
########################################################################

library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)
library(data.table)
library(gt)
library(webshot2)  # For PNG export - install if needed: install.packages("webshot2")

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

########################################################################
# COMMON DATA PREPARATION FUNCTIONS
########################################################################

PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2,
                                   lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>% dplyr::select(lc_k2, rc_k2) %>% rowMeans(na.rm=T)
  avg_m_k2 <- tmp_df %>% dplyr::select(lm_k2, rm_k2) %>% rowMeans(na.rm=T)
  avg_c_f <- tmp_df %>% dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  avg_m_f <- tmp_df %>% dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  avg_m_k2_f <- avg_m_k2 / avg_m_f
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f')
  return(results)
}

# Function to get SGLT2 info and filter
get_sglt2_filtered_data <- function(dat_results) {
  dat2 <- dat_results
  dat2$group2 <- NA
  need_med_info <- dat2 %>% filter(is.na(group2))
  
  RH <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RENALHEIR-SGLT2.csv')
  names(RH) <- c('Subject', 'rep_instr', 'rep_inst', 'SGLT2')
  RH2 <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RenalHEIRitage-SGLT2Use.csv')
  names(RH2) <- c('Subject', 'event', 'rep_instr', 'rep_inst', 'mrn', 'SGLT2', 'SGLT2_ever')
  RH2 <- RH2 %>% filter(!is.na(mrn))
  improve <- data.table::fread('C:/Users/netio/Downloads/IMPROVET2D-SGLT2i_DATA_LABELS_2025-08-25_0938.csv')
  names(improve)[5] <- 'SGLT2'
  names(improve)[1] <- 'record_id'
  improve <- improve %>% filter(!is.na(SGLT2)) %>% filter(SGLT2 != '')
  
  improve_small <- improve %>% filter(record_id %in% need_med_info$record_id)
  RH_small <- RH %>% filter(Subject %in% need_med_info$record_id)
  RH2_small <- RH2 %>% filter(mrn %in% need_med_info$mrn)
  
  for(i in c(1:nrow(RH_small))){
    if(nrow(RH_small) == 0) next
    if(RH_small$SGLT2[i] == 'No'){
      dat2$group2[which(dat2$record_id == RH_small$Subject[i])] <- 'T2D-No SGLTi2'
      dat2$epic_sglti2_1[which(dat2$record_id == RH_small$Subject[i])] <- 'No'
    }else if(RH_small$SGLT2[i] == 'Yes'){
      dat2$group2[which(dat2$record_id == RH_small$Subject[i])] <- 'T2D-SGLTi2'
      dat2$epic_sglti2_1[which(dat2$record_id == RH_small$Subject[i])] <- 'Yes'
    }
  }
  
  for(i in c(1:nrow(RH2_small))){
    if(nrow(RH2_small) == 0) next
    if(RH2_small$SGLT2[i] == 'No'){
      dat2$group2[which(dat2$mrn == RH2_small$mrn[i])] <- 'T2D-No SGLTi2'
      dat2$epic_sglti2_1[which(dat2$mrn == RH2_small$mrn[i])] <- 'No'
    }else if(RH2_small$SGLT2[i] == 'Yes'){
      dat2$group2[which(dat2$mrn == RH2_small$mrn[i])] <- 'T2D-SGLTi2'
      dat2$epic_sglti2_1[which(dat2$mrn == RH2_small$mrn[i])] <- 'Yes'
    }
  }
  
  for(i in c(1:nrow(improve_small))){
    if(nrow(improve_small) == 0) next
    if(improve_small$SGLT2[i] == 'No'){
      dat2$group2[which(dat2$record_id == improve_small$record_id[i])] <- 'T2D-No SGLTi2'
      dat2$epic_sglti2_1[which(dat2$record_id == improve_small$record_id[i])] <- 'No'
    }else if(improve_small$SGLT2[i] == 'Yes'){
      dat2$group2[which(dat2$record_id == improve_small$record_id[i])] <- 'T2D-SGLTi2'
      dat2$epic_sglti2_1[which(dat2$record_id == improve_small$record_id[i])] <- 'Yes'
    }
  }
  
  dat2$epic_sglti2_1[which(dat2$group == 'Lean Control')] <- 'No'
  dat2 <- dat2 %>% filter(epic_sglti2_1 != 'Yes')
  
  return(dat2)
}

########################################################################
# TABLE 1: scRNAseq Cohort (from Step 4 in your script)
########################################################################

cat("\n=== TABLE 1: scRNAseq Cohort ===\n")

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

# Get scRNAseq participants from Seurat metadata
load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')
# Get scRNAseq participants from Seurat metadata
load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')
tmp_meta <- so_kpmp_sc@meta.data
scrnaseq_ids <- tmp_meta %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

# Remove co-enrolled participant (keep only first biopsy)
scrnaseq_ids <- scrnaseq_ids %>% filter(record_id != "RH2-14-T")

remove(so_kpmp_sc)

dat_scrnaseq <- dat %>% 
  semi_join(scrnaseq_ids, by='record_id') %>% 
  filter(visit == 'baseline') %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

cat("scRNAseq cohort N:", nrow(dat_scrnaseq), "\n")

# Create Table 1
table1 <- dat_scrnaseq %>%
  mutate(
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    acr_u = as.numeric(acr_u),
    microalbuminuria = factor(ifelse(acr_u >= 30, "Yes", "No"), levels = c("No", "Yes")),
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    group = as.factor(group)
  ) %>%
  dplyr::select(age, sex, race_ethnicity, bmi, hba1c, acr_u, microalbuminuria, group) %>%
  tbl_summary(
    by = group,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      acr_u ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(age ~ 1, bmi ~ 1, hba1c ~ 2, acr_u ~ 1, all_categorical() ~ c(0, 1)),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      acr_u ~ "UACR, mg/g",
      microalbuminuria ~ "Microalbuminuria (UACR ≥30)"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "fisher.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD); Median (Q1, Q3); n (%)")

########################################################################
# TABLE 2: PET/CT Cohort (from Step 1 - Aim 1)
########################################################################

cat("\n=== TABLE 2: PET/CT Cohort ===\n")

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

# Remove existing PET columns to avoid duplicates
dat <- dat %>% dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f')))

tmp_results <- PET_avg(dat)
dat_results <- dat %>% bind_cols(tmp_results)
dat_results <- dat_results %>% filter(!is.na(avg_c_k2))
dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

dat2 <- get_sglt2_filtered_data(dat_results)

cat("PET/CT cohort N:", nrow(dat2), "\n")
print(table(dat2$group))

# Create Table 2
table2 <- dat2 %>%
  mutate(
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    acr_u = as.numeric(acr_u),
    microalbuminuria = factor(ifelse(acr_u >= 30, "Yes", "No"), levels = c("No", "Yes")),
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    group = as.factor(group)
  ) %>%
  dplyr::select(age, sex, race_ethnicity, bmi, hba1c, acr_u, microalbuminuria, group) %>%
  tbl_summary(
    by = group,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      acr_u ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(age ~ 1, bmi ~ 1, hba1c ~ 2, acr_u ~ 1, all_categorical() ~ c(0, 1)),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      acr_u ~ "UACR, mg/g",
      microalbuminuria ~ "Microalbuminuria (UACR ≥30)"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "fisher.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD); Median (Q1, Q3); n (%)")

########################################################################
# TABLE 3: UACR-PET Analysis Cohort (from Step 3)
########################################################################

cat("\n=== TABLE 3: UACR-PET Cohort ===\n")

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

# Remove existing PET columns to avoid duplicates
dat <- dat %>% dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f')))

tmp_results <- PET_avg(dat)
dat_results <- dat %>% bind_cols(tmp_results)
dat_results <- dat_results %>% filter(!is.na(avg_c_k2))
dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes'))

cat("UACR-PET cohort N:", nrow(dat_results), "\n")
print(table(dat_results$group))

# Create Table 3
table3 <- dat_results %>%
  mutate(
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    acr_u = as.numeric(acr_u),
    microalbuminuria = factor(ifelse(acr_u >= 30, "Yes", "No"), levels = c("No", "Yes")),
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    group = as.factor(group)
  ) %>%
  dplyr::select(age, sex, race_ethnicity, bmi, hba1c, acr_u, microalbuminuria, group) %>%
  tbl_summary(
    by = group,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      acr_u ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(age ~ 1, bmi ~ 1, hba1c ~ 2, acr_u ~ 1, all_categorical() ~ c(0, 1)),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      acr_u ~ "UACR, mg/g",
      microalbuminuria ~ "Microalbuminuria (UACR ≥30)"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ "fisher.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD); Median (Q1, Q3); n (%)")

########################################################################
# TABLE 4 & 5: GBM and Arteriosclerosis (from Step 2)
########################################################################

cat("\n=== TABLE 4 & 5: GBM/Arteriosclerosis Cohort ===\n")

# Load clinical data with GBM and arteriosclerosis
dat_clinical <- data.table::fread("C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
dat_clinical <- dat_clinical[-stringr::str_which(dat_clinical$record_id, '-O')]

dat_clinical <- dat_clinical %>% 
  dplyr::select(-any_of(c("age", "sex", "race_ethnicity", "bmi", "hba1c", "acr_u", "group")))

# Reload harmonized data to get demographics
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

# Filter to baseline visits before joining
dat_baseline <- dat %>% filter(visit == "baseline")

# Check what record_ids overlap
cat("record_ids in dat_clinical:", length(unique(dat_clinical$record_id)), "\n")
cat("record_ids in dat_baseline:", length(unique(dat_baseline$record_id)), "\n")
cat("Overlapping record_ids:", length(intersect(dat_clinical$record_id, dat_baseline$record_id)), "\n")

# Merge clinical data with demographics from harmonized data (including acr_u)
dat_clinical_merged <- dat_clinical %>%
  left_join(dat_baseline %>% dplyr::select(record_id, age, sex, race_ethnicity, bmi, hba1c, acr_u, group), 
            by = "record_id")

cat("GBM/Arteriosclerosis cohort N:", nrow(dat_clinical_merged), "\n")
cat("acr_u available:", sum(!is.na(dat_clinical_merged$acr_u)), "out of", nrow(dat_clinical_merged), "\n")
# Prepare data for Table 4 (GBM)
dat_gbm <- dat_clinical_merged %>%
  filter(!is.na(`GBM thickness`) & `GBM thickness` != '') %>%
  mutate(
    GBM_Status = factor(ifelse(`GBM thickness` == "yes", "Yes", "No"), levels = c("No", "Yes")),
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    acr_u = as.numeric(acr_u),
    microalbuminuria = factor(ifelse(acr_u >= 30, "Yes", "No"), levels = c("No", "Yes")),
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    group = as.factor(group),
    arteriosclerosis = as.factor(arteriosclerosis)
  )

cat("Table 4 (GBM) N:", nrow(dat_gbm), "\n")

# Create Table 4
table4 <- dat_gbm %>%
  dplyr::select(age, sex, race_ethnicity, bmi, hba1c, acr_u, microalbuminuria, group, 
                arteriosclerosis, GBM_Status) %>%
  tbl_summary(
    by = GBM_Status,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      acr_u ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(age ~ 1, bmi ~ 1, hba1c ~ 2, acr_u ~ 1, all_categorical() ~ c(0, 1)),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      acr_u ~ "UACR, mg/g",
      microalbuminuria ~ "Microalbuminuria (UACR ≥30)",
      group ~ "Group",
      arteriosclerosis ~ "Arteriosclerosis"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "fisher.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD); Median (Q1, Q3); n (%)")

# Prepare data for Table 5 (Arteriosclerosis)
dat_arterio <- dat_clinical_merged %>%
  filter(!is.na(arteriosclerosis) & arteriosclerosis != '') %>%
  mutate(
    Arterio_Status = factor(ifelse(arteriosclerosis == "yes", "Yes", "No"), levels = c("No", "Yes")),
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    acr_u = as.numeric(acr_u),
    microalbuminuria = factor(ifelse(acr_u >= 30, "Yes", "No"), levels = c("No", "Yes")),
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    group = as.factor(group),
    `GBM thickness` = as.factor(`GBM thickness`)
  )

cat("Table 5 (Arteriosclerosis) N:", nrow(dat_arterio), "\n")

# Create Table 5
table5 <- dat_arterio %>%
  dplyr::select(age, sex, race_ethnicity, bmi, hba1c, acr_u, microalbuminuria, group,
                `GBM thickness`, Arterio_Status) %>%
  tbl_summary(
    by = Arterio_Status,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      acr_u ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(age ~ 1, bmi ~ 1, hba1c ~ 2, acr_u ~ 1, all_categorical() ~ c(0, 1)),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      acr_u ~ "UACR, mg/g",
      microalbuminuria ~ "Microalbuminuria (UACR ≥30)",
      group ~ "Group",
      `GBM thickness` ~ "GBM Thickening"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "fisher.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD); Median (Q1, Q3); n (%)")

########################################################################
# EXPORT TABLES AS PNG AND HTML
########################################################################

# Function to save table as PNG and HTML
save_table_files <- function(tbl, filename_base, title) {
  
  # Convert to gt table with formatting
  gt_tbl <- tbl %>%
    as_gt() %>%
    tab_header(title = title) %>%
    tab_options(
      table.font.size = px(12),
      heading.title.font.size = px(14),
      heading.title.font.weight = "bold",
      column_labels.font.weight = "bold",
      table.width = pct(100)
    )
  
  # Save as HTML
  gtsave(gt_tbl, filename = paste0(base_path, filename_base, ".html"))
  cat(paste0(filename_base, ".html saved!\n"))
  
  # Save as PNG
  gtsave(gt_tbl, filename = paste0(base_path, filename_base, ".png"), vwidth = 800, vheight = 900)
  cat(paste0(filename_base, ".png saved!\n"))
}

# Save all tables
save_table_files(table1, "Table1_scRNAseq_Demographics", 
                 "Table 1. Participant Characteristics for Single-Cell RNA Sequencing Cohort")

save_table_files(table2, "Table2_PET_CT_Demographics",
                 "Table 2. Participant Characteristics for PET/CT Cohorts")

save_table_files(table3, "Table3_UACR_PET_Demographics",
                 "Table 3. Participant Characteristics for UACR-PET Analysis Cohort")

save_table_files(table4, "Table4_GBM_Demographics",
                 "Table 4. Participant Characteristics Stratified by GBM Thickening Status")

save_table_files(table5, "Table5_Arteriosclerosis_Demographics",
                 "Table 5. Participant Characteristics Stratified by Arteriosclerosis Status")

cat("\n=== ALL PNG AND HTML FILES SAVED ===\n")
cat("You can now insert the PNG images into your Word document.\n")
cat("Or open the HTML files in a browser, select all, and paste into Word.\n")






















