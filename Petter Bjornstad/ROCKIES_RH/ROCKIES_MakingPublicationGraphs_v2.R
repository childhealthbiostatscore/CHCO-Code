########################################################################
# ROCKIES Publication - Demographics Tables for Word Documents
########################################################################

library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

########################################################################
# DATA PREPARATION (from your existing code)
########################################################################

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

# Load clinical data for GBM/Arteriosclerosis tables
dat_clinical <- fread("C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
dat_clinical <- dat_clinical[-str_which(dat_clinical$record_id, '-O')]

########################################################################
# HELPER FUNCTION: Create publication-ready table
########################################################################

create_demo_table <- function(data, stratify_var, title, 
                               vars_to_include = NULL,
                               var_labels = NULL) {
  
  # Default variables if not specified
  if (is.null(vars_to_include)) {
    vars_to_include <- c("age", "sex", "race", "ethnicity", "bmi", 
                         "hba1c", "egfr_ckd_epi", "acr_u", 
                         "diabetes_duration", "sbp", "dbp")
  }
  
  # Filter to available variables
  available_vars <- vars_to_include[vars_to_include %in% names(data)]
  
  # Default labels
  default_labels <- list(
    age = "Age (years)",
    sex = "Sex",
    race = "Race",
    ethnicity = "Ethnicity", 
    bmi = "BMI (kg/m²)",
    hba1c = "HbA1c (%)",
    egfr_ckd_epi = "eGFR (mL/min/1.73m²)",
    acr_u = "UACR (mg/g)",
    diabetes_duration = "Diabetes Duration (years)",
    sbp = "Systolic BP (mmHg)",
    dbp = "Diastolic BP (mmHg)",
    gfr_bsa_plasma = "mGFR (mL/min/1.73m²)",
    avg_c_k2 = "Cortical K2",
    avg_c_f = "Cortical F",
    avg_c_k2_f = "Cortical K2/F"
  )
  
  if (!is.null(var_labels)) {
    default_labels <- modifyList(default_labels, var_labels)
  }
  
  # Create the table
  tbl <- data %>%
    dplyr::select(all_of(c(stratify_var, available_vars))) %>%
    tbl_summary(
      by = all_of(stratify_var),
      missing = "ifany",
      missing_text = "Missing",
      label = default_labels[available_vars],
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = list(
        all_continuous() ~ 1,
        all_categorical() ~ c(0, 1)
      )
    ) %>%
    add_p(
      test = list(
        all_continuous() ~ "wilcox.test",
        all_categorical() ~ "fisher.test"
      )
    ) %>%
    add_overall() %>%
    modify_header(label = "**Characteristic**") %>%
    modify_spanning_header(c("stat_1", "stat_2") ~ paste0("**", stratify_var, "**")) %>%
    bold_labels()
  
  return(tbl)
}

########################################################################
# TABLE 1: scRNAseq Cohort (32 participants with kidney biopsy)
########################################################################

# Filter for scRNAseq cohort - adjust filter criteria based on your data
# Assuming you have a variable indicating scRNAseq participants
dat_scrnaseq <- dat %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes')) %>%
  # Add any additional filters for scRNAseq cohort
  mutate(
    T2D_Status = ifelse(group == "Type 2 Diabetes", "T2D", "Control")
  )

# Create Table 1
table1 <- dat_scrnaseq %>%
  create_demo_table(
    stratify_var = "T2D_Status",
    title = "Table 1. Participant Characteristics for Single-Cell RNA Sequencing Cohort",
    vars_to_include = c("age", "sex", "race", "ethnicity", "bmi", 
                        "hba1c", "egfr_ckd_epi", "acr_u", "diabetes_duration")
  )

########################################################################
# TABLE 2: PET/CT Cohorts (18 T2D + 11 HC)
########################################################################

# Use dat2 from your existing code (after SGLT2i filtering)
# Recreate dat2 or use the filtered version
PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2, lc_f, rc_f, lm_f, rm_f)
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

dat_results <- dat %>% 
  dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f'))) %>%
  bind_cols(PET_avg(.)) %>%
  filter(!is.na(avg_c_k2)) %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

dat_pet <- dat_results %>%
  mutate(
    Cohort = case_when(
      group == "Type 2 Diabetes" ~ "T2D (RENAL-HEIRitage)",
      group == "Lean Control" ~ "Healthy Control (CROCODILE)"
    )
  )

table2 <- dat_pet %>%
  create_demo_table(
    stratify_var = "Cohort",
    title = "Table 2. Participant Characteristics for PET/CT Cohorts",
    vars_to_include = c("age", "sex", "race", "ethnicity", "bmi", 
                        "hba1c", "egfr_ckd_epi", "acr_u", "sbp", "dbp",
                        "avg_c_k2", "avg_c_f", "avg_c_k2_f")
  )

########################################################################
# TABLE 3: UACR-PET Analysis Cohort (40 participants)
########################################################################

dat_uacr <- dat %>%
  dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f'))) %>%
  bind_cols(PET_avg(.)) %>%
  filter(!is.na(avg_c_k2)) %>%
  filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes')) %>%
  filter(record_id != 'CRC-55') %>%
  mutate(
    Group = case_when(
      group == "Type 2 Diabetes" ~ "T2D",
      group == "Lean Control" ~ "Lean Control",
      group == "Obese Control" ~ "Obese Control"
    )
  )

table3 <- dat_uacr %>%
  create_demo_table(
    stratify_var = "Group",
    title = "Table 3. Participant Characteristics for UACR-PET Analysis Cohort",
    vars_to_include = c("age", "sex", "race", "ethnicity", "bmi", 
                        "hba1c", "egfr_ckd_epi", "acr_u",
                        "avg_c_k2", "avg_c_f", "avg_c_k2_f")
  )

########################################################################
# TABLE 4: GBM Thickening Status
########################################################################

dat_gbm <- dat_clinical %>%
  filter(`GBM thickness` != '' & !is.na(`GBM thickness`)) %>%
  mutate(
    GBM_Status = ifelse(`GBM thickness` == "yes", "GBM Thickening", "No GBM Thickening")
  )

# Merge with full dataset for demographics
dat_gbm_full <- dat_gbm %>%
  left_join(dat %>% dplyr::select(record_id, age, sex, race, ethnicity, bmi, 
                                   hba1c, egfr_ckd_epi, acr_u, diabetes_duration),
            by = "record_id")

table4 <- dat_gbm_full %>%
  create_demo_table(
    stratify_var = "GBM_Status",
    title = "Table 4. Participant Characteristics Stratified by GBM Thickening Status",
    vars_to_include = c("age", "sex", "race", "bmi", "hba1c", 
                        "egfr_ckd_epi", "acr_u", "avg_c_k2", "avg_c_f", "avg_c_k2_f")
  )

########################################################################
# TABLE 5: Arteriosclerosis Status
########################################################################

dat_arterio <- dat_clinical %>%
  filter(arteriosclerosis != '' & !is.na(arteriosclerosis)) %>%
  mutate(
    Arterio_Status = ifelse(arteriosclerosis == "yes", "Arteriosclerosis", "No Arteriosclerosis")
  )

dat_arterio_full <- dat_arterio %>%
  left_join(dat %>% dplyr::select(record_id, age, sex, race, ethnicity, bmi, 
                                   hba1c, egfr_ckd_epi, acr_u, diabetes_duration),
            by = "record_id")

table5 <- dat_arterio_full %>%
  create_demo_table(
    stratify_var = "Arterio_Status",
    title = "Table 5. Participant Characteristics Stratified by Arteriosclerosis Status",
    vars_to_include = c("age", "sex", "race", "bmi", "hba1c", 
                        "egfr_ckd_epi", "acr_u", "avg_c_k2", "avg_c_f", "avg_c_k2_f")
  )

########################################################################
# EXPORT TO WORD DOCUMENTS
########################################################################

# Function to save table as Word document
save_table_to_word <- function(tbl, filename, title) {
  
  # Convert to flextable
  ft <- tbl %>%
    as_flex_table() %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    autofit() %>%
    set_caption(caption = title)
  
  # Create Word document
  doc <- read_docx() %>%
    body_add_par(title, style = "heading 1") %>%
    body_add_flextable(ft)
  
  # Save
  print(doc, target = paste0(base_path, filename))
  print(paste0(filename, " saved!"))
}

# Save individual tables
save_table_to_word(table1, "Table1_scRNAseq_Demographics.docx", 
                   "Table 1. Participant Characteristics for Single-Cell RNA Sequencing Cohort")

save_table_to_word(table2, "Table2_PET_CT_Demographics.docx",
                   "Table 2. Participant Characteristics for PET/CT Cohorts")

save_table_to_word(table3, "Table3_UACR_PET_Demographics.docx",
                   "Table 3. Participant Characteristics for UACR-PET Analysis Cohort")

save_table_to_word(table4, "Table4_GBM_Demographics.docx",
                   "Table 4. Participant Characteristics Stratified by GBM Thickening Status")

save_table_to_word(table5, "Table5_Arteriosclerosis_Demographics.docx",
                   "Table 5. Participant Characteristics Stratified by Arteriosclerosis Status")

########################################################################
# COMBINED WORD DOCUMENT WITH ALL TABLES
########################################################################

# Create combined document
combined_doc <- read_docx() %>%
  body_add_par("ROCKIES Publication - Demographics Tables", style = "heading 1") %>%
  body_add_par("") %>%
  
  # Table 1
  body_add_par("Table 1. Participant Characteristics for Single-Cell RNA Sequencing Cohort", 
               style = "heading 2") %>%
  body_add_flextable(table1 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit()) %>%
  body_add_par("") %>%
  body_add_break() %>%
  
  # Table 2
  body_add_par("Table 2. Participant Characteristics for PET/CT Cohorts", 
               style = "heading 2") %>%
  body_add_flextable(table2 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit()) %>%
  body_add_par("") %>%
  body_add_break() %>%
  
  # Table 3
  body_add_par("Table 3. Participant Characteristics for UACR-PET Analysis Cohort", 
               style = "heading 2") %>%
  body_add_flextable(table3 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit()) %>%
  body_add_par("") %>%
  body_add_break() %>%
  
  # Table 4
  body_add_par("Table 4. Participant Characteristics Stratified by GBM Thickening Status", 
               style = "heading 2") %>%
  body_add_flextable(table4 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit()) %>%
  body_add_par("") %>%
  body_add_break() %>%
  
  # Table 5
  body_add_par("Table 5. Participant Characteristics Stratified by Arteriosclerosis Status", 
               style = "heading 2") %>%
  body_add_flextable(table5 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit())

# Save combined document
print(combined_doc, target = paste0(base_path, "All_Demographics_Tables.docx"))
print("Combined demographics tables saved!")

print("All demographics tables complete!")






















########################################################################
# ROCKIES Publication - Demographics Tables for Word Documents
########################################################################

library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

########################################################################
# DATA PREPARATION (from your existing code)
########################################################################

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

# Load clinical data for GBM/Arteriosclerosis tables
dat_clinical <- data.table::fread("C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
dat_clinical <- dat_clinical[-str_which(dat_clinical$record_id, '-O')]

########################################################################
# HELPER FUNCTION: Create publication-ready table
########################################################################

create_demo_table <- function(data, stratify_var, title, 
                              vars_to_include = NULL,
                              var_labels = NULL) {
  
  # Default variables if not specified
  if (is.null(vars_to_include)) {
    vars_to_include <- c("age", "sex", "race", "ethnicity", "bmi", 
                         "hba1c", "egfr_ckd_epi", "acr_u", 
                         "diabetes_duration", "sbp", "dbp")
  }
  
  # Filter to available variables
  available_vars <- vars_to_include[vars_to_include %in% names(data)]
  
  # Default labels
  default_labels <- list(
    age = "Age (years)",
    sex = "Sex",
    race = "Race",
    ethnicity = "Ethnicity", 
    bmi = "BMI (kg/m²)",
    hba1c = "HbA1c (%)",
    egfr_ckd_epi = "eGFR (mL/min/1.73m²)",
    acr_u = "UACR (mg/g)",
    diabetes_duration = "Diabetes Duration (years)",
    sbp = "Systolic BP (mmHg)",
    dbp = "Diastolic BP (mmHg)",
    gfr_bsa_plasma = "mGFR (mL/min/1.73m²)",
    avg_c_k2 = "Cortical K2",
    avg_c_f = "Cortical F",
    avg_c_k2_f = "Cortical K2/F"
  )
  
  if (!is.null(var_labels)) {
    default_labels <- modifyList(default_labels, var_labels)
  }
  
  # Create the table
  tbl <- data %>%
    dplyr::select(all_of(c(stratify_var, available_vars))) %>%
    tbl_summary(
      by = all_of(stratify_var),
      missing = "ifany",
      missing_text = "Missing",
      label = default_labels[available_vars],
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = list(
        all_continuous() ~ 1,
        all_categorical() ~ c(0, 1)
      )
    ) %>%
    add_p(
      test = list(
        all_continuous() ~ "wilcox.test",
        all_categorical() ~ "fisher.test"
      )
    ) %>%
    add_overall() %>%
    modify_header(label = "**Characteristic**") %>%
    modify_spanning_header(c("stat_1", "stat_2") ~ paste0("**", stratify_var, "**")) %>%
    bold_labels()
  
  return(tbl)
}

########################################################################
# TABLE 1: scRNAseq Cohort (32 participants with kidney biopsy)
########################################################################

# Filter for scRNAseq cohort - adjust filter criteria based on your data
# Assuming you have a variable indicating scRNAseq participants
dat_scrnaseq <- dat %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes')) %>%
  # Add any additional filters for scRNAseq cohort
  mutate(
    T2D_Status = ifelse(group == "Type 2 Diabetes", "T2D", "Control")
  )

# Create Table 1
table1 <- dat_scrnaseq %>%
  create_demo_table(
    stratify_var = "T2D_Status",
    title = "Table 1. Participant Characteristics for Single-Cell RNA Sequencing Cohort",
    vars_to_include = c("age", "sex", "race", "ethnicity", "bmi", 
                        "hba1c", "egfr_ckd_epi", "acr_u", "diabetes_duration")
  )

########################################################################
# TABLE 2: PET/CT Cohorts (18 T2D + 11 HC)
########################################################################

# Use dat2 from your existing code (after SGLT2i filtering)
# Recreate dat2 or use the filtered version
PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2, lc_f, rc_f, lm_f, rm_f)
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

dat_results <- dat %>% 
  dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f'))) %>%
  bind_cols(PET_avg(.)) %>%
  filter(!is.na(avg_c_k2)) %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

dat_pet <- dat_results %>%
  mutate(
    Cohort = case_when(
      group == "Type 2 Diabetes" ~ "T2D (RENAL-HEIRitage)",
      group == "Lean Control" ~ "Healthy Control (CROCODILE)"
    )
  )

table2 <- dat_pet %>%
  create_demo_table(
    stratify_var = "Cohort",
    title = "Table 2. Participant Characteristics for PET/CT Cohorts",
    vars_to_include = c("age", "sex", "race", "ethnicity", "bmi", 
                        "hba1c", "egfr_ckd_epi", "acr_u", "sbp", "dbp",
                        "avg_c_k2", "avg_c_f", "avg_c_k2_f")
  )

########################################################################
# TABLE 3: UACR-PET Analysis Cohort (40 participants)
########################################################################

dat_uacr <- dat %>%
  dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f'))) %>%
  bind_cols(PET_avg(.)) %>%
  filter(!is.na(avg_c_k2)) %>%
  filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes')) %>%
  filter(record_id != 'CRC-55') %>%
  mutate(
    Group = case_when(
      group == "Type 2 Diabetes" ~ "T2D",
      group == "Lean Control" ~ "Lean Control",
      group == "Obese Control" ~ "Obese Control"
    )
  )

table3 <- dat_uacr %>%
  create_demo_table(
    stratify_var = "Group",
    title = "Table 3. Participant Characteristics for UACR-PET Analysis Cohort",
    vars_to_include = c("age", "sex", "race", "ethnicity", "bmi", 
                        "hba1c", "egfr_ckd_epi", "acr_u",
                        "avg_c_k2", "avg_c_f", "avg_c_k2_f")
  )

########################################################################
# TABLE 4: GBM Thickening Status
########################################################################

dat_gbm <- dat_clinical %>%
  filter(`GBM thickness` != '' & !is.na(`GBM thickness`)) %>%
  mutate(
    GBM_Status = ifelse(`GBM thickness` == "yes", "GBM Thickening", "No GBM Thickening")
  )

# Merge with full dataset for demographics
dat_gbm_full <- dat_gbm %>%
  left_join(dat %>% dplyr::select(record_id, age, sex, race, ethnicity, bmi, 
                                  hba1c, egfr_ckd_epi, acr_u, diabetes_duration),
            by = "record_id")

table4 <- dat_gbm_full %>%
  create_demo_table(
    stratify_var = "GBM_Status",
    title = "Table 4. Participant Characteristics Stratified by GBM Thickening Status",
    vars_to_include = c("age", "sex", "race", "bmi", "hba1c", 
                        "egfr_ckd_epi", "acr_u", "avg_c_k2", "avg_c_f", "avg_c_k2_f")
  )

########################################################################
# TABLE 5: Arteriosclerosis Status
########################################################################

dat_arterio <- dat_clinical %>%
  filter(arteriosclerosis != '' & !is.na(arteriosclerosis)) %>%
  mutate(
    Arterio_Status = ifelse(arteriosclerosis == "yes", "Arteriosclerosis", "No Arteriosclerosis")
  )

dat_arterio_full <- dat_arterio %>%
  left_join(dat %>% dplyr::select(record_id, age, sex, race, ethnicity, bmi, 
                                  hba1c, egfr_ckd_epi, acr_u, diabetes_duration),
            by = "record_id")

table5 <- dat_arterio_full %>%
  create_demo_table(
    stratify_var = "Arterio_Status",
    title = "Table 5. Participant Characteristics Stratified by Arteriosclerosis Status",
    vars_to_include = c("age", "sex", "race", "bmi", "hba1c", 
                        "egfr_ckd_epi", "acr_u", "avg_c_k2", "avg_c_f", "avg_c_k2_f")
  )

########################################################################
# EXPORT TO WORD DOCUMENTS
########################################################################

# Function to save table as Word document
save_table_to_word <- function(tbl, filename, title) {
  
  # Convert to flextable
  ft <- tbl %>%
    as_flex_table() %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    autofit() %>%
    set_caption(caption = title)
  
  # Create Word document
  doc <- read_docx() %>%
    body_add_par(title, style = "heading 1") %>%
    body_add_flextable(ft)
  
  # Save
  print(doc, target = paste0(base_path, filename))
  print(paste0(filename, " saved!"))
}

# Save individual tables
save_table_to_word(table1, "Table1_scRNAseq_Demographics.docx", 
                   "Table 1. Participant Characteristics for Single-Cell RNA Sequencing Cohort")

save_table_to_word(table2, "Table2_PET_CT_Demographics.docx",
                   "Table 2. Participant Characteristics for PET/CT Cohorts")

save_table_to_word(table3, "Table3_UACR_PET_Demographics.docx",
                   "Table 3. Participant Characteristics for UACR-PET Analysis Cohort")

save_table_to_word(table4, "Table4_GBM_Demographics.docx",
                   "Table 4. Participant Characteristics Stratified by GBM Thickening Status")

save_table_to_word(table5, "Table5_Arteriosclerosis_Demographics.docx",
                   "Table 5. Participant Characteristics Stratified by Arteriosclerosis Status")

########################################################################
# COMBINED WORD DOCUMENT WITH ALL TABLES
########################################################################

# Create combined document
combined_doc <- read_docx() %>%
  body_add_par("ROCKIES Publication - Demographics Tables", style = "heading 1") %>%
  body_add_par("") %>%
  
  # Table 1
  body_add_par("Table 1. Participant Characteristics for Single-Cell RNA Sequencing Cohort", 
               style = "heading 2") %>%
  body_add_flextable(table1 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit()) %>%
  body_add_par("") %>%
  body_add_break() %>%
  
  # Table 2
  body_add_par("Table 2. Participant Characteristics for PET/CT Cohorts", 
               style = "heading 2") %>%
  body_add_flextable(table2 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit()) %>%
  body_add_par("") %>%
  body_add_break() %>%
  
  # Table 3
  body_add_par("Table 3. Participant Characteristics for UACR-PET Analysis Cohort", 
               style = "heading 2") %>%
  body_add_flextable(table3 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit()) %>%
  body_add_par("") %>%
  body_add_break() %>%
  
  # Table 4
  body_add_par("Table 4. Participant Characteristics Stratified by GBM Thickening Status", 
               style = "heading 2") %>%
  body_add_flextable(table4 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit()) %>%
  body_add_par("") %>%
  body_add_break() %>%
  
  # Table 5
  body_add_par("Table 5. Participant Characteristics Stratified by Arteriosclerosis Status", 
               style = "heading 2") %>%
  body_add_flextable(table5 %>% as_flex_table() %>% fontsize(size = 9, part = "all") %>% 
                       font(fontname = "Times New Roman", part = "all") %>% autofit())

# Save combined document
print(combined_doc, target = paste0(base_path, "All_Demographics_Tables.docx"))
print("Combined demographics tables saved!")

print("All demographics tables complete!")

