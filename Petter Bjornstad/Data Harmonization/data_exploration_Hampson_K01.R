#Cohorts for K01 
library(tidyverse)
# SUMMARY ----
#Read in harmonized dataset
dir.dat <- c("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Data Harmonization/Data Clean")
meta <- read.csv(fs::path(dir.dat,"soma_harmonized_dataset.csv"))
# meta_PANTHER <- meta %>% 
#   filter(study=="PANTHER")
#Collaps 
meta_collaps <- meta %>%
  filter(study=="RENAL-HEIR" | study=="IMPROVE" | study=="RENAL-HEIRitage" | study=="PANTHER" | study=="CROCODILE") %>% 
  filter(visit=="baseline") %>%
  arrange(screen_date) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(mrn, visit)) 
length(unique(meta_collaps$mrn))#306
length(unique(meta_collaps$record_id))#306


#Filter to studies of interest: PANTHER, RH/RH2, CROCODILE, IMPROVE
meta_collaps <- meta_collaps %>% 
  filter(study=="RENAL-HEIR" | study=="IMPROVE" | study=="RENAL-HEIRitage" | study=="PANTHER" | study=="CROCODILE")
length(unique(meta_collaps$mrn))#307
length(unique(meta_collaps$record_id))#353

#"mrn","record_id","procedure",
prot.names <- c(colnames(meta_collaps)[which(grepl("seq.",colnames(meta_collaps)) & !grepl("_urine_cradj",colnames(meta_collaps)) & !grepl("_urine",colnames(meta_collaps)))])
met.names <- c(colnames(meta_collaps)[which(startsWith(colnames(meta_collaps),"c0") | startsWith(colnames(meta_collaps),"c1"))])

meta_omics <- meta_collaps %>% 
  dplyr::select(all_of(c("mrn","record_id","procedure","study","age","sex","bmi","race","ethnicity","group",met.names,prot.names)))

first_all_na <- which(sapply(meta_omics, function(x) all(is.na(x))))[1] #7484
# If such a column exists, keep only columns before it
meta_omics <- meta_omics %>%
  dplyr::select(1:(first_all_na - 1))

#Remake after filtering out missing proteins
prot.names <- c(colnames(meta_omics)[which(grepl("seq.",colnames(meta_omics)) & !grepl("_urine_cradj",colnames(meta_omics)) & !grepl("_urine",colnames(meta_omics)))])
met.names <- c(colnames(meta_omics)[which(startsWith(colnames(meta_omics),"c0") | startsWith(colnames(meta_omics),"c1"))])


meta_summary <- meta_omics %>% 
  group_by(study) %>% 
  summarize(
    # Sample size
    n_samples = n(),
    
    # Missing value counts for omics data
    n_missing_prot = sum(rowSums(is.na(select(cur_data(), all_of(prot.names))))),
    n_missing_met = sum(rowSums(is.na(select(cur_data(), all_of(met.names))))),
    
    # Proportion of missing values
    prop_missing_prot = n_missing_prot / (n() * length(prot.names)),
    prop_missing_met = n_missing_met / (n() * length(met.names)),
    
    # Age summary
    mean_age = mean(age, na.rm = TRUE),
    sd_age = sd(age, na.rm = TRUE),
    median_age = median(age, na.rm = TRUE),
    n_missing_age = sum(is.na(age)),
    
    # Sex distribution (assuming binary coding or factor)
    n_female = sum(sex == "Female" | sex == "F" | sex == 1, na.rm = TRUE),
    n_male = sum(sex == "Male" | sex == "M" | sex == 0, na.rm = TRUE),
    prop_female = n_female / (n_female + n_male),
    n_missing_sex = sum(is.na(sex)),
    
    # Race/ethnicity (most common categories)
    n_missing_race = sum(is.na(race)),
    n_missing_ethnicity = sum(is.na(ethnicity)),
    
    # BMI summary
    mean_bmi = mean(bmi, na.rm = TRUE),
    sd_bmi = sd(bmi, na.rm = TRUE),
    median_bmi = median(bmi, na.rm = TRUE),
    n_missing_bmi = sum(is.na(bmi)),
    
    # Group information (assuming this is a treatment/control variable)
    n_missing_group = sum(is.na(group)),
    
    .groups = "drop"
  )

# If you want to see the most common race/ethnicity categories by study:
race_ethnicity_summary <- meta_omics %>%
  group_by(study) %>%
  summarize(
    most_common_race = names(sort(table(race), decreasing = TRUE))[1],
    most_common_ethnicity = names(sort(table(ethnicity), decreasing = TRUE))[1],
    race_categories = paste(names(table(race)), collapse = ", "),
    ethnicity_categories = paste(names(table(ethnicity)), collapse = ", "),
    .groups = "drop"
  )

# Combine the summaries
meta_summary_complete <- meta_summary %>%
  left_join(race_ethnicity_summary, by = "study")
# write.csv(meta_summary_complete,"/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/4_Grant Applications/K01 Grant/Summary_of_Cohorts.csv")
# meta_summary_complete <- meta_summary_complete %>% 
#   dplyr::select(-c("n_missing_prot","n_missing_met"))


# CROCODILE ----
CROCODILE <- meta_omics %>% 
  filter(study=="CROCODILE")
CROC_N <- length(unique(CROCODILE$mrn))
summary(CROCODILE$age)  #18.75 - 30.15
length(which(CROCODILE$group=="Type 1 Diabetes")) #39 type 1
length(which(CROCODILE$group=="Lean Control")) #21 lean controls
summary(CROCODILE$bmi) #19.29 - 31.93 
# Calculate participants with at least one non-missing metabolomics value
CROCODILE_metabolomics <- CROCODILE %>%
  select(all_of(met.names)) %>%
  rowwise() %>%
  mutate(has_metabolomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_metabolomics = sum(has_metabolomics)) %>%
  pull(has_metabolomics)
CROCODILE_metabolomics #50 have metabolomics
CROCODILE_missing_Metabolomics <- CROC_N - CROCODILE_metabolomics
CROCODILE_missing_Metabolomics #10 missing metabolomics

# Calculate participants with at least one non-missing proteomics value
CROCODILE_proteomics <- CROCODILE %>%
  select(all_of(prot.names)) %>%
  rowwise() %>%
  mutate(has_proteomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_proteomics = sum(has_proteomics)) %>%
  pull(has_proteomics)
CROCODILE_proteomics #48 with proteomics
CROCODILE_missing_Proteomics <- CROC_N - CROCODILE_proteomics
CROCODILE_missing_Proteomics #12 missing metabolomics
  
# IMPROVE ----
IMPROVE <- meta_omics %>% 
  filter(study=="IMPROVE")
summary(IMPROVE$age) 
# Calculate participants with at least one non-missing metabolomics value
IMPROVE_metabolomics <- IMPROVE %>%
  select(all_of(met.names)) %>%
  rowwise() %>%
  mutate(has_metabolomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_metabolomics = sum(has_metabolomics)) %>%
  pull(has_metabolomics)
IMPROVE_metabolomics #0 have metabolomics, 11 missing metabolomics

# Calculate participants with at least one non-missing proteomics value
IMPROVE_proteomics <- IMPROVE %>%
  select(all_of(prot.names)) %>%
  rowwise() %>%
  mutate(has_proteomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_proteomics = sum(has_proteomics)) %>%
  pull(has_proteomics)
IMPROVE_proteomics #7 with proteomics, 4 missing proteomics

# PANTHER ----
PANTHER <- meta_omics %>% 
  filter(study=="PANTHER")
PANTHER_N <- length(unique(PANTHER$mrn))
summary(PANTHER$age)  #18.75 - 30.15
unique(PANTHER$group) #"Type 2 Diabetes" "Obese Control"   "Lean Control"
length(which(PANTHER$group=="Type 1 Diabetes")) #39 type 1
length(which(PANTHER$group=="Lean Control")) #21 lean controls
summary(PANTHER$bmi) #19.29 - 31.93 
# Calculate participants with at least one non-missing metabolomics value
PANTHER_metabolomics <- PANTHER %>%
  select(all_of(met.names)) %>%
  rowwise() %>%
  mutate(has_metabolomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_metabolomics = sum(has_metabolomics)) %>%
  pull(has_metabolomics)
PANTHER_metabolomics #0 have metabolomics
PANTHER_missing_Metabolomics <- PANTHER_N - PANTHER_metabolomics
PANTHER_missing_Metabolomics #117 missing metabolomics

# Calculate participants with at least one non-missing proteomics value
PANTHER_proteomics <- PANTHER %>%
  select(all_of(prot.names)) %>%
  rowwise() %>%
  mutate(has_proteomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_proteomics = sum(has_proteomics)) %>%
  pull(has_proteomics)
PANTHER_proteomics #89 with proteomics
PANTHER_missing_Proteomics <- PANTHER_N - PANTHER_proteomics
PANTHER_missing_Proteomics #28 missing metabolomics

# RENAL-HEIR ----
RH <- meta_omics %>% 
  filter(study=="RENAL-HEIR")
RH_N <- length(unique(RH$mrn))
summary(RH$age)  #18.75 - 30.15
unique(RH$group) #"Type 2 Diabetes" "Obese Control"   "Lean Control"
length(which(RH$group=="Type 2 Diabetes")) #40 type 2
length(which(RH$group=="Lean Control")) #11 lean controls
length(which(RH$group=="Obese Control")) #10 lean controls
summary(RH$bmi) #19.29 - 31.93 
# Calculate participants with at least one non-missing metabolomics value
RH_metabolomics <- RH %>%
  select(all_of(met.names)) %>%
  rowwise() %>%
  mutate(has_metabolomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_metabolomics = sum(has_metabolomics)) %>%
  pull(has_metabolomics)
RH_metabolomics #46 have metabolomics
RH_missing_Metabolomics <- RH_N - RH_metabolomics
RH_missing_Metabolomics #15 missing metabolomics

# Calculate participants with at least one non-missing proteomics value
RH_proteomics <- RH %>%
  select(all_of(prot.names)) %>%
  rowwise() %>%
  mutate(has_proteomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_proteomics = sum(has_proteomics)) %>%
  pull(has_proteomics)
RH_proteomics #48 with proteomics
RH_missing_Proteomics <- RH_N - RH_proteomics
RH_missing_Proteomics #13 missing metabolomics

# RENAL-HEIRitage ----
RH2 <- meta_omics %>% 
  filter(study=="RENAL-HEIRitage")
RH2_N <- length(unique(RH2$mrn))
summary(RH2$age)  # 11.59  -   76.19
unique(RH2$group) #"Type 2 Diabetes" "Obese Control"   "Lean Control"
length(which(RH2$group=="Type 2 Diabetes")) #33 type 2
length(which(RH2$group=="Lean Control")) #9 lean controls
length(which(RH2$group=="Obese Control")) #15 lean controls
summary(RH2$bmi) #17.54- 52.79 
# Calculate participants with at least one non-missing metabolomics value
RH2_metabolomics <- RH2 %>%
  select(all_of(met.names)) %>%
  rowwise() %>%
  mutate(has_metabolomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_metabolomics = sum(has_metabolomics)) %>%
  pull(has_metabolomics)
RH2_metabolomics #46 have metabolomics
RH2_missing_Metabolomics <- RH2_N - RH2_metabolomics
RH2_missing_Metabolomics #11 missing metabolomics

# Calculate participants with at least one non-missing proteomics value
RH2_proteomics <- RH2 %>%
  select(all_of(prot.names)) %>%
  rowwise() %>%
  mutate(has_proteomics = any(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  summarise(has_proteomics = sum(has_proteomics)) %>%
  pull(has_proteomics)
RH2_proteomics #48 with proteomics
RH2_missing_Proteomics <- RH2_N - RH2_proteomics
RH2_missing_Proteomics #9 missing metabolomics
mean(RH2$age)
hist(RH2$age)
