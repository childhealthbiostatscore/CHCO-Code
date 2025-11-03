library(dplyr)
library(tidyr)
library(arsenal)
library(purrr)
library(Hmisc)
library(growthcleanr)
library(stringr)

dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

# withdrew/LTFU
exclude <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/PANTHER/Data_Cleaned/panther_withdrew_ltfu_list.csv")$record_id

dat <- dat %>%
  filter(study == "PANTHER") %>%
  mutate(visit = case_when(visit %in% c("baseline", "screening") ~ "baseline")) %>%
  dplyr::summarise(dplyr::across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   dplyr::across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  fill(hba1c) %>% fill(acr_u) %>% ungroup() %>%
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other"),
                race = case_when(race == "Black/African American & White" ~ "More Than One", 
                                 race == "Hawaiian or Pacific Islander & White" ~ "More Than One",
                                 T ~ race),
                status = case_when(record_id %in% exclude ~ "Removed",
                                   T ~ "Currently enrolled"),
                age = case_when(record_id == "PAN-18-O" ~ 10, T ~ age), # 10 at 2025, manually editing because no date of consent?
                tanner_stage_comp = as.character(coalesce(tan_fgd, tan_fph, tan_tveq, tan_mgd, tan_mph, breast_tanner)),
                age_mo = (age * 12),
                sex = case_when(sex == "Male" ~ "male",
                                sex == "Female" ~ "female"),
                male_ind = case_when(sex == "male" ~ 1, sex == "female" ~ 0)) %>% 
  filter(visit == "baseline")
  
bmi_percentile = ext_bmiz(data = subset(dat, 
                                        select = c("record_id", "sex", "age_mo", "weight", "height", "bmi")), 
                          age = "age_mo", 
                          wt = "weight", 
                          ht = "height", 
                          bmi = "bmi", 
                          adjust.integer.age = F) %>% 
  dplyr:: select(record_id, bmip, bmiz) %>%
  filter(!is.na(bmip))
dat <- left_join(dat, bmi_percentile, by = "record_id") %>%
  dplyr::mutate(sex_group_risk = case_when(sex == "male" & group_risk == "High" ~ "M, High",
                                           sex == "male" & group_risk == "Low" ~ "M, Low",
                                           sex == "female" & group_risk == "High" ~ "F, High",
                                           sex == "female" & group_risk == "Low" ~ "F, Low"))

table(subset(dat, record_id %in% exclude)$group_risk)
ier <- dat %>%
  # filter(record_id %nin% exclude) %>%
  dplyr::select(race, ethnicity, sex, age, status) %>%
  dplyr::rename(Race = race,
                Ethnicity = ethnicity,
                Gender = sex,
                Age = age,
                Status = status) %>%
  dplyr::mutate(Age = floor(Age), 
                "Age Unit" = "Years",
                Gender = str_to_sentence(Gender)) %>%
  arrange(Status)

summary(tableby(group_risk ~ age + sex + bmi  + bmip + tanner_stage_comp + race_ethnicity_condensed + total_kidney_volume_ml + mm_ir + mm_di + gfr_raw_plasma + erpf_raw_plasma , 
                data = subset(dat, record_id %nin% exclude), total = T, test = F,  numeric.simplify = T, numeric.stats = c("meansd"), digits = 1))

summary(tableby(group_risk ~ kwt(mm_di, "Nmiss", "medianq1q3", "range"), data = subset(dat,record_id %nin% exclude & mm_di < 10000), total = T, test = F))
summary(tableby(group_risk ~ kwt(mm_ir, "Nmiss", "medianq1q3", "range"), data = subset(dat,record_id %nin% exclude & mm_ir < 100), total = T, test = F))


write.csv(ier, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/PANTHER/Data_Cleaned/2025_panther_ier_2.csv", row.names = F)

summary_table <- ier %>%
  group_by(Race, Ethnicity, Gender) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = c(Ethnicity, Gender),
    values_from = count,
    values_fill = 0
  ) %>%
  # reorder columns
  select(
    Race,
    `Not Hispanic or Latino_Female`,
    `Not Hispanic or Latino_Male`,
    `Hispanic or Latino_Female`,
    `Hispanic or Latino_Male`
  ) %>%
  # reorder rows
  mutate(Race = factor(Race, levels = c(
    "American Indian or Alaskan Native",
    "Asian",
    "Hawaiian or Pacific Islander",
    "Black or African American",
    "White",
    "More Than One",
    "Unknown"
  ))) %>%
  arrange(Race)
summary_table


# Excluding LTFU
exc_ier <- dat %>%
  filter(record_id %nin% exclude) %>%
  dplyr::select(race, ethnicity, sex, age, status) %>%
  dplyr::rename(Race = race,
                Ethnicity = ethnicity,
                Gender = sex,
                Age = age,
                Status = status) %>%
  dplyr::mutate(Age = floor(Age), 
                "Age Unit" = "Years",
                Gender = str_to_sentence(Gender)) %>%
  arrange(Status)

exc_summary_table <- exc_ier %>%
  group_by(Race, Ethnicity, Gender) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = c(Ethnicity, Gender),
    values_from = count,
    values_fill = 0
  ) %>%
  # reorder columns
  select(
    Race,
    `Not Hispanic or Latino_Female`,
    `Not Hispanic or Latino_Male`,
    `Hispanic or Latino_Female`,
    `Hispanic or Latino_Male`
  ) %>%
  # reorder rows
  mutate(Race = factor(Race, levels = c(
    "American Indian or Alaskan Native",
    "Asian",
    "Hawaiian or Pacific Islander",
    "Black or African American",
    "White",
    "More Than One",
    "Unknown"
  ))) %>%
  arrange(Race)
exc_summary_table


# Projected IER based on CROC
# 40 T1D, (50% M/50% F); age & sex distribution to PANTHER, race/eth distribution to CROC

# 20 HC, (50% M/50% F); age & sex distribution to PANTHER, race/eth distribution to CROC

dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

croc_dat <- dat %>%
  filter(study == "CROCODILE") %>%
  mutate(visit = case_when(visit %in% c("baseline", "screening") ~ "baseline")) %>%
  dplyr::summarise(dplyr::across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   dplyr::across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  fill(hba1c) %>% fill(acr_u) %>% ungroup() %>%
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other"),
                race = case_when(race == "Black/African American & White" ~ "More Than One", 
                                 race == "Hawaiian or Pacific Islander & White" ~ "More Than One",
                                 T ~ race),
                status = case_when(record_id %in% exclude ~ "Removed",
                                   T ~ "Currently enrolled"),
                age = case_when(record_id == "PAN-18-O" ~ 10, T ~ age), # 10 at 2025, manually editing because no date of consent?
                tanner_stage_comp = as.character(coalesce(tan_fgd, tan_fph, tan_tveq, tan_mgd, tan_mph, breast_tanner)),
                age_mo = (age * 12),
                sex = case_when(sex == "Male" ~ "male",
                                sex == "Female" ~ "female"),
                male_ind = case_when(sex == "male" ~ 1, sex == "female" ~ 0)) %>% 
  filter(visit == "baseline")

croc_dat %>%
  # filter(record_id %nin% exclude) %>%
  dplyr::select(race, ethnicity, sex, age, status) %>%
  dplyr::rename(Race = race,
                Ethnicity = ethnicity,
                Gender = sex,
                Age = age,
                Status = status) %>%
  dplyr::mutate(Age = floor(Age), 
                "Age Unit" = "Years",
                Gender = str_to_sentence(Gender)) %>%
  arrange(Status) %>%
  group_by(Race, Ethnicity, Gender) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = c(Ethnicity, Gender),
    values_from = count,
    values_fill = 0
  ) %>%
  # reorder columns
  select(
    Race,
    `Not Hispanic or Latino_Female`,
    `Not Hispanic or Latino_Male`,
    `Hispanic or Latino_Female`,
    `Hispanic or Latino_Male`
  ) %>%
  # reorder rows
  mutate(Race = factor(Race, levels = c(
    "American Indian or Alaskan Native",
    "Asian",
    "Hawaiian or Pacific Islander",
    "Black or African American",
    "White",
    "More Than One",
    "Unknown"
  ))) %>%
  arrange(Race)

