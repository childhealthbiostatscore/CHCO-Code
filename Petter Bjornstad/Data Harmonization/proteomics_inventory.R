library(tidyverse)

user <- Sys.info()[["user"]]

root_path <- "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW"



harm_dat <- read.csv(
  file.path(root_path, "Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_olink_harmonized_dataset.csv"),
  na.strings = ""
)
## Pull in proteomics


proteomics <- harm_dat %>% 
  dplyr::select(record_id, mrn, visit, group, study, age, bmi, seq.10000.28, seq.10000.28_urine) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(!is.na(seq.10000.28) | !is.na(seq.10000.28_urine))

# RH/RH2
#-----RH/RH2 T2D-----

rh_rh2_proteomics_t2d <- proteomics %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage")) %>%
  filter(group == "Type 2 Diabetes") %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("RH/RH2 - T2D:")
nrow(rh_rh2_proteomics_t2d)

visit_counts <- rh_rh2_proteomics_t2d %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- rh_rh2_proteomics_t2d %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- rh_rh2_proteomics_t2d %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")


cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")
#-----RH/RH2 LEAN-----

rh_rh2_proteomics_l <- proteomics %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage")) %>%
  filter(group == "Lean Control") %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("RH/RH2 - Lean")
nrow(rh_rh2_proteomics_l)

visit_counts <- rh_rh2_proteomics_l %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- rh_rh2_proteomics_l %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- rh_rh2_proteomics_l %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")# N = 0
#-----RH/RH2 OBESE-----

rh_rh2_proteomics_o <- proteomics %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage")) %>%
  filter(group == "Obese Control") %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("RH/RH2 - Obese")
nrow(rh_rh2_proteomics_o)

visit_counts <- rh_rh2_proteomics_o %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- rh_rh2_proteomics_o %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- rh_rh2_proteomics_o %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")# N = 0
# TODAY
if (user == "shivaniramesh") {
  root_path <- "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/"
} 

load(file.path(root_path, "/TODAY subaward/Somalogic data raw/soma.Rdata"))
nrow(soma)
# N = 644
## Unique
today_unique <- soma %>% distinct(releaseid, .keep_all = T)
nrow(today_unique)
# N = 379

# TEEN-LABS

load(file.path(root_path, "Teen Labs/Data_Cleaned/analysis_dataset.RData"))
tl_proteomics <- df %>%
  filter(!is.na(seq.10000.28))
nrow(tl_proteomics)
# N = 324
## Unique
tl_unique <- tl_proteomics %>% distinct(ID, .keep_all = T)
nrow(tl_unique)
# N = 64

# RH/RH2/TODAY/TEEN-LABS: N = 644+61+324 = 1029

# CASPER, CROCODILE, ATTEMPT
#-----CASPER T1D-----

casper_proteomics <- proteomics %>%
  filter(study %in% c("CASPER")) %>%
  filter(group == "Type 1 Diabetes" | is.na(group)) %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("CASPER - T1D")
nrow(casper_proteomics)

visit_counts <- casper_proteomics %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- casper_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- casper_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")# N = 0
#-----IMPROVE T2D-----

improve_proteomics <- proteomics %>%
  filter(study %in% c("IMPROVE")) %>%
  filter(group == "Type 2 Diabetes" | is.na(group)) %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("IMPROVE - T2D")
nrow(improve_proteomics)

visit_counts <- improve_proteomics %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- improve_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- improve_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

#Participants that have missing visit 1 or more samples (printing raw record_id, visit values), some studies have multiple visits per participant
cat("Participants that have missing visit 1 or more samples:\n", 
    paste0(improve_proteomics$record_id[improve_proteomics$record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]], 
           " - ", improve_proteomics$visit[improve_proteomics$record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]], ")", collapse = ", \n"), "\n")

cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")# N = 30 
#-----CROC T1D-----
croc_proteomics <- proteomics %>%
  filter(study %in% c("CROCODILE")) %>%
  filter(group == "Type 1 Diabetes" | is.na(group)) %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("CROCODILE - T1D")
nrow(croc_proteomics)

visit_counts <- croc_proteomics %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- croc_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- croc_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")# N = 30 for CROC
#-----CROC LEAN-----
croc_proteomics_l <- proteomics %>%
  filter(study %in% c("CROCODILE")) %>%
  filter(group == "Lean Control" | is.na(group)) %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("CROCODILE - Lean")
nrow(croc_proteomics_l)

visit_counts <- croc_proteomics_l %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- croc_proteomics_l %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- croc_proteomics_l %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")
# N = 18 for CROC
#-----ATTEMPT BLOOD-----
attempt_blood_proteomics <- proteomics %>%
  filter(study %in% c("ATTEMPT")) %>%
  filter(group == "Type 1 Diabetes" | is.na(group)) %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("ATTEMPT (Blood) - T1D")
nrow(attempt_blood_proteomics)

visit_counts <- attempt_blood_proteomics %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- attempt_blood_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- attempt_blood_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

cat("Participants that have missing visit 1 or more samples:\n", 
    paste0(attempt_blood_proteomics$record_id[attempt_blood_proteomics$record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]], 
           " - ", attempt_blood_proteomics$visit[attempt_blood_proteomics$record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]], ")", collapse = ", \n"), "\n")


cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")
# N = 188

#-----ATTEMPT URINE-----

attempt_urine_proteomics <- proteomics %>%
  filter(study %in% c("ATTEMPT")) %>%
  filter(group == "Type 1 Diabetes" | is.na(group)) %>%
  distinct(mrn, seq.10000.28_urine, .keep_all = T)
cat("ATTEMPT (urine) - T1D")
nrow(attempt_urine_proteomics)

visit_counts <- attempt_urine_proteomics %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- attempt_urine_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- attempt_urine_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

cat("Participants that have missing visit 1 or more samples:\n", 
    paste0(attempt_urine_proteomics$record_id[attempt_urine_proteomics$record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]], 
           " - ", attempt_urine_proteomics$visit[attempt_urine_proteomics$record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]], ")", collapse = ", \n"), "\n")


cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")

#-----PENGUIN-----

penguin_proteomics <- proteomics %>%
  filter(study %in% c("PENGUIN")) %>%
  filter(group == "PKD" | is.na(group)) %>%
  distinct(mrn, seq.10000.28, .keep_all = T)
cat("PENGUIN - PKD")
nrow(penguin_proteomics)

visit_counts <- penguin_proteomics %>%
  group_by(record_id) %>%
  summarise(n_visits = n(), .groups = "drop")

# Summaries
n_total <- n_distinct(visit_counts$record_id)
n_one_visit <- sum(visit_counts$n_visits == 1)
n_multiple_visits <- sum(visit_counts$n_visits > 1)

samples_one_visit <- penguin_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits == 1]) %>%
  nrow()

samples_multiple_visits <- penguin_proteomics %>%
  filter(record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]) %>%
  nrow()

cat("Samples in 1-visit participants:", samples_one_visit, "\n")
cat("Samples in >1-visit participants:", samples_multiple_visits, "\n")

#Participants that have missing visit 1 or more samples (printing raw record_id, visit values), some studies have multiple visits per participant
cat("Participants that have missing visit 1 or more samples:\n", 
    paste0(penguin_proteomics$record_id[penguin_proteomics$record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]], 
           " - ", penguin_proteomics$visit[penguin_proteomics$record_id %in% visit_counts$record_id[visit_counts$n_visits > 1]], collapse = ", \n"), "\n")

cat("Total participants:", n_total, "\n")
cat("Participants with only 1 visit:", n_one_visit, "\n")
cat("Participants with >1 visit:", n_multiple_visits, "\n")# N = 30 

#Table: Rows: record_id (from all studies), study, group; Columns: all unique visit values across studies; Fill in existing proteomics data with blood or urine (urine on in the sections above marked as usrine, assume all others are blood)
#all_proteomics <- proteomics %>%
#  mutate(proteomics_type = case_when(!is.na(seq.10000.28) ~ "Blood",
#                                    !is.na(seq.10000.28_urine) ~ "Urine",
#                                    TRUE ~ NA_character_)) %>%
#  select(record_id, study, group, visit, proteomics_type) %>%
#  filter(!is.na(proteomics_type)) %>%
#  distinct(record_id, study, group, visit, proteomics_type, .keep_all = T) %>%
#  pivot_wider(names_from = visit, values_from = proteomics_type, values_fill = NA)


  

