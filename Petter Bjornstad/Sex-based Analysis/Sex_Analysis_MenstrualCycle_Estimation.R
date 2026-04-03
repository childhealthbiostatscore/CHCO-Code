####Sex Differences. Estimating menstrual cycle 






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





harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))




# Assumptions
cycle_length <- 28
follicular_days <- 1:13
ovulation_days <- 14
luteal_days <- 15:28
menstruation_days <- 1:5

# Function to classify cycle day into phase
get_phase <- function(cycle_day) {
  day <- ((cycle_day - 1) %% cycle_length) + 1
  case_when(
    day %in% menstruation_days ~ "menstruation",
    day %in% setdiff(follicular_days, menstruation_days) ~ "follicular",
    day == ovulation_days ~ "ovulation",
    day %in% luteal_days ~ "luteal",
    TRUE ~ NA_character_
  )
}

# Function to calculate phase probabilities for a given day difference
calculate_phase_probabilities <- function(days_difference) {
  possible_clamp_days <- follicular_days
  n_possibilities <- length(possible_clamp_days)
  
  biopsy_days <- possible_clamp_days + days_difference
  biopsy_phases <- sapply(biopsy_days, get_phase)
  
  phase_counts <- table(biopsy_phases)
  phase_probs <- phase_counts / n_possibilities
  
  all_phases <- c("menstruation", "follicular", "ovulation", "luteal")
  result <- setNames(rep(0, 4), all_phases)
  result[names(phase_probs)] <- as.numeric(phase_probs)
  
  return(as.list(result))
}

# ============================================================
# Extract relevant procedures from harmonized_data
# Filter to women aged 10-60
# ============================================================

# First filter to eligible participants
eligible_participants <- harmonized_data %>%
  filter(sex == "Female" | sex == "F",
         age >= 10 & age <= 60) %>%
  dplyr::select(mrn) %>%
  distinct()

# Get clamp dates (reference point - known follicular phase)
clamp_dates <- harmonized_data %>%
  semi_join(eligible_participants, by = "mrn") %>%
  filter(procedure %in% c("clamp", "croc_clamp")) %>%
  dplyr::select(mrn, clamp_date = screen_date, visit_clamp = visit) %>%
  group_by(mrn) %>%
  slice_min(clamp_date, n = 1, with_ties = FALSE) %>%
  ungroup()

# Get kidney biopsy dates
biopsy_dates <- harmonized_data %>%
  semi_join(eligible_participants, by = "mrn") %>%
  filter(procedure == "kidney_biopsy") %>%
  dplyr::select(mrn, biopsy_date = screen_date, visit_biopsy = visit)

# ============================================================
# Join and calculate phase probabilities
# ============================================================

analysis_data <- biopsy_dates %>%
  inner_join(clamp_dates, by = "mrn") %>%
  mutate(
    days_difference = as.numeric(difftime(biopsy_date, clamp_date, units = "days"))
  )

# Check the data
cat("Number of participants with both procedures:", nrow(analysis_data), "\n")
cat("Days difference range:", range(analysis_data$days_difference), "\n\n")

# Calculate phase probabilities for each participant
phase_results <- analysis_data %>%
  rowwise() %>%
  mutate(probs = list(calculate_phase_probabilities(days_difference))) %>%
  unnest_wider(probs) %>%
  ungroup()

# Add most likely phase
phase_results <- phase_results %>%
  mutate(
    most_likely_phase = case_when(
      menstruation == pmax(menstruation, follicular, ovulation, luteal) ~ "menstruation",
      follicular == pmax(menstruation, follicular, ovulation, luteal) ~ "follicular",
      ovulation == pmax(menstruation, follicular, ovulation, luteal) ~ "ovulation",
      TRUE ~ "luteal"
    ),
    max_probability = pmax(menstruation, follicular, ovulation, luteal)
  )

# ============================================================
# View results
# ============================================================

print(phase_results %>% 
        dplyr::select(mrn, clamp_date, biopsy_date, days_difference,
                      menstruation, follicular, ovulation, luteal, 
                      most_likely_phase, max_probability))

# Summary
cat("\n--- Summary across all participants ---\n")
phase_results %>%
  summarise(
    n = n(),
    avg_days_diff = mean(days_difference),
    across(c(menstruation, follicular, ovulation, luteal), mean, .names = "avg_{.col}")
  ) %>%
  print()

# Distribution of most likely phases
cat("\n--- Most likely phase distribution ---\n")
table(phase_results$most_likely_phase) %>% print()