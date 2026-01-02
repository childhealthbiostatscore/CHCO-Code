# ---
# CKM Outcomes Wrangling Script
# Author: Darwin Del Castillo
# Date: `r lubridate::today()`
# ---

# Importing dataset
CKM_outcomes <- read_csv("./ckm/CKM stage summary.csv")

# Fixing variable names
CKM_outcomes <- CKM_outcomes |> janitor::clean_names()

###########################
# Plasma Data Preparation #
###########################

# take only the plasma baseline samples
plasma$Date.Drawn <- as.Date(plasma$Date.Drawn, format = "%m/%d/%Y")
baseline_plasma <- plasma |>
  arrange(releaseid, Date.Drawn) |>
  group_by(releaseid) |>
  filter(row_number() == 1)
baseline_plasma <- baseline_plasma |> arrange(Date.Drawn)

# merging ckm outcomes with plasma baseline
cleaned_baseline_plasma_data <- merge(baseline_plasma,
                                      CKM_outcomes,
                                      by = "releaseid",
                                      all.x = TRUE,
                                      all.y = FALSE)

# Note: CKM_outcomes already includes covariates

# identify columns corresponding to metabolites (plasma: .in.uM)
is_seq_plasma <- function(.x) grepl("\\.in\\.uM$", .x)
seq_plasma <- is_seq_plasma(names(cleaned_baseline_plasma_data))

# convert to numeric
cleaned_baseline_plasma_data[, seq_plasma] <-
  apply(cleaned_baseline_plasma_data[, seq_plasma], 2, as.numeric)

# are there metabolites with low variability?
no_var <- caret::nearZeroVar(cleaned_baseline_plasma_data[, seq_plasma])
if(length(no_var) > 0){
  cleaned_baseline_plasma_data <- cleaned_baseline_plasma_data[, -no_var]
  seq_plasma <- is_seq_plasma(names(cleaned_baseline_plasma_data))
}

# log transform
cleaned_log_baseline_plasma_data <- cleaned_baseline_plasma_data |>
  mutate(across(all_of(names(cleaned_baseline_plasma_data)[seq_plasma]),
                ~log(.x)))

# prepare CKM stages for analysis as categorical factor
cleaned_log_baseline_plasma_data <- cleaned_log_baseline_plasma_data |>
  mutate(ckm_syn_base_factor = factor(ckm_syn_base,
                                      levels = c("Stage 2",
                                                 "Stage 2+",
                                                 "Stage 3")))

##########################
# Urine Data Preparation #
##########################

# take only the urine baseline samples
urine$Date.Drawn <- as.Date(urine$Date.Drawn, format = "%m/%d/%Y")
baseline_urine <- urine |>
  arrange(releaseid, Date.Drawn) |>
  group_by(releaseid) |>
  filter(row_number() == 1)
baseline_urine <- baseline_urine |> arrange(Date.Drawn)

# merging ckm outcomes with urine baseline
cleaned_baseline_urine_data <- merge(baseline_urine,
                                     CKM_outcomes,
                                     by = "releaseid",
                                     all.x = TRUE,
                                     all.y = FALSE)

# identify creatinine-adjusted metabolite columns (uM/mM.Creatinine only)
is_seq_urine <- function(.x) grepl("\\.in\\.uM/mM\\.Creatinine$", .x)
seq_urine <- is_seq_urine(names(cleaned_baseline_urine_data))

# convert to numeric
cleaned_baseline_urine_data[, seq_urine] <-
  apply(cleaned_baseline_urine_data[, seq_urine], 2, as.numeric)

# are there metabolites with low variability?
no_var_urine <- caret::nearZeroVar(cleaned_baseline_urine_data[, seq_urine])
if(length(no_var_urine) > 0){
  cleaned_baseline_urine_data <- cleaned_baseline_urine_data[, -no_var_urine]
  seq_urine <- is_seq_urine(names(cleaned_baseline_urine_data))
}

# log transform
cleaned_log_baseline_urine_data <- cleaned_baseline_urine_data |>
  mutate(across(all_of(names(cleaned_baseline_urine_data)[seq_urine]),
                ~log(.x)))

# prepare CKM stages for analysis as categorical factor
cleaned_log_baseline_urine_data <- cleaned_log_baseline_urine_data |>
  mutate(ckm_syn_base_factor = factor(ckm_syn_base,
                                      levels = c("Stage 2",
                                                 "Stage 2+",
                                                 "Stage 3")))