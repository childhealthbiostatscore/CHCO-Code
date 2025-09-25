# ---
# CKM Outcomes Wrangling Script
# Author: Darwin Del Castillo
# Date: `r lubridate::today()`
# ---

# Importing dataset
CKM_outcomes <- read_csv("./CKM_outcomes/CKM stage summary.csv")

# Fixing variable names
CKM_outcomes <- CKM_outcomes |> janitor::clean_names()

# Transforming CKM stages to factors
cleaned_log <- cleaned_log |> 
  mutate(ckm_syn_base_factor = case_when(ckm_syn_base == "Stage 2" ~ 0,
                                         ckm_syn_base == "Stage 2+" ~ 1,
                                         ckm_syn_base == "Stage 3" ~ 2,
                                         .default	= NA_real_))
#only ckm stages 2, 2+ and 3
ckm_stages <- cleaned_log$ckm_syn_base_factor
ckm_stages <- cbind(rep(1,nrow(cleaned_log)),ckm_stages)