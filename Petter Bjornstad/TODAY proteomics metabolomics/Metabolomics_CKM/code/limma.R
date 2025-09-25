# ---
# limma baseline analysis (unadjusted and adjusted)
# Author: Darwin Del Castillo
# Date: `r lubridate::today()`
# ---

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

# Note: CKM_outcomes already includes covariates like MAC, HbA1c, log_trig, sbp, si_1_ins0

# identify columns corresponding to proteins
is_seq <- function(.x) grepl("uM", .x)
seq <- is_seq(names(cleaned_baseline_plasma_data))

# convert to numeric
cleaned_baseline_plasma_data[,seq] <- apply(cleaned_baseline_plasma_data[, seq],
                                            2, 
                                            as.numeric)

# log transform
cleaned_log_baseline_plasma_data <- cleaned_baseline_plasma_data |>
  modify_if(is_seq(names(cleaned_baseline_plasma_data)), log)

# prepare CKM stages for analysis
cleaned_log_baseline_plasma_data <- cleaned_log_baseline_plasma_data |>
  mutate(ckm_syn_base_factor = case_when(ckm_syn_base == "Stage 2" ~ 0,
                                         ckm_syn_base == "Stage 2+" ~ 1,
                                         ckm_syn_base == "Stage 3" ~ 2,
                                         .default = NA_real_))

# create design matrix for unadjusted analysis
ckm_stages <- cleaned_log_baseline_plasma_data$ckm_syn_base_factor
ckm_stages <- cbind(rep(1, nrow(cleaned_log_baseline_plasma_data)), ckm_stages)

################################
# Linear models for CKM stages #
################################

# unadjusted linear models 
## moderated t-tests with and without FDR adjustment
ymat <- t(cleaned_log_baseline_plasma_data[, seq])
fit <- lmFit(ymat, ckm_stages)
fit <- eBayes(fit)
results_ckm <- topTable(fit, coef = 2, number = nrow(ymat))
results_ckm <- results_ckm[order(results_ckm$P.Value), ]

# keeping observations with nonmissing values covariates
cleaned_baseline_plasma_data_adj <- cleaned_baseline_plasma_data |>
  filter(!is.na(hb_a1c) & !is.na(log_trig) & !is.na(sbp) & !is.na(si_1_ins0)) |>
  mutate(ckm_syn_base_factor = case_when(ckm_syn_base == "Stage 2" ~ 0,
                                         ckm_syn_base == "Stage 2+" ~ 1,
                                         ckm_syn_base == "Stage 3" ~ 2,
                                         .default = NA_real_))

# adjusted linear models
## make design matrix
design_adj <- model.matrix(~ ckm_syn_base_factor + hb_a1c + log_trig + sbp + si_1_ins0,
                           data = cleaned_baseline_plasma_data_adj)

## log transform the adjusted dataset
cleaned_log_baseline_plasma_data_adj <- cleaned_baseline_plasma_data_adj |>
  modify_if(is_seq(names(cleaned_baseline_plasma_data_adj)), log)

## create adjusted analysis matrices
ymat_adj <- t(cleaned_log_baseline_plasma_data_adj[, seq])

## fit adjusted model
fit_adj <- lmFit(ymat_adj, design_adj)
fit_adj <- eBayes(fit_adj)
results_ckm_adj <- topTable(fit_adj, coef = 2, number = nrow(ymat_adj))
results_ckm_adj <- results_ckm_adj[order(results_ckm_adj$P.Value), ]