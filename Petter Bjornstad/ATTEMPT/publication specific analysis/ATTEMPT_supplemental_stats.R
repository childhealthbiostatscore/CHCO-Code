# =============================================================================
# ATTEMPT MS - Supplemental statistics tables
# -----------------------------------------------------------------------------
# Builds a single xlsx file with one sheet per clinical/quantitative figure
# or table in "ATTEMPT MS figures and tables.qmd".
#
# Each sheet contains:
#   - Descriptive statistics (N, Mean, SD, Median, Q1, Q3) by treatment x visit
#     (or category) for the variable shown on the figure.
#   - Model-based estimates (paired t-tests, DiD, LMM emmeans pairwise contrasts)
#     that produced the stars on the figure, with 95% CIs and p-values.
#
# Output: <publication folder>/ATTEMPT_supplemental_stats.xlsx
#
# Run this script in the SAME R environment used for the qmd (i.e. with the
# Bjornstad keys.json available so s3readRDS() works). It does not modify the
# qmd or any source data.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
  library(aws.s3)
  library(openxlsx)
  library(arsenal)
  library(lmerTest)
  library(emmeans)
  library(broom)
  library(broom.mixed)
  library(growthcleanr)
})

# -----------------------------------------------------------------------------
# 0. Paths & S3 setup (mirrors the qmd)
# -----------------------------------------------------------------------------
keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
Sys.setenv(
  "AWS_ACCESS_KEY_ID"     = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION"    = "",
  "AWS_REGION"            = "",
  "AWS_S3_ENDPOINT"       = "s3.kopah.uw.edu"
)

user <- Sys.info()[["user"]]
if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path  <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path  <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root_path / git_path for this user.")
}

out_xlsx <- file.path(git_path, "ATTEMPT", "publication specific analysis",
                      "ATTEMPT_supplemental_stats.xlsx")

# -----------------------------------------------------------------------------
# 1. Load data (same as qmd) and replicate the data prep
# -----------------------------------------------------------------------------
attempt_dat <- s3readRDS(object = "cleaned_data/attempt_dat.rds",
                         bucket = "attempt", region = "")
delta_df    <- s3readRDS(object = "cleaned_data/attempt_delta_df.rds",
                         bucket = "attempt", region = "")
load(file = file.path(root_path, "/CROCODILE/Data_Cleaned/croc_data.RData"))
croc_dat <- dat

# Unit conversions (match qmd)
attempt_dat$creatinine_s_mgdl      <- attempt_dat$creatinine_s / 88.4
attempt_dat$emu_urine_acr_mean_mgg <- attempt_dat$emu_urine_acr_mean * 8.84

# Carry baseline (visit -4) values forward to visit 0 for MRI/mGFR vars (match qmd)
vars_to_copy <- c('mgfr_jodal_bsa', 'mgfr_jodal', 'mri_r2_cortex_l',
                  'mri_r2_cortex_r', 'avg_c_r2', 'mri_r2_kidney_l',
                  'mri_r2_kidney_r', 'avg_k_r2', 'avg_m_r2')
v4_vals <- attempt_dat[attempt_dat$visit == -4, c('subject_id', vars_to_copy)]
for (sid in unique(v4_vals$subject_id)) {
  patient_vals <- v4_vals[v4_vals$subject_id == sid, vars_to_copy]
  for (v in vars_to_copy) {
    attempt_dat[attempt_dat$subject_id == sid & attempt_dat$visit == 0, v] <-
      patient_vals[[v]]
  }
}

# BMI percentile (used in the supplemental cohort-compare table)
attempt_bmip <- tryCatch({
  subs <- attempt_dat %>%
    dplyr::mutate(sex  = case_when(sex == "Male" ~ 0, sex == "Female" ~ 1),
                  agem = age * 12) %>%
    dplyr::select(subject_id, visit, agem, wt = weight, ht = height, bmi, sex)
  growthcleanr::ext_bmiz(data = subs) %>%
    dplyr::select(subject_id, visit, bmip)
}, error = function(e) {
  message("growthcleanr::ext_bmiz failed: ", conditionMessage(e),
          " - bmip will be NA for TS2 cohort-compare table.")
  NULL
})
if (!is.null(attempt_bmip)) {
  attempt_dat <- attempt_dat %>%
    dplyr::select(-any_of("bmip")) %>%
    left_join(attempt_bmip, by = c("subject_id", "visit"))
} else {
  attempt_dat$bmip <- NA_real_
}

# Histology factor levels (match qmd)
attempt_dat$`Glomeruli sclerosed` <-
  factor(attempt_dat$`Glomeruli sclerosed`, levels = c("None", "Global", "Segmental"))
attempt_dat$`GBM thickening`      <-
  factor(attempt_dat$`GBM thickening`, levels = c("None", "Mild", "Moderate", "Severe"))
attempt_dat$`Vessel pathology`    <-
  factor(attempt_dat$`Vessel pathology`, levels = c("None", "Arteriolosclerosis", "Arteriolohyalinosis"))

# Cohort-compare grouping
attempt_dat$biopsy_cohort <- paste0(
  attempt_dat$treatment, ":",
  ifelse(attempt_dat$site %in% c("London", "Toronto"),
         "Non-biopsy cohort", "Biopsy-eligible cohort")
)

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------

# Descriptive stats for a continuous variable
desc_stats <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(data.frame(N = 0, Mean = NA, SD = NA,
                      Median = NA, Q1 = NA, Q3 = NA))
  }
  data.frame(
    N      = length(x),
    Mean   = mean(x),
    SD     = stats::sd(x),
    Median = stats::median(x),
    Q1     = unname(stats::quantile(x, 0.25)),
    Q3     = unname(stats::quantile(x, 0.75))
  )
}

# Descriptive stats by group(s)
desc_by <- function(df, value, ...) {
  by_vars <- rlang::ensyms(...)
  df %>%
    group_by(!!!by_vars) %>%
    summarise(
      N      = sum(!is.na(.data[[value]])),
      Mean   = mean(.data[[value]], na.rm = TRUE),
      SD     = stats::sd(.data[[value]], na.rm = TRUE),
      Median = stats::median(.data[[value]], na.rm = TRUE),
      Q1     = unname(stats::quantile(.data[[value]], 0.25, na.rm = TRUE)),
      Q3     = unname(stats::quantile(.data[[value]], 0.75, na.rm = TRUE)),
      .groups = "drop"
    )
}

# Categorical N(%) by group (preserves factor levels with 0 counts)
cat_by <- function(df, var, group) {
  df %>%
    filter(!is.na(.data[[var]])) %>%
    dplyr::count(!!sym(group), !!sym(var), name = "N", .drop = FALSE) %>%
    group_by(!!sym(group)) %>%
    mutate(Pct = round(100 * N / sum(N), 1)) %>%
    ungroup() %>%
    rename(Variable_level = !!sym(var)) %>%
    mutate(Variable = var, .before = 1)
}

# Run a treatment x visit LMM on the given outcome (typically a delta from
# baseline pre-computed in delta_df) and extract per-visit between-arm
# contrasts (Dapa - Placebo) via emmeans. Covariates match the qmd captions.
# If the outcome at the baseline visit is 0 for every subject (true deltas),
# the function still fits but the baseline contrast will be 0 by construction.
fit_lmm_treat_visit <- function(data, y, visits, baseline_visit,
                                include_site = FALSE,
                                bonferroni_adjust = NULL,
                                covariates = c("age", "sex",
                                               "diabetes_dx_duration", "bmi")) {
  d <- data %>%
    dplyr::filter(visit %in% visits, !is.na(.data[[y]])) %>%
    dplyr::mutate(
      treatment = factor(treatment, levels = c("Placebo", "Dapagliflozin 5mg")),
      visit_f   = factor(visit, levels = visits)
    )
  # Keep only covariates that exist in the data (delta_df may be slimmer)
  cov_use <- covariates[covariates %in% names(d)]
  if (include_site && "site" %in% names(d)) cov_use <- c(cov_use, "site")
  rhs_cov <- if (length(cov_use)) paste("+", paste(cov_use, collapse = " + ")) else ""
  base_fml <- paste0(y, " ~ treatment * visit_f ", rhs_cov,
                     " + (1 | subject_id)")
  fit <- tryCatch(
    lmerTest::lmer(as.formula(base_fml), data = d, REML = TRUE),
    error = function(e) {
      message("LMM failed for ", y, ": ", conditionMessage(e),
              " - retrying without covariates")
      tryCatch(
        lmerTest::lmer(as.formula(paste0(y, " ~ treatment * visit_f + (1 | subject_id)")),
                       data = d, REML = TRUE),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(fit)) return(NULL)

  if (is.null(bonferroni_adjust)) {
    adj <- if (length(visits) > 2) "bonferroni" else "none"
  } else adj <- bonferroni_adjust

  # Between-arm at each visit (Dapagliflozin - Placebo)
  em_between <- emmeans::emmeans(fit, ~ treatment | visit_f)
  con_between <- pairs(em_between, reverse = TRUE, adjust = adj) %>%
    summary(infer = TRUE) %>% as.data.frame() %>%
    transmute(Visit = visit_f, Contrast = "Dapagliflozin - Placebo",
              Estimate = estimate, SE, df, Lower_CI = lower.CL, Upper_CI = upper.CL,
              p_value = p.value)

  # Within-arm change (post - baseline)
  em_within <- emmeans::emmeans(fit, ~ visit_f | treatment)
  con_within <- pairs(em_within, reverse = TRUE, adjust = adj) %>%
    summary(infer = TRUE) %>% as.data.frame() %>%
    transmute(Treatment = treatment,
              Contrast = as.character(contrast),
              Estimate = estimate, SE, df,
              Lower_CI = lower.CL, Upper_CI = upper.CL,
              p_value = p.value)

  list(model_formula = base_fml,
       between_visits = con_between,
       within_treatment = con_within,
       n_used = nobs(fit),
       adjustment = adj)
}

# Ensure required covariates / metadata exist in delta_df by left-joining
# subject-level fields from attempt_dat (baseline visit 0).
ensure_delta_meta <- function(delta_df, attempt_dat,
                              vars = c("age", "sex", "diabetes_dx_duration",
                                       "bmi", "site")) {
  need <- vars[!vars %in% names(delta_df)]
  if (length(need) == 0) return(delta_df)
  base <- attempt_dat %>% dplyr::filter(visit == 0) %>%
    dplyr::select(any_of(c("subject_id", need))) %>%
    dplyr::distinct(subject_id, .keep_all = TRUE)
  delta_df %>% dplyr::left_join(base, by = "subject_id")
}

# Add a section row to a sheet-data list
push <- function(lst, title, df) {
  c(lst, list(list(title = title, df = df)))
}

# Write a list of (title, df) blocks to one sheet with section headers
write_sheet_sections <- function(wb, sheet, blocks) {
  if (!sheet %in% names(wb)) addWorksheet(wb, sheet)
  hdr_style <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                           border = "TopBottomLeftRight")
  sec_style <- createStyle(textDecoration = "bold", fontSize = 12,
                           fgFill = "#4472C4", fontColour = "#FFFFFF")
  row <- 1
  for (b in blocks) {
    writeData(wb, sheet, x = b$title, startRow = row, startCol = 1)
    addStyle(wb, sheet, sec_style, rows = row, cols = 1, gridExpand = TRUE)
    row <- row + 1
    writeData(wb, sheet, x = b$df, startRow = row, startCol = 1,
              headerStyle = hdr_style)
    row <- row + nrow(b$df) + 2
  }
  setColWidths(wb, sheet, cols = 1:20, widths = "auto")
}

# -----------------------------------------------------------------------------
# 3. Build workbook
# -----------------------------------------------------------------------------
wb <- createWorkbook()

# ---- README sheet ---------------------------------------------------------
addWorksheet(wb, "README")
readme_df <- data.frame(
  Sheet = c(
    "T1_Biopsy_baseline",
    "T2_Histology",
    "F1A_mGFR_by_baseline_cat",
    "F1B_Kidney_R2_full_cohort",
    "F1C_CROC_Medulla_R2",
    "F1D_mGFR_Denver",
    "F1E_HbA1c_Denver",
    "F1F_TIR_Denver",
    "F1G_Kidney_R2_Denver",
    "F1H_Medulla_R2_Denver",
    "F1I_Cortex_R2_Denver",
    "TS1_Full_baseline",
    "TS2_Cohort_compare"
  ),
  Figure = c(
    "Table 1 (tbl-co-baseline)",
    "Table 2 (tbl-histology)",
    "Fig 1A (mgfr_jodal_bsa_category_bars)",
    "Fig 1B (delta_avg_k_r2_overtime)",
    "Fig 1C (croc_r2_p / fig-croc-m-r2)",
    "Fig 1D (delta_mgfr_jodal_bsa_denver)",
    "Fig 1E (delta_hba1c_overtime_2_denver)",
    "Fig 1F (delta_tir_overtime_denver)",
    "Fig 1G (delta_avg_k_r2_overtime_denver)",
    "Fig 1H (delta_avg_m_r2_overtime_denver)",
    "Fig 1I (delta_avg_c_r2_overtime_denver)",
    "Table S1 (tbl-full-baseline)",
    "Table S2 (tbl-full-baseline-cohort-compare)"
  ),
  Description = c(
    "Baseline characteristics of biopsy-eligible cohort (Denver, visit 0)",
    "Qualitative histology of baseline biopsies (Denver, biopsy_yn = Yes)",
    "Mean delta mGFR (V-4 -> V16) by baseline mGFR category, full cohort",
    "Delta kidney R2* from baseline by treatment x visit, full cohort",
    "Medulla R2*: within-arm deltas (Placebo, Dapa), DiD, and T1D - HC baseline",
    "Delta mGFR from baseline by treatment x visit, biopsy-eligible cohort (Denver)",
    "Delta HbA1c from baseline by treatment x visit, biopsy-eligible cohort (Denver)",
    "Delta TIR from baseline by treatment x visit, biopsy-eligible cohort (Denver)",
    "Delta kidney R2* from baseline by treatment x visit, biopsy-eligible cohort (Denver)",
    "Delta medulla R2* from baseline by treatment x visit, biopsy-eligible cohort (Denver)",
    "Delta cortex R2* from baseline by treatment x visit, biopsy-eligible cohort (Denver)",
    "Baseline characteristics of full cohort (visit 0)",
    "Baseline characteristics: full vs biopsy-eligible cohort (visit 0)"
  )
)
writeData(wb, "README", readme_df, headerStyle = createStyle(textDecoration = "bold",
                                                             fgFill = "#4472C4",
                                                             fontColour = "#FFFFFF"))
setColWidths(wb, "README", cols = 1:3, widths = c(30, 40, 90))

# Continuous baseline variables used in Tables 1, S1, S2
baseline_vars <- c("age", "weight", "height", "bmi", "hba1c",
                   "emu_urine_acr_mean_mgg", "creatinine_s_mgdl",
                   "cystatin_c_s", "albumin_serum_gl", "cgm_tir",
                   "diabetes_dx_duration", "sbp", "dbp",
                   "cholesterol_serum_mmoll", "triglycerides_serum_mmoll",
                   "pulse", "ldl_serum_mmoll", "hdl_serum_mmoll",
                   "mgfr_jodal_bsa", "mgfr_jodal",
                   "avg_c_r2", "avg_k_r2", "avg_m_r2", "tkv")

# Build per-figure data prep frame for cohort filters
biopsy_baseline <- subset(attempt_dat, visit == 0 & site == "Denver")
full_baseline   <- subset(attempt_dat, visit == 0)

# ---- Table 1: Biopsy baseline -------------------------------------------
t1_cont_stats <- lapply(intersect(baseline_vars, names(biopsy_baseline)), function(v) {
  desc_by(biopsy_baseline, v, treatment) %>% mutate(Variable = v, .before = 1)
}) %>% bind_rows()
t1_sex <- cat_by(biopsy_baseline, "sex", "treatment")
write_sheet_sections(wb, "T1_Biopsy_baseline",
  list(list(title = "Continuous variables - by treatment (Denver, visit 0)",
            df = t1_cont_stats),
       list(title = "Categorical: sex - N (Pct) by treatment",
            df = t1_sex)))

# ---- Table 2: Histology -------------------------------------------------
histo_df <- subset(attempt_dat, visit == 0 & site == "Denver" & biopsy_yn == "Yes")
histo_vars <- c("Glomeruli number", "Glomeruli sclerosed", "GBM thickening",
                "Tubular atrophy", "Vessel pathology")
t2_blocks <- list()
for (v in histo_vars) {
  if (!v %in% names(histo_df)) next
  x <- histo_df[[v]]
  # Treat as continuous if numeric, categorical otherwise
  if (is.numeric(x)) {
    t2_blocks <- push(t2_blocks,
      paste0("Continuous: ", v, " - by treatment"),
      desc_by(histo_df, v, treatment) %>% mutate(Variable = v, .before = 1))
  } else {
    t2_blocks <- push(t2_blocks,
      paste0("Categorical: ", v, " - N (Pct) by treatment"),
      cat_by(histo_df, v, "treatment"))
  }
}
write_sheet_sections(wb, "T2_Histology", t2_blocks)

# ---- Figure 1A: mGFR by baseline mGFR category --------------------------
# Categorize baseline mGFR (visit -4) into <90, 90-110, >=110
f1a_dat <- attempt_dat %>%
  dplyr::mutate(treatment = case_when(treatment == "Dapagliflozin 5mg" ~ "Dapagliflozin",
                                      TRUE ~ treatment)) %>%
  dplyr::filter(visit %in% c(-4, 16), !is.na(mgfr_jodal_bsa))
# Baseline category from visit -4
f1a_base <- f1a_dat %>% dplyr::filter(visit == -4) %>%
  transmute(subject_id, baseline_mgfr = mgfr_jodal_bsa,
            mgfr_cat = cut(mgfr_jodal_bsa, breaks = c(-Inf, 90, 110, Inf),
                           labels = c("<90", "90-110", ">=110"), right = FALSE))
f1a_long <- f1a_dat %>%
  inner_join(f1a_base %>% select(subject_id, mgfr_cat), by = "subject_id") %>%
  arrange(subject_id, visit) %>%
  group_by(subject_id, treatment, mgfr_cat) %>%
  dplyr::filter(all(c(-4, 16) %in% visit)) %>%
  summarise(delta = mgfr_jodal_bsa[visit == 16] - mgfr_jodal_bsa[visit == -4],
            .groups = "drop") %>%
  dplyr::filter(!is.na(delta))

f1a_desc <- desc_by(f1a_long, "delta", treatment, mgfr_cat) %>%
  rename(Treatment = treatment, mGFR_category = mgfr_cat) %>%
  arrange(mGFR_category, Treatment)

# Within-arm paired t-tests vs 0
f1a_within <- f1a_long %>%
  group_by(treatment, mgfr_cat) %>%
  summarise(
    N        = n(),
    Mean     = mean(delta, na.rm = TRUE),
    SE       = stats::sd(delta, na.rm = TRUE) / sqrt(N),
    Lower_CI = Mean + qt(0.025, df = N - 1) * SE,
    Upper_CI = Mean + qt(0.975, df = N - 1) * SE,
    p_value  = if (N > 1) t.test(delta, mu = 0)$p.value else NA_real_,
    .groups  = "drop"
  ) %>%
  rename(Treatment = treatment, mGFR_category = mgfr_cat)

# DiD via lm with treatment x category interaction (explicit factor levels)
f1a_did <- tryCatch({
  d <- f1a_long %>%
    dplyr::mutate(treatment = factor(treatment, levels = c("Placebo", "Dapagliflozin")))
  fit <- lm(delta ~ treatment * mgfr_cat, data = d)
  em  <- emmeans::emmeans(fit, ~ treatment | mgfr_cat)
  pairs(em, reverse = TRUE) %>% summary(infer = TRUE) %>% as.data.frame() %>%
    transmute(mGFR_category = mgfr_cat,
              Contrast = "Dapagliflozin - Placebo",
              Estimate = estimate, SE, df,
              Lower_CI = lower.CL, Upper_CI = upper.CL,
              p_value = p.value)
}, error = function(e) data.frame())

write_sheet_sections(wb, "F1A_mGFR_by_baseline_cat", list(
  list(title = "Descriptives: Delta mGFR (V16 - V-4) by treatment x baseline mGFR category",
       df = f1a_desc),
  list(title = "Within-arm paired t-tests vs 0 (delta mGFR)",
       df = f1a_within),
  list(title = "Difference-in-differences (lm: delta ~ treatment * category; emmeans Dapagliflozin - Placebo)",
       df = f1a_did)
))

# ---- Figure 1B: Kidney R2* over time, full cohort (DELTAS from delta_df)
# delta_df stores within-subject changes from baseline under the original
# variable name; baseline visit values are 0 for every subject.
delta_df_meta <- ensure_delta_meta(delta_df, attempt_dat)

f1b_dat  <- delta_df_meta %>% dplyr::filter(visit %in% c(-4, 16))
f1b_desc <- desc_by(f1b_dat, "avg_k_r2", treatment, visit) %>%
  arrange(visit, treatment)
f1b_mod  <- fit_lmm_treat_visit(delta_df_meta, "avg_k_r2",
                                visits = c(-4, 16), baseline_visit = -4,
                                include_site = TRUE)
f1b_blocks <- list(list(title = "Descriptives: delta Kidney R2* (from baseline) by treatment x visit",
                        df = f1b_desc))
if (!is.null(f1b_mod)) {
  f1b_blocks <- push(f1b_blocks,
    paste0("LMM on deltas (", f1b_mod$model_formula, "); n=", f1b_mod$n_used,
           "; adjustment=", f1b_mod$adjustment),
    data.frame(Info = "See contrasts below. Baseline contrast is 0 by construction."))
  f1b_blocks <- push(f1b_blocks, "Between-arm at each visit (Dapagliflozin - Placebo)",
                     f1b_mod$between_visits)
  f1b_blocks <- push(f1b_blocks, "Within-arm change (post - baseline)",
                     f1b_mod$within_treatment)
}
write_sheet_sections(wb, "F1B_Kidney_R2_full_cohort", f1b_blocks)

# ---- Figure 1C: CROCODILE reversal medulla R2* -------------------------
# Coerce visit to character so the two cohorts can be combined cleanly
croc_part <- croc_dat %>%
  dplyr::select(subject_id = record_id, visit, group, avg_m_r2) %>%
  dplyr::mutate(visit = as.character(visit))
attempt_part <- attempt_dat %>% dplyr::filter(site == "Denver") %>%
  dplyr::select(subject_id, visit, group = treatment, avg_m_r2) %>%
  dplyr::mutate(visit = as.character(visit),
                subject_id = as.character(subject_id))
croc_comp <- bind_rows(croc_part, attempt_part) %>%
  dplyr::filter(!is.na(avg_m_r2))

# Descriptives per group/visit
f1c_desc <- desc_by(croc_comp, "avg_m_r2", group, visit) %>% arrange(group, visit)

# Within-arm (ATTEMPT) paired changes - require both visits
paired_changes <- croc_comp %>%
  dplyr::filter(group %in% c("Placebo", "Dapagliflozin 5mg"),
                visit %in% c("-4", "16")) %>%
  dplyr::mutate(visit_num = as.numeric(visit)) %>%
  group_by(subject_id, group) %>%
  dplyr::filter(n_distinct(visit_num) == 2) %>%
  summarise(delta = avg_m_r2[visit_num == 16] - avg_m_r2[visit_num == -4],
            .groups = "drop")
f1c_within <- paired_changes %>% group_by(group) %>%
  summarise(N = n(),
            Mean_delta = mean(delta, na.rm = TRUE),
            SD_delta   = stats::sd(delta, na.rm = TRUE),
            SE         = SD_delta / sqrt(N),
            Lower_CI   = Mean_delta + qt(0.025, df = N - 1) * SE,
            Upper_CI   = Mean_delta + qt(0.975, df = N - 1) * SE,
            p_value    = if (N > 1) t.test(delta, mu = 0)$p.value else NA_real_,
            .groups    = "drop") %>%
  rename(Group = group)

# DiD (Dapa - Placebo) of within-arm deltas
paired_changes$group <- factor(paired_changes$group,
                               levels = c("Placebo", "Dapagliflozin 5mg"))
tt_did <- t.test(delta ~ group, data = paired_changes)
did_mean <- unname(diff(tt_did$estimate))         # Dapa - Placebo
did_ci   <- -rev(tt_did$conf.int)
f1c_did <- data.frame(
  Contrast = "Dapagliflozin - Placebo (DiD of delta medulla R2*)",
  N        = nrow(paired_changes),
  Estimate = did_mean,
  Lower_CI = did_ci[1], Upper_CI = did_ci[2],
  p_value  = tt_did$p.value
)

# T1D - HC baseline (CROCODILE)
base_cmp <- croc_comp %>%
  dplyr::filter(group %in% c("Lean Control", "Type 1 Diabetes"),
                visit == "baseline") %>%
  dplyr::mutate(group = factor(group, levels = c("Lean Control", "Type 1 Diabetes")))
tt_t1d <- t.test(avg_m_r2 ~ group, data = base_cmp)
t1d_minus_lean_mean <- unname(diff(tt_t1d$estimate))
t1d_minus_lean_ci   <- -rev(tt_t1d$conf.int)
f1c_t1d <- data.frame(
  Contrast = "Type 1 Diabetes - Lean Control (baseline medulla R2*)",
  N_T1D    = sum(base_cmp$group == "Type 1 Diabetes"),
  N_HC     = sum(base_cmp$group == "Lean Control"),
  Estimate = t1d_minus_lean_mean,
  Lower_CI = t1d_minus_lean_ci[1], Upper_CI = t1d_minus_lean_ci[2],
  p_value  = tt_t1d$p.value
)

write_sheet_sections(wb, "F1C_CROC_Medulla_R2", list(
  list(title = "Descriptives: Medulla R2* by group x visit (ATTEMPT Denver + CROCODILE)",
       df = f1c_desc),
  list(title = "Within-arm paired t-tests (delta = V16 - V-4) vs 0 - ATTEMPT Denver",
       df = f1c_within),
  list(title = "DiD on within-arm deltas (Dapagliflozin - Placebo)", df = f1c_did),
  list(title = "Baseline T1D - HC contrast (Welch's t-test)",        df = f1c_t1d)
))

# ---- Figures 1D-I: Denver timelines (DELTAS from delta_df) -------------
# Stats here match what the figures plot: within-subject deltas from baseline.
# Baseline visit deltas are 0 for every subject in both arms; the between-arm
# contrast is only meaningful at the post-baseline visits.
denver_delta <- subset(delta_df_meta, site == "Denver")

# Resolve a column name in `denver_delta` from a list of candidates
resolve_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NULL else hit[1]
}

timeline_specs <- list(
  list(sheet = "F1D_mGFR_Denver",       y_cands = c("mgfr_jodal_bsa"),  visits = c(-4, 16),   base = -4),
  list(sheet = "F1E_HbA1c_Denver",      y_cands = c("hba1c"),           visits = c(0, 4, 16), base = 0),
  list(sheet = "F1F_TIR_Denver",        y_cands = c("tir", "cgm_tir"),  visits = c(0, 16),    base = 0),
  list(sheet = "F1G_Kidney_R2_Denver",  y_cands = c("avg_k_r2"),        visits = c(-4, 16),   base = -4),
  list(sheet = "F1H_Medulla_R2_Denver", y_cands = c("avg_m_r2"),        visits = c(-4, 16),   base = -4),
  list(sheet = "F1I_Cortex_R2_Denver",  y_cands = c("avg_c_r2"),        visits = c(-4, 16),   base = -4)
)

for (spec in timeline_specs) {
  y_var <- resolve_col(denver_delta, spec$y_cands)
  if (is.null(y_var)) {
    message("Skipping ", spec$sheet, ": none of ",
            paste(spec$y_cands, collapse = "/"), " in delta_df (Denver).")
    next
  }
  d_sub <- denver_delta %>% dplyr::filter(visit %in% spec$visits)
  desc  <- desc_by(d_sub, y_var, treatment, visit) %>% arrange(visit, treatment)
  mod   <- fit_lmm_treat_visit(denver_delta, y_var,
                               visits = spec$visits, baseline_visit = spec$base,
                               include_site = FALSE)
  blocks <- list(list(title = paste0("Descriptives: delta ", y_var,
                                     " (from baseline) by treatment x visit"),
                      df = desc))
  if (!is.null(mod)) {
    blocks <- push(blocks,
      paste0("LMM on deltas (", mod$model_formula, "); n=", mod$n_used,
             "; adjustment=", mod$adjustment),
      data.frame(Info = "See contrasts below. Baseline contrast is 0 by construction."))
    blocks <- push(blocks, "Between-arm at each visit (Dapagliflozin - Placebo)",
                   mod$between_visits)
    blocks <- push(blocks, "Within-arm change (post - baseline)",
                   mod$within_treatment)
  }
  write_sheet_sections(wb, spec$sheet, blocks)
}

# ---- Table S1: Full cohort baseline ------------------------------------
ts1_cont <- lapply(intersect(baseline_vars, names(full_baseline)), function(v) {
  desc_by(full_baseline, v, treatment) %>% mutate(Variable = v, .before = 1)
}) %>% bind_rows()
ts1_sex <- cat_by(full_baseline, "sex", "treatment")
write_sheet_sections(wb, "TS1_Full_baseline", list(
  list(title = "Continuous variables - by treatment (visit 0, full cohort)",
       df = ts1_cont),
  list(title = "Categorical: sex - N (Pct) by treatment",
       df = ts1_sex)
))

# ---- Table S2: Full vs Biopsy cohort comparison ------------------------
ts2_vars <- c(baseline_vars, "bmip")
ts2_cont <- lapply(intersect(ts2_vars, names(full_baseline)), function(v) {
  desc_by(full_baseline, v, biopsy_cohort) %>% mutate(Variable = v, .before = 1)
}) %>% bind_rows()
ts2_sex <- cat_by(full_baseline, "sex", "biopsy_cohort")
write_sheet_sections(wb, "TS2_Cohort_compare", list(
  list(title = "Continuous variables - by treatment x cohort (visit 0)",
       df = ts2_cont),
  list(title = "Categorical: sex - N (Pct) by treatment x cohort",
       df = ts2_sex)
))

# -----------------------------------------------------------------------------
# 4. Save
# -----------------------------------------------------------------------------
saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
message("Wrote: ", out_xlsx)
message("Sheets: ", paste(names(wb), collapse = ", "))



