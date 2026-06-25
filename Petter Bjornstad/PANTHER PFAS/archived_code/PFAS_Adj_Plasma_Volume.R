## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
#load libraries
library(ggplot2)
library(tidyverse)
library(janitor)
library(lme4)
library(lmerTest)
library(broom)
library(broom.mixed)
library(emmeans)
library(ggeffects)
library(patchwork)
library(gWQS)
library(corrplot)
library(gt)
library(table1)
library(lme4)
library(ggpubr)
library(pheatmap)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
#Directories
#computer <- "mac studio"
computer <- "mac laptop"
if (computer == "mac studio") {
  user <- Sys.info()[["user"]]
  if (user == "sleidholt") {
    dir.dat <- c("/Users/Shared/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
    dir.results <- c("/Users/Shared/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Savanah Leidholt/PFAS/Results/total_plasma_volume_adjusted_results")
    git_path <- "/Users/sleidholt/Documents/GitHub/CHCO-Code/Petter Bjornstad/PANTHER PFAS/"
  }
} else {
  user <- Sys.info()[["user"]]
  if (user == "savanahleidholt") {
    dir.dat <- c("/Users/savanahleidholt/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
    dir.results <- c("/Users/savanahleidholt/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Savanah Leidholt/PFAS/Results/total_plasma_volume_adjusted_results")
    git_path <- "/Users/savanahleidholt/Desktop/CHCO-Code/Petter Bjornstad/PANTHER PFAS"
  } 
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
source(file.path(git_path, "Functions.R"))

pfas_imputed <- read_rds(file = file.path(dir.dat, "Savanah Leidholt", "PFAS", "Data", "pfas_imputed.rds"))
baseline_acru <- read_rds(file = file.path(dir.dat, "Savanah Leidholt", "PFAS", "Data", "baseline_acru.rds"))
meta_nopfas <- read_rds(file = file.path(dir.dat, "Savanah Leidholt", "PFAS", "Data", "meta_nopfas.rds"))
dat_baseline <- read_rds(file = file.path(dir.dat, "Savanah Leidholt", "PFAS", "Data", "dat_baseline.rds"))
dat_long <- read_rds(file = file.path(dir.dat, "Savanah Leidholt", "PFAS", "Data", "dat_long.rds"))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------

pfas_all <- c(
  "N.EtFOSAA","N.MeFOSAA","PFBA","PFDA","PFDoA","PFHpA","PFHxA","PFNA","PFOA",
  "PFBS","PFPeA","PFTeDA","PFTrDA","PFUnA","PFDoS","PFHps","PFHxS","PFNS",
  "PFOS","PFOSA","PFPeAS"
)

pfas_vars <- c(
  "N.MeFOSAA","PFDA","PFHpA","PFNA","PFBS",
  "PFHps","PFHxS","PFOA","PFOS","PFPeAS"
)

baseline_model_outcomes <- c(
  "log_acr_u",
 # "log_microalbumin_u",
  "eGFR_CKiD_U25_avg",
  "erpf_bsa_plasma",
  "avg_c_r2",
  "avg_k_r2",
  "bmi",
  "pah_clear_bsa",
  "dexa_lean_kg",
  "dexa_fat_kg",
  "igf_1",
  "bl_dheas",
  "mm_si",
  "mm_ir",
  "mm_di"
)

group_time_outcomes <- c(
  "log_acr_u",
  #"log_microalbumin_u",
  "eGFR_CKiD_U25_avg",
  "erpf_bsa_plasma",
  "pah_clear_bsa",
  "avg_c_r2",
  "avg_k_r2",
  "bmi",
  "dexa_lean_kg",
  "dexa_fat_kg",
  "igf_1",
  "bl_dheas"
)


outcome_label_map <- c(
  #kidney function
  acr_u = "uACR",
  log_acr_u = "log1p(uACR)",
 # microalbumin_u = "Microalbuminuria",
 # log_microalbumin_u = "log1p(Microalbuminuria)",
  eGFR_CKiD_U25_avg = "eGFR",
  erpf_bsa_plasma = "ERPF",
  pah_clear_bsa = "PAH clearance",
  #kidney oxygenation
  avg_c_r2 = "Cortical R2",
  avg_k_r2= "Kidney R2",
  #metabolism
  bmi = "Body Mass Index",
  dexa_lean_kg = "DEXA Lean Mass (kg)",
  dexa_fat_kg = "DEXA Fat Mass (kg)",
  mm_si = "Insulin Sensitivity",
  mm_ir = "Insulin Resistance",
  mm_di = "Disposition Index",
  #hormones
  igf_1 = "IGF-1",
  bl_dheas = "DHEA-S",
  tan_mgd = "Tanner Stage",
  tan_fgd = "Tanner Stage"
)

pfas_label_map <- c(
  log2_N.MeFOSAA_bl = "log2(N-MeFOSAA)",
  log2_PFDA_bl = "log2(PFDA)",
  log2_PFHpA_bl = "log2(PFHpA)",
  log2_PFNA_bl = "log2(PFNA)",
  log2_PFBS_bl = "log2(PFBS)",
  log2_PFHps_bl = "log2(PFHpS)",
  log2_PFHxS_bl = "log2(PFHxS)",
  log2_PFOA_bl = "log2(PFOA)",
  log2_PFOS_bl = "log2(PFOS)",
  log2_PFPeAS_bl = "log2(PFPeAS)"
)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------
fmt_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

get_pfas_label <- function(x) {
  dplyr::recode(x, !!!pfas_label_map, .default = x)
}

get_outcome_label <- function(outcome) {
  dplyr::recode(outcome, !!!outcome_label_map, .default = outcome)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# Outcome cleaning / log transforms

dat_baseline <- dat_baseline %>%
  mutate(
    microalbumin_u = ifelse(microalbumin_u > 3000, NA, microalbumin_u),
    acr_u = ifelse(acr_u > 2367, NA, acr_u),
    log_acr_u = log1p(acr_u),
    log_microalbumin_u = log1p(microalbumin_u),
    log_acr_u_bl = log1p(acr_u_bl)
  )

dat_long <- dat_long %>%
  mutate(
    microalbumin_u = ifelse(microalbumin_u > 3000, NA, microalbumin_u),
    acr_u = ifelse(acr_u > 2367, NA, acr_u),
    log_acr_u = log1p(acr_u),
    log_microalbumin_u = log1p(microalbumin_u),
    log_acr_u_bl = log1p(acr_u_bl)
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define baseline PFAS variable names
pfas_bl_vars <- paste0(pfas_vars, "_bl")

# names for plasma-total PFAS columns
pfas_total_ng_bl_vars <- paste0(pfas_bl_vars, "_total_ng")


# Add blood/plasma volume and plasma-total PFAS

dat_baseline <- dat_baseline %>%
  mutate(
    blood_vol_ml = case_when(
      sex == "Male" ~ (weight * (70/sqrt(bmi/22))),
      sex == "Female" ~ (weight * (60/sqrt(bmi/22)))),
    hct_prop = hct / 100,
    plasma_vol_ml = blood_vol_ml * (1 - hct_prop)
  ) %>%
  mutate(
    across(
      all_of(pfas_bl_vars),
      ~ .x * plasma_vol_ml,
      .names = "{.col}_total_ng"
    )
  )

dat_long <- dat_long %>%
  mutate(
    blood_vol_ml = case_when(
      sex == "Male" ~ (weight * (70/sqrt(bmi/22))),
      sex == "Female" ~ (weight * (60/sqrt(bmi/22)))),
    hct_prop = hct / 100,
    plasma_vol_ml = blood_vol_ml * (1 - hct_prop)
  ) %>%
  mutate(
    across(
      all_of(pfas_bl_vars),
      ~ .x * plasma_vol_ml,
      .names = "{.col}_total_ng"
    )
  )



## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create summed plasma-total PFAS burden
# and ln-transform summed burden

dat_baseline <- dat_baseline %>%
  mutate(
    n_pfas_total_nonmiss = rowSums(!is.na(across(all_of(pfas_total_ng_bl_vars)))),
    total_plasma_pfas_bl = ifelse(
      n_pfas_total_nonmiss == 0,
      NA_real_,
      rowSums(across(all_of(pfas_total_ng_bl_vars)), na.rm = TRUE)
    ),
    ln_total_plasma_pfas_bl = ifelse(
      total_plasma_pfas_bl > 0,
      log(total_plasma_pfas_bl),
      NA_real_
    )
  ) %>%
  select(-n_pfas_total_nonmiss)

dat_long <- dat_long %>%
  mutate(
    n_pfas_total_nonmiss = rowSums(!is.na(across(all_of(pfas_total_ng_bl_vars)))),
    total_plasma_pfas_bl = ifelse(
      n_pfas_total_nonmiss == 0,
      NA_real_,
      rowSums(across(all_of(pfas_total_ng_bl_vars)), na.rm = TRUE)
    ),
    ln_total_plasma_pfas_bl = ifelse(
      total_plasma_pfas_bl > 0,
      log(total_plasma_pfas_bl),
      NA_real_
    )
  ) %>%
  select(-n_pfas_total_nonmiss)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
dat_baseline <- dat_baseline %>%
  mutate(
    across(
      all_of(pfas_total_ng_bl_vars),
      ~ ifelse(.x > 0, log(.x), NA_real_),
      .names = "ln_{.col}"
    )
  )

dat_long <- dat_long %>%
  mutate(
    across(
      all_of(pfas_total_ng_bl_vars),
      ~ ifelse(.x > 0, log(.x), NA_real_),
      .names = "ln_{.col}"
    )
  )

# vector of log2 plasma-total PFAS variables for modeling
ln_pfas_total_bl_vars <- paste0("ln_", pfas_total_ng_bl_vars)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
dat_baseline2 <- dat_long %>% 
  select(record_id, visit, weight, height, hct, sex, group, bmi, age, sex) %>%
  filter(group == "Lean Control") %>%
  filter(visit == "baseline_screening") %>%
  mutate(
    height = height/100
  ) 

dat_baseline2 <- dat_baseline2 %>%
  mutate(
    # Calculate Nadler blood volume (L) by sex
    nadler_bv_ml = case_when(
      sex == "Male"   ~ ((height^3 * 0.3669) + (0.03219 * weight) + 0.6041) * 1000,
      sex == "Female" ~ ((height^3 * 0.3561) + (0.03208 * weight) + 0.1833) * 1000),
    LBB_bv_ml = case_when(
      sex == "Male" ~ (weight * (70/sqrt(bmi/22))),
      sex == "Female" ~ (weight * (60/sqrt(bmi/22)))
    )
    )

dat_baseline2 <- dat_baseline2 %>%
  mutate(
    hct_prop = hct / 100,
    LBB_plasma_ml  = LBB_bv_ml * (1 - hct_prop),
    nadler_plasma_ml  = nadler_bv_ml * (1 - hct_prop),
    pv_diff = LBB_plasma_ml - nadler_plasma_ml
  )





## --------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf(file.path(dir.results, "Nadler_v_LBB.pdf"), width = 10, height = 8)
# Keep only the two plasma volume estimates and pivot longer
plot_dat <- dat_baseline2 %>%
  pivot_longer(
    cols = c(LBB_plasma_ml, nadler_plasma_ml),
    names_to = "equation_type",
    values_to = "plasma_volume"
  ) %>%
  mutate(
    equation_type = recode(
      equation_type,
      LBB_plasma_ml = "Modified LBB",
      nadler_plasma_ml = "Nadler")
  )

# ANCOVA: test if mean difference ≠ 0, adjusting for age and sex
model <- lm(pv_diff ~ age + sex, data = dat_baseline2)
summary(model)

# Or if you want a formal intercept test (is the mean diff ≠ 0?):
emmeans(model, ~ 1)  # intercept = adjusted mean difference

#boxplot
 ggplot(plot_dat, aes(x = equation_type, y = plasma_volume, fill = equation_type)) +
  scale_fill_viridis_d(option = "C") +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(shape = sex), width = 0.1, alpha = 0.7) +
  stat_compare_means(paired = F, method = "t.test") +
  theme_bw() +
  labs(
    title = "Plasma volume estimates in Lean Controls by equation type",
    x = "Equation type",
    y = "Estimated plasma volume",
    size = "Sex"
  )
 
 dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# data for correlation matrix
pfas_corr_df <- dat_baseline %>%
  dplyr::select(all_of(ln_pfas_total_bl_vars)) %>%
  dplyr::rename_with(~ sapply(.x, get_pfas_label))

# correlation matrix using pairwise complete observations
pfas_corr_mat <- cor(
  pfas_corr_df,
  use = "pairwise.complete.obs",
  method = "spearman"   # use "pearson" if you prefer linear correlation
)

pfas_corr_mat


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
pfas_cor_tests <- expand.grid(
  var1 = ln_pfas_total_bl_vars,
  var2 = ln_pfas_total_bl_vars,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(var1 < var2) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    test = list(
      cor.test(
        dat_baseline[[var1]],
        dat_baseline[[var2]],
        method = "spearman",
        use = "pairwise.complete.obs",
        exact = FALSE
      )
    ),
    rho = test$estimate,
    p_value = test$p.value
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    var1_label = sapply(var1, get_pfas_label),
    var2_label = sapply(var2, get_pfas_label),
    p_fdr = p.adjust(p_value, method = "fdr")
  ) %>%
  dplyr::select(var1_label, var2_label, rho, p_value, p_fdr)

pfas_cor_tests

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
cor_mtest <- function(mat, method = "spearman") {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  colnames(p.mat) <- colnames(mat)
  rownames(p.mat) <- colnames(mat)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- suppressWarnings(
        cor.test(mat[, i], mat[, j], method = method, exact = FALSE)
      )
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  p.mat
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
pfas_corr_df <- dat_baseline %>%
  dplyr::select(all_of(ln_pfas_total_bl_vars)) %>%
  dplyr::rename_with(~ sapply(.x, get_pfas_label))

pfas_corr_mat <- cor(
  pfas_corr_df,
  use = "pairwise.complete.obs",
  method = "spearman"
)

pfas_p_mat <- cor_mtest(pfas_corr_df, method = "spearman")


## ----fig.width=16, fig.height=12-------------------------------------------------------------------------------------------------------------------------------
pdf(file = file.path(dir.results, "pfas_correlation_plot.pdf"))

corrplot::corrplot(
  pfas_corr_mat,
  method = "color",
  type = "upper",
  order = "hclust",
  tl.col = "black",
  tl.srt = 45,
  col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
  p.mat = pfas_p_mat,
  sig.level = c(0.001, 0.01, 0.05),
  insig = "label_sig",
  pch.col = "black",
  pch.cex = 1.2,
  mar = c(0, 0, 1, 0),
  title = "Spearman correlation matrix of baseline log2-transformed PFAS"
)

dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
run_baseline_group_model <- function(data, outcome, covars = c("age", "sex")) {
  
  rhs <- c("group", covars)
  m0<- as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + ")))
  
  model_data <- data %>%
    dplyr::select(all.vars(m0)) %>%
    tidyr::drop_na()
  
  m1 <- lm(m0, data = model_data)
  
  # omnibus test for group
  group_rows <- broom::tidy(stats::anova(m1)) %>%
    dplyr::filter(term == "group") %>%
    dplyr::mutate(
      outcome = outcome,
      n = nrow(model_data)
    )
  
  # pairwise Tukey-adjusted comparisons
  emm <- emmeans::emmeans(m1, ~ group)
  pairwise_tbl <- summary(
    pairs(emm, adjust = "tukey"),
    infer = TRUE
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      outcome = outcome,
      n = nrow(model_data)
    )
  
  list(
    model = m1,
    omnibus = group_rows,
    emmeans = as.data.frame(summary(emm)),
    pairwise = pairwise_tbl
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
pfas_group_vars <- c(ln_pfas_total_bl_vars, "ln_total_plasma_pfas_bl")

baseline_group_pfas_results <- lapply(
  pfas_group_vars,
  function(outcome) {
    run_baseline_group_model(
      data = dat_baseline,
      outcome = outcome,
      covars = c("age", "sex")
    )
  }
)
names(baseline_group_pfas_results) <- pfas_group_vars

baseline_group_pfas_omnibus <- bind_rows(
  lapply(baseline_group_pfas_results, `[[`, "omnibus")
) %>%
  dplyr::mutate(
    outcome_label = sapply(outcome, get_pfas_label),
    p_fdr = p.adjust(p.value, method = "fdr"),
    p_fmt = fmt_p(p.value),
    p_fdr_fmt = fmt_p(p_fdr)
  )

baseline_group_pfas_pairwise <- bind_rows(
  lapply(baseline_group_pfas_results, `[[`, "pairwise")
) %>%
  dplyr::mutate(
    outcome_label = sapply(outcome, get_pfas_label),
    p_fmt = fmt_p(p.value)
  )

baseline_group_pfas_omnibus
baseline_group_pfas_pairwise

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_baseline_pfas_group_boxplot <- function(data,
                                             exposure,
                                             results_list) {
  
  plot_data <- data %>%
    dplyr::select(group, sex, age, all_of(exposure)) %>%
    tidyr::drop_na()
  
  pair_tbl <- results_list[[exposure]]$pairwise %>%
    dplyr::mutate(
      group1 = stringr::str_trim(stringr::word(contrast, 1, sep = " - ")),
      group2 = stringr::str_trim(stringr::word(contrast, 2, sep = " - "))
    )
  
  y_max <- max(plot_data[[exposure]], na.rm = TRUE)
  y_min <- min(plot_data[[exposure]], na.rm = TRUE)
  y_step <- 0.08 * (y_max - y_min)
  
  pair_tbl <- pair_tbl %>%
    dplyr::mutate(
      y.position = y_max + seq_len(n()) * y_step,
      p_label = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  omnibus_p <- results_list[[exposure]]$omnibus$p.value[1]
  
  ggplot(plot_data, aes(x = group, y = .data[[exposure]], fill = group)) +
    scale_fill_viridis_d(option = "C") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 2) +
    ggpubr::stat_pvalue_manual(
      pair_tbl,
      label = "p_label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    ) +
    theme_bw() +
    labs(
      title = paste("Baseline", get_pfas_label(exposure), "by group"),
      subtitle = paste0("Omnibus group p = ", fmt_p(omnibus_p), "; pairwise p adjusted by Tukey"),
      x = "Group",
      y = get_pfas_label(exposure)
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------

plots_by_group <- lapply(
  pfas_group_vars,
  function(pfas) {
    plot_baseline_pfas_group_boxplot(
      data = dat_baseline,
      exposure = pfas,
      results_list = baseline_group_pfas_results
    )
  }
)

plots_by_group

combined_plot2 <- wrap_plots(plots_by_group, ncol = 3)
combined_plot2

pdf(file.path(dir.results, "pfas_by_group.pdf"), width = 18, height = 16)
print(combined_plot2)
dev.off()



## --------------------------------------------------------------------------------------------------------------------------------------------------------------
run_baseline_pfas_lm <- function(data, outcome, exposure, adjust_acr = TRUE) {
  
  rhs <- c(exposure, "age", "sex")
  
  m0 <- as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + ")))
  
  model_data <- data %>%
    dplyr::select(all.vars(m0)) %>%
    tidyr::drop_na()
  
  m1 <- lm(m0, data = model_data)
  
  list(
    model = m1,
    tidy = broom::tidy(m1) %>%
      dplyr::mutate(
        analysis = "baseline_cross_sectional",
        outcome = outcome,
        exposure = exposure,
        model_formula = paste(deparse(m0), collapse = " ")
      )
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
baseline_results_list <- lapply(baseline_model_outcomes, function(outcome) {
  lapply(pfas_group_vars, function(exposure) {
    run_baseline_pfas_lm(dat_baseline, outcome, exposure)
  })
})

baseline_results <- dplyr::bind_rows(
  lapply(baseline_results_list, function(x) dplyr::bind_rows(lapply(x, `[[`, "tidy")))
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
baseline_pfas_only <- baseline_results %>%
  dplyr::mutate(
    exposure_clean = stringr::str_replace_all(exposure, "`", "")
  ) %>%
  dplyr::filter(term %in% c(exposure, exposure_clean)) %>%
  dplyr::group_by(analysis, outcome) %>%
  dplyr::mutate(
    p_fdr = p.adjust(p.value, method = "fdr"),
    sig_fdr = p_fdr < 0.05,
    p_fmt = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    p_fdr_fmt = ifelse(p_fdr < 0.001, "<0.001", sprintf("%.3f", p_fdr))
  ) %>%
  dplyr::ungroup()

baseline_results_final <- baseline_results %>%
  dplyr::left_join(
    baseline_pfas_only %>%
      dplyr::select(analysis, outcome, exposure, term, p_fdr, sig_fdr, p_fdr_fmt),
    by = c("analysis", "outcome", "exposure", "term")
  )

baseline_table <- baseline_results_final %>%
  dplyr::filter(term %in% c(exposure, stringr::str_replace_all(exposure, "`", ""))) %>%
  dplyr::mutate(
    estimate_ci = paste0(
      round(estimate, 3), " (",
      round(estimate - 1.96 * std.error, 3), ", ",
      round(estimate + 1.96 * std.error, 3), ")"
    ),
    p_fmt = fmt_p(p.value)
  ) %>%
  dplyr::select(
    analysis, outcome, exposure, estimate, std.error, statistic, p.value,
    p_fdr, sig_fdr, estimate_ci, p_fmt, p_fdr_fmt
  )

baseline_table


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
forest_df <- baseline_table %>%
  dplyr::mutate(
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_list <- list()


for (outcome in unique(forest_df$outcome)) {
  df <- forest_df %>%
    dplyr::filter(outcome == !!outcome)
  p <- ggplot(df, aes(x = estimate, y = exposure)) +
  geom_point(aes(colour = sig_fdr), size = 2.5) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "plum")
) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = paste("PFAS associations with", get_outcome_label(outcome)),
    x = "Effect size (β)",
    y = "PFAS",
    shape = "FDR <0.05"
  )
plot_list[[outcome]] <- p
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_list

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
run_baseline_pfas_acr_lm <- function(data, outcome, exposure, adjust_acr = TRUE) {
  
  rhs <- c(exposure, "age", "sex")
  if (adjust_acr && !outcome %in% c("acr_u", "log_acr_u")) {
   rhs <- c(rhs, "log_acr_u_bl")
  }
  
  m0 <- as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + ")))
  
  model_data <- data %>%
    dplyr::select(all.vars(m0)) %>%
    tidyr::drop_na()
  
  m1 <- lm(m0, data = model_data)
  
  list(
    model = m1,
    tidy = broom::tidy(m1) %>%
      dplyr::mutate(
        analysis = "baseline_acr_cross_sectional",
        outcome = outcome,
        exposure = exposure,
        model_formula = paste(deparse(m0), collapse = " ")
      )
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
baseline_acr_results_list <- lapply(baseline_model_outcomes, function(outcome) {
  lapply(pfas_group_vars, function(exposure) {
    run_baseline_pfas_acr_lm(dat_baseline, outcome, exposure)
  })
})

baseline_acr_results <- dplyr::bind_rows(
  lapply(baseline_acr_results_list, function(x) dplyr::bind_rows(lapply(x, `[[`, "tidy")))
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
baseline_acr_pfas_only <- baseline_acr_results %>%
  dplyr::mutate(
    exposure_clean = stringr::str_replace_all(exposure, "`", "")
  ) %>%
  dplyr::filter(term %in% c(exposure, exposure_clean)) %>%
  dplyr::group_by(analysis, outcome) %>%
  dplyr::mutate(
    p_fdr = p.adjust(p.value, method = "fdr"),
    sig_fdr = p_fdr < 0.05,
    p_fmt = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    p_fdr_fmt = ifelse(p_fdr < 0.001, "<0.001", sprintf("%.3f", p_fdr))
  ) %>%
  dplyr::ungroup()

baseline_acr_results_final <- baseline_acr_results %>%
  dplyr::left_join(
    baseline_acr_pfas_only %>%
      dplyr::select(analysis, outcome, exposure, term, p_fdr, sig_fdr, p_fdr_fmt),
    by = c("analysis", "outcome", "exposure", "term")
  )

baseline_acr_table <- baseline_acr_results_final %>%
  dplyr::filter(term %in% c(exposure, stringr::str_replace_all(exposure, "`", ""))) %>%
  dplyr::mutate(
    estimate_ci = paste0(
      round(estimate, 3), " (",
      round(estimate - 1.96 * std.error, 3), ", ",
      round(estimate + 1.96 * std.error, 3), ")"
    ),
    p_fmt = fmt_p(p.value)
  ) %>%
  dplyr::select(
    analysis, outcome, exposure, estimate, std.error, statistic, p.value,
    p_fdr, sig_fdr, estimate_ci, p_fmt, p_fdr_fmt
  )

baseline_acr_table


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
forest_baseline_acr_df <- baseline_acr_table %>%
  dplyr::mutate(
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_list <- list()


for (outcome in unique(forest_baseline_acr_df$outcome)) {
  df <- forest_baseline_acr_df %>%
    dplyr::filter(outcome == !!outcome)
  p <- ggplot(df, aes(x = estimate, y = exposure)) +
  geom_point(aes(colour = sig_fdr), size = 2.5) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "plum")
) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = paste("PFAS associations with", get_outcome_label(outcome)),
    x = "Effect size (β)",
    y = "PFAS",
    shape = "FDR <0.05"
  )
plot_list[[outcome]] <- p
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_list

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
run_baseline_pfas_tanner_lm <- function(data, outcome, exposure, adjust_acr = TRUE) {
  
  rhs <- c(exposure, "age", "sex", "tan_stage")
  
  
  m0 <- as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + ")))
  
  model_data <- data %>%
    dplyr::select(all.vars(m0)) %>%
    tidyr::drop_na()
  
  m1 <- lm(m0, data = model_data)
  
  list(
    model = m1,
    tidy = broom::tidy(m1) %>%
      dplyr::mutate(
        analysis = "baseline_cross_sectional",
        outcome = outcome,
        exposure = exposure,
        model_formula = paste(deparse(m0), collapse = " ")
      )
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
baseline_tanner_results_list <- lapply(baseline_model_outcomes, function(outcome) {
  lapply(pfas_group_vars, function(exposure) {
    run_baseline_pfas_tanner_lm(dat_baseline, outcome, exposure)
  })
})

baseline_tanner_results <- dplyr::bind_rows(
  lapply(baseline_tanner_results_list, function(x) dplyr::bind_rows(lapply(x, `[[`, "tidy")))
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
baseline_tanner_pfas_only <- baseline_tanner_results %>%
  dplyr::mutate(
    exposure_clean = stringr::str_replace_all(exposure, "`", "")
  ) %>%
  dplyr::filter(term %in% c(exposure, exposure_clean)) %>%
  dplyr::group_by(analysis, outcome) %>%
  dplyr::mutate(
    p_fdr = p.adjust(p.value, method = "fdr"),
    sig_fdr = p_fdr < 0.05,
    p_fmt = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    p_fdr_fmt = ifelse(p_fdr < 0.001, "<0.001", sprintf("%.3f", p_fdr))
  ) %>%
  dplyr::ungroup()

baseline_tanner_results_final <- baseline_tanner_results %>%
  dplyr::left_join(
    baseline_tanner_pfas_only %>%
      dplyr::select(analysis, outcome, exposure, term, p_fdr, sig_fdr, p_fdr_fmt),
    by = c("analysis", "outcome", "exposure", "term")
  )

baseline_tanner_table <- baseline_tanner_results_final %>%
  dplyr::filter(term %in% c(exposure, stringr::str_replace_all(exposure, "`", ""))) %>%
  dplyr::mutate(
    estimate_ci = paste0(
      round(estimate, 3), " (",
      round(estimate - 1.96 * std.error, 3), ", ",
      round(estimate + 1.96 * std.error, 3), ")"
    ),
    p_fmt = fmt_p(p.value)
  ) %>%
  dplyr::select(
    analysis, outcome, exposure, estimate, std.error, statistic, p.value,
    p_fdr, sig_fdr, estimate_ci, p_fmt, p_fdr_fmt
  )

baseline_tanner_table


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
forest_baseline_tanner_df <- baseline_tanner_table %>%
  dplyr::mutate(
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_list <- list()


for (outcome in unique(forest_baseline_tanner_df$outcome)) {
  df <- forest_baseline_tanner_df %>%
    dplyr::filter(outcome == !!outcome)
  p <- ggplot(df, aes(x = estimate, y = exposure)) +
  geom_point(aes(colour = sig_fdr), size = 2.5) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "plum")
) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = paste("PFAS associations with", get_outcome_label(outcome)),
    x = "Effect size (β)",
    y = "PFAS",
    shape = "FDR <0.05"
  )
plot_list[[outcome]] <- p
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_list


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
#trying to mimic hailey's code for longitudinal data
# Standardize outcome variables so betas are comparable across outcomes
dat_long <- dat_long %>%
  mutate(across(all_of(group_time_outcomes), ~ scale(.)[,1],
                .names = "{.col}_z"))

# Update outcome_vars_long to use standardized versions
group_time_outcomes2 <- paste0(group_time_outcomes, "_z")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
model_3a <- function(data, outcome, exposure) {
  
  # -----------------------------
  # age, sex, and group as covariates
  # -----------------------------
  covars <- c("age", "sex", "group")
  
  # -----------------------------
  # Keep needed columns only
  # -----------------------------
  needed_vars <- c("record_id", outcome, exposure, "visit_num", covars)
  
  model_data <- data %>%
    dplyr::select(dplyr::all_of(needed_vars)) %>%
    tidyr::drop_na()
  
  # -----------------------------
  # Outlier removal per outcome
  # -----------------------------
  m_out <- mean(model_data[[outcome]], na.rm = TRUE)
  s_out <- sd(model_data[[outcome]], na.rm = TRUE)
  
  removed_ids <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) > 3 * s_out) %>%
    dplyr::pull(record_id) %>%
    unique()
  
  model_data <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) <= 3 * s_out)
  
  # -----------------------------
  # Sample size info
  # -----------------------------
  n_subj <- dplyr::n_distinct(model_data$record_id)
  n_obs  <- nrow(model_data)
  
  # -----------------------------
  # Formula
  # -----------------------------
  rhs <- c(exposure, "visit_num", paste0(exposure, ":visit_num"), covars)
  
  m0 <- as.formula(paste0(
    outcome, " ~ ",
    paste(rhs, collapse = " + "),
    " + (1 + visit_num | record_id)"
  ))
  
  # -----------------------------
  # Fit model with fallback
  # -----------------------------
  m1 <- tryCatch(
    lmerTest::lmer(m0, data = model_data, REML = TRUE),
    error = function(e) {
      fallback_formula <- as.formula(
        paste0(
          outcome, " ~ ",
          paste(rhs, collapse = " + "),
          " + (1 | record_id)"
        )
      )
      lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
    }
  )
  
  if (lme4::isSingular(m1, tol = 1e-4)) {
    fallback_formula <- as.formula(
      paste0(
        outcome, " ~ ",
        paste(rhs, collapse = " + "),
        " + (1 | record_id)"
      )
    )
    m1 <- lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
  }
  
  # -----------------------------
  # Tidy output
  # -----------------------------
  tidy_tbl <- broom.mixed::tidy(m1, effects = "fixed") %>%
    dplyr::mutate(
      analysis = "baseline_pfas_longitudinal_base",
      outcome = outcome,
      exposure = exposure,
      n_subj = n_subj,
      n_obs = n_obs,
      removed_ids = if (length(removed_ids) > 0) paste(removed_ids, collapse = ", ") else NA_character_,
      model_formula = paste(deparse(formula(m1)), collapse = " ")
    )
  
  list(
    model = m1,
    tidy = tidy_tbl
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod3a_list <- lapply(group_time_outcomes2, function(outcome) {
  lapply(pfas_group_vars, function(exposure) {
    model_3a(dat_long, outcome, exposure)
  })
})
names(mod3a_list) <- group_time_outcomes2

mod3a_results <- dplyr::bind_rows(
  lapply(mod3a_list, function(x) {
    dplyr::bind_rows(lapply(x, `[[`, "tidy"))
  })
)

mod3a_results


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod3a_results_fdr <- mod3a_results %>%
  dplyr::mutate(
    term_type = dplyr::case_when(
      term == exposure ~ "main",
      term == paste0(exposure, ":visit_num") ~ "interaction",
      term == paste0("visit_num:", exposure) ~ "interaction",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::group_by(term_type) %>%
  dplyr::mutate(
    p_fdr = ifelse(!is.na(term_type), p.adjust(p.value, method = "fdr"), NA_real_),
    sig_fdr = dplyr::case_when(
      !is.na(p_fdr) & p_fdr < 0.001 ~ "***",
      !is.na(p_fdr) & p_fdr < 0.01  ~ "**",
      !is.na(p_fdr) & p_fdr < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::ungroup()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
hm_mod3a <- mod3a_results_fdr %>%
  dplyr::group_by(outcome, exposure) %>%
  dplyr::summarise(
    beta_main = estimate[term == exposure][1],
    pval_main = p.value[term == exposure][1],
    pval_main_fdr = p_fdr[term == exposure][1],
    beta_int = estimate[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    pval_int = p.value[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    pval_int_fdr = p_fdr[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    beta_visit0 = beta_main,
    beta_visit1 = beta_main + beta_int,
    beta_visit2 = beta_main + 2 * beta_int
  ) %>%
  dplyr::select(
    outcome, exposure,
    beta_visit0, beta_visit1, beta_visit2,
    pval_main, pval_main_fdr,
    pval_int, pval_int_fdr
  ) %>%
  tidyr::pivot_longer(
    cols = c(beta_visit0, beta_visit1, beta_visit2),
    names_to = "visit",
    values_to = "beta"
  ) %>%
  dplyr::mutate(
    visit = factor(
      visit,
      levels = c("beta_visit0", "beta_visit1", "beta_visit2"),
      labels = c("Baseline", "Year 1", "Year 2")
    ),
    sig = dplyr::case_when(
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.001 ~ "***",
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.01  ~ "**",
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.05  ~ "*",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.001 ~ "***",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.01  ~ "**",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.05  ~ "*",
      TRUE ~ ""
    ),
    exposure = factor(
      exposure,
      levels = pfas_group_vars,
      labels = sapply(pfas_group_vars, get_pfas_label)
    ),
    outcome = factor(
      outcome,
      levels = group_time_outcomes2,
      labels = sapply(group_time_outcomes2, get_outcome_label)
    )
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
beta_mod3a <- max(abs(hm_mod3a$beta), na.rm = TRUE)

plot_model_3a<- ggplot(
  hm_mod3a,
  aes(x = exposure, y = outcome, fill = beta)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), size = 3.5, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-beta_mod3a, beta_mod3a),
    name = "Beta"
  ) +
  facet_wrap(~ visit, nrow = 1) +
  labs(
    title = "Baseline PFAS associations with longitudinal outcomes",
    subtitle = paste0(
      "Base model (Age+Sex)\n",
      "* FDR<0.05  ** FDR<0.01  *** FDR<0.001"
    ),
    x = "Baseline PFAS",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

plot_model_3a

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
model_3b <- function(data, outcome, exposure) {
  
  # -----------------------------
  # Build covariate list
  # -----------------------------
  covars <- c("age", "sex", "group")
  
  if (!outcome %in% c("acr_u", "log_acr_u")) {
    covars <- c(covars, "log_acr_u_bl")
  }
  
  # -----------------------------
  # Keep needed columns only
  # -----------------------------
  needed_vars <- c("record_id", outcome, exposure, "visit_num", covars)
  
  model_data <- data %>%
    dplyr::select(dplyr::all_of(needed_vars)) %>%
    tidyr::drop_na()
  
  # -----------------------------
  # Outlier removal per outcome
  # -----------------------------
  m_out <- mean(model_data[[outcome]], na.rm = TRUE)
  s_out <- sd(model_data[[outcome]], na.rm = TRUE)
  
  removed_ids <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) > 3 * s_out) %>%
    dplyr::pull(record_id) %>%
    unique()
  
  model_data <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) <= 3 * s_out)
  
  # -----------------------------
  # Sample size info
  # -----------------------------
  n_subj <- dplyr::n_distinct(model_data$record_id)
  n_obs  <- nrow(model_data)
  
  # -----------------------------
  # Formula
  # -----------------------------
  rhs <- c(exposure, "visit_num", paste0(exposure, ":visit_num"), covars)
  
  m0 <- as.formula(paste0(
    outcome, " ~ ",
    paste(rhs, collapse = " + "),
    " + (1 + visit_num | record_id)"
  ))
  
  # -----------------------------
  # Fit model with fallback
  # -----------------------------
  m1 <- tryCatch(
    lmerTest::lmer(m0, data = model_data, REML = TRUE),
    error = function(e) {
      fallback_formula <- as.formula(
        paste0(
          outcome, " ~ ",
          paste(rhs, collapse = " + "),
          " + (1 | record_id)"
        )
      )
      lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
    }
  )
  
  if (lme4::isSingular(m1, tol = 1e-4)) {
    fallback_formula <- as.formula(
      paste0(
        outcome, " ~ ",
        paste(rhs, collapse = " + "),
        " + (1 | record_id)"
      )
    )
    m1 <- lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
  }
  
  # -----------------------------
  # Tidy output
  # -----------------------------
  tidy_tbl <- broom.mixed::tidy(m1, effects = "fixed") %>%
    dplyr::mutate(
      analysis = "baseline_pfas_longitudinal_acr",
      outcome = outcome,
      exposure = exposure,
      n_subj = n_subj,
      n_obs = n_obs,
      removed_ids = if (length(removed_ids) > 0) paste(removed_ids, collapse = ", ") else NA_character_,
      model_formula = paste(deparse(formula(m1)), collapse = " ")
    )
  
  list(
    model = m1,
    tidy = tidy_tbl
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod3b_list <- lapply(group_time_outcomes2, function(outcome) {
  lapply(pfas_group_vars, function(exposure) {
    model_3b(dat_long, outcome, exposure)
  })
})
names(mod3b_list) <- group_time_outcomes2

mod_3b_results <- dplyr::bind_rows(
  lapply(mod3b_list, function(x) {
    dplyr::bind_rows(lapply(x, `[[`, "tidy"))
  })
)

mod_3b_results


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod3b_results_fdr <- mod_3b_results %>%
  dplyr::mutate(
    term_type = dplyr::case_when(
      term == exposure ~ "main",
      term == paste0(exposure, ":visit_num") ~ "interaction",
      term == paste0("visit_num:", exposure) ~ "interaction",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::group_by(term_type) %>%
  dplyr::mutate(
    p_fdr = ifelse(!is.na(term_type), p.adjust(p.value, method = "fdr"), NA_real_),
    sig_fdr = dplyr::case_when(
      !is.na(p_fdr) & p_fdr < 0.001 ~ "***",
      !is.na(p_fdr) & p_fdr < 0.01  ~ "**",
      !is.na(p_fdr) & p_fdr < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::ungroup()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
hm_mod3b <- mod3b_results_fdr %>%
  dplyr::group_by(outcome, exposure) %>%
  dplyr::summarise(
    beta_main = estimate[term == exposure][1],
    pval_main = p.value[term == exposure][1],
    pval_main_fdr = p_fdr[term == exposure][1],
    beta_int = estimate[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    pval_int = p.value[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    pval_int_fdr = p_fdr[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    beta_visit0 = beta_main,
    beta_visit1 = beta_main + beta_int,
    beta_visit2 = beta_main + 2 * beta_int
  ) %>%
  dplyr::select(
    outcome, exposure,
    beta_visit0, beta_visit1, beta_visit2,
    pval_main, pval_main_fdr,
    pval_int, pval_int_fdr
  ) %>%
  tidyr::pivot_longer(
    cols = c(beta_visit0, beta_visit1, beta_visit2),
    names_to = "visit",
    values_to = "beta"
  ) %>%
  dplyr::mutate(
    visit = factor(
      visit,
      levels = c("beta_visit0", "beta_visit1", "beta_visit2"),
      labels = c("Baseline", "Year 1", "Year 2")
    ),
    sig = dplyr::case_when(
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.001 ~ "***",
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.01  ~ "**",
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.05  ~ "*",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.001 ~ "***",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.01  ~ "**",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.05  ~ "*",
      TRUE ~ ""
    ),
    exposure = factor(
      exposure,
      levels = pfas_group_vars,
      labels = sapply(pfas_group_vars, get_pfas_label)
    ),
    outcome = factor(
      outcome,
      levels = group_time_outcomes2,
      labels = sapply(group_time_outcomes2, get_outcome_label)
    )
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
beta_mod3b <- max(abs(hm_mod3b$beta), na.rm = TRUE)

plot_mod3b <- ggplot(
  hm_mod3b,
  aes(x = exposure, y = outcome, fill = beta)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), size = 3.5, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-beta_mod3b, beta_mod3b),
    name = "Beta"
  ) +
  facet_wrap(~ visit, nrow = 1) +
  labs(
    title = "Baseline PFAS associations with longitudinal outcomes",
    subtitle = paste0(
      "Baseline uACR-adjusted mixed models.\n",
      "* FDR<0.05  ** FDR<0.01  *** FDR<0.001"
    ),
    x = "Baseline PFAS",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

plot_mod3b

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
model_3c <- function(data, outcome, exposure) {
  
  # -----------------------------
  # Build covariate list
  # -----------------------------
  covars <- c("age", "sex", "group", "tan_stage")
  
  # -----------------------------
  # Keep needed columns only
  # -----------------------------
  needed_vars <- c("record_id", outcome, exposure, "visit_num", covars)
  
  model_data <- data %>%
    dplyr::select(dplyr::all_of(needed_vars)) %>%
    tidyr::drop_na()
  
  # -----------------------------
  # Outlier removal per outcome
  # -----------------------------
  m_out <- mean(model_data[[outcome]], na.rm = TRUE)
  s_out <- sd(model_data[[outcome]], na.rm = TRUE)
  
  removed_ids <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) > 3 * s_out) %>%
    dplyr::pull(record_id) %>%
    unique()
  
  model_data <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) <= 3 * s_out)
  
  # -----------------------------
  # Sample size info
  # -----------------------------
  n_subj <- dplyr::n_distinct(model_data$record_id)
  n_obs  <- nrow(model_data)
  
  # -----------------------------
  # Formula
  # -----------------------------
  rhs <- c(exposure, "visit_num", paste0(exposure, ":visit_num"), covars)
  
  m0 <- as.formula(paste0(
    outcome, " ~ ",
    paste(rhs, collapse = " + "),
    " + (1 + visit_num | record_id)"
  ))
  
  
  # -----------------------------
  # Fit model with fallback
  # -----------------------------
  m1 <- tryCatch(
    lmerTest::lmer(m0, data = model_data, REML = TRUE),
    error = function(e) {
      fallback_formula <- as.formula(
        paste0(
          outcome, " ~ ",
          paste(rhs, collapse = " + "),
          " + (1 | record_id)"
        )
      )
      lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
    }
  )
  
  if (lme4::isSingular(m1, tol = 1e-4)) {
    fallback_formula <- as.formula(
      paste0(
        outcome, " ~ ",
        paste(rhs, collapse = " + "),
        " + (1 | record_id)"
      )
    )
    m1 <- lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
  }
  
  # -----------------------------
  # Tidy output
  # -----------------------------
  tidy_tbl <- broom.mixed::tidy(m1, effects = "fixed") %>%
    dplyr::mutate(
      analysis = "baseline_pfas_longitudinal_tanner",
      outcome = outcome,
      exposure = exposure,
      n_subj = n_subj,
      n_obs = n_obs,
      removed_ids = if (length(removed_ids) > 0) paste(removed_ids, collapse = ", ") else NA_character_,
      model_formula = paste(deparse(formula(m1)), collapse = " ")
    )
  
  list(
    model = m1,
    tidy = tidy_tbl
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod3c_list<- lapply(group_time_outcomes2, function(outcome) {
  lapply(pfas_group_vars, function(exposure) {
    model_3c(dat_long, outcome, exposure)
  })
})
names(mod3c_list) <- group_time_outcomes2

mod3c_results <- dplyr::bind_rows(
  lapply(mod3c_list, function(x) {
    dplyr::bind_rows(lapply(x, `[[`, "tidy"))
  })
)

mod3c_results


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod3c_results_fdr <- mod3c_results %>%
  dplyr::mutate(
    term_type = dplyr::case_when(
      term == exposure ~ "main",
      term == paste0(exposure, ":visit_num") ~ "interaction",
      term == paste0("visit_num:", exposure) ~ "interaction",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::group_by(term_type) %>%
  dplyr::mutate(
    p_fdr = ifelse(!is.na(term_type), p.adjust(p.value, method = "fdr"), NA_real_),
    sig_fdr = dplyr::case_when(
      !is.na(p_fdr) & p_fdr < 0.001 ~ "***",
      !is.na(p_fdr) & p_fdr < 0.01  ~ "**",
      !is.na(p_fdr) & p_fdr < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::ungroup()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------
# Heatmap of Tanner-adjusted longitudinal outcomes
#-----------------------------------
hm_mod3c <- mod3c_results_fdr %>%
  dplyr::group_by(outcome, exposure) %>%
  dplyr::summarise(
    beta_main = estimate[term == exposure][1],
    pval_main = p.value[term == exposure][1],
    pval_main_fdr = p_fdr[term == exposure][1],
    beta_int = estimate[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    pval_int = p.value[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    pval_int_fdr = p_fdr[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    beta_visit0 = beta_main,
    beta_visit1 = beta_main + beta_int,
    beta_visit2 = beta_main + 2 * beta_int
  ) %>%
  dplyr::select(
    outcome, exposure,
    beta_visit0, beta_visit1, beta_visit2,
    pval_main, pval_main_fdr,
    pval_int, pval_int_fdr
  ) %>%
  tidyr::pivot_longer(
    cols = c(beta_visit0, beta_visit1, beta_visit2),
    names_to = "visit",
    values_to = "beta"
  ) %>%
  dplyr::mutate(
    visit = factor(
      visit,
      levels = c("beta_visit0", "beta_visit1", "beta_visit2"),
      labels = c("Baseline", "Year 1", "Year 2")
    ),
    sig = dplyr::case_when(
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.001 ~ "***",
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.01  ~ "**",
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.05  ~ "*",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.001 ~ "***",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.01  ~ "**",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.05  ~ "*",
      TRUE ~ ""
    ),
    exposure = factor(
      exposure,
      levels = pfas_group_vars,
      labels = sapply(pfas_group_vars, get_pfas_label)
    ),
    outcome = factor(
      outcome,
      levels = group_time_outcomes2,
      labels = sapply(group_time_outcomes2, get_outcome_label)
    )
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
beta_mod3c <- max(abs(hm_mod3c$beta), na.rm = TRUE)

plot_mod3c <- ggplot(
  hm_mod3c,
  aes(x = exposure, y = outcome, fill = beta)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), size = 3.5, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-beta_mod3c, beta_mod3c),
    name = "Beta"
  ) +
  facet_wrap(~ visit, nrow = 1) +
  labs(
    title = "Baseline PFAS associations with longitudinal outcomes",
    subtitle = paste0(
      "Tanner-adjusted mixed models.\n",
      "* FDR<0.05  ** FDR<0.01  *** FDR<0.001"
    ),
    x = "Baseline PFAS",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

plot_mod3c


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------
# baseline PFAS predicting longitudinal outcomes,
# adjusted for age, sex, group, baseline uACR, and Tanner stage
#------------------------------------------------------------------
model_3d <- function(data, outcome, exposure) {
  
  # -----------------------------
  # Build covariate list
  # -----------------------------
  covars <- c("age", "sex", "group", "tan_stage")
  
  # Add baseline uACR except when it's the outcome
  if (!outcome %in% c("acr_u", "log_acr_u")) {
    covars <- c(covars, "log_acr_u_bl")
  }
  
  # -----------------------------
  # Keep needed columns only
  # -----------------------------
  needed_vars <- c("record_id", outcome, exposure, "visit_num", covars)
  
  model_data <- data %>%
    dplyr::select(dplyr::all_of(needed_vars)) %>%
    tidyr::drop_na()
  
  # -----------------------------
  # Outlier removal per outcome
  # -----------------------------
  m_out <- mean(model_data[[outcome]], na.rm = TRUE)
  s_out <- sd(model_data[[outcome]], na.rm = TRUE)
  
  removed_ids <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) > 3 * s_out) %>%
    dplyr::pull(record_id) %>%
    unique()
  
  model_data <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) <= 3 * s_out)
  
  # -----------------------------
  # Sample size info
  # -----------------------------
  n_subj <- dplyr::n_distinct(model_data$record_id)
  n_obs  <- nrow(model_data)
  
  # -----------------------------
  # Formula
  # -----------------------------
  rhs <- c(exposure, "visit_num", paste0(exposure, ":visit_num"), covars)
  
  m0 <- as.formula(paste0(
    outcome, " ~ ",
    paste(rhs, collapse = " + "),
    " + (1 + visit_num | record_id)"
  ))
  
  
  # -----------------------------
  # Fit model with fallback
  # -----------------------------
  m1 <- tryCatch(
    lmerTest::lmer(m0, data = model_data, REML = TRUE),
    error = function(e) {
      fallback_formula <- as.formula(
        paste0(
          outcome, " ~ ",
          paste(rhs, collapse = " + "),
          " + (1 | record_id)"
        )
      )
      lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
    }
  )
  
  if (lme4::isSingular(m1, tol = 1e-4)) {
    fallback_formula <- as.formula(
      paste0(
        outcome, " ~ ",
        paste(rhs, collapse = " + "),
        " + (1 | record_id)"
      )
    )
    m1 <- lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
  }
  
  # -----------------------------
  # Tidy output
  # -----------------------------
  tidy_tbl <- broom.mixed::tidy(m1, effects = "fixed") %>%
    dplyr::mutate(
      analysis = "baseline_pfas_longitudinal_acr_tanner",
      outcome = outcome,
      exposure = exposure,
      n_subj = n_subj,
      n_obs = n_obs,
      removed_ids = if (length(removed_ids) > 0) paste(removed_ids, collapse = ", ") else NA_character_,
      model_formula = paste(deparse(formula(m1)), collapse = " ")
    )
  
  list(
    model = m1,
    tidy = tidy_tbl
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod3d_list<- lapply(group_time_outcomes2, function(outcome) {
  lapply(pfas_group_vars, function(exposure) {
    model_3d(dat_long, outcome, exposure)
  })
})
names(mod3d_list) <- group_time_outcomes2

mod3d_results<- dplyr::bind_rows(
  lapply(mod3d_list, function(x) {
    dplyr::bind_rows(lapply(x, `[[`, "tidy"))
  })
)

mod3d_results


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod3d_results_fdr<- mod3d_results %>%
  dplyr::mutate(
    term_type = dplyr::case_when(
      term == exposure ~ "main",
      term == paste0(exposure, ":visit_num") ~ "interaction",
      term == paste0("visit_num:", exposure) ~ "interaction",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::group_by(term_type) %>%
  dplyr::mutate(
    p_fdr = ifelse(!is.na(term_type), p.adjust(p.value, method = "fdr"), NA_real_),
    sig_fdr = dplyr::case_when(
      !is.na(p_fdr) & p_fdr < 0.001 ~ "***",
      !is.na(p_fdr) & p_fdr < 0.01  ~ "**",
      !is.na(p_fdr) & p_fdr < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::ungroup()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
hm_mod3d <- mod3d_results_fdr %>%
  dplyr::group_by(outcome, exposure) %>%
  dplyr::summarise(
    
    # Main PFAS effect (baseline)
    beta_main = estimate[term == exposure][1],
    pval_main = p.value[term == exposure][1],
    pval_main_fdr = p_fdr[term == exposure][1],
    
    # Interaction (change over time)
    beta_int = estimate[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    
    pval_int = p.value[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    
    pval_int_fdr = p_fdr[
      term == paste0(exposure, ":visit_num") |
      term == paste0("visit_num:", exposure)
    ][1],
    
    .groups = "drop"
  ) %>%
  
  # Convert to visit-specific effects
  dplyr::mutate(
    beta_visit0 = beta_main,
    beta_visit1 = beta_main + beta_int,
    beta_visit2 = beta_main + 2 * beta_int
  ) %>%
  
  dplyr::select(
    outcome, exposure,
    beta_visit0, beta_visit1, beta_visit2,
    pval_main, pval_main_fdr,
    pval_int, pval_int_fdr
  ) %>%
  
  # Long format for plotting
  tidyr::pivot_longer(
    cols = c(beta_visit0, beta_visit1, beta_visit2),
    names_to = "visit",
    values_to = "beta"
  ) %>%
  
  dplyr::mutate(
    
    visit = factor(
      visit,
      levels = c("beta_visit0", "beta_visit1", "beta_visit2"),
      labels = c("Baseline", "Year 1", "Year 2")
    ),
    
    # Significance stars
    sig = dplyr::case_when(
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.001 ~ "***",
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.01  ~ "**",
      visit == "Baseline" & !is.na(pval_main_fdr) & pval_main_fdr < 0.05  ~ "*",
      
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.001 ~ "***",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.01  ~ "**",
      visit %in% c("Year 1", "Year 2") & !is.na(pval_int_fdr) & pval_int_fdr < 0.05  ~ "*",
      
      TRUE ~ ""
    ),
    
    exposure = factor(
      exposure,
      levels = pfas_group_vars,
      labels = sapply(pfas_group_vars, get_pfas_label)
    ),
    
    outcome = factor(
      outcome,
      levels = group_time_outcomes2,
      labels = sapply(group_time_outcomes2, get_outcome_label)
    )
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------

beta_mod3d <- max(abs(hm_mod3d$beta), na.rm = TRUE)

plot_mod3d <- ggplot(
  hm_mod3d,
  aes(x = exposure, y = outcome, fill = beta)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), size = 3.5, color = "black") +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-beta_mod3d, beta_mod3d),
    name = "Beta"
  ) +
  
  facet_wrap(~ visit, nrow = 1) +
  
  labs(
    title = "Baseline PFAS associations with longitudinal outcomes",
    subtitle = paste0(
      "uACR + Tanner-adjusted mixed models\n",
      "* FDR<0.05  ** FDR<0.01  *** FDR<0.001"
    ),
    x = "Baseline PFAS",
    y = NULL
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

plot_mod3d






