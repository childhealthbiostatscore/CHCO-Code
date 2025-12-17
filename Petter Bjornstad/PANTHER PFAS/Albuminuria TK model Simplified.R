# Initial Analysis

source(fs::path(here::here("!libraries.R")))
source(fs::path(here::here("!directories.R")))
source(fs::path(here::here("!functions.R")))
library(broom)
library(tidyverse)

# PFOA constants -----
vd_pfoa = 0.2 #Vd (PFOA)	0.2
half_life_pfoa = 3.2 #t1/2 (PFOA)	3.2

# 1. Renal-HEIR data ----------
# Read in data
rhdf <- readxl::read_xlsx(fs::path(dir_cleaned_data,
                                   "RH_Complete_data_imputed_lg2_V5.xlsx")) |>
  filter(!is.na(pfas_t_pfda))

rhdf <- rhdf |> 
  rename(year = date_year)

# Select key variables from analysis
dat <- rhdf |>
  tidylog::select(-contains("met_"), 
                  -contains("prot_"), 
                  -contains("lipid_")) 

## 1.1 Identifycalculate key variables -------

# 1.1. Body surface area (m^2): 
dat$bsa_dubois

## 1.2. Unindexed eGFR (mL/min) ----
#  (two options- direct vs. indirectly measured eGFR): 
dat$eGFRunindexed = dat$bsa_dubois * dat$eGFR_CKD_epi/1.73
dat$mGFRunindexed = dat$gfr_raw_plasma 
# dat$mGFRunindexed = dat$gfr_bsa_plasma

## 1.3. uACR (mg/g, already calculated) ----- 
dat$acr_u

## 1.4. Daily creatinine excretion (units out: ) ----
# creatinine_s: mg/dL
dat$cr_daily = 14.4 * dat$mGFRunindexed * dat$creatinine_s 

## 1.5. Albumin excretion rate (units out: mg/day) ----
# aer_4_coltime units: mcg * mL / min
# dat$aer_4_mg_ml_min =  /1000
dat$aer_mg_day = dat$aer_4_coltime * 1.44

## 1.6. Fraction bound of PFAS ----
#assume 99.92 (based on low variance in NHANES data): 
# K_alb_w = 10^4.48
# dat$fbound = 1 / (1 + 1/K_alb_w * 1/dat$))
fbound = 0.9992  

## 1.7. Fraction unbound of PFAS: ---------
funbound = 1 - fbound

## 1.8. Baseline clearance: ----
dat$CLbase_half_life_vd = (log(2)*vd_pfoa*dat$weight)/(half_life_pfoa*365) * 1000
dat$CLbase_gfr = dat$eGFRunindexed * funbound * 1400

## 1.9. Estimate serum albumin concentration ----
dat$Calb = dat$tot_protein/1.65 #(From NHANES PMID: 40661684)

## 1.10. Albumin clearance: CLalb = AER/Calb × funbound ----
dat$CLalb = (dat$aer_mg_day/dat$Calb) * funbound * 1000

# CLtotal = CLbase + Clalb
dat$CLtotal = dat$CLalb + dat$CLbase_gfr

# % albumin = CLalb / CLbase × 100
dat$pct_alb = 100*dat$CLalb/dat$CLbase_gfr


## 1.11. Clean results --------
# Transform vars and create names similar to NHANES
dat <- dat |> 
  mutate(
    log10_cl_alb = log10(CLalb+0.001), 
    log10_cl_base = log10(CLbase_gfr+ 0.001)
  )



## 1.12. Test associations ------------
mod <- dat |> 
  tidylog::filter(!is.na(pfas_t_pfoa), 
                  !is.na(log10_cl_alb)) %>%
  lm(log(pfas_t_pfoa) ~ log10_cl_alb + log10_cl_base + age, .)
summary(mod)


# 2. NHANES data --------
## 2.1. Read in data -------
nhanes <- readxl::read_excel("/Users/jagoodri/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/Goodrich Projects Active/1_proj_dir/NHANES/0_data/data_from_Fabian_Fischer/NHANES data analyzed_cleaned.xlsx") |>
  janitor::clean_names()

nhanes <- nhanes |> rename(log10_cl_alb = log10_cl_alb)
nhanes$log10_cl_base = log10(nhanes$cl_base_m_l_day_from_e_gfr)

nhanes$hispanic = if_else(nhanes$race %in% c("Mexican American", "Other Hispanic"), "Hispanic", "Not Hispanic")

nhanes$year = nhanes$numerical

nhanes_filtered <- nhanes |> 
  tidylog::filter(
    !is.na(log10_cl_alb), 
    !is.na(log10_cl_base))

mean(nhanes_filtered$log10_cl_alb)


## 2.2 Run the models -----
covars = "year + age"

# Model with albumin
mod_full <- lm(
  as.formula(paste("log(pfoa_ng_m_l) ~ log10_cl_alb + log10_cl_base + ", covars)),
  data = nhanes)

summary(mod_full)

# Model without albumin
mod_no_alb <- lm(
  as.formula(paste("log(pfoa_ng_m_l) ~ log10_cl_base + ", covars)),
  # as.formula(paste("log(pfoa_ng_m_l) ~ log10_cl_base")), 
  data = nhanes)


anova(mod_full, mod_no_alb)
summary(mod_full)
summary(mod_no_alb)

# 3. Estimate Renal-HEIR PFOA based on NHANES TK model ----
# Get prediction dataset
pred_dat <- dat |> 
  tidylog::filter(!is.na(log10_cl_alb), !is.na(log10_cl_base))

# Play around with how year is modeled- for this, decided to just use the raw data
# pred_dat$year = mean(dat$year, na.rm=T)
# pred_dat$year = 2020

## 3.1. Predictions -----
# Main Model (without race)
pred_dat$log_pfoa_est_w_alb = predict(mod_full, newdata = pred_dat) |> as.numeric()

# Model without albumin
pred_dat$log_pfoa_est_no_alb = predict(mod_no_alb, newdata = pred_dat) |> as.numeric()

# # Model with race (performs worse than model with race) 
# dat$log_pfoa_est_w_race = predict(mod_race, newdata = as.data.frame(dat)) |> as.numeric()
# # Model with Hispanic ethnicity (performs worse than both other models) 
# dat$log_pfoa_est_w_hisp = predict(mod_hispanic, newdata = as.data.frame(dat)) |> 
#   as.numeric()

## 3.2. Compare predicted vs. measured PFOA ------
# Remove people without any albuminuria
pred_datf <- pred_dat |> 
  tidylog::filter(!is.na(log_pfoa_est_w_alb),
                  aer_mg_day > 0)

# Correlation test- predicted vs. measured
(ct_alb    <- cor.test(pred_datf$log_pfoa_est_w_alb,  log(pred_datf$pfas_t_pfoa)))
(ct_no_alb <- cor.test(pred_datf$log_pfoa_est_no_alb, log(pred_datf$pfas_t_pfoa)))

ct_alb$conf.int[1:2]
ct_no_alb$conf.int[1:2]

# Create nice label
(lab_alb <- sprintf("R² = %.2f, p = %s", unname(ct_alb$estimate)^2, 
                   format.pval(ct_alb$p.value, digits = 3, eps = 1e-3)))
(lab_no_alb <- sprintf("R² = %.2f, p = %s", unname(ct_no_alb$estimate)^2, 
                      format.pval(ct_no_alb$p.value, digits = 3, eps = 1e-3)))


## 3.3. Plot estimated vs. measured PFOA -------
llim = 0.2
ulim = 1.4

(panel_a <- pred_datf |> 
    ggplot(aes(x = exp(log_pfoa_est_w_alb),
               y = pfas_t_pfoa)) + 
    # geom_abline(slope = 1, color = "grey50", linetype = 2) + 
    stat_smooth(method = "lm", color = "grey30") + 
    geom_point(alpha = 0.7) + 
    scale_y_log10(limits = c(llim, ulim)) +
    scale_x_log10(limits = c(llim, ulim)) +
    annotate("text", 
             x = 0.3, 
             y = llim, 
             label = lab_alb,
             hjust = 0,
             vjust = 0.5, 
             size = 6) +
    ylab("Measured Blood PFAS (µg/L)") + 
    xlab("TK Estimated Blood PFAS (µg/L)") +
    cowplot::theme_cowplot() )


(panel_b <- pred_datf |> 
    ggplot(aes(x = exp(log_pfoa_est_no_alb), y = pfas_t_pfoa)) + 
    # geom_abline(slope = 1, color = "grey50", linetype = 2) + 
    stat_smooth(method = "lm", color = "grey30") + 
    geom_point(alpha = 0.7) + 
    scale_y_log10(limits = c(llim, ulim)) +
    scale_x_log10(limits = c(llim, ulim)) +
    annotate("text", x = 0.3, 
             y = llim, 
             label = lab_no_alb,
             hjust = 0,
             vjust = 0.5, 
             size = 6) +
    ylab("Measured Blood PFAS (µg/L)") + 
    xlab("TK Estimated Blood PFAS (µg/L)") +
    cowplot::theme_cowplot() )


## 3.4. Combine panels --------
# two null (blank) panels
blank_panel <- ggplot() + theme_void()

top_row <- cowplot::plot_grid(
  blank_panel, blank_panel,
  labels = c("A. Model Accounting for Albumin Clearance",
             "B. Model Not Accounting for Albumin Clearance"),
  label_size = 12,
  ncol = 2,
  label_x = 0.02,  # far upper left of each blank panel
  label_y = 0.98,
  hjust   = 0,
  vjust   = 1
)

bottom_row <- cowplot::plot_grid(
  panel_a, panel_b,
  ncol = 2
)

cowplot::plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.06, 1)  # adjust height of label row as needed
)


ggsave(filename = fs::path(dir_fig, "Renal_HEIR_tk_model_prediction comparison.jpeg"), 
       width = 10, height = 5, units = "in", dpi = 300)

## 3.5. Compare estimated geometric means -------
pred_datf |> 
  summarise(
    # pfoa_est_w_alb = exp(mean(log_pfoa_est_w_alb)), 
    # pfoa_est_w_no_alb = exp(mean(log_pfoa_est_no_alb)),
    # pfoa_est_w_race_gm = exp(mean(log_pfoa_est_no_alb)), 
    # pfas_t_pfoa_gm = exp(mean(log(pfas_t_pfoa))), 
    RMSE_w_alb2 = sqrt(mean((log_pfoa_est_w_alb-log(pfas_t_pfoa))^2)), 
    RMSE_no_alb = sqrt(mean((log_pfoa_est_no_alb-log(pfas_t_pfoa))^2)), 
    # bias_w_alb  = mean(exp(log_pfoa_est_w_alb) - (pfas_t_pfoa)),
    # bias_no_alb = mean(exp(log_pfoa_est_no_alb) - (pfas_t_pfoa))
  )


# calculate timed urine duration
summary(pred_datf$aer_4_coltime) 

