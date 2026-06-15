### PANTHER multivariable models ###

source("/Users/kristenmiller/Documents/GitHub/CHCO-Code/Petter Bjornstad/PANTHER/panther_baseline_manuscript.qmd")
### outcome variables:

#Table 2:

#add p-value from linear model to this table:
tab.2<-kable(summary(tableby(group_risk ~  total_kidney_volume_ml+ht_adj_tkv+avg_c_r2+avg_k_r2+avg_c_t1+avg_pcascl+avg_c_adc+
                        mm_si+mm_ir+mm_bcell+mm_di+gfr_raw_plasma+gfr_bsa_plasma+erpf_raw_plasma+erpf_bsa_plasma +
                        ra+re+anova(rvr, digits = 3)+glomerular_pressure, 
                      data = panther, 
                      control = table_controls,
                      test = T), digits = 1))

#fit linear models:
outcomes<-c("total_kidney_volume_ml","ht_adj_tkv","avg_c_r2","avg_k_r2","avg_c_t1","avg_pcascl","avg_c_adc","mm_si",
            "mm_ir","mm_bcell","mm_di","gfr_raw_plasma","gfr_bsa_plasma","erpf_raw_plasma","erpf_bsa_plasma",
              "ra","re","rvr","glomerular_pressure",
            "log_lh", "log_estrad", "log_tot_test", "log_free_test", "bl_dheas", "igf_1", "igf1_z_score_calc")

run_adj_lm <- function(outcomes, exposure, data) {
  
  results <- lapply(outcomes, function(outcome) {
    
    formula <- as.formula(paste(outcome, "~ sex +", exposure))
    lm_adj <- lm(formula, data = data)
    
    # Overall F test p-value for the exposure
    anova_p <- round(anova(lm_adj)[exposure, "Pr(>F)"], 3)
    anova_p <- ifelse(anova_p == 0, "<0.001", as.character(anova_p))
    
    data.frame(
      outcome  = outcome,
      exposure = exposure,
      pvalue   = anova_p
    )
  })
  
  do.call(rbind, results)
}

results_group_risk <- run_adj_lm(outcomes = outcomes,
                      exposure = "group_risk",
                      data = panther)

results_tanner_stage <- run_adj_lm(outcomes = outcomes,
                                 exposure = "tanner_stage_comp_panther_cat",
                                 data = panther)


#Table s1: 

### predictors:
sex
total_kidney_volume_ml
sbp
metformin_timepoint
sglti_timepoint
race_ethnicity_condensed
age
hba1c
group #diabetes status