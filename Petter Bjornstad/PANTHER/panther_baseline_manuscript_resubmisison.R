### PANTHER multivariable models ###
#source("/Users/kristenmiller/Documents/GitHub/CHCO-Code/Petter Bjornstad/PANTHER/panther_baseline_manuscript.qmd")

output_dir<-"/Users/kristenmiller/Library/CloudStorage/OneDrive-SharedLibraries-UW/bjornstad_pyle_tommerdahl_lab - Documents/Dissemination/Manuscripts/P031_PANTHER_baseline_yc/Drafts"
setwd(output_dir)


### outcome variables:
outcomes_table2<-c("total_kidney_volume_ml", "ht_adj_tkv","avg_c_r2","avg_k_r2","avg_c_t1","avg_pcascl","avg_c_adc",
                   "gfr_raw_plasma","gfr_bsa_plasma","erpf_raw_plasma","erpf_bsa_plasma",
                   "ra","re","rvr","glomerular_pressure")
outcomes_tableS2<- c("lh","estrad","tot_test","free_test","bl_dheas","igf_1","igf1_z_score_calc","cpeptide",
                     "dexa_lean_mass","dexa_body_fat","dexa_fat_kg","dexa_trunk_mass",
                     "mm_si", "mm_ir","mm_bcell","mm_di")
run_adj_lm <- function(outcomes, exposure, data) {
  
  # Helper to extract F-test p-value, estimate, and 95% CI for exposure
  get_model_stats <- function(lm_fit, exposure) {
    
    p <- round(anova(lm_fit)[exposure, "Pr(>F)"], 3)
    p <- ifelse(p == 0, "<0.001", as.character(p))
    
    coef_rows <- grepl(exposure, rownames(summary(lm_fit)$coefficients))
    est <- round(coef(lm_fit)[coef_rows], 3)
    ci  <- round(confint(lm_fit)[coef_rows, , drop = FALSE], 3)
    
    est_ci <- paste0(est, " (", ci[, 1], ", ", ci[, 2], ")")
    est_ci <- paste(est_ci, collapse = "; ")
    
    list(p = p, est_ci = est_ci)
  }
  
  results <- lapply(outcomes, function(outcome) {
    
    # Model 0: unadjusted
    lm0 <- lm(as.formula(paste(outcome, "~", exposure)), data = data)
    s0  <- get_model_stats(lm0, exposure)
    
    # Model 1: sex adjusted
    lm1 <- lm(as.formula(paste(outcome, "~ sex +", exposure)), data = data)
    s1  <- get_model_stats(lm1, exposure)
    
    # Model 2: sex + TKV adjusted
    lm2 <- lm(as.formula(paste(outcome, "~ sex + total_kidney_volume_ml +", exposure)), data = data)
    s2  <- get_model_stats(lm2, exposure)
    
    # Model 3: sex + TKV + SBP adjusted
    lm3 <- lm(as.formula(paste(outcome, "~ sex + total_kidney_volume_ml + sbp +", exposure)), data = data)
    s3  <- get_model_stats(lm3, exposure)
    
    # Use Hmisc label if available, otherwise fall back to variable name
    outcome_label <- label(data[[outcome]])
    if (is.null(outcome_label) || outcome_label == "") outcome_label <- outcome
    
    # Mean +/- SD by exposure group
    mean_sd <- tapply(data[[outcome]], data[[exposure]], function(x) {
      paste0(round(mean(x, na.rm = TRUE), 3), " (", round(sd(x, na.rm = TRUE), 3),")")
    })
    mean_sd_wide <- as.data.frame(t(mean_sd))
    colnames(mean_sd_wide) <- paste0(exposure, "_", names(mean_sd))
    
    data.frame(
      outcome               = outcome_label,
      mean_sd_wide,
      est_ci_unadj          = s0$est_ci,
      p_unadj               = s0$p,
      est_ci_sex            = s1$est_ci,
      p_sex                 = s1$p,
      est_ci_sex_tkv        = s2$est_ci,
      p_sex_tkv             = s2$p,
      est_ci_sex_tkv_sbp    = s3$est_ci,
      p_sex_tkv_sbp         = s3$p
    )
  })
  
  do.call(rbind, results)
}

# fit series of adjusted linear models based on covariates suggested by reviewers
table2_resub <- run_adj_lm(outcomes = outcomes_table2,
                                     exposure = "group_risk",
                                     data = panther)
tableS2_resub <- run_adj_lm(outcomes = outcomes_tableS2,
                            exposure = "group_risk",
                            data = panther)
tableS1a_resub <- run_adj_lm(outcomes = outcomes_table2,
                                          exposure = "tanner_stage_comp_panther_cat",
                                          data = panther)
tableS1a_resub<-tableS1a_resub[,-which(colnames(tableS1a_resub)%in%c("est_ci_unadj","est_ci_sex","est_ci_sex_tkv","est_ci_sex_tkv_sbp"))]
tableS1b_resub <- run_adj_lm(outcomes = outcomes_tableS2,
                              exposure = "tanner_stage_comp_panther_cat",
                              data = panther)
tableS1b_resub<-tableS1b_resub[,-which(colnames(tableS1b_resub)%in%c("est_ci_unadj","est_ci_sex","est_ci_sex_tkv","est_ci_sex_tkv_sbp"))]

#write.csv(table2_resub,"table2_resub.csv",row.names=F)
#write.csv(tableS2_resub,"tableS2_resub.csv",row.names=F)
#write.csv(tableS1a_resub,"tableS1a_resub.csv",row.names=F)
#write.csv(tableS1b_resub,"tableS1b_resub.csv",row.names=F)

# sensitivity analysis: no diabetes 
panther_no_diabetes<-subset(panther,panther$group!="Type 2 Diabetes")
table2_resub_no_diabetes <- run_adj_lm(outcomes = outcomes_table2,
                           exposure = "group_risk",
                           data = panther_no_diabetes)
tableS2_resub_no_diabetes <- run_adj_lm(outcomes = outcomes_tableS2,
                            exposure = "group_risk",
                            data = panther_no_diabetes)
write.csv(table2_resub_no_diabetes,"table2_resub_no_diabetes.csv",row.names=F)
write.csv(tableS2_resub_no_diabetes,"tableS2_resub_no_diabetes.csv",row.names=F)
### predictors:
# sex
# total_kidney_volume_ml
# sbp
# metformin_timepoint
# sglti_timepoint
# race_ethnicity_condensed
# age
# hba1c
# group #diabetes status