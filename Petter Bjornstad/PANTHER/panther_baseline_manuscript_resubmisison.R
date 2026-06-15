### PANTHER multivariable models ###
source("/Users/kristenmiller/Documents/GitHub/CHCO-Code/Petter Bjornstad/PANTHER/panther_baseline_manuscript.qmd")

output_dir<-"/Users/kristenmiller/Library/CloudStorage/OneDrive-SharedLibraries-UW/bjornstad_pyle_tommerdahl_lab - Documents/Dissemination/Manuscripts/P031_PANTHER_baseline_yc/Drafts"
setwd(output_dir)
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
outcomes_table2<-c("total_kidney_volume_ml","ht_adj_tkv","avg_c_r2","avg_k_r2","avg_c_t1","avg_pcascl","avg_c_adc",
                   
                   "gfr_raw_plasma","gfr_bsa_plasma","erpf_raw_plasma","erpf_bsa_plasma",
              "ra","re","rvr","glomerular_pressure")

run_sex_adj_lm <- function(outcomes, exposure, data) {
  
  results <- lapply(outcomes, function(outcome) {
    
    # outcome<-"total_kidney_volume_ml"
    # exposure<-"group_risk"
    # data<-panther
    formula <- as.formula(paste(outcome, "~ sex +", exposure))
    lm_adj <- lm(formula, data = data)

    # Overall F test p-value for the exposure
    anova_p <- round(anova(lm_adj)[exposure, "Pr(>F)"], 3)
    anova_p <- ifelse(anova_p == 0, "<0.001", as.character(anova_p))
    
    outcome_label <- label(data[[outcome]])
    if (is.null(outcome_label) || outcome_label == "") outcome_label <- outcome

    # Estimated marginal means by exposure level
    emm <- as.data.frame(emmeans(lm_adj, specs = exposure))
    emm$emmean <- round(emm$emmean, 3)
    emm$SE     <- round(emm$SE, 3)
    emm$emm_se <- paste0(emm$emmean, " (", emm$SE,")")
    emm_wide <- as.data.frame(t(emm$emm_se))
    colnames(emm_wide) <- paste0(exposure, "_", emm[[exposure]])
    
    data.frame(
      outcome  = outcome_label,
      #exposure = exposure,
      emm_wide,
      pvalue   = anova_p
    )
  })
  
  do.call(rbind, results)
}
run_sex_tkv_adj_lm <- function(outcomes, exposure, data) {
  
  results <- lapply(outcomes, function(outcome) {
    
    # outcome<-"gfr_raw_plasma"
    # exposure<-"group_risk"
    # data<-panther
    formula <- as.formula(paste(outcome, "~ sex + total_kidney_volume_ml + ", exposure))
    lm_adj <- lm(formula, data = data)
    
    # Overall F test p-value for the exposure
    anova_p <- round(anova(lm_adj)[exposure, "Pr(>F)"], 3)
    anova_p <- ifelse(anova_p == 0, "<0.001", as.character(anova_p))
    
    outcome_label <- label(data[[outcome]])
    if (is.null(outcome_label) || outcome_label == "") outcome_label <- outcome
    
    # Estimated marginal means by exposure level
    emm <- as.data.frame(emmeans(lm_adj, specs = exposure))
    emm$emmean <- round(emm$emmean, 3)
    emm$SE     <- round(emm$SE, 3)
    emm$emm_se <- paste0(emm$emmean, " (", emm$SE,")")
    emm_wide <- as.data.frame(t(emm$emm_se))
    colnames(emm_wide) <- paste0(exposure, "_", emm[[exposure]])
    
    data.frame(
      outcome  = paste0(outcome_label,"*"),
      #exposure = exposure,
      emm_wide,
      pvalue   = anova_p
    )
  })
  
  do.call(rbind, results)
}

#Table 2: 
adj_pvalues_sex_table2 <- run_sex_adj_lm(outcomes = outcomes_table2,
                                     exposure = "group_risk",
                                     data = panther)
adj_pvalues_sex_tkv_table2 <- run_sex_tkv_adj_lm(outcomes = c("gfr_raw_plasma","erpf_raw_plasma"),
                                     exposure = "group_risk",
                                     data = panther)
table2_resub<-rbind(adj_pvalues_sex_table2,adj_pvalues_sex_tkv_table2)
write.csv(table2_resub,"table2_resub.csv",row.names=F)

#Table S2:
outcomes_tableS2<- c("lh","estrad","tot_test","free_test","bl_dheas","igf_1","igf1_z_score_calc","cpeptide",
                     "dexa_lean_mass","dexa_body_fat","dexa_fat_kg","dexa_trunk_mass",
                     "mm_si", "mm_ir","mm_bcell","mm_di")
tableS2_resub <- run_adj_lm(outcomes = outcomes_tableS2,
                                 exposure = "group_risk",
                                 data = panther)
write.csv(tableS2_resub,"tableS2_resub.csv",row.names=F)

#Table S1: split into matching tables of risk group:
adj_pvalues_sex_tableS1 <- run_sex_adj_lm(outcomes = outcomes_table2,
                                         exposure = "tanner_stage_comp_panther_cat",
                                         data = panther)
adj_pvalues_sex_tkv_tableS1 <- run_sex_tkv_adj_lm(outcomes = c("gfr_raw_plasma","erpf_raw_plasma"),
                                                 exposure = "tanner_stage_comp_panther_cat",
                                                 data = panther)
tableS1_a_resub<-rbind(adj_pvalues_sex_tableS1,adj_pvalues_sex_tkv_tableS1)
write.csv(tableS1_a_resub,"tableS1_a_resub.csv",row.names=F)

#Table S2:
tableS1_b_resub <- run_adj_lm(outcomes = outcomes_tableS2,
                            exposure = "tanner_stage_comp_panther_cat",
                            data = panther)
write.csv(tableS1_b_resub,"tableS2_b_resub.csv",row.names=F)


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