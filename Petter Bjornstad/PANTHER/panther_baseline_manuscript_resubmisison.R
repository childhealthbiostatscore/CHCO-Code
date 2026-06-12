### PANTHER multivariable models ###


### outcome variables:

#Table 2:
lm(avg_c_r2~sex+group_risk,data=panther)
lm(avg_c_r2~sex+tanner_stage_comp_panther_cat,data=panther)

kable(summary(tableby(group_risk ~  total_kidney_volume_ml + ht_adj_tkv + avg_c_r2 + avg_k_r2 + avg_c_t1 + avg_pcascl + avg_c_adc + 
                        mm_si + mm_ir + mm_bcell + mm_di + gfr_raw_plasma + gfr_bsa_plasma + erpf_raw_plasma + erpf_bsa_plasma +
                        ra + re + anova(rvr, digits = 3) + glomerular_pressure, 
                      data = panther, 
                      control = table_controls,
                      test = T), digits = 1))
#Table s1: 
"log_lh", "log_estrad", "log_tot_test", "log_free_test", "bl_dheas", "igf_1", "igf1_z_score_calc"

### predictors:
group_risk
tanner_stage_comp_panther_cat
sex
total_kidney_volume_ml
sbp
metformin_timepoint
sglti_timepoint
race_ethnicity_condensed
age
hba1c
group #diabetes status