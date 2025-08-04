"""
This code is designed to pull data from multiple REDCap projects and harmonize
the data in a single dataset. Some studies are cross-sectional but include
measures at multiple visits, and some studies are longitudinal. So, this code
outputs data in a "semi-long" format with one row per study procedure, and a
visit column for longitudinal clustering. The data cleaning process for each 
individual dataset is a separate function, and this function puts them together 
and performs some formatting tweaks.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def harmonize_data():
    # Libraries
    import os
    import sys
    sys.path.insert(0, os.path.expanduser('~') +
                    "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    import pandas as pd
    import numpy as np
    from natsort import natsorted, ns
    from casper import clean_casper
    from coffee import clean_coffee
    from crocodile import clean_crocodile
    from improve import clean_improve
    from penguin import clean_penguin
    from renal_heir import clean_renal_heir
    from renal_heiritage import clean_renal_heiritage
    from panther import clean_panther
    from panda import clean_panda
    from attempt import clean_attempt
    from harmonization_functions import calc_egfr
    # Use individual data functions to import cleaned DFs
    casper = clean_casper()
    coffee = clean_coffee()
    crocodile = clean_crocodile()
    improve = clean_improve()
    penguin = clean_penguin()
    renal_heir = clean_renal_heir()
    renal_heiritage = clean_renal_heiritage()
    panther = clean_panther()
    panda = clean_panda()
    attempt = clean_attempt()
    # Merge
    harmonized = pd.concat([casper, coffee], join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, crocodile],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, improve],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, penguin],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, renal_heir],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, renal_heiritage],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, panther],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, panda],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, attempt],
                           join='outer', ignore_index=True)
                           
    # Fix levels of categorical variables
    harmonized["visit"] = \
        pd.Categorical(harmonized["visit"],
                       categories=['screening', 'baseline', 'pre_surgery',
                                   '3_months_post_surgery', '4_months_post', '12_months_post_surgery',
                                   'year_1', 'year_2', 'year_3', 'year_4'],
                       ordered=True)
    harmonized["race"].replace(
        ["American Indian or Alaskan Native & White",
         "Black or African American & White",
         'American Indian or Alaskan Native & Black or African American',
         'Asian & White'], "More Than One", inplace=True)
    harmonized["race"].replace(
        {'Black/African American': 'Black or African American',
         "": "Unknown"}, inplace=True)
    harmonized["ethnicity"].replace({"": "Unknown"}, inplace=True)
    race_ethnicity = harmonized["race"] + \
        ", " + harmonized["ethnicity"]
    harmonized = pd.concat([harmonized, race_ethnicity], axis=1)
    harmonized.rename({0: "race_ethnicity"}, axis=1, inplace=True)
    # Replace blanks with missing
    harmonized.replace("", np.nan, inplace=True)
    # Convert to numeric
    harmonized = harmonized.apply(
        pd.to_numeric, errors='ignore')    
    # Date variables
    dates = ["dob", "date", "diabetes_dx_date"]
    harmonized[dates] = \
        harmonized[dates].apply(pd.to_datetime, errors='coerce')
    # Calculated variables
    # Age
    age = round((harmonized["date"] - harmonized["dob"]).dt.days / 365.25, 2)
    harmonized = pd.concat([harmonized, age], axis=1)
    harmonized.rename({0: "age"}, axis=1, inplace=True)
    # BMI
    harmonized["bmi"] = pd.to_numeric(harmonized["weight"])/((pd.to_numeric(harmonized["height"])/100)**2)
    # Diabetes duration
    disease_duration = \
        round((harmonized["date"] -
              harmonized["diabetes_dx_date"]).dt.days / 365.25, 2)
    harmonized = pd.concat([harmonized, disease_duration], axis=1)
    harmonized.rename({0: "diabetes_duration"}, axis=1, inplace=True)
    # eGFR
    harmonized["age"] = pd.to_numeric(harmonized["age"], errors="coerce")
    harmonized["creatinine_s"] = pd.to_numeric(harmonized["creatinine_s"], errors="coerce")
    harmonized["cystatin_c_s"] = pd.to_numeric(harmonized["cystatin_c_s"], errors="coerce")
    harmonized["bun"] = pd.to_numeric(harmonized["bun"], errors="coerce")
    harmonized["height"] = pd.to_numeric(harmonized["height"], errors="coerce")
    harmonized = calc_egfr(harmonized, age="age",
                           serum_creatinine="creatinine_s", cystatin_c="cystatin_c_s",
                           bun="bun", height="height", sex="sex", male="Male", female="Female", alpha=0.5)
    # Kidney volume
    harmonized["total_kidney_volume_ml"] = \
        harmonized["left_kidney_volume_ml"] + \
        harmonized["right_kidney_volume_ml"]
    harmonized["total_kidney_volume_ml_manual"] = harmonized.apply(lambda row: row["volume_left_manual"] + row["volume_right_manual"], axis=1)
    harmonized = harmonized.assign(
        ht_adj_tkv = harmonized["total_kidney_volume_ml"] / (harmonized.groupby("record_id")["height"].transform("mean") / 100))
    harmonized = harmonized.assign(
        ht_adj_tkv_manual = harmonized["total_kidney_volume_ml_manual"] / (harmonized.groupby("record_id")["height"].transform("mean") / 100))
    # PCASL
    harmonized["pcasl3d_left"] = pd.to_numeric(harmonized["pcasl3d_left"], errors='coerce')
    harmonized["pcasl3d_right"] = pd.to_numeric(harmonized["pcasl3d_right"], errors='coerce')
    harmonized["avg_pcascl"]= \
        harmonized[["pcasl3d_left", "pcasl3d_right"]].apply(lambda x: x.mean(), axis=1)
    # Average R2*
    harmonized["avg_k_r2"]= \
        harmonized[["bold_l_bl_kidney", "bold_r_bl_kidney"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_c_r2"]= \
        harmonized[["bold_l_bl_cortex", "bold_r_bl_cortex"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_r2"]= \
        harmonized[["bold_l_bl_medulla", "bold_r_bl_medulla"]].apply(lambda x: x.mean(), axis=1)   
    # Average T1
    harmonized["avg_k_t1"]= \
        harmonized[["bold_l_t1_kidney", "bold_r_t1_kidney"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_c_t1"]= \
        harmonized[["bold_l_t1_cortex", "bold_r_t1_cortex"]].apply(lambda x: x.mean(), axis=1)
    # Average ADC
    harmonized["avg_c_adc"] = \
        harmonized[["adc_left", "adc_right"]].apply(lambda x: x.mean(), axis=1)
    # # Average voxelwise
    harmonized["avg_m_k2_wo_cyst_vw"] = \
        harmonized[["lm_k2_wo_cyst_vw", "rm_k2_wo_cyst_vw"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_c_k2_wo_cyst_vw"] = \
        harmonized[["lc_k2_wo_cyst_vw", "rc_k2_wo_cyst_vw"]].apply(lambda x: x.mean(), axis=1)

    # Average K1
    harmonized["avg_c_k1"]= \
        harmonized[["lc_k1", "rc_k1"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_k1"]= \
        harmonized[["lm_k1", "rm_k1"]].apply(lambda x: x.mean(), axis=1)   
    # Average k2
    harmonized["avg_c_k2"]= \
        harmonized[["lc_k2", "rc_k2"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_k2"]= \
        harmonized[["lm_k2", "rm_k2"]].apply(lambda x: x.mean(), axis=1) 
    # Average F
    harmonized["avg_c_f"]= \
        harmonized[["lc_f", "rc_f"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_k2"]= \
        harmonized[["lm_f", "rm_f"]].apply(lambda x: x.mean(), axis=1)         
        
        
    # Calculate FSOC = bl_bold - pf_bold
    cols = [c for c in harmonized.columns if "_bl_" in c] + \
        [c for c in harmonized.columns if "_pf_" in c]
    harmonized[cols] = harmonized[cols].apply(
        pd.to_numeric, errors='coerce', axis=1)
    harmonized = harmonized.assign(
        fsoc_r_cortex=harmonized["bold_r_bl_cortex"] -
        harmonized["bold_r_pf_cortex"],
        fsoc_r_medulla=harmonized["bold_r_bl_medulla"] -
        harmonized["bold_r_pf_medulla"],
        fsoc_r_kidney=harmonized["bold_r_bl_kidney"] -
        harmonized["bold_r_pf_kidney"],
        fsoc_l_cortex=harmonized["bold_l_bl_cortex"] -
        harmonized["bold_l_pf_cortex"],
        fsoc_l_medulla=harmonized["bold_l_bl_medulla"] -
        harmonized["bold_l_pf_medulla"],
        fsoc_l_kidney=harmonized["bold_l_bl_kidney"] -
        harmonized["bold_l_pf_kidney"])
    # Average FSOC
    harmonized["avg_k_fsoc"]= \
        harmonized[["fsoc_l_kidney", "fsoc_r_kidney"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_c_fsoc"]= \
        harmonized[["fsoc_l_cortex", "fsoc_r_cortex"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_fsoc"]= \
        harmonized[["fsoc_l_medulla", "fsoc_r_medulla"]].apply(lambda x: x.mean(), axis=1)        
    # UACR
    harmonized["acr_u"] = \
        pd.to_numeric(harmonized["microalbumin_u"], errors="coerce") * 100 / \
        pd.to_numeric(harmonized["creatinine_u"], errors="coerce")
    # Albuminuria
    alb = []
    for a in harmonized["acr_u"]:
        if a < 30:
            alb.append("A1")
        elif a >= 30 and a <= 300:
            alb.append("A2")
        elif a > 300:
            alb.append("A3")
        else:
            alb.append(np.nan)
    harmonized["albuminuria_cat"] = alb
    harmonized["elevated_albuminuria"] = pd.cut(
        harmonized["acr_u"], [-float("inf"), 30, float("inf")], right=False, labels=["No", "Yes"])
    # FFA suppression negative to be 0
    harmonized["ffa_suppression"] = np.where(
        harmonized["ffa_suppression"] < 0, 0, harmonized["ffa_suppression"])
    harmonized["p1_ffa_suppression"] = np.where(
        harmonized["p1_ffa_suppression"] < 0, 0, harmonized["p1_ffa_suppression"])
    harmonized["p2_ffa_suppression"] = np.where(
        harmonized["p2_ffa_suppression"] < 0, 0, harmonized["p2_ffa_suppression"])
    # FFA suppression combined
    harmonized = \
        harmonized.assign(ffa_suppression_combined=harmonized["ffa_suppression"].where(
            harmonized["ffa_suppression"].notnull(), harmonized["p2_ffa_suppression"]))
    # Fasting Insulin
    harmonized[["insulin_minus_20", "insulin_minus_10", "insulin_minus_5", "insulin_0"]] = harmonized[[
        "insulin_minus_20", "insulin_minus_10", "insulin_minus_5", "insulin_0"]].apply(pd.to_numeric)
    harmonized["fasting_insulin"] = \
        harmonized[["insulin_minus_20", "insulin_minus_10",
                    "insulin_minus_5", "insulin_0"]].apply(lambda x: x.mean(), axis=1)
    # Fasting FFA
    harmonized["fasting_ffa"] = \
        harmonized[["ffa_minus_20", "ffa_minus_10", "ffa_minus_5", "ffa_0"]].apply(
            lambda x: x.mean(), axis=1)
    # Fasting Glucose
    harmonized["fbg"] = \
        harmonized[["glucose_minus_120", "glucose_minus_90","glucose_minus_20",
        "glucose_minus_10", "glucose_minus_5", "glucose_0", "fbg"]].apply(
            lambda x: x.mean(), axis=1)
    # HOMA-IR (https://link.springer.com/article/10.1007/BF00280883), FBG entered as mg/dL, converting to mmol/L (18 mg/dL = 1 mmol/L)
    harmonized = harmonized.assign(
        homa_ir=(harmonized["fasting_insulin"] * (harmonized["fbg"]/18))/22.5)
    # Adipose IR (fasting_insulin * fasting_ffa)
    harmonized = harmonized.assign(
        adipose_ir=harmonized["fasting_ffa"] * harmonized["fasting_insulin"])
    # SEARCH IS score eIS (estimated insulin sensitivity) (https://academic.oup.com/jcem/article/101/2/686/2811091#81467317)
    harmonized = harmonized.assign(
        search_eis = np.exp((4.64725 - (0.02032 * harmonized.groupby("record_id")["waistcm"].transform("mean")) - 
        (0.09779 * harmonized.groupby("record_id")["hba1c"].transform("mean")) - 
        (0.00235 * harmonized.groupby("record_id")["triglycerides"].transform("mean")))))
    # Co-enroll IDs
    casper_mrns = harmonized.loc[harmonized['study'] == 'CASPER', ['mrn', 'record_id']]
    casper_id_map = dict(zip(casper_mrns['mrn'], casper_mrns['record_id']))
    harmonized['casper_id'] = harmonized.apply(lambda row: casper_id_map[row['mrn']] if row['mrn'] in casper_id_map else '', axis=1)
    coffee_mrns = harmonized.loc[harmonized['study'] == 'COFFEE', ['mrn', 'record_id']]
    coffee_id_map = dict(zip(coffee_mrns['mrn'], coffee_mrns['record_id']))
    harmonized['coffee_id'] = harmonized.apply(lambda row: coffee_id_map[row['mrn']] if row['mrn'] in coffee_id_map else '', axis=1)
    croc_mrns = harmonized.loc[harmonized['study'] == 'CROCODILE', ['mrn', 'record_id']]
    croc_id_map = dict(zip(croc_mrns['mrn'], croc_mrns['record_id']))
    harmonized['croc_id'] = harmonized.apply(lambda row: croc_id_map[row['mrn']] if row['mrn'] in croc_id_map else '', axis=1)
    improve_mrns = harmonized.loc[harmonized['study'] == 'IMPROVE', ['mrn', 'record_id']]
    improve_id_map = dict(zip(improve_mrns['mrn'], improve_mrns['record_id']))
    harmonized['improve_id'] = harmonized.apply(lambda row: improve_id_map[row['mrn']] if row['mrn'] in improve_id_map else '', axis=1)
    penguin_mrns = harmonized.loc[harmonized['study'] == 'PENGUIN', ['mrn', 'record_id']]
    penguin_id_map = dict(zip(penguin_mrns['mrn'], penguin_mrns['record_id']))
    harmonized['penguin_id'] = harmonized.apply(lambda row: penguin_id_map[row['mrn']] if row['mrn'] in penguin_id_map else '', axis=1)
    rh_mrns = harmonized.loc[harmonized['study'] == 'RENAL-HEIR', ['mrn', 'record_id']]
    rh_id_map = dict(zip(rh_mrns['mrn'], rh_mrns['record_id']))
    harmonized['rh_id'] = harmonized.apply(lambda row: rh_id_map[row['mrn']] if row['mrn'] in rh_id_map else '', axis=1)
    rh2_mrns = harmonized.loc[harmonized['study'] == 'RENAL-HEIRitage', ['mrn', 'record_id']]
    rh2_id_map = dict(zip(rh2_mrns['mrn'], rh2_mrns['record_id']))
    harmonized['rh2_id'] = harmonized.apply(lambda row: rh2_id_map[row['mrn']] if row['mrn'] in rh2_id_map else '', axis=1)
    panther_mrns = harmonized.loc[harmonized['study'] == 'PANTHER', ['mrn', 'record_id']]
    panther_id_map = dict(zip(panther_mrns['mrn'], panther_mrns['record_id']))
    harmonized['panther_id'] = harmonized.apply(lambda row: panther_id_map[row['mrn']] if row['mrn'] in panther_id_map else '', axis=1)
    panda_mrns = harmonized.loc[harmonized['study'] == 'PANDA', ['mrn', 'record_id']]
    panda_id_map = dict(zip(panda_mrns['mrn'], panda_mrns['record_id']))
    harmonized['panda_id'] = harmonized.apply(lambda row: panda_id_map[row['mrn']] if row['mrn'] in panda_id_map else '', axis=1)
    attempt_mrns = harmonized.loc[harmonized['study'] == 'ATTEMPT', ['mrn', 'record_id']]
    attempt_id_map = dict(zip(attempt_mrns['mrn'], attempt_mrns['record_id']))
    harmonized['attempt_id'] = harmonized.apply(lambda row: attempt_id_map[row['mrn']] if row['mrn'] in attempt_id_map else '', axis=1)
    
    # Sort columns
    id_cols = ["record_id", "casper_id", "coffee_id", "croc_id", "improve_id", 
                "penguin_id", "rh_id", "rh2_id", "panther_id", "panda_id", "attempt_id",
                "mrn", "co_enroll_id", "study", "dob", "diabetes_dx_date",
               "sex", "race", "ethnicity", "visit", "procedure", "date", "group"]
    other_cols = harmonized.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    harmonized = harmonized[id_cols + other_cols]
    # Sort rows
    harmonized.sort_values(
        ["study", "record_id", "visit", "procedure", "date"], inplace=True)
    # Format dates nicely
    harmonized[dates] = harmonized[dates].apply(
        lambda x: x.dt.strftime('%Y-%m-%d'))
    # Return
    harmonized = harmonized.astype(object)
    return harmonized
