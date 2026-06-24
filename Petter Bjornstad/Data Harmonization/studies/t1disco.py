"""
This code is designed to pull data from the T1-DISCO REDCap project and output data in a long format with one row per study procedure per visit.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"
# =====================================================================
# WHAT THIS FILE DOES
# Cleans the T1-DISCO REDCap project into a harmonized DataFrame
# (long format, one row per study procedure per visit).
# Called by: data_harmonization.py via clean_t1disco()
#
# INPUTS:  REDCap API (token from api_tokens.csv), data_dictionary_master.csv
# OUTPUT:  returns a pandas DataFrame (not written to disk here)
# DEPENDS: harmonization_functions.combine_checkboxes, pfas_data_merge.PFAS
#          (PFAS is imported but currently unused in this function)
#
# SECTIONS BELOW: Demographics, Medical History (commented out), Vitals,
#                 Labs, Clamp, Renal Clearance, Missingness, Merge datasets
# =====================================================================


def clean_t1disco():
    # Libraries
    import os
    import sys
    sys.path.insert(0, os.path.expanduser('~') +
                    "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    import redcap
    import pandas as pd
    import numpy as np
    from datetime import timedelta
    from natsort import natsorted, ns
    # __pipeline_path_bootstrap__ (files moved into subfolders; resolve repo root from this file)
    import os as _os, sys as _sys
    _ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
    for _sub in ("calculations", "exports", "studies",):
        _p = _os.path.join(_ROOT, _sub)
        if _p not in _sys.path:
            _sys.path.insert(0, _p)
    from harmonization_functions import combine_checkboxes
    from pfas_data_merge import PFAS
    # REDCap project variables
    import getpass
    user = getpass.getuser() 
    
    if user == "choiyej":
        base_data_path = "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/"
        git_path = "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
    elif user == "pylell":
        base_data_path = "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
        git_path = "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    elif user == "kristenmiller":
        base_data_path = "/Users/kristenmiller/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/"
        git_path = "/Users/kristenmiller/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    elif user == "shivaniramesh":
        base_data_path = os.path.expanduser("~/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
        git_path = "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    else:
        sys.exit(f"Unknown user: please specify root path for this user. (Detected user: {user})")

    # Look up this study's API token and open the REDCap project connection
    tokens = pd.read_csv(base_data_path + "/Data Harmonization/api_tokens.csv")
        #"/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "T1-DISCO", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Sentinel values (and their string forms plus blank) treated as missing
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]

    # Load the master data dictionary
    dictionary = pd.read_csv(base_data_path + "Data Harmonization/data_dictionary_master.csv")
    redcap_cols = ["redcap_event_name",
                   "redcap_repeat_instrument", "redcap_repeat_instance"]
    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    var =  [v for v in meta.loc[meta["form_name"] == "demographics", "field_name"]]
    # # Export
    demo = pd.DataFrame(proj.export_records(fields=var))
    demo.replace(rep, np.nan, inplace=True)    
    # T1-DISCO is a single-group T1D study
    demo["group"] = "Type 1 Diabetes"
    # Drop REDCap bookkeeping and identifying/free-text columns
    demo.drop(redcap_cols + ["name_first", "name_last", "phone", "race_other"], axis=1, inplace=True)
    #Race columns combined into one
    demo = combine_checkboxes(demo, base_name="race", levels=[
        "American Indian or Alaskan Native", "Asian",
        "Hawaiian or Pacific Islander", "Black or African American",
        "White", "Unknown", "Other"])
    # Same for ethnicity
    demo = combine_checkboxes(demo,
                            base_name="ethnicity",
                            levels=["Hispanic or Latino",
                                    "Not Hispanic or Latino",
                                    "Unknown/Not Reported"])
    demo["sex"] = demo["sex"].replace({"1": "Male", "2": "Female", "3": "Other"})
    demo = demo.drop(columns = ["age_consent"])
    # --------------------------------------------------------------------------
    # Medical History
    # --------------------------------------------------------------------------
    # var =  [v for v in meta.loc[meta["form_name"] == "medical_history", "field_name"]]
    # med_hist = pd.DataFrame(proj.export_records(fields=var))
    # med_hist.replace(rep, np.nan, inplace=True)
    # med_hist["visit"] = med_hist["redcap_event_name"].replace({"visit_1__screen_arm_1": "Screening", "visit_2__baseline_arm_1": "Baseline", 
    #                                                "visit_3__1_month_arm_1": "1_month", "visit_4__2_month_arm_1": "2_month", 
    #                                                "visit_5__4_month_arm_1": "4_month", "visit_6__6_month_arm_1": "6_month",
    #                                                "visit_7__posttreat_arm_1": "Post-Treatment", "visit_8__9_month_arm_1": "9_month"})
    # med_hist.drop(redcap_cols + ["htn_med___6", "cholesterol_med___3", "cgm_model", "htn_med_other", "abnormal_cholesterol_yn", 
    #                              "hx_cv_positive___4", "hx_cv_other", "hx_cv_meds", "hx_cv_medlist", "hx_met", 
    #                              "hx_met_positive___1", "hx_met_positive___2", "hx_met_positive___3", "hx_met_positive___5", "hx_met_positive___6", 
    #                              "hx_met_positive___7", "hx_met_other", "steroids", "steroid_type___1", "steroid_type___2", "steroid_type___3", 
    #                              "steroid_type___4", "steroid_type___5", "antipsychotics", "which_glp1_ra", "when_was_glp1_ra_stopped", "thyroid_med",
    #                              "antipsychotic_type___1", "antipsychotic_type___2", "thyroid_hx", "recent_dka", "recent_hh",
    #                              "bleeding_dx", "clotting_dx", "allergies", "mri", "tobacco_amount", "vape_years", "vape_amount", "regular_cycles_yn",
    #                              "sexually_active_yn"], axis=1, inplace=True)

    # med_hist = med_hist.rename({"cgm_use_yn": "cgm_yn", "cgm_brand": "cgm_type","recent_severe_hypo": "hypoglycemia", "hx_htn": "hypertension", 
                                
    #                             "htn_med___1": "ace_inhibitor", "htn_med___2": "angiotensin_receptor_blocker", 
    #                             "htn_med___3": "beta_blocker", "htn_med___4": "calcium_channel_blocker", "htn_med___5": "diuretic", 
    #                             "cholesterol_med___1": "statin", "cholesterol_med___2": "fibrates", "hx_met_medlist": "met_meds_yes", 
    #                             "glp1ra_hx": "epic_ever_glp1ra_1", "actively_taking_glp1_ra": "epic_glp1ra_1"}, axis = 1)

    # med_hist["cgm_type"] = med_hist["cgm_type"].replace({"1": "Dexcom", "2": "Medtronic", "3": "Abbott", "4": "Other"})
    # med_hist = med_hist.rename({
    #     "hx_cv_positive___1": "congestive_heart_failure", "hx_cv_positive___2": "coronary_artery_disease", 
    #     "hx_cv_positive___3": "peripheral_vascular_disease"}, axis=1)
    # med_hist = med_hist.rename({
    #     "hx_met_positive___4": "apnea"}, axis = 1)

    # meds_yn = ["cgm_yn", "ace_inhibitor", "angiotensin_receptor_blocker", "beta_blocker", "calcium_channel_blocker", "diuretic", "hypertension", 
    #            "hx_htn_control", "statin", "fibrates", "hypoglycemia", "hx_cv", "congestive_heart_failure", "coronary_artery_disease", 
    #            "peripheral_vascular_disease", "apnea", "hx_met_meds", "epic_ever_glp1ra_1", "epic_glp1ra_1", "hx_anemia", "tobacco", "vape"]
    # med_hist[meds_yn] = med_hist[meds_yn].apply(lambda col: col.map(
    #      lambda x: "Yes" if x == "1" or x is True else ("No" if pd.notna(x) and x != "" else x)
    #      ))
    # ----------------------------------------------------------------------------
    # Vitals
    # ----------------------------------------------------------------------------
    var = ["record_id", "visit_type", "anthro_date", "anthro_height", "anthro_weight", "anthro_bmi"]
    vitals = pd.DataFrame(proj.export_records(fields=var))
    vitals.drop(redcap_cols, axis=1, inplace=True)
    vitals["visit_type"] = vitals["visit_type"].replace({"1": "screening", "2": "baseline", "3": "2_month", "4": "4_month", "5": "6_month", "6": "9_month"})
    vitals = vitals.rename({"visit_type": "visit", "anthro_date": "date", "anthro_height": "height", "anthro_weight": "weight", "anthro_bmi": "bmi"}, axis=1)
    vitals.replace(rep, np.nan, inplace=True)
    vitals["procedure"] = "labs"

    # ----------------------------------------------------------------------------
    #Labs
    # ----------------------------------------------------------------------------
    var = ["record_id", "visit_type", "labs_date", "bun", "screen_creat_s"]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs.drop(redcap_cols, axis=1, inplace=True)
    labs["visit_type"] = labs["visit_type"].replace({"1": "screening", "2": "baseline", "3": "2_month", "4": "4_month", "5": "6_month", "6": "9_month"})
    labs = labs.rename({"visit_type": "visit", "labs_date": "date", "screen_creat_s": "creatinine_s"}, axis=1)
    labs.replace(rep, np.nan, inplace=True)
    labs["procedure"] = "labs"

    # ----------------------------------------------------------------------------
    # Clamp
    # ----------------------------------------------------------------------------
    var = ["record_id", "clamp_visit", "clamp_date", "clamp_wt", "clamp_ht", "clamp_bmi"]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    clamp.drop(redcap_cols, axis=1, inplace=True)
    clamp["clamp_visit"] = clamp["clamp_visit"].replace({"1": "screening", "2": "baseline", "3": "2_month", "4": "4_month", "5": "6_month", "6": "9_month"})
    clamp = clamp.rename({"clamp_visit": "visit", "clamp_date": "date", "clamp_wt": "weight", "clamp_ht": "height", "clamp_bmi": "bmi"}, axis=1)
    clamp.replace(rep, np.nan, inplace=True)
    clamp["procedure"] = "clamp"

    # ----------------------------------------------------------------------------
    # Renal Clearance
    # ----------------------------------------------------------------------------
    var = ["record_id", "renal_visit", "urine_microalbumin_0", "urine_creatinine_0"]
    rc = pd.DataFrame(proj.export_records(fields=var))
    print("DEBUG rc raw shape:", rc.shape, "| columns:", rc.columns.tolist())
    rc.drop(redcap_cols, axis=1, inplace=True)
    rc["renal_visit"] = rc["renal_visit"].replace({"1": "baseline", "2": "post_treatment"})
    rc = rc.rename({"renal_visit": "visit", "urine_microalbumin_0": "microalbumin_u", "urine_creatinine_0": "creatinine_u"}, axis=1)
    rc.replace(rep, np.nan, inplace=True)
    rc["procedure"] = "renal_clearance_testing"
    print("DEBUG rc after processing:", rc.shape, "| microalbumin_u non-null:", rc['microalbumin_u'].notna().sum())

    # ----------------------------------------------------------------------------
    # Missingness
    # ----------------------------------------------------------------------------
    # Drop rows that are mostly empty (keep rows with at least thresh non-null values)
    demo.dropna(thresh=5, axis=0, inplace=True)
    vitals.dropna(thresh=3, axis=0, inplace=True)
    labs.dropna(thresh=3, axis=0, inplace=True)
    clamp.dropna(thresh=4, axis=0, inplace=True)
    rc.dropna(thresh=3, axis=0, inplace=True)

    # ----------------------------------------------------------------------------
    # Merge datasets
    # ----------------------------------------------------------------------------
    # Stack each procedure as its own rows, then attach demographics by record
    df = pd.concat([vitals, labs], join='outer', ignore_index=True)
    df = pd.concat([df, clamp], join='outer', ignore_index=True)
    df = pd.concat([df, rc], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.copy()

    df["study"] = "T1-DISCO"
    # Prefix record IDs with the study tag unless already present
    df["record_id"] = df["record_id"].astype(str)
    mask = ~df["record_id"].str.startswith("DSCO-")
    df.loc[mask, "record_id"] = "DSCO-" + df.loc[mask, "record_id"]
    # Put identifier/visit columns first, then natural-sort the remaining measures
    dem_cols = ["mrn", "dob", "group", "race", "ethnicity"]
    id_cols = ["record_id", "study"] + \
    dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # Impose chronological visit ordering
    df["visit"] = pd.Categorical(df["visit"],
                                 categories=["screening", "baseline", "2_month", "4_month", "6_month", "post_treatment", "9_month"],
                                 ordered=True)
    return df
