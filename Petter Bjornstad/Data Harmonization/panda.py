"""
This code is designed to pull data from the PANDA REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_panda():
    # Libraries
    import os
    import sys
    import re
    sys.path.insert(0, os.path.expanduser('~') +
                    "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    import redcap
    import pandas as pd
    import numpy as np
    from natsort import natsorted, ns
    from harmonization_functions import combine_checkboxes
    # REDCap project variables
    import getpass
    user = getpass.getuser()  # safer than os.getlogin(), works in more environments

    if user == "choiyej":
        base_data_path = "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/"
        git_path = "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
    elif user == "pylell":
        base_data_path = "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
        git_path = "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    elif user == "shivaniramesh":
        base_data_path = os.path.expanduser("~/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
        git_path = "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    else:
        sys.exit(f"Unknown user: please specify root path for this user. (Detected user: {user})")

    tokens = pd.read_csv(base_data_path + "/Data Harmonization/api_tokens.csv")        #"/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "PANDA", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]

    dictionary = pd.read_csv(base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv")

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["record_id", "crc_substudy", "dob", "diabetes_dx_date", "sex", "race", "ethnicity", "participation_status","mrn"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols, events=["screening_arm_1", "screening__annual_arm_2"]))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo["co_enroll_id"] = demo["crc_substudy"]
    demo["co_enroll_id"] = demo["co_enroll_id"].apply(lambda x: "CRC-" + x if isinstance(x, str) and x.strip() else x)
    dem_cols.remove("crc_substudy")
    demo.drop(["crc_substudy"], axis=1, inplace=True)
    # Race columns combined into one
    demo = combine_checkboxes(demo, base_name="race", levels=[
        'American Indian/Alaska Native', 'Asian', 'Hawaiian/Pacific Islander', 'Black/African American', 'White', 'Other', 'Unknown'])
    # Same for ethnicity
    demo = combine_checkboxes(demo,
                              base_name="ethnicity",
                              levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
    # Relevel sex and group
    demo["sex"].replace({1: "Male", 2: "Female", 3: "Other",
                        "1": "Male", "2": "Female", "3": "Other"}, inplace=True)
    demo["group"]="Type 1 Diabetes"
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)
    demo.drop(["redcap_event_name"], axis = 1, inplace=True)

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    med = med[med['redcap_event_name'].str.contains('screen', na=False)]
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    med_list = {'diabetes_tx___1': "insulin_pump_timepoint",
                'diabetes_tx___2': "insulin_injections_timepoint",
                "diabetes_meds_other___1": "metformin_timepoint",
                "diabetes_meds_other___2": "tzd_timepoint",
                "diabetes_meds_other___3": "glp1_agonist_timepoint",
                "diabetes_meds_other___4": "sglti_timepoint",
                "diabetes_meds_other___5": "other_diabetes_med_timepoint",
                "htn_med___1": "ace_inhibitor",
                "htn_med___2": "angiotensin_receptor_blocker",
                "htn_med___3": "beta_blocker",
                "htn_med___4": "ca_channel_blocker",
                "htn_med___5": "diuretic",
                "htn_med___6": "statin",
                "pump_basal_rate": "pump_basal_rate", 
                "cgm_yn": "cgm_yn"}
    og_names = list(med_list.keys())
    hx_list = [col for col in med.columns if col.startswith('hx_')]
    med = med[["record_id"] + og_names + hx_list]
    med.rename(med_list, axis=1, inplace=True)
    # RAASi
    med = med.assign(raasi_timepoint=np.maximum(pd.to_numeric(
        med["ace_inhibitor"]), pd.to_numeric(med["angiotensin_receptor_blocker"])))
    # Metformin
    med.rename({"diabetes_med_other___1": "metformin_timepoint"},
               axis=1, inplace=True)
    # Insulin
    med = med.assign(insulin_med_timepoint=np.maximum(pd.to_numeric(
        med["insulin_pump_timepoint"]), pd.to_numeric(med["insulin_injections_timepoint"])))
    # Replace 0/1 values with yes/no
    med.iloc[:, 1:] = med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"})
    med["procedure"] = "medications"
    med["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # EPIC Medications
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "epic_meds", "field_name"]]
    epic_med = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    epic_med.replace(rep, np.nan, inplace=True)
    # Replace 0/1 values with yes/no
    epic_med.iloc[:, 1:] = epic_med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 2: "No", "2": "No", 1: "Yes", "1": "Yes"})
    epic_med["procedure"] = "epic_medications"
    epic_med["visit"] = "baseline"
    
    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "vitals_anthropometric_measurements", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.columns = phys.columns.str.replace(r"phys_", "", regex=True)
    phys.rename({"sysbp": "sbp", "diasbp": "dbp"}, inplace=True, axis=1)
    phys["visit"] = phys["redcap_event_name"].apply(lambda x: re.search(r"annual_visit_(\d+)", x))
    phys["visit"] = phys["visit"].apply(lambda x: f"year_{x.group(1)}" if x else "baseline")
    phys.drop(["redcap_event_name"], axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    screen = pd.DataFrame(proj.export_records(forms=["screening_labs", "annual_labs"]))
    screen = screen[screen['redcap_event_name'].str.contains('screen', na=False)]
    # Replace missing values
    screen.replace(rep, np.nan, inplace=True)
    screen.rename({"screen_a1c": "hba1c"}, axis=1, inplace=True)
    screen_num = list(screen.columns.values.tolist())
    screen_num = [s for s in screen_num if any(s.startswith(x) for x in ['screen', 'ann', 'yr'])]
    screen[screen_num] = screen[screen_num].apply(
        pd.to_numeric, errors='coerce')
    screen["creatinine_s"] = screen[['screen_creat_s', 'yr_creat_s']].mean(axis=1)
    screen["creatinine_u"] = screen[['screen_creat_u', 'yr_creat_u']].mean(axis=1)
    screen["microalbumin_u"] = screen[['screen_microalbumin_u', 'yr_microalbumin_u']].mean(axis=1)
    screen["acr_u"] = screen[['screen_uacr', 'ann_uacr']].mean(axis=1)
    screen["date"] = screen.labs_date.fillna(screen.visit_date_yr)
    screen.drop(["redcap_event_name", "screening_labs_complete", "annual_labs_complete",      
                "prescreen_a1c", "prescreen_a1c_date",
                "yr_upt", "screen_upt", "annual_pi_s", "annual_pi_p", "yr_morning_u",
                "yr_morning_urine1", "yr_morning_urine2", "annual_pi_u", 
                'screen_creat_s', 'yr_creat_s', 'screen_creat_u', 'yr_creat_u',
                'screen_microalbumin_u', 'yr_microalbumin_u',
                'screen_uacr', 'ann_uacr', 'labs_date', 'visit_date_yr'], axis=1, inplace=True)    
    screen["visit"] = "baseline"
    screen["procedure"] = "screening"
    # Assume medication review done at screening
    med = pd.merge(med, screen[['record_id', 'date']], on='record_id', how='left')

    # --------------------------------------------------------------------------
    # Labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "annual_labs", "field_name"]]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs = labs.loc[~labs.redcap_event_name.str.contains("screen|annual_visit_1_arm_1")]
    # Replace missing values
    labs.replace(rep, np.nan, inplace=True)
    labs.columns = labs.columns.str.replace(
        r"visit_|bl_|_yr|yr_", "", regex=True)
    labs.rename({"uacr": "acr_u", "a1c": "hba1c",
                "na_u": "sodium_u", "na_s": "sodium_s"}, axis=1, inplace=True)
    labs["procedure"] = "clamp"
    labs["visit"] = labs["redcap_event_name"].apply(lambda x: re.search(r"annual_visit_(\d+)", x))
    labs["visit"] = labs["visit"].apply(lambda x: f"year_{x.group(1)}" if x else "baseline")
    labs.drop(["redcap_event_name"], axis=1, inplace=True)
    
    # --------------------------------------------------------------------------
    # CROC Labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"] 
                                               == "crc_substudy_baseline_vitalslabs", "field_name"]]
    croc_labs = pd.DataFrame(proj.export_records(fields=var))
    croc_labs = croc_labs.loc[~croc_labs.redcap_event_name.str.contains("screen|annual_visit_1_arm_1")]
    # Replace missing values
    croc_labs.drop(["study_id_0", "baseline_vitals_clamp"], axis=1, inplace=True)    
    croc_labs.replace(rep, np.nan, inplace=True)
    croc_labs.columns = croc_labs.columns.str.replace(
        r"visit_|bl_|_yr|yr_|_clamp", "", regex=True)
    croc_labs.rename({"uacr": "acr_u", "a1c": "hba1c",
                "na_u": "sodium_u", "na_s": "sodium_s",
                "s_creat": "creatinine_s", "u_creat": "creatinine_u",
                "u_mab": "urine_mab_baseline", "u_glucose": "urine_glucose",
                "u_na": "na_u", "triglyceride_clamp": "triglycerides", "protein": "tot_protein"}, axis=1, inplace=True)
    croc_labs["procedure"] = "croc_clamp"
    croc_labs["visit"] = croc_labs["redcap_event_name"].apply(lambda x: re.search(r"annual_visit_(\d+)", x))
    croc_labs["visit"] = croc_labs["visit"].apply(lambda x: f"year_{x.group(1)}" if x else "baseline")
    croc_labs.drop(["redcap_event_name"], axis=1, inplace=True)


    # --------------------------------------------------------------------------
    # BOLD/ASL MRI
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "main_study_visit_mri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    mri.replace(rep, np.nan, inplace=True)
    mri.columns = mri.columns.str.replace(
        r"mri_", "", regex=True)
    mri.columns = mri.columns.str.replace(
        r"_pre_r", "_r_bl", regex=True)
    mri.columns = mri.columns.str.replace(
        r"_pre_l", "_l_bl", regex=True)
    mri.columns = mri.columns.str.replace(
        r"_post_r", "_r_pf", regex=True)
    mri.columns = mri.columns.str.replace(
        r"_post_l", "_l_pf", regex=True)        
    mri.rename({"volume_right": "right_kidney_volume_ml",
                "volume_left": "left_kidney_volume_ml"},
               axis=1, inplace=True)
    mri["procedure"] = "mri"
    mri["visit"] = mri["redcap_event_name"].apply(lambda x: re.search(r"annual_visit_(\d+)", x))
    mri["visit"] = mri["visit"].apply(lambda x: f"year_{x.group(1)}" if x else "baseline")
    mri.drop(["redcap_event_name"], axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # DXA Scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "main_study_visit_dxa_scan", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    dxa.replace(rep, np.nan, inplace=True)
    dxa.columns = dxa.columns.str.replace(
        r"dxa_|_percent", "", regex=True)
    dxa.rename({"dxa_date": "date", "bodyfat": "body_fat", 
                "leanmass": "lean_mass",
                "trunkmass": "trunk_mass", "fatmass_kg": "fat_kg",
                "leanmass_kg": "lean_kg", "trunkmass_kg": "trunk_kg",
                "bmd": "bone_mineral_density"}, axis=1, inplace=True)
    dxa_cols = dxa.columns[3:].to_list()
    dxa.rename(dict(zip(dxa_cols, ["dexa_" + d for d in dxa_cols])),
               axis=1, inplace=True)
    dxa["procedure"] = "dxa"
    dxa["visit"] = dxa["redcap_event_name"].apply(lambda x: re.search(r"annual_visit_(\d+)", x))
    dxa["visit"] = dxa["visit"].apply(lambda x: f"year_{x.group(1)}" if x else "baseline")
    dxa.drop(["redcap_event_name"], axis=1, inplace=True)
    dictionary.loc[dictionary['variable_name'].isin(dxa.columns), 'form_name'] = 'dxa'


    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "renal_biopsy", "field_name"]]
    var = var + ["gloms", "gloms_gs", "ifta", "vessels_other", "fia",
                 "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins",
                 "glom_volume_con", "mes_matrix_area",
                 "mes_index", "mes_volume_weibel", "mes_volume_wiggins",
                 "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima",
                 "art_media", "pod_nuc_density", "pod_cell_volume"]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    biopsy.replace(rep, np.nan, inplace=True)
    biopsy.drop([col for col in biopsy.columns if '_yn' in col] +
                [col for col in biopsy.columns if 'procedure_' in col] +
                [col for col in biopsy.columns if 'phys_' in col],
                axis=1, inplace=True)
    biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"labs_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"vitals_", "", regex=True)
    biopsy.rename({"hg": "hemoglobin"}, inplace=True, axis=1)
    biopsy["procedure"] = "kidney_biopsy"
    biopsy["visit"] = biopsy["redcap_event_name"].apply(lambda x: re.search(r"annual_visit_(\d+)", x))
    biopsy["visit"] = biopsy["visit"].apply(lambda x: f"year_{x.group(1)}" if x else "baseline")
    biopsy.drop(["redcap_event_name"], axis=1, inplace=True)

    
    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "crc_substudy_clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    clamp.replace(rep, np.nan, inplace=True)
    # Format
    clamp.drop(["clamp_yn", "clamp_ffa",
                "clamp_insulin", "hct_yn", "clamp_bg"], axis=1, inplace=True)
    clamp.rename({"clamp_wt": "weight", "clamp_ht": "height",
                  "cystatin_c": "cystatin_c_s", "hct_210": "hematocrit_210",
                  "acr_baseline": "acr_u", "acr_250": "acr_u_pm",
                  "clamp_d20": "d20_infusion"},
                 inplace=True, axis=1)
    clamp.columns = clamp.columns.str.replace(r"clamp_", "", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"insulin_minus", "insulin_minus_", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"ffa_minus", "ffa_minus_", regex=True)
    clamp.columns = clamp.columns.str.replace(r"bg_", "glucose_", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"glucose_minus", "glucose_minus_", regex=True)
    clamp["procedure"] = "clamp"
    clamp["visit"] = clamp["redcap_event_name"].apply(lambda x: re.search(r"annual_visit_(\d+)", x))
    clamp["visit"] = clamp["visit"].apply(lambda x: f"year_{x.group(1)}" if x else "baseline")
    clamp.drop(["redcap_event_name"], axis=1, inplace=True)
    clamp["insulin_sensitivity_method"] = "hyperinsulinemic_euglycemic_clamp"

    num_vars = ["d20_infusion", "weight"]
    clamp[num_vars] = clamp[num_vars].apply(
        pd.to_numeric, errors='coerce')
    
    clamp["gir_190"] = (clamp["d20_infusion"] * 190 / 60) / clamp["weight"] # previously M-value
    clamp["gir_200"] = (clamp["d20_infusion"] * 200 / 60) / clamp["weight"]
    
     # FFA
    ffa = [c for c in clamp.columns if "ffa_" in c]
    clamp[ffa] = clamp[ffa].apply(
        pd.to_numeric, errors='coerce')
    clamp["baseline_ffa"] = \
        clamp[['ffa_minus_20', 'ffa_minus_10', 'ffa_0']].mean(axis=1)
    clamp["p1_steady_state_ffa"] = \
        clamp[['ffa_70', 'ffa_80', 'ffa_90']].mean(axis=1)
    clamp["p2_steady_state_ffa"] = \
        clamp[['ffa_250', 'ffa_260', 'ffa_270']].mean(axis=1)
    clamp["p1_ffa_suppression"] = (
        (clamp["baseline_ffa"] - clamp["p1_steady_state_ffa"]) / clamp["baseline_ffa"]) * 100
    clamp["p2_ffa_suppression"] = (
        (clamp["baseline_ffa"] - clamp["p2_steady_state_ffa"]) / clamp["baseline_ffa"]) * 100
        
    # Insulin
    insulin = [c for c in clamp.columns if "insulin_" in c]
    clamp[insulin] = clamp[insulin].apply(
        pd.to_numeric, errors='coerce')
    clamp["baseline_insulin"] = \
        clamp[['insulin_minus_20', 'insulin_minus_10', 'insulin_0']].mean(
            axis=1)
    clamp["p1_steady_state_insulin"] = \
        clamp[['insulin_70', 'insulin_80', 'insulin_90']].mean(axis=1)
    clamp["p2_steady_state_insulin"] = \
        clamp[['insulin_250', 'insulin_260', 'insulin_270']].mean(axis=1)
    clamp["ffa_method"] = "hyperinsulinemic_euglycemic_clamp"
    clamp["ffa_method"] = "hyperinsulinemic_euglycemic_clamp"

    # --------------------------------------------------------------------------
    # Renal Clearance Testing
    # --------------------------------------------------------------------------

    var = ["record_id"] + ["hct_210"] + ["visit_map"] + ["phys_map"] + [v for v in meta.loc[meta["form_name"]
                                               == "crc_substudy_renal_clearance_testing", "field_name"]]
    rct = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    rct.replace(rep, np.nan, inplace=True)
    rename = {"gfr_raw": "gfr_raw_plasma_urine", "gfr_bsa": "gfr_bsa_plasma_urine",
              "erpf_raw": "erpf_raw_plasma_urine", "erpf": "erpf_bsa_plasma_urine",
              "gfr_15mgmin": "gfr_raw_plasma", "gfrbsa": "gfr_bsa_plasma",
              "erpf_pah_85": "erpf_raw_plasma", "erpfbsa": "erpf_bsa_plasma"}
    rct.rename(rename, axis=1, inplace=True)
    
    # Calculate variables
    rct_vars = ["gfr_raw_plasma", "erpf_raw_plasma", "visit_map", "phys_map", "hct_210"]
    rct[rct_vars] = rct[rct_vars].apply(pd.to_numeric, errors='coerce')
    rct["map"] = rct[["visit_map", "phys_map"]].mean(axis=1)
    rct["erpf_raw_plasma_seconds"] = rct["erpf_raw_plasma"]/60
    rct["gfr_raw_plasma_seconds"] = rct["gfr_raw_plasma"]/60
    # Filtration Fraction
    rct["ff"] = rct["gfr_raw_plasma"]/rct["erpf_raw_plasma"] 
    # Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
    rct["kfg"] = 0.1012
    # Filtration pressure across glomerular capillaries
    rct["deltapf"] = (rct["gfr_raw_plasma"]/60)/rct["kfg"] 
    
    # No tot_protein variable collected
    # Plasma protein mean concentration
    # rct["cm"] = (rct["bl_tot_protein"]/rct["ff"])*np.log(1/(1-rct["ff"])) 
    # # Pi G (Oncotic pressure)
    # rct["pg"] = 5*(rct["cm"]-2)
    # # Glomerular Pressure
    # rct["glomerular_pressure"] = rct["pg"] + rct["deltapf"] + 10
    # Afferent Arteriolar Resistance
    # rct["ra"] = ((rct["map"] - rct["glomerular_pressure"]) / rct["rbf_seconds"]) * 1328    
    # rct.loc[~(rct['ra'] > 0), 'ra']=np.nan    
    
    # Renal Blood Flow
    rct["rbf"] = (rct["erpf_raw_plasma"]) / (1 - rct["hct_210"]/100)
    rct["rbf_seconds"] = (rct["erpf_raw_plasma_seconds"]) / (1 - rct["hct_210"]/100)
    # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    rct["rvr"] = rct["map"] / rct["rbf"]
    # Efferent Arteriorlar Resistance 
    rct["re"] = (rct["gfr_raw_plasma_seconds"]) / (rct["kfg"] * (rct["rbf_seconds"] - (rct["gfr_raw_plasma_seconds"]))) * 1328

    # Reduce rct dataset
    rct = rct[["record_id", "ff", "kfg", "deltapf", "rbf", "rvr", "re", 
               "pah_raw", "pah_sd", "pah_cv", "pahcl_12_8mgmin"] + list(rename.values())] 
    rct["procedure"] = "clamp"
    rct["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Astrazeneca urine metabolomics
    # --------------------------------------------------------------------------
    
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "astrazeneca_urine_metabolomics", "field_name"]]
    az_u_metab = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    az_u_metab.replace(rep, np.nan, inplace=True)
    az_u_metab["procedure"] = "az_u_metab"
    az_u_metab["visit"] = "baseline"
    az_u_metab["date"] = labs["date"]
    
    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=5, axis=0, inplace=True)
    epic_med.dropna(thresh=6, axis=0, inplace=True)
    phys.dropna(thresh=5, axis=0, inplace=True)
    screen.dropna(thresh=4, axis=0, inplace=True)
    labs.dropna(thresh=5, axis=0, inplace=True)
    mri.dropna(thresh=4, axis=0, inplace=True)
    dxa.dropna(thresh=5, axis=0, inplace=True)
    clamp.dropna(thresh=7, axis=0, inplace=True)
    rct.dropna(thresh=5, axis=0, inplace=True)
    biopsy.dropna(thresh=20, axis=0, inplace=True)
    az_u_metab.dropna(thresh=6, axis=0, inplace=True)
    croc_labs.dropna(thresh=5, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------
    # Procedure = clamp
    clamp_merge = pd.merge(clamp, labs, how="outer")
    clamp_merge = pd.merge(clamp_merge, rct,  how="outer")
    # Everything else
    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp_merge], join='outer', ignore_index=True)
    df = pd.concat([df, croc_labs], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, epic_med], join='outer', ignore_index=True)
    df = pd.concat([df, az_u_metab], join='outer', ignore_index=True)
    df = pd.merge(df, demo, on='record_id', how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df["study"] = "PANDA"
    id_cols = ["record_id", "co_enroll_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # Sort
    df.sort_values(["record_id", "date", "procedure"], inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Print final data
    return df
