"""
This code is designed to pull data from the RPC2 REDCap project and output data in a long format with one row per study procedure per visit.
"""
__author__ = ["Ye Ji Choi"]
__credits__ = ["Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Ye Ji Choi"
__email__ = "yejichoi@uw.edu"
__status__ = "Dev"


# Function to clean and structure REDCap data
def clean_rpc2_redcap():
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
    from harmonization_functions import combine_checkboxes
    import re
    
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
            #"/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "RPC2", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Columns to drop
    redcap_cols = ["redcap_event_name",
                   "redcap_repeat_instrument", "redcap_repeat_instance"]
    # Get metadata
    meta = pd.DataFrame(proj.metadata)
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]
    invalid = ["", " ", np.nan]
    
    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["subject_id", "date_of_consent", "mr_number", "dob", "gender", "race", "ethnicity", "participation_status", "diabetes_hx_type"]
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    print(demo[redcap_cols])
    demo.replace(rep, np.nan, inplace=True)
    demo.rename({"mr_number": "mrn", "gender": "sex",  "date_of_consent":"consent_date"}, axis=1, inplace=True)
    demo = combine_checkboxes(demo, base_name="race", levels=["American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
    demo = combine_checkboxes(demo, base_name="ethnicity", levels=["Hispanic", "Non-Hispanic", "Unknown/Not Reported"])
    demo["sex"].replace({1: "Male", 0: "Female", 2: "Other",
                        "1": "Male", "0": "Female", "2": "Other"}, inplace=True)
    demo["group"] = demo["diabetes_hx_type"].replace({"1": "Type 1 Diabetes", "2": "Type 2 Diabetes", 1: "Type 1 Diabetes", 2: "Type 2 Diabetes"})
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)
    demo.drop(["diabetes_hx_type", "redcap_event_name", "redcap_repeat_instance", "redcap_repeat_instrument"], axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "medical_history", "field_name"]]    
    
    med = pd.DataFrame(proj.export_records(fields=var))
    med.replace(rep, np.nan, inplace=True)
    med_list = {"diabetes_med___1": "metformin_timepoint",
                "diabetes_med___2": "insulin_timepoint",
                "diabetes_med_other___1": "tzd_timepoint",
                "diabetes_med_other___2": "glp1_agonist_timepoint",
                "diabetes_med_other___3": "sglti_timepoint",
                "diabetes_med_other___4": "other_diabetes_med_timepoint",
                "htn_med_type___1": "ace_inhibitor",
                "htn_med_type___2": "angiotensin_receptor_blocker",
                "htn_med_type___3": "beta_blocker",
                "htn_med_type___4": "ca_channel_blocker",
                "htn_med_type___5": "diuretic",
                "addl_hld_meds___1": "statin",
                "addl_hld_meds___2": "fibrates",
                "addl_hld_meds___3": "niacin",
                "meds_weight_type___1": "topiramate",
                "meds_weight_type___2": "phentermine",
                "mra_med": "mra",
                "uric_acid_med": "uric_acid_med"
                }

    og_names = list(med_list.keys())
    med = med[["subject_id"] + og_names]
    med.rename(med_list, axis=1, inplace=True)    
    
    # RAASi
    med = med.assign(raasi_timepoint=np.maximum(pd.to_numeric(
        med["ace_inhibitor"]), pd.to_numeric(med["angiotensin_receptor_blocker"])))
    # Replace 0/1 values with yes/no
    med.iloc[:, 1:] = med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 2: "No", "2": "No", 1: "Yes", "1": "Yes"})

    # --------------------------------------------------------------------------
    # Med Dispense
    # --------------------------------------------------------------------------
    var = ["subject_id","study_med_disp_date"]
    disp = pd.DataFrame(proj.export_records(fields=var))
    disp.replace(rep, np.nan, inplace=True)
    disp["procedure"] = "med_dispense"
    
    disp["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
        inplace=True)
    disp.rename({"study_med_disp_date": "date", "redcap_event_name": "visit"}, axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Physical exam/vitals
    # --------------------------------------------------------------------------
    var = ["subject_id", "vitals_date", "weight", "height", "bmi", "sys_bp", "dys_bp", "pulse"]
    
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "vitals", "field_name"]]  
                                               
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys.replace(rep, np.nan, inplace=True)
    phys.rename({"sys_bp": "sbp", "dys_bp": "dbp", "vitals_date": "date"}, axis=1, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.drop(["bmi_percentile", "bp_arm", "bp_position", "bp_rest_yn", "sys_bp_1", "sys_bp_2", "sys_bp_3", 
                "dys_bp_1","dys_bp_2", "dys_bp_3", "pulse_1", "pulse_2", "pulse_3", "vitals_norm", "vitals_no"], axis=1, inplace=True)

    phys["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
        inplace=True)
    phys["procedure"] = "physical_exam"
    # set to kidney_biopsy for specific events
    phys.loc[
        phys["redcap_event_name"].isin(["v3_arm_1", "v8_arm_1"]),
        "procedure"
    ] = "kidney_biopsy"
    # then rename
    phys = phys.rename(columns={"redcap_event_name": "visit"})
    #phys.drop(["redcap_event_name"], axis=1, inplace=True)
    
    # --------------------------------------------------------------------------
    # Screening Lab results
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "screening_labs", "field_name"]]    
    screen = pd.DataFrame(proj.export_records(fields=var))
    screen.replace(rep, np.nan, inplace=True)
    screen.rename({ "date_of_screen": "date"}, axis=1, inplace=True)
    screen["procedure"] = "labs"
    screen["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
        inplace=True)
    screen = screen.rename(columns={"redcap_event_name": "visit", "screen_urine_acr": "acr_u"})
    
    screen.drop(["time_of_screen_blood", "screen_egfr", "time_of_screen_urine", "screen_pregnant"], axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Lab results
    # --------------------------------------------------------------------------
    var = ["subject_id", "phys_date", "vitals_date"] + [v for v in meta.loc[meta["form_name"]
                                               == "labs", "field_name"]]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs.replace(rep, np.nan, inplace=True)
    labs['date'] = labs['phys_date'].fillna(labs['vitals_date'])
    labs.rename({ "labs_hg": "hemoglobin", "labs_pltct": "pltct", "serum_sodium": "sodium_s", 
                  "serum_creatinine": "creatinine_s", "totprot_base": "tot_prot", "spotalbumin": "microalbumin_u", 
                  "spotcreatinine": "creatinine_u", "labs_bun": "bun"}, axis=1, inplace=True)
    labs.drop(["wbc", "rbc", "mcv", "mch", "rdw", "mpv", "immgran", "chol_fractions", "ua_color", "ua_norm", 
              "ua_specgrav", "ua_ph", "ua_leukest", "ua_nitrate", "ua_protein", "ua_glucose", "ua_ketone", 
              "ua_urobilinogen", "ua_bilirubin", "ua_blood", "urineculture", "uric_acid", "phys_date", "vitals_date"], axis=1, inplace=True)            
                  
    labs["procedure"] = "labs"
    labs["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
        inplace=True)
    labs = labs.rename(columns={"redcap_event_name": "visit"})
    
    # --------------------------------------------------------------------------
    # Kidney Hemodynamic Outcomes
    # --------------------------------------------------------------------------

    # var = ["subject_id", "phys_date", "vitals_date"] + [v for v in meta.loc[meta["form_name"]
    #                                            == "kidney_hemodynamic_outcomes", "field_name"]]
    # kidney_outcomes = pd.DataFrame(proj.export_records(fields=var))
    # kidney_outcomes.replace(rep, np.nan, inplace=True)
    # kidney_outcomes['date'] = kidney_outcomes['phys_date'].fillna(kidney_outcomes['vitals_date'])
    # #kidney_outcomes.rename({"serum_creatinine": "creatinine_s", "urine_albumin": "screen_urine_acr", "creatinine_u": "creatinine_u"}, axis=1, inplace=True)
    # kidney_outcomes["procedure"] = "renal_clearance_testing"
    # # kidney_outcomes["redcap_event_name"].replace(
    # #     {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         # "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         # "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         # "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
    # #     inplace=True)
    # kidney_outcomes.drop(["phys_date", "vitals_date"], axis=1, inplace=True)            
    # kidney_outcomes = kidney_outcomes.rename(columns={"redcap_event_name": "visit"})
    # 
    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "kidney_biopsy", "field_name"]]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    biopsy.replace(rep, np.nan, inplace=True)
    biopsy.drop([col for col in biopsy.columns if '_yn' in col] +
                [col for col in biopsy.columns if 'procedure_' in col] +
                [col for col in biopsy.columns if '_image_' in col],
                axis=1, inplace=True)
    biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)

    biopsy["procedure"] = "kidney_biopsy"
    biopsy["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
        inplace=True)
    biopsy = biopsy.rename(columns={"redcap_event_name": "visit"})
    
    # --------------------------------------------------------------------------
    # MRI Outcomes
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_boldasl_mri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri.replace(rep, np.nan, inplace=True)
    mri.columns = mri.columns.str.replace(
        r"mri_", "", regex=True)
    mri.rename({"volume_right": "right_kidney_volume_ml",
                "volume_left": "left_kidney_volume_ml",
                "asl_left": "pcasl3d_left","asl_right": "pcasl3d_right", "mri_date": "date"},
               axis=1, inplace=True)
    mri.drop([col for col in mri.columns if '_outcomes' in col],
              axis=1, inplace=True)

    mri["procedure"] = "bold_mri"
    mri["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
        inplace=True)
    mri = mri.rename(columns={"redcap_event_name": "visit"})
    #print(mri)

    # --------------------------------------------------------------------------
    # Renal Clearance testing
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_renal_clearance_testing", "field_name"]]
    rct = pd.DataFrame(proj.export_records(fields=var))
    rct.replace(rep, np.nan, inplace=True)
    rename = {"gfr_raw": "gfr_raw_plasma_urine", "gfr_bsa": "gfr_bsa_plasma_urine",
              "erpf_raw": "erpf_raw_plasma_urine", "erpf": "erpf_bsa_plasma_urine",
              "gfr_15mgmin": "gfr_raw_plasma", "gfrbsa": "gfr_bsa_plasma",
              "erpf_pah_85": "erpf_raw_plasma", "erpfbsa": "erpf_bsa_plasma",
              "pah_bsa": "pah_bsa_plasma_urine", "pahbsa": "pah_clear_bsa",
              "pahcl_12_8mgmin": "pah_clear_abs"}
    rct.rename(rename, axis=1, inplace=True)


    rct.drop([col for col in rct.columns if '_time' in col] +
              [col for col in rct.columns if '_yn' in col],
              axis=1, inplace=True)

    # # Calculate variables
    # rct_vars = ["gfr_raw_plasma", "erpf_raw_plasma",
    #             "totprot_base", "map"]
    # rct[rct_vars] = rct[rct_vars].apply(pd.to_numeric, errors='coerce')
    # rct["erpf_raw_plasma_seconds"] = rct["erpf_raw_plasma"] / 60
    # rct["gfr_raw_plasma_seconds"] = rct["gfr_raw_plasma"] / 60
    # # Filtration Fraction
    # rct["ff"] = rct["gfr_raw_plasma"] / rct["erpf_raw_plasma"]
    # # Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
    # rct["kfg"] = np.select(
    #     [rct["diabetes_hx"].eq("1"), rct["diabetes_hx"].eq("0")], [0.1012, 0.1733])
    # # Filtration pressure across glomerular capillaries
    # rct["deltapf"] = (rct["gfr_raw_plasma"] / 60) / rct["kfg"]
    # # Plasma protein mean concentration
    # rct["cm"] = (rct["totprot_base"] / rct["ff"]) * \
    #     np.log(1 / (1 - rct["ff"]))
    # # Pi G (Oncotic pressure)
    # rct["pg"] = 5 * (rct["cm"] - 2)
    # # Glomerular Pressure
    # rct["glomerular_pressure"] = rct["pg"] + rct["deltapf"] + 10
    # # Renal Blood Flow
    # rct["rbf"] = (rct["erpf_raw_plasma"]) / (1 - rct["hct_210"] / 100)
    # rct["rbf_seconds"] = (rct["erpf_raw_plasma_seconds"]
    #                       ) / (1 - rct["hct_210"] / 100)
    # # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    # rct["rvr"] = rct["map"] / rct["rbf"]
    # # Efferent Arteriolar Resistance
    # rct["re"] = (rct["gfr_raw_plasma_seconds"]) / (rct["kfg"] *
    #                                                (rct["rbf_seconds"] - (rct["gfr_raw_plasma_seconds"]))) * 1328
    # # Afferent Arteriolar Resistance
    # rct["ra"] = ((rct["map"] - rct["glomerular_pressure"]) /
    #              rct["rbf_seconds"]) * 1328
    # rct.loc[~(rct['ra'] > 0), 'ra'] = np.nan
    # Reduce rct dataset
    # rct = rct[["subject_id", "ff", "kfg", "deltapf", "cm", "pg",
    #            "glomerular_pressure",
    #            "pah_raw", "pah_sd", "pah_cv"] + list(rename.values())]
    rct["procedure"] = "renal_clearance_testing"
    
    rct["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
        inplace=True)
    rct = rct.rename(columns={"redcap_event_name": "visit"})
    #print(rct)

    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    demo.dropna(thresh=5, axis=0, inplace=True)
    med.dropna(thresh=3, axis=0, inplace=True)
    disp.dropna(thresh=4, axis=0, inplace=True)
    phys.dropna(thresh=5, axis=0, inplace=True)
    labs.dropna(thresh=5, axis=0, inplace=True)
    # kidney_outcomes.dropna(thresh=5, axis=0, inplace=True)
    biopsy.dropna(thresh=5, axis=0, inplace=True)
    mri.dropna(thresh=5, axis=0, inplace=True)
    rct.dropna(thresh=4, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    # df = pd.concat([med, phys], join='outer', ignore_index=True)
    df = pd.concat([phys, labs], join='outer', ignore_index=True)
    # df = pd.concat([df, kidney_outcomes], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.concat([df, disp], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.concat([df, rct], join='outer', ignore_index=True)
    df = pd.merge(df, demo, on='subject_id', how="outer")
    df = pd.merge(df, med, on='subject_id', how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    print("UNIQUE VISITS: ", df["visit"].unique())

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    #df.rename({"study_visit": "visit"}, axis=1, inplace=True)
    df["study"] = "RPC2"
    id_cols = ["subject_id", "mrn", "study", "dob", "sex", "race", "ethnicity", "visit", "procedure", "group"]
    cols = ['dob', 'group', 'sex', 'race', 'ethnicity', 'mrn']

    df[cols] = (
        df
        .sort_values(['subject_id'])   # important for reproducibility
        .groupby('subject_id')[cols]
        .transform(lambda x: x.ffill().bfill())
    )

    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # Sort
    df.sort_values(["subject_id", "visit", "procedure"], inplace=True)
    # Rename subject identifier
    df.rename({"subject_id": "record_id"}, axis=1, inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data

    datecheck = df.copy()
    fields = ["record_id", "visit", "procedure", "date"]
    mask_fsoc = datecheck["date"].replace("", np.nan).notna()
    sub_fsoc = datecheck[mask_fsoc][fields]
    # list of record_ids that have failed screening
    drop_ids = ["RPC-10", "RPC-12", "RPC-14", "RPC-19", "RPC-22", "RPC-23", "RPC-24", 
                "RPC-30", "RPC-32", "RPC-33", "RPC-34", "RPC-35", "RPC-36", "RPC-39", 
                "RPC-40", "RPC-42"]  
    df = df[~df['record_id'].isin(drop_ids)]


    # sub_fsoc.to_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/rpc2.csv", index=False)
    return df


