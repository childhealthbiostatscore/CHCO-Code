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
        base_data_path = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/"
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

    dem_cols = ["subject_id", "age_current", "date_of_consent", "mr_number", "dob", "gender", "race", "ethnicity", "participation_status", "diabetes_hx_type"]
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    demo.replace(rep, np.nan, inplace=True)
    demo.rename({"mr_number": "mrn", "gender": "sex", "age_current": "age"}, axis=1, inplace=True)
    demo = combine_checkboxes(demo, base_name="race", levels=["American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
    demo = combine_checkboxes(demo, base_name="ethnicity", levels=["Hispanic", "Non-Hispanic", "Unknown/Not Reported"])
    demo["sex"].replace({1: "Male", 0: "Female", 2: "Other",
                        "1": "Male", "0": "Female", "2": "Other"}, inplace=True)
    demo["group"] = demo["diabetes_hx_type"].replace({"1": "Type 1 Diabetes", "2": "Type 2 Diabetes"})
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)
    demo["age"] = pd.to_numeric(demo["age"], downcast='integer', errors='coerce')
    

    print(demo)

    

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------
    var = ["subject_id", "vitals_date", "diabetes_med", "diabetes_med_other", "htn_med_type", "addl_hld_meds", "insulin_med", "kidneybx_1", "kidneybx_2", "study_med_date"]
    med = pd.DataFrame(proj.export_records(fields=var))
    med.replace(rep, np.nan, inplace=True)
    med.rename({"diabetes_med_other": "sglt2i_timepoint", "insulin_med": "insulin_med_timepoint", "study_med_date" : "date"}, axis=1, inplace=True)
    med.iloc[:, 1:] = med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 2: "No", "2": "No", 1: "Yes", "1": "Yes"})
    med["procedure"] = "medications"
    med = med.rename(columns={"redcap_event_name": "visit"})
    med["visit"] = "baseline"
    print(med)
    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------
    var = ["subject_id", "phys_date", "weight", "height", "bmi", "sys_bp", "dys_bp", "pulse"]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys.replace(rep, np.nan, inplace=True)
    phys.rename({"bp_systolic": "sbp", "bp_diastolic": "dbp", "phys_date": "date"}, axis=1, inplace=True)
    phys["procedure"] = "physical_exam"
    phys["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "baseline", "v2_gfr_mri_arm_1": "baseline", "v4_arm_1": "Post-biopsy vist", "v7_gfr_mri_arm_1": "Post-biopsy GFR MRI", "v61_med_dispense_arm_1": "Med Dispense 1",  "v62_med_dispense_arm_1": "Med Dispense 2"}, 
        inplace=True)
    phys = phys.rename(columns={"redcap_event_name": "visit"})
    print(phys)
    #phys.drop(["redcap_event_name"], axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Vitals
    # --------------------------------------------------------------------------
    var = ["subject_id", "vitals_date", "weight", "height", "bmi"]
    vit = pd.DataFrame(proj.export_records(fields=var))
    vit.replace(rep, np.nan, inplace=True)
    vit.rename({"vitals_date" : "date"}, axis=1, inplace=True)
    vit["procedure"] = "physical_exam"
    # --------------------------------------------------------------------------
    # Lab results
    # --------------------------------------------------------------------------
    var = ["subject_id", "date_of_screen", "serum_creatinine", "screen_urine_acr", "microalbumin_u", "creatinine_u", "creatinine_s", "cystatin_c_s","hba1c", "chol_base", "bun", "hdl", "ldl", "triglycerides"]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs.replace(rep, np.nan, inplace=True)
    labs.rename({ "date_of_screen": "date"}, axis=1, inplace=True)
    labs["procedure"] = "labs"
    labs["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "baseline", "v2_gfr_mri_arm_1": "baseline", "v4_arm_1": "post_biopsy", "v7_gfr_mri_arm_1": "Post-biopsy GFR MRI", "v61_med_dispense_arm_1": "Med Dispense 1",  "v62_med_dispense_arm_1": "Med Dispense 2"}, 
        inplace=True)
    labs = labs.rename(columns={"redcap_event_name": "visit"})
    print(labs)
    
    # --------------------------------------------------------------------------
    # Kidney Hemodynamic Outcomes
    # --------------------------------------------------------------------------

    var = ["subject_id", "vitals_date", "ff", "glomerular_pressure", "rbf", "aff_arteriolar_resistance", "eff_arteriolar_resistance"]#, "gfr_raw", "gfr_ecv_percent", "gfr_ecv_std", "gfr_cv", "gfr_bsa", "gfr_15mgmin", "gfrbsa", "erpf_raw", "erpf", "erpfbsa"]
    kidney_outcomes = pd.DataFrame(proj.export_records(fields=var))
    kidney_outcomes.replace(rep, np.nan, inplace=True)
    #kidney_outcomes.rename({"serum_creatinine": "creatinine_s", "urine_albumin": "screen_urine_acr", "creatinine_u": "creatinine_u"}, axis=1, inplace=True)
    kidney_outcomes["procedure"] = "kidney_hemodynamic_outcomes"
    kidney_outcomes["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "baseline", "v2_gfr_mri_arm_1": "baseline", "v4_arm_1": "post_biopsy", "v7_gfr_mri_arm_1": "Post-biopsy GFR MRI", "v61_med_dispense_arm_1": "Med Dispense 1",  "v62_med_dispense_arm_1": "Med Dispense 2"}, 
        inplace=True)
    kidney_outcomes = kidney_outcomes.rename(columns={"redcap_event_name": "visit"})

    print(kidney_outcomes)
    
    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------
    var = ["subject_id", "bx_date", "gloms", "gloms_gs", "ifta", "mes_index", "pod_nuc_density", "bx_kit_id"]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    biopsy.replace(rep, np.nan, inplace=True)
    biopsy.rename({"bx_date": "date", "bx_kit_id":"kit_id"}, axis=1, inplace=True)
    biopsy["procedure"] = "kidney_biopsy"
    biopsy["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "baseline", "v2_gfr_mri_arm_1": "baseline", "v4_arm_1": "post_biopsy", "v7_gfr_mri_arm_1": "Post-biopsy GFR MRI", "v61_med_dispense_arm_1": "Med Dispense 1",  "v62_med_dispense_arm_1": "Med Dispense 2"}, 
        inplace=True)
    biopsy = biopsy.rename(columns={"redcap_event_name": "visit"})
    print(biopsy)
    
    # --------------------------------------------------------------------------
    # MRI Outcomes
    # --------------------------------------------------------------------------
    var = ["subject_id","mri_date", "asl_right", "asl_left", "adc_right", "adc_left", "bold_r_bl_cortex", "bold_l_bl_cortex"]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri.replace(rep, np.nan, inplace=True)
    mri.rename({"asl_left": "pcasl3d_left","asl_right": "pcasl3d_right", "mri_date": "date"}, axis=1, inplace=True)

    mri["procedure"] = "bold_mri"
    mri["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "baseline", "v2_gfr_mri_arm_1": "baseline", "v4_arm_1": "post_biopsy", "v7_gfr_mri_arm_1": "Post-biopsy GFR MRI", "v61_med_dispense_arm_1": "Med Dispense 1",  "v62_med_dispense_arm_1": "Med Dispense 2"}, 
        inplace=True)
    mri = mri.rename(columns={"redcap_event_name": "visit"})
    print(mri)


    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([demo, med], join='outer', ignore_index=True)
    df = pd.concat([df, phys], join='outer', ignore_index=True)
    df = pd.concat([df, labs], join='outer', ignore_index=True)
    df = pd.concat([df, kidney_outcomes], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    #df.rename({"study_visit": "visit"}, axis=1, inplace=True)
    df["study"] = "RPC2"
    id_cols = ["subject_id", "mrn", "study", "dob", "sex", "race", "ethnicity", "visit", "procedure", "group"]
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
    return df


