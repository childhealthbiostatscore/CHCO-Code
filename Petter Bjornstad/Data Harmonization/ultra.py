"""
This code is designed to pull data from the ULTRA-T2D REDCap project and output data in a long format with one row per study procedure per visit.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_ultra():
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
    token = tokens.loc[tokens["Study"] == "ULTRA-T2D", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Columns to drop
    redcap_cols = ["redcap_event_name",
                   "redcap_repeat_instrument", "redcap_repeat_instance"]
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    var =  [v for v in meta.loc[meta["form_name"] == "demographics", "field_name"]]
    # # Export
    demo = pd.DataFrame(proj.export_records(fields=var))
    # print(type(demo))       # Should be <class 'list'> or <class 'pandas.DataFrame'> depending on your wrapper
    # print(len(demo))        # Number of records exported
    # print(demo[:3]) 
    # Replace missing values
    if demo is None:
        print("Warning: No records returned. Skipping demo processing.")
    else:
        demo.replace(rep, np.nan, inplace=True)    
        demo["group"] = "Type 2 Diabetes"
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
        demo["participation_status"] = demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"})#, inplace=True)
        demo["sex"] = "Male"
    # --------------------------------------------------------------------------
    # Medical History
    # --------------------------------------------------------------------------

    var =  [v for v in meta.loc[meta["form_name"] == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    med.drop(redcap_cols, axis=1, inplace=True)
    # Replace missing values
    med = med.replace(rep, np.nan)

    med = med.rename(
          {"diabetes_meds___1": "metformin_timepoint", "diabetes_meds___2": "insulin_timepoint", "diabetes_meds___3": "tzd_timepoint", 
           "diabetes_meds___4": "glp1_agonist_timepoint", "diabetes_meds___5": "sglti_timepoint"}, axis=1)#, inplace=True)
    #TO DROP: diabetes_meds___6, diabetes_meds___7, insulin_type, other_cvd, cvd_other, cardio_meds_yes, medhx_met, met_hx___3, met_hx___5, met_other, t2d_less_21
    med = med.drop(columns = ["diabetes_meds___6", "diabetes_meds___7", "insulin_type", "other_cvd", "cvd_other", "cardio_meds_yes", "medhx_met", "met_hx___3", 
                              "met_hx___5", "met_other", "t2d_less_21", "insulin_type___1", "insulin_type___2", "insulin_type___3", "diabetes_med_other", 
                              "insulin_other", "hypertension_med_type___7", "cvd_type___5", "hypertension_other", "cardio_meds", "ua_med_prev_when",
                              "ua_med_prev_dur"], errors='ignore')
    med = med.drop(columns = [c for c in med.columns if c.startswith("incl")], errors='ignore')
    #Metformin
    # med.rename({"diabetes_med___1": "metformin_timepoint"},
    #            axis=1, inplace=True)
    #Insulin
    # med["cvd_type"].replace(
    #     {1: "Congestive Heart Failure", 2: "Myocardial Infarction", 3: "Coronary Artery Disease", 4: "Peripheral Vascular Disease", 5: "Other"}, inplace=True)
    # Hypertension Y/N
    # med["met_hx"].replace(
    #     {1: "Hyperlipidemia", 3: "Obesity", 4: "Sleep Apnea", 5: "Other"}, inplace=True)
    med["med_hx_hypertension"] = med["med_hx_hypertension"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"})#, inplace=True)
    med["hypertension_meds"] = med["hypertension_meds"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"})#, inplace=True)
    
    # med = med.rename({"med_hx_hypertension": "hypertension", "diabetes_diag":"diabetes_dx_date"})
    med = med.rename(
        {"hypertension_med_type___1" : "ace_inhibitor", "hypertension_med_type___2": "angiotensin_receptor_blocker", "hypertension_med_type___3": "beta_blocker", 
         "hypertension_med_type___4": "calcium_channel_blocker", "hypertension_med_type___5": "diuretic", "hypertension_med_type___6": "statin"}, axis=1)
    med["medhx_cvd"] = med["medhx_cvd"].replace({0: "No", "0": "No", 1: "Yes", "1": "Yes"})#, inplace=True)
    med = med.rename({
        "cvd_type___1": "congestive_heart_failure", "cvd_type___2": "myocardial_infarction", "cvd_type___3": "coronary_artery_disease", "cvd_type___4": "peripheral_vascular_disease"}, axis=1)
    med = med.rename(
        {
            "met_hx___1": "hyperlipidemia", "met_hx___4": "apnea"}, axis = 1)
    med = med.rename({ "metdx_meds": "hx_met_meds"}, axis = 1)
    med["procedure"] = "screening"
    meds_yn = ["metformin_timepoint", "insulin_timepoint", "tzd_timepoint", "glp1_agonist_timepoint", "sglti_timepoint", "ace_inhibitor", "angiotensin_receptor_blocker", "beta_blocker", 
               "calcium_channel_blocker", "diuretic", "statin", "congestive_heart_failure", "myocardial_infarction", "coronary_artery_disease", "peripheral_vascular_disease", "medhx_cvd", 
               "hyperlipidemia", "apnea", "hx_met_meds", "met_meds_yes"]
    med[meds_yn] = med[meds_yn].apply(lambda col: col.map(
         lambda x: "Yes" if x == "1" or x is True else ("No" if pd.notna(x) and x != "" else x)
         ))
    
    
    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var =  [v for v in meta.loc[meta["form_name"] == "physical_exam", "field_name"]]

    phys = pd.DataFrame(proj.export_records(fields=var))
    phys = phys.drop(redcap_cols, axis=1)
    phys = phys.drop(columns = ["pe_normal", "pe_abnormal", "physical_exam_complete"], errors='ignore')

    # Replace missing values
    phys = phys.replace(rep, np.nan)
    phys = phys.rename({"pe_sbp": "sbp", "pe_dbp": "dbp", "pe_date":"date", "pe_height":"height", "pe_weight":"weight", "pe_bmi":"bmi", "pe_waist":"waistcm", "pe_hip": "hipcm"}, axis=1)
    phys["procedure"] = "physical_exam"

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    
    screen = screen.replace(rep, np.nan)#, inplace=True)  # Replace missing values
    screen.drop(redcap_cols + ['prescreen_a1c_date'],
                axis=1, inplace=True)
    screen["g6pd"] = screen["g6pd"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"})#, inplace=True)
    screen.rename({"a1c": "hba1c", "screening_labs_date":"date_of_screen"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"

    # --------------------------------------------------------------------------
    # Vitals
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "study_visit_vitalslabs", "field_name"]]
    vital = pd.DataFrame(proj.export_records(fields=var))
    if not vital.empty and "visit" in vital.columns:
        vital = vital.replace(rep, np.nan)#, inplace=True)  # Replace missing values
    vital.drop(redcap_cols + ["sv_vitals_yn", "pilabs_yn"],
                axis=1, inplace=True)
    vital.rename({"studyvisit_type":"visit", "sv_date":"date", "labs_s_creatinine": "creatinine_s", "labs_u_creatinine":"creatinine_u", "labs_u_microalbumin":"microalbumin_u", "labs_cystatinc":"cystatin_c_s"},
                  axis=1, inplace=True)
    vital.columns = vital.columns.str.replace(
        r"vitals_", "", regex=True)
    vital.columns = vital.columns.str.replace(
        r"labs_", "", regex=True)
    vital.columns = vital.columns.str.replace(
        r"pilabs_", "", regex=True)
    print("studyvisit_type unique values:", vital["visit"].unique())
    print("studyvisit_type value counts:\n",
          vital["visit"].value_counts(dropna=False))
    vital["procedure"] = "labs"

    # --------------------------------------------------------------------------
    # Imaging
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "imaging", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri.replace(rep, np.nan, inplace=True)  # Replace missing values
    mri.drop(redcap_cols, axis=1, inplace=True)#+ ["mri_cardio", "mri_abdo",
                            #"mri_aortic", "study_visit_mri"],
    print([c for c in mri.columns if 'date' in c.lower() or 'imaging' in c.lower()]) 
    mri.columns = mri.columns.str.replace(
        r"mri_|visit_", "", regex=True)
    mri.rename({"lvsv": "lv_stroke_volume", "rvsv" : "rv_stroke_volume", "rvco": "rv_cardiac_output",
                "lvco": "lv_cardiac_output", "myo_mass_dias" : "lved_mass", "myo_mass_syst": "lves_mass",
                "lv_myo_mass_dias" : "lv_myo_mass_diast",
                "imaging_hr": "lv_hr", "imaging_date": "date",
                "af_pwv_xcor3": "af_pwv",
                "radial_peak": "grs", "circum_peak": "gcs", "long_peak": "gls"}, axis=1, inplace=True)
    print("af_pwv PULLED?")                                                       
    print("af_pwv" in mri.columns)                            
    print("TOTAL NOT NULLS IN date:")                                           
    print(mri["date"].notna().sum())
     
                
    mri["procedure"] = "cardio_abdominal_mri"

    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=3, axis=0, inplace=True)
    vital.dropna(thresh=3, axis=0, inplace=True)
    phys.dropna(thresh=3, axis=0, inplace=True)
    screen.dropna(thresh=3, axis=0, inplace=True)
    demo.dropna(thresh=3, axis=0, inplace=True)
    mri.dropna(thresh=4, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([med, vital], join='outer', ignore_index=True)
    #df = pd.concat([df, vital], join='outer', ignore_index=True)
    df = pd.concat([df, phys], join='outer', ignore_index=True)
    df = pd.concat([df, screen], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)

    df = pd.merge(df, demo, how="outer")
#     df = df.loc[:, ~df.columns.str.startswith('redcap_')]
#     df = df.copy()

#     # --------------------------------------------------------------------------
#     # Reorganize
#     # --------------------------------------------------------------------------

    #df.rename({"study_visit": "visit"}, axis=1, inplace=True)
    df["study"] = "ULTRA"
    id_cols = ["record_id", "study"] + \
        [ "dob", "participation_status", "mrn"] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # Change study visit names
    # df["visit"].replace({np.nan: "baseline", '1': "baseline",
    #                      '2': "3_months_post_surgery", '3': "12_months_post_surgery"}, inplace=True)
    df["visit"] = df["visit"].replace({
        np.nan: "baseline",
        "": "baseline",
        1: "baseline",
        '1': "baseline",
        2: "1_week",
        '2': "1_week",
    })
    df["visit"] = pd.Categorical(df["visit"],
                                 categories=["baseline", "1_week"],
                                 ordered=True)
    print("visit unique values:", df["visit"].unique())
    print("visit value counts:\n",
          df["visit"].value_counts(dropna=False))
    # Fix subject IDs
    #df["subject_id"] = df["subject_id"].str.replace(r"2D-", "_", regex=True)
    # Sort
    df.sort_values(["record_id", "visit", "procedure"], inplace=True)
    # Rename subject identifier
    #df.rename({"subject_id": "record_id"}, axis=1, inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data
    return df
