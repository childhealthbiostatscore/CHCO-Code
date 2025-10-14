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

    dem_cols = ["record_id", "dob", "race", "ethnicity", "participation_status", "mrn"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo["group"] = "Type 2 Diabetes"
    demo.drop(redcap_cols, axis=1, inplace=True)
    demo.rename({"gender": "sex", "mr_number": "mrn"},
                inplace=True, axis=1)
    # Race columns combined into one
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

    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)

    # --------------------------------------------------------------------------
    # Medical History
    # --------------------------------------------------------------------------

    var = ["record_id", "diabetes_diag", "med_hx_hypertension"] #"insulin_type",met_hx
           #"cvd_type", "met_hx", "med_hx_hypertension"] #"diabetes_meds"]
    med = pd.DataFrame(proj.export_records(fields=var))
    med.drop(redcap_cols, axis=1, inplace=True)
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)

    # med["diabetes_meds"].replace(
    #      {"diabetes_meds___1": "Metformin", "diabetes_meds___2": "Insulin", "diabetes_meds___3": "Thiazolinediones (TZDs)", "diabetes_meds___4": "GLP-1 agonists", "diabetes_meds___5": "SGLT-1/2 inhibitors", "diabetes_meds___6": "Other", "diabetes_meds___7":"None"}, inplace=True)
    # Metformin
    # med["insuline_type"].replace(
    #     {1: "Long acting", 2: "Short acting", 3: "Other"}, inplace=True)
    # med.rename({"diabetes_med___1": "metformin_timepoint"},
    #            axis=1, inplace=True)
    # Insulin
    # med["cvd_type"].replace(
    #     {1: "Congestive Heart Failure", 2: "Myocardial Infarction", 3: "Coronary Artery Disease", 4: "Peripheral Vascular Disease", 5: "Other"}, inplace=True)
    # # Hypertension Y/N
    # med["met_hx"].replace(
    #     {1: "Hyperlipidemia", 3: "Obesity", 4: "Sleep Apnea", 5: "Other"}, inplace=True)
    med["med_hx_hypertension"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"med_hx_hypertension": "hypertension", "diabetes_diag":"diabetes_dx_date"})


    med["procedure"] = "screening"
    
    
    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["record_id", "pe_date", "pe_height", "pe_weight",
           "pe_bmi", "pe_sbp", "pe_dbp", "pe_waist", "pe_hip"]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys.drop(redcap_cols, axis=1, inplace=True)
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"pe_sbp": "sbp", "pe_dbp": "dbp", "pe_date":"date", "pe_height":"height", "pe_weight":"weight", "pe_bmi":"bmi", "pe_waist":"waistcm", "pe_hip": "hipcm"}, inplace=True, axis=1)
    phys["procedure"] = "physical_exam"

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    screen.replace(rep, np.nan, inplace=True)  # Replace missing values
    screen.drop(redcap_cols + ['prescreen_a1c', 'prescreen_a1c_date'],#, "screening_labs_complete"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|_of_screen", "", regex=True)
    screen.rename({"a1c": "hba1c"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"

    # --------------------------------------------------------------------------
    # Vitals
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "study_visit_vitalslabs", "field_name"]]
    vital = pd.DataFrame(proj.export_records(fields=var))
    vital.replace(rep, np.nan, inplace=True)  # Replace missing values
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
    
    vital["procedure"] = "labs"

    # --------------------------------------------------------------------------
    # Imaging
    # --------------------------------------------------------------------------

    # var = ["record_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
    #                                                            == "imaging", "field_name"]]
    # mri = pd.DataFrame(proj.export_records(fields=var))
    # mri.replace(rep, np.nan, inplace=True)  # Replace missing values
    # mri.drop(redcap_cols + ["mri_cardio", "mri_abdo",
    #                         "mri_aortic", "study_visit_mri"],
    #          axis=1, inplace=True)
    # mri.columns = mri.columns.str.replace(
    #     r"mri_|visit_", "", regex=True)
    # mri["procedure"] = "cardio_abdominal_mri"

    
    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=3, axis=0, inplace=True)
    vital.dropna(thresh=3, axis=0, inplace=True)
    phys.dropna(thresh=3, axis=0, inplace=True)
    screen.dropna(thresh=3, axis=0, inplace=True)
    demo.dropna(thresh=3, axis=0, inplace=True)
    

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([med, vital], join='outer', ignore_index=True)
    #df = pd.concat([df, vital], join='outer', ignore_index=True)
    df = pd.concat([df, phys], join='outer', ignore_index=True)
    df = pd.concat([df, screen], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    #df.rename({"study_visit": "visit"}, axis=1, inplace=True)
    df["study"] = "ULTRA"
    id_cols = ["record_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # Change study visit names
    df["visit"].replace({np.nan: "baseline", '1': "baseline",
                         '2': "3_months_post_surgery", '3': "12_months_post_surgery"}, inplace=True)
    df["visit"] = pd.Categorical(df["visit"],
                                 categories=["baseline", "3_months_post_surgery",
                                 "12_months_post_surgery"],
                                 ordered=True)
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
