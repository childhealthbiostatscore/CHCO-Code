"""
This code is designed to pull data from the SWEETHEART REDCap project and output data in a long format with one row per study procedure per visit.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_sweetheart():
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
        git_path = "/Users/shivaniramesh/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    else:
        sys.exit(f"Unknown user: please specify root path for this user. (Detected user: {user})")

    tokens = pd.read_csv(base_data_path + "/Data Harmonization/api_tokens.csv")        #"/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "SWEETHEART", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Columns to drop
    #redcap_cols = ["redcap_event_name",
                   #"redcap_repeat_instrument", "redcap_repeat_instance"]
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["record_id", "co_enroll_id", "group", "dob", "sex", 
                 "race", "ethnicity", "participation_status", "mrn"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo["group"].replace({"2": "Type 2 Diabetes", "3": "Obese Control",
                        "4": "Lean Control"}, inplace=True)
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
    # Relevel sex and group
    demo["sex"].replace({1: "Male", 2: "Female", 3: "Other",
                        "1": "Male", "2": "Female", "3": "Other"}, inplace=True)

    demo["group_risk"] = np.where(demo.group.str.contains("lean", case=False), "Low", "High")
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)

    # --------------------------------------------------------------------------
    # Medical History
    # --------------------------------------------------------------------------

    var = ["record_id", "diabetes_dx_date", "diabetes_duration", "diabetes_med",
           "diabetes_med_other", "htn_med_type", "addl_hld_meds", "hypertension"]
    med = pd.DataFrame(proj.export_records(fields=var))

    #med.drop(redcap_cols, axis=1, inplace=True)
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    # SGLT2i (diabetes_med_other___4), RAASi (htn_med_type___1, htn_med_type___2), Metformin (diabetes_med_other___1)
    med = med[["record_id", "diabetes_med_other___3", "htn_med_type___1",
               "htn_med_type___2", "htn_med_type___3", "htn_med_type___5", "diabetes_med___1", "diabetes_med___2", "addl_hld_meds___1", "hypertension"]]
    # SGLT2i
    med["diabetes_med_other___3"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med_other___3": "sglti_timepoint"},
               axis=1, inplace=True)
    # RAASi
    med = med.assign(raasi_timepoint=np.maximum(pd.to_numeric(
        med["htn_med_type___1"]), pd.to_numeric(med["htn_med_type___2"])))
    med["raasi_timepoint"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    # Metformin
    med["diabetes_med___1"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med___1": "metformin_timepoint"},
               axis=1, inplace=True)
    # Insulin
    med["diabetes_med___2"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med___2": "insulin_med_timepoint"},
               axis=1, inplace=True)
    # Hypertension Y/N
    med["hypertension"].replace(
        {2: "No", "2": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med["htn_med_type___1"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"htn_med_type___1": "ace_inhibitor"},
               axis=1, inplace=True)
    med["htn_med_type___3"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"htn_med_type___3": "beta_blocker"},
               axis=1, inplace=True)
    med["htn_med_type___5"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"htn_med_type___5": "diuretic"},
               axis=1, inplace=True)

    # Statin
    med["addl_hld_meds___1"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"addl_hld_meds___1": "statin"},
           axis=1, inplace=True)
    med.rename({"med_date": "date"},
           axis=1, inplace=True)
    med["procedure"] = "medications"
    

    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(
        fields=var))
    
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys.drop(["phys_normal", "phys_abnormal"], axis=1, inplace=True)
    
    phys.rename({"phys_sysbp": "sbp", "phys_diasbp": "dbp"
                 }, inplace=True, axis=1)
    non_nan_count = phys['sbp'].notna().sum()
    #print(f"NONNAN SBP COUNT:", non_nan_count)
    phys.columns = phys.columns.str.replace(r"phys_", "", regex=True)
    phys["procedure"] = "physical_exam"

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))

    screen.replace(rep, np.nan, inplace=True)  # Replace missing values
    screen.drop( ['prescreen_a1c_date', 'prescreen_a1c'],
                axis=1, inplace=True)
    
    screen.rename({ "screen_uacr": "acr_u",
                    "screen_a1c":"hba1c", "screen_creat_s": "creatinine_s", "screen_creat_u":"creatinine_u", "labs_date":"date_of_screen"},
                  axis=1, inplace=True)
    print("NONNAN ACRU:", screen['acr_u'].notna().sum())
    screen.columns = screen.columns.str.replace(r"screen_|_of_screen", "", regex=True)
    screen["procedure"] = "screening"

    # --------------------------------------------------------------------------
    # Imaging
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "imaging", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))

    screen.replace(rep, np.nan, inplace=True)  # Replace missing values
    
    screen.columns = screen.columns.str.replace(
        r"screen_|_of_screen", "", regex=True)
    
    screen["procedure"] = "imaging"




    
    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=3, axis=0, inplace=True)
    demo.dropna(thresh=3, axis=0, inplace=True)
    phys.dropna(thresh=3, axis=0, inplace=True)
    screen.dropna(thresh=3, axis=0, inplace=True)
    

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([phys, med], join='outer', ignore_index=True)
    df = pd.concat([df, screen], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    #df.rename({"study_visit": "visit"}, axis=1, inplace=True)
    df["study"] = "SWEETHEART"
    id_cols = ["record_id", "study"] + \
        dem_cols[1:] + [ "procedure"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    df["record_id"] = "SWHT_" + df["record_id"].astype(str)

    # Change study visit names
    #df["visit"].replace({np.nan: "baseline", '1': "baseline",
                         #'2': "3_months_post_surgery", '3': "12_months_post_surgery"}, inplace=True)
    #df["visit"] = pd.Categorical(df["visit"],
                                 #
    # Fix subject IDs
    #df["record_id"] = df["record_id"].str.replace(r"2D-", "_", regex=True)
    # Sort
    df.sort_values(["record_id", "procedure"], inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data
    return df
