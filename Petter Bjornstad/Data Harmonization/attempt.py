"""
This code is designed to pull data from the ATTEMPT REDCap project and output data in a long format with one row per study procedure per visit.
"""
__author__ = ["Ye Ji Choi"]
__credits__ = ["Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Ye Ji Choi"
__email__ = "yejichoi@uw.edu"
__status__ = "Dev"




# Function to clean and structure REDCap data
def clean_attempt():
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
    tokens = pd.read_csv(
        "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "ATTEMPT", "Token"].iloc[0]
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
    
    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["subject_id", "consent_date", "mrn", "dob", "sex", "race", "race_other", "ethnicity", "participation_status", "t1d", "t1d_date"]    
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    demo.replace(rep, np.nan, inplace=True)
    
    # Map categorical variables
    demo["sex"].replace({"1": "Female", "2": "Male", "3": "Other"}, inplace=True)
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)
    
    # Combine checkbox variables
    demo = combine_checkboxes(demo, base_name="race", levels=["American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
    demo = combine_checkboxes(demo, base_name="ethnicity", levels=["Hispanic", "Non-Hispanic", "Unknown/Not Reported"])
    
    # Create group variable based on diabetes status
    demo["group"] = demo["t1d"].replace({"1": "Type 1 Diabetes"})
    demo.rename({"t1d_date": "diabetes_dx_date"}, axis=1, inplace=True)
    demo.drop("redcap_event_name", axis=1, inplace=True)
    dem_cols = ["subject_id", "consent_date", "mrn", "dob", "sex", "race", "race_other", "ethnicity", "participation_status", "group", "diabetes_dx_date"]    
    
    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "epic_meds", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    med.replace(rep, np.nan, inplace=True)
    med.replace({"0": "No", "1": "Yes"}, inplace=True)
    med.rename({"study_visit_date": "date", "redcap_event_name": "visit"}, axis=1, inplace=True)
    med["procedure"] = "medications"
    med["visit"].replace({"screening_visit_arm_1": "baseline", "visit_2_arm_1": "baseline", "visit_3_arm_1": "4_months_post"}, inplace=True)

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------
    var = ["subject_id", "date_clamp", "sbp_clamp", "dbp_clamp", "clamp_height", "weight_clamp"]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    clamp.replace(rep, np.nan, inplace=True)
    clamp.rename({"sbp_clamp": "sbp", "dbp_clamp": "dbp", "clamp_height": "height", "weight_clamp": "weight", "date_clamp": "date", "redcap_event_name": "visit"}, axis=1, inplace=True)
    clamp["procedure"] = "clamp"
    clamp["visit"].replace({"screening_visit_arm_1": "baseline", "visit_2_arm_1": "baseline", "visit_3_arm_1": "4_months_post"}, inplace=True)

    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "kidney_biopsy", "field_name"]]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    biopsy.replace(rep, np.nan, inplace=True)
    biopsy.drop([col for col in biopsy.columns if '_yn' in col] +
                [col for col in biopsy.columns if 'procedure_' in col] +
                ["core_diagnostic", "core_hypo_cryo", "core_oct", "core_rna"],
                axis=1, inplace=True)
    biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"labs_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"vitals_", "", regex=True)
    biopsy.rename({"redcap_event_name": "visit"}, axis=1, inplace=True)
    biopsy["procedure"] = "kidney_biopsy"
    biopsy["visit"].replace({"screening_visit_arm_1": "baseline", "visit_2_arm_1": "baseline", "visit_3_arm_1": "4_months_post"}, inplace=True)
    
    # --------------------------------------------------------------------------
    # MRI
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "study_visit_boldasl_mri", "field_name"]]
    bold_mri = pd.DataFrame(proj.export_records(fields=var))
    bold_mri.replace(rep, np.nan, inplace=True)
    bold_mri.rename({"mri_lab_date": "date", "redcap_event_name": "visit"}, axis=1, inplace=True)
    bold_mri["procedure"] = "bold_mri"
    bold_mri["visit"].replace({"screening_visit_arm_1": "baseline", "visit_2_arm_1": "baseline", "visit_3_arm_1": "4_months_post"}, inplace=True)


    # --------------------------------------------------------------------------
    # Safety Labs
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "safety_labs", "field_name"]]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs.replace(rep, np.nan, inplace=True)
    labs.rename({"hgb": "hemoglobin", "hct": "hematocrit", "plt": "pltct", "creat": "creat_s", "redcap_event_name": "visit"}, axis=1, inplace=True)
    labs["visit"].replace({"screening_visit_arm_1": "baseline", "visit_2_arm_1": "baseline", "visit_3_arm_1": "4_months_post"}, inplace=True)
    
    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    demo.dropna(thresh=9, axis=0, inplace=True)
    med.dropna(thresh=9, axis=0, inplace=True)
    clamp.dropna(thresh=5, axis=0, inplace=True)
    bold_mri.dropna(thresh=4, axis=0, inplace=True)
    labs.dropna(thresh=5, axis=0, inplace=True)
    biopsy.dropna(thresh=12, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([clamp, med], join='outer', ignore_index=True)
    df = pd.concat([df, bold_mri], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.concat([df, labs], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df["study"] = "ATTEMPT"
    id_cols = ["subject_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
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
