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
    
    # REDCap project variables
    tokens = pd.read_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
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
    
    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["subject_id", "date_of_consent", "mr_number", "dob", "gender", "race", "ethnicity", "participation_status", "diabetes_hx_type"]
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    demo.replace(rep, np.nan, inplace=True)
    demo.rename({"mr_number": "mrn", "gender": "sex"}, axis=1, inplace=True)
    demo = combine_checkboxes(demo, base_name="race", levels=["American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
    demo = combine_checkboxes(demo, base_name="ethnicity", levels=["Hispanic", "Non-Hispanic", "Unknown/Not Reported"])
    demo["sex"].replace({1: "Male", 0: "Female", 2: "Other",
                        "1": "Male", "0": "Female", "2": "Other"}, inplace=True)
    demo["group"] = demo["diabetes_hx_type"].replace({"1": "Type 1 Diabetes", "2": "Type 2 Diabetes"})
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)
    
    # Map Yes/No/Other values
    yes_no_map = {"0": "No", "1": "Yes", "2": "Unknown"}

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------
    var = ["subject_id", "visit", "diabetes_med", "diabetes_med_other", "htn_med_type", "addl_hld_meds", "insulin_med"]
    med = pd.DataFrame(proj.export_records(fields=var))
    med.replace(rep, np.nan, inplace=True)
    med.rename({"diabetes_med_other": "sglt2i_timepoint", "insulin_med": "insulin_med_timepoint"}, axis=1, inplace=True)
    med["procedure"] = "medications"
    
    # Physical exam
    var = ["subject_id", "visit", "weight", "height", "bmi", "bp_systolic", "bp_diastolic", "pulse"]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys.replace(rep, np.nan, inplace=True)
    phys.rename({"bp_systolic": "sbp", "bp_diastolic": "dbp"}, axis=1, inplace=True)
    phys["procedure"] = "physical_exam"
    
    # Lab results
    var = ["subject_id", "visit", "serum_creatinine", "urine_albumin", "urine_creatinine", "hba1c", "chol_base", "hdl", "ldl", "triglycerides"]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs.replace(rep, np.nan, inplace=True)
    labs.rename({"serum_creatinine": "creatinine_s", "urine_albumin": "acr_u", "urine_creatinine": "creatinine_u"}, axis=1, inplace=True)
    labs["procedure"] = "labs"
    
    # Kidney Hemodynamic Outcomes
    var = ["subject_id", "visit", "gfr", "rpf", "filtration_fraction", "glomerular_pressure", "rbf", "aff_arteriolar_resistance", "eff_arteriolar_resistance"]
    kidney_outcomes = pd.DataFrame(proj.export_records(fields=var))
    kidney_outcomes.replace(rep, np.nan, inplace=True)
    kidney_outcomes["procedure"] = "kidney_hemodynamic_outcomes"
    
    # Kidney Biopsy
    var = ["subject_id", "visit", "bx_date", "gloms", "gloms_gs", "ifta", "mes_index", "pod_nuc_density"]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    biopsy.replace(rep, np.nan, inplace=True)
    biopsy["procedure"] = "kidney_biopsy"
    
    # MRI Outcomes
    var = ["subject_id", "visit", "mri_date", "asl_right", "asl_left", "adc_right", "adc_left", "bold_r_bl_cortex", "bold_l_bl_cortex"]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri.replace(rep, np.nan, inplace=True)
    mri["procedure"] = "mri_outcomes"
    
    # Merge all datasets
    df = pd.concat([med_hist, med, phys, labs, kidney_outcomes, biopsy, mri], ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df.sort_values(["subject_id", "visit", "procedure"], inplace=True)
    
    return df
