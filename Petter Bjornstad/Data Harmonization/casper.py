"""
This code is designed to pull data from the CASPER REDCap project and output 
data in a "semi-long" format with one row per study procedure, and a visit 
column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_casper():
    # Libraries
    import os
    import sys
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
    token = tokens.loc[tokens["Study"] == "CASPER", "Token"].iloc[0]
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

    dem_cols = ["subject_id", "dob", "diagnosis", 
                "gender", "race", "ethnicity", "participation_status","mr_number"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo["co_enroll_id"] = ""
    demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date", "mr_number": "mrn"},
                inplace=True, axis=1)
    dem_cols[2] = "diabetes_dx_date"
    dem_cols[3] = "sex"
    dem_cols[7] = "mrn"
    # Race columns combined into one
    demo = combine_checkboxes(demo, base_name="race", levels=[
        "American Indian or Alaskan Native",
        "Asian",
        "Hawaiian or Pacific Islander",
        "Black or African American",
        "White",
        "Unknown",
        "Other"])
    # Same for ethnicity
    demo = combine_checkboxes(demo,
                              base_name="ethnicity",
                              levels=["Hispanic or Latino",
                                      "Not Hispanic or Latino",
                                      "Unknown/Not Reported"])
    # Relevel sex and group
    demo["sex"].replace({1: "Male", 0: "Female", 3: "Other",
                        "1": "Male", "0": "Female", "3": "Other"}, inplace=True)
    demo["group"] = "Type 1 Diabetes"
    demo["group_risk"] = np.where(demo.group.str.contains("lean", case=False), "Low", "High")
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    # SGLT2i (diabetes_med_other___4), RAASi (htn_med_type___1, htn_med_type___2), Metformin (diabetes_med_other___1)
    med = med[["subject_id", "diabetes_med_other___4", "htn_med_type___1",
               "htn_med_type___2", "diabetes_med_other___1", "diabetes_med___1",
               "diabetes_med___2"]]
    # SGLT2i
    med["diabetes_med_other___4"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med_other___4": "sglti_timepoint"},
               axis=1, inplace=True)
    # RAASi
    med = med.assign(raasi_timepoint=np.maximum(pd.to_numeric(
        med["htn_med_type___1"]), pd.to_numeric(med["htn_med_type___2"])))
    med.drop(med[['htn_med_type___1', 'htn_med_type___2']],
             axis=1, inplace=True)
    med["raasi_timepoint"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    # Metformin
    med["diabetes_med_other___1"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med_other___1": "metformin_timepoint"},
               axis=1, inplace=True)
    # Insulin
    med = med.assign(insulin_med_timepoint=np.maximum(pd.to_numeric(
        med["diabetes_med___1"]), pd.to_numeric(med["diabetes_med___2"])))
    med.drop(med[['diabetes_med___1', 'diabetes_med___2']],
             axis=1, inplace=True)
    med["insulin_med_timepoint"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med["procedure"] = "medications"
    med["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.drop(["phys_norm", "phys_no", "lmp", "screen_bmi_percentile",
               "male_activity_factor", "fem_activity_factor", "schofield_male",
               "schofield_female"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sys_bp": "sbp", "dys_bp": "dbp",
                 "waist_circumference": "waistcm",
                 "hip_circumference": "hipcm"}, inplace=True, axis=1)
    phys["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    screen.replace(rep, np.nan, inplace=True)
    screen.drop(['a1c_pre', 'a1c_pre_date', "screen_pregnant"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|_of_screen", "", regex=True)
    screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
                   "urine_cre": "creatinine_u", "urine_mab": "microalbumin_u"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"
    # Assume medication review done at screening
    med["date"] = screen["date"]
    screen["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    clamp.replace(rep, np.nan, inplace=True)
    # Format
    clamp.drop(["baseline", "fasting_labs", "bg_labs", "urine_labs", "hct_lab",
                "a1c_clamp_time", "clamp_a1c", "clamp_a1c_date"],
               axis=1, inplace=True)
    clamp.columns = clamp.columns.str.replace(
        r"clamp_", "", regex=True)
    clamp.rename({"cystatin_c": "cystatin_c_s", 
                  "serum_creatinine": "creatinine_s",
                  "urine_mab_baseline": "microalbumin_u",
                  "urine_cre_baseline": "creatinine_u", 
                  "pls": "pulse",
                  "urine_sodium": "sodium_u", 
                  "serum_sodium": "sodium_s",
                  "acr_baseline": "acr_u", 
                  "acr_250": "acr_u_pm", 
                  "total_protein": "tot_protein",
                  "glucose_bl": "urine_glucose_bl"},
                 inplace=True, axis=1)
    clamp["procedure"] = "clamp"
    clamp["visit"] = "baseline"
    clamp["insulin_sensitivity_method"] = "hyperglycemic_clamp"
    clamp["ffa_method"] = "hyperglycemic_clamp"
    # No insulin, c peptide, or FFA
    num_vars = ["d20_infusion", "weight", "hematocrit_minus_5", "hematocrit_90", "hematocrit_120"]
    clamp[num_vars] = clamp[num_vars].apply(
        pd.to_numeric, errors='coerce')
    clamp["gir_190"] = (clamp["d20_infusion"] * 190 / 60) / clamp["weight"] # previously M-value
    clamp["gir_200"] = (clamp["d20_infusion"] * 200 / 60) / clamp["weight"]
    # Hematocrit average
    clamp["hematocrit_avg"] = clamp[["hematocrit_minus_5", "hematocrit_90", "hematocrit_120"]].mean(axis=1)

    # --------------------------------------------------------------------------
    # DXA Scan
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "body_composition_dxa", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    dxa.replace(rep, np.nan, inplace=True)
    dxa.columns = dxa.columns.str.replace(
        r"dexa_", "", regex=True)
    dxa_cols = dxa.columns[2:].to_list()
    dxa.rename(dict(zip(dxa_cols, ["dexa_" + d for d in dxa_cols])),
               axis=1, inplace=True)
    dxa["procedure"] = "dxa"
    dxa["visit"] = "baseline"
    dictionary.loc[dictionary['variable_name'].isin(dxa.columns), 'form_name'] = 'dxa'


    # --------------------------------------------------------------------------
    # Outcomes
    # --------------------------------------------------------------------------

    var = ["subject_id"] + ["hematocrit_minus_5"] + ["hematocrit_90"] + ["hematocrit_120"] + ["map"] + ["clamp_map"] + ["total_protein"] + [v for v in meta.loc[meta["form_name"]
                                                == "outcomes", "field_name"]]
    out = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    out.replace(rep, np.nan, inplace=True)
    out.drop(["kidney_outcomes", "egfr", "metab_outcomes",
              "asl_outcomes", "bold_outcomes"],
             axis=1, inplace=True)
    # Kidney outcomes like GFR, etc. were collected with the clamp, not
    # necessarily the day of the MRI
    mri_cols = [c for c in out.columns if ("bold_" in c) or ("asl_" in c)]
    mri = out[["subject_id", "mri_date"] + mri_cols].copy()
    mri.rename({"mri_date": "date",
                "asl_left": "pcasl3d_left",
                "asl_right": "pcasl3d_right"}, axis=1, inplace=True)
    out = out[list(set(out.columns).difference(mri_cols))]
    rename = {"gfr": "gfr_raw_plasma", "gfr_bsa": "gfr_bsa_plasma",
              "rpf": "erpf_raw_plasma", "rpf_bsa": "erpf_bsa_plasma",
              "pah_bsa": "pah_clear_bsa", "abs_pah": "pah_clear_abs",
              "gfr_ecv": "gfr_ecv_percent", "gfr_standard": "gfr_ecv_std"}
    out.rename(rename, axis=1, inplace=True)
    # Calculate variables
    out_vars = ["gfr_raw_plasma", "erpf_raw_plasma", "total_protein", "map", "clamp_map", 
                "hematocrit_minus_5", "hematocrit_90", "hematocrit_120"]
    out[out_vars] = out[out_vars].apply(pd.to_numeric, errors='coerce')
    out["map"] = out[["map", "clamp_map"]].mean(axis=1)
    out["hematocrit_avg"] = out[["hematocrit_minus_5", "hematocrit_90", "hematocrit_120"]].mean(axis=1)
    out["erpf_raw_plasma_seconds"] = out["erpf_raw_plasma"]/60
    out["gfr_raw_plasma_seconds"] = out["gfr_raw_plasma"]/60
    # Filtration Fraction
    out["ff"] = out["gfr_raw_plasma"]/out["erpf_raw_plasma"] 
    # Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
    out["kfg"] = 0.1012
    # Filtration pressure across glomerular capillaries
    out["deltapf"] = (out["gfr_raw_plasma"]/60)/out["kfg"] 
    # Plasma protein mean concentration
    out["cm"] = (out["total_protein"]/out["ff"])*np.log(1/(1-out["ff"])) 
    # Pi G (Oncotic pressure)
    out["pg"] = 5*(out["cm"]-2)
    # Glomerular Pressure
    out["glomerular_pressure"] = out["pg"] + out["deltapf"] + 10
    # Renal Blood Flow
    out["rbf"] = (out["erpf_raw_plasma"]) / (1 - out["hematocrit_avg"]/100)
    out["rbf_seconds"] = (out["erpf_raw_plasma_seconds"]) / (1 - out["hematocrit_avg"]/100)
    # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    out["rvr"] = out["map"] / out["rbf"]
    # Efferent Arteriolar Resistance 
    out["re"] = (out["gfr_raw_plasma_seconds"]) / (out["kfg"] * (out["rbf_seconds"] - (out["gfr_raw_plasma_seconds"]))) * 1328
    # Afferent Arteriolar Resistance
    out["ra"] = ((out["map"] - out["glomerular_pressure"]) / out["rbf_seconds"]) * 1328    
    out.loc[~(out['ra'] > 0), 'ra']=np.nan    
    out.drop(["gfr_raw_plasma_seconds", "rbf_seconds", "gfr_raw_plasma_seconds", "erpf_raw_plasma_seconds", 
              "total_protein", "map", "clamp_map", "hematocrit_minus_5", "hematocrit_90", "hematocrit_120", "hematocrit_avg"],
             axis=1, inplace=True)
    out["date"] = clamp["date"]
    out["procedure"] = "clamp"
    out["visit"] = "baseline"
    mri["procedure"] = "bold_mri"
    mri["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=4, axis=0, inplace=True)
    phys.dropna(thresh=3, axis=0, inplace=True)
    screen.dropna(thresh=3, axis=0, inplace=True)
    mri.dropna(thresh=4, axis=0, inplace=True)
    dxa.dropna(thresh=4, axis=0, inplace=True)
    clamp.dropna(thresh=6, axis=0, inplace=True)
    out.dropna(thresh=4, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp], join='outer', ignore_index=True)
    df = pd.merge(df, out, how='outer')
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df["study"] = "CASPER"
    id_cols = ["subject_id", "co_enroll_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # Sort
    df.sort_values(["subject_id", "procedure"], inplace=True)
    # Rename subject identifier
    df.rename({"subject_id": "record_id"}, axis=1, inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data
    return df
