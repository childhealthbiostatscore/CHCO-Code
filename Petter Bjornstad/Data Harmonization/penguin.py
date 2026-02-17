"""
This code is designed to pull data from the PENGUIN REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_penguin():
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
    token = tokens.loc[tokens["Study"] == "PENGUIN", "Token"].iloc[0]
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

    dem_cols = ["record_id", "dob", "group", "sex", "race", "ethnicity", "participation_status", "mrn"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
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
                              levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
    # Relevel sex and group
    demo["sex"].replace({2: "Male", 1: "Female", 3: "Other",
                        "2": "Male", "1": "Female", "3": "Other"}, inplace=True)
    demo["diabetes_dx_date"] = np.nan
    demo["co_enroll_id"] = np.nan
    demo["group"] = "PKD"
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    med_list = {"htn_med___1": "ace_inhibitor",
                "htn_med___2": "angiotensin_receptor_blocker",
                "htn_med___3": "beta_blocker",
                "htn_med___4": "ca_channel_blocker",
                "htn_med___5": "diuretic",
                "htn_med___6": "statin"}
    og_names = list(med_list.keys())
    med = med[["record_id"] + og_names]
    med.rename(med_list, axis=1, inplace=True)
    # Insulin med (no one with diabetes/insulin)
    med["insulin_med_timepoint"] = 0
    # RAASi
    med = med.assign(raasi_timepoint=np.maximum(pd.to_numeric(
        med["ace_inhibitor"]), pd.to_numeric(med["angiotensin_receptor_blocker"])))
    # Replace 0/1 values with yes/no
    med.iloc[:, 1:] = med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"})
    med["procedure"] = "medications"
    med["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys.drop(["phys_age", "phys_normal", "phys_abnormal"],
              axis=1, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sysbp": "sbp", "diasbp": "dbp"}, inplace=True, axis=1)
    phys["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    screen.replace(rep, np.nan, inplace=True)
    screen.drop(["screen_creat_lab", "screen_upt", "screen_menstrual"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|labs_", "", regex=True)
    screen.rename({"creat_s": "creatinine_s"}, axis=1, inplace=True)
    screen["procedure"] = "screening"
    screen["visit"] = "baseline"
    med["date"] = phys["date"]

    # --------------------------------------------------------------------------
    # Baseline labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_baseline_vitalslabs", "field_name"]]
    labs = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    labs.replace(rep, np.nan, inplace=True)
    labs.drop(["baseline_vitals", "visit_upt", "visit_uptresult", "visit_weight", "visit_height", 
               "baseline_labs", "u24_labs", "pilabs_yn", "metabolomics_yn", "kim_yn"], axis=1, inplace=True)
    labs.columns = labs.columns.str.replace(
        r"visit_|bl_", "", regex=True)
    labs.rename({"a1c": "hba1c", 
                "uacr": "acr_u",
                "glucose_u": "urine_glucose_bl"
                }, axis=1, inplace=True)
    labs["procedure"] = "clamp"
    labs["visit"] = "baseline"

    dictionary.loc[dictionary['variable_name'].isin(labs.columns), 'form_name'] = 'labs'

    # --------------------------------------------------------------------------
    # DXA Scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_dxa_scan", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    dxa.replace(rep, np.nan, inplace=True)
    dxa.rename({"bodyfat_percent": "body_fat",
                "leanmass_percent": "lean_mass", "fatmass_kg": "fat_kg",
                "leanmass_kg": "lean_kg", "trunkmass_kg": "trunk_kg",
                "trunkmass_percent": "trunk_mass",
                "bmd": "bone_mineral_density"}, axis=True, inplace=True)
    dxa.columns = dxa.columns.str.replace(
        r"dxa_", "", regex=True)
    dxa_cols = dxa.columns[2:].to_list()
    dxa.rename(dict(zip(dxa_cols, ["dexa_" + d for d in dxa_cols])),
               axis=1, inplace=True)
    dxa["procedure"] = "dxa"
    dxa["visit"] = "baseline"
    dictionary.loc[dictionary['variable_name'].isin(dxa.columns), 'form_name'] = 'dxa'


    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_he_clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    clamp.replace(rep, np.nan, inplace=True)
    # Format
    clamp.drop(["clamp_yn", "clamp_ffa",
                "clamp_insulin", "hct_yn", "clamp_bg"], axis=1, inplace=True)
    clamp.rename({"clamp_wt": "weight", "clamp_ht": "height", "hct_210": "hematocrit_avg"},
                 inplace=True, axis=1)
    clamp.columns = clamp.columns.str.replace(r"clamp_", "", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"insulin_minus", "insulin_minus_", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"ffa_minus", "ffa_minus_", regex=True)
    clamp.columns = clamp.columns.str.replace(r"bg_", "glucose_", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"glucose_minus", "glucose_minus_", regex=True)
    clamp.rename({"d20": "d20_infusion"},
                 inplace=True, axis=1)  
    num_vars = ["d20_infusion", "weight"]
    clamp[num_vars] = clamp[num_vars].apply(
        pd.to_numeric, errors='coerce')
    
    clamp["gir_190"] = (clamp["d20_infusion"] * 190 / 60) / clamp["weight"] # previously M-value
    clamp["gir_200"] = (clamp["d20_infusion"] * 200 / 60) / clamp["weight"]
    
    clamp["procedure"] = "clamp"
    clamp["visit"] = "baseline"
    clamp["insulin_sensitivity_method"] = "hyperinsulinemic_euglycemic_clamp"
    # FFA
    # See /home/timvigers/Work/CHCO/Petter Bjornstad/IHD/Background/Renal Heir Equations.docx
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
        clamp[['insulin_minus_20', 'insulin_minus_10', 'insulin_0']].mean(axis=1)
    clamp["p1_steady_state_insulin"] = \
        clamp[['insulin_70', 'insulin_80', 'insulin_90']].mean(axis=1)
    clamp["p2_steady_state_insulin"] = \
        clamp[['insulin_250', 'insulin_260', 'insulin_270']].mean(axis=1)
    clamp["ffa_method"] = "hyperinsulinemic_euglycemic_clamp"

    dictionary.loc[dictionary['variable_name'].isin(clamp.columns), 'form_name'] = 'clamp'


    # --------------------------------------------------------------------------
    # Renal Clearance Testing
    # --------------------------------------------------------------------------

    var = ["record_id"] + ["group"] + ["bl_tot_protein"] + ["hct_210"] + ["visit_map"] + ["phys_map"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_renal_clearance_testing", "field_name"]]
    rct = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    rct.replace(rep, np.nan, inplace=True)
    rename = {"gfr": "gfr_raw_plasma_urine", "gfr_bsa": "gfr_bsa_plasma_urine",
              "erpf": "erpf_raw_plasma_urine", "erpf_bsa": "erpf_bsa_plasma_urine",
              "gfr_15mgmin": "gfr_raw_plasma", "gfrbsa": "gfr_bsa_plasma",
              "erpf_pah_85": "erpf_raw_plasma", "erpfbsa": "erpf_bsa_plasma",
              "pah_bsa": "pah_bsa_plasma_urine"}
    rct.rename(rename, axis=1, inplace=True)
    # Calculate variables
    rct_vars = ["gfr_raw_plasma", "erpf_raw_plasma", "bl_tot_protein", "visit_map", "phys_map", "hct_210"]
    rct[rct_vars] = rct[rct_vars].apply(pd.to_numeric, errors='coerce')
    rct["map"] = rct[["visit_map", "phys_map"]].mean(axis=1)
    rct["erpf_raw_plasma_seconds"] = rct["erpf_raw_plasma"]/60
    rct["gfr_raw_plasma_seconds"] = rct["gfr_raw_plasma"]/60
    # Filtration Fraction
    rct["ff"] = rct["gfr_raw_plasma"]/rct["erpf_raw_plasma"] 
    # Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
    rct["kfg"] = np.select([rct["group"].eq("1"), rct["group"].eq("2")], [0.1012, 0.1733]) 
    # Filtration pressure across glomerular capillaries
    rct["deltapf"] = (rct["gfr_raw_plasma"]/60)/rct["kfg"] 
    # Plasma protein mean concentration
    rct["cm"] = (rct["bl_tot_protein"]/rct["ff"])*np.log(1/(1-rct["ff"])) 
    # Pi G (Oncotic pressure)
    rct["pg"] = 5*(rct["cm"]-2)
    # Glomerular Pressure
    rct["glomerular_pressure"] = rct["pg"] + rct["deltapf"] + 10
    # Renal Blood Flow
    rct["rbf"] = (rct["erpf_raw_plasma"]) / (1 - rct["hct_210"]/100)
    rct["rbf_seconds"] = (rct["erpf_raw_plasma_seconds"]) / (1 - rct["hct_210"]/100)
    # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    rct["rvr"] = rct["map"] / rct["rbf"]
    # Efferent Arteriorlar Resistance 
    rct["re"] = (rct["gfr_raw_plasma_seconds"]) / (rct["kfg"] * (rct["rbf_seconds"] - (rct["gfr_raw_plasma_seconds"]))) * 1328
    # Afferent Arteriolar Resistance
    rct["ra"] = ((rct["map"] - rct["glomerular_pressure"]) / rct["rbf_seconds"]) * 1328    
    rct.loc[~(rct['ra'] > 0), 'ra']=np.nan    
    # Reduce rct dataset
    rct = rct[["record_id", "ff", "kfg", "deltapf", "cm", "pg", "glomerular_pressure", "rbf", "rvr", "ra", "re"] + list(rename.values())] 
    rct["procedure"] = "clamp"
    rct["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # PET scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "pet_scan", "field_name"]]
    pet = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    pet.replace(rep, np.nan, inplace=True)
    pet.drop(["petcom_yn"], axis=1, inplace=True)
    pet.columns = pet.columns.str.replace(r"pet_", "", regex=True)
    pet["procedure"] = "pet_scan"
    pet["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Voxelwise
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "voxelwise", "field_name"]]
    voxelwise = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    voxelwise.replace(rep, np.nan, inplace=True)
    voxelwise["procedure"] = "pet_scan"
    voxelwise["visit"] = "baseline"
    
    # --------------------------------------------------------------------------
    # fMRI
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "fmri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    mri.replace(rep, np.nan, inplace=True)
    mri["procedure"] = "mri"
    mri["visit"] = "baseline"
    
    # --------------------------------------------------------------------------
    # Brain biomarkers
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "biomarkers", "field_name"]]
    brain = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    brain.replace(rep, np.nan, inplace=True)
    brain["procedure"] = "brain_biomarkers"
    brain["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=4, axis=0, inplace=True)
    phys.dropna(thresh=4, axis=0, inplace=True)
    screen.dropna(thresh=4, axis=0, inplace=True)
    labs.dropna(thresh=4, axis=0, inplace=True)
    mri.dropna(thresh=4, axis=0, inplace=True)
    dxa.dropna(thresh=4, axis=0, inplace=True)
    clamp.dropna(thresh=6, axis=0, inplace=True)
    rct.dropna(thresh=4, axis=0, inplace=True)
    pet.dropna(thresh=4, axis=0, inplace=True)
    voxelwise.dropna(thresh=4, axis=0, inplace=True)
    brain.dropna(thresh=2, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------
    # Procedure = clamp
    clamp_merge = pd.merge(clamp, labs, how="outer")
    clamp_merge = pd.merge(clamp_merge, rct,  how="outer")
    # Everything else
    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp_merge], join='outer', ignore_index=True)
    pet = pd.merge(pet, voxelwise, how = 'outer')
    df = pd.concat([df, pet], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.concat([df, brain], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df["study"] = "PENGUIN"
    id_cols = ["record_id", "co_enroll_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # SORT
    df.sort_values(["record_id", "procedure"], inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data

    tocsv_path = base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv"
    dictionary.to_csv(tocsv_path, index=False) 
    return df
