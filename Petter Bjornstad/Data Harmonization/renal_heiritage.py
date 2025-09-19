"""
This code is designed to pull data from the RENAL HEIRitage REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_renal_heiritage():
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

    tokens = pd.read_csv(base_data_path + "/Data Harmonization/api_tokens.csv")
        #"/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "Renal-HEIRitage", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["record_id", "dob", "group_rh2", "sex", "race", "ethnicity", "sglt2i", 
                "participation_status", "rh_id", "diabetes_dx_date", "mrn"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    demo = demo.loc[demo["redcap_event_name"].str.startswith('screen', na=False)].copy()
    demo.drop(["redcap_event_name"], inplace=True, axis=1)
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo.rename({"sglt2i": "sglt2i_ever",
                  "rh_id": "co_enroll_id",
                  "group_rh2": "group"},
                inplace=True, axis=1)
    dem_cols[2] = "group"
    dem_cols[6] = "sglt2i_ever"
    dem_cols[8] = "co_enroll_id"
    # Race columns combined into one
    demo = combine_checkboxes(demo, base_name="race", levels=[
        "American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
    # Same for ethnicity
    demo = combine_checkboxes(demo,
                              base_name="ethnicity",
                              levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
    # Relevel sex and group and participation status
    demo["sex"].replace({1: "Female", 2: "Male", 3: "Other",
                        "1": "Female", "2": "Male", "3": "Other"}, inplace=True)
    demo["group"].replace({1: "Type 2 Diabetes", 2: "Obese Control", 3: "Lean Control",
                           "1": "Type 2 Diabetes", "2": "Obese Control",
                           "3": "Lean Control"}, inplace=True)
    demo["group_risk"] = np.where(demo.group.str.contains("lean", case=False), "Low", "High")
    demo["sglt2i_ever"].replace({1: "Yes", 0: "No", "1": "Yes", "0": "No", np.nan: "No"},
                       inplace=True)
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)

    # --------------------------------------------------------------------------
    # Medical History
    # --------------------------------------------------------------------------

    var = ["record_id"] + ["mra_med"] + ["screen_a1c"] + ["insulin_inj"] + ["htn_med_type"] +  \
    ["sglt2i_med"] + ["diabetes_med_other"] + ["addl_hld_meds"] + ["meds_weight_type"] + \
    ["uric_acid_med"] + ["diabetes_med"]
    med = pd.DataFrame(proj.export_records(fields=var)) 
    med = med.loc[med["redcap_event_name"].str.startswith('screen', na=False)]
    med.drop(["redcap_event_name"], inplace=True, axis=1)
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    med.rename({"screen_a1c": "hba1c",
                "sglt2i_med": "sglti_timepoint",
                "diabetes_med_other___3": "glp1_agonist_timepoint",
                "htn_med_type___5": "diuretic",
                "htn_med_type___4": "ca_channel_blocker",
                "htn_med_type___3": "beta_blocker",
                "addl_hld_meds___1": "statin",
                "addl_hld_meds___2": "fibrates",
                "meds_weight_type___1": "topiramate",
                "meds_weight_type___2": "phentermine",
                "uric_acid_med": "uric_acid_med"},
                inplace=True, axis=1)
    med['insulin_med_timepoint'] = med.apply(lambda row: "1" if row['insulin_inj'] == "1" or row['diabetes_med___2'] == "1" else 0, axis=1)
    med['raasi_timepoint'] = med.apply(lambda row: "1" if row['htn_med_type___1'] == "1" or row['htn_med_type___2'] == "1" else 0, axis=1)
    # Replace 1 with yes 0 and 2 with no
    meds = ["sglti_timepoint", "glp1_agonist_timepoint", "diuretic", "ca_channel_blocker",\
    "beta_blocker", "statin", "fibrates", "topiramate", "phentermine", "uric_acid_med", "insulin_med_timepoint", "raasi_timepoint", "mra_med"]
    med[meds] = med[meds].applymap(lambda x: "Yes" if x == "1" or x is True else ("No" if pd.notna(x) and x != "" else x))
    med["procedure"] = "screening"
    med["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # EPIC Medications
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "epic_meds", "field_name"]]
    epic_med = pd.DataFrame(proj.export_records(fields=var))
    epic_med = epic_med.loc[epic_med["redcap_event_name"]. str.startswith('annual_blood_and_u_arm_1', na=False)]
    epic_med.drop(["redcap_event_name"], axis=1, inplace=True)
    epic_med.drop(["redcap_repeat_instance"], axis=1, inplace=True)
    epic_med.drop(["redcap_repeat_instrument"], axis=1, inplace=True)
    
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
                                                == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys = phys.loc[phys["redcap_event_name"].str.startswith('screen', na=False)]
    phys.drop(["redcap_event_name"], inplace=True, axis=1)
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.drop(["male_activity_factor", "fem_activity_factor", "schofield_male",
               "schofield_female", "phys_norm", "screen_bmi_percentile", "phys_age", "insulin_inj"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sysbp": "sbp", "diasbp": "dbp",
                "waist_circumference": "waistcm", 
                "hip_circumference": "hipcm"}, inplace=True, axis=1)
    phys["visit"] = "baseline"
    med["date"] = phys["date"]

    # --------------------------------------------------------------------------
    # Annual Labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "annual_labs", "field_name"]]
    annual_labs = pd.DataFrame(proj.export_records(fields=var))
    annual_labs = annual_labs.loc[annual_labs["redcap_event_name"].str.startswith('annual', na=False)]
    annual_labs.drop(["redcap_event_name"], inplace=True, axis=1)
    # Replace missing values
    annual_labs.replace(rep, np.nan, inplace=True)
    # Format
    annual_labs["procedure"] = "labs"
    annual_labs["visit"] = "baseline"
    annual_labs.drop(["annual_labs", "redcap_repeat_instance", "redcap_repeat_instrument", "study_visit_annual"], axis=1, inplace=True)
    annual_labs.columns = annual_labs.columns.str.replace(
        r"an_|annual_lab_", "", regex=True)
    annual_labs.rename({"uacr": "acr_u"},
                 inplace=True, axis=1)

    # --------------------------------------------------------------------------
    # Renal Clearance Testing
    # --------------------------------------------------------------------------

    var = ["record_id"] + ["group_rh2"] + ["phys_map"] + [v for v in meta.loc[meta["form_name"] == 
    "study_visit_renal_clearance_testing", "field_name"]] +[v for v in meta.loc[meta["form_name"] == 
    "renal_clearance_baseline_labs", "field_name"]]
    rct = pd.DataFrame(proj.export_records(fields=var))
    numeric_cols = [col for col in rct.columns if col not in ['record_id', 'redcap_event_name'] 
                    and not col.startswith('tm_') and not col.endswith(('_date', '_com', '_start'))]
    rct[numeric_cols] = rct[numeric_cols].apply(pd.to_numeric, errors = 'ignore')
    rct = rct.groupby('record_id', as_index=False).max()
    rct.drop(["redcap_event_name"], inplace=True, axis=1)
    # Replace missing values
    rct.replace(rep, np.nan, inplace=True)
    rename = {"gfr_raw": "gfr_raw_plasma_urine", "gfr_bsa": "gfr_bsa_plasma_urine",
              "erpf_raw": "erpf_raw_plasma_urine", "erpf": "erpf_bsa_plasma_urine",
              "gfr_15mgmin": "gfr_raw_plasma", "gfrbsa": "gfr_bsa_plasma",
              "erpf_pah_85": "erpf_raw_plasma", "erpfbsa": "erpf_bsa_plasma",
              "phys_map": "map", "pah_bsa": "pah_bsa_plasma_urine", "pahbsa": "pah_clear_bsa",
              "pahcl_12_8mgmin": "pah_clear_abs", 
              "rc_lab_date": "date"}
    rct.rename(rename, axis=1, inplace=True)
    rct.columns = rct.columns.str.replace(
        r"bl_", "", regex=True)
    # Calculate variables
    rct_vars = ["gfr_raw_plasma", "erpf_raw_plasma", "tot_protein", "map", "hct"]
    rct[rct_vars] = rct[rct_vars].apply(pd.to_numeric, errors='coerce')
    rct["erpf_raw_plasma_seconds"] = rct["erpf_raw_plasma"]/60
    rct["gfr_raw_plasma_seconds"] = rct["gfr_raw_plasma"]/60
    # Filtration Fraction
    rct["ff"] = rct["gfr_raw_plasma"]/rct["erpf_raw_plasma"] 
    # Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
    rct["kfg"] = np.select([rct["group_rh2"].eq(1), rct["group_rh2"].eq(2), rct["group_rh2"].eq(3)], [0.1012, 0.1733, 0.1733]) 
    # Filtration pressure across glomerular capillaries
    rct["deltapf"] = (rct["gfr_raw_plasma"]/60)/rct["kfg"] 
    # Plasma protein mean concentration
    rct["cm"] = (rct["tot_protein"]/rct["ff"])*np.log(1/(1-rct["ff"])) 
    # Pi G (Oncotic pressure)
    rct["pg"] = 5*(rct["cm"]-2)
    # Glomerular Pressure
    rct["glomerular_pressure"] = rct["pg"] + rct["deltapf"] + 10
    # Renal Blood Flow
    rct["rbf"] = (rct["erpf_raw_plasma"]) / (1 - rct["hct"]/100)
    rct["rbf_seconds"] = (rct["erpf_raw_plasma_seconds"]) / (1 - rct["hct"]/100)
    # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    rct["rvr"] = rct["map"] / rct["rbf"]
    # Efferent Arteriolar Resistance 
    rct["re"] = (rct["gfr_raw_plasma_seconds"]) / (rct["kfg"] * (rct["rbf_seconds"] - (rct["gfr_raw_plasma_seconds"]))) * 1328
    # Afferent Arteriolar Resistance
    rct["ra"] = ((rct["map"] - rct["glomerular_pressure"]) / rct["rbf_seconds"]) * 1328  
    rct.loc[~(rct['ra'] > 0), 'ra']=np.nan    
    # Reduce rct dataset
    rct.drop(["rbf_seconds", "erpf_raw_plasma_seconds", "redcap_repeat_instrument", "map", "group_rh2"], axis=1, inplace=True)
    rct["procedure"] = "renal_clearance_testing"
    rct["visit"] = "baseline"
    rct = rct[rct['date'].notna()]
    
    # --------------------------------------------------------------------------
    # Dextran
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "dextran_data", "field_name"]]
    dextran = pd.DataFrame(proj.export_records(fields=var))
    dextran = dextran.loc[dextran["redcap_event_name"].str.startswith('renal', na=False)]
    dextran.drop(["redcap_event_name"], inplace=True, axis=1)
    # Replace missing values
    dextran.replace(rep, np.nan, inplace=True)

    # --------------------------------------------------------------------------
    # Outcomes
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "study_visit_boldasl_mri", "field_name"]]
    out = pd.DataFrame(proj.export_records(fields=var))
    out = out.loc[out["redcap_event_name"].str.startswith("kidney_and", na=False)].copy()
    out.drop(["redcap_event_name", "adc_outcomes"], axis=1, inplace=True)
    # Replace missing values
    out.replace(rep, np.nan, inplace=True)
    # Kidney outcomes like GFR, etc. were collected with the clamp, not
    # necessarily the day of the MRI
    bold_mri_cols = [c for c in out.columns if ("bold_" in c) or ("asl_" in c)]
    bold_mri = out[["record_id"] + bold_mri_cols].copy()
    out = out[list(set(out.columns).difference(bold_mri_cols))]
    rename = {"volume_left": "left_kidney_volume_ml",
              "volume_right": "right_kidney_volume_ml",
              "mri_lab_date": "date"}
    out.rename(rename, axis=1, inplace=True)
    out["procedure"] = "clamp"
    out["visit"] = "baseline"
    out = out[out['date'].notna()]
    bold_mri["procedure"] = "bold_mri"
    bold_mri["visit"] = "baseline"
    bold_mri["date"] = out["date"]
  
    # --------------------------------------------------------------------------
    # PET scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "pet_scan", "field_name"]]
    pet = pd.DataFrame(proj.export_records(fields=var))
    pet = pet.loc[pet["redcap_event_name"].str.startswith('pet', na=False)]
    pet.drop(["redcap_event_name"], axis=1, inplace=True)
    # Replace missing values
    pet.replace(rep, np.nan, inplace=True)
    pet.drop(["petcom_yn"], axis=1, inplace=True)
    pet.columns = pet.columns.str.replace(r"pet_", "", regex=True)
    pet["procedure"] = "pet_scan"
    pet["visit"] = "baseline"
    
    # --------------------------------------------------------------------------
    # Liver PET scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "liver_pet_scan", "field_name"]]
    liver_pet = pd.DataFrame(proj.export_records(fields=var))
    liver_pet = liver_pet.loc[liver_pet["redcap_event_name"].str.startswith('pet_scan', na=False)]
    liver_pet.drop(["redcap_event_name"], axis=1, inplace=True)
    # Replace missing values
    liver_pet.replace(rep, np.nan, inplace=True)
    liver_pet["procedure"] = "liver_pet_scan"
    liver_pet["visit"] = "baseline"
    
    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "kidney_biopsy", "field_name"]]
    var = var + ["cortex_dx_other", "gloms", "gloms_gs", "ifta", "vessels_other", "fia",
                 "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins",
                 "glom_volume_con", "mes_matrix_area",
                 "mes_index", "mes_volume_weibel", "mes_volume_wiggins",
                 "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima",
                 "art_media", "pod_nuc_density", "pod_cell_volume"]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    biopsy = biopsy.loc[biopsy["redcap_event_name"].str.startswith("kidney_bio", na=False)]
    # Replace missing values
    biopsy.replace(rep, np.nan, inplace=True)
    biopsy.drop([col for col in biopsy.columns if '_yn' in col] +
                [col for col in biopsy.columns if 'procedure_' in col] +
                ["core_diagnostic", "core_hypo_cryo", "core_oct", "core_rna"],
                axis=1, inplace=True)
    biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"labs_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"vitals_", "", regex=True)
    biopsy.rename({"hg": "hemoglobin"}, axis=1, inplace=True)
    biopsy["procedure"] = "kidney_biopsy"
    biopsy["visit"] = "baseline"
    
    # --------------------------------------------------------------------------
    # Voxelwise
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "voxelwise", "field_name"]]
    voxelwise = pd.DataFrame(proj.export_records(fields=var))
    voxelwise.drop(["redcap_event_name", "redcap_repeat_instrument", "redcap_repeat_instance"], axis=1, inplace=True)
    # Replace missing values
    voxelwise.replace(rep, np.nan, inplace=True)
    voxelwise["procedure"] = "pet_scan"
    voxelwise["visit"] = "baseline"
    
    # --------------------------------------------------------------------------
    # Neurocognitive Tracking
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "neurocognitive_tracking", "field_name"]]
    neuro = pd.DataFrame(proj.export_records(fields=var))
    neuro = neuro.loc[neuro["redcap_event_name"].str.startswith('renal', na=False)]
    neuro.drop(["redcap_event_name"], axis=1, inplace=True)
    # Replace missing values
    neuro.replace(rep, np.nan, inplace=True)
    neuro["procedure"] = "neurocognitive_tracking"
    neuro["visit"] = "baseline"
    neuro = neuro.loc[neuro.filter(regex='^procedr_').sum(axis=1) != 0]
    
    # --------------------------------------------------------------------------
    # Astrazeneca urine metabolomics
    # --------------------------------------------------------------------------
    
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "astrazeneca_urine_metabolomics", "field_name"]]
    az_u_metab = pd.DataFrame(proj.export_records(fields=var))
    az_u_metab.drop(["redcap_event_name"], axis=1, inplace=True)
    az_u_metab.drop(["redcap_repeat_instance"], axis=1, inplace=True)
    az_u_metab.drop(["redcap_repeat_instrument"], axis=1, inplace=True)
    # Replace missing values
    az_u_metab.replace(rep, np.nan, inplace=True)
    az_u_metab["procedure"] = "az_u_metab"
    az_u_metab["visit"] = "baseline"
    az_u_metab["date"] = annual_labs["date"]

    # --------------------------------------------------------------------------
    # Plasma metabolomics
    # --------------------------------------------------------------------------
    
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "metabolomics_blood_raw", "field_name"]]
    plasma_metab = pd.DataFrame(proj.export_records(fields=var))
    plasma_metab.drop(["redcap_event_name"], axis=1, inplace=True)
    plasma_metab.drop(["redcap_repeat_instance"], axis=1, inplace=True)
    plasma_metab.drop(["redcap_repeat_instrument"], axis=1, inplace=True)

    # Replace missing values
    plasma_metab.replace(rep, np.nan, inplace=True)
    plasma_metab["procedure"] = "plasma_metab"
    plasma_metab["date"] = annual_labs["date"]
    plasma_metab["visit"] = "baseline"
    
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
    epic_med.dropna(thresh=5, axis=0, inplace=True)
    phys.dropna(thresh=4, axis=0, inplace=True)
    annual_labs.dropna(thresh=4, axis=0, inplace=True)
    rct.dropna(thresh=4, axis=0, inplace=True)
    dextran.dropna(thresh=4, axis=0, inplace=True)
    out.dropna(thresh=4, axis=0, inplace=True)
    bold_mri.dropna(thresh=4, axis=0, inplace=True)
    pet.dropna(thresh=6, axis=0, inplace=True)
    liver_pet.dropna(thresh=5, axis=0, inplace=True)
    brain.dropna(thresh=2, axis=0, inplace=True)
    biopsy.dropna(thresh=7, axis=0, inplace=True)
    neuro.dropna(thresh=4, axis=0, inplace=True)
    voxelwise.dropna(thresh=4, axis=0, inplace=True)
    az_u_metab.dropna(thresh=5, axis=0, inplace=True)
    plasma_metab.dropna(thresh=10, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([phys, annual_labs], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, epic_med], join='outer', ignore_index=True)
    df = pd.concat([df, rct], join='outer', ignore_index=True)
    df = pd.concat([df, dextran], join='outer', ignore_index=True)
    df = pd.merge(df, out, how='outer')
    df = pd.concat([df, bold_mri], join='outer', ignore_index=True)
    pet = pd.merge(pet, voxelwise, how = 'outer')
    df = pd.concat([df, pet], join='outer', ignore_index=True)
    df = pd.concat([df, liver_pet], join='outer', ignore_index=True)
    df = pd.concat([df, brain], join='outer', ignore_index=True)
    df = pd.concat([df, neuro], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.concat([df, az_u_metab], join='outer', ignore_index=True)
    df = pd.concat([df, plasma_metab], join='outer', ignore_index=True)
    df = pd.merge(df, demo, on='record_id', how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df["study"] = "RENAL-HEIRitage"
    id_cols = ["record_id", "study"] + \
        dem_cols[2:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    df = df[id_cols + other_cols]
    # Sort
    df.sort_values(["record_id", "procedure"], inplace=True)
    # Rename subject identifier
    df.rename({"record_id": "record_id"}, axis=1, inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data
    return df
