"""
This code is designed to pull data from the RENAL HEIR REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_renal_heir():
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

    tokens = pd.read_csv(base_data_path + "/Data Harmonization/api_tokens.csv")        #"/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "Renal-HEIR", "Token"].iloc[0]
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

    dem_cols = ["subject_id", "co_enroll_id", "dob", "diagnosis", 
                "group", "gender", "race", "ethnicity", "sglt2i", "participation_status","mr_number"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date", "sglt2i": "sglt2i_ever", "mr_number": "mrn"},
                inplace=True, axis=1)
    dem_cols[3] = "diabetes_dx_date"
    dem_cols[5] = "sex"
    dem_cols[8] = "sglt2i_ever"
    dem_cols[10] = "mrn"
    # Race columns combined into one
    demo = combine_checkboxes(demo, base_name="race", levels=[
        "American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
    # Same for ethnicity
    demo = combine_checkboxes(demo,
                              base_name="ethnicity",
                              levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
    # Relevel sex and group and participation status
    demo["sex"].replace({1: "Male", 0: "Female", 2: "Other",
                        "1": "Male", "0": "Female", "2": "Other"}, inplace=True)
    demo["group"].replace({2: "Type 2 Diabetes", 3: "Obese Control",
                           4: "Lean Control",
                           "2": "Type 2 Diabetes", "3": "Obese Control",
                           "4": "Lean Control"}, inplace=True)
    demo["group_risk"] = np.where(demo.group.str.contains("lean", case=False), "Low", "High")
    demo["sglt2i_ever"].replace({1: "Yes", 0: "No", "1": "Yes", "0": "No", np.nan: "No"},
                       inplace=True)
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                             == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    # SGLT2i
    med = med[["subject_id", "mra_med", "diabetes_med___1", "diabetes_med_other___1", "addl_hld_meds___2", "diabetes_med_other___2",
              "diabetes_med_other___3", "diabetes_med___3", "addl_hld_meds___1", "htn_med_type___1",
              "htn_med_type___2", "htn_med_type___3", "htn_med_type___4", "htn_med_type___5"]]

    # RAASi
    med = med.assign(raasi_timepoint=np.maximum(pd.to_numeric(
        med["htn_med_type___1"]), pd.to_numeric(med["htn_med_type___2"])))
    med_list = {"raasi_timepoint": "raasi_timepoint",
                "diabetes_med_other___3": "sglti_timepoint",
                "diabetes_med_other___1": "tzd_timepoint",
                "addl_hld_meds___2": "fibrates_timepoint",
                "diabetes_med___1": "metformin_timepoint",
                "diabetes_med___3": "insulin_med_timepoint",
                "addl_hld_meds___1": "statin",
                "diabetes_med_other___2": "glp1_agonist_timepoint",
                "htn_med_type___1": "ace_inhibitor",
                "htn_med_type___2": "angiotensin_receptor_blocker",
                "htn_med_type___3": "beta_blocker",
                "htn_med_type___4": "ca_channel_blocker",
                "htn_med_type___5": "diuretic"}
    og_names = list(med_list.keys())
    med.rename(med_list, axis=1, inplace=True)
    med.iloc[:, 1:2] = med.iloc[:, 1:2].replace(
        {2: "No", "2": "No", 1: "Yes", "1": "Yes"})
    med.iloc[:, 2:] = med.iloc[:, 2:].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"})

    med["procedure"] = "medications"
    med["visit"] = "baseline"
    
    # --------------------------------------------------------------------------
    # EPIC Medications
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "epic_meds", "field_name"]]
    epic_med = pd.DataFrame(proj.export_records(fields=var))
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

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.drop(["male_activity_factor", "fem_activity_factor", "schofield_male",
               "schofield_female", "phys_norm", "phys_no", "lmp", "screen_bmi_percentile"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sys_bp": "sbp", "dys_bp": "dbp",
                "waist_circumference": "waistcm", "hip_circumference": "hipcm"}, inplace=True, axis=1)
    phys["visit"] = "baseline"
    med["date"] = phys["date"]

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    screen.replace(rep, np.nan, inplace=True)
    screen.drop(["a1c_pre", "a1c_pre_date", "screen_pregnant"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"_of_screen|screen_", "", regex=True)
    screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
                   "urine_mab": "microalbumin_u", "urine_cre": "creatinine_u"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"
    screen["visit"] = "baseline"

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
    dxa_cols = dxa.columns[4:].to_list()
    dxa.rename(dict(zip(dxa_cols, ["dexa_" + d for d in dxa_cols])),
               axis=1, inplace=True)
               
    dxa["procedure"] = "dxa"
    dxa["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    clamp.replace(rep, np.nan, inplace=True)
    # PB suspects error in labs for blood samples at T=-10 for RH-72-T, omit from analysis (05/11/23)
    clamp.loc[clamp['subject_id']=='RH-72-T', 'cpeptide_minus_10'] = np.nan
    clamp.loc[clamp['subject_id']=='RH-72-T', 'insulin_minus_10'] = np.nan
    clamp.loc[clamp['subject_id']=='RH-72-T', 'ffa_minus_10'] = np.nan
    clamp.loc[clamp['subject_id']=='RH-72-T', 'glucose_minus_10'] = np.nan
    # Format
    clamp.drop(["baseline", "fasting_labs", "urine_labs", "hct_lab",
                "bg_labs", "ffa_lab", "cpep_lab", "insulin_labs"], axis=1, inplace=True)
    clamp.columns = clamp.columns.str.replace(
        r"clamp_", "", regex=True)
    clamp.rename({"serum_creatinine": "creatinine_s",
                  "serum_sodium": "sodium_s",
                  "urine_sodium": "sodium_u",
                  "cystatin_c": "cystatin_c_s",
                  "urine_mab_baseline": "microalbumin_u",
                  "urine_cre_baseline": "creatinine_u",
                  "pls": "pulse", 
                  "acr_baseline": "acr_u", 
                  "acr_250": "acr_u_pm",
                  "total_protein": "tot_protein",
                  "glucose_bl": "urine_glucose_bl"},
                 inplace=True, axis=1)
    clamp.columns = clamp.columns.str.replace(r"clamp_", "", regex=True)
    clamp["procedure"] = "clamp"
    clamp["visit"] = "baseline"
    clamp["insulin_sensitivity_method"] = "hyperglycemic_clamp"
    # M
    num_vars = ["d20_infusion", "weight"]
    clamp[num_vars] = clamp[num_vars].apply(
        pd.to_numeric, errors='coerce')
    clamp["gir_190"] = (clamp["d20_infusion"] * 190 / 60) / clamp["weight"] # previously M-value
    clamp["gir_200"] = (clamp["d20_infusion"] * 200 / 60) / clamp["weight"]
    # FFA
    # See /home/timvigers/Work/CHCO/Petter Bjornstad/IHD/Background/Renal Heir Equations.docx
    ffa = [c for c in clamp.columns if "ffa_" in c]
    clamp[ffa] = clamp[ffa].apply(
        pd.to_numeric, errors='coerce')
    clamp["baseline_ffa"] = \
        clamp[['ffa_minus_10', 'ffa_minus_5']].mean(axis=1)
    clamp["steady_state_ffa"] = \
        clamp[['ffa_220', 'ffa_230', 'ffa_240', 'ffa_250']].mean(axis=1)
    clamp["ffa_suppression"] = (
        (clamp["baseline_ffa"] - clamp["steady_state_ffa"]) / clamp["baseline_ffa"]) * 100
    clamp["ffa_method"] = "hyperglycemic_clamp"
    # Insulin
    ins = ['insulin_minus_10', 'insulin_minus_5', 'insulin_2', 'insulin_4',
           'insulin_6', 'insulin_8', 'insulin_10', 'insulin_120',
           'insulin_220', 'insulin_230', 'insulin_240', 'insulin_250']
    clamp[ins] = clamp[ins].apply(
        pd.to_numeric, errors='coerce')
    clamp["steady_state_insulin"] = clamp[['insulin_220',
                                           'insulin_230', 'insulin_240', 'insulin_250']].mean(axis=1) * 6
    # C peptide
    cpep = [c for c in clamp.columns if "cpeptide" in c]
    clamp[cpep] = clamp[cpep].apply(
        pd.to_numeric, errors='coerce')
    clamp["steady_state_cpeptide"] = clamp[['cpeptide_220',
                                            'cpeptide_230', 'cpeptide_240', 'cpeptide_250']].mean(axis=1)
    # ACPRg
    clamp["acprg"] = clamp[['cpeptide_2', 'cpeptide_4',
                            'cpeptide_6', 'cpeptide_8', 'cpeptide_10']].mean(axis=1) - clamp[['cpeptide_minus_10', 'cpeptide_minus_5']].mean(axis=1)
    # Negative ACPRg to 0.01
    clamp.loc[clamp["acprg"] < 0, "acprg"] = 0.01
    # AIRg
    clamp["airg"] = clamp[['insulin_2', 'insulin_4',
                           'insulin_6', 'insulin_8', 'insulin_10']].mean(axis=1) * 6 - clamp[['insulin_minus_10', 'insulin_minus_5']].mean(axis=1) * 6
    # Negative AIRg to 0.01
    clamp.loc[clamp["airg"] < 0, "airg"] = 0.01    
    # DI
    clamp["di"] = \
        (clamp["gir_190"] / clamp["steady_state_insulin"]) * clamp["airg"]

    # Hematocrit average
    hematocrit_vars = ["hematocrit_90", "hematocrit_120"]
    clamp[hematocrit_vars] = clamp[hematocrit_vars].apply(
        pd.to_numeric, errors='coerce')
    clamp["hct"] = clamp[["hematocrit_90", "hematocrit_120"]].mean(axis=1)

    # --------------------------------------------------------------------------
    # Outcomes
    # --------------------------------------------------------------------------

    var = ["subject_id"] + ["group"] + ["hematocrit_90"] + ["hematocrit_120"] + ["map"] + ["clamp_map"] + ["total_protein"] + [v for v in meta.loc[meta["form_name"]
                                                == "outcomes", "field_name"]]
    out = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    out.replace(rep, np.nan, inplace=True)
    out.drop(["kidney_outcomes", "egfr", "metab_outcomes", "asl_outcomes", "adc_outcomes", "tkv_outcomes"],
             axis=1, inplace=True)
    # Kidney outcomes like GFR, etc. were collected with the clamp, not
    # necessarily the day of the MRI
    bold_mri_cols = [c for c in out.columns if ("bold_" in c) or ("asl_" in c)]
    bold_mri = out[["subject_id", "mri_date"] + bold_mri_cols].copy()
    bold_mri.rename({"mri_date": "date",
                     "asl_left": "pcasl3d_left",
                     "asl_right": "pcasl3d_right"}, axis=1, inplace=True)
    bold_mri.drop(["bold_outcomes"], axis=1, inplace=True)
    out = out[list(set(out.columns).difference(bold_mri_cols))]
    rename = {"gfr": "gfr_raw_plasma", "gfr_bsa": "gfr_bsa_plasma",
              "rpf": "erpf_raw_plasma", "rpf_bsa": "erpf_bsa_plasma",
              "volume_right": "right_kidney_volume_ml",
              "volume_left": "left_kidney_volume_ml"}
    out.rename(rename, axis=1, inplace=True)
    
    # Calculate variables
    out_vars = ["gfr_raw_plasma", "erpf_raw_plasma", "total_protein", "map", "clamp_map", "hematocrit_90", "hematocrit_120"]
    out[out_vars] = out[out_vars].apply(pd.to_numeric, errors='coerce')
    out["map"] = out[["map", "clamp_map"]].mean(axis=1)
    out["hct"] = out[["hematocrit_90", "hematocrit_120"]].mean(axis=1)
    out["erpf_raw_plasma_seconds"] = out["erpf_raw_plasma"]/60
    out["gfr_raw_plasma_seconds"] = out["gfr_raw_plasma"]/60
    # Filtration Fraction
    out["ff"] = out["gfr_raw_plasma"]/out["erpf_raw_plasma"] 
    # Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
    out["kfg"] = np.select([out["group"].eq("2"), out["group"].eq("3"), out["group"].eq("4")], [0.1012, 0.1733, 0.1733]) 
    # Filtration pressure across glomerular capillaries
    out["deltapf"] = (out["gfr_raw_plasma"]/60)/out["kfg"] 
    # Plasma protein mean concentration
    out["cm"] = (out["total_protein"]/out["ff"])*np.log(1/(1-out["ff"])) 
    # Pi G (Oncotic pressure)
    out["pg"] = 5*(out["cm"]-2)
    # Glomerular Pressure
    out["glomerular_pressure"] = out["pg"] + out["deltapf"] + 10
    # Renal Blood Flow
    out["rbf"] = (out["erpf_raw_plasma"]) / (1 - out["hct"]/100)
    out["rbf_seconds"] = (out["erpf_raw_plasma_seconds"]) / (1 - out["hct"]/100)
    # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    out["rvr"] = out["map"] / out["rbf"]
    # Efferent Arteriolar Resistance 
    out["re"] = (out["gfr_raw_plasma_seconds"]) / (out["kfg"] * (out["rbf_seconds"] - (out["gfr_raw_plasma_seconds"]))) * 1328
    # Afferent Arteriolar Resistance
    out["ra"] = ((out["map"] - out["glomerular_pressure"]) / out["rbf_seconds"]) * 1328    
    out.loc[~(out['ra'] > 0), 'ra']=np.nan    
    out.drop(["gfr_raw_plasma_seconds", "rbf_seconds", "gfr_raw_plasma_seconds", "erpf_raw_plasma_seconds", 
            "total_protein", "map", "clamp_map", "hematocrit_90", "hematocrit_120", "hct", "group"],
             axis=1, inplace=True)
    out.rename({"mri_date": "date"}, axis=1, inplace=True)
    out["procedure"] = "clamp"
    out["visit"] = "baseline"
    bold_mri["procedure"] = "bold_mri"
    bold_mri["visit"] = "baseline"
    outcome_vars = out.columns.difference(["subject_id", "procedure", "visit", "date"])
    dictionary.loc[dictionary['variable_name'].isin(outcome_vars), 'form_name'] = 'outcomes'


    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------

    var = ["subject_id", ] + [v for v in meta.loc[meta["form_name"]
                                                  == "kidney_biopsy", "field_name"]]
    var = var + ["cortex_dx_other", "gloms", "gloms_gs", "ifta", "vessels_other", "fia",
                 "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins",
                 "glom_volume_con", "mes_matrix_area",
                 "mes_index", "mes_volume_weibel", "mes_volume_wiggins",
                 "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima",
                 "art_media", "pod_nuc_density", "pod_cell_volume",
                 "gbm_thick_artmean", "gbm_thick_harmmean"]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
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
    # Astrazeneca urine metabolomics
    # --------------------------------------------------------------------------
    
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "astrazeneca_urine_metabolomics", "field_name"]]
    az_u_metab = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    az_u_metab.replace(rep, np.nan, inplace=True)
    az_u_metab["procedure"] = "az_u_metab"
    az_u_metab["visit"] = "baseline"
    az_u_metab["date"] = phys["date"]

    # --------------------------------------------------------------------------
    # Plasma metabolomics
    # --------------------------------------------------------------------------
    
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "metabolomics_blood_raw", "field_name"]]
    plasma_metab = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    plasma_metab.replace(rep, np.nan, inplace=True)
    plasma_metab["procedure"] = "plasma_metab"
    plasma_metab["visit"] = "baseline"
    plasma_metab["date"] = clamp["date"]
    
    # --------------------------------------------------------------------------
    # Brain biomarkers
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "biomarkers", "field_name"]]
    brain = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    brain.replace(rep, np.nan, inplace=True)
    brain["procedure"] = "brain_biomarkers"
    brain["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Pavel Labs
    # --------------------------------------------------------------------------
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "pavel_labs", "field_name"]]
    var = var + ["pavel_wbc", "pavel_neutrophil", "pavel_lymphocyte", "pavel_monocyte","pavel_rbc", "pavel_hgb", "pavel_hct", "pavel_plt", "pavel_mpv"]
    pavel = pd.DataFrame(proj.export_records(fields=var))
    pavel.replace(rep, np.nan, inplace=True)
    pavel["procedure"] = "labs"
    pavel["visit"] = "baseline"

    
    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=4, axis=0, inplace=True)
    epic_med.dropna(thresh=5, axis=0, inplace=True)
    phys.dropna(thresh=4, axis=0, inplace=True)
    screen.dropna(thresh=4, axis=0, inplace=True)
    dxa.dropna(thresh=4, axis=0, inplace=True)
    clamp.dropna(thresh=6, axis=0, inplace=True)
    out.dropna(thresh=4, axis=0, inplace=True)
    bold_mri.dropna(thresh=4, axis=0, inplace=True)
    brain.dropna(thresh=2, axis=0, inplace=True)
    biopsy.dropna(thresh=4, axis=0, inplace=True)
    az_u_metab.dropna(thresh=5, axis=0, inplace=True)
    plasma_metab.dropna(thresh=10, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, epic_med], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp], join='outer', ignore_index=True)
    df = pd.merge(df, out, how='outer')
    df = pd.concat([df, bold_mri], join='outer', ignore_index=True)
    df = pd.concat([df, brain], join='outer', ignore_index=True)
    df = pd.concat([df, pavel], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.concat([df, az_u_metab], join='outer', ignore_index=True)
    df = pd.concat([df, plasma_metab], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df["study"] = "RENAL-HEIR"
    id_cols = ["subject_id", "co_enroll_id", "study"] + \
        dem_cols[2:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    df = df[id_cols + other_cols]
    # Sort
    df.sort_values(["subject_id", "procedure"], inplace=True)
    # Rename subject identifier
    df.rename({"subject_id": "record_id"}, axis=1, inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data
    return df
