"""
This code is designed to pull data from the IMPROVE REDCap project and output data in a long format with one row per study procedure per visit.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_improve():
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
    token = tokens.loc[tokens["Study"] == "IMPROVE", "Token"].iloc[0]
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

    dem_cols = ["subject_id", "co_enroll_id", "dob", "diagnosis", 
                "gender", "race", "ethnicity", "sglt2i", "participation_status", "mr_number"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols,
                                            events=["screening_arm_1"]))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo["group"] = "Type 2 Diabetes"
    demo.drop(redcap_cols, axis=1, inplace=True)
    demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date", "mr_number": "mrn"},
                inplace=True, axis=1)
    dem_cols[3] = "diabetes_dx_date"
    dem_cols[4] = "sex"
    dem_cols[7] = "sglt2i_ever"
    dem_cols[9] = "mrn"
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
    demo["sex"].replace({1: "Male", 0: "Female", 2: "Other",
                        "1": "Male", "0": "Female", "2": "Other"}, inplace=True)
    demo["sglt2i"].replace({1: "Yes", 0: "No", "1": "Yes", "0": "No", np.NaN: "No"},
                           inplace=True)
    demo["group_risk"] = np.where(demo.group.str.contains("lean", case=False), "Low", "High")
    demo.rename({"sglt2i": "sglt2i_ever"}, axis=1, inplace=True)
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit", "med_date", "diabetes_med",
           "diabetes_med_other", "htn_med_type", "addl_hld_meds", "hypertension"]
    med = pd.DataFrame(proj.export_records(fields=var))
    med.drop(redcap_cols, axis=1, inplace=True)
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    # SGLT2i (diabetes_med_other___4), RAASi (htn_med_type___1, htn_med_type___2), Metformin (diabetes_med_other___1)
    med = med[["subject_id", "study_visit", "med_date", "diabetes_med_other___3", "htn_med_type___1",
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
    # EPIC Medications
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "epic_meds", "field_name"]]
    epic_med = pd.DataFrame(proj.export_records(fields=var))
    epic_med.drop(redcap_cols, axis=1, inplace=True)
    # Replace missing values
    epic_med.replace(rep, np.nan, inplace=True)
    # Replace 0/1 values with yes/no
    epic_med.iloc[:, 1:] = epic_med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 2: "No", "2": "No", 1: "Yes", "1": "Yes"})
    epic_med["procedure"] = "epic_medications"
    
    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(
        fields=var, events=["screening_arm_1"]))
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys.drop(redcap_cols + ["phys_norm", "phys_no", "breast_tanner",
                             "testicular_volume", "lmp", "screen_bmi_percentile",
                             "activity_factor_male", "activity_factor_female",
                             "schofield_male", "schofield_female"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sys_bp": "sbp", "dys_bp": "dbp",
                 "waist_circumference": "waistcm",
                 "hip_circumference": "hipcm"}, inplace=True, axis=1)
    phys["procedure"] = "physical_exam"

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var,
                                              events=["screening_arm_1"]))
    screen.replace(rep, np.nan, inplace=True)  # Replace missing values
    screen.drop(redcap_cols + ['a1c_pre', 'a1c_pre_date', "screen_pregnant"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|_of_screen", "", regex=True)
    screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
                   "urine_cre": "creatinine_u", "urine_mab": "microalbumin_u"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"

    # --------------------------------------------------------------------------
    # Accelerometry
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "accelerometry", "field_name"]]
    accel = pd.DataFrame(proj.export_records(fields=var))
    accel.replace(rep, np.nan, inplace=True)  # Replace missing values
    accel = accel.loc[accel["acc_wear_percent"] != ""]
    accel.drop(redcap_cols + ["study_visit_accel"], axis=1, inplace=True)
    accel.columns = accel.columns.str.replace(
        r"acc_|accel_", "", regex=True)
    accel["procedure"] = "accelerometry"

    # --------------------------------------------------------------------------
    # Cardio/Abdominal MRI
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "cardioabdominal_mri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri.replace(rep, np.nan, inplace=True)  # Replace missing values
    mri.drop(redcap_cols + ["mri_cardio", "mri_abdo",
                            "mri_aortic", "study_visit_mri"],
             axis=1, inplace=True)
    mri.columns = mri.columns.str.replace(
        r"mri_|visit_", "", regex=True)
    mri["procedure"] = "cardio_abdominal_mri"

    # --------------------------------------------------------------------------
    # MMTT + Metabolic Cart
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "mmtt_metabolic_cart", "field_name"]]
    mmtt = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    mmtt.replace(rep, np.nan, inplace=True)
    # Drop unnecessary columns
    mmtt.drop(redcap_cols + ["study_visit_mttt", "mmtt_vitals", "mmtt_pregnant",
                             "mmtt_lmp", "mmtt_brmr", "mmtt_60rmr", "mmtt_base_labs", "mmtt_ffa_labs", "mmtt_insulin",
                             "mmtt_glp1", "mmtt_cpep", "mmtt_yy", "mmtt_glucagon",
                             "mmtt_gluc"],
              axis=1, inplace=True)
    mmtt.columns = mmtt.columns.str.replace(
        r"mmtt_", "", regex=True)
    mmtt.columns = mmtt.columns.str.replace(
        r"_neg_", "_minus_", regex=True)
    mmtt.columns = mmtt.columns.str.replace(r"bg_", "glucose_", regex=True)
    mmtt.columns = mmtt.columns.str.replace(r"cpep_", "cpeptide_", regex=True)
    mmtt.rename({"wt": "weight", "ht": "height", "waist": "waistcm",
                "hip": "hipcm", "hr": "pulse", "sys_bp": "sbp",
                 "dia_bp": "dbp", "hba1c_base": "hba1c", "na_base": "sodium_base",
                 "totprot_base": "tot_protein", "glucose_base": "glucose_bl",
                 "gluc_base": "glucose_bl_cmp"},
                inplace=True, axis=1)
    mmtt["procedure"] = "mmtt"
    # FFA
    # See /home/timvigers/Work/CHCO/Petter Bjornstad/IHD/Background/Renal Heir Equations.docx
    ffa = [c for c in mmtt.columns if "ffa_" in c]
    mmtt[ffa] = mmtt[ffa].apply(
        pd.to_numeric, errors='coerce')
    mmtt["baseline_ffa"] = mmtt[['ffa_minus_10', 'ffa_0']].mean(axis=1)
    mmtt["steady_state_ffa"] = mmtt['ffa_120']
    mmtt["ffa_suppression"] = (
        (mmtt["baseline_ffa"] - mmtt["steady_state_ffa"]) / mmtt["baseline_ffa"]) * 100
    mmtt["ffa_method"] = "mmtt"

    # --------------------------------------------------------------------------
    # DXA
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "body_composition_dxa_bod_pod", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    dxa.replace(rep, np.nan, inplace=True)
    dxa.drop(redcap_cols + ["study_visit_bodycomp", "dxa_complete",
                            "bodpod_complete"],
             axis=1, inplace=True)
    dxa.columns = dxa.columns.str.replace(
        r"dxa_", "dexa_", regex=True)
    dxa.columns = dxa.columns.str.replace(
        r"bp_", "bod_pod_", regex=True)
    dxa.rename({"dexa_bmd": "dexa_bone_mineral_density",
                "bod_pod_fat_mass": "bod_pod_fat_kg",
                "bodcomp_date": "date"}, axis=1, inplace=True)
    dxa["procedure"] = "dxa"

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    clamp.replace(rep, np.nan, inplace=True)
    # Format
    clamp.drop(redcap_cols + ["study_visit_clamp", "baseline", "fasting_labs",
                              "bg_labs", "ns_bolus", "urine_labs", "cpep_lab",
                              "iohexol_bolus", "pah_bolus", "hct_lab", "insulin_labs"],
               axis=1, inplace=True)
    clamp.columns = clamp.columns.str.replace(
        r"clamp_", "", regex=True)
    clamp.rename({"cystatin_c": "cystatin_c_s", 
                  "urine_mab": "microalbumin_u",
                  "serum_creatinine": "creatinine_s", 
                  "acr_baseline": "acr_u",
                  "urine_mab_baseline": "microalbumin_u",
                  "urine_cre_baseline": "creatinine_u",
                  "serum_sodium": "sodium_s", 
                  "urine_sodium": "sodium_u",
                  "pls": "pulse", 
                  "acr_250": "acr_u_pm", 
                  "total_protein": "tot_protein",
                  "glucose_bl": "urine_glucose_bl"
                  }, inplace=True, axis=1)
    clamp["procedure"] = "clamp"
    clamp["insulin_sensitivity_method"] = "hyperglycemic_clamp"
    # M
    num_vars = ["d20_infusion", "weight"]
    clamp[num_vars] = clamp[num_vars].apply(pd.to_numeric, errors='coerce')
    clamp["gir_190"] = (clamp["d20_infusion"] * 190 / 60) / clamp["weight"] # previously M-value
    clamp["gir_200"] = (clamp["d20_infusion"] * 200 / 60) / clamp["weight"]
    # No FFA
    # Insulin
    ins = ['insulin_minus_10', 'insulin_minus_5', 'insulin_2', 'insulin_4',
           'insulin_6', 'insulin_8', 'insulin_10', 'insulin_120', 'insulin_220',
           'insulin_230', 'insulin_240', 'insulin_245', 'insulin_249',
           'insulin_250', 'insulin_252', 'insulin_253', 'insulin_254',
           'insulin_255']
    clamp[ins] = clamp[ins].apply(
        pd.to_numeric, errors='coerce')
    clamp["steady_state_insulin"] = clamp[['insulin_220',
                                           'insulin_230', 'insulin_240', 'insulin_245', 'insulin_249',
                                           'insulin_250', 'insulin_252', 'insulin_253', 'insulin_254',
                                           'insulin_255']].mean(axis=1) * 6
    # C peptide
    cpep = [c for c in clamp.columns if "cpeptide" in c]
    clamp[cpep] = clamp[cpep].apply(
        pd.to_numeric, errors='coerce')
    clamp["steady_state_cpeptide"] = clamp[['cpeptide_220', 'cpeptide_230',
                                            'cpeptide_240', 'cpeptide_245',
                                            'cpeptide_249', 'cpeptide_250',
                                            'cpeptide_252', 'cpeptide_253',
                                            'cpeptide_254', 'cpeptide_255']].mean(axis=1)
    # ACPRg
    clamp["acprg"] = clamp[['cpeptide_2', 'cpeptide_4',
                            'cpeptide_6', 'cpeptide_8', 'cpeptide_10']].mean(axis=1) - clamp[['cpeptide_minus_10', 'cpeptide_minus_5']].mean(axis=1)
    # Negative ACPRg to 0.01
    clamp.loc[clamp["acprg"] < 0, "acprg"] = 0.01
    #AIRg
    clamp["airg"] = clamp[['insulin_2', 'insulin_4',
                           'insulin_6', 'insulin_8', 'insulin_10']].mean(axis=1) * 6 - clamp[['insulin_minus_10', 'insulin_minus_5']].mean(axis=1) * 6
    # Negative AIRg to 0.01
    clamp.loc[clamp["airg"] < 0, "airg"] = 0.01
    # DI
    clamp["di"] = (clamp["gir_190"] /
                   clamp["steady_state_insulin"]) * clamp["airg"]
    # Accelerometry done 1 day before clamp
    accel["date"] = [pd.to_datetime(d) - timedelta(days=1)
                     for d in clamp["date"]]
    accel["date"] = accel["date"].dt.strftime('%m/%d/%Y')
    # Hematocrit average
    hematocrit_vars = ["hematocrit_90", "hematocrit_120"]
    clamp[hematocrit_vars] = clamp[hematocrit_vars].apply(
        pd.to_numeric, errors='coerce')
    clamp["hct"] = clamp[["hematocrit_90", "hematocrit_120"]].mean(axis=1)

    # --------------------------------------------------------------------------
    # Outcomes
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + ["hematocrit_90"] + ["hematocrit_120"] + ["map"] + ["clamp_map"] + ["total_protein"] + [v for v in meta.loc[meta["form_name"]
                                                               == "outcomes", "field_name"]]
    out = pd.DataFrame(proj.export_records(fields=var))
    out.replace(rep, np.nan, inplace=True)  # Replace missing values
    out.drop(redcap_cols + ["kidney_outcomes", "egfr", "metab_outcomes", "study_visit_outcomes",
                            "asl_outcomes", "bold_outcomes", "adc_outcomes"],
             axis=1, inplace=True)
    # Kidney outcomes like GFR, etc. were collected with the clamp, not
    # necessarily the day of the MRI
    bold_mri_cols = [c for c in out.columns if ("bold_" in c) or ("asl_" in c)]
    bold_mri = out[["subject_id", "study_visit",
                    "mri_date"] + bold_mri_cols].copy()
    bold_mri.rename({"mri_date": "date",
                     "asl_left": "pcasl3d_left",
                     "asl_right": "pcasl3d_right"}, axis=1, inplace=True)
    out = out[list(set(out.columns).difference(bold_mri_cols))]
    rename = {"gfr": "gfr_raw_plasma", "gfr_bsa": "gfr_bsa_plasma",
              "rpf": "erpf_raw_plasma", "erpf_bsa": "erpf_bsa_plasma",
              "pah_bsa": "pah_clear_bsa", "abs_pah": "pah_clear_abs"}
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
    out["rbf"] = (out["erpf_raw_plasma"]) / (1 - out["hct"]/100)
    out["rbf_seconds"] = (out["erpf_raw_plasma_seconds"]) / (1 - out["hct"]/100)
    # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    out["rvr"] = out["map"] / out["rbf"]
    # Efferent Arteriolar Resistance 
    out["re"] = (out["gfr_raw_plasma_seconds"]) / (out["kfg"] * (out["rbf_seconds"] - (out["gfr_raw_plasma_seconds"]))) * 1328
    # Afferent Arteriolar Resistance
    out["ra"] = ((out["map"] - out["glomerular_pressure"]) / out["rbf_seconds"]) * 1328    
    out.loc[~(out['ra'] > 0), 'ra']=np.nan    
    out.drop(["gfr_raw_plasma_seconds", "rbf_seconds", "erpf_raw_plasma_seconds",
              "hematocrit_90" , "hematocrit_120" , "map" , "clamp_map" , "total_protein" , "hct"],
             axis=1, inplace=True)
    out.rename({"mri_date": "date"}, axis=1, inplace=True)
    out["procedure"] = "clamp"
    bold_mri["procedure"] = "bold_mri"

    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "kidney_biopsy", "field_name"]]
    var = var + ["gloms", "gloms_gs", "ifta", "vessels_other", "fia",
                 "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins",
                 "glom_volume_con", "mes_matrix_area",
                 "mes_index", "mes_volume_weibel", "mes_volume_wiggins",
                 "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima",
                 "art_media", "pod_nuc_density", "pod_cell_volume"]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    biopsy.replace(rep, np.nan, inplace=True)
    biopsy.drop(redcap_cols + [col for col in biopsy.columns if '_yn' in col] +
                [col for col in biopsy.columns if 'procedure_' in col] +
                ["core_diagnostic", "core_hypo_cryo", "core_oct", "core_rna"],
                axis=1, inplace=True)
    biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"labs_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"vitals_", "", regex=True)
    biopsy.rename({"hg": "hemoglobin"}, inplace=True, axis=1)
    biopsy["procedure"] = "kidney_biopsy"

    # --------------------------------------------------------------------------
    # Astrazeneca urine metabolomics
    # --------------------------------------------------------------------------
    
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "astrazeneca_urine_metabolomics", "field_name"]]
    az_u_metab = pd.DataFrame(proj.export_records(fields=var))
    az_u_metab.drop(redcap_cols, axis=1, inplace=True)
    # Replace missing values
    az_u_metab.replace(rep, np.nan, inplace=True)
    az_u_metab["procedure"] = "az_u_metab"
    az_u_metab["date"] = clamp["date"]
    
    # --------------------------------------------------------------------------
    # Plasma metabolomics
    # --------------------------------------------------------------------------
    
    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "metabolomics_blood_raw", "field_name"]]
    plasma_metab = pd.DataFrame(proj.export_records(fields=var))
    plasma_metab.drop(redcap_cols, axis=1, inplace=True)
    # Replace missing values
    plasma_metab.replace(rep, np.nan, inplace=True)
    plasma_metab["procedure"] = "plasma_metab"
    plasma_metab["date"] = clamp["date"]
    
    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=9, axis=0, inplace=True)
    epic_med.dropna(thresh=5, axis=0, inplace=True)
    phys.dropna(thresh=3, axis=0, inplace=True)
    screen.dropna(thresh=3, axis=0, inplace=True)
    accel.dropna(thresh=6, axis=0, inplace=True)
    mri.dropna(thresh=4, axis=0, inplace=True)
    mmtt.dropna(thresh=5, axis=0, inplace=True)
    dxa.dropna(thresh=4, axis=0, inplace=True)
    clamp.dropna(thresh=5, axis=0, inplace=True)
    out.dropna(thresh=10, axis=0, inplace=True)
    bold_mri.dropna(thresh=4, axis=0, inplace=True)
    biopsy.dropna(thresh=12, axis=0, inplace=True)
    az_u_metab.dropna(thresh=10, axis=0, inplace=True)
    plasma_metab.dropna(thresh=10, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([accel, med], join='outer', ignore_index=True)
    df = pd.concat([df, epic_med], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.concat([df, mmtt], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp], join='outer', ignore_index=True)
    df = pd.concat([df, bold_mri], join='outer', ignore_index=True)
    df = pd.merge(df, out, how='outer')
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.concat([df, screen], join='outer', ignore_index=True)
    df = pd.concat([df, phys], join='outer', ignore_index=True)
    df = pd.concat([df, az_u_metab], join='outer', ignore_index=True)
    df = pd.concat([df, plasma_metab], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df.rename({"study_visit": "visit"}, axis=1, inplace=True)
    df["study"] = "IMPROVE"
    id_cols = ["subject_id", "study"] + \
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
    df["subject_id"] = df["subject_id"].str.replace(r"2D-", "_", regex=True)
    # Sort
    df.sort_values(["subject_id", "visit", "procedure"], inplace=True)
    # Rename subject identifier
    df.rename({"subject_id": "record_id"}, axis=1, inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data
    return df
