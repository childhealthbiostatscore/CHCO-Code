"""
This code is designed to pull data from the CROCODILE REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_crocodile():
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
        git_path = "/Users/shivaniramesh/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    else:
        sys.exit(f"Unknown user: please specify root path for this user. (Detected user: {user})")

    tokens = pd.read_csv(base_data_path + "/Data Harmonization/api_tokens.csv")        #"/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "CROCODILE", "Token"].iloc[0]
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

    dem_cols = ["record_id", "dob", "diabetes_dx_date",
                "group", "sex", "race", "ethnicity", "participation_status",  "mrn"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo["co_enroll_id"] = ""
    # Race columns combined into one
    demo = combine_checkboxes(demo, base_name="race", levels=[
        'American Indian/Alaska Native', 'Asian', 'Hawaiian/Pacific Islander', 'Black/African American', 'White', 'Other', 'Unknown'])
    # Same for ethnicity
    demo = combine_checkboxes(demo,
                              base_name="ethnicity",
                              levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
    # Relevel sex and group
    demo["sex"].replace({1: "Male", 2: "Female", 3: "Other",
                        "1": "Male", "2": "Female", "3": "Other"}, inplace=True)
    demo["group"].replace({1: "Type 1 Diabetes", 2: "Lean Control",
                           "1": "Type 1 Diabetes", "2": "Lean Control"}, inplace=True)
    demo["group_risk"] = np.where(
        demo.group.str.contains("lean", case=False), "Low", "High")
    demo["participation_status"].replace(
        {"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)
    
    dictionary.loc[dictionary['variable_name'].isin(demo.columns), 'form_name'] = 'demographics'


    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    med_list = {'diabetes_tx___1': "insulin_pump_timepoint",
                'diabetes_tx___2': "insulin_injections_timepoint",
                "diabetes_meds_other___1": "metformin_timepoint",
                "diabetes_meds_other___2": "tzd_timepoint",
                "diabetes_meds_other___3": "glp1_agonist_timepoint",
                "diabetes_meds_other___4": "sglti_timepoint",
                "diabetes_meds_other___5": "other_diabetes_med_timepoint",
                "htn_med___1": "ace_inhibitor",
                "htn_med___2": "angiotensin_receptor_blocker",
                "htn_med___3": "beta_blocker",
                "htn_med___4": "ca_channel_blocker",
                "htn_med___5": "diuretic",
                "htn_med___6": "statin",
                "pump_basal_rate": "pump_basal_rate",
                "cgm_yn": "cgm_yn",
                "mra_med": "mra",
                "fibrates_med": "fibrates",
                "meds_weight_type___1": "topiramate",
                "meds_weight_type___2": "phentermine",
                "uric_acid_med": "uric_acid_med"
                }
    dictionary.loc[dictionary['variable_name'] == 'ace_inhibitor', 'form_name'] = 'medical_history'

    og_names = list(med_list.keys())
    hx_list = [col for col in med.columns if col.startswith('hx_')]
    med = med[["record_id"] + og_names + hx_list]
    med.rename(med_list, axis=1, inplace=True)
    # RAASi
    med = med.assign(raasi_timepoint=np.maximum(pd.to_numeric(
        med["ace_inhibitor"]), pd.to_numeric(med["angiotensin_receptor_blocker"])))
    # Metformin
    med.rename({"diabetes_med_other___1": "metformin_timepoint"},
               axis=1, inplace=True)
    # Insulin
    med = med.assign(insulin_med_timepoint=np.maximum(pd.to_numeric(
        med["insulin_pump_timepoint"]), pd.to_numeric(med["insulin_injections_timepoint"])))
    # Replace 0/1 values with yes/no
    med.iloc[:, 1:] = med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 2: "No", "2": "No", 1: "Yes", "1": "Yes"})
    med["procedure"] = "medications"
    med["visit"] = "baseline"
    dictionary.loc[dictionary['variable_name'].isin(med.columns), 'form_name'] = 'medical_history'


    # --------------------------------------------------------------------------
    # EPIC Medications
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "epic_meds", "field_name"]]
    epic_med = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    epic_med.replace(rep, np.nan, inplace=True)
    # Replace 0/1 values with yes/no
    print("First 5 valid values in epic_ever_sglt2i_1 before replacement:")
    print(epic_med["epic_ever_sglt2i_1"].dropna().head())
    print("Data type of first epic_ever_sglt2i_1 value:")   
    print(epic_med["epic_ever_sglt2i_1"].dropna().head().apply(type))
    epic_med.iloc[:, 1:] = epic_med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 2: "No", "2": "No", 1: "Yes", "1": "Yes"})
    print("First 5 valid values in epic_ever_sglt2i_1 after replacement:")
    print(epic_med["epic_ever_sglt2i_1"].dropna().head())
    epic_med["procedure"] = "epic_medications"
    epic_med["visit"] = "baseline"

    
    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.drop(["phys_normal", "phys_abnormal"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_", "", regex=True)
    phys.rename({"sysbp": "sbp", "diasbp": "dbp"}, inplace=True, axis=1)
    phys["visit"] = "baseline"

    dictionary.loc[dictionary['variable_name'].isin(phys.columns), 'form_name'] = 'physical_exam'


    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    screen.replace(rep, np.nan, inplace=True)
    screen.drop(["prescreen_a1c", "prescreen_a1c_date",
                "screen_menstrual", "screen_upt"], axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"labs_|screen_", "", regex=True)
    screen.rename({"creat_s": "creatinine_s", "uacr": "acr_u", "a1c": "hba1c",
                   "creat_u": "creatinine_u", "hg": "hemoglobin"}, axis=1, inplace=True)
    screen["procedure"] = "screening"
    screen["visit"] = "baseline"
    # Assume medication review done at screening
    med["date"] = screen["date"]

    dictionary.loc[dictionary['variable_name'].isin(screen.columns), 'form_name'] = 'labs'


    # --------------------------------------------------------------------------
    # Labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_baseline_vitalslabs", "field_name"]]
    labs = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    labs.replace(rep, np.nan, inplace=True)
    labs.drop(["baseline_vitals", "visit_upt", "visit_weight", "visit_height",
               "visit_uptresult", "baseline_labs", "pilabs_yn", "pi_copeptin",
               "pi_renin", "pi_angiotensin2", "pi_osmo_s", "pi_osmo_u", "pi_lithium_s",
               "pi_lithium_u", "metabolomics_yn", "kim_yn", "pi_kim_ykl40", "pi_kim_ngal",
               "pi_kim_kim1", "pi_kim_il18", "pi_kim_tnfr1", "pi_kim_tnfr2"], axis=1, inplace=True)
    labs.columns = labs.columns.str.replace(
        r"visit_|bl_", "", regex=True)
    labs.rename({"uacr": "acr_u",
                "a1c": "hba1c",
                 "na_u": "sodium_u",
                 "na_s": "sodium_s",
                 "glucose_u": "urine_glucose_bl"}, axis=1, inplace=True)
    labs["procedure"] = "clamp"
    labs["visit"] = "baseline"

    dictionary.loc[dictionary['variable_name'].isin(labs.columns), 'form_name'] = 'labs'


    # --------------------------------------------------------------------------
    # BOLD/ASL MRI
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_boldasl_mri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    mri.replace(rep, np.nan, inplace=True)
    mri.columns = mri.columns.str.replace(
        r"mri_", "", regex=True)
    mri.rename({"volume_right": "right_kidney_volume_ml",
                "volume_left": "left_kidney_volume_ml"},
               axis=1, inplace=True)
    mri["date"] = labs["date"]
    mri["procedure"] = "bold_mri"
    mri["visit"] = "baseline"

    dictionary.loc[dictionary['variable_name'].isin(mri.columns), 'form_name'] = 'bold_mri'


    # --------------------------------------------------------------------------
    # DXA Scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_dxa_scan", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    dxa.replace(rep, np.nan, inplace=True)
    dxa.columns = dxa.columns.str.replace(
        r"dxa_|_percent", "", regex=True)
    dxa.rename({"bodyfat": "body_fat", "leanmass": "lean_mass",
                "trunkmass": "trunk_mass", "fatmass_kg": "fat_kg",
                "leanmass_kg": "lean_kg", "trunkmass_kg": "trunk_kg",
                "bmd": "bone_mineral_density"}, axis=1, inplace=True)
    dxa_cols = dxa.columns[4:].to_list()
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
    clamp.rename({"clamp_wt": "weight",
                  "clamp_ht": "height",
                  "cystatin_c": "cystatin_c_s",
                  "hct_210": "hematocrit_210",
                  "acr_baseline": "acr_u",
                  "acr_250": "acr_u_pm",
                  "clamp_d20": "d20_infusion"},
                 inplace=True, axis=1)
    clamp.columns = clamp.columns.str.replace(r"clamp_", "", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"insulin_minus", "insulin_minus_", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"ffa_minus", "ffa_minus_", regex=True)
    clamp.columns = clamp.columns.str.replace(r"bg_", "glucose_", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"glucose_minus", "glucose_minus_", regex=True)
    clamp["procedure"] = "clamp"
    clamp["visit"] = "baseline"
    clamp["insulin_sensitivity_method"] = "hyperinsulinemic_euglycemic_clamp"
    
    num_vars = ["d20_infusion", "weight"]
    clamp[num_vars] = clamp[num_vars].apply(
        pd.to_numeric, errors='coerce')
    
    clamp["gir_190"] = (clamp["d20_infusion"] * 190 / 60) / clamp["weight"] # previously M-value
    clamp["gir_200"] = (clamp["d20_infusion"] * 200 / 60) / clamp["weight"]

    # FFA
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
        clamp[['insulin_minus_20', 'insulin_minus_10', 'insulin_0']].mean(
            axis=1)
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
    rename = {"gfr_raw": "gfr_raw_plasma_urine", "gfr_bsa": "gfr_bsa_plasma_urine",
              "erpf_raw": "erpf_raw_plasma_urine", "erpf": "erpf_bsa_plasma_urine",
              "gfr_15mgmin": "gfr_raw_plasma", "gfrbsa": "gfr_bsa_plasma",
              "erpf_pah_85": "erpf_raw_plasma", "erpfbsa": "erpf_bsa_plasma",
              "pah_bsa": "pah_bsa_plasma_urine", "pahbsa": "pah_clear_bsa",
              "pahcl_12_8mgmin": "pah_clear_abs"}
    rct.rename(rename, axis=1, inplace=True)

    # Calculate variables
    rct_vars = ["gfr_raw_plasma", "erpf_raw_plasma",
                "bl_tot_protein", "visit_map", "phys_map", "hct_210"]
    rct[rct_vars] = rct[rct_vars].apply(pd.to_numeric, errors='coerce')
    rct["map"] = rct[["visit_map", "phys_map"]].mean(axis=1)
    rct["erpf_raw_plasma_seconds"] = rct["erpf_raw_plasma"] / 60
    rct["gfr_raw_plasma_seconds"] = rct["gfr_raw_plasma"] / 60
    # Filtration Fraction
    rct["ff"] = rct["gfr_raw_plasma"] / rct["erpf_raw_plasma"]
    # Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
    rct["kfg"] = np.select(
        [rct["group"].eq("1"), rct["group"].eq("2")], [0.1012, 0.1733])
    # Filtration pressure across glomerular capillaries
    rct["deltapf"] = (rct["gfr_raw_plasma"] / 60) / rct["kfg"]
    # Plasma protein mean concentration
    rct["cm"] = (rct["bl_tot_protein"] / rct["ff"]) * \
        np.log(1 / (1 - rct["ff"]))
    # Pi G (Oncotic pressure)
    rct["pg"] = 5 * (rct["cm"] - 2)
    # Glomerular Pressure
    rct["glomerular_pressure"] = rct["pg"] + rct["deltapf"] + 10
    # Renal Blood Flow
    rct["rbf"] = (rct["erpf_raw_plasma"]) / (1 - rct["hct_210"] / 100)
    rct["rbf_seconds"] = (rct["erpf_raw_plasma_seconds"]
                          ) / (1 - rct["hct_210"] / 100)
    # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    rct["rvr"] = rct["map"] / rct["rbf"]
    # Efferent Arteriolar Resistance
    rct["re"] = (rct["gfr_raw_plasma_seconds"]) / (rct["kfg"] *
                                                   (rct["rbf_seconds"] - (rct["gfr_raw_plasma_seconds"]))) * 1328
    # Afferent Arteriolar Resistance
    rct["ra"] = ((rct["map"] - rct["glomerular_pressure"]) /
                 rct["rbf_seconds"]) * 1328
    rct.loc[~(rct['ra'] > 0), 'ra'] = np.nan
    # Reduce rct dataset
    rct = rct[["record_id", "ff", "kfg", "deltapf", "cm", "pg",
               "glomerular_pressure", "rbf", "rvr", "ra", "re", 
               "pah_raw", "pah_sd", "pah_cv"] + list(rename.values())]
    rct["procedure"] = "renal_clearance_testing"
    rct["visit"] = "baseline"

    dictionary.loc[dictionary['variable_name'].isin(rct.columns), 'form_name'] = 'renal_clearance'


    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "optional_kidney_biopsy_56ba", "field_name"]]
    var = var + ["gloms", "gloms_gs", "ifta", "vessels_other", "fia",
                 "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins",
                 "glom_volume_con", "mes_matrix_area",
                 "mes_index", "mes_volume_weibel", "mes_volume_wiggins",
                 "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima",
                 "art_media", "pod_nuc_density", "pod_cell_volume",
                 "gbm_thick_artmean", "gbm_thick_harmmean",
                 "cortex_total_area", "cortex_analyzed_area", 
                 "cortex_percentage", "pt_total_number", 
                 "pt_total_area", "pt_epithelium_area", "pt_lumen_area", "pt_nuclear_count", 
                 "pt_total_nuc_area", "pt_density", "pt_avg_area", 
                 "pt_epithelium_avg_area", "pt_lumen_avg_area", 
                 "fractional_pt_total_area", "fractional_pt_epithelium_area", 
                 "fractional_pt_lumen_area", "pt_nuc_density_number_tubule", 
                 "pt_nuc_area", "pt_nuc_density_number_cortex", "pt_nuc_density_area_cortex"]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    biopsy.replace(rep, np.nan, inplace=True)
    biopsy.drop([col for col in biopsy.columns if '_yn' in col] +
                [col for col in biopsy.columns if 'procedure_' in col] +
                [col for col in biopsy.columns if '_image_' in col],
                axis=1, inplace=True)
    biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"labs_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"vitals_", "", regex=True)
    biopsy.rename({"hg": "hemoglobin"}, inplace=True, axis=1)
    biopsy["procedure"] = "kidney_biopsy"
    biopsy["visit"] = "baseline"
    epic_med["date"] = biopsy["date"]

    dictionary.loc[dictionary['variable_name'].isin(biopsy.columns), 'form_name'] = 'kidney_biopsy'


    # --------------------------------------------------------------------------
    # PET scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "optional_pet_scan", "field_name"]]
    pet = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    pet.replace(rep, np.nan, inplace=True)
    pet.drop(["petcon_yn"], axis=1, inplace=True)
    pet.columns = pet.columns.str.replace(r"pet_", "", regex=True)
    pet["procedure"] = "pet_scan"
    pet["visit"] = "baseline"

    dictionary.loc[dictionary['variable_name'].isin(pet.columns), 'form_name'] = 'pet_scan'


    
    # --------------------------------------------------------------------------
    # Liver PET scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "liver_pet_scan", "field_name"]]
    liver_pet = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    liver_pet.replace(rep, np.nan, inplace=True)
    liver_pet["procedure"] = "liver_pet_scan"
    liver_pet["visit"] = "baseline"

    dictionary.loc[dictionary['variable_name'].isin(liver_pet.columns), 'form_name'] = 'liver_pet_scan'

    
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
    
    dictionary.loc[dictionary['variable_name'].isin(voxelwise.columns), 'form_name'] = 'voxelwise'

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
    brain["date"] = screen["date"]
    
    # --------------------------------------------------------------------------
    # Pavel Labs
    # --------------------------------------------------------------------------
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "pavel_labs", "field_name"]]
    var = var + ["pavel_wbc", "pavel_neutrophils", "pavel_lymphocyte", "pavel_monocyte","pavel_rbc", "pavel_hgb", "pavel_hct", "pavel_plt", "pavel_mpv"]
    pavel = pd.DataFrame(proj.export_records(fields=var))
    pavel.replace(rep, np.nan, inplace=True)
    pavel["procedure"] = "labs"
    pavel["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Metabolomics (Blood and Tissue)
    # --------------------------------------------------------------------------
    
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"].isin(["metabolomics", "metabolomics_blood_raw"]), "field_name"]]
    metabolomics_blood = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    metabolomics_blood.replace(rep, np.nan, inplace=True)
    metabolomics_blood["procedure"] = "plasma_metab"
    metabolomics_blood["visit"] = "baseline"
    metabolomics_blood["date"] = screen["date"]
    
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "metabolomics_tissue", "field_name"]]
    metabolomics_tissue = pd.DataFrame(proj.export_records(fields=var))
    
    tissue_met = [c for c in metabolomics_tissue.columns if "_tissue" in c]
    metabolomics_tissue[tissue_met] = metabolomics_tissue[tissue_met].apply(
        pd.to_numeric, errors='coerce')
    metabolomics_tissue = metabolomics_tissue.groupby(by=["record_id"]).agg("mean", numeric_only = True).reset_index()
    # Replace missing values
    metabolomics_tissue.replace(rep, np.nan, inplace=True)
    metabolomics_tissue["procedure"] = "metabolomics_tissue"
    metabolomics_tissue["visit"] = "baseline"
    metabolomics_tissue["date"] = screen["date"]

    
    # --------------------------------------------------------------------------
    # Astrazeneca urine metabolomics
    # --------------------------------------------------------------------------
    
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "astrazeneca_urine_metabolomics", "field_name"]]
    az_u_metab = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    az_u_metab.replace(rep, np.nan, inplace=True)
    az_u_metab["procedure"] = "az_u_metab"
    az_u_metab["visit"] = "baseline"
    az_u_metab["date"] = labs["date"]
    # --------------------------------------------------------------------------
    # Lipidomics
    # --------------------------------------------------------------------------
    
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                  == "lipidomics", "field_name"]]
    lip = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    lip.replace(rep, np.nan, inplace=True)
    lip["procedure"] = "lipidomics"
    lip["visit"] = "baseline"
    dictionary.loc[dictionary['variable_name'].isin(lip.columns), 'form_name'] = 'lipidomics'

    
    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=5, axis=0, inplace=True)
    epic_med.dropna(thresh=5, axis=0, inplace=True)
    phys.dropna(thresh=4, axis=0, inplace=True)
    screen.dropna(thresh=4, axis=0, inplace=True)
    labs.dropna(thresh=4, axis=0, inplace=True)
    mri.dropna(thresh=4, axis=0, inplace=True)
    dxa.dropna(thresh=4, axis=0, inplace=True)
    clamp.dropna(thresh=6, axis=0, inplace=True)
    rct.dropna(thresh=4, axis=0, inplace=True)
    biopsy.dropna(thresh=4, axis=0, inplace=True)
    pet.dropna(thresh=4, axis=0, inplace=True)
    liver_pet.dropna(thresh=5, axis=0, inplace=True)
    brain.dropna(thresh=2, axis=0, inplace=True)
    voxelwise.dropna(thresh=4, axis=0, inplace=True)
    metabolomics_blood.dropna(thresh=4, axis=0, inplace=True)
    metabolomics_tissue.dropna(thresh=2, axis=0, inplace=True)
    az_u_metab.dropna(thresh=5, axis=0, inplace=True)
    lip.dropna(thresh=5, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------
    
    # Procedure = clamp
    clamp_merge = pd.merge(clamp, labs, how="outer")
    clamp_merge = pd.merge(clamp_merge, rct, how="outer")
    # Everything else
    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True) 
    df = pd.concat([df, epic_med], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp_merge], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    pet = pd.merge(pet, voxelwise, how = 'outer')
    df = pd.concat([df, pet], join='outer', ignore_index=True)
    df = pd.concat([df, liver_pet], join='outer', ignore_index=True)
    df = pd.concat([df, brain], join='outer', ignore_index=True)
    df = pd.concat([df, pavel], join='outer', ignore_index=True)

    df = pd.concat([df, metabolomics_blood], join='outer', ignore_index=True)
    df = pd.concat([df, metabolomics_tissue], join='outer', ignore_index=True)
    df = pd.concat([df, az_u_metab], join='outer', ignore_index=True)
    df = pd.concat([df, lip], join='outer', ignore_index=True)
    df = pd.merge(df, demo, on='record_id', how="outer")
    df = df.loc[:, ~df.columns.str.startswith('redcap_')]
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df["study"] = "CROCODILE"
    id_cols = ["record_id", "co_enroll_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # Sort
    df.sort_values(["record_id", "date", "procedure"], inplace=True)
    # Rename IDs
    df["record_id"] = ["CRC-" + str(i).zfill(2) for i in df["record_id"]]
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Print final data
    tocsv_path = base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv"
    dictionary.to_csv(tocsv_path, index=False)
    return df
