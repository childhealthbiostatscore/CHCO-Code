"""
This code is designed to pull data from the PANTHER REDCap project and output data in a long format with one row per study procedure per visit.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_panther():
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
    token = tokens.loc[tokens["Study"] == "PANTHER", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999]
    rep = rep + [str(r) for r in rep] + [""]

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["record_id", "group", "dob", "t2d_date", 
                "sex", "race", "ethnicity", "sglt2i", "participation_status","mrn"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols,
                                            events=["screening_arm_1"]))
    # Replace missing values
    demo.replace(rep, np.nan, inplace=True)
    demo.drop(["redcap_event_name"], inplace=True, axis=1)
    demo.rename({"t2d_date": "diabetes_dx_date"},
                inplace=True, axis=1)
    dem_cols[3] = "diabetes_dx_date"
    dem_cols[7] = "sglt2i_ever"
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
    demo["sex"].replace({1: "Female", 2: "Male", 3: "Other",
                        "1": "Female", "2": "Male", "3": "Other"}, inplace=True)
    demo["sglt2i"].replace({1: "Yes", 0: "No", "1": "Yes", "0": "No", np.NaN: "No"},
                           inplace=True)
    demo["group"].replace({"1": "Type 2 Diabetes", "2": "Obese Control", "3": "Lean Control"}, inplace=True)
    demo["group_risk"] = np.where(demo.group.str.contains("lean", case=False), "Low", "High")
    demo.rename({"sglt2i": "sglt2i_ever"}, axis=1, inplace=True)
    demo["participation_status"].replace({"1": "Participated", "2": "Removed", "3": "Will Participate"}, inplace=True)

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["record_id", "screen_creat_s", "screen_a1c", "consent_date"]
    screen = pd.DataFrame(proj.export_records(fields=var, events=["screening_arm_1"]))
    # Replace missing values
    screen.replace(rep, np.nan, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|_of_screen|consent_", "", regex=True)
    screen.rename({"creat_s": "creatinine_s", "a1c": "hba1c"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"
    screen["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["record_id", "dkd_meds_dose", "t2d_meds_dose", "htn_med",\
    "hx_cv_medlist", "meds_weight_type", "glp1ra_med", "mra_med", \
    "fibrates_med", "uric_acid_med", "screen_birthweight"]
    med = pd.DataFrame(proj.export_records(fields=var, events=["screening_arm_1"]))
    # Replace missing values
    med.replace(rep, np.nan, inplace=True)
    # Visit ID
    med["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"}, 
        inplace=True)
    med = med.rename(columns={"redcap_event_name": "visit"})
    med.rename({"screen_birthweight": "birthweight"}, inplace=True, axis=1)
    # Meds by regex
    med["sglti_timepoint"]=med.apply(lambda x: x.str.contains('sgl|canag|dapagl|empag|floz|sglt2i_med', case=False).any(), axis=1)
    med["ace_inhibitor"]=med.apply(lambda x: x.str.contains('lisino|liso', case=False).any(), axis=1)
    med["raasi_timepoint"]=np.where(np.logical_or.reduce((med.ace_inhibitor==True, med.htn_med___1=="1", med.htn_med___2=="1")), True, False)
    med["metformin_timepoint"]=med.apply(lambda x: x.str.contains('met|axpin|diage|gluci|glucophage|metabe', case=False).any(), axis=1)
    med["insulin_med_timepoint"]=med.apply(lambda x: x.str.contains('insul|lantus', case=False).any(), axis=1)
    med.rename({"glp1ra_med": "glp1_agonist_timepoint",
                "mra_med": "mra",
                "htn_med___5": "diuretic",
                "htn_med___4": "ca_channel_blocker",
                "htn_med___3": "beta_blocker",
                "htn_med___6": "statin",
                "fibrates_med": "fibrates",
                "meds_weight_type___1": "topiramate",
                "meds_weight_type___2": "phentermine",
                "uric_acid_med": "uric_acid_med"
                },
                inplace=True, axis=1)
    meds=["sglti_timepoint", "raasi_timepoint", "metformin_timepoint", "insulin_med_timepoint", "glp1_agonist_timepoint", "mra", "diuretic", "ca_channel_blocker", "beta_blocker", "statin", "fibrates", "topiramate", "phentermine", "uric_acid_med"]
    med[meds] = med[meds].applymap(lambda x: "Yes" if x == "1" or x is True else ("No" if pd.notna(x) and x != "" else x))
    # SGLT2i, RAASi, Metformin
    med = med[["record_id", "visit", "sglti_timepoint", "raasi_timepoint", "metformin_timepoint", "insulin_med_timepoint",\
    "glp1_agonist_timepoint", "mra", "diuretic", "ca_channel_blocker", "beta_blocker", "statin", "fibrates", "topiramate", "phentermine", "uric_acid_med", "birthweight"]].copy()
    med["procedure"] = "medications"
    med["date"] = screen["date"]

    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(
        fields=var))
    phys["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"}, 
        inplace=True)
    phys = phys.rename(columns={"redcap_event_name": "visit"})
    # Replace missing values
    phys.replace(rep, np.nan, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sysbp": "sbp", "diasbp": "dbp"}, inplace=True, axis=1)
    phys.drop(["norm", "abnormal", "age"], axis=1, inplace=True)
    #phys=phys.loc[:, ~phys.columns.str.startswith("tan_")].copy()
    phys["acan"].replace({"1": "Yes", "0": "No"}, inplace=True)
    phys["acan_sev"].replace({"1": "Mild", "2": "Moderate", "3": "Severe"}, inplace=True)
    phys["procedure"] = "physical_exam"

    # --------------------------------------------------------------------------
    # IVGTT
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"].str.startswith("ivgtt"), "field_name"]] + [v for v in meta.loc[meta["form_name"].str.startswith("minmod"), "field_name"]]
    ivgtt = pd.DataFrame(proj.export_records(fields=var))
    ivgtt["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"}, 
        inplace=True)
    ivgtt = ivgtt.rename(columns={"redcap_event_name": "visit"})    
    ivgtt.replace(rep, np.nan, inplace=True)  # Replace missing values
    ivgtt.drop({"bl_lipid", "visit_upt", "visit_uptresult", "visit_npo", "visit_bgl_250",
                "visit_24hrurine", "ivgtt_pe", "ivgtt_wc", "baseline_vitals", "ivgtt_bl_labs",
                "bls_d25_end", "bls_ins_end"}, axis=1, inplace=True)
    ivgtt.columns = ivgtt.columns.str.replace(
        r"visit_|_iv", "", regex=True)
    ivgtt.columns = ivgtt.columns.str.replace(
        r"neg", "minus_", regex=True)
    ivgtt.columns = ivgtt.columns.str.replace(
        r"ins_", "insulin_", regex=True)        
    ivgtt.rename({"am_bg": "fbg", 
                  "bl_adipo": "adipo_base", 
                  "bl_leptin": "leptin_base",
                  "bl_glucose": "glucose_bl", 
                  "bl_insulin": "fasting_insulin",
                  "bl_cpep": "cpeptide", 
                  "ast": "gotast_base", 
                  "alt": "gptalt_base",
                  "bl_cholesterol": "chol_base", 
                  "bl_triglycerides": "triglycerides", 
                  "bl_hdl": "hdl_base", 
                  "bl_ldl": "ldl_base", 
                  "bl_nonhdl": "nonhdl_base",
                  "bls_d25_st": "d25_bolus_time", 
                  "bls_ins_st": "insulin_bolus_time",
                  "bl_estrad": "estrad",
                  "bl_tot_test" : "tot_test",
                  "bl_free_test" : "free_test"},
              axis=1, inplace=True)
    # Insulin
    # ins=list(ivgtt.loc[:, ivgtt.columns.str.startswith("insulin_")].columns.values)
    # ivgtt[ins] = ivgtt[ins].apply(
    #     pd.to_numeric, errors='coerce')
    # ivgtt["steady_state_insulin"] = ivgtt[['insulin_220',
    #                                        'insulin_230', 'insulin_240', 'insulin_245', 'insulin_249',
    #                                        'insulin_250', 'insulin_252', 'insulin_253', 'insulin_254',
    #                                        'insulin_255']].mean(axis=1) * 6
    # # C peptide
    # cpep = [c for c in ivgtt.columns if "cpeptide" in c]
    # ivgtt[cpep] = ivgtt[cpep].apply(
    #     pd.to_numeric, errors='coerce')
    # ivgtt["steady_state_cpeptide"] = ivgtt[['cpeptide_220', 'cpeptide_230',
    #                                         'cpeptide_240', 'cpeptide_245',
    #                                         'cpeptide_249', 'cpeptide_250',
    #                                         'cpeptide_252', 'cpeptide_253',
    #                                         'cpeptide_254', 'cpeptide_255']].mean(axis=1)
    # # ACPRg
    # ivgtt["acprg"] = ivgtt[['cpeptide_2', 'cpeptide_4',
    #                         'cpeptide_6', 'cpeptide_8', 'cpeptide_10']].mean(axis=1) - ivgtt[['cpeptide_minus_10', 'cpeptide_minus_5']].mean(axis=1)
    # # AIRg (USE FROM MINMOD)
    # ivgtt["airg"] = ivgtt[['insulin_2', 'insulin_4',
    #                        'insulin_6', 'insulin_8', 'insulin_10']].mean(axis=1) * 6 - ivgtt[['insulin_minus_10', 'insulin_minus_5']].mean(axis=1) * 6
    # # DI (USE FROM MINMOD)
    # ivgtt["di"] = (ivgtt["gir_190"] /
    #                ivgtt["steady_state_insulin"]) * ivgtt["airg"]
              
    ivgtt["procedure"] = "ivgtt"
    # --------------------------------------------------------------------------
    # DXA
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "study_visit_dxa_scan", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    dxa["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"}, 
        inplace=True)
    dxa = dxa.rename(columns={"redcap_event_name": "visit"})       
    # Replace missing values
    dxa.replace(rep, np.nan, inplace=True)
    dxa.columns = dxa.columns.str.replace(
        r"dxa_", "dexa_", regex=True)
    # DEXA Redcap form flipped weight/height vars
    dxa.rename({"bmd": "dexa_bone_mineral_density", "dexa_date": "date",
                "dexa_height": "weight", "dexa_weight": "height", 
                "bodyfat_percent": "dexa_body_fat", "leanmass_percent": "dexa_lean_mass",
                "trunkmass_percent": "dexa_trunk_mass", "fatmass_kg": "dexa_fat_kg", 
                "leanmass_kg": "dexa_lean_kg", "trunkmass_kg": "dexa_trunk_kg"}, axis=1, inplace=True)
    dxa.drop(["dexa_age"], axis=1, inplace=True)
    dxa["procedure"] = "dexa"

    # --------------------------------------------------------------------------
    # Renal Clearance Testing
    # --------------------------------------------------------------------------

    var = ["record_id"] + ["group"] + ["phys_map"] + ["dxa_height"] + ["dxa_weight"] + [v for v in meta.loc[meta["form_name"] == 
    "study_visit_renal_clearance_testing", "field_name"]] +[v for v in meta.loc[meta["form_name"] == 
    "renal_clearance_baseline_labs", "field_name"]]
    rct = pd.DataFrame(proj.export_records(fields=var))
    # Replace missing values
    rct.replace(rep, np.nan, inplace=True)
    rct["group"] = rct.groupby(["record_id"])["group"].ffill()
    rct = rct.loc[rct["redcap_event_name"] != "screening_arm_1"].copy()
    rct["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"},  
        inplace=True)
    rct = rct.rename(columns={"redcap_event_name": "visit"})       
    rct["procedure"] = "renal_clearance_testing"
    rename = {"gfr_raw": "gfr_raw_plasma_urine", "gfr_bsa": "gfr_bsa_plasma_urine",
              "erpf_raw": "erpf_raw_plasma_urine", "erpf": "erpf_bsa_plasma_urine",
              "gfr_15mgmin": "gfr_raw_plasma", "gfrbsa": "gfr_bsa_plasma",
              "erpf_pah_85": "erpf_raw_plasma", "erpfbsa": "erpf_bsa_plasma",
              "pahbsa": "pah_clear_bsa", "phys_map": "map", "bolus_ioh": "iohexol_vol", 
              "ioh_bol_com": "iohexol_time", "infusion_ioh": "iohexol_infusion_vol",
              "iohexol_yn": "iohexol_bolus", "pah_yn": "pah_bolus", "bolus_pah": "pah_vol",
              "infusion_pah": "pah_infusion_vol", "pah_bol_com": "pah_time",
              "dxa_height": "weight", "dxa_weight": "height"}
    rct.rename(rename, axis=1, inplace=True)
    rct = rct[rct.columns[~rct.columns.str.endswith('_yn')]].copy()
    rct.columns = rct.columns.str.replace(
        r"bl_|pi_|kim_", "", regex=True)
    # Calculate variables
    rct_vars = ["gfr_raw_plasma", "erpf_raw_plasma", "tot_protein", "map", "hct"]
    rct[rct_vars] = rct[rct_vars].apply(pd.to_numeric, errors='coerce')
    rct["erpf_raw_plasma_seconds"] = rct["erpf_raw_plasma"]/60
    rct["gfr_raw_plasma_seconds"] = rct["gfr_raw_plasma"]/60
    # Filtration Fraction
    rct["ff"] = rct["gfr_raw_plasma"]/rct["erpf_raw_plasma"] 
    # Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
    rct["kfg"] = np.select([rct["group"].eq("1"), rct["group"].eq("2"), rct["group"].eq("3")], [0.1012, 0.1733, 0.1733]) 
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
    rct.rename(rename, axis=1, inplace=True)
    rct.drop(["rbf_seconds", "erpf_raw_plasma_seconds", "group", "rc_labs", "map", "kfg", 
              "gfr_raw_plasma_seconds"], axis=1, inplace=True)
    rct["date"] = ivgtt["date"]

    # --------------------------------------------------------------------------
    # CGM
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "cgm_data", "field_name"]]
    cgm = pd.DataFrame(proj.export_records(fields=var))
    cgm["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"}, 
        inplace=True)
    cgm = cgm.rename(columns={"redcap_event_name": "visit"})   
    cgm.replace(rep, np.nan, inplace=True)  # Replace missing values
    cgm["cgm_yn"].replace({"1":"Yes", "0":"No"}, inplace=True)
    cgm["procedure"] = "cgm"
    cgm["date"] = phys["date"]

    # --------------------------------------------------------------------------
    # Outcomes
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "study_visit_boldasl_mri", "field_name"]]
    out = pd.DataFrame(proj.export_records(fields=var))
    out["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"},  
        inplace=True)
    out = out.rename(columns={"redcap_event_name": "visit"})       
    # Replace missing values
    out.replace(rep, np.nan, inplace=True)
    # Kidney outcomes like GFR, etc. were collected with the clamp, not
    # necessarily the day of the MRI
    bold_mri_cols = [c for c in out.columns if ("bold_" in c) or ("asl_" in c)]
    bold_mri = out[["record_id", "visit"] + bold_mri_cols].copy()
    out = out[list(set(out.columns).difference(bold_mri_cols))].copy()
    rename = {"volume_left": "left_kidney_volume_ml",
              "volume_right": "right_kidney_volume_ml",
              "mri_date": "date"}
    out.rename(rename, axis=1, inplace=True)
    out.drop(["mri_perf_data", "mri_diff_data", "mri_kid_scan", "mri_staff", "mri_comp", "mri_read_date"], axis=1, inplace=True)
    bold_mri.drop(["mri_bold_data"], axis=1, inplace=True)
    out["procedure"] = "clamp"
    bold_mri["procedure"] = "bold_mri"
    bold_mri["date"] = out["date"]

    # --------------------------------------------------------------------------
    # Pc MRI 2d
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                                               == "pc_mri_2d", "field_name"]]
    pc_mri = pd.DataFrame(proj.export_records(fields=var))
    pc_mri["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"}, inplace=True)
    pc_mri = pc_mri.rename(columns={"redcap_event_name": "visit"})       
    pc_mri.replace(rep, np.nan, inplace=True)  # Replace missing values
    pc_mri.rename({"hr_mri": "hr_mri_ra"}, axis=1, inplace=True)
    pc_mri["procedure"] = "pc_mri_2d"
    pc_mri["date"] = out["date"]

    # --------------------------------------------------------------------------
    # Missingness
    # --------------------------------------------------------------------------

    med.dropna(thresh=3, axis=0, inplace=True)
    screen.dropna(thresh=5, axis=0, inplace=True)
    phys.dropna(thresh=4, axis=0, inplace=True)
    ivgtt.dropna(thresh=4, axis=0, inplace=True)
    dxa.dropna(thresh=4, axis=0, inplace=True)
    rct.dropna(thresh=5, axis=0, inplace=True)
    cgm.dropna(thresh=4, axis=0, inplace=True)
    out.dropna(thresh=4, axis=0, inplace=True)
    bold_mri.dropna(thresh=4, axis=0, inplace=True)
    pc_mri.dropna(thresh=4, axis=0, inplace=True)

    # --------------------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------------------

    df = pd.concat([med, phys], join='outer', ignore_index=True)
    df = pd.concat([df, screen], join='outer', ignore_index=True)
    df = pd.concat([df, ivgtt], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, rct], join='outer', ignore_index=True)
    df = pd.concat([df, cgm], join='outer', ignore_index=True)
    df = pd.concat([df, bold_mri], join='outer', ignore_index=True)
    df = pd.concat([df, pc_mri], join='outer', ignore_index=True)
    df = pd.merge(df, out, how='outer')
    df = pd.merge(df, demo, how="outer")
    df = df.copy()

    # --------------------------------------------------------------------------
    # Reorganize
    # --------------------------------------------------------------------------

    df["study"] = "PANTHER"
    id_cols = ["record_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # Sort
    df.sort_values(["record_id", "visit", "procedure"], inplace=True)
    # Drop empty columns
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data
    return df
