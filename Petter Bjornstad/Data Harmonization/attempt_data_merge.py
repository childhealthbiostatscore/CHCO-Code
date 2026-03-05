"""
Merges attempt_dat.rds (from Kopah S3) into the harmonized dataset.
Called from data_harmonization.py after harmonize_data() is assembled.

"""


def merge_attempt_dat(harmonized):
    import os
    import sys
    import json
    import getpass
    import tempfile
    import boto3
    import pyreadr
    import pandas as pd
    import numpy as np
    from botocore.client import Config

    user = getpass.getuser()

    if user == "choiyej":
        base_data_path = "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/"
        keys_path = "/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json"
    elif user == "pylell":
        base_data_path = "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
        keys_path = ""
    elif user == "shivaniramesh":
        base_data_path = os.path.expanduser("~/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
        keys_path = "/Users/shivaniramesh/Desktop/keys.json"
    else:
        sys.exit(f"Unknown user: please specify root path for this user. (Detected user: {user})")

    # ------------------------------------------------------------------
    # Read attempt_dat.rds from Kopah S3
    # ------------------------------------------------------------------
    with open(keys_path, "r") as f:
        keys = json.load(f)

    s3 = boto3.client(
        "s3",
        endpoint_url="https://s3.kopah.uw.edu",
        aws_access_key_id=keys["MY_ACCESS_KEY"],
        aws_secret_access_key=keys["MY_SECRET_KEY"],
        config=Config(signature_version="s3v4"),
    )

    with tempfile.NamedTemporaryFile(suffix=".rds", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        s3.download_file(Bucket="attempt", Key="cleaned_data/attempt_dat.rds", Filename=tmp_path)
        result = pyreadr.read_r(tmp_path)
        attempt_dat = result[None]
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    # ------------------------------------------------------------------
    # Prep for merge: rename ID, recode visit to match harmonized labels
    # ------------------------------------------------------------------
    attempt_dat = attempt_dat.rename(columns={"subject_id": "record_id"})
    # subject_id is float64 (e.g. 10001.0) — convert via int to avoid '10001.0' string
    attempt_dat["record_id"] = attempt_dat["record_id"].astype(float).astype(int).astype(str)

    visit_map = {-4: "screening", 0: "baseline", 4: "4_weeks_post", 16: "4_months_post", 18: "follow_up"}
    attempt_dat["visit"] = pd.to_numeric(attempt_dat["visit"], errors="coerce").map(visit_map)

    # ------------------------------------------------------------------
    # Unit conversions and renames before merge
    # ------------------------------------------------------------------

    #  Anthropometrics — convert m → cm, rename to harmonized names
    attempt_dat["height"]          = pd.to_numeric(attempt_dat["height_m"],       errors="coerce") * 100
    attempt_dat["waistcm"]         = pd.to_numeric(attempt_dat["waist_m"],        errors="coerce") * 100
    attempt_dat["hipcm"]           = pd.to_numeric(attempt_dat["hip_m"],          errors="coerce") * 100
    attempt_dat["waist_hip_ratio"] = attempt_dat["waist_to_hip_ratio"]

    #  Blood labs — central
    attempt_dat["hba1c"]         = pd.to_numeric(attempt_dat["hba1c_percent"],            errors="coerce")           # % → % (same units)
    attempt_dat["cholesterol"]   = pd.to_numeric(attempt_dat["cholesterol_serum_mmoll"],  errors="coerce")           # mmol/L → mmol/L (same units)
    attempt_dat["hdl"]           = pd.to_numeric(attempt_dat["hdl_serum_mmoll"],          errors="coerce") * 38.67   # mmol/L → mg/dL
    attempt_dat["ldl"]           = pd.to_numeric(attempt_dat["ldl_serum_mmoll"],          errors="coerce") * 38.67   # mmol/L → mg/dL
    attempt_dat["triglycerides"] = pd.to_numeric(attempt_dat["triglycerides_serum_mmoll"],errors="coerce") * 88.5745 # mmol/L → mg/dL
    attempt_dat["creatinine_s"]  = pd.to_numeric(attempt_dat["creatinine_serum_umoll"],   errors="coerce") / 88.4194 # µmol/L → mg/dL
    attempt_dat["ca_base"]       = pd.to_numeric(attempt_dat["calcium_serum_mmoll"],      errors="coerce") * 4.0078  # mmol/L → mg/dL
    attempt_dat["c00009"]        = pd.to_numeric(attempt_dat["phosphate_serum_mmoll"],    errors="coerce")           # mmol/L (NOTE: c00009 units undocumented in harmonized)

    #  Blood labs — local
    attempt_dat["hct"]       = pd.to_numeric(attempt_dat["hct_local"],          errors="coerce")           # L/L → L/L (same units as ATTEMPT rows in harmonized)
    attempt_dat["pltct"]     = pd.to_numeric(attempt_dat["plt_local"],          errors="coerce")           # ×10⁹/L → same scale as harmonized pltct
    attempt_dat["sodium_s"]  = pd.to_numeric(attempt_dat["sodium_blood_local"], errors="coerce")           # mmol/L → same scale as harmonized sodium_s
    attempt_dat["bun"]       = pd.to_numeric(attempt_dat["bun_local"],          errors="coerce") * 2.8011  # mmol/L → mg/dL

    #  BOLD MRI — rename cortex/kidney to harmonized names; keep medulla as-is
    attempt_dat["bold_l_bl_cortex"] = pd.to_numeric(attempt_dat["mri_r2_cortex_l"], errors="coerce")
    attempt_dat["bold_r_bl_cortex"] = pd.to_numeric(attempt_dat["mri_r2_cortex_r"], errors="coerce")
    attempt_dat["bold_l_bl_kidney"] = pd.to_numeric(attempt_dat["mri_r2_kidney_l"], errors="coerce")
    attempt_dat["bold_r_bl_kidney"] = pd.to_numeric(attempt_dat["mri_r2_kidney_r"], errors="coerce")

    #  Glucose monitoring
    attempt_dat["cgm_type"] = attempt_dat["bgm_cgm_device"]
    # bgm_type skipped (numeric code only, no key); bgm_cgm_duration skipped

    # EMU — unit conversions and renames
    # Mean ACR: mg/mmol × 8.84 → mg/g → acr_u
    attempt_dat["acr_u"] = pd.to_numeric(attempt_dat["emu_urine_acr_mean"], errors="coerce") * 8.84

    # Mean albumin (mg/L = mcg/mL, direct) → microalbumin_u
    attempt_dat["microalbumin_u"] = (
        pd.to_numeric(attempt_dat["emu_1_albumin_mgl"], errors="coerce")
        .add(pd.to_numeric(attempt_dat["emu_2_albumin_mgl"], errors="coerce"), fill_value=np.nan)
        .add(pd.to_numeric(attempt_dat["emu_3_albumin_mgl"], errors="coerce"), fill_value=np.nan)
    )
    _emu_alb_n = (
        attempt_dat[["emu_1_albumin_mgl", "emu_2_albumin_mgl", "emu_3_albumin_mgl"]]
        .apply(pd.to_numeric, errors="coerce").notna().sum(axis=1)
    )
    attempt_dat["microalbumin_u"] = attempt_dat["microalbumin_u"] / _emu_alb_n.replace(0, np.nan)

    # Mean creatinine: µmol/L ÷ 88.4194 → mg/dL → creatinine_u
    _emu_cr = (
        attempt_dat[["emu_1_creatinine_umoll", "emu_2_creatinine_umoll", "emu_3_creatinine_umoll"]]
        .apply(pd.to_numeric, errors="coerce")
    )
    attempt_dat["creatinine_u"] = _emu_cr.mean(axis=1) / 88.4194

    for i in [1, 2, 3]:
        attempt_dat[f"emu_{i}_acr_u"]          = pd.to_numeric(attempt_dat[f"emu_{i}_acr_mgmmol"],        errors="coerce") * 8.84
        attempt_dat[f"emu_{i}_microalbumin_u"]  = pd.to_numeric(attempt_dat[f"emu_{i}_albumin_mgl"],       errors="coerce")        # mg/L = mcg/mL
        attempt_dat[f"emu_{i}_creatinine_u"]    = pd.to_numeric(attempt_dat[f"emu_{i}_creatinine_umoll"],  errors="coerce") / 88.4194
        attempt_dat[f"emu_{i}_date"]            = attempt_dat[f"date_emu_{i}"]

    # 24-hour urine — unit conversions
    attempt_dat["u24_na"] = (
        pd.to_numeric(attempt_dat["sodium_urine24h"], errors="coerce") * 1_000_000  # mmol/L → nmol/L
    )
    attempt_dat["u24_vl"] = (
        pd.to_numeric(attempt_dat["urine24h_volume_24h_litres"], errors="coerce") * 1000  # L → mL
    )
    attempt_dat["u24_mab"] = (
        pd.to_numeric(attempt_dat["microalbumin_urine24h"], errors="coerce") *
        pd.to_numeric(attempt_dat["urine24h_volume_24h_litres"], errors="coerce") *
        1_000_000  # g/L × L × 1e6 → mcg/24hr
    )

    #  Medical / family / smoking history — family history renames
    attempt_dat["famhx_t1d"]      = attempt_dat["fam_hx_t1d"]
    attempt_dat["fam_t1d"]        = attempt_dat["fam_hx_t1d_family"]
    attempt_dat["famhx_htn"]      = attempt_dat["fam_hx_hypertension"]
    attempt_dat["fam_htn"]        = attempt_dat["fam_hx_htn_family"]
    attempt_dat["famhx_hyperlipid"] = attempt_dat["fam_hx_dyslipidemia"]
    attempt_dat["fam_hyperlipid"]   = attempt_dat["fam_hx_dyslipidemia_family"]
    # fam_hx_ckd kept as-is (CKD != kidney failure/transplant; no exact harmonized match)
    # med_hx_1-5_* kept as-is (ATTEMPT-specific free-text medical history entries)
    # smoking_history kept as-is (no harmonized equivalent)

    #  Demographics
    attempt_dat["sex"] = attempt_dat["gender"]
    # dob_month/dob_year skipped — harmonized already has full dob for ATTEMPT
    # ethnicity_categorized skipped — already mapped to string labels in harmonized ethnicity

    #  Diabetes management
    # date_t1d_diagnosis comes as object/string ('%m/%d/%y') from pyreadr — convert to datetime
    attempt_dat["diabetes_dx_date"] = pd.to_datetime(attempt_dat["date_t1d_diagnosis"], format="%m/%d/%y", errors="coerce")
    # insulin_type/pump_insulin not in attempt_dat; all other diabetes mgmt vars skipped

    # Convert all other date columns from string ('%m/%d/%y') to datetime so they are
    # compatible with harmonized datetime columns and write correctly to CSV.
    _date_cols = [
        "date_emu_1", "date_emu_2", "date_emu_3",
        "date_urine24h_testing", "urine24h_start_date", "urine24h_stop_date",
    ]
    for _dc in _date_cols:
        if _dc in attempt_dat.columns:
            attempt_dat[_dc] = pd.to_datetime(attempt_dat[_dc], format="%m/%d/%y", errors="coerce")

    # ------------------------------------------------------------------
    # Variables approved for merge (fill in per section review)
    # ------------------------------------------------------------------

    # Anthropometrics
    # height_m, waist_m, hip_m: convert to cm and fill into existing harmonized cols
    # waist_to_hip_ratio: rename to waist_hip_ratio
    # bmi_z/bmi_percentile/waist_to_hip_ratio already exist in harmonized — skip
    anthro_vars = [
        "height", "waistcm", "hipcm",   # converted from m → cm
        "waist_hip_ratio",               # renamed from waist_to_hip_ratio
        "bsa_haycock_m2",
        "height_z", "height_percentile",
        "weight_z", "weight_percentile",
    ]

    #  Blood labs — central
    # Converted: hba1c, cholesterol, hdl, ldl, triglycerides, creatinine_s, ca_base, c00009
    # Added as-is: glucose_serum_mmoll, albumin_serum_gl, magnesium_serum_mmoll,
    #              uricacid_serum_umoll, bhb1_serum_mmoll, tsh_serum_miul, pth_serum_pmoll
    blood_central_vars = [
        "hba1c", "cholesterol", "hdl", "ldl", "triglycerides", "creatinine_s",
        "ca_base", "c00009",
        "glucose_serum_mmoll", "albumin_serum_gl", "magnesium_serum_mmoll",
        "uricacid_serum_umoll", "bhb1_serum_mmoll", "tsh_serum_miul", "pth_serum_pmoll",
    ]

    #  Blood labs — local
    # Converted: hct, pltct, sodium_s, bun
    # Added as-is: hgb_local, wbc_local, rbc_local, nrbc_local, mcv_local, mch_local,
    #              mchc_local, rdwcv_local, mpv_local, alp_local, alt_local (→ alt),
    #              bilirubin_local, creatine_kinase_local, potassium_blood_local,
    #              chloride_blood_local, bicarbonate_local, lipase_local
    # Skipped: serology_cs_yn, serology_cs_labs, serology_lab, date_collection_local,
    #          time_collection_local, comments
    blood_local_vars = [
        "hct", "pltct", "sodium_s", "bun",
        "hgb_local", "wbc_local", "rbc_local", "nrbc_local",
        "mcv_local", "mch_local", "mchc_local", "rdwcv_local", "mpv_local",
        "alt_local", "alp_local", "bilirubin_local", "creatine_kinase_local",
        "potassium_blood_local", "chloride_blood_local", "bicarbonate_local", "lipase_local",
    ]

    #  BOLD MRI
    # Renamed to harmonized names: bold_l/r_bl_cortex, bold_l/r_bl_kidney
    # Added as-is: mri_r2_medulla_l/r (different pipeline from harmonized bold_*_bl_medulla)
    # Skipped: mri_time (time only, no date), mri_technician, mri_radiologist,
    #          mri_reading_date, bold_mri, mri_kidneys_scanned
    bold_mri_vars = [
        "bold_l_bl_cortex", "bold_r_bl_cortex",
        "bold_l_bl_kidney", "bold_r_bl_kidney",
        "mri_r2_medulla_l", "mri_r2_medulla_r",
    ]

    #  Demographics
    # gender → sex; dob_month/dob_year/ethnicity_categorized skipped (redundant with harmonized)
    demo_vars = [
        "sex",
        "ethnicity_other",
    ]

    # Section 6: Diabetes management
    # date_t1d_diagnosis renamed to diabetes_dx_date; all other vars skipped
    diabetes_mgmt_vars = [
        "diabetes_dx_date",
    ]

    #  Drug compliance — all skipped
    drug_compliance_vars = []

    #  eGFR
    # Added: egfr_ckdepi40_cr (new equation), egfr_creatinine_serum_umoll_ageadjusted (new derived)
    # Skipped: egfr_ckidu25_cr/cysc (will be computed by data_harmonization.py as eGFR_CKiD_U25_*)
    #          egfr_creatinine_serum_mgdl (redundant with creatinine_s in Section 2)
    #          egfr_cystatin_c_serum_mgl (redundant with cystatin_c_s)
    #          egfr_k_cr, egfr_k_cysc, egfr_age, egfr_sex, egfr_height_m (intermediates/redundant)
    # NOTE: creatinine_s in harmonized has mixed units for ATTEMPT rows — investigate separately
    egfr_vars = [
        "egfr_ckdepi40_cr",
        "egfr_creatinine_serum_umoll_ageadjusted",
    ]

    #  Glucose monitoring
    # bgm_cgm_device renamed to cgm_type (matches harmonized dict label "If Yes, what brand?")
    # bgm_type skipped (numeric code only); bgm_cgm_duration skipped
    glucose_monitoring_vars = [
        "cgm_type",
    ]

    #  Medical / family / smoking history
    # Renamed to harmonized family_medical_history names: famhx_t1d/fam_t1d, famhx_htn/fam_htn,
    #   famhx_hyperlipid/fam_hyperlipid
    # Added as-is: fam_hx_ckd (no exact harmonized match), fam_hx_liver, fam_hx_other,
    #   fam_hx_other_details, fam_hx_other_family, med_hx_yn, medical_hx_number,
    #   med_hx_1-3_description/dx_year/active/treatment, smoking_history
    #   (med_hx_4/5_* dropped — too sparse)
    med_fam_history_vars = [
        # Family history — renamed
        "famhx_t1d", "fam_t1d",
        "famhx_htn", "fam_htn",
        "famhx_hyperlipid", "fam_hyperlipid",
        # Family history — ATTEMPT-specific (no harmonized match)
        "fam_hx_ckd", "fam_hx_ckd_family",
        "fam_hx_liver",
        "fam_hx_other", "fam_hx_other_details", "fam_hx_other_family",
        # Medical history — free-text ATTEMPT-specific entries
        "med_hx_yn", "medical_hx_number",
        "med_hx_1_description", "med_hx_1_dx_year", "med_hc_1_active", "med_hc_1_treatment",
        "med_hx_2_description", "med_hx_2_dx_year", "med_hx_2_active", "med_hx_2_treatment",
        "med_hx_3_description", "med_hx_3_dx_year", "med_hx_3_active", "med_hx_3_treatment",
        # med_hx_4/5_* dropped — too sparse (n≤2)
        # Smoking
        "smoking_history",
    ]

    # mGFR
    # All kept as-is (no renaming); ATTEMPT-specific calculation methods (BM adult/child/combined, Jodal)
    # mgfr_si/mgfr_si_adjusted are raw/BSA-adjusted plasma clearance but kept under original names
    mgfr_vars = [
        "mgfr_si", "mgfr_si_adjusted",
        "mgfr_bm_adult", "mgfr_bm_adult_adjusted", "mgfr_bm_adult_adjusted_cl1",
        "mgfr_bm_child", "mgfr_bm_child_adjusted", "mgfr_bm_child_adjusted_cl1",
        "mgfr_bm_combined", "mgfr_bm_combined_adjusted", "mgfr_bm_combined_adjusted_cl1",
        "mgfr_jodal", "mgfr_jodal_adjusted", "mgfr_jodal_bsa",
        "mgfr_jodal_f", "mgfr_jodal_f_bsa",
        "mgfr_ketone_iv", "mgfr_ketones_pre", "mgfr_ketones_t240",
        "mgfr_bsa_haycock",
    ]

    #  24-hour urine
    # Converted and renamed:
    #   sodium_urine24h (mmol/L) × 1e6 → u24_na (nmol/L)
    #   urine24h_volume_24h_litres (L) × 1000 → u24_vl (mL)
    #   microalbumin_urine24h (g/L) × volume (L) × 1e6 → u24_mab (mcg/24hr)
    # Added as-is: creatinine_urine24h (labeled mmol/L, likely µmol/L), creatinine_urine24h_gL
    #   (no mg/dL conversion available), glucose_urine24h, urea_urine24h,
    #   urine24h_yn, date_urine24h_testing, urine24h_start_date, urine24h_start_time,
    #   urine24h_stop_date, urine24h_stop_time, urine24h_time, urine_24hour_complete
    # Dropped: urine24h_volume_litres (redundant — u24_vl is converted version),
    #   microalbumin_urine24h_mgL (redundant — u24_mab is converted version),
    #   sample_id_urine24h (admin barcode), screen_urine_acr (all NA)
    # Skipped: tdid_basal_u_kg, tdid_bolus_u_kg, tdid_u_kg (insulin dose vars, not urine)
    urine_24h_vars = [
        # Renamed/converted
        "u24_na", "u24_vl", "u24_mab",
        # Added as-is
        "creatinine_urine24h", "creatinine_urine24h_gL",
        "glucose_urine24h", "urea_urine24h",
        "urine24h_yn", "date_urine24h_testing",
        "urine24h_start_date", "urine24h_start_time",
        "urine24h_stop_date", "urine24h_stop_time", "urine24h_time",
        "urine_24hour_complete",
    ]

    #  Early morning urine (EMU)
    # Means → harmonized names:
    #   emu_urine_acr_mean (mg/mmol) × 8.84 → acr_u (mg/g)
    #   mean(emu_1/2/3_albumin_mgl, mg/L=mcg/mL) → microalbumin_u
    #   mean(emu_1/2/3_creatinine_umoll) ÷ 88.4194 → creatinine_u (mg/dL)
    # Individual days → emu_N_acr_u (×8.84), emu_N_microalbumin_u (1:1), emu_N_creatinine_u (÷88.4194)
    # Dates: date_emu_N → emu_N_date
    # Added as-is: emu_urine_acr_mean_pooled, emu_acr_microalbuminuria, emu_collection_details,
    #   bhcg_urine_test, bhcg_urine_result, spot_urine_date, spot_urine_time
    # Dropped: emu_N_sample_id (admin barcodes)
    emu_vars = [
        # Means → harmonized
        "acr_u", "microalbumin_u", "creatinine_u",
        # Individual days
        "emu_1_acr_u", "emu_2_acr_u", "emu_3_acr_u",
        "emu_1_microalbumin_u", "emu_2_microalbumin_u", "emu_3_microalbumin_u",
        "emu_1_creatinine_u", "emu_2_creatinine_u", "emu_3_creatinine_u",
        "emu_1_date", "emu_2_date", "emu_3_date",
        # As-is
        "emu_urine_acr_mean_pooled", "emu_acr_microalbuminuria", "emu_collection_details",
        "bhcg_urine_test", "bhcg_urine_result",
        "spot_urine_date", "spot_urine_time",
    ]

    #  Local urinalysis
    # All kept as-is — dipstick is qualitative (0/1), CHC lab is semi-quantitative;
    # harmonized spot_* vars have no data and use a different collection method
    urinalysis_vars = [
        # Dipstick
        "bilirubin_urine_dip", "blood_urine_dip", "glucose_urine_dip",
        "ketones_urine_dip", "leukocytes_urine_dip", "nitrite_urine_dip",
        "ph_urine_dip", "protein_urine_dip", "spgravity_urine_dip", "urobilinogen_urine_dip",
        # CHC lab (semi-quantitative)
        "glucose_urinalysis_chc", "glucose_urinalysis_chc_value",
        "ketones_urinalysis_chc", "ketones_urinalysis_chc_value",
        "protein_urinalysis_chc", "protein_urinalysis_chc_value",
        "ph_urinalysis_chc", "spgravity_urinalysis_chc",
        "leukocyte_esterase_urinalysis_chc", "nitrite_urinalysis_chc", "blood_urinalysis_chc",
        # All-site composite
        "glucose_urinalysis_all", "ketones_urinalysis_all", "protein_urinalysis_all",
        "ph_urinalysis_all", "spgravity_urinalysis_all", "leukocytes_urinalysis_all",
        "nitrites_urinalysis_all", "hgb_blood_urinalysis_all",
        # Admin — dropped: urinalysis_lab, urinalysis_cs_yn, urinalysis_cs_labs
    ]

    #  Vitals (partial), CGM derived, biopsy, TKV
    # Checked against data: cardiac MRI (all n=0), long/circum/radial morphometry (all n=0),
    #   sys_bp/dys_bp/vitals_date/imaging_hr (n=0), protein (n=0) — all skipped (no data)
    cardiac_mri_vars = [
        # Vitals — only PWV has data; BP/date/HR all empty
        "pwv_cf_1", "pwv_cf_2",
        # CGM derived
        "cgm_tir", "cgm_window", "avg_ketones",
        # Biopsy (protein skipped — no data)
        "biopsy_yn", "GBM thickening", "Glomeruli number", "Glomeruli sclerosed",
        "Tubular atrophy", "Vessel pathology", "Pathology_report_ID", "LN2_ID",
        # Kidney morphometry — only TKV has data; long/circum/radial series all empty
        "tkv",
    ]

    vars_to_merge = (
        anthro_vars + blood_central_vars + blood_local_vars +
        bold_mri_vars + demo_vars + diabetes_mgmt_vars +
        drug_compliance_vars + egfr_vars + glucose_monitoring_vars +
        med_fam_history_vars + mgfr_vars + urine_24h_vars +
        emu_vars + urinalysis_vars + cardiac_mri_vars
    )

    

    # ------------------------------------------------------------------
    # Subset and merge — only update ATTEMPT rows, only new columns
    # ------------------------------------------------------------------
    if not vars_to_merge:
        print("merge_attempt_dat: no variables approved yet, returning harmonized unchanged.")
        return harmonized
    #sanity check for size of harmonized[study == "ATTEMPT"]
    attempt_count = (harmonized["study"] == "ATTEMPT").sum()
    print(f"merge_attempt_dat: {attempt_count} ATTEMPT rows in harmonized before merge.")
    
    keep_cols = ["record_id", "visit"] + [v for v in vars_to_merge if v in attempt_dat.columns]
    attempt_subset = attempt_dat[keep_cols].drop_duplicates(subset=["record_id", "visit"])

    # Only merge into ATTEMPT rows; new columns will be NaN for all other studies
    new_cols = [c for c in attempt_subset.columns if c not in harmonized.columns]
    attempt_rows = harmonized["study"] == "ATTEMPT"

    harmonized = harmonized.merge(
        attempt_subset[["record_id", "visit"] + new_cols],
        on=["record_id", "visit"],
        how="left",
    )

    # For columns that already exist in harmonized, fill ATTEMPT gaps with attempt_dat values
    existing_cols = [c for c in vars_to_merge if c in harmonized.columns and c not in new_cols]
    if existing_cols:
        attempt_fill = attempt_subset[["record_id", "visit"] + existing_cols].copy()
        attempt_fill.columns = ["record_id", "visit"] + [f"_attempt_{c}" for c in existing_cols]
        harmonized = harmonized.merge(attempt_fill, on=["record_id", "visit"], how="left")
        for c in existing_cols:
            harmonized[c] = harmonized[c].combine_first(harmonized[f"_attempt_{c}"])
            harmonized.drop(columns=[f"_attempt_{c}"], inplace=True)

    print(f"merge_attempt_dat: merged {len(vars_to_merge)} variables into harmonized dataset.")
    return harmonized
