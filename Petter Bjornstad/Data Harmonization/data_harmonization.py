"""
This code is designed to pull data from multiple REDCap projects and harmonize
the data in a single dataset. Some studies are cross-sectional but include
measures at multiple visits, and some studies are longitudinal. So, this code
outputs data in a "semi-long" format with one row per study procedure, and a
visit column for longitudinal clustering. The data cleaning process for each 
individual dataset is a separate function, and this function puts them together 
and performs some formatting tweaks.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def harmonize_data():
    # Libraries
    import os
    import sys
    sys.path.insert(0, os.path.expanduser('~') +
                    "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    import pandas as pd
    import numpy as np
    from natsort import natsorted, ns
    from casper import clean_casper
    from coffee import clean_coffee
    from crocodile import clean_crocodile
    from improve import clean_improve
    from penguin import clean_penguin
    from renal_heir import clean_renal_heir
    from renal_heiritage import clean_renal_heiritage
    from panther import clean_panther
    from panda import clean_panda
    from attempt import clean_attempt
    from rpc2 import clean_rpc2_redcap
    from ultra import clean_ultra
    from sweetheart import clean_sweetheart
    from harmonization_functions import calc_egfr, create_study_id_columns, biopsy_merge
    import getpass
    user = getpass.getuser() 
    # Use individual data functions to import cleaned DFs
    casper = clean_casper()
    coffee = clean_coffee()
    crocodile = clean_crocodile()
    improve = clean_improve()
    penguin = clean_penguin()
    renal_heir = clean_renal_heir()
    renal_heiritage = clean_renal_heiritage()
    panther = clean_panther()
    panda = clean_panda()
    attempt = clean_attempt()
    rpc2 = clean_rpc2_redcap()
    ultra = clean_ultra()
    sweetheart = clean_sweetheart()
    # Merge
    harmonized = pd.concat([casper, coffee], join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, crocodile],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, improve],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, penguin],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, renal_heir],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, renal_heiritage],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, panther],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, panda],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, attempt],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, rpc2],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, ultra],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, sweetheart],
                           join='outer', ignore_index=True)
    

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
    dictionary = pd.read_csv(base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv")

    

                           
    # Fix levels of categorical variables
    harmonized["visit"] = \
        pd.Categorical(harmonized["visit"],
                       categories=['screening', 'baseline', 'pre_surgery',
                                   '3_months_post_surgery', '4_months_post', '12_months_post_surgery',
                                   'year_1', 'year_2', 'year_3', 'year_4' , "post_biopsy","Post-biopsy GFR MRI","Med Dispense 1", "Med Dispense 2"],
                       ordered=True)
    harmonized["race"].replace(
        ["American Indian or Alaskan Native & White",
         "Black or African American & White",
         'American Indian or Alaskan Native & Black or African American',
         'Asian & White'], "More Than One", inplace=True)
    harmonized["race"].replace(
        {'Black/African American': 'Black or African American',
         "": "Unknown"}, inplace=True)
    harmonized["ethnicity"].replace({"": "Unknown"}, inplace=True)
    race_ethnicity = harmonized["race"] + \
        ", " + harmonized["ethnicity"]
    harmonized = pd.concat([harmonized, race_ethnicity], axis=1)
    harmonized.rename({0: "race_ethnicity"}, axis=1, inplace=True)

    # Replace blanks with missing
    harmonized.replace("", np.nan, inplace=True)

    operational = ['study', 'visit', 'record_id', 'procedure', 'date']
    dictionary.loc[dictionary['variable_name'].isin(operational), 'form_name'] = 'operational'
    dictionary.loc[dictionary['variable_name']== "race_ethnicity", 'form_name'] = 'demographics'



    # Date variables


    dates = ["dob", "date", "diabetes_dx_date"]
    for col in harmonized.columns:
        try:
            harmonized[col] = pd.to_numeric(harmonized[col], errors="ignore")
        except Exception as e:
            print(f"Skipping numeric conversion for {col}: {e}")

    def normalize_date(date_str):
        """
        Normalize date strings into YYYY-MM-DD.
        Handles dash or slash separators, missing day/month, and 2-digit years.
        """
        if not isinstance(date_str, str) or not date_str.strip():
            return None

        s = date_str.strip()

        # Handle common placeholders
        if s.lower() in {"na", "none", "null", "missing", "nan"}:
            return None

        # Replace slashes with dashes
        s = s.replace("/", "-")

        parts = s.split("-")
        nums = []
        for p in parts:
            try:
                nums.append(int(p))
            except ValueError:
                return None

        year = month = day = None

        # If only 1 number: assume year
        if len(nums) == 1:
            year = nums[0]
            month = 1
            day = 1
        else:
            for n in nums:
                if n > 31 and year is None:
                    year = n
                elif 1 <= n <= 12 and month is None:
                    month = n
                elif 1 <= n <= 31 and day is None:
                    day = n

            # Fill missing month/day
            if year is None:
                # fallback: pick largest number as year
                year = max(nums)
            if month is None:
                month = 1
            if day is None:
                day = 1

        # Handle 2-digit years (assume 1900s if > current year?)
        if year < 100:
            year += 2000 if year <= 25 else 1900  # adjust based on your data

        return f"{year:04d}-{month:02d}-{day:02d}"
    harmonized.loc[harmonized['record_id'] == 'RH2-22-T', 'diabetes_dx_date'] = '2013-06-26'
    harmonized.loc[harmonized['record_id'] == 'RH2-42-T', 'diabetes_dx_date'] = '2004-10-06'
    harmonized.loc[harmonized['record_id'] == 'RH2-13-O', 'diabetes_dx_date'] = '2017-03-24'


    
    for col in dates:
        if col in harmonized.columns:
            harmonized[col] = harmonized[col].apply(normalize_date)

    for col in dates:
        if col in harmonized.columns:
            harmonized[col] = pd.to_datetime(harmonized[col], errors="coerce")
            print(f"{col}: {harmonized[col].notna().sum()} non-null datetime values")
    
    # harmonized[dates] = \
    #     harmonized[dates].apply(pd.to_datetime, errors='coerce')

    # ----------------------
    # Calculated variables
    # ----------------------


    # Age
    print("\nStudies with dob:")
    print(harmonized.groupby('study')['dob'].apply(lambda x: x.notna().sum()))

    print("\nStudies with date:")
    print(harmonized.groupby('study')['date'].apply(lambda x: x.notna().sum())) 
    harmonized = harmonized.drop(columns=['age'], errors='ignore')
    age = round((harmonized["date"] - harmonized["dob"]).dt.days / 365.25, 2)
    harmonized = pd.concat([harmonized, age], axis=1)
    harmonized.rename({0: "age"}, axis=1, inplace=True)

    print("\nStudies with age:")
    print(harmonized.groupby('study')['age'].apply(lambda x: x.notna().sum()))
    # set column x to 'a' where column y == 'b'
    dictionary.loc[dictionary['variable_name'] == 'age', 'form_name'] = 'demographics'
    # BMI
    harmonized["bmi"] = pd.to_numeric(harmonized["weight"])/((pd.to_numeric(harmonized["height"])/100)**2)
    dictionary.loc[dictionary['variable_name'] == 'bmi', 'form_name'] = 'physical_exam'

    # Diabetes duration
    disease_duration = \
        round((harmonized["date"] -
              harmonized["diabetes_dx_date"]).dt.days / 365.25, 2)

    harmonized = pd.concat([harmonized, disease_duration], axis=1)
    harmonized.rename({0: "diabetes_duration"}, axis=1, inplace=True)
    dictionary.loc[dictionary['variable_name'] == 'diabetes_duration', 'form_name'] = 'medical_history'
    # eGFR
    if "age" in harmonized.columns:
        harmonized = harmonized.loc[:, ~harmonized.columns.duplicated()]
        harmonized["age"] = pd.to_numeric(harmonized["age"], errors="coerce")
    else:
        print("WARNING: 'age' column not found in harmonized data.")
    harmonized["age"] = pd.to_numeric(harmonized["age"], errors="coerce")
    harmonized["creatinine_s"] = pd.to_numeric(harmonized["creatinine_s"], errors="coerce")
    harmonized["cystatin_c_s"] = pd.to_numeric(harmonized["cystatin_c_s"], errors="coerce")
    harmonized["bun"] = pd.to_numeric(harmonized["bun"], errors="coerce")
    harmonized["height"] = pd.to_numeric(harmonized["height"], errors="coerce")
    harmonized = calc_egfr(harmonized, age="age",
                           serum_creatinine="creatinine_s", cystatin_c="cystatin_c_s",
                           bun="bun", height="height", sex="sex", male="Male", female="Female", alpha=0.5)
    egfr = ['eGFR_Schwartz', 'eGFR_bedside_Schwartz', 'eGFR_Zap', 'eGFR_fas_cr',
         'eGFR_fas_cr_cysc', 'eGFR_CKD_epi',
         'eGFR_CKiD_U25_Creat', 'eGFR_CKiD_U25_CystatinC', 'eGFR_CKiD_U25_avg']
    dictionary.loc[dictionary['variable_name'].isin(egfr), 'form_name'] = 'eGFR'

    # Kidney volume
    harmonized["total_kidney_volume_ml"] = \
        harmonized["left_kidney_volume_ml"] + \
        harmonized["right_kidney_volume_ml"]
    harmonized["total_kidney_volume_ml_manual"] = harmonized.apply(lambda row: row["volume_left_manual"] + row["volume_right_manual"], axis=1)
    harmonized = harmonized.assign(
        ht_adj_tkv = harmonized["total_kidney_volume_ml"] / (harmonized.groupby("record_id")["height"].transform("mean") / 100))
    harmonized = harmonized.assign(
        ht_adj_tkv_manual = harmonized["total_kidney_volume_ml_manual"] / (harmonized.groupby("record_id")["height"].transform("mean") / 100))
    
    kidney_vol_vars = ['total_kidney_volume_ml', 'total_kidney_volume_ml_manual', 'left_kidney_volume_ml', 'right_kidney_volume_ml']
    
    dictionary.loc[dictionary['variable_name'].isin(kidney_vol_vars), 'form_name'] = 'fmri'

    # PCASL
    harmonized["pcasl3d_left"] = pd.to_numeric(harmonized["pcasl3d_left"], errors='coerce')
    harmonized["pcasl3d_right"] = pd.to_numeric(harmonized["pcasl3d_right"], errors='coerce')
    harmonized["avg_pcascl"]= \
        harmonized[["pcasl3d_left", "pcasl3d_right"]].apply(lambda x: x.mean(), axis=1)
    pcasl_vars = ['pcasl3d_left', 'pcasl3d_right', 'avg_pcascl']
    dictionary.loc[dictionary['variable_name'].isin(pcasl_vars), 'form_name'] = 'study_visit_boldasl_mri'

    # Average R2*
    harmonized["avg_k_r2"]= \
        harmonized[["bold_l_bl_kidney", "bold_r_bl_kidney"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_c_r2"]= \
        harmonized[["bold_l_bl_cortex", "bold_r_bl_cortex"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_r2"]= \
        harmonized[["bold_l_bl_medulla", "bold_r_bl_medulla"]].apply(lambda x: x.mean(), axis=1) 
    outcomes_vars = ['avg_k_r2', 'avg_c_r2', 'avg_m_r2']
    dictionary.loc[dictionary['variable_name'].isin(outcomes_vars), 'form_name'] = 'outcomes'
    
    # Average T1
    harmonized["avg_k_t1"]= \
        harmonized[["bold_r_t1_kidney", "bold_l_t1_kidney"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_c_t1"]= \
        harmonized[["bold_l_t1_cortex", "bold_r_t1_cortex"]].apply(lambda x: x.mean(), axis=1)
    
    boldasl_mri_vars = ['avg_k_t1', 'avg_c_t1', 'bold_l_t1_kidney', 'bold_r_t1_kidney', 'bold_r_t1_cortex', 'bold_l_t1_cortex']
    dictionary.loc[dictionary['variable_name'].isin(boldasl_mri_vars), 'form_name'] = 'study_visit_boldasl_mri'
    
    # Average ADC
    harmonized["avg_c_adc"] = \
        harmonized[["adc_left", "adc_right"]].apply(lambda x: x.mean(), axis=1)
    boldasl_mri_vars = ['avg_c_adc']
    dictionary.loc[dictionary['variable_name'].isin(boldasl_mri_vars), 'form_name'] = 'study_visit_boldasl_mri'
    
    # # Average voxelwise
    harmonized["avg_m_k2_wo_cyst_vw"] = \
        harmonized[["lm_k2_wo_cyst_vw", "rm_k2_wo_cyst_vw"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_c_k2_wo_cyst_vw"] = \
        harmonized[["lc_k2_wo_cyst_vw", "rc_k2_wo_cyst_vw"]].apply(lambda x: x.mean(), axis=1)
    voxelwise_vars = ['avg_m_k2_wo_cyst_vw', 'avg_c_k2_wo_cyst_vw']
    dictionary.loc[dictionary['variable_name'].isin(voxelwise_vars), 'form_name'] = 'voxelwise'
    
    # Average K1
    harmonized["avg_c_k1"]= \
        harmonized[["lc_k1", "rc_k1"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_k1"]= \
        harmonized[["lm_k1", "rm_k1"]].apply(lambda x: x.mean(), axis=1)  
    optional_pet_vars = ['avg_c_k1', 'avg_m_k1']
    dictionary.loc[dictionary['variable_name'].isin(optional_pet_vars), 'form_name'] = 'optional_pet_scan'
     
    # Average k2
    harmonized["avg_c_k2"]= \
        harmonized[["lc_k2", "rc_k2"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_k2"]= \
        harmonized[["lm_k2", "rm_k2"]].apply(lambda x: x.mean(), axis=1) 
    optional_pet_vars = ['avg_c_k2', 'avg_m_k2']
    dictionary.loc[dictionary['variable_name'].isin(optional_pet_vars), 'form_name'] = 'optional_pet_scan'
    
    # Average F
    harmonized["avg_c_f"]= \
        harmonized[["lc_f", "rc_f"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_f"]= \
        harmonized[["lm_f", "rm_f"]].apply(lambda x: x.mean(), axis=1)  
    optional_pet_vars = ['avg_c_f', 'avg_m_f']
    dictionary.loc[dictionary['variable_name'].isin(optional_pet_vars), 'form_name'] = 'optional_pet_scan'       
    
        
    # Calculate FSOC = bl_bold - pf_bold
    cols = [c for c in harmonized.columns if "_bl_" in c] + \
        [c for c in harmonized.columns if "_pf_" in c]
    harmonized[cols] = harmonized[cols].apply(
        pd.to_numeric, errors='coerce', axis=1)
    harmonized = harmonized.assign(
        fsoc_r_cortex=harmonized["bold_r_bl_cortex"] -
        harmonized["bold_r_pf_cortex"],
        fsoc_r_medulla=harmonized["bold_r_bl_medulla"] -
        harmonized["bold_r_pf_medulla"],
        fsoc_r_kidney=harmonized["bold_r_bl_kidney"] -
        harmonized["bold_r_pf_kidney"],
        fsoc_l_cortex=harmonized["bold_l_bl_cortex"] -
        harmonized["bold_l_pf_cortex"],
        fsoc_l_medulla=harmonized["bold_l_bl_medulla"] -
        harmonized["bold_l_pf_medulla"],
        fsoc_l_kidney=harmonized["bold_l_bl_kidney"] -
        harmonized["bold_l_pf_kidney"])
    fsoc_vars = ['fsoc_r_cortex', 'fsoc_r_medulla', 'fsoc_r_kidney', 'fsoc_l_cortex', 'fsoc_l_medulla', 'fsoc_l_kidney']
    dictionary.loc[dictionary['variable_name'].isin(fsoc_vars), 'form_name'] = 'FSOC'
   
    # Average FSOC
    harmonized["avg_k_fsoc"]= \
        harmonized[["fsoc_l_kidney", "fsoc_r_kidney"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_c_fsoc"]= \
        harmonized[["fsoc_l_cortex", "fsoc_r_cortex"]].apply(lambda x: x.mean(), axis=1)
    harmonized["avg_m_fsoc"]= \
        harmonized[["fsoc_l_medulla", "fsoc_r_medulla"]].apply(lambda x: x.mean(), axis=1)
    fsoc_vars = ['avg_k_fsoc', 'avg_c_fsoc', 'avg_m_fsoc']
    dictionary.loc[dictionary['variable_name'].isin(fsoc_vars), 'form_name'] = 'FSOC'
            
    # UACR
    harmonized["acr_u"] = \
        pd.to_numeric(harmonized["microalbumin_u"], errors="coerce") * 100 / \
        pd.to_numeric(harmonized["creatinine_u"], errors="coerce")
    uacr_vars = ['acr_u']
    dictionary.loc[dictionary['variable_name'].isin(uacr_vars), 'form_name'] = 'UACR'
    
    
    # Albuminuria
    alb = []
    for a in harmonized["acr_u"]:
        if a < 30:
            alb.append("A1")
        elif a >= 30 and a <= 300:
            alb.append("A2")
        elif a > 300:
            alb.append("A3")
        else:
            alb.append(np.nan)
    harmonized["albuminuria_cat"] = alb
    harmonized["elevated_albuminuria"] = pd.cut(
        harmonized["acr_u"], [-float("inf"), 30, float("inf")], right=False, labels=["No", "Yes"])
    albuminuria_vars = ['albuminuria_cat', 'elevated_albuminuria']
    dictionary.loc[dictionary['variable_name'].isin(albuminuria_vars), 'form_name'] = 'Albuminuria'
    

    # FFA suppression negative to be 0
    harmonized["ffa_suppression"] = np.where(
        harmonized["ffa_suppression"] < 0, 0, harmonized["ffa_suppression"])
    harmonized["p1_ffa_suppression"] = np.where(
        harmonized["p1_ffa_suppression"] < 0, 0, harmonized["p1_ffa_suppression"])
    harmonized["p2_ffa_suppression"] = np.where(
        harmonized["p2_ffa_suppression"] < 0, 0, harmonized["p2_ffa_suppression"])
    # FFA suppression combined
    harmonized = \
        harmonized.assign(ffa_suppression_combined=harmonized["ffa_suppression"].where(
            harmonized["ffa_suppression"].notnull(), harmonized["p2_ffa_suppression"]))
    ffa_vars = ['ffa_suppression', 'p1_ffa_suppression', 'p2_ffa_suppression', 'ffa_suppression_combined']
    dictionary.loc[dictionary['variable_name'].isin(ffa_vars), 'form_name'] = 'FFA'
    

    # Fasting Insulin
    harmonized[["insulin_minus_20", "insulin_minus_10", "insulin_minus_5", "insulin_0"]] = harmonized[[
        "insulin_minus_20", "insulin_minus_10", "insulin_minus_5", "insulin_0"]].apply(pd.to_numeric)
    harmonized["fasting_insulin"] = \
        harmonized[["insulin_minus_20", "insulin_minus_10",
                    "insulin_minus_5", "insulin_0"]].apply(lambda x: x.mean(), axis=1)
    # Fasting FFA
    harmonized["fasting_ffa"] = \
        harmonized[["ffa_minus_20", "ffa_minus_10", "ffa_minus_5", "ffa_0"]].apply(
            lambda x: x.mean(), axis=1)
    # Fasting Glucose
    harmonized["fbg"] = \
        harmonized[["glucose_minus_120", "glucose_minus_90","glucose_minus_20",
        "glucose_minus_10", "glucose_minus_5", "glucose_0", "fbg"]].apply(
            lambda x: x.mean(), axis=1)
    clamp_vars = ['fasting_insulin', 'fasting_ffa', 'fbg']
    dictionary.loc[dictionary['variable_name'].isin(clamp_vars), 'form_name'] = 'clamp'
    dictionary.loc[dictionary['variable_name'] == 'mm_si', 'units'] = '(mu/I)^-1.min^-1'

    harmonized["m_i_p2_raw_lean"] = pd.to_numeric(harmonized["p2_raw_leanm"])/pd.to_numeric(harmonized["p2_steady_state_insulin"]) #CROC, PENGUIN
    harmonized["m_i_gir_190"] = pd.to_numeric(harmonized["gir_190"])/pd.to_numeric(harmonized["steady_state_insulin"]) #CROC?, IMPR, RH
    mi_vars = ['m_i_p2_raw_lean', 'm_i_gir_190', "m_i"]

    dictionary.loc[dictionary['variable_name'].isin(mi_vars), 'form_name'] = 'clamp'
    dictionary.loc[dictionary['variable_name'] == 'm_value', 'units'] = 'mg/kg lean/min'
    

    # HOMA-IR (https://link.springer.com/article/10.1007/BF00280883), FBG entered as mg/dL, converting to mmol/L (18 mg/dL = 1 mmol/L)
    harmonized = harmonized.assign(
        homa_ir=(harmonized["fasting_insulin"] * (harmonized["fbg"]/18))/22.5)
    # Adipose IR (fasting_insulin * fasting_ffa)
    harmonized = harmonized.assign(
        adipose_ir=harmonized["fasting_ffa"] * harmonized["fasting_insulin"])
    # SEARCH IS score eIS (estimated insulin sensitivity) (https://academic.oup.com/jcem/article/101/2/686/2811091#81467317)
    harmonized = harmonized.assign(
        search_eis = np.exp((4.64725 - (0.02032 * harmonized.groupby("record_id")["waistcm"].transform("mean")) - 
        (0.09779 * harmonized.groupby("record_id")["hba1c"].transform("mean")) - 
        (0.00235 * harmonized.groupby("record_id")["triglycerides"].transform("mean")))))
    clamp_vars = ['homa_ir', 'adipose_ir', 'search_eis']
    dictionary.loc[dictionary['variable_name'].isin(clamp_vars), 'form_name'] = 'clamp'
    
    # Merge copeptin values (CASPER, RH2, CROCODILE)
    # Step 1: Replace empty strings with NaN in both columns
    harmonized['copeptin'].replace("", np.nan, inplace=True)
    harmonized['pi_copeptin'].replace("", np.nan, inplace=True)

    # Step 2: Convert both columns to numeric (coerce errors to NaN)
    harmonized['copeptin'] = pd.to_numeric(harmonized['copeptin'], errors='coerce')
    harmonized['pi_copeptin'] = pd.to_numeric(harmonized['pi_copeptin'], errors='coerce')

    # Step 3: Move pi_copeptin value to copeptin if copeptin is missing but pi_copeptin is present
    mask_move = harmonized['copeptin'].isna() & harmonized['pi_copeptin'].notna()
    harmonized.loc[mask_move, 'copeptin'] = harmonized.loc[mask_move, 'pi_copeptin']
    harmonized.loc[mask_move, 'pi_copeptin'] = np.nan  # clear pi_copeptin after move

    # Step 4: If both have values but differ, keep copeptin and clear pi_copeptin
    mask_diff = harmonized['copeptin'].notna() & harmonized['pi_copeptin'].notna() & (harmonized['copeptin'] != harmonized['pi_copeptin'])
    harmonized.loc[mask_diff, 'pi_copeptin'] = np.nan

    # Step 5: If both have values and are equal, clear pi_copeptin (optional cleanup)
    mask_equal = harmonized['copeptin'].notna() & harmonized['pi_copeptin'].notna() & (harmonized['copeptin'] == harmonized['pi_copeptin'])
    harmonized.loc[mask_equal, 'pi_copeptin'] = np.nan

    # Step 6: Drop pi_copeptin column (all handled)
    harmonized.drop(columns=['pi_copeptin'], inplace=True)

    #aer_24 = (u24_mab * u24_vl) / 1440
    #bsa_dubois = 0.007184 × W^0.425 × H^0.725
    # aer_4_coltime = (urine_mab_250 * urine_vol) / 250
    harmonized["aer_24"] = (pd.to_numeric(harmonized["u24_mab"], errors="coerce") * pd.to_numeric(harmonized["u24_vl"], errors="coerce")) / 1440
    harmonized["aer_4_coltime"] = (pd.to_numeric(harmonized["urine_mab_250"], errors="coerce") * pd.to_numeric(harmonized["urine_vol"], errors="coerce")) / 250
    harmonized["bsa_dubois"] = 0.007184 * (pd.to_numeric(harmonized["weight"], errors="coerce")**0.425) * (pd.to_numeric(harmonized["height"], errors="coerce")**0.725)
   
    dictionary.loc[dictionary['variable_name'] == 'aer_24', 'units'] = 'mcg * mL / min'
    dictionary.loc[dictionary['variable_name'] == 'aer_4_coltime', 'units'] = 'mcg * mL / min'
    dictionary.loc[dictionary['variable_name'] == 'bsa_dubois', 'units'] = 'm^2'
    dictionary.loc[dictionary['variable_name'] == 'aer_24', 'label'] = 'AER over 24 hrs, mcg * mL / min'
    dictionary.loc[dictionary['variable_name'] == 'aer_4_coltime', 'label'] = 'AER over 250 mins, mcg * mL / min'
    dictionary.loc[dictionary['variable_name'] == 'bsa_dubois', 'label'] = 'BSA Calculated with Dubois formula, m^2'


    dictionary.loc[dictionary['variable_name'] == 'ace_inhibitor', 'form_name'] = 'medical_history'
    dictionary.loc[dictionary['variable_name'] == 'acprg', 'form_name'] = 'clamp'
    dictionary.loc[dictionary['variable_name'] == 'acr_u_pm', 'form_name'] = 'clamp'
    dictionary.loc[dictionary['variable_name'] == 'copeptin', 'form_name'] = 'copeptin'



    harmonized = biopsy_merge(harmonized)




    # Co-enroll IDs
    # casper_mrns = harmonized.loc[harmonized['study'] == 'CASPER', ['mrn', 'record_id']]
    # casper_id_map = dict(zip(casper_mrns['mrn'], casper_mrns['record_id']))
    # harmonized['casper_id'] = harmonized.apply(lambda row: casper_id_map[row['mrn']] if row['mrn'] in casper_id_map else '', axis=1)
    # coffee_mrns = harmonized.loc[harmonized['study'] == 'COFFEE', ['mrn', 'record_id']]
    # coffee_id_map = dict(zip(coffee_mrns['mrn'], coffee_mrns['record_id']))
    # harmonized['coffee_id'] = harmonized.apply(lambda row: coffee_id_map[row['mrn']] if row['mrn'] in coffee_id_map else '', axis=1)
    # croc_mrns = harmonized.loc[harmonized['study'] == 'CROCODILE', ['mrn', 'record_id']]
    # croc_id_map = dict(zip(croc_mrns['mrn'], croc_mrns['record_id']))
    # harmonized['croc_id'] = harmonized.apply(lambda row: croc_id_map[row['mrn']] if row['mrn'] in croc_id_map else '', axis=1)
    # improve_mrns = harmonized.loc[harmonized['study'] == 'IMPROVE', ['mrn', 'record_id']]
    # improve_id_map = dict(zip(improve_mrns['mrn'], improve_mrns['record_id']))
    # harmonized['improve_id'] = harmonized.apply(lambda row: improve_id_map[row['mrn']] if row['mrn'] in improve_id_map else '', axis=1)
    # penguin_mrns = harmonized.loc[harmonized['study'] == 'PENGUIN', ['mrn', 'record_id']]
    # penguin_id_map = dict(zip(penguin_mrns['mrn'], penguin_mrns['record_id']))
    # harmonized['penguin_id'] = harmonized.apply(lambda row: penguin_id_map[row['mrn']] if row['mrn'] in penguin_id_map else '', axis=1)
    # rh_mrns = harmonized.loc[harmonized['study'] == 'RENAL-HEIR', ['mrn', 'record_id']]
    # rh_id_map = dict(zip(rh_mrns['mrn'], rh_mrns['record_id']))
    # harmonized['rh_id'] = harmonized.apply(lambda row: rh_id_map[row['mrn']] if row['mrn'] in rh_id_map else '', axis=1)
    # rh2_mrns = harmonized.loc[harmonized['study'] == 'RENAL-HEIRitage', ['mrn', 'record_id']]
    # rh2_id_map = dict(zip(rh2_mrns['mrn'], rh2_mrns['record_id']))
    # harmonized['rh2_id'] = harmonized.apply(lambda row: rh2_id_map[row['mrn']] if row['mrn'] in rh2_id_map else '', axis=1)
    # panther_mrns = harmonized.loc[harmonized['study'] == 'PANTHER', ['mrn', 'record_id']]
    # panther_id_map = dict(zip(panther_mrns['mrn'], panther_mrns['record_id']))
    # harmonized['panther_id'] = harmonized.apply(lambda row: panther_id_map[row['mrn']] if row['mrn'] in panther_id_map else '', axis=1)
    # panda_mrns = harmonized.loc[harmonized['study'] == 'PANDA', ['mrn', 'record_id']]
    # panda_id_map = dict(zip(panda_mrns['mrn'], panda_mrns['record_id']))
    # harmonized['panda_id'] = harmonized.apply(lambda row: panda_id_map[row['mrn']] if row['mrn'] in panda_id_map else '', axis=1)
    # attempt_mrns = harmonized.loc[harmonized['study'] == 'ATTEMPT', ['mrn', 'record_id']]
    # attempt_id_map = dict(zip(attempt_mrns['mrn'], attempt_mrns['record_id']))
    # harmonized['attempt_id'] = harmonized.apply(lambda row: attempt_id_map[row['mrn']] if row['mrn'] in attempt_id_map else '', axis=1)
    
    create_study_id_columns(harmonized)

    # Sort columns
    id_cols = ["record_id", "casper_id", "coffee_id", "croc_id", "improve_id", 
                "penguin_id", "rh_id", "rh2_id", "panther_id", "panda_id", "attempt_id", "rpc2_id", "swth_id",
                "mrn", "co_enroll_id", "study", "dob", "diabetes_dx_date",
               "sex", "race", "ethnicity", "visit", "procedure", "date", "group"]
    other_cols = harmonized.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    harmonized = harmonized[id_cols + other_cols]

    # Sort rows
    harmonized.sort_values(
        ["study", "record_id", "visit", "procedure", "date"], inplace=True)
    # Format dates nicely
    harmonized[dates] = harmonized[dates].apply(
        lambda x: x.dt.strftime('%Y-%m-%d'))
    # Return
    harmonized = harmonized.astype(object)
    tocsv_path = base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv"
    dictionary.to_csv(tocsv_path, index=False)
    return harmonized
