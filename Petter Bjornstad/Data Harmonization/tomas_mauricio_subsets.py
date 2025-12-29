import pandas as pd
import numpy as np
import sys
import os

# Add the path to import the RPC2 cleaning function
sys.path.insert(0, os.path.expanduser('~') + "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from rpc2 import clean_rpc2_redcap

# Pull RPC2 data directly from REDCap (before harmonization loses the visits)
rpc2_raw = clean_rpc2_redcap()

print("=== DIAGNOSTIC: Raw RPC2 data from REDCap ===")
print(f"Total rows: {len(rpc2_raw)}")
print(f"Visits:")
print(rpc2_raw['visit'].value_counts(dropna=False))
print()

# Select fields of interest - adjust based on what's available
fields_wanted = ["record_id", "dob", "date", "sex", "age", "group", "visit", "procedure", 
                 "creatinine_s", "bmi", "height", "weight", "acr_u", "screen_urine_acr",
                 "microalbumin_u", "creatinine_u", "hba1c", "cystatin_c_s", "bun"]

# Keep only columns that exist in the data
available_fields = [f for f in fields_wanted if f in rpc2_raw.columns]
missing_fields = [f for f in fields_wanted if f not in rpc2_raw.columns]
if missing_fields:
    print(f"Note: These fields are not in the data: {missing_fields}")

rpc2 = rpc2_raw[available_fields].copy()

# Calculate acr_u if not present but components are available
if 'acr_u' not in rpc2.columns or rpc2['acr_u'].isna().all():
    if 'microalbumin_u' in rpc2.columns and 'creatinine_u' in rpc2.columns:
        rpc2['acr_u'] = (pd.to_numeric(rpc2['microalbumin_u'], errors='coerce') * 100 / 
                        pd.to_numeric(rpc2['creatinine_u'], errors='coerce'))
        print("Calculated acr_u from microalbumin_u and creatinine_u")

# Remove study failed participants
study_failed = ["RPC-10", "RPC-12", "RPC-14", "RPC-19", "RPC-22", "RPC-23", "RPC-24", 
                "RPC-30", "RPC-32", "RPC-33", "RPC-34", "RPC-35", "RPC-36", "RPC-39", 
                "RPC-40", "RPC-42"]
rpc2 = rpc2[~rpc2['record_id'].isin(study_failed)]

print(f"After removing failed participants: {len(rpc2)} rows")
print(f"Visits after filtering:")
print(rpc2['visit'].value_counts(dropna=False))
print()

# Drop procedure column if present (will aggregate across procedures)
if 'procedure' in rpc2.columns:
    rpc2 = rpc2.drop(columns=['procedure'])

# Drop helper columns used for acr_u calculation
cols_to_drop = ['microalbumin_u', 'creatinine_u', 'screen_urine_acr']
rpc2 = rpc2.drop(columns=[c for c in cols_to_drop if c in rpc2.columns], errors='ignore')

# Convert date columns
if 'date' in rpc2.columns:
    rpc2['date'] = pd.to_datetime(rpc2['date'], errors='coerce')
if 'dob' in rpc2.columns:
    rpc2['dob'] = pd.to_datetime(rpc2['dob'], errors='coerce')

# Convert age to numeric
if 'age' in rpc2.columns:
    rpc2['age'] = pd.to_numeric(rpc2['age'], errors='coerce')
    
    # Calculate age from date of birth where missing
    if 'dob' in rpc2.columns and 'date' in rpc2.columns:
        calculated_age = round((rpc2["date"] - rpc2["dob"]).dt.days / 365.25, 2)
        rpc2['age'] = rpc2['age'].fillna(calculated_age)

# Create age_screening column - age at screening visit (before aggregation)
if 'age' in rpc2.columns:
    screening_ages = rpc2[rpc2['visit'] == 'screening'].groupby('record_id')['age'].first()
    rpc2['age_screening'] = rpc2['record_id'].map(screening_ages)

# Define categorical and numeric columns based on what's available
all_categorical = ['dob', 'date', 'sex', 'group']
all_numeric = ['age', 'age_screening', 'creatinine_s', 'bmi', 'height', 'weight', 
               'acr_u', 'hba1c', 'cystatin_c_s', 'bun']

categorical_cols = [c for c in all_categorical if c in rpc2.columns]
numeric_cols = [c for c in all_numeric if c in rpc2.columns]

# Convert numeric columns to actual numeric types (they may be strings/objects)
for col in numeric_cols:
    rpc2[col] = pd.to_numeric(rpc2[col], errors='coerce')

# Create aggregation dictionary: first for categorical, mean for numeric
agg_dict = {col: 'first' for col in categorical_cols}
agg_dict.update({col: 'mean' for col in numeric_cols})

# Group by record_id and visit - one row per pair
rpc2_filtered = rpc2.groupby(['record_id', 'visit'], as_index=False).agg(agg_dict)

# Calculate eGFR (before converting age to Int64 to avoid type issues)
def calc_egfr(df, age="age", age_screening="age_screening", serum_creatinine="creatinine_s", 
              sex="sex", male="Male", female="Female"):
    data = df.copy()
    
    if serum_creatinine not in data.columns or sex not in data.columns:
        print("Warning: Cannot calculate eGFR - missing required columns")
        return data
    
    serum_creatinine_vals = pd.to_numeric(data[serum_creatinine], errors="coerce")
    sex_data = data[sex].copy()
    sex_data = sex_data.replace({male: "M", female: "F", "Other": np.nan, "": np.nan})
    
    age_data = pd.to_numeric(data[age], errors="coerce") if age in data.columns else pd.Series([np.nan]*len(data))
    
    # Use age_screening as fallback when age is missing
    if age_screening in data.columns:
        age_screening_data = pd.to_numeric(data[age_screening], errors="coerce")
        age_data = age_data.fillna(age_screening_data)
    
    f = sex_data.replace({"M": 0, "F": 1})
    a = sex_data.replace({"M": -0.302, "F": -0.241})
    k = sex_data.replace({"M": 0.9, "F": 0.7})
    
    eGFR_CKD_epi = 142 * (np.minimum(serum_creatinine_vals / k, 1)**a) * \
        (np.maximum(serum_creatinine_vals / k, 1)**-1.200) * \
        (0.9938**age_data) * (1.012 * f + (1 - f))
    
    data["eGFR_CKD_epi"] = eGFR_CKD_epi
    return data

rpc2_filtered = calc_egfr(rpc2_filtered, age="age", age_screening="age_screening", 
                          serum_creatinine="creatinine_s", sex="sex",
                          male="Male", female="Female")

# Round and convert age to proper data type (after eGFR calculation)
if 'age' in rpc2_filtered.columns:
    rpc2_filtered['age'] = rpc2_filtered['age'].round().astype('Int64')

# Drop rows where ALL of the specified clinical variables are empty
clinical_vars = ['creatinine_s', 'bmi', 'height', 'weight', 'acr_u', 'hba1c', 'cystatin_c_s', 'bun', 'eGFR_CKD_epi']
clinical_vars_present = [c for c in clinical_vars if c in rpc2_filtered.columns]
rpc2_filtered = rpc2_filtered[rpc2_filtered[clinical_vars_present].notna().any(axis=1)].copy()

print(f"After dropping rows with all clinical vars empty: {len(rpc2_filtered)} rows")

# Define visit order and sort
visit_order = ["v1_screening_arm_1",  "v2_gfr_mri_arm_1", 
                                   "v3_arm_1", "v4_arm_1",
                                   "p5_phone_visit_arm_1",
                                   "v61_med_dispense_arm_1",
                                   "v62_med_dispense_arm_1",
                                   "v7_gfr_mri_arm_1",
                                   "v8_arm_1"]
rpc2_filtered['visit'] = pd.Categorical(rpc2_filtered['visit'], categories=visit_order, ordered=True)
rpc2_filtered = rpc2_filtered.sort_values(['record_id', 'visit']).reset_index(drop=True)

# Remove dob and age_screening columns if present
cols_to_drop = [c for c in ['dob', 'age_screening'] if c in rpc2_filtered.columns]
if cols_to_drop:
    rpc2_filtered = rpc2_filtered.drop(columns=cols_to_drop)

# Final output
print("=== FINAL DATASET ===")
print(f"Total rows: {len(rpc2_filtered)}")
print(f"Unique record_ids: {rpc2_filtered['record_id'].nunique()}")
print(f"\nVisit counts:")
print(rpc2_filtered['visit'].value_counts())
print(f"\nRows per record_id:")
print(rpc2_filtered.groupby('record_id').size().value_counts().sort_index())

print(rpc2_filtered)

# Save
output_path = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/rpc2_update.csv"
rpc2_filtered.to_csv(output_path, index=False)
print(f"\nSaved to: {output_path}")