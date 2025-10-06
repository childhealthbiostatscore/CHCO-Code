import pandas as pd
import numpy as np

harmonized = pd.read_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv")

fields = ["record_id", "dob", "date", "sex", "age", "group", "visit", "procedure", "creatinine_s", "bmi", "height", "acr_u", "hba1c", "cystatin_c_s", "bun"]

# Initial filtering - ONLY by study and group
rpc2 = harmonized[harmonized['study'] == "RPC2"][fields].drop_duplicates().copy()

# Store information about which rows originally had date values
rpc2['had_original_date'] = rpc2['date'].notna()

# Make sure record_id is never marked as duplicate (even if there are multiple)
cols = rpc2.columns
duplicates = cols.duplicated()
duplicates = duplicates & (cols != "record_id")

rpc2 = rpc2.loc[:, ~duplicates]

rpc2['group'] = rpc2.groupby('record_id')['group'].transform('first')
rpc2['sex'] = rpc2.groupby('record_id')['sex'].transform('first')
rpc2['dob'] = rpc2.groupby('record_id')['dob'].transform('first')

mean_heights = rpc2.groupby(['record_id', 'visit'])['height'].transform('mean')
rpc2['height'] = rpc2['height'].fillna(mean_heights)
mean_bmi = rpc2.groupby(['record_id', 'visit'])['bmi'].transform('mean')
rpc2['bmi'] = rpc2['bmi'].fillna(mean_bmi)

# Convert date columns BEFORE filling missing dates
rpc2['date'] = pd.to_datetime(rpc2['date'])
rpc2['dob'] = pd.to_datetime(rpc2['dob'])

# NEW: Fill missing dates from adjacent rows with same visit or procedure
# Sort by record_id to ensure we're looking at the right adjacent rows
rpc2 = rpc2.sort_index().reset_index(drop=True)

# To handle potential non-numeric values, first coerce them to NaN.
rpc2['age'] = pd.to_numeric(rpc2['age'], errors='coerce')

# Calculate age from date of birth using your formula
calculated_age = round((rpc2["date"] - rpc2["dob"]).dt.days / 365.25, 2)

# Fill missing age values with calculated age (preserves existing non-null values)
rpc2['age'] = rpc2['age'].fillna(calculated_age)

# Round and convert to proper data type
rpc2['age'] = rpc2['age'].round().astype('Int64')

# NEW: Create age_screening column - age at screening visit (lowercase 'screening')
# Get age at screening for each record_id (where visit == 'screening')
screening_ages = rpc2[rpc2['visit'] == 'screening'].groupby('record_id')['age'].first()
rpc2['age_screening'] = rpc2['record_id'].map(screening_ages)

def calc_egfr(df, age="age", age_screening="age_screening", serum_creatinine="creatinine_s", 
              sex="sex", male="Male", female="Female"):
    # Make a copy of the dataframe
    data = df.copy()
    # Format input
    serum_creatinine = pd.to_numeric(data[serum_creatinine], errors="coerce")
    sex_data = data[sex].copy()
    sex_data.replace({male: "M", female: "F", "Other": np.nan, "": np.nan}, inplace=True)
    age_data = pd.to_numeric(data[age], errors="coerce")
    
    # NEW: Use age_screening as fallback when age is missing
    age_screening_data = pd.to_numeric(data[age_screening], errors="coerce")
    age_data = age_data.fillna(age_screening_data)
    
    f = sex_data.replace({"M": 0, "F": 1})
    a = sex_data.replace({"M": -0.302, "F": -0.241})
    k = sex_data.replace({"M": 0.9, "F": 0.7})
    eGFR_CKD_epi = 142 * (np.minimum(serum_creatinine / k, 1)**a) * \
        (np.maximum(serum_creatinine / k, 1)**-1.200) * \
        (0.9938**age_data) * (1.012 * f + (1 - f))
    egfr = pd.concat([eGFR_CKD_epi], axis=1)
    egfr.columns = ["eGFR_CKD_epi"]
    data = pd.concat([data, egfr], axis=1)
    return data

rpc2 = calc_egfr(rpc2, age="age", age_screening="age_screening", 
                 serum_creatinine="creatinine_s", sex="sex",
                 male="Male", female="Female")

# Modified filtering logic: Keep all rows that originally had dates, 
# OR rows that have at least one of the clinical variables
clinical_vars = ['creatinine_s','acr_u', 'hba1c', 'cystatin_c_s', 'bun', 'eGFR_CKD_epi']

# Create mask for rows to keep:
# 1. Rows that originally had date values, OR
# 2. Rows that have at least one non-null clinical variable
keep_mask = (
    rpc2['had_original_date'] | 
    rpc2[clinical_vars].notna().any(axis=1)
)

# Apply the filter
rpc2_filtered = rpc2[keep_mask].copy()

# Drop the helper column
rpc2_filtered = rpc2_filtered.drop('had_original_date', axis=1)

print(f"Original RPC2 rows: {len(rpc2)}")
print(f"Rows with original dates: {rpc2['had_original_date'].sum()}")
print(f"Rows with dates after filling from adjacent rows: {rpc2['date'].notna().sum()}")
print(f"After filtering (keeping date rows + clinical data rows): {len(rpc2_filtered)} rows")

print(rpc2_filtered)

rpc2_filtered.to_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/tomas_m_rpc2_subset.csv", index=False)