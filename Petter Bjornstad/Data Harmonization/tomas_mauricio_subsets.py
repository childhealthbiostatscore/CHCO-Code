import pandas as pd
import numpy as np

harmonized = pd.read_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv")

fields = ["record_id", "dob", "date", "sex", "age", "group", "visit", "procedure", "creatinine_s", "bmi", "height", "weight", "acr_u", "hba1c", "cystatin_c_s", "bun"]

# Initial filtering - ONLY by study and group
rpc2 = harmonized[harmonized['study'] == "RPC2"][fields].drop_duplicates().copy()

# Make sure record_id is never marked as duplicate (even if there are multiple)
cols = rpc2.columns
duplicates = cols.duplicated()
duplicates = duplicates & (cols != "record_id")
rpc2 = rpc2.loc[:, ~duplicates]

rpc2['group'] = rpc2.groupby('record_id')['group'].transform('first')
rpc2['sex'] = rpc2.groupby('record_id')['sex'].transform('first')
rpc2['dob'] = rpc2.groupby('record_id')['dob'].transform('first')

# Fill missing values with mean by record_id and visit
mean_heights = rpc2.groupby(['record_id', 'visit'])['height'].transform('mean')
rpc2['height'] = rpc2['height'].fillna(mean_heights)
mean_weight = rpc2.groupby(['record_id', 'visit'])['weight'].transform('mean')
rpc2['weight'] = rpc2['weight'].fillna(mean_weight)
mean_bmi = rpc2.groupby(['record_id', 'visit'])['bmi'].transform('mean')
rpc2['bmi'] = rpc2['bmi'].fillna(mean_bmi)
mean_creatinine = rpc2.groupby(['record_id', 'visit'])['creatinine_s'].transform('mean')
rpc2['creatinine_s'] = rpc2['creatinine_s'].fillna(mean_creatinine)
mean_uacr = rpc2.groupby(['record_id', 'visit'])['acr_u'].transform('mean')
rpc2['acr_u'] = rpc2['acr_u'].fillna(mean_uacr)
mean_hba1c = rpc2.groupby(['record_id', 'visit'])['hba1c'].transform('mean')
rpc2['hba1c'] = rpc2['hba1c'].fillna(mean_hba1c)
mean_cystatin = rpc2.groupby(['record_id', 'visit'])['cystatin_c_s'].transform('mean')
rpc2['cystatin_c_s'] = rpc2['cystatin_c_s'].fillna(mean_cystatin)
mean_bun = rpc2.groupby(['record_id', 'visit'])['bun'].transform('mean')
rpc2['bun'] = rpc2['bun'].fillna(mean_bun)

# Convert date columns BEFORE filling missing dates
rpc2['date'] = pd.to_datetime(rpc2['date'])
rpc2['dob'] = pd.to_datetime(rpc2['dob'])

# Fill missing dates with first date value for same record_id and visit
rpc2['date'] = rpc2.groupby(['record_id', 'visit'])['date'].transform('first')

# Sort by record_id to ensure we're looking at the right adjacent rows
rpc2 = rpc2.sort_index().reset_index(drop=True)

# To handle potential non-numeric values, first coerce them to NaN
rpc2['age'] = pd.to_numeric(rpc2['age'], errors='coerce')

# Calculate age from date of birth
calculated_age = round((rpc2["date"] - rpc2["dob"]).dt.days / 365.25, 2)

# Fill missing age values with calculated age
rpc2['age'] = rpc2['age'].fillna(calculated_age)

# Round and convert to proper data type
rpc2['age'] = rpc2['age'].round().astype('Int64')

# Create age_screening column - age at screening visit
screening_ages = rpc2[rpc2['visit'] == 'screening'].groupby('record_id')['age'].first()
rpc2['age_screening'] = rpc2['record_id'].map(screening_ages)

def calc_egfr(df, age="age", age_screening="age_screening", serum_creatinine="creatinine_s", 
              sex="sex", male="Male", female="Female"):
    data = df.copy()
    serum_creatinine = pd.to_numeric(data[serum_creatinine], errors="coerce")
    sex_data = data[sex].copy()
    sex_data.replace({male: "M", female: "F", "Other": np.nan, "": np.nan}, inplace=True)
    age_data = pd.to_numeric(data[age], errors="coerce")
    
    # Use age_screening as fallback when age is missing
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

# Drop rows that are missing ALL of the specified clinical variables
clinical_vars = ['creatinine_s', 'bmi', 'height', 'weight', 'acr_u', 'hba1c', 'cystatin_c_s', 'bun']
rpc2_filtered = rpc2[rpc2[clinical_vars].notna().any(axis=1)].copy()

# Remove study failed participants
study_failed = ['RPC-06','RPC-09','RPC-10','RPC-11', 'RPC-12', 'RPC-14', 'RPC-15', 'RPC-18','RPC-19', 'RPC-20','RPC-22', 'RPC-23', 'RPC-24', 'RPC-28','RPC-30', 'RPC-31',
                'RPC-32', 'RPC-33', 'RPC-34', 'RPC-35', 'RPC-36', 'RPC-39', 'RPC-40', 'RPC-42', 'RPC-43', 'RPC-45']
rpc2_filtered = rpc2_filtered[~rpc2_filtered['record_id'].isin(study_failed)]

# Drop procedure and age_screening columns before dropping duplicates
rpc2_filtered = rpc2_filtered.drop(columns=['procedure', 'age_screening'])

# Drop duplicates
rpc2_filtered.drop_duplicates(inplace=True)

# Group by record_id and visit, taking mean of numerical variables
# First, identify non-numeric columns to keep using 'first'
non_numeric_cols = ['record_id', 'visit', 'dob', 'date', 'sex', 'group']
numeric_cols = [col for col in rpc2_filtered.columns if col not in non_numeric_cols]

# Create aggregation dictionary
agg_dict = {col: 'first' for col in non_numeric_cols if col not in ['record_id', 'visit']}
agg_dict.update({col: 'mean' for col in numeric_cols})

# Group by record_id and visit
rpc2_filtered = rpc2_filtered.groupby(['record_id', 'visit'], as_index=False).agg(agg_dict)

# Remove dob column
rpc2_filtered = rpc2_filtered.drop(columns=['dob'])

print(f"Original RPC2 rows: {len(rpc2)}")
print(f"After filtering and grouping: {len(rpc2_filtered)} rows")
print(f"Rows with dates: {rpc2_filtered['date'].notna().sum()}")

print(rpc2_filtered)

rpc2_filtered.to_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/rpc2_update.csv", index=False)