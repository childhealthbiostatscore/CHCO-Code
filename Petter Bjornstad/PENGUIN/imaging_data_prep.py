import pandas as pd
import numpy as np

# --- CONFIGURATION ---
INPUT_FILE = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"
OUTPUT_FILE = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/PENGUIN/Data_Cleaned/imaging_data.csv"

# Define all columns to extract from harmonized dataset
ID_COLS = ['record_id', 'visit']

CATEGORICAL_COLS = [
    'study', 'procedure', 'group', 'sex', 'adpkd_classification'
]

NUMERIC_COLS = [
    'age', 'bmi',
    'eGFR_CKD_epi', 'eGFR_fas_cr', 'eGFR_fas_cr_cysc', 'acr_u',
    'gfr_raw_plasma', 'gfr_bsa_plasma', 'erpf_raw_plasma', 
    'erpf_bsa_plasma', 'glomerular_pressure', 'ra', 're',

    # PET
    'rc_f','rc_k2', 'rc_vb', 'rc_k1', 
    'rm_f', 'rm_k2', 'rm_vb', 'rm_k1', 
    'lc_f', 'lc_k2', 'lc_vb', 'lc_k1', 
    'lm_f', 'lm_k2', 'lm_vb', 'lm_k1',

    'avg_c_f',	'avg_c_k1',	'avg_c_k2',
    'avg_m_f', 'avg_m_k1', 'avg_m_k2',

    # fMRI
    'total_cyst_volume_ml',
    'total_number_of_cysts',
    'cyst_parenchyma_sa',

    # voxelwise
    'rtot_k1_w_cyst', 'rtot_k1_wo_cyst',
    'ltot_k1_w_cyst', 'ltot_k1_wo_cyst',
    'rtot_k2_w_cyst', 'rtot_k2_wo_cyst',
    'ltot_k2_w_cyst', 'ltot_k2_wo_cyst',
    'rtot_k1_w_cyst_percent_avg', 'rtot_k1_wo_cyst_percent_avg',
    'ltot_k1_w_cyst_percent_avg', 'ltot_k1_wo_cyst_percent_avg',
    'rtot_k2_w_cyst_percent_avg', 'rtot_k2_wo_cyst_percent_avg',
    'ltot_k2_w_cyst_percent_avg', 'ltot_k2_wo_cyst_percent_avg',
    'rc_k1_w_cyst_vw', 'rc_k1_wo_cyst_vw',
    'lc_k1_w_cyst_vw', 'lc_k1_wo_cyst_vw',
    'rc_k1_w_cyst_percent_avg', 'rc_k1_wo_cyst_percent_avg',
    'lc_k1_w_cyst_percent_avg', 'lc_k1_wo_cyst_percent_avg',
    'rc_k2_w_cyst_percent_avg', 'rc_k2_wo_cyst_percent_avg',
    'lc_k2_w_cyst_percent_avg', 'lc_k2_wo_cyst_percent_avg',

    'avg_c_k2_wo_cyst_vw', 'avg_m_k2_wo_cyst_vw',

    'ht_adj_tkv'
]


# Define columns that should be filled with within-visit means
# (exclude identifiers and categorical variables)
COLUMNS_TO_FILL = [
    'eGFR_CKD_epi', 'eGFR_fas_cr', 'eGFR_fas_cr_cysc', 'acr_u',
    'gfr_raw_plasma', 'gfr_bsa_plasma', 'erpf_raw_plasma', 
    'erpf_bsa_plasma', 'glomerular_pressure', 'ra', 're',
    #PET
    'rc_f','rc_k2', 'rc_vb', 'rc_k1',
    'rm_f', 'rm_k2', 'rm_vb', 'rm_k1',
    'lc_f', 'lc_k2', 'lc_vb', 'lc_k1',
    'lm_f', 'lm_k2', 'lm_vb', 'lm_k1',

    #fMRI
    'total_cyst_volume_ml', 
    'total_number_of_cysts', 
    'adpkd_classification', 
    'cyst_parenchyma_sa',
    #Voxelwise
    'rtot_k1_w_cyst', 'rtot_k1_wo_cyst', 
    'ltot_k1_w_cyst', 'ltot_k1_wo_cyst',

    'rtot_k2_w_cyst', 'rtot_k2_wo_cyst',
    'ltot_k2_w_cyst', 'ltot_k2_wo_cyst',

    'rtot_k1_w_cyst_percent_avg', 'rtot_k1_wo_cyst_percent_avg',
    'ltot_k1_w_cyst_percent_avg', 'ltot_k1_wo_cyst_percent_avg',

    'rtot_k2_w_cyst_percent_avg', 'rtot_k2_wo_cyst_percent_avg',
    'ltot_k2_w_cyst_percent_avg', 'ltot_k2_wo_cyst_percent_avg',

    'rc_k1_w_cyst_vw', 'rc_k1_wo_cyst_vw',
    'lc_k1_w_cyst_vw', 'lc_k1_wo_cyst_vw',

    'rc_k1_w_cyst_percent_avg', 'rc_k1_wo_cyst_percent_avg',
    'lc_k1_w_cyst_percent_avg', 'lc_k1_wo_cyst_percent_avg',

    'rc_k2_w_cyst_percent_avg', 'rc_k2_wo_cyst_percent_avg',
    'lc_k2_w_cyst_percent_avg', 'lc_k2_wo_cyst_percent_avg',

    'ht_adj_tkv'
]

# Categorical/identifier columns to forward-fill within record_id
DEMOGRAPHIC_COLUMNS = ['group', 'sex']

columns = ID_COLS + CATEGORICAL_COLS + NUMERIC_COLS

import pandas as pd
import numpy as np
# Read the harmonized dataset (suppress dtype warning)
harmonized = pd.read_csv("Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", low_memory=False)
# Define studies and columns to keep
studies = ['PENGUIN', 'CROCODILE']
groups = ['Lean Control', 'PKD']

# Filter for selected studies and columns
df = harmonized.loc[harmonized['study'].isin(studies), columns].copy()
df = df[df['group'].isin(groups)].copy()
# Replace 'NA' strings with NaN
df.replace('NA', np.nan, inplace=True)


print("=== ORIGINAL DATA ===")
print(f"Total rows: {len(df)}")
print(f"Unique visits: {df['visit'].unique()}")
print(f"Baseline rows: {(df['visit'] == 'baseline').sum()}")
# Define columns to conserve from all procedure types
conserve_cols = COLUMNS_TO_FILL
# First, fill age from any row with same record_id (age shouldn't vary)

df['age'] = df['age'].fillna(df.groupby('record_id')['age'].transform('first'))

    
# For each combination of record_id and visit, calculate mean of conserve_cols
# from all procedure types, then fill NaN values with the group mean
for col in conserve_cols:
# Calculate mean for each record_id + visit combination across all procedures
    group_means = df.groupby(['record_id', 'visit'])[col].transform('mean')
# Fill NaN values with the group mean
    df[col] = df[col].fillna(group_means)
# Apply transformations for bmi
df['bmi'] = df['bmi'].fillna(
df.groupby(['record_id', 'visit'])['bmi'].transform('first')
)
# Filter to keep only baseline visits
df_baseline = df[df['visit'] == 'baseline'].copy()
# Drop columns that are entirely NaN and remove duplicates
df_baseline.dropna(how='all', axis=1, inplace=True)
df_baseline.drop_duplicates(inplace=True)

df_baseline.drop(columns=['visit'], inplace=True)
# Display summary statistics
print(f"\n=== FINAL SUMMARY ===")
print(f"Final baseline rows: {len(df_baseline)}")
print(f"Unique procedures in baseline: {df_baseline['procedure'].unique()}")
print(f"\nRows per procedure type:")
print(df_baseline['procedure'].value_counts())
print(f"\nMissing values per column:")
missing_summary = df_baseline.isnull().sum()
print(missing_summary)


improve_mask = df_baseline['study'].isin(studies)
cardio_mri_mask = df_baseline['procedure'] == 'pet_scan'
df_baseline = df_baseline[~improve_mask | cardio_mri_mask].copy()
print(f"Final total rows: {len(df_baseline)}")



# Save to CSV
df_baseline.to_csv(OUTPUT_FILE, index=False)