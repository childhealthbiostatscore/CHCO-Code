import pandas as pd
import numpy as np

# --- CONFIGURATION ---
INPUT_FILE = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"
OUTPUT_FILE = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/cgm_data_inventory.csv"

ID_COLS = ['record_id', 'visit']

# Numeric CGM data columns (excludes cgm_yn which is a yes/no flag).
# Only PANTHER collects these via procedure='cgm'; other studies only have cgm_yn
# in their medications rows and are excluded by the procedure filter.
CGM_COLS = [
    'cgm_avg_glucose',
    'cgm_glucose_var',
    'cgm_high',
    'cgm_low',
    'cgm_target',
    'cgm_v_high',
    'cgm_v_low',
]

# --- LOAD ---
harmonized = pd.read_csv(INPUT_FILE, low_memory=False)
harmonized.replace('NA', np.nan, inplace=True)

existing_cgm_cols = [c for c in CGM_COLS if c in harmonized.columns]
missing_cgm_cols = sorted(set(CGM_COLS) - set(existing_cgm_cols))
if missing_cgm_cols:
    print(f"Warning: CGM columns not found in dataset ({len(missing_cgm_cols)}): {missing_cgm_cols}")

# --- FILTER TO CGM PROCEDURE ROWS ONLY ---
# Only studies with procedure='cgm' have meaningful CGM numeric data
cols_to_load = [c for c in ID_COLS + ['study', 'procedure', 'cgm_yn'] + existing_cgm_cols if c in harmonized.columns]
df = harmonized[cols_to_load].copy()
df = df[(df['procedure'] == 'cgm') & (df['visit'] != 'screening')].copy()

# --- COLLAPSE TO ONE ROW PER record_id x visit ---
# Numeric CGM cols: take mean; cgm_yn/study: take first value
study_lookup = df.groupby(ID_COLS)[['study', 'cgm_yn']].first().reset_index()
agg_rules = {c: 'mean' for c in existing_cgm_cols}
df_condensed = df.groupby(ID_COLS, as_index=False).agg(agg_rules)
df_condensed = df_condensed.merge(study_lookup, on=ID_COLS, how='left')

total = len(existing_cgm_cols)

# --- BUILD INVENTORY ---
def compute_status(row):
    if row['cgm_yn'] == 'No':
        return 'cgm not worn'
    values = row[existing_cgm_cols]
    n_missing = values.isna().sum()
    if n_missing == 0:
        return 'complete'
    elif n_missing == total:
        return 'all missing'
    else:
        return f'{n_missing} of {total} missing'

# data_in_redcap = TRUE if cgm_yn was filled in (someone recorded whether CGM was worn)
df_condensed['data_in_redcap'] = df_condensed['cgm_yn'].notna()
df_condensed['status'] = df_condensed.apply(compute_status, axis=1)

inventory = df_condensed[['record_id', 'study', 'visit', 'data_in_redcap', 'status']].copy()

# --- SAVE ---
inventory.to_csv(OUTPUT_FILE, index=False)
print(f"Saved CGM inventory to: {OUTPUT_FILE}")

# --- SUMMARY ---
print(f"\n=== SUMMARY ===")
print(f"Total record_id × visit rows: {len(inventory)}")
print(f"Unique record_ids: {inventory['record_id'].nunique()}")
print(f"\nStatus breakdown:")
print(inventory['status'].value_counts().to_string())
print(f"\ndata_in_redcap breakdown:")
print(inventory['data_in_redcap'].value_counts().to_string())
