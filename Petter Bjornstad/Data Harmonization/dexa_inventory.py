import pandas as pd
import numpy as np

# --- CONFIGURATION ---
INPUT_FILE = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"

ID_COLS = ['record_id', 'visit','procedure']

CATEGORICAL_COLS = [
    'study', 'group', 'sex'
]

NUMERIC_COLS = [
    'age', 'bmi',
    'bod_pod_body_fat',
    'bod_pod_est_rmr',
    'bod_pod_fat_kg',
    'bod_pod_lean_kg',
    'bod_pod_lean_mass',
    'dexa_ag_ratio',
    'dexa_body_fat',
    'dexa_bone_mineral_density',
    'dexa_est_vat',
    'dexa_fat_kg',
    'dexa_lean_kg',
    'dexa_lean_mass',
    'dexa_trunk_kg',
    'dexa_trunk_mass',
    'height',
    'weight',
    'bia_complete',
    'bia_body_fat',
    'bia_lean_mass',
    'bia_fat_mass_kg',
    'bia_lean_mass_kg',
    'tanita_scale_complete',
    'tanita_body_fat',
    'tanita_lean_mass',
    'tanita_fat_mass_kg',
    'tanita_lean_mass_kg',
]

# Read the dataset
harmonized = pd.read_csv(INPUT_FILE, low_memory=False)
existing_cols = [c for c in ID_COLS + CATEGORICAL_COLS + NUMERIC_COLS if c in harmonized.columns]
missing_cols = sorted(set(ID_COLS + CATEGORICAL_COLS + NUMERIC_COLS) - set(existing_cols))
print(f"Missing columns ({len(missing_cols)}):")
print(missing_cols)

df = harmonized[existing_cols].copy()
df.replace('NA', np.nan, inplace=True)

# Fill numeric columns with mean per record_id + visit
for col in NUMERIC_COLS:
    if col in df.columns:
        df[col] = df.groupby(['record_id', 'visit'])[col].transform('mean')

# Fill categorical columns with first value per record_id + visit
for col in CATEGORICAL_COLS:
    if col in df.columns:
        df[col] = df.groupby(['record_id', 'visit'])[col].transform('first')

# Drop irrelevant studies
irrelevant_studies = ['ATTEMPT', 'COFFEE', 'RENAL-HEIRitage', 'RPC2', 'SWEETHEART', 'ULTRA']
if 'study' in df.columns:
    df = df[~df['study'].isin(irrelevant_studies)].copy()

# Collapse to ONE row per record_id × visit (numeric=mean, categorical=first)
numeric_existing = [c for c in NUMERIC_COLS if c in df.columns]
categorical_existing = [c for c in CATEGORICAL_COLS if c in df.columns]

agg_rules = {c: 'mean' for c in numeric_existing}
agg_rules.update({c: 'first' for c in categorical_existing})

df_condensed = df.groupby(ID_COLS, as_index=False).agg(agg_rules)

# Drop columns that are entirely NaN (after collapse)
df_condensed.dropna(how='all', axis=1, inplace=True)

# Save the condensed subset
out_file_subset = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/dexa_all_subset.csv"
df_condensed.to_csv(out_file_subset, index=False)
print(f"Saved subsetted dataset to: {out_file_subset}")

# Count numeric validity per row: non-NaN and non-zero
valid_numeric_count = (df_condensed[numeric_existing].notna() & (df_condensed[numeric_existing] != 0)).sum(axis=1)
df_condensed = df_condensed.loc[valid_numeric_count >= 2]

# Build summary: number of unique record_id per study × group × visit
all_results = []
for visit, df_visit in df_condensed.groupby('visit'):
    for (study, group), subset in df_visit.groupby(['study', 'group']):
        num_records = subset['record_id'].nunique()
        all_results.append({
            'study': study,
            'visit': visit,
            'group': group,
            'num_records': num_records
        })

# Convert summary to DataFrame and save
summary_df = pd.DataFrame(all_results)
out_file_summary = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/bodycomp_inventory_all_visits.csv"
summary_df.to_csv(out_file_summary, index=False)
print(f"Saved combined inventory to: {out_file_summary}")

print("\n=== FINAL SUMMARY ===")
print("Number of unique record_id with at least 3 valid numeric body composition measures:")
print(df_condensed['record_id'].nunique())
# Count numeric validity per row: non-NaN and non-zero
valid_numeric_count = (df_condensed[numeric_existing].notna() & (df_condensed[numeric_existing] != 0)).sum(axis=1)
df_condensed['valid_numeric_count'] = valid_numeric_count

# Split record_ids into two lists per study × group
for study in df_condensed['study'].dropna().unique():
    for group in df_condensed['group'].dropna().unique():
        subset = df_condensed[(df_condensed['study'] == study) & (df_condensed['group'] == group)]
        record_ids_minimal = subset.loc[subset['valid_numeric_count'] <= 2, 'record_id'].unique()
        record_ids_sufficient = subset.loc[subset['valid_numeric_count'] > 2, 'record_id'].unique()

        print(f"\nStudy: {study}, Group: {group}")
        print(f"Record_ids with more than 2 valid measures ({len(record_ids_sufficient)}):")
        print(record_ids_sufficient.tolist())
        print(f"Record_ids with 2 or fewer valid measures ({len(record_ids_minimal)}):")
        print(record_ids_minimal.tolist())


