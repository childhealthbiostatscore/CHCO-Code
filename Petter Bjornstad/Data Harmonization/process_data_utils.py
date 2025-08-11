import pandas as pd
import numpy as np

harmonized = pd.read_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv")  

def subset_harmonized(data, 
                      vars=None, 
                      study_subset=None, 
                      numeric_fn=np.mean, 
                      non_numeric_fn=lambda x: x.dropna().iloc[-1] if not x.dropna().empty else np.nan,
                      group_by_vars=["record_id"],
                      summarize=True):
    # Required columns to always keep
    necessary_cols = ["study"] + group_by_vars
    
    # Check variables exist
    if vars is not None:
        missing_vars = set(vars) - set(data.columns)
        if missing_vars:
            raise ValueError(f"These variables are missing from dataset: {', '.join(missing_vars)}")
        cols_to_select = list(dict.fromkeys(necessary_cols + vars))  # preserve order, unique
        data = data.loc[:, cols_to_select]
    else:
        data = data.copy()

    # Filter by study if requested
    if study_subset is not None:
        data = data[data["study"].isin(study_subset)]

    # Replace empty strings or whitespace-only strings with NaN for all columns
    data = data.applymap(lambda x: np.nan if (isinstance(x, str) and x.strip() == "") else x)

    # Impute numeric missing with numeric_fn
    numeric_cols = data.select_dtypes(include="number").columns
    for col in numeric_cols:
        missing_mask = data[col].isna()
        if missing_mask.any():
            impute_value = numeric_fn(data.loc[~missing_mask, col])
            data.loc[missing_mask, col] = impute_value

    if summarize:
        # Group and summarize
        def summarize_group(group):
            result = {}
            for col in group.columns:
                if col in group_by_vars:
                    result[col] = group.name if isinstance(group.name, str) else group.name[0]
                elif pd.api.types.is_numeric_dtype(group[col]):
                    result[col] = numeric_fn(group[col].dropna())
                else:
                    # non_numeric_fn on non-numeric, default last non-na
                    result[col] = non_numeric_fn(group[col])
            return pd.Series(result)

        summarized = data.groupby(group_by_vars).apply(summarize_group).reset_index(drop=True)

        # Always keep necessary cols (study, group_by_vars)
        # Drop columns that are all NA except necessary columns
        cols_to_keep = [c for c in summarized.columns if (c in necessary_cols) or (summarized[c].notna().any())]
        summarized = summarized.loc[:, cols_to_keep]

        print("NA counts per column BEFORE dropping all-NA columns:")
        print(summarized.isna().sum())

        return summarized

    else:
        # Drop columns all NA except necessary columns
        cols_to_keep = [c for c in data.columns if (c in necessary_cols) or (data[c].notna().any())]
        data = data.loc[:, cols_to_keep]
        return data



vars_needed = [
    "copeptin", "hba1c", "sbp", "dbp",
    "mm_si", "baseline_ffa", "eGFR_CKD_epi", "gfr_raw_plasma", "gfr_bsa_plasma",
    "group", "age", "sex", "bmi", "diabetes_duration"
]

studies_to_include = ["CASPER", "CROCODILE", "RENAL-HEIRitage"]

subset_data = subset_harmonized(
    data=harmonized,
    vars=vars_needed,
    study_subset=studies_to_include,
    numeric_fn=np.mean,  # you can use np.median, np.min, etc.
    non_numeric_fn=lambda x: x.dropna().iloc[-1] if not x.dropna().empty else np.nan,
    group_by_vars=["record_id"],
    summarize=True
)

print(subset_data.head())


subset_data.to_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/Daniel_Casillas_Subset.csv", index=False)
