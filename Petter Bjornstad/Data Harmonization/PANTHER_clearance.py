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
from pfas_data_merge import PFAS
import getpass


# REDCap variable names -> harmonized dataset names
REDCAP_TO_HARMONIZED = {
    "gfr_15mgmin": "gfr_raw_plasma",
    "gfrbsa":      "gfr_bsa_plasma",
    "erpf_pah_85": "erpf_raw_plasma",
    "erpfbsa":     "erpf_bsa_plasma",
}
REDCAP_KEY_VARS = list(REDCAP_TO_HARMONIZED.keys())


def get_paths_and_tokens():
    user = getpass.getuser()
    if user == "choiyej":
        base_data_path = "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/"
        git_path = "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
    elif user == "pylell":
        base_data_path = "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
        git_path = "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    elif user == "shivaniramesh":
        base_data_path = os.path.expanduser("~/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
        git_path = "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    else:
        sys.exit(f"Unknown user: please specify root path for this user. (Detected user: {user})")
    tokens = pd.read_csv(base_data_path + "/Data Harmonization/api_tokens.csv")
    return base_data_path, tokens


def panther_redcap_data(base_data_path, tokens):
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "PANTHER", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    meta = pd.DataFrame(proj.metadata)
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]
    var = list(dict.fromkeys(["record_id", "group", "phys_map", "dxa_height", "dxa_weight"] + REDCAP_KEY_VARS +
               [v for v in meta.loc[meta["form_name"] == "study_visit_renal_clearance_testing", "field_name"]]))
    rct = pd.DataFrame(proj.export_records(fields=var))
    rct.replace(rep, np.nan, inplace=True)
    rct["group"] = rct.groupby(["record_id"])["group"].ffill()
    rct = rct.loc[rct["redcap_event_name"] != "screening_arm_1"].copy()
    rct["redcap_event_name"].replace(
        {"screening_arm_1": "screening", "baseline_arm_1": "baseline", "year_1_arm_1": "year_1", "year_2_arm_1": "year_2",  "year_3_arm_1": "year_3",  "year_4_arm_1": "year_4"},
        inplace=True)
    rct = rct.rename(columns={"redcap_event_name": "visit"})

    baseline_include = [
        "PAN-01-T", "PAN-02-T", "PAN-04-O", "PAN-05-T", "PAN-08-T", "PAN-09-C",
        "PAN-11-C", "PAN-12-O", "PAN-13-C", "PAN-14-C", "PAN-15-C", "PAN-17-O",
        "PAN-20-T", "PAN-22-C", "PAN-23-C", "PAN-24-T", "PAN-25-O", "PAN-26-O",
        "PAN-27-C", "PAN-28-C", "PAN-34-C", "PAN-35-C", "PAN-36-O", "PAN-37-C",
        "PAN-38-O", "PAN-39-O", "PAN-40-C", "PAN-41-C", "PAN-42-O", "PAN-44-O",
        "PAN-45-O", "PAN-46-C", "PAN-47-C", "PAN-48-C", "PAN-49-O", "PAN-51-O",
        "PAN-52-C", "PAN-53-C", "PAN-54-T", "PAN-55-O", "PAN-56-C", "PAN-57-C",
        "PAN-59-O", "PAN-60-O", "PAN-61-O", "PAN-63-C", "PAN-64-C", "PAN-65-O",
        "PAN-66-C", "PAN-68-C", "PAN-71-O", "PAN-73-C", "PAN-74-C", "PAN-75-O",
        "PAN-76-O", "PAN-77-T", "PAN-78-O", "PAN-79-C", "PAN-80-O", "PAN-83-O",
        "PAN-84-O", "PAN-85-C", "PAN-86-C", "PAN-87-C", "PAN-88-O", "PAN-89-C",
        "PAN-90-C", "PAN-91-C", "PAN-92-O", "PAN-93-O", "PAN-94-O", "PAN-95-O",
        "PAN-96-O", "PAN-97-C", "PAN-98-O", "PAN-99-O", "PAN-100-T", "PAN-101-O",
        "PAN-102-O", "PAN-103-O", "PAN-104-T", "PAN-105-O", "PAN-106-O", "PAN-107-O",
        "PAN-108-O", "PAN-109-C", "PAN-110-C",
    ]

    year_1_include = [
        "PAN-01-T", "PAN-02-T", "PAN-04-O", "PAN-05-T", "PAN-08-T", "PAN-09-C",
        "PAN-11-C", "PAN-12-O", "PAN-13-C", "PAN-14-C", "PAN-15-C", "PAN-16-O",
        "PAN-17-O", "PAN-20-T", "PAN-22-C", "PAN-23-C", "PAN-24-T", "PAN-25-O",
        "PAN-26-O", "PAN-27-C", "PAN-28-C", "PAN-34-C", "PAN-35-C", "PAN-36-O",
        "PAN-37-C", "PAN-38-O", "PAN-39-O", "PAN-40-C", "PAN-41-C", "PAN-42-O",
        "PAN-44-O", "PAN-45-O", "PAN-46-C", "PAN-47-C", "PAN-48-C", "PAN-49-O",
        "PAN-51-O", "PAN-52-C", "PAN-53-C", "PAN-54-T", "PAN-55-O", "PAN-56-C",
        "PAN-57-C", "PAN-59-O", "PAN-60-O", "PAN-61-O", "PAN-63-C", "PAN-64-C",
        "PAN-65-O", "PAN-66-C", "PAN-68-C", "PAN-71-O", "PAN-73-C", "PAN-74-C",
        "PAN-75-O", "PAN-76-O", "PAN-77-T", "PAN-78-O", "PAN-79-C", "PAN-80-O",
        "PAN-83-O", "PAN-85-C", "PAN-86-C", "PAN-87-C", "PAN-88-O", "PAN-89-C",
        "PAN-90-C", "PAN-91-C", "PAN-92-O", "PAN-93-O", "PAN-94-O", "PAN-95-O",
        "PAN-99-O", "PAN-102-O", "PAN-103-O",
    ]

    year_2_include = [
        "PAN-02-T", "PAN-04-O", "PAN-05-T", "PAN-08-T", "PAN-09-C", "PAN-11-C",
        "PAN-12-O", "PAN-13-C", "PAN-14-C", "PAN-15-C", "PAN-17-O", "PAN-20-T",
        "PAN-22-C", "PAN-23-C", "PAN-24-T", "PAN-25-O", "PAN-26-O", "PAN-27-C",
        "PAN-28-C", "PAN-34-C", "PAN-35-C", "PAN-36-O", "PAN-37-C", "PAN-40-C",
        "PAN-41-C", "PAN-45-O", "PAN-46-C", "PAN-47-C", "PAN-48-C", "PAN-49-O",
        "PAN-52-C", "PAN-53-C", "PAN-54-T", "PAN-55-O", "PAN-56-C", "PAN-57-C",
    ]

    visits = ["baseline", "year_1", "year_2"]
    rct_filtered = rct[rct["visit"].isin(visits)].copy()

    rows = []
    for rid in rct["record_id"].unique():
        for v in visits:
            if v == "baseline" and rid not in baseline_include:
                continue
            if v == "year_1" and rid not in year_1_include:
                continue
            if v == "year_2" and rid not in year_2_include:
                continue
            subset = rct_filtered[(rct_filtered["record_id"] == rid) & (rct_filtered["visit"] == v)]
            if subset.empty:
                rows.append({"record_id": rid, "visit": v, **{h: np.nan for h in REDCAP_TO_HARMONIZED.values()}})
            else:
                row = subset.iloc[0]
                rows.append({"record_id": rid, "visit": v, **{h: row[r] for r, h in REDCAP_TO_HARMONIZED.items()}})

    return pd.DataFrame(rows)

def crocodile_redcap_data(base_data_path, tokens):
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "CROCODILE", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]
    dictionary = pd.read_csv(base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv")

    var = list(dict.fromkeys(["record_id", "group", "bl_tot_protein", "hct_210", "visit_map", "phys_map"] + REDCAP_KEY_VARS +
               [v for v in meta.loc[meta["form_name"] == "study_visit_renal_clearance_testing", "field_name"]]))
    rct = pd.DataFrame(proj.export_records(fields=var))
    rct.replace(rep, np.nan, inplace=True)
    rct["visit"] = "baseline"
    rct_filtered = rct[rct["visit"] == "baseline"]

    rows = []
    for rid in rct["record_id"].unique():
        subset = rct_filtered[(rct_filtered["record_id"] == rid) & (rct_filtered["visit"] == "baseline")]
        if subset.empty:
            rows.append({"record_id": rid, "visit": "baseline", **{h: np.nan for h in REDCAP_TO_HARMONIZED.values()}})
        else:
            row = subset.iloc[0]
            rows.append({"record_id": rid, "visit": "baseline", **{h: row[r] for r, h in REDCAP_TO_HARMONIZED.items()}})

    return pd.DataFrame(rows)

def penguin_redcap_data(base_data_path, tokens):
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "PENGUIN", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]
    dictionary = pd.read_csv(base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv")
    var = list(dict.fromkeys(["record_id", "group", "bl_tot_protein", "hct_210", "visit_map", "phys_map"] + REDCAP_KEY_VARS +
               [v for v in meta.loc[meta["form_name"] == "study_visit_renal_clearance_testing", "field_name"]]))
    rct = pd.DataFrame(proj.export_records(fields=var))
    rct.replace(rep, np.nan, inplace=True)
    rct["visit"] = "baseline"
    rct_filtered = rct[rct["visit"] == "baseline"]
    rows = []
    for rid in rct["record_id"].unique():
        subset = rct_filtered[(rct_filtered["record_id"] == rid) & (rct_filtered["visit"] == "baseline")]
        if subset.empty:
            rows.append({"record_id": rid, "visit": "baseline", **{h: np.nan for h in REDCAP_TO_HARMONIZED.values()}})
        else:
            row = subset.iloc[0]
            rows.append({"record_id": rid, "visit": "baseline", **{h: row[r] for r, h in REDCAP_TO_HARMONIZED.items()}})
    return pd.DataFrame(rows)



def rpc2_redcap_data(base_data_path, tokens):
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "RPC2", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
    rep = rep + [str(r) for r in rep] + [""]

    var = list(dict.fromkeys(["study_id"] + REDCAP_KEY_VARS +
               [v for v in meta.loc[meta["form_name"] == "study_visit_renal_clearance_testing", "field_name"]]))
    rct = pd.DataFrame(proj.export_records(fields=var))
    rct.replace(rep, np.nan, inplace=True)
    rct.rename({"subject_id": "record_id"}, axis=1, inplace=True)

    rct["redcap_event_name"].replace(
        {"v1_screening_arm_1": "screening", "p5_phone_visit_arm_1": "treatment_period_2",
         "v2_gfr_mri_arm_1": "baseline", "v3_arm_1": "baseline", "v4_arm_1": "treatment_period_1",
         "v7_gfr_mri_arm_1": "post_treatment", "v61_med_dispense_arm_1": "treatment_period_3",
         "v62_med_dispense_arm_1": "treatment_period_4", "v8_arm_1": "post_treatment"},
        inplace=True)
    rct = rct.rename(columns={"redcap_event_name": "visit"})

    exclude_ids = [
        "RPC-06", "RPC-09", "RPC-10", "RPC-11", "RPC-12", "RPC-14", "RPC-15",
        "RPC-18", "RPC-19", "RPC-20", "RPC-22", "RPC-23", "RPC-24", "RPC-28",
        "RPC-30", "RPC-31", "RPC-32", "RPC-33", "RPC-34", "RPC-35", "RPC-36",
        "RPC-39", "RPC-40", "RPC-42", "RPC-43", "RPC-44",
    ]
    rct = rct[~rct["record_id"].isin(exclude_ids)]

    visits = ["baseline", "post_treatment"]
    rct_filtered = rct[rct["visit"].isin(visits)].copy()
    rows = []
    for rid in rct["record_id"].unique():
        for v in visits:
            subset = rct_filtered[(rct_filtered["record_id"] == rid) & (rct_filtered["visit"] == v)]
            if subset.empty:
                rows.append({"record_id": rid, "visit": v, **{h: np.nan for h in REDCAP_TO_HARMONIZED.values()}})
            else:
                row = subset.iloc[0]
                rows.append({"record_id": rid, "visit": v, **{h: row[r] for r, h in REDCAP_TO_HARMONIZED.items()}})

    return pd.DataFrame(rows)



def harmonized_clearance_data(base_data_path):
    # Load data dictionary and get renal clearance variable names
    dd = pd.read_csv(base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv")
    renal_vars = dd.loc[dd["form_name"] == "renal_clearance", "variable_name"].tolist()

    # Load harmonized dataset, filter to renal clearance columns
    harm = pd.read_csv(base_data_path + "Data Harmonization/Data Clean/harmonized_dataset.csv", low_memory=False)
    renal_cols_present = [v for v in renal_vars if v in harm.columns]
    harm_sub = harm[["record_id", "visit"] + renal_cols_present].copy()

    # Collapse to one row per (record_id, visit): mean for numeric, first non-missing for categorical
    num_cols = [c for c in renal_cols_present if pd.api.types.is_numeric_dtype(harm_sub[c])]
    cat_cols = [c for c in renal_cols_present if c not in num_cols]
    agg_dict = {**{c: "mean" for c in num_cols}, **{c: "first" for c in cat_cols}}
    harm_agg = harm_sub.groupby(["record_id", "visit"], as_index=False).agg(agg_dict)

    harm_key_vars = ["gfr_raw_plasma", "gfr_bsa_plasma", "erpf_raw_plasma", "erpf_bsa_plasma"]
    return harm_agg[["record_id", "visit"] + harm_key_vars]


if __name__ == "__main__":
    base_data_path, tokens = get_paths_and_tokens()

    panther = panther_redcap_data(base_data_path, tokens)
    panther["study"] = "PANTHER"

    crocodile = crocodile_redcap_data(base_data_path, tokens)
    crocodile["study"] = "CROCODILE"

    penguin = penguin_redcap_data(base_data_path, tokens)
    penguin["study"] = "PENGUIN"

    rpc2 = rpc2_redcap_data(base_data_path, tokens)
    rpc2["study"] = "RPC2"

    output = pd.concat([panther, penguin, rpc2], ignore_index=True)

    key_vars = ["gfr_raw_plasma", "gfr_bsa_plasma", "erpf_raw_plasma", "erpf_bsa_plasma"]
    n_missing = output[key_vars].isna().sum(axis=1)
    output["status"] = n_missing.map(lambda n: "complete" if n == 0 else ("all missing" if n == 4 else f"{n} of 4 missing"))

    output = output[["study", "record_id", "visit"] + key_vars + ["status"]]

    status_order = {"all missing": 0, "1 of 4 missing": 1, "2 of 4 missing": 2, "3 of 4 missing": 3, "complete": 4}
    output["status_order"] = output["status"].map(status_order)
    output = output.sort_values(by=["status_order", "study", "record_id"]).drop(columns="status_order").reset_index(drop=True)

    print(f"Total record_id x visit combinations: {len(output)}")
    output.to_csv(base_data_path + "Data Harmonization/Data Clean/clearance_inventory.csv", index=False)
