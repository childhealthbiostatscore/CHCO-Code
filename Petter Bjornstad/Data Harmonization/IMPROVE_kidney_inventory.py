import os
import sys
import getpass
import numpy as np
import pandas as pd
import redcap
from datetime import date


user = getpass.getuser()
base_paths = {
    "shivaniramesh": os.path.expanduser("~/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/"),
    "choiyej": "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/",
}
if user not in base_paths:
    sys.exit(f"Unknown user: {user}. Add your path to base_paths.")
BASE = base_paths[user]
OUTPUT_DIR = os.path.join(BASE, "IMPROVE T2D", "Data_Cleaned")
os.makedirs(OUTPUT_DIR, exist_ok=True)

rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999, -9999.0]
rep = rep + [str(r) for r in rep] + [""]

VISITS = ["screening", "baseline", "3_months_post_surgery", "12_months_post_surgery"]
VISIT_MAP = {
    np.nan: "baseline", "nan": "baseline",
    "1": "baseline", "2": "3_months_post_surgery", "3": "12_months_post_surgery",
}


FIELD_MAP = [
    ("uric_acid_med",                       "patient_medications",     "direct"),
    ("mra_med",                             "patient_medications",     "direct"),
    ("screen_serum_creatinine",             "screening_labs",          "direct"),
    ("screen_urine_mab",                    "screening_labs",          "direct"),
    ("screen_urine_cre",                    "screening_labs",          "direct"),
    ("screen_urine_acr",                    "screening_labs",          "direct"),
    ("mmtt_bun_base",                       "mmtt_metabolic_cart",     "direct"),
    ("mmtt_creat_base",                     "mmtt_metabolic_cart",     "direct"),
    ("cystatin_c",                          "clamp",                   "direct"),
    ("serum_creatinine",                    "clamp",                   "direct"),
    ("clamp_urine_mab_baseline",            "clamp",                   "direct"),
    ("clamp_urine_cre_baseline",            "clamp",                   "direct"),
    ("clamp_acr_baseline",                  "clamp",                   "direct"),
    ("clamp_urine_sodium",                  "clamp",                   "direct"),
    ("clamp_glucose_bl",                    "clamp",                   "direct"),
    ("urine_glucose",                       "clamp",                   "direct"),
    ("clamp_urine_mab_250",                 "clamp",                   "direct"),
    ("clamp_urine_cre_250",                 "clamp",                   "direct"),
    ("clamp_acr_250",                       "clamp",                   "direct"),
    ("clamp_urine_vol",                     "clamp",                   "direct"),
    ("iohexol_bolus",                       "clamp",                   "direct"),
    ("iohexol_time",                        "clamp",                   "direct"),
    ("iohexol_vol",                         "clamp",                   "direct"),
    ("iohexol_120",                         "clamp",                   "direct"),
    ("iohexol_150",                         "clamp",                   "direct"),
    ("iohexol_180",                         "clamp",                   "direct"),
    ("iohexol_210",                         "clamp",                   "direct"),
    ("iohexol_240",                         "clamp",                   "direct"),
    ("pah_bolus",                           "clamp",                   "direct"),
    ("pah_time",                            "clamp",                   "direct"),
    ("pah_vol",                             "clamp",                   "direct"),
    ("pah_minus_10",                        "clamp",                   "direct"),
    ("pah_90",                              "clamp",                   "direct"),
    ("pah_120",                             "clamp",                   "direct"),
    ("gfr",                                 "outcomes",                "direct"),
    ("gfr_bsa",                             "outcomes",                "direct"),
    ("abs_pah",                             "outcomes",                "direct"),
    ("pah_bsa",                             "outcomes",                "direct"),
    ("rpf",                                 "outcomes",                "direct"),
    ("erpf_bsa",                            "outcomes",                "direct"),
    ("ecv",                                 "outcomes",                "direct"),
    ("gfr_ecv_percent",                     "outcomes",                "direct"),
    ("gfr_ecv_std",                         "outcomes",                "direct"),
    ("asl_right",                           "outcomes",                "direct"),
    ("asl_left",                            "outcomes",                "direct"),
    ("bold_r_bl_cortex",                    "outcomes",                "direct"),
    ("bold_r_bl_medulla",                   "outcomes",                "direct"),
    ("bold_r_bl_kidney",                    "outcomes",                "direct"),
    ("bold_r_pf_cortex",                    "outcomes",                "direct"),
    ("bold_r_pf_medulla",                   "outcomes",                "direct"),
    ("bold_r_pf_kidney",                    "outcomes",                "direct"),
    ("bold_l_bl_cortex",                    "outcomes",                "direct"),
    ("bold_l_bl_medulla",                   "outcomes",                "direct"),
    ("bold_l_bl_kidney",                    "outcomes",                "direct"),
    ("bold_l_pf_cortex",                    "outcomes",                "direct"),
    ("bold_l_pf_medulla",                   "outcomes",                "direct"),
    ("bold_l_pf_kidney",                    "outcomes",                "direct"),
    ("adc_right",                           "outcomes",                "direct"),
    ("adc_left",                            "outcomes",                "direct"),
    ("volume_right",                        "outcomes",                "direct"),
    ("volume_left",                         "outcomes",                "direct"),
    ("volume_right_manual",                 "outcomes",                "direct"),
    ("volume_left_manual",                  "outcomes",                "direct"),
    ("gloms",                               "biopsy_results_michigan", "direct"),
    ("glom_enlarge___1",                    "biopsy_results_michigan", "direct"),
    ("glom_enlarge___2",                    "biopsy_results_michigan", "direct"),
    ("glom_enlarge___3",                    "biopsy_results_michigan", "direct"),
    ("glom_enlarge___4",                    "biopsy_results_michigan", "direct"),
    ("glom_enlarge___5",                    "biopsy_results_michigan", "direct"),
    ("glom_enlarge___6",                    "biopsy_results_michigan", "direct"),
    ("glom_enlarge___7",                    "biopsy_results_michigan", "direct"),
    ("gloms_gs",                            "biopsy_results_michigan", "direct"),
    ("ifta",                                "biopsy_results_michigan", "direct"),
    ("fia",                                 "biopsy_results_michigan", "direct"),
    ("glom_tuft_area",                      "kidney_biopsy",           "external"),
    ("glom_volume_weibel",                  "kidney_biopsy",           "external"),
    ("glom_volume_wiggins",                 "kidney_biopsy",           "external"),
    ("glom_volume_con",                     "kidney_biopsy",           "external"),
    ("mes_matrix_area",                     "kidney_biopsy",           "external"),
    ("mes_index",                           "kidney_biopsy",           "external"),
    ("mes_volume_weibel",                   "kidney_biopsy",           "external"),
    ("mes_volume_wiggins",                  "kidney_biopsy",           "external"),
    ("mes_volume_con",                      "kidney_biopsy",           "external"),
    ("gbm_thick_artmean",                   "kidney_biopsy",           "external"),
    ("gbm_thick_harmmean",                  "kidney_biopsy",           "external"),
    ("az_creatine_p",                       "az_urine_metabolites",    "direct"),
    ("az_hippuric_acid_p",                  "az_urine_metabolites",    "direct"),
    ("az_kynurenic_acid_n",                 "az_urine_metabolites",    "direct"),
    ("az_kynurenine_p",                     "az_urine_metabolites",    "direct"),
    ("az_oxalic_acid_n",                    "az_urine_metabolites",    "direct"),
    ("az_taurine_p",                        "az_urine_metabolites",    "direct"),
    ("az_trimethylamine_n_oxide_tmao_p",    "az_urine_metabolites",    "direct"),
]

BASELINE_ONLY_FORMS = {
    "az_urine_metabolites", "medications", "patient_medications",
    "kidney_biopsy", "biopsy_results_michigan",
}
SCREENING_ONLY_FORMS = {"screening_labs"}

VISITS_NOT_DONE = {
    ("IT_01", "3_months_post_surgery"),
    ("IT_01", "12_months_post_surgery"),
    ("IT_04", "3_months_post_surgery"),
    ("IT_04", "12_months_post_surgery"),
    ("IT_07", "3_months_post_surgery"),
    ("IT_09", "3_months_post_surgery"),
    ("IT_13", "12_months_post_surgery"),
    ("IT_17", "3_months_post_surgery"),
    ("IT_17", "12_months_post_surgery"),
    ("IT_18", "baseline"),
    ("IT_18", "3_months_post_surgery"),
    ("IT_18", "12_months_post_surgery"),
    ("IT_20", "3_months_post_surgery"),
    ("IT_20", "12_months_post_surgery"),
    ("IT_22", "3_months_post_surgery"),
    ("IT_23", "3_months_post_surgery"),
    ("IT_23", "12_months_post_surgery"),
}

VISIT_NOTES = {
    ("IT_07", "12_months_post_surgery"): "12mo clamp row in harmonized is phantom (demographics only, no clamp data)",
    ("IT_09", "12_months_post_surgery"): "Partial urine labs only (ACR, creatinine_u, microalbumin_u); no GFR/ERPF",
    ("IT_10", "12_months_post_surgery"): "12mo clamp + outcomes absent in REDCap — ask Tyler if visit was done",
    ("IT_19", "3_months_post_surgery"):  "3mo clamp + outcomes absent in REDCap — ask Tyler if visit was done",
    ("IT_19", "12_months_post_surgery"): "12mo clamp + outcomes absent in REDCap — ask Tyler if visit was done",
    ("IT_20", "baseline"):               "Baseline outcomes not collected (physical constraints); phantom GFR/ERPF row in harmonized — ask Tyler to remove from REDCap",
    ("IT_22", "12_months_post_surgery"): "12mo outcomes not entered in REDCap — ask Tyler",
}

BIOPSY_VISITS = {
    ("IT_07", "baseline"), ("IT_07", "12_months_post_surgery"),
    ("IT_08", "baseline"), ("IT_08", "12_months_post_surgery"),
    ("IT_09", "baseline"),
    ("IT_10", "baseline"), ("IT_10", "12_months_post_surgery"),
    ("IT_11", "baseline"), ("IT_11", "12_months_post_surgery"),
    ("IT_12", "baseline"), ("IT_12", "12_months_post_surgery"),
    ("IT_13", "baseline"),
    ("IT_14", "12_months_post_surgery"),
    ("IT_19", "baseline"),
}

MRI_VISITS = {
    ("IT_23", "baseline"),
    ("IT_01", "baseline"),
    ("IT_02", "12_months_post_surgery"), ("IT_02", "3_months_post_surgery"), ("IT_02", "baseline"),
    ("IT_03", "12_months_post_surgery"), ("IT_03", "3_months_post_surgery"), ("IT_03", "baseline"),
    ("IT_04", "baseline"),
    ("IT_05", "3_months_post_surgery"),  ("IT_05", "baseline"),
    ("IT_06", "12_months_post_surgery"), ("IT_06", "3_months_post_surgery"), ("IT_06", "baseline"),
    ("IT_07", "baseline"),               ("IT_07", "12_months_post_surgery"),
    ("IT_08", "12_months_post_surgery"), ("IT_08", "3_months_post_surgery"), ("IT_08", "baseline"),
    ("IT_09", "baseline"),               ("IT_09", "12_months_post_surgery"),
    ("IT_10", "12_months_post_surgery"), ("IT_10", "3_months_post_surgery"), ("IT_10", "baseline"),
    ("IT_11", "12_months_post_surgery"), ("IT_11", "3_months_post_surgery"), ("IT_11", "baseline"),
    ("IT_12", "12_months_post_surgery"), ("IT_12", "3_months_post_surgery"), ("IT_12", "baseline"),
    ("IT_13", "3_months_post_surgery"),  ("IT_13", "baseline"),
    ("IT_14", "12_months_post_surgery"), ("IT_14", "3_months_post_surgery"), ("IT_14", "baseline"),
    ("IT_15", "12_months_post_surgery"), ("IT_15", "3_months_post_surgery"), ("IT_15", "baseline"),
    ("IT_16", "12_months_post_surgery"), ("IT_16", "3_months_post_surgery"), ("IT_16", "baseline"),
    ("IT_17", "baseline"),
    ("IT_18", "baseline"), ("IT_18", "3_months_post_surgery"), ("IT_18", "12_months_post_surgery"),
    ("IT_19", "baseline"), ("IT_19", "3_months_post_surgery"), ("IT_19", "12_months_post_surgery"),
    ("IT_20", "baseline"), ("IT_20", "3_months_post_surgery"), ("IT_20", "12_months_post_surgery"),
    ("IT_22", "baseline"), ("IT_22", "3_months_post_surgery"), ("IT_22", "12_months_post_surgery"),
}

SECTIONS = {
    "med": [
        "uric_acid_med", "mra_med",
    ],
    "screening": [
        "screen_serum_creatinine", "screen_urine_mab", "screen_urine_cre", "screen_urine_acr",
    ],
    "labs": [
        "cystatin_c", "serum_creatinine", "clamp_urine_mab_baseline", "clamp_urine_cre_baseline",
        "clamp_acr_baseline", "clamp_urine_sodium", "clamp_glucose_bl", "urine_glucose",
        "clamp_urine_mab_250", "clamp_urine_cre_250", "clamp_acr_250", "clamp_urine_vol",
    ],
    "clamp": [
        "iohexol_bolus", "iohexol_time", "iohexol_vol", "iohexol_120", "iohexol_150",
        "iohexol_180", "iohexol_210", "iohexol_240", "pah_bolus", "pah_time", "pah_vol",
        "pah_minus_10", "pah_90", "pah_120",
    ],
    "clearance": [
        "gfr", "gfr_bsa", "abs_pah", "pah_bsa", "rpf", "erpf_bsa", "ecv",
        "gfr_ecv_percent", "gfr_ecv_std",
    ],
    "mri": [
        "asl_right", "asl_left", "bold_r_bl_cortex", "bold_r_bl_medulla", "bold_r_bl_kidney",
        "bold_r_pf_cortex", "bold_r_pf_medulla", "bold_r_pf_kidney", "bold_l_bl_cortex",
        "bold_l_bl_medulla", "bold_l_bl_kidney", "bold_l_pf_cortex", "bold_l_pf_medulla",
        "bold_l_pf_kidney", "adc_right", "adc_left", "volume_right", "volume_left",
        "volume_right_manual", "volume_left_manual",
    ],
    "biopsy": [
        "gloms", "glom_enlarge___1", "glom_enlarge___2", "glom_enlarge___3", "glom_enlarge___4",
        "glom_enlarge___5", "glom_enlarge___6", "glom_enlarge___7", "gloms_gs", "ifta", "fia",
        "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins", "glom_volume_con",
        "mes_matrix_area", "mes_index", "mes_volume_weibel", "mes_volume_wiggins", "mes_volume_con",
        "gbm_thick_artmean", "gbm_thick_harmmean",
    ],
    "mmtt": [
        "mmtt_bun_base", "mmtt_creat_base",
    ],
    "urine_metab": [
        "az_creatine_p", "az_hippuric_acid_p", "az_kynurenic_acid_n", "az_kynurenine_p",
        "az_oxalic_acid_n", "az_taurine_p", "az_trimethylamine_n_oxide_tmao_p",
    ],
}

SKIP_VALUES = {"", "nan", "--"}


def connect_improve(tokens):
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "IMPROVE", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    meta = pd.DataFrame(proj.metadata)
    return proj, meta


def pull_form(proj, meta, form_name):
    base = ["subject_id", "study_visit"]
    form_fields = meta.loc[meta["form_name"] == form_name, "field_name"].tolist()
    var = list(dict.fromkeys(base + form_fields))
    df = pd.DataFrame(proj.export_records(fields=var))
    df.replace(rep, np.nan, inplace=True)
    for col in ["redcap_event_name", "redcap_repeat_instrument", "redcap_repeat_instance"]:
        if col in df.columns:
            df.drop(columns=col, inplace=True)
    if form_name in BASELINE_ONLY_FORMS:
        df["visit"] = "baseline"
    elif form_name in SCREENING_ONLY_FORMS:
        df["visit"] = "screening"
    else:
        sv = df["study_visit"].astype(str).str.strip()
        df["visit"] = sv.map(lambda v: VISIT_MAP.get(v, "baseline"))
    df.rename(columns={"subject_id": "record_id"}, inplace=True)
    return df


def build_redcap_lookup(proj, meta):
    form_to_fields = {}
    for rc_field, form, kind in FIELD_MAP:
        if kind != "direct":
            continue
        form_to_fields.setdefault(form, set()).add(rc_field)

    lookup = {}
    for form, fields in form_to_fields.items():
        print(f"  Pulling: {form} ...")
        try:
            df = pull_form(proj, meta, form)
        except Exception as e:
            print(f"    ERROR: {e}")
            continue
        for rc_field in fields:
            if rc_field not in df.columns:
                print(f"    WARNING: {rc_field!r} not in {form} export")
                continue
            for _, row in df.iterrows():
                key = (str(row["record_id"]), row["visit"], rc_field)
                if key not in lookup or pd.isna(lookup[key]):
                    lookup[key] = row[rc_field]
    return lookup


def build_harmonized_lookup(harm_improve):
    external_cols = list(dict.fromkeys(
        rc_field for rc_field, form, kind in FIELD_MAP if kind == "external"
    ))
    lookup = {}
    for (rid, visit), grp in harm_improve.groupby(["record_id", "visit"]):
        for col in external_cols:
            if col not in harm_improve.columns:
                val = np.nan
            else:
                non_null = grp[col].dropna()
                val = non_null.iloc[0] if len(non_null) else np.nan
            lookup[(str(rid), visit, col)] = val
    return lookup


def build_completeness(output):
    rows = []
    for _, row in output.iterrows():
        rid, visit = row["record_id"], row["visit"]
        r = {"record_id": rid, "visit": visit, "notes": row.get("notes", "")}
        total_present = 0
        total_vars = 0

        for section, cols in SECTIONS.items():
            present_cols = [c for c in cols if c in output.columns]
            n_total = len(present_cols)

            if visit == "screening" and section != "screening":
                r[section] = "--"
                continue
            if section == "biopsy" and (rid, visit) not in BIOPSY_VISITS:
                r[section] = "--"
                continue
            if section == "mri" and (rid, visit) not in MRI_VISITS:
                r[section] = "--"
                continue

            n_present = sum(
                1 for c in present_cols if str(row[c]).strip() not in SKIP_VALUES
            )
            r[section] = f"{n_present} of {n_total}"
            total_present += n_present
            total_vars += n_total

        r["_sort"] = total_present / total_vars if total_vars else 0
        rows.append(r)

    return (
        pd.DataFrame(rows)
        .sort_values("_sort")
        .drop(columns="_sort")
        .reset_index(drop=True)
    )



if __name__ == "__main__":
    today = date.today().strftime("%Y%m%d")
    tokens = pd.read_csv(BASE + "Data Harmonization/api_tokens.csv")

    print("Connecting to IMPROVE REDCap...")
    proj, meta = connect_improve(tokens)

    print("Pulling REDCap form data...")
    rc_lookup = build_redcap_lookup(proj, meta)

    print("Loading harmonized dataset...")
    harm = pd.read_csv(
        BASE + "Data Harmonization/Data Clean/harmonized_dataset.csv",
        dtype={"record_id": str},
        low_memory=False,
    )
    harm_improve = harm[harm["study"] == "IMPROVE"].copy()
    record_ids = sorted(harm_improve["record_id"].unique())
    print(f"IMPROVE: {len(record_ids)} participants, {len(harm_improve)} rows")

    h_lookup = build_harmonized_lookup(harm_improve)

    fields_ordered = list(dict.fromkeys(rc_field for rc_field, form, kind in FIELD_MAP))
    kind_map = {rc_field: kind for rc_field, form, kind in FIELD_MAP}

    def get_value(rid, visit, rc_field):
        if kind_map[rc_field] == "external":
            val = h_lookup.get((str(rid), visit, rc_field), np.nan)
        else:
            val = rc_lookup.get((str(rid), visit, rc_field), np.nan)
        return "" if pd.isna(val) else val

    # Build inventory
    rows = []
    for rid in record_ids:
        for visit in VISITS:
            row = {"record_id": rid, "visit": visit}
            for rc_field in fields_ordered:
                row[rc_field] = get_value(rid, visit, rc_field)
            rows.append(row)

    output = pd.DataFrame(rows)

    # Drop visits not done
    output = output[~output.apply(
        lambda r: (r["record_id"], r["visit"]) in VISITS_NOT_DONE, axis=1
    )].reset_index(drop=True)

    output["notes"] = output.apply(
        lambda r: VISIT_NOTES.get((r["record_id"], r["visit"]), ""), axis=1
    )
    output = output[["record_id", "visit", "notes"] + fields_ordered]

    mri_cols  = [c for c in SECTIONS["mri"]    if c in output.columns]
    biopsy_cols = [c for c in SECTIONS["biopsy"] if c in output.columns]
    for i, row in output.iterrows():
        if (row["record_id"], row["visit"]) not in MRI_VISITS:
            for c in mri_cols:
                if str(row[c]).strip() in ("", "nan"):
                    output.at[i, c] = "--"
        if (row["record_id"], row["visit"]) not in BIOPSY_VISITS:
            for c in biopsy_cols:
                if str(row[c]).strip() in ("", "nan"):
                    output.at[i, c] = "--"

    inv_path = os.path.join(OUTPUT_DIR, f"IMPROVE_clamp_outcomes_inventory_{today}.csv")
    output.to_csv(inv_path, index=False)
    print(f"\nWrote inventory: {inv_path}  ({len(output)} rows × {len(output.columns)} columns)")

    completeness = build_completeness(output)
    comp_path = os.path.join(OUTPUT_DIR, f"IMPROVE_completeness_by_section_{today}.csv")
    completeness.to_csv(comp_path, index=False)
    print(f"Wrote completeness: {comp_path}  ({len(completeness)} rows × {len(completeness.columns)} columns)")
