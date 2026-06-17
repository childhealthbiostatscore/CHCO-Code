"""
This code is designed pulls our harmonized data and selects the variables Fadhl
is using in his pre-/post-surgery analysis.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"

import os
os.chdir(
    "C:/Users/timbv/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
import pandas as pd
import numpy as np
# __pipeline_path_bootstrap__ (files moved into subfolders; resolve repo root from this file)
import os as _os, sys as _sys
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for _sub in ("compile", "studies", "calculations", "exports",):
    _p = _os.path.join(_ROOT, _sub)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)
from data_harmonization import harmonize_data
from natsort import natsorted, ns
# Get dataset
df = harmonize_data()
# Variable lists
vars = ["record_id", "co_enroll_id", "study", "kit_id", "group", "dob",
        "diabetes_dx_date", "race", "ethnicity", "race_ethnicity", "visit", "procedure", "date", "diabetes_duration", "age", "sex",
        "sglt2i_ever", "sglti_timepoint", "bmi", "sbp", "dbp", "creatinine_s", "cystatin_c_s", "bun",
        "gfr", "gfr_bsa", "eGFR_fas_cr",
        "acr_u", "hba1c", "dexa_body_fat", "dexa_fat_kg", "dexa_trunk_kg", "bod_pod_body_fat", "bod_pod_fat_kg",
        "bold_r_bl_cortex", "bold_r_bl_medulla", "bold_r_bl_kidney",
        "bold_r_pf_cortex", "bold_r_pf_medulla", "bold_r_pf_kidney",
        "bold_l_bl_cortex", "bold_l_bl_medulla", "bold_l_bl_kidney",
        "bold_l_pf_cortex", "bold_l_pf_medulla", "bold_l_pf_kidney",
        "gloms", "gloms_gs", "ifta", "vessels_other", "fia", "glom_tuft_area",
        "glom_volume_weibel", "glom_volume_wiggins", "glom_volume_con",
        "mes_matrix_area", "mes_index", "mes_volume_weibel", "mes_volume_wiggins", "mes_volume_con",
        "glom_nuc_count", "mes_nuc_count", "art_intima", "art_media", "pod_nuc_density", "pod_cell_volume", "he_clamp",
        "raw_m", "steady_state_insulin", "steady_state_cpeptide", "acprg", "airg", "di", "p1_raw_m", "p1_raw_leanm", "p1_gc_m", "p1_gc_leanm", "p2_raw_m", "p2_raw_leanm", "p2_gc_m", "p2_gc_leanm"] + \
    ["baseline_ffa", "steady_state_ffa", "ffa_supression"]
fadhl = df[vars]
fadhl = fadhl[fadhl["visit"] != "3_months_post_surgery"]
fadhl.to_csv("~/fadhl_pre_aggregation.csv", index=False)
# Group rows by visit, get non-missing values
fadhl["visit"] = fadhl["visit"].replace({"pre_surgery": "baseline"})
fadhl = fadhl.groupby(by=["record_id", "visit"]).agg("last").reset_index()
fadhl = fadhl.loc[~pd.isna(fadhl["study"])]
# Add Michigan ID
id = fadhl["record_id"].copy()
id.loc[(fadhl["study"] == "IMPROVE") & (fadhl["visit"] == "baseline")] += "_BL"
id.loc[(fadhl["study"] == "IMPROVE") & (
    fadhl["visit"] == "12_months_post_surgery")] += "_12M"
fadhl["michigan_id"] = id
fadhl.set_index(["michigan_id"], drop=True, inplace=True)
# Write
fadhl = fadhl.astype(object)
fadhl.to_csv("Z:/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Kidney Morphometrics/ADA 2023/Data_Cleaned/harmonized_latest_value.csv", index=False)
