"""
This code is designed pulls our harmonized data and selects the variables Fadhl
is using in his pre-/post-surgery analysis.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"
import os
import sys

sys.path.insert(
    0, os.path.expanduser("~") + "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization"
)
# __pipeline_path_bootstrap__ (files moved into subfolders; resolve repo root from this file)
import os as _os, sys as _sys
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for _sub in ("compile", "studies", "calculations", "exports",):
    _p = _os.path.join(_ROOT, _sub)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)
from data_harmonization import harmonize_data

# Get dataset
df = harmonize_data()
# Select
df = df[
    [
        "record_id",
        "co_enroll_id",
        "study",
        "group",
        "diabetes_dx_date",
        "age",
        "height",
        "sex",
        "race",
        "ethnicity",
        "visit",
        "bmi",
        "creatinine_s",
        "cystatin_c_s",
        "acr_u",
        "gfr_bsa_plasma",
        "gfr_bsa_plasma_urine",
        "gfr_raw_plasma",
        "gfr_raw_plasma_urine",
        "eGFR_bedside_Schwartz",
        "eGFR_CKD_epi",
        "eGFR_fas_cr",
        "eGFR_fas_cr_cysc",
        "eGFR_Zap",
        "eGFR_Schwartz",
    ]
]
# Group rows by visit, get non-missing values
# Also, limit to baseline visits
df["visit"] = df["visit"].replace({"pre_surgery": "baseline"})
df = df[df["visit"] == "baseline"]
df = df.groupby(by=["record_id"]).agg("last").reset_index()
df.to_csv(
    "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/nathan_bekelman.csv",
    index=False,
)
