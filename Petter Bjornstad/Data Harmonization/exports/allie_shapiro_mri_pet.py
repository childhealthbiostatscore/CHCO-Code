"""
This code is designed pulls our harmonized data and selects the variables Allie Shapiro is using in her neuroimaging analysis.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"

import os
os.chdir(
    "C:/Users/timbv/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
import pandas as pd
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
# Filter and select
mri = df[df["record_id"].isin(
    ["CRC-10", "CRC-12", "CRC-13", "CRC-18", "CRC-20", "CRC-24", "CRC-28", "CRC-31", "CRC-33", "CRC-34", "CRC-35"])]
mri = mri.groupby(by=["record_id"]).agg("last").reset_index()
# Select variables of interest
cols = ["record_id", "group", "age", "sex", "bmi",
        'adc_left', 'adc_right', 'pasl2d_left', 'pasl2d_right', 'pcasl3d_left', 'pcasl3d_right', "volume_left", "volume_right"] + \
    [c for c in mri.columns if "_bl_" in c] + \
    [c for c in mri.columns if "_pf_" in c] + \
    [c for c in mri.columns if "fsoc" in c] + \
    ["rc_k1", "rc_k2", "rc_f", "rc_vb", "rm_k1", "rm_k2", "rm_f", "rm_vb",
     "lc_k1", "lc_k2", "lc_f", "lc_vb", "lm_k1", "lm_k2", "lm_f", "lm_vb"]
mri = mri[cols]
# Write CSV
mri.to_csv("~/mri_and_pet.csv", index=False)
