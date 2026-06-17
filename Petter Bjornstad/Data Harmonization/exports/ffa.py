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
os.chdir(os.path.expanduser('~'))
os.chdir("GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
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
# Get dataset
df = harmonize_data()
# Subset columns
df = df[["record_id", "co_enroll_id", "study", "visit", "procedure", "insulin_sensitivity_method",
         "ffa_method", 'ffa_minus_20', 'ffa_minus_10', 'ffa_minus_5', 'ffa_0',
         'ffa_2', 'ffa_4', 'ffa_6', 'ffa_8', 'ffa_10', 'ffa_30', 'ffa_60',
         'ffa_70', 'ffa_80', 'ffa_90', 'ffa_120', 'ffa_180', 'ffa_220',
         'ffa_230', 'ffa_240', 'ffa_250', 'ffa_260', 'ffa_270', 'baseline_ffa',
         'steady_state_ffa', 'p1_steady_state_ffa', 'p2_steady_state_ffa',
         'ffa_suppression', 'p1_ffa_suppression', 'p2_ffa_suppression']]
# Subset rows
df = df.loc[df["procedure"].isin(["clamp", "mmtt"])]
df = df.loc[~df["study"].isin(["CASPER", "COFFEE"])]
# Write
df.to_csv("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/ffa_dataset.csv", index=False)
