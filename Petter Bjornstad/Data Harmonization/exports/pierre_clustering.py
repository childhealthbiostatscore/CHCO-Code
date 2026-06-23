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
import pandas as pd

# Get dataset
df = harmonize_data()
# Variable list from Pierre
dict = pd.read_csv(
    "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/variables_for_Pierre_clustering.csv"
)
dict = dict[dict["Unnamed: 2"] == "x"]
cols = [c for c in dict["variable_name"]] + ["hematocrit_210"] + ["tot_protein"]
# Select columns
df = df[cols]
# Select rows
df = df[df["study"].isin(["CROCODILE", "PANDA", "RENAL-HEIR", "RENAL-HEIRitage"])]
df = df[df["procedure"].isin(["screening", "bold_mri", "clamp", "kidney_biopsy", "medications", "renal_clearance_testing"])]
# Write
df.to_csv(
    "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/pierre_clustering.csv",
    index=False,
)
