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
sys.path.insert(0, os.path.expanduser('~') +
                "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
# __pipeline_path_bootstrap__ (files moved into subfolders; resolve repo root from this file)
import os as _os, sys as _sys
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for _sub in ("compile", "studies", "calculations", "exports",):
    _p = _os.path.join(_ROOT, _sub)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)
from data_harmonization import harmonize_data
from datetime import datetime
import pandas as pd
# Get dataset
df = harmonize_data()
# Write
df.to_csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", index=False)
# Update for Michigan
time = datetime.now().strftime('%Y_%m_%d_%I%M%p')
# Re-save dictionary
dictionary = pd.read_csv(
    "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/data_dictionary_master.csv")
dictionary.to_csv("/Users/timvigers/Library/CloudStorage/Dropbox/Work/Shared/Michigan/CHCO_Sample_IDs_Clinical/chco_data_dictionary_" +
                  time + ".csv", index=False)
# Save cleaned data (de-identified)
df = df.drop(["dob"], axis=1)
df.to_csv("/Users/timvigers/Library/CloudStorage/Dropbox/Work/Shared/Michigan/CHCO_Sample_IDs_Clinical/chco_harmonized_dataset_" +
          time + ".csv", index=False)
