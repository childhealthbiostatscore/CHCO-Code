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
from data_harmonization import harmonize_data
from datetime import datetime
import pandas as pd
# Get dataset
df = harmonize_data()
# Select
df = df[["record_id", "co_enroll_id", "study", "dob", "diabetes_dx_date", "sex",
        "race", "ethnicity", "visit", "gfr_bsa_plasma", "gfr_bsa_plasma_urine",
         "gfr_raw_plasma", "gfr_raw_plasma_urine", "eGFR_bedside_Schwartz",
         "eGFR_CKD_epi", "eGFR_fas_cr"]]
# Group rows by visit, get non-missing values
# Also, limit to baseline visits
df["visit"] = df["visit"].replace({"pre_surgery": "baseline"})
df = df[df["visit"] == "baseline"]
df = df.groupby(by=["record_id"]).agg("last").reset_index()
df.to_csv("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/nathan_bekelman.csv", index=False)