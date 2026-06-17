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
import sys
import pandas as pd
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
# Get dataset
df = harmonize_data()
# Filter and select
df = df[df["study"] == "CROCODILE"]
# Import neuro markers
ab40 = pd.read_excel(
    "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Cleaned/neuro_markers.xlsx",
    sheet_name="Ab40 Results")
ab40 = ab40[["record_id", 'AEB Rep1', 'AEB Rep2', 'Ave. AEB', 'CV',
             'Conc. Rep1 (pg/mL)', 'Conc. Rep2 (pg/mL)', 'Ave. Conc. (pg/mL)',
             'CV.1', 'Dilution Corrected Conc. (pg/mL)']]
ab40.columns = [c + " AB40" if c != "record_id" else c for c in ab40.columns]
ab42 = pd.read_excel(
    "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Cleaned/neuro_markers.xlsx",
    sheet_name="Ab42 Results")
ab42 = ab42[["record_id", 'AEB Rep1', 'AEB Rep2', 'Ave. AEB', 'CV',
             'Conc. Rep1 (pg/mL)', 'Conc. Rep2 (pg/mL)', 'Ave. Conc. (pg/mL)',
             'CV.1', 'Dilution Corrected Conc. (pg/mL)']]
ab42.columns = [c + " AB42" if c != "record_id" else c for c in ab42.columns]
tau = pd.read_excel(
    "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Cleaned/neuro_markers.xlsx",
    sheet_name="Tau Results")
tau = tau[["record_id", 'AEB Rep1', 'AEB Rep2', 'Ave. AEB', 'CV',
           'Conc. Rep1 (pg/mL)', 'Conc. Rep2 (pg/mL)', 'Ave. Conc. (pg/mL)',
           'CV.1', 'Dilution Corrected Conc. (pg/mL)']]
tau.columns = [c + " Tau" if c != "record_id" else c for c in tau.columns]
nfl = pd.read_excel(
    "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Cleaned/neuro_markers.xlsx",
    sheet_name="NFL Results")
nfl = nfl[["record_id", 'AEB Rep1', 'AEB Rep2', 'Ave. AEB', 'CV',
           'Conc. Rep1 (pg/mL)', 'Conc. Rep2 (pg/mL)', 'Ave. Conc. (pg/mL)',
           'CV.1', 'Dilution Corrected Conc. (pg/mL)']]
nfl.columns = [c + " NFL" if c != "record_id" else c for c in nfl.columns]
gfap = pd.read_excel(
    "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Cleaned/neuro_markers.xlsx",
    sheet_name="GFAP Results")
gfap = gfap[["record_id", 'AEB Rep1', 'AEB Rep2', 'Ave. AEB', 'CV',
             'Conc. Rep1 (pg/mL)', 'Conc. Rep2 (pg/mL)', 'Ave. Conc. (pg/mL)',
             'CV.1', 'Dilution Corrected Conc. (pg/mL)']]
gfap.columns = [c + " GFAP" if c != "record_id" else c for c in gfap.columns]
ptau181 = pd.read_excel(
    "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Cleaned/neuro_markers.xlsx",
    sheet_name="pTau 181 Results")
ptau181 = ptau181[["record_id", 'AEB Rep1', 'AEB Rep2', 'Ave. AEB', 'CV',
                   'Conc. Rep1 (pg/mL)', 'Conc. Rep2 (pg/mL)', 'Ave. Conc. (pg/mL)',
                   'CV.1', 'Dilution Corrected Conc. (pg/mL)']]
ptau181.columns = [c + " pTau 181" if c !=
                   "record_id" else c for c in ptau181.columns]
# Select variables of interest
cols = ["record_id", "group", "age", "diabetes_duration", "sex",
        "bmi", "sbp", "dbp", "insulin_pump_timepoint",
        "microalbumin_u", "p1_ffa_suppression", "p1_gc_leanm",
        "p1_gc_m", "p1_raw_leanm", "p1_raw_m", "p2_gc_leanm", "p2_gc_m",
        "p2_raw_leanm", "p2_raw_m"]
df = df[cols]
df = df.groupby(by=["record_id"]).agg("last").reset_index()
# Merge
df = pd.merge(df, ab40)
df = pd.merge(df, ab42)
df = pd.merge(df, tau)
df = pd.merge(df, nfl)
df = pd.merge(df, gfap)
df = pd.merge(df, ptau181)
# Write CSV
df.to_csv("~/crc_neuro_markers.csv", index=False)
