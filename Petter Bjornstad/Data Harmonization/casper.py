"""
This code is designed to pull data from the CASPER REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"

# Libraries
import os
import redcap
import pandas as pd
import numpy as np
os.chdir("/Users/timvigers/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from harmonization_functions import combine_checkboxes
from harmonization_functions import find_duplicate_columns

# REDCap project variables
tokens = pd.read_csv(
    "~/Dropbox/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"
token = tokens.loc[tokens["Study"] == "CASPER", "Token"].iloc[0]
proj = redcap.Project(url=uri, token=token)
# Get project metadata
meta = pd.DataFrame(proj.metadata)

# ------------------------------------------------------------------------------
# Demographics
# ------------------------------------------------------------------------------

dem_cols = ["subject_id", "dob", "diagnosis", "gender", "race", "ethnicity"]
# Export
demo = pd.DataFrame(proj.export_records(fields=dem_cols))
demo["co_enroll_id"] = np.nan
demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date"},
            inplace=True, axis=1)
dem_cols[2] = "diabetes_dx_date"
dem_cols[3] = "sex"
# Race columns combined into one
demo = combine_checkboxes(demo, base_name="race", levels=[
                          "American Indian or Alaskan Native",
                          "Asian",
                          "Hawaiian or Pacific Islander",
                          "Black or African American",
                          "White",
                          "Unknown",
                          "Other"])
# Same for ethnicity
demo = combine_checkboxes(demo,
                          base_name="ethnicity",
                          levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
# Relevel sex and group
demo["sex"].replace({1: "Male", 0: "Female", 3: "Other",
                     "1": "Male", "0": "Female", "3": "Other"}, inplace=True)
demo["group"] = "Type 1 Diabetes"

# ------------------------------------------------------------------------------
# Medications
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Physical exam
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "physical_exam", "field_name"]]
phys = pd.DataFrame(proj.export_records(fields=var))
phys["procedure"] = "physical_exam"
phys.drop(["phys_norm", "phys_no", "breast_tanner",
           "testicular_volume", "lmp", "screen_bmi_percentile", "male_activity_factor", "fem_activity_factor", "schofield_male", "schofield_female"], axis=1, inplace=True)
phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
phys.rename({"sys_bp": "sbp", "dys_bp": "dbp", "waist_circumference": "waistcm",
            "hip_circumference": "hipcm"}, inplace=True, axis=1)

# ------------------------------------------------------------------------------
# Screening labs
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "screening_labs", "field_name"]]
screen = pd.DataFrame(proj.export_records(fields=var))
screen.drop(['a1c_pre', 'a1c_pre_date', "screen_pregnant"],
            axis=1, inplace=True)
screen.columns = screen.columns.str.replace(
    r"screen_|_of_screen", "", regex=True)
screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
              "urine_cre": "creatinine_u", "urine_mab": "mab_u"},
              axis=1, inplace=True)
screen["procedure"] = "screening"

# ------------------------------------------------------------------------------
# Clamp
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "clamp", "field_name"]]
clamp = pd.DataFrame(proj.export_records(fields=var))
clamp.drop(["baseline", "fasting_labs", "bg_labs", "urine_labs", "hct_lab",
            "a1c_clamp_time", "clamp_a1c", "clamp_a1c_date"],
           axis=1, inplace=True)
clamp = clamp.loc[clamp["clamp_date"] != ""]
clamp.columns = clamp.columns.str.replace(
    r"clamp_", "", regex=True)
clamp["procedure"] = "clamp"

# ------------------------------------------------------------------------------
# DXA Scan
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "body_composition_dxa", "field_name"]]
dxa = pd.DataFrame(proj.export_records(fields=var))
dxa.columns = dxa.columns.str.replace(
    r"dexa_", "", regex=True)
dxa["procedure"] = "dxa"

# ------------------------------------------------------------------------------
# Outcomes
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "outcomes", "field_name"]]
out = pd.DataFrame(proj.export_records(fields=var))
out.drop(["kidney_outcomes", "egfr", "metab_outcomes",
          "asl_outcomes", "bold_outcomes"],
         axis=1, inplace=True)
out = out.loc[out["mri_date"] != ""]
out.columns = out.columns.str.replace(
    r"mri_", "", regex=True)
out["procedure"] = "kidney_outcomes"
out.to_csv("~/out.csv")

# MERGE
casper = pd.merge(phys, screen, how="outer")
casper = pd.merge(casper, dxa, how="outer")
casper = pd.merge(casper, clamp, how="outer")
casper = pd.merge(casper, demo, how="outer")
# REORGANIZE
casper["visit"] = "baseline"
casper["study"] = "CASPER"
id_cols = ["subject_id", "co_enroll_id", "study"] + \
    dem_cols[1:] + ["visit", "procedure", "date"]
other_cols = casper.columns.difference(id_cols).tolist()
casper = casper[id_cols + other_cols]
# SORT
casper.sort_values(["subject_id", "date", "procedure"], inplace=True)
# Check for duplicated column names
dups = find_duplicate_columns(casper)
dups.to_csv("~/croc_duplicate_columns.csv", index=False)
# Print final data
casper.to_csv("~/casper.csv", index=False)