"""
This code is designed to pull data from the RENAL HEIR REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
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
token = tokens.loc[tokens["Study"] == "Renal-HEIR", "Token"].iloc[0]
proj = redcap.Project(url=uri, token=token)
# Get project metadata
meta = pd.DataFrame(proj.metadata)

# ------------------------------------------------------------------------------
# Demographics
# ------------------------------------------------------------------------------

dem_cols = ["subject_id", "co_enroll_id", "dob", "diagnosis",
            "group", "gender", "race", "ethnicity"]
# Export
demo = pd.DataFrame(proj.export_records(fields=dem_cols))
demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date"},
            inplace=True, axis=1)
dem_cols[3] = "diabetes_dx_date"
dem_cols[5] = "sex"
# Race columns combined into one
demo = combine_checkboxes(demo, base_name="race", levels=[
                          "American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
# Same for ethnicity
demo = combine_checkboxes(demo,
                          base_name="ethnicity",
                          levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
# Relevel sex and group
demo["sex"].replace({1: "Male", 0: "Female", 2: "Other",
                     "1": "Male", "0": "Female", "2": "Other"}, inplace=True)
demo["group"].replace({2: "Type 2 Diabetes", 3: "Obese Control",
                       4: "Lean Control",
                       "2": "Type 2 Diabetes", "3": "Obese Control",
                       "4": "Lean Control"}, inplace=True)

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
phys.drop(["male_activity_factor", "fem_activity_factor", "schofield_male",
           "schofield_female", "phys_norm", "phys_no", "breast_tanner", "testicular_volume", "lmp", "screen_bmi_percentile"], axis=1, inplace=True)
phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
phys.rename({"sysbp": "sbp", "diasbp": "dbp",
            "waist_circumference": "waistcm", "hip_circumference": "hipcm"}, inplace=True, axis=1)

# ------------------------------------------------------------------------------
# Screening labs
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "screening_labs", "field_name"]]
screen = pd.DataFrame(proj.export_records(fields=var))
screen.drop(["a1c_pre", "a1c_pre_date", "screen_pregnant"],
            axis=1, inplace=True)
screen.columns = screen.columns.str.replace(
    r"_of_screen|screen_", "", regex=True)
screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
              "urine_mab": "mab_u", "urine_cre": "creatinine_u"},
              axis=1, inplace=True)
screen["procedure"] = "screening"

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
# Clamp
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "clamp", "field_name"]]
clamp = pd.DataFrame(proj.export_records(fields=var))
clamp.drop(["baseline", "fasting_labs", "urine_labs", "hct_lab",
           "bg_labs", "ffa_lab", "cpep_lab", "insulin_labs"], axis=1, inplace=True)
clamp.columns = clamp.columns.str.replace(
    r"clamp_", "", regex=True)
clamp.rename({"serum_creatinine": "creatinine_s", "serum_sodium": "sodium_s"},
             inplace=True, axis=1)
clamp.columns = clamp.columns.str.replace(r"clamp_", "", regex=True)
clamp["procedure"] = "clamp"

# ------------------------------------------------------------------------------
# Outcomes
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "outcomes", "field_name"]]
out = pd.DataFrame(proj.export_records(fields=var))
out.drop(["kidney_outcomes", "egfr", "metab_outcomes", "asl_outcomes"],
         axis=1, inplace=True)
out.columns = out.columns.str.replace(
    r"mri_", "", regex=True)
out["procedure"] = "kidney_outcomes"

# ------------------------------------------------------------------------------
# Kidney Biopsy
# ------------------------------------------------------------------------------

var = ["subject_id", ] + [v for v in meta.loc[meta["form_name"]
                                              == "kidney_biopsy", "field_name"]]
var = var + ["gloms", "gloms_gs", "ifta", "vessels_other", "fia",
             "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins",
             "glom_volume_con", "mes_matrix_area",
             "mes_index", "mes_volume_weibel", "mes_volume_wiggins",
             "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima",
             "art_media", "pod_nuc_density", "pod_cell_volume"]
biopsy = pd.DataFrame(proj.export_records(fields=var))
biopsy = biopsy.loc[biopsy["bx_date"] != ""]
biopsy.drop([col for col in biopsy.columns if '_yn' in col] +
            [col for col in biopsy.columns if 'procedure_' in col],
            axis=1, inplace=True)
biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)
biopsy.columns = biopsy.columns.str.replace(r"labs_", "", regex=True)
biopsy.columns = biopsy.columns.str.replace(r"vitals_", "", regex=True)
biopsy["procedure"] = "kidney_biopsy"

# MERGE
renal_heir = pd.merge(phys, screen, how="outer")
renal_heir = pd.merge(renal_heir, dxa, how="outer")
renal_heir = pd.merge(renal_heir, clamp, how="outer")
renal_heir = pd.merge(renal_heir, out, how="outer")
renal_heir = pd.merge(renal_heir, biopsy, how="outer")
renal_heir = pd.merge(renal_heir, demo, how="outer")
# REORGANIZE
renal_heir["visit"] = "baseline"
renal_heir["study"] = "RENAL-HEIR"
id_cols = ["subject_id", "co_enroll_id", "study"] + \
    dem_cols[1:] + ["visit", "procedure", "date"]
other_cols = renal_heir.columns.difference(id_cols).tolist()
renal_heir = renal_heir[id_cols + other_cols]
# SORT
renal_heir.sort_values(["subject_id", "date", "procedure"], inplace=True)