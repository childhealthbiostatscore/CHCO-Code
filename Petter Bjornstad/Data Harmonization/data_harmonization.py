"""
This code is designed to pull data from multiple REDCap projects and harmonize
the data in a single dataset. Some studies are cross-sectional but include
measures at multiple visits, and some studies are longitudinal. So, this code
outputs data in a "semi-long" format with one row per study procedure, and a
visit column for longitudinal clustering.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"

# Libraries
import sys
import re
import redcap
import pandas as pd
import numpy as np
sys.path.insert(1,
                "/Users/timvigers/GitHub/shared-resources/Data Cleaning/Data Manipulation in Python")
from combine_redcap_checkboxes import combine_redcap_checkboxes
# REDCap import variables
tokens = pd.read_csv(
    "~/Dropbox/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"
croc_token = tokens.loc[tokens["Study"] == "CROCODILE", "Token"].iloc[0]
croc = redcap.Project(url=uri, token=croc_token)

# -----------------------------------------------------------------------------
# For each cross-sectional study, import procedure forms separately then stack
# longwise for "semi-long" data.
# -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# CROCODILE
# ------------------------------------------------------------------------------
croc_meta = pd.DataFrame(croc.metadata)
croc_df = croc.export_records(format_type="df")
# DEMOGRAPHICS
dem_cols = ["record_id", "dob", "diabetes_dx_date",
            "group", "sex", "race", "ethnicity"]
# Export
croc_demo = croc.export_records(format_type="df", fields=dem_cols)
# Race columns combined into one
croc_demo = combine_redcap_checkboxes(croc_demo, base_name="race",
                                      levels=[
                                          'American Indian/Alaska Native', 'Asian', 'Hawaiian/Pacific Islander', 'Black/African American', 'White', 'Other', 'Unknown'])
# Same for ethnicity
croc_demo = combine_redcap_checkboxes(croc_demo,
                                      base_name="ethnicity",
                                      levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
# Relevel sex and group
croc_demo["sex"].replace({1: "Male", 2: "Female", 3: "Other"}, inplace=True)
croc_demo["group"].replace({1: "Type 1 Diabetes", 2: "Control"}, inplace=True)
# PHYSICAL EXAM
var = [v for v in croc_meta.loc[croc_meta["form_name"]
                                == "physical_exam", "field_name"]]
phys = croc_df.filter(var, axis=1)
phys["procedure"] = "physical_exam"
phys.drop(["phys_normal", "phys_abnormal"], axis=1, inplace=True)
phys.columns = phys.columns.str.replace(r"phys_", "")
# SCREENING
var = [v for v in croc_meta.loc[croc_meta["form_name"]
                                == "screening_labs", "field_name"]]
screen = croc_df.filter(var, axis=1)
screen.drop(["prescreen_a1c", "prescreen_a1c_date",
            "screen_menstrual", "screen_upt"], axis=1, inplace=True)
screen.columns = screen.columns.str.replace(
    r"labs_|screen_", "", regex=True)
screen.rename({"creat_s": "creatinine_s",
              "creat_u": "creatinine_u"}, axis=1, inplace=True)
screen["procedure"] = "screening"


# MERGE
crocodile = pd.merge(phys, screen, on=["record_id", "date", "procedure"],
                     how="outer")
# REORGANIZE
crocodile["visit"] = "baseline"
crocodile["study"] = "CROCODILE"
id_cols = ["study", "visit", "procedure", "date"]
other_cols = crocodile.columns.difference(id_cols, sort=False).tolist()
crocodile = crocodile[id_cols + other_cols]
# SORT
crocodile.sort_values(["record_id", "date"], inplace=True)
# Print final data
crocodile.to_csv("~/crocodile.csv")