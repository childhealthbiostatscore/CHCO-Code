#def update_dict():
# Libraries
import os
import sys
sys.path.insert(0, os.path.expanduser('~') +
                "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
import pandas as pd
import numpy as np
from natsort import natsorted, ns
from harmonization_functions import calc_egfr, create_study_id_columns
import redcap  # Add this import for the redcap module
# Use individual data functions to import cleaned DFs
import getpass
user = getpass.getuser() 
if user == "choiyej":
    base_data_path = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/"
    git_path = "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
elif user == "pylell":
    base_data_path = "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
    git_path = "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
elif user == "shivaniramesh":
    base_data_path = os.path.expanduser("~/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
    git_path = "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
else:
    sys.exit(f"Unknown user: please specify root path for this user. (Detected user: {user})")
dictionary = pd.read_csv(base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv")

tokens = pd.read_csv(base_data_path + "Data Harmonization/api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"

studies = ["CASPER", "COFFEE", "CROCODILE", "IMPROVE", "PENGUIN",
            "RENAL_HEIR", "RENAL_HEIRITAGE", "PANTHER", "PANDA", "ATTEMPT"]

    # Connect to REDCap projects and get metadata
casper_token = tokens.loc[tokens["Study"] == "CASPER", "Token"].iloc[0]
casper = redcap.Project(url=uri, token=casper_token)
print(casper)
casper_metadata = pd.DataFrame(casper.metadata)
print(casper_metadata)
forms = casper_metadata["form_name"].unique()


coffee_token = tokens.loc[tokens["Study"] == "COFFEE", "Token"].iloc[0]
coffee = redcap.Project(url=uri, token=coffee_token)
coffee_metadata = pd.DataFrame(coffee.metadata)

crocodile_token = tokens.loc[tokens["Study"] == "CROCODILE", "Token"].iloc[0]
crocodile = redcap.Project(url=uri, token=crocodile_token)
crocodile_metadata = pd.DataFrame(crocodile.metadata)

improve_token = tokens.loc[tokens["Study"] == "IMPROVE", "Token"].iloc[0]
improve = redcap.Project(url=uri, token=improve_token)
improve_metadata = pd.DataFrame(improve.metadata)

penguin_token = tokens.loc[tokens["Study"] == "PENGUIN", "Token"].iloc[0]
penguin = redcap.Project(url=uri, token=penguin_token)
penguin_metadata = pd.DataFrame(penguin.metadata)

renal_heir_token = tokens.loc[tokens["Study"] == "Renal-HEIR", "Token"].iloc[0]
renal_heir = redcap.Project(url=uri, token=renal_heir_token)
renal_heir_metadata = pd.DataFrame(renal_heir.metadata)

renal_heiritage_token = tokens.loc[tokens["Study"] == "Renal-HEIRitage", "Token"].iloc[0]
renal_heiritage = redcap.Project(url=uri, token=renal_heiritage_token)
renal_heiritage_metadata = pd.DataFrame(renal_heiritage.metadata)

panther_token = tokens.loc[tokens["Study"] == "PANTHER", "Token"].iloc[0]
panther = redcap.Project(url=uri, token=panther_token)
panther_metadata = pd.DataFrame(panther.metadata)

panda_token = tokens.loc[tokens["Study"] == "PANDA", "Token"].iloc[0]
panda = redcap.Project(url=uri, token=panda_token)
panda_metadata = pd.DataFrame(panda.metadata)

attempt_token = tokens.loc[tokens["Study"] == "ATTEMPT", "Token"].iloc[0]
attempt = redcap.Project(url=uri, token=attempt_token)
attempt_metadata = pd.DataFrame(attempt.metadata)




studies = [casper_metadata, coffee_metadata, crocodile_metadata, improve_metadata, penguin_metadata,
            renal_heir_metadata, renal_heiritage_metadata, panther_metadata, panda_metadata, attempt_metadata]
variables = dictionary['variable_name'].tolist()

for study in studies:
    for variable_name, label in zip(study['field_name'], study['field_label']):

        if variable_name.lower() != variable_name or " " in variable_name:
            print(f" Suspicious variable_name: {variable_name} (label: {label})")

        if variable_name not in variables:
            new_row = pd.DataFrame({'variable_name': [variable_name], 'label': [label]})
            dictionary = pd.concat([dictionary, new_row], ignore_index=True)
            variables.append(variable_name)

        
dictionary['form_name'] = ''  # Initialize the new column
dictionary['field_type'] = ''  # Initialize the new column
all_metadata = pd.concat(studies, ignore_index=True)
for index, row in dictionary.iterrows():
    matching_row = all_metadata[all_metadata['field_name'] == row['variable_name']]
    if not matching_row.empty:
        # Get the form_name for the first match and assign it
        dictionary.at[index, 'form_name'] = matching_row['form_name'].iloc[0]
        dictionary.at[index, 'field_type'] = matching_row['field_type'].iloc[0]

dictionary = dictionary.drop_duplicates(subset=['variable_name', 'label'])

# Drop invalid labels (numeric 0, string "0", empty string, NaN)
invalid_labels = [0, '0', '']  

dictionary = dictionary[~dictionary['label'].isin(invalid_labels)]
dictionary = dictionary.dropna(subset=['label'])
      
dictionary = dictionary.drop_duplicates(subset=['variable_name', 'label'])

harmonized = pd.read_csv(base_data_path + "Data Harmonization/Data Clean/harmonized_dataset.csv")
tocsv_path = base_data_path + "Data Harmonization/Data Clean/data_dictionary_master.csv"
dictionary.to_csv(tocsv_path, index=False)

        
