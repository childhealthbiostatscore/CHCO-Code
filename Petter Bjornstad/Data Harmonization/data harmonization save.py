# Libraries
import os
import sys
sys.path.insert(0, os.path.expanduser('~') +
                  "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from data_harmonization import harmonize_data

clean = harmonize_data()
clean.to_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW//Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", index=False)
