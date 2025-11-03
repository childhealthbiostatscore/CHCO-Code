import pandas as pd
import numpy as np


def PFAS():
    import os
    import sys
    import pandas as pd
    sys.path.insert(0, os.path.expanduser('~') +
                    "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
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

    df = pd.read_excel(
    base_data_path + "PANTHER/Data_Raw/PFAS Emory 2025/8- PFAS Quantification/CLU0126-PFAS_Concentrations_Final_Report.xlsx",
    header=None,           # don't auto-detect headers
    skiprows=3,            # skip rows above header
    usecols="B:Z",         # only these columns
    nrows=95               # stop after 40 rows of data
    )
    df.columns = df.iloc[0]
    df = df.drop(0).reset_index(drop=True)

    df.rename({"N-ethylperfluoro-1-octanesulfonamidoacetic_acid_(N-EtFOSAA)": "N-EtFOSAA",	"N-methylperfluoro-1-octanesulfonamidoacetic_acid_(N-MeFOSAA)":"N-MeFOSAA",
                "Perfluoro-n-butanoic_acid-Deg_(PFBA)": "PFBA", "Perfluoro-n-decanoic_acid_(PFDA)": "PFDA", "Perfluoro-n-dodecanoic_acid_(PFDoA)":"PFDoA", 
                "Perfluoro-n-heptanoic_acid_(PFHpA)":"PFHpA", "Perfluoro-n-hexanoic_acid-Deg_(PFHxA)":"PFHxA",	"Perfluoro-n-nonanoic_acid_(PFNA)":"PFNA", 
                "Perfluoro-n-octanoic_acid_(PFOA)":"PFOA",	"Perfluoro-n-pentanoic_acid-Deg_(PFPeA)":"PFPeA", "Perfluoro-n-tetradecanoic_acid_(PFTeDA)":"PFTeDA",
                "Perfluoro-n-tridecanoic_acid_(PFTrDA)":"PFTrDA",	"Perfluoro-n-undecanoic_acid_(PFUnA)":"PFUnA",	"Perfluorobutane-1-sulfonic_acid_(PFBS)":"PFBS",
            	"Perfluorodecane-1-sulfonic_acid_(PFDoS)":"PFDoS",	"Perfluoroheptanesulfonic_acid_(PFHps)":"PFHps", "Perfluorohexane-1-sulfonic_acid_(PFHxS)":"PFHxS",
            	"Perfluorononanesulfonic_acid_(PFNS)":"PFNS", "Perfluorooctane-1-sulfonic_acid_(PFOS)":"PFOS",	"Perfluorooctane_sulfonamide_(PFOSA)":"PFOSA",
            	"Perfluoropentanesulfonic_acid_(PFPeAS)":"PFPeAS"}, axis=1, inplace=True)
    df["Sample.ID"] = df["Sample.ID"].str.removesuffix("_2")
    df["Sample.ID"] = df["Sample.ID"].str.removesuffix("-W")
    df["Sample.ID"] = df["Sample.ID"].str.removesuffix("_W_2")    

    df["visit"] = "baseline"
    df["procedure"] = "pfas"
    
    # Define replacements
    replacements = {
        "PAN-110-O": "PAN-110-C",
        "PAN-61-C": "PAN-61-O",
        "PAN-71-C": "PAN-71-O",
        "PAN-16-O-W": "PAN-16-O",
        "PAN-08-T-W": "PAN-08-T",
        "PAN-56-C-W": "PAN-56-C",
        "PAN-57-C-W": "PAN-57-C"
    }
    
    # Replace values
    df["Sample.ID"] = df["Sample.ID"].replace(replacements)
    
    df.drop(columns=['Run_order', 'Batch', 'Sample_type'], inplace=True)

    return df
