"""
Helper functions for data harmonization.
"""


def combine_checkboxes(df, base_name="", levels=[], sep=" & ",
                       drop=True):
    # Libraries
    import pandas as pd
    import numpy as np
    # Make a copy of the dataframe
    data = df.copy()
    # A list of integer strings (one for each level)
    num = [str(s) for s in range(1, len(levels) + 1)]
    # Get column names
    cols = [base_name + "___" + s for s in num]
    # Which columns are checked
    values = data[cols].apply(lambda x: " & ".join(
        [str(s + 1) for s in np.where((x == 1) | (x == "1") | (x == "Checked"))[0]]), axis=1)
    # Dictionary for number to level translation
    column_translation = dict(zip(num, levels))
    # Replace numbers with levels
    for word, replacement in column_translation.items():
        values = [s.replace(word, replacement) for s in values]
    # Add to dataframe
    data[base_name] = values
    # Drop old columns
    if drop:
        data.drop(cols, axis=1, inplace=True)
    # Return data
    return data


def calc_egfr(df, age="age", serum_creatinine="creatinine_s",
              cystatin_c="cystatin_c_s", bun="bun", height="height", sex="sex",
              male="Male", female="Female", alpha=0.5):
    import pandas as pd
    import numpy as np
    # Make a copy of the dataframe
    data = df.copy()
    # Format input
    sex_data = data[sex].copy()
    sex_data.replace({male: "M", female: "F", "Other": np.nan, "": np.nan}, inplace=True)
    qcr = np.floor(data[age])
    serum_creatinine = pd.to_numeric(data[serum_creatinine], errors="coerce")
    cystatin_c = pd.to_numeric(data[cystatin_c], errors="coerce")
    height = pd.to_numeric(data[height], errors="coerce")
    bun = pd.to_numeric(data[bun], errors="coerce")
    age_data = pd.to_numeric(data[age], errors="coerce")

    # Younger participants
    qcr.replace({8: 0.46, 9: 0.49, 10: 0.51, 11: 0.53,
                12: 0.57, 13: 0.59, 14: 0.61}, inplace=True)
    # Females
    qcr.loc[(qcr == 15) & (sex_data == "F")] = 0.64
    qcr.loc[(qcr == 16) & (sex_data == "F")] = 0.67
    qcr.loc[(qcr == 17) & (sex_data == "F")] = 0.69
    qcr.loc[(qcr == 18) & (sex_data == "F")] = 0.69
    qcr.loc[(qcr >= 19) & (sex_data == "F")] = 0.70
    # Males
    qcr.loc[(qcr == 15) & (sex_data == "M")] = 0.72
    qcr.loc[(qcr == 16) & (sex_data == "M")] = 0.78
    qcr.loc[(qcr == 17) & (sex_data == "M")] = 0.82
    qcr.loc[(qcr == 18) & (sex_data == "M")] = 0.85
    qcr.loc[(qcr == 19) & (sex_data == "M")] = 0.88
    qcr.loc[(qcr > 19) & (sex_data == "M")] = 0.90

    # FAS and other equations
    eGFR_fas_cr = 107.3 / (serum_creatinine / qcr)
    f1 = serum_creatinine / qcr
    f2 = 1 - alpha
    f3 = cystatin_c / 0.82
    eGFR_fas_cr_cysc = 107.3 / ((0.5 * f1) + (f2 * f3))
    eGFR_Zap = (507.76 * np.exp((0.003 * height)) /
                ((cystatin_c ** 0.635) * ((serum_creatinine * 88.4) ** 0.547)))
    m = pd.to_numeric(sex_data.replace({"M": 1, "F": 0}))
    
    eGFR_Schwartz = 39.1 * ((height / serum_creatinine) ** 0.516) * ((1.8 / cystatin_c) ** 0.294) * ((30 / bun) ** 0.169) * \
        (1.099 ** m) * ((height / 1.4) ** 0.188)
    eGFR_bedside_Schwartz = (41.3 * (height / 100)) / serum_creatinine
    f = pd.to_numeric(sex_data.replace({"M": 0, "F": 1}))
    a = pd.to_numeric(sex_data.replace({"M": -0.302, "F": -0.241}))
    k = pd.to_numeric(sex_data.replace({"M": 0.9, "F": 0.7}))
    eGFR_CKD_epi = 142 * (np.minimum(serum_creatinine / k, 1)**a) * \
        (np.maximum(serum_creatinine / k, 1)**-1.200) * \
        (0.9938**age_data) * (1.012 * f + (1 - f))

    # CKiD U25 eGFR (Pierce CB, Muñoz A, Ng DK, Warady BA, Furth SL, Schwartz GJ. “Age- and sex-dependent clinical equations to estimate glomerular filtration rates in children and young adults with chronic kidney disease.” Kidney International. 2021;99(4):948-956.)
    def kappa_creat(age, sex):
        if sex == "F":
            if age < 12:
                return 36.1 * (1.008**(age - 12))
            elif age < 18:
                return 36.1 * (1.023**(age - 12))
            else:
                return 41.4
        else:
            if age < 12:
                return 39.0 * (1.008**(age - 12))
            elif age < 18:
                return 39.0 * (1.045**(age - 12))
            else:
                return 50.8

    def kappa_cys(age, sex):
        if sex == "F":
            if age < 12:
                return 79.9 * (1.004**(age - 12))
            elif age < 18:
                return 79.9 * (0.974**(age - 12))
            else:
                return 68.3
        else:
            if age < 15:
                return 87.2 * (1.011**(age - 15))
            elif age < 18:
                return 87.2 * (0.960**(age - 15))
            else:
                return 77.1

    kappa_cr = [kappa_creat(a, s) for a, s in zip(age_data, sex_data)]
    kappa_cysc = [kappa_cys(a, s) for a, s in zip(age_data, sex_data)]
    kappa_cr = pd.Series(kappa_cr, index=df.index)
    kappa_cysc = pd.Series(kappa_cysc, index=df.index)
    eGFR_CKiD_U25_Creat = kappa_cr * ((height/100) / serum_creatinine)
    eGFR_CKiD_U25_CystatinC = kappa_cysc * (1 / cystatin_c)
    eGFR_CKiD_U25_avg = (eGFR_CKiD_U25_Creat + eGFR_CKiD_U25_CystatinC) / 2

    # Combine results
    egfr = pd.concat(
        [eGFR_Schwartz, eGFR_bedside_Schwartz, eGFR_Zap, eGFR_fas_cr,
         eGFR_fas_cr_cysc, eGFR_CKD_epi,
         eGFR_CKiD_U25_Creat, eGFR_CKiD_U25_CystatinC, eGFR_CKiD_U25_avg], axis=1)

    egfr.columns = ["eGFR_Schwartz", "eGFR_bedside_Schwartz", "eGFR_Zap",
                    "eGFR_fas_cr", "eGFR_fas_cr_cysc", "eGFR_CKD_epi",
                    "eGFR_CKiD_U25_Creat", "eGFR_CKiD_U25_CystatinC", "eGFR_CKiD_U25_avg"]

    data = pd.concat([data, egfr], axis=1)
    
    # Return
    return data
  
def add_id_column(df, study_name):
    mrns = df.loc[df['study'] == study_name, ['mrn', 'record_id']].drop_duplicates('mrn')
    id_map = dict(zip(mrns['mrn'], mrns['record_id']))
    df[f'{study_name.lower().replace("-", "")}_id'] = df['mrn'].map(id_map)   
  

def create_study_id_columns(harmonized):
    import os
    study_list = [('CASPER', 'casper_id'), ('COFFEE', 'coffee_id'), ('CROCODILE','croc_id'), ('IMPROVE','improve_id'), ('PENGUIN','penguin_id'), ('RENAL-HEIR','rh_id'), ('RENAL-HEIRitage','rh2_id'), ('PANTHER','panther_id'), ('PANDA','panda_id'), ('ATTEMPT','attempt_id'), ('RPC2','rpc2_id'), ('SWEETHEART', 'swth_id')]

    for study, id in study_list: 
        study_mrns = harmonized.loc[harmonized['study'] == study, ['mrn', 'record_id']]
        study_id_map = dict(zip(study_mrns['mrn'], study_mrns['record_id']))
        harmonized[id] = harmonized.apply(lambda row: study_id_map[row['mrn']] if row['mrn'] in study_id_map else '', axis=1)

    # ---- Create and Save ID Linkage Matrix ----
    linkage_df = (
        harmonized
        .dropna(subset=['mrn', 'study', 'record_id'])
        .drop_duplicates(subset=['mrn', 'study'])
        .pivot(index='mrn', columns='study', values='record_id')
        .fillna('')
        .reset_index()
    )
    import getpass
    user = getpass.getuser()  # safer than os.getlogin(), works in more environments

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
    output_path = (base_data_path + "Data Harmonization/id_linkage_matrix.csv")
    linkage_df.to_csv(output_path, index=False)
    return harmonized