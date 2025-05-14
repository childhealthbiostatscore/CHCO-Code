# Functions for Gomez calculations

data["erpf_raw_plasma_seconds"] = data["erpf_raw_plasma"]/60
data["gfr_raw_plasma_seconds"] = data["gfr_raw_plasma"]/60

# Filtration Fraction
data["ff"] = data["gfr_raw_plasma"]/data["erpf_raw_plasma"] 

# Kfg for group (T1D/T2D kfg: 0.1012, Control kfg: 0.1733)
data["kfg"] = np.select([data["group"].eq("1"), data["group"].eq("2")], [0.1012, 0.1733]) 

# Filtration pressure across glomerular capillaries
data["deltapf"] = (data["gfr_raw_plasma"]/60)/data["kfg"] 

# Plasma protein mean concentration
data["cm"] = (data["bl_tot_protein"]/data["ff"])*np.log(1/(1-data["ff"])) 

# Pi G (Oncotic pressure)
data["pg"] = 5*(data["cm"]-2)

# Glomerular Pressure
data["glomerular_pressure"] = data["pg"] + data["deltapf"] + 10

# Renal Blood Flow
data["rbf"] = (data["erpf_raw_plasma"]) / (1 - data["hct_210"]/100)
data["rbf_seconds"] = (data["erpf_raw_plasma_seconds"]) / (1 - data["hct_210"]/100)

# Renal Vascular Resistance (mmHg*l^-1*min^-1)
data["rvr"] = data["map"] / data["rbf"]

# Efferent Arteriorlar Resistance 
data["re"] = (data["gfr_raw_plasma_seconds"]) / (data["kfg"] * (data["rbf_seconds"] - (data["gfr_raw_plasma_seconds"]))) * 1328

# Afferent Arteriolar Resistance
data["ra"] = ((data["map"] - data["glomerular_pressure"]) / data["rbf_seconds"]) * 1328
