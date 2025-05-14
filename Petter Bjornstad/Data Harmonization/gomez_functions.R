library(dplyr)

data <- data %>%
  mutate(
    erpf_raw_plasma_seconds = erpf_raw_plasma / 60,
    gfr_raw_plasma_seconds = gfr_raw_plasma / 60,
    
    # Filtration Fraction
    ff = gfr_raw_plasma / erpf_raw_plasma,
    
    # Kfg for group (T1D/T2D: 0.1012, Control: 0.1733)
    kfg = case_when(
      group == "1" ~ 0.1012,
      group == "2" ~ 0.1733,
      TRUE ~ NA_real_
    ),
    
    # Filtration pressure across glomerular capillaries
    deltapf = (gfr_raw_plasma / 60) / kfg,
    
    # Plasma protein mean concentration
    cm = (bl_tot_protein / ff) * log(1 / (1 - ff)),
    
    # Pi G (Oncotic pressure)
    pg = 5 * (cm - 2),
    
    # Glomerular Pressure
    glomerular_pressure = pg + deltapf + 10,
    
    # Renal Blood Flow
    rbf = erpf_raw_plasma / (1 - hct_210 / 100),
    rbf_seconds = erpf_raw_plasma_seconds / (1 - hct_210 / 100),
    
    # Renal Vascular Resistance (mmHg*l^-1*min^-1)
    rvr = map / rbf,
    
    # Efferent Arteriolar Resistance
    re = (gfr_raw_plasma_seconds / (kfg * (rbf_seconds - gfr_raw_plasma_seconds))) * 1328,
    
    # Afferent Arteriolar Resistance
    ra = ((map - glomerular_pressure) / rbf_seconds) * 1328
  )