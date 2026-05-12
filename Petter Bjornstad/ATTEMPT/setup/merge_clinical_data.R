# READ IN AND MERGE ATTEMPT DATA SENT BY ANTOINE
# AUTHOR: LAURA PYLE

library(readxl)
library(dplyr)

# DEMOGRAPHICS
demo <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_Demographics")

# ANTHROPOMETRICS
anthro <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_Anthropometrics")
anthro <- anthro %>% select(-c(site, date_visit))

# FAMILY HISTORY
famhx <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_MedFamSmokingHx")
famhx <- famhx %>% select(-c(site, date_visit))

# DIABETES MANAGEMENT
dm <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_DiabetesManagement")
dm <- dm %>% select(-c(site, date_visit))

# GLUCOSE MANAGEMENT
gm <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_GlucoseMonitoring")
gm <- gm %>% select(-c(site, date_visit))

# URINE LABS
urinelabs <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_LocalUrineLabs")
urinelabs <- urinelabs %>% select(-c(site, date_visit))

# URINE 24 H
urine24 <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_Urine24h")
urine24 <- urine24 %>% select(-c(site))

# URINE EMU
urineemu <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_UrineEMU")
urineemu <- urineemu %>% select(-c(site, date_visit))

# LOCAL BLOOD LABS
localblood <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_LocalBloodLabs")
localblood <- localblood %>% select(-c(site, date_visit))

# CENTRAL BLOOD LABS
centralblood <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_CentralBloodLabs")
centralblood <- centralblood %>% select(-c(site, date_visit))

# COMPLIANCE
compliance <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_Compliance")
compliance <- compliance %>% select(-c(site, date_visit))

# EGFR
egfr <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_eGFR")
egfr <- egfr %>% select(-c(site, date_visit))

# MGFR
mgfr <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_mGFR")
mgfr <- mgfr %>% select(-c(site, date_visit))

# BOLD MRI
boldmri <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Raw/ATTEMPT_AnalyticalDatasets_Denver.xlsx',
                  sheet = "ATTEMPT_BoldMRI")
boldmri <- boldmri %>% select(-c(site, date_visit))

# merge by ID and visit
attempt_clinical <- full_join(demo, anthro, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, anthro, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, famhx, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, dm, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, gm, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, urinelabs, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, urine24, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, urineemu, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, localblood, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, centralblood, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, compliance, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, egfr, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, mgfr, by = c("subject_id", "visit"))
attempt_clinical <- full_join(attempt_clinical, boldmri, by = c("subject_id", "visit"))

attempt_clinical <- attempt_clinical %>% arrange(subject_id, visit)
