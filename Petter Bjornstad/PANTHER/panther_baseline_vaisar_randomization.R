library(togolab)
library(readxl)
library(dplyr)

#keys <- read.csv("/Users/lpyle/Library/CloudStorage/OneDrive-UW/lpyle@uw.edu_accessKeys.csv")
#togo_setup_s3(keys = keys) 
togo_paths()

# read in samples that were sent to Tomas
samples <- read_xlsx('/Users/lpyle/Library/CloudStorage/OneDrive-SharedLibraries-UW/bjornstad_pyle_tommerdahl_lab - Documents/Data Science/Projects/P031_PANTHER_baseline_yc/resources/PANTHER Samples for Tomas Vaisar Lab.xlsx')
samples <- samples[,1]
samples <- samples[1:97,]

# read in harmonized data to get clinical characteristics
harm <- togo_load_harmonized()
panther <- harm %>% filter(study == "PANTHER")
panther_baseline <- panther %>% filter(visit == "baseline")
panther_screening <- panther %>% filter(visit == "screening")
panther_keep <- panther_screening %>% select(record_id, sex, group_risk, tan_fgd, tan_fph, tan_tveq, tan_mgd, tan_mph, breast_tanner)
panther_keep <- panther_keep %>% dplyr::mutate(
              tanner_stage_comp = coalesce(tan_fgd, tan_fph, tan_tveq, tan_mgd, tan_mph, breast_tanner),
              tanner_stage_comp_panther = case_when(tanner_stage_comp > 3 ~ 4, T~ tanner_stage_comp))
