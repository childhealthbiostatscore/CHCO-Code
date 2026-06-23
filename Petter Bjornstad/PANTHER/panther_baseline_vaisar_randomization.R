library(togolab)
library(readxl)

togo_paths(setup_s3 = FALSE)

# read in samples that were sent to Tomas
samples <- read_xlsx('/Users/lpyle/Library/CloudStorage/OneDrive-SharedLibraries-UW/bjornstad_pyle_tommerdahl_lab - Documents/Data Science/Projects/P031_PANTHER_baseline_yc/resources/PANTHER Samples for Tomas Vaisar Lab.xlsx')
samples <- samples[,1]
samples <- samples[1:97,]

# read in harmonized data to get clinical characteristics