# ATTEMPT request from Jaime (JXMM) on 9/4/25

# Could I also ask you to send the clinical information of the donors and the blood SomaScan expression? Sometimes, we are interested in correlations of proteins or correlations with clinical information. Given the rare population and the new version of proteins, the dataset is unique.
# 
# Please, make sure the data is fully anonymised. We just need to map the clinical information to the individual protein expression, so please create pseudo labels to remove the donor ID.

library(dplyr)
library(purrr)
library(jsonlite)
library(aws.s3)
library(digest)

user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

## Create an S3 client
keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# SOMA data pull
harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/soma_harmonized_dataset.csv"), na.strings = "")

# ATTEMPT clinical data pull
attempt_dat <- s3readRDS(object = "cleaned_data/attempt_dat.rds", bucket = "attempt", region = "")

# anonymize IDs
set.seed(123)
unique_ids <- unique(c(harm_dat$record_id, attempt_dat$subject_id))
mapping <- setNames(sample(seq_along(unique_ids)), unique_ids)
harm_dat$anon_id <- mapping[as.character(harm_dat$record_id)]
attempt_dat$anon_id <- mapping[as.character(attempt_dat$subject_id)]

soma_subset <- harm_dat %>%
  filter(study == "ATTEMPT") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(urine_creat_proteomics > 1) %>%
  filter(!is.na(seq.10000.28)) %>%
  dplyr::select(anon_id, visit, matches("^seq\\..*\\d$")) 

View(soma_subset)

# SOMA analyte info data pull
load(file.path(root_path, "Data Harmonization/Combined SomaScan/analytes_2.Rdata"))
View(analytes_attempt)

attempt_dat_subset <- attempt_dat %>%
  mutate(site = case_when(site == "Denver" ~ "A", site == "Toronto" ~ "B", site == "London" ~ "C"),
         age = floor(age)) %>%
  dplyr::select(anon_id, visit, treatment, site, age, sex, mgfr_jodal, mgfr_jodal_bsa, 
                hba1c, emu_urine_acr_mean, weight, height, bmi, sbp, dbp,
                cgm_tir, diabetes_dx_duration,
                cholesterol_serum_mmoll, triglycerides_serum_mmoll, 
                ldl_serum_mmoll, hdl_serum_mmoll) %>%
  filter(anon_id %in% soma_subset$anon_id)

# save all files
write.csv(soma_subset, file.path(root_path, "ATTEMPT/Data Clean/Plasma proteomics subset/JXMM_091225_plasma_soma.csv"), row.names = F)
write.csv(analytes_attempt, file.path(root_path, "ATTEMPT/Data Clean/Plasma proteomics subset/JXMM_091225_plasma_soma_analytes.csv"), row.names = F)
write.csv(attempt_dat_subset, file.path(root_path, "ATTEMPT/Data Clean/Plasma proteomics subset/JXMM_091225_plasma_soma_clin.csv"), row.names = F)
