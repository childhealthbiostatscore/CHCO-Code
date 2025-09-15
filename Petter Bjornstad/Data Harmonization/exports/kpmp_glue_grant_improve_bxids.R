# PB90 IMPROVE biopsy IDs for EO (for sequencing info for the KPMP glue grant)
library(dplyr)
library(purrr)
library(jsonlite)
library(aws.s3)

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

# pull in PB90 meta data
pb90_meta <- s3readRDS(bucket = "scrna", "data_clean/pb90_meta.rds", region = "")

pb90_improve <- pb90_meta %>% 
  filter(cohort == "IMPROVE") %>%
  distinct(kit_id, .keep_all = T)

improve_ids <- pb90_improve %>%
  select(record_id, visit, cryostor_id, kit_id)

write.csv(improve_ids, file.path(root_path, "Data Harmonization/Data Exports/kpmp_gg_improveIDs.csv"), row.names = F)

