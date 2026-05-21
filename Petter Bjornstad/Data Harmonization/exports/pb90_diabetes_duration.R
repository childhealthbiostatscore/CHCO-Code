# obtain diabetes duration for PB90 per Tomas's request
library(jsonlite)
library(dplyr)
library(kableExtra)
library(knitr)
library(purrr)
library(tidyr)
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else {
  stop("Unknown user: please specify root path for this user.")
}
Sys.setenv(
  AWS_ACCESS_KEY_ID     = keys$MY_ACCESS_KEY,
  AWS_SECRET_ACCESS_KEY = keys$MY_SECRET_KEY,
  AWS_DEFAULT_REGION    = "",
  AWS_REGION            = "",
  AWS_S3_ENDPOINT       = "s3.kopah.uw.edu"
)

harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na = "") 

dat <- harm_dat %>%
  dplyr::filter(procedure == "kidney_biopsy") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit))

pb90_meta <- s3readRDS("data_clean/pb90_meta.rds", "scrna",
          region = "")

pb90_meta <- pb90_meta %>%
  left_join(subset(dat, select = c(record_id, visit, diabetes_duration)))

pb90_meta_subset <- pb90_meta %>%
  distinct(record_id, visit, group, diabetes_duration)
write.csv(pb90_meta_subset, file.path(root_path, "Data Harmonization/Data Exports/pb90_diabetes_duration.csv"), na = "", row.names = F)
