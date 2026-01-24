
library(dplyr)
library(purrr)
library(jsonlite)
library(aws.s3)
library(digest)
library(readxl)

user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
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

# data pull
harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na.strings = "")

dat_subset <- harm_dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit))


# read in biopsy master spreadsheet (updated real time)
biopsy_master <- read_excel("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/Biopsy master tracker.xlsx", skip = 2) 
colnames(biopsy_master) <- gsub(" ", "_", colnames(biopsy_master))  # Replace spaces with underscores
colnames(biopsy_master) <- gsub("/", "_", colnames(biopsy_master))  # Replace slashes with underscores
colnames(biopsy_master) <- gsub("\\(", "_", colnames(biopsy_master))  # Replace opening parentheses with underscores
colnames(biopsy_master) <- gsub("\\)", "", colnames(biopsy_master))  # Remove closing parentheses
biopsy_master <- biopsy_master %>%
  dplyr::rename(record_id = Study_ID,
                visit = Visit_ID) %>%
  mutate(visit = case_when(Study == "ATTEMPT" ~ visit,
                           Study == "REMODEL" ~ visit, 
                           Study == "IMPROVE" & visit == "4M" ~ "3_months_post_surgery",
                           visit == "12M" ~ "12_months_post_surgery",
                           visit == "Follow-up" ~ "follow-up",
                           T ~ "baseline"),
         record_id = case_when(startsWith(Coenroll_ID__Same_visit_only, "RPC") & !is.na(Coenroll_ID__Same_visit_only) ~ Coenroll_ID__Same_visit_only,
                               T ~ record_id),
         Study = case_when(startsWith(Coenroll_ID__Same_visit_only, "RPC") & !is.na(Coenroll_ID__Same_visit_only) ~ "RPC2",
                           T ~ Study))
biopsy_completed <- biopsy_master %>%
  filter(scRNA_status == "Complete")


# sphyg + scRNA
sphyg_dat <- dat_subset %>%
  filter(pwv_succ == 1) %>%
  filter(record_id %in% biopsy_completed$record_id) %>%
  mutate(pwv = "sphygmocor") %>%
  dplyr::select(record_id, visit, pwv)

# CMRI + scRNA
cmri_dat <- dat_subset %>%
  filter(!is.na(lvedv)) %>%
  filter(paste0(record_id,visit) %in% paste0(biopsy_completed$record_id, biopsy_completed$visit)) %>%
  mutate(pwv = "cmr") %>%
  dplyr::select(record_id, visit, pwv)

pwv_scrna <- rbind(sphyg_dat, cmri_dat)
