# load libraries
library(tidyverse)
library(dplyr)
library(REDCapR)
library(readxl)
library(jsonlite)
library(Hmisc)
library(aws.s3)
# Set up environment for Kopah
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
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

## Create an S3 client
keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

dat <- read.csv(file.path(root_path, "/Data Harmonization/Data Clean/harmonized_dataset.csv"), na.strings = c(" ", "", "-9999",-9999))

# read in biopsy master spreadsheet (updated real time)
biopsy_master <- read_excel("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/Biopsy master tracker.xlsx", skip = 2) 
colnames(biopsy_master) <- gsub(" ", "_", colnames(biopsy_master))  # Replace spaces with underscores
colnames(biopsy_master) <- gsub("/", "_", colnames(biopsy_master))  # Replace slashes with underscores
colnames(biopsy_master) <- gsub("\\(", "_", colnames(biopsy_master))  # Replace opening parentheses with underscores
colnames(biopsy_master) <- gsub("\\)", "", colnames(biopsy_master))  # Remove closing parentheses
biopsy_master <- biopsy_master %>%
  dplyr::rename(record_id =Study_ID ,
                visit = Visit_ID) %>%
  mutate(visit = case_when(Study == "ATTEMPT" ~ visit,
                           Study == "REMODEL" ~ visit, 
                           Study == "IMPROVE" & visit == "4M" ~ "3_months_post_surgery",
                           visit == "12M" ~ "12_months_post_surgery",
                           visit == "Follow-up" ~ "follow-up",
                           T ~ "baseline"),
         record_id = case_when(startsWith(record_id, "RPC") & !is.na(Coenroll_ID__Same_visit_only) ~ Coenroll_ID__Same_visit_only,
                               T ~ record_id)) # RH2-60-T and RH2-48-T coenrolled into RPC2 as of 02/03/25 data pull

# subset biopsy master to PANDA
biopsy_master_panda <- biopsy_master %>%
  filter(Study == "PANDA") %>%
  filter(!is.na(record_id))

# subset harm dat to PANDA
dat_panda <- dat %>%
  # filter(procedure == "kidney_biopsy") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(mrn, date)) %>%
  filter(!is.na(panda_id)) %>%
  filter(study != "ATTEMPT")
  
# MRI in PANDA
dat_panda_mri <- dat_panda %>%
  filter(
    if_all(
      c(avg_pcascl, total_kidney_volume_ml, avg_c_r2, avg_k_r2),
      ~ is.na(.x)
    )
  ) %>%
  filter(study == "PANDA") %>%
  dplyr::select(record_id, panda_id, croc_id, visit, avg_pcascl, total_kidney_volume_ml, avg_c_r2, avg_k_r2)

# record ID in master spreadsheet but not in REDCap
biopsy_master_panda$record_id[biopsy_master_panda$record_id %nin% dat_panda$record_id]

# record ID in REDCap but not master spreadsheet
dat_panda$record_id[dat_panda$record_id %nin% biopsy_master_panda$record_id & dat_panda$record_id %nin% biopsy_master_panda$Coenroll_ID__Same_visit_only]

# list of PB90
pb90_meta <- s3readRDS("data_clean/pb90_meta.rds", bucket = "scrna", region = "") 

pb90_meta_unique <- pb90_meta %>%
  distinct(kit_id, cryostor_id, .keep_all = T) %>%
  mutate(kit_id = gsub("KI", "KL", kit_id),
         kit_id = gsub("Kl", "KL", kit_id))

# list of ATTEMPT seurat
attempt_meta <- s3readRDS("cleaned_data/attempt_clean_so_metadata.rds", bucket = "attempt", region = "")

attempt_meta_unique <- attempt_meta %>%
  distinct(bx_kit_id, Cryostor_ID, .keep_all = T) %>%
  mutate(kit_id = gsub("KI", "KL", bx_kit_id),
         kit_id = gsub("Kl", "KL", kit_id),
         cryostor_id = Cryostor_ID)

# list of PANDA
dat_panda_ids <- dat_panda %>%
  dplyr::select(ends_with("_id"), study) %>%
  mutate(kit_id = gsub("KI", "KL", kit_id),
         kit_id = gsub("Kl", "KL", kit_id))

length(unique(dat_panda$cryostor_id))
length(unique(dat_panda$kit_id))

dat_panda_ids <- dat_panda_ids %>%
  mutate(run_yn = case_when(kit_id %in% pb90_meta_unique$kit_id | cryostor_id %in% pb90_meta_unique$cryostor_id~ "pb90", 
                            kit_id %in% attempt_meta_unique$kit_id | cryostor_id %in% attempt_meta_unique$cryostor_id~ "attempt", 
                            T ~ "n")) %>%
  filter(!is.na(cryostor_id)) %>%
  filter(kit_id != "KL-0029530") 

table(dat_panda_ids$run_yn)

dat_panda_ids %>%
  filter(study == "PANDA" & !is.na(croc_id))

dat_panda_ids %>% # Sample collected for diagnostics only
  filter(study == "PANDA", run_yn == "n")


write.csv(dat_panda_ids, file.path(root_path, "PANDA/Data_Cleaned/panda_outstanding_scrna.csv"), row.names = F, na = "")


# IDs to run
dat_panda_ids_run <- dat_panda_ids %>%
  filter(run_yn == "n" | (study == "PANDA" & run_yn == "pb90"))

write.csv(dat_panda_ids_run, file.path(root_path, "PANDA/Data_Cleaned/panda_outstanding_scrna_torun_010625.csv"), row.names = F, na = "")
