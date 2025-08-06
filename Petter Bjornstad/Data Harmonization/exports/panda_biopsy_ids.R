# PANDA Biopsy N (with unique biopsy, not coenrolled into other studies like ATTEMPT/CROC)
# Some have repeat biopsies that are unique to PANDA to capture

library(reticulate)
library(dplyr)
library(tidyr)
library(purrr)
library(REDCapR)
library(readxl)

harm_dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

panda <- harm_dat %>% filter(!is.na(panda_id)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  arrange(mrn) %>%
  filter(!is.na(kit_id))

panda_unique_bx <- panda |>
  add_count(kit_id)|>
    filter(n == 1) |>
    distinct() |>
    filter(study == "PANDA") %>%
    mutate(site = case_when(grepl("^PNDA[- ]1\\d{2}$", record_id) ~ "CU",
                            grepl("^PNDA[- ]2\\d{2}$", record_id) ~ "UW"),
           panda_id = gsub(" ", "-", panda_id)) %>%
  dplyr::select(record_id, mrn, visit, ends_with("_id"), study)

# panda_followup <- panda %>%
#   dplyr::select(record_id, panda_id, p1_raw_m, p2_raw_m) %>%
#   filter(!is.na(p1_raw_m))

# pull in master biopsy spreadsheet

biopsy_master <- read_excel("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/Biopsy master tracker.xlsx", skip = 2) 
colnames(biopsy_master) <- gsub(" ", "_", colnames(biopsy_master))  # Replace spaces with underscores
colnames(biopsy_master) <- gsub("/", "_", colnames(biopsy_master))  # Replace slashes with underscores
colnames(biopsy_master) <- gsub("\\(", "_", colnames(biopsy_master))  # Replace opening parentheses with underscores
colnames(biopsy_master) <- gsub("\\)", "", colnames(biopsy_master))  # Remove closing parentheses

biopsy_master <- biopsy_master %>%
  dplyr::rename(record_id=Study_ID,
         visit=Visit_ID) %>%
  mutate(visit = case_when(Study == "ATTEMPT" ~ visit,
                           Study == "REMODEL" ~ visit, 
                           Study == "IMPROVE" & visit == "4M" ~ "3_months_post_surgery",
                           visit == "12M" ~ "12_months_post_surgery",
                           visit == "Follow-up" ~ "follow-up",
                           T ~ "baseline"),
         record_id = case_when(startsWith(record_id, "RPC") & !is.na(Coenroll_ID__Same_visit_only) ~ Coenroll_ID__Same_visit_only,
                               T ~ record_id)) # RH2-60-T and RH2-48-T coenrolled into RPC2 as of 02/03/25 data pull

biopsy_master_panda <- biopsy_master %>%
  dplyr::select(record_id, Study, Biopsy_date, visit, Shipped_Y_N, scRNA_status, ends_with("_ID")) %>%
  filter(!is.na(record_id)) %>%
  # filter(Shipped_Y_N == "Yes") %>%
  # filter(scRNA_status != "No sample") %>%
  filter(Study == "PANDA") %>%
  mutate(sequenced = case_when(scRNA_status == "Complete" ~ "Yes",
                               T~ "No"),
         panda_id = gsub(" ", "-", record_id))

# biopsy_master_panda$Kit_ID[biopsy_master_panda$Kit_ID %nin% panda_unique_bx$kit_id]
# panda_unique_bx$kit_id[panda_unique_bx$kit_id %nin% biopsy_master_panda$Kit_ID]

# PNDA 205 missing cryostor

master_panda_sub <- biopsy_master_panda %>%
  dplyr::select(panda_id, Kit_ID, Cryostor_ID) %>%
  dplyr::rename(kit_id = Kit_ID,
                cryostor_id = Cryostor_ID) %>%
  mutate(croc_id = NA,
         attempt_id = NA)

harm_panda_sub <- panda_unique_bx %>%
  dplyr::select(panda_id, croc_id, attempt_id, kit_id, cryostor_id)

combined_panda <- rbind(master_panda_sub, harm_panda_sub) %>%
  group_by(panda_id) %>%
  fill(attempt_id, croc_id, kit_id, cryostor_id, .direction = "updown") %>%
  filter(is.na(attempt_id)) %>%
  distinct(panda_id, .keep_all = T) %>%
  mutate(site = case_when(grepl("^PNDA[- ]1\\d{2}$", panda_id) ~ "Colorado",
                          grepl("^PNDA[- ]2\\d{2}$", panda_id) ~ "Washington"),
         # cryostor_id = case_when(panda_id == "PNDA-219" ~ "S-2504-01205",
         #                         T ~ cryostor_id),
         # kit_id = case_when(panda_id == "PNDA-204" ~ "KL-00325650",
         #                    T ~ kit_id)
         ) %>%
  # filter(!is.na(cryostor_id)) %>%
  dplyr::select(-attempt_id, croc_id, site)

# pnda14 <- combined_panda %>%
#   filter(startsWith(panda_id, "PNDA-2"))

# pnda17 <- combined_panda %>%
#   filter(startsWith(panda_id, "PNDA-2"))
# 
# pnda17$panda_id[pnda17$panda_id %nin% pnda14$panda_id]

write.csv(combined_panda, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/PANDA/Data_Cleaned/panda_unique_biopsy_ids.csv", row.names = F, na = "")

# crocodile Kit IDs


croc <- harm_dat %>% filter(record_id %in% combined_panda$croc_id[!is.na(combined_panda$croc_id)]) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  arrange(mrn) %>%
  filter(!is.na(kit_id))

croc_unique_bx <- croc |>
  mutate(kit_id_croc = kit_id, cryostor_id_croc = cryostor_id, croc_id = record_id) |> 
  dplyr::select(croc_id, kit_id_croc, cryostor_id_croc) 

croc_panda_combined <- left_join(combined_panda, croc_unique_bx)

write.csv(croc_panda_combined, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/PANDA/Data_Cleaned/panda_unique_biopsy_ids.csv", row.names = F, na = "")
