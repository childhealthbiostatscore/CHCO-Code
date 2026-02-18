# Run this script to pull IDs of people with MRI in the BIC
# specify user for paths
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "shivaniramesh") { # hyak version
  root_path <- "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive"
  git_path <- "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

library(tidyverse)
library(purrr)
library(Hmisc)
library(stringr)

get_filtered_subfolders <- function(root_path, keywords = c("bjornstad", "nadeau", 
                                                            "kendrick", "kelsey")) {
  #' Get top-level folder, first-level subfolders matching keywords, and their subfolders
  #' 
  #' @param root_path Path to the root folder to scan
  #' @param keywords Character vector of keywords to filter first-level subfolders
  #' @return A dataframe with columns 'top_folder', 'level1_folder', 'level2_folder', and paths
  
  # Get the top folder name
  top_folder_name <- basename(root_path)
  
  # Get all items in the root folder
  all_items <- list.files(root_path, full.names = TRUE)
  
  # Filter to only directories
  all_dirs <- all_items[dir.exists(all_items)]
  
  # Get just the folder names (not full paths)
  folder_names <- basename(all_dirs)
  
  # Create pattern for grep - combine keywords with OR operator
  pattern <- paste(keywords, collapse = "|")
  
  # Filter folders that contain any of the keywords (case-insensitive)
  matching_indices <- grep(pattern, folder_names, ignore.case = TRUE)
  
  # Initialize results list
  results_list <- list()
  
  if (length(matching_indices) > 0) {
    # For each matching first-level folder
    for (i in matching_indices) {
      level1_folder <- folder_names[i]
      level1_path <- all_dirs[i]
      
      # Get all subfolders in this first-level folder
      level2_items <- list.files(level1_path, full.names = TRUE)
      level2_dirs <- level2_items[dir.exists(level2_items)]
      
      if (length(level2_dirs) > 0) {
        # If there are subfolders, add each as a row
        for (level2_path in level2_dirs) {
          level2_folder <- basename(level2_path)
          results_list[[length(results_list) + 1]] <- data.frame(
            top_folder = top_folder_name,
            level1_folder = level1_folder,
            level2_folder = level2_folder,
            level1_path = level1_path,
            level2_path = level2_path,
            stringsAsFactors = FALSE
          )
        }
      } else {
        # If no subfolders, still include the level1 folder with NA for level2
        results_list[[length(results_list) + 1]] <- data.frame(
          top_folder = top_folder_name,
          level1_folder = level1_folder,
          level2_folder = NA,
          level1_path = level1_path,
          level2_path = NA,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Combine all results into one dataframe
    df <- bind_rows(results_list)
    
  } else {
    # Return empty dataframe with proper structure if no matches
    df <- data.frame(
      top_folder = character(),
      level1_folder = character(),
      level2_folder = character(),
      level1_path = character(),
      level2_path = character(),
      stringsAsFactors = FALSE
    )
    message("No subfolders found matching the keywords")
  }
  
  return(df)
}

# Set your folder path
folder_path <- "/Volumes/bic-server"

extract_num <- function(x) {
  hits <- stringr::str_extract_all(x, rx)[[1]]
  if (length(hits) == 0) return(NA_character_)
  hits <- hits[nchar(hits) == max(nchar(hits))]
  dplyr::last(hits)  # tie-breaker
}

# Get two levels (keyword-matched folders and their immediate subfolders)
bic_folders <- get_filtered_subfolders(folder_path)

bic_folders_clean <- bic_folders %>%
  filter(level1_folder %nin% c("Bjornstad_GZBU_231527", "Kendrick_BC_161572", "Kelsey_HIP_070988")) %>%
  dplyr::mutate(bic_id = gsub("ATTEMPT[_.]|ATTEMPTR.|_192947|_171874|_170820|_170802|_161752|_191282|_180704|_18-704|_200277|_212999|_1611752|_213019", 
                              "", 
                              level2_folder)) %>%
  filter(!grepl("[._]twix", bic_id, ignore.case = T)) %>%
  dplyr::mutate(study = case_when(bic_id == "RH2.04.O" ~ "RENAL-HEIRitage",
                                  grepl("ATTEMPT", level1_path) ~ "ATTEMPT",
                                  grepl("CFE", level1_path) ~ "COFFEE",
                                  grepl("CRC", level1_path) ~ "CROCODILE",
                                  grepl("CS", level1_path) ~ "CASPER",
                                  grepl("Bjornstad_HEIR", level1_path) ~ "RENAL-HEIR",
                                  grepl("IMPROVE", level1_path) ~ "IMPROVE",
                                  grepl("PENGUIN", level1_path) ~ "PENGUIN",
                                  grepl("RH2", level1_path) ~ "RENAL-HEIRitage",
                                  grepl("PANTHER", level1_path) ~ "PANTHER",
                                  grepl("RPC2", level1_path) ~ "RPC2",
                                  grepl("PANDA", level1_path) ~ "PANDA",
  ),
  
  record_number = case_when(study == "ATTEMPT" ~ sub("^([0-9]+)[._].*", "\\1", bic_id),
                            study == "PANDA" ~ as.character(as.numeric(str_extract(bic_id, "(?<=[._-])[0-9]+(?=$|[._-])")) + 100),
                            T ~ sapply(stringr::str_extract_all(bic_id, "(?<=[._-])[0-9]+(?![0-9])"),
                                       function(h) if (length(h)) tail(h[nchar(h) == max(nchar(h))], 1) else NA_character_)
  ),
  record_number = as.numeric(record_number),
  record_id = case_when(study == "ATTEMPT" ~ sub("^([0-9]+)[._].*", "\\1", bic_id),
                        study %in% c("COFFEE", "CROCODILE", "CASPER") ~ gsub("[_.]", "-", bic_id),
                        study == "RENAL-HEIR" ~ sub("-0$", "-O", gsub("[_.]", "-", bic_id)),
                        study == "IMPROVE" ~ sub("\\_[123]$", "", gsub("IT2D", "IT", gsub("[.]", "_", bic_id))),
                        study == "PENGUIN" ~ sub("PENGUIN[._]", "PEN-", bic_id),
                        study == "RENAL-HEIRitage" ~ sub("-0$", "-O", gsub("[_.]", "-", gsub("RH.2", "RH2", gsub("\\.[123]$|\\_ABD$|\\_BRAIN$", "", bic_id)))),
                        study == "PANTHER" ~ gsub("PANTHER", "PAN", gsub("[._]", "-", gsub("\\.[123]$", "", bic_id))),
                        study == "RPC2" ~ ifelse(grepl("\\.[123]$|\\.BL$", bic_id), sub("RPC2", "RPC", gsub("\\.", "-", sub("\\.[123]$|\\.BL$", "", bic_id))), sub("RPC2[.]", "RPC-", bic_id)),
                        study == "PANDA" ~ sub("PANDA.", "PNDA-1", bic_id)
  ),
  visit_id = case_when(bic_id %in% c("IT2D_14_2", "IT2D_15_2", "IT2D.16.2") ~ "3",
                       study %in% c("ATTEMPT", "RPC2", "PANTHER", "IMPROVE") ~ ifelse(grepl("[._-][1-3]$", bic_id), sub(".*[._-]([1-3])$", "\\1", bic_id), "1"),
                       study %in% c("COFFEE", "CROCODILE", "CASPER", "PANDA", "RENAL-HEIRitage", "PENGUIN", "RENAL-HEIR") ~ "1"
  ),
  visit_id = as.numeric(visit_id),
  ) %>%
  filter(!is.na(study))  %>%
  filter(!is.na(record_number)) %>%
  filter(record_number != 0) %>%
  dplyr::select(bic_id, record_id_bic = record_id, record_number, visit_id, study)

write.csv(bic_folders_clean, file.path(root_path, "Data Harmonization/Data Clean/MRI/bic_folders_list.csv"), row.names = F)

# pull in data from harmonized dataset
bic_folders_clean <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/MRI/bic_folders_list.csv"))
harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na = "")

mri_dat <- harm_dat %>%
  dplyr::mutate(visit = case_when(study == "ATTEMPT" & visit == "screening" ~ "baseline", T ~ visit)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  dplyr::mutate(record_number = case_when(study == "ATTEMPT" ~ as.character(sub("^([0-9]+)[._].*", "\\1", record_id)),
                                          T ~ stringr::str_extract(record_id, "(?<=[._\\- ])[0-9]+(?=$|[._\\- ])")),
                record_number = as.numeric(record_number),
                visit_id = case_when(visit == "v2_gfr_mri_arm_1" ~ 1,
                                     visit == "v7_gfr_mri_arm_1" ~ 2,
                                     visit == "baseline" ~ 1,
                                     visit == "4_months_post" ~ 2,
                                     visit == "3_months_post_surgery" ~ 2,
                                     visit == "12_months_post_surgery" ~ 3,
                                     visit == "year_1" ~ 2,
                                     visit == "year_2" ~ 3,
                                     visit == "post_treatment" ~ 2),
                data_in_redcap = case_when(dplyr::if_all(
                  c(avg_c_r2, avg_k_r2, avg_m_r2,
                    total_kidney_volume_ml,
                    avg_c_t1, avg_k_t1,
                    avg_c_adc, avg_pcascl), is.na) ~ FALSE,
                  TRUE ~ TRUE
                )
  ) %>%
  filter(!is.na(visit_id)) %>%
  dplyr::select(record_id, 
                record_number, study, procedure, visit, visit_id, data_in_redcap,
                avg_c_r2, avg_k_r2, avg_m_r2, total_kidney_volume_ml, avg_c_t1, avg_k_t1, avg_c_adc, avg_pcascl) 


panda_uw_mri_ids <- harm_dat %>%
  filter(procedure == "mri" & (grepl("PNDA 2", record_id) | grepl("PNDA-2", record_id) | grepl("PNDA_2", record_id))) %>%
  filter(is.na(avg_c_r2)) %>%
  filter(!is.na(date)) %>%
  pull(record_id)

panther_manual_ids <- c("PAN_202_T",
                        "PAN_203_O", # for 1 year follow up
                        "PAN_204_T",
                        "PAN_205_O",
                        "PAN_207_O",
                        "PAN_208_O",
                        "PAN_209_O",
                        "PAN_210_O",
                        "PAN_211_O",
                        "PAN_208_O"
)

mri_dat_combined <- full_join(mri_dat, bic_folders_clean, 
                              by = join_by(record_number, study, visit_id)) %>%
  filter(!(study == "PANDA" & visit != "baseline")) %>%
  mutate(data_in_bic = case_when(record_id == "PAN_203_O" & visit == "baseline" ~ FALSE,
                                 record_id == "PAN_203_O" & visit == "year_1" ~ TRUE,
                                 record_id %in% panther_manual_ids & visit == "baseline" ~ TRUE,
                                 record_id %in% panda_uw_mri_ids ~ TRUE,
                                 is.na(bic_id) ~ FALSE, 
                                 bic_id == "IT2D_14" ~ FALSE, # comment in REDCap - pt unable to be scanned
                                 T ~ TRUE),
         record_id = coalesce(record_id, record_id_bic),
         data_in_redcap = case_when(is.na(data_in_redcap) ~ FALSE, 
                                    T ~ data_in_redcap),
         status = case_when(data_in_bic & !data_in_redcap ~ "Awaiting analysis",
                            data_in_bic & data_in_redcap ~ "Complete",
                            data_in_redcap & !data_in_bic ~ "Complete",
                            !data_in_bic ~ "N/A")) %>%
  dplyr::select(record_id, bic_id, visit_id, study, visit, data_in_redcap, data_in_bic, status) %>%
  arrange(status)

write.csv(mri_dat_combined, file.path(root_path, "Data Harmonization/Data Clean/MRI/mri_data_compiled.csv"), row.names = F, na = "")

table(subset(mri_dat_combined, status == "Awaiting analysis")$study)
