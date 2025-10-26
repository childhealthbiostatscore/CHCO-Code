# Run this script to pull IDs of people with MRI in the BIC
# specify user for paths
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

library(tidyverse)
library(purrr)
library(Hmisc)

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

# Get two levels (keyword-matched folders and their immediate subfolders)
bic_folders <- get_filtered_subfolders(folder_path)

bic_folders_clean <- bic_folders %>%
       filter(level1_folder %nin% c("Bjornstad_GZBU_231527", "Kendrick_BC_161572", "Kelsey_HIP_070988")) %>%
       mutate(record_id = gsub("ATTEMPT[_.]|ATTEMPTR.|_192947|_171874|_170820|_170802|_161752|_191282|_180704|_18-704|_200277|_212999|_1611752|_213019", 
                               "", 
                               level2_folder)) %>%
  filter(!grepl("[._]twix", record_id, ignore.case = T)) %>%
  dplyr::select(record_id, level1_path)

write.csv(bic_folders_clean, file.path(root_path, "Data Harmonization/Data Clean/bic_folders_list.csv"), row.names = F)
