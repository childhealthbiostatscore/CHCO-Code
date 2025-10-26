library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# specify user for paths
user <- Sys.info()[["user"]]
if (user == "laurapyle") {
  data_path <- "/Users/laurapyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/PBMCs"
  harm_path <- '/Users/laurapyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean'
} else if (user == "lpyle") {
  data_path <- "/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/PBMCs"
  harm_path <- '/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean'
} else if (user == "pylell") {
  data_path <- '/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/PBMCs'
  harm_path <- '/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean'
 } else {
  stop("Unknown user: please specify root path for this user.")
}

setwd(data_path)

convert_cell_count_to_numeric <- function(x) {
  # Vectorized function to handle NA and empty strings
  sapply(x, function(val) {
    # Return NA for missing values
    if (is.na(val) || val == "") return(NA_real_)
    
    # Convert to lowercase for easier matching
    val_lower <- tolower(val)
    
    # Handle "K" notation (e.g., "600K" = 600,000)
    if (str_detect(val_lower, "k")) {
      num <- as.numeric(str_extract(val_lower, "[0-9.]+"))
      return(num * 1000)
    }
    
    # Handle scientific notation with "e" (e.g., "1.25e6" = 1,250,000)
    if (str_detect(val_lower, "e[0-9]+")) {
      num <- as.numeric(str_extract(val_lower, "[0-9.]+e[0-9]+"))
      return(num)
    }
    
    # Handle "x 10^6" or "x10^6" notation
    if (str_detect(val_lower, "10\\^")) {
      # Extract the coefficient (number before "x")
      coef <- as.numeric(str_extract(val_lower, "[0-9.]+"))
      
      # Extract the exponent (number after "^")
      exponent <- as.numeric(str_extract(val_lower, "(?<=\\^)[0-9]+"))
      
      # Handle cases like "1.6 x 10^" where exponent is missing
      if (is.na(exponent)) return(NA_real_)
      
      return(coef * 10^exponent)
    }
    
    # If none of the above patterns match, try direct numeric conversion
    num <- as.numeric(str_extract(val_lower, "[0-9.]+"))
    return(num)
  })
}

# Function to reformat IDs
reformat_ids <- function(ids) {
  # Remove any trailing asterisks or extra whitespace
  ids <- trimws(gsub("\\*$", "", ids))
  
  # Apply reformatting rules
  reformatted <- sapply(ids, function(id) {
    if (grepl("^CRC", id)) {
      # Extract the number part
      num <- gsub("^CRC-?", "", id)
      return(paste0("CRC-", num))
      
    } else if (grepl("^IT2D", id)) {
      # Extract the number part
      num <- gsub("^IT2D-?", "", id)
      return(paste0("IT_", num))
      
    } else if (grepl("^RH", id)) {
      # Extract the number part (remove RH and any existing -T or T suffix)
      num <- gsub("^RH-?", "", id)
      num <- gsub("-?T$", "", num)
      return(paste0("RH-", num, "-T"))
      
    } else if (grepl("^RPC", id)) {
      # Extract the number part
      num <- gsub("^RPC-?", "", id)
      return(paste0("RPC-", num))
      
    } else {
      # If doesn't match any pattern, return as is
      return(id)
    }
  }, USE.NAMES = FALSE)
  
  return(reformatted)
}

assign_group_dplyr <- function(df_with_ids, df_lookup, id_col = "ID") {
  # Create a unified lookup table from all ID columns
  lookup_unified <- df_lookup %>%
    select(croc_id, improve_id, rh_id, rpc2_id, group) %>%
    pivot_longer(cols = c(croc_id, improve_id, rh_id, rpc2_id),
                 names_to = "id_type",
                 values_to = "id") %>%
    filter(!is.na(id)) %>%
    distinct(id, .keep_all = TRUE) %>%
    select(id, group)
  
  # Join with the main dataframe
  df_with_ids <- df_with_ids %>%
    left_join(lookup_unified, by = setNames("id", id_col))
  
  return(df_with_ids)
}

pbmc <- readxl::read_xlsx("./Bjornstad PBMC Manifest_17Jan2025-2.xlsx")
pbmc$ID_reformatted <- reformat_ids(pbmc$`Cell Line`)
pbmc_with_cell_count <- pbmc %>% filter(!is.na(`Alterations/Notes`))
pbmc_vial_count <- as.data.frame(table(pbmc_with_cell_count$`Cell Line`))
names(pbmc_vial_count) <- c("Cell Line", "Vial Count")
pbmc_with_cell_count <- full_join(pbmc_with_cell_count, pbmc_vial_count, by = "Cell Line")
pbmc_with_cell_count <- pbmc_with_cell_count %>%
  mutate(cell_count_numeric = convert_cell_count_to_numeric(`Alterations/Notes`))

pbmc_3x10e6 <- pbmc_with_cell_count %>% filter(cell_count_numeric > 3000000)
pbmc_3x10e6_gt1 <- pbmc_3x10e6 %>% filter(`Vial Count` > 1)


# now merge w/ harmonized data to get disease groups
filename <- "harmonized_dataset.csv"
full_path <- file.path(harm_path, filename)
harm_data <- read.csv(full_path)
#Sort by screening date
harm_data <- harm_data %>%  arrange(date_of_screen) 
collapsed_data <- harm_data	%>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))
pbmc_3x10e6_gt1 <- assign_group_dplyr(pbmc_3x10e6_gt1, collapsed_data, id_col = "ID_reformatted")

