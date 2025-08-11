
library(tidyverse)
library(dplyr)


harmonized <- read.csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW(2)/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv")  
dict <- read.csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW(2)/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/data_dictionary_master.csv")


processed_subset_harmonized <- function(data,
                               vars = NULL,                     # variables to include (NULL = all)
                               study_subset = NULL,              # vector of study names
                               numeric_fn = mean,                # numeric summarization/imputation function
                               non_numeric_fn = last,            # summarization function for non-numeric
                               group_by_vars = "record_id") {    # grouping variables
  
  # 1. Subset to variables + required columns
  required_cols <- c("study", group_by_vars)
  if (!is.null(vars)) {
    missing_vars <- setdiff(vars, colnames(data))
    if (length(missing_vars) > 0) {
      stop(paste("These variables are missing from dataset:", paste(missing_vars, collapse = ", ")))
    }
    data <- data %>% select(all_of(required_cols), all_of(vars))
  } else {
    data <- data %>% select(all_of(required_cols), everything())
  }
  
  # 2. Filter by study if requested
  if (!is.null(study_subset)) {
    data <- data %>% filter(study %in% study_subset)
  }
  
  # 3. Normalize missing values (convert "" or " " to NA)
  data <- data %>%
    mutate(across(everything(), ~ na_if(trimws(.), "")))
  
  # 4. Impute missing values before summarizing
  data <- data %>%
    mutate(
      across(where(is.numeric),
             ~ ifelse(is.na(.x), numeric_fn(.x, na.rm = TRUE), .x)),
      across(where(negate(is.numeric)),
             ~ {
               colname <- cur_column()
               if (grepl("_id$", colname)) {
                 .x  # keep NA
               } else {
                 ifelse(is.na(.x), "not recorded", .x)
               }
             })
    )
  
  # 5. Summarize
  summarized <- data %>%
    summarise(
      across(where(is.numeric),
             ~ numeric_fn(.x, na.rm = TRUE)),
      across(where(negate(is.numeric)),
             ~ non_numeric_fn(.x)),
      .by = all_of(group_by_vars)
    )
  
  return(summarized)
}


#result <- processed_subset_harmonized(
#  data = harmonized,
#  vars = c("casper_id", "croc_id", "mrn"),                 # NULL for all vars
#  study_subset = c("CASPER", "CROCODILE"),          # filter studies
#  numeric_fn = mean,                              # mean, median, min, max, first, last, etc.
#  non_numeric_fn = last,                            # summarization & imputation for non-numeric
#  group_by_vars = c("record_id")           # one or more grouping variables
#)

#print(result)
