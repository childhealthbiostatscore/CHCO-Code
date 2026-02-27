# PANTHER labs compile for UACR confirmation
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
library(janitor)
library(readr)
library(lubridate)

lab_path <- file.path(root_path, "PANTHER/Data_Raw/Labs")

files <- list.files(lab_path, pattern = "\\.csv$", full.names = TRUE)

labs_df <- files %>%
  set_names(basename(.)) %>%
  map_dfr(
    ~ read_csv(
      .x,
      col_types = cols(
        MRN = col_character(),
        `Collection Date` = col_character(),
        `Collection Date...12` = col_character(),
        `Collection Date...13` = col_character(),
        `Order Comments` = col_character(),
        .default = col_guess()
      ),
      name_repair = "unique",
      show_col_types = FALSE
    ),
    .id = "source_file"
  ) %>%
  clean_names()

# Data cleaning
urine_albumin_tests <- c("MALBUR", "MALBUT")
urine_creatinine_tests <- c("Creatinine Urine Spot", "Creatinine-Urine")

labs_df_clean <- labs_df %>%
  dplyr::mutate(
    collection_date_combined = coalesce(as_date(mdy_hm(collection_date)), 
                                        as_date(mdy(collection_date_12)), 
                                        as_date(mdy(collection_date_13))),
    collection_month = month(collection_date_combined),
    collection_day   = day(collection_date_combined),
    collection_year  = year(collection_date_combined),
    test_name = case_when(test_name %in% urine_creatinine_tests ~ "urine_cr",
                          T ~ test_name),
    lloq = str_extract(result_text, "(?<=<)\\s*\\d+\\.?\\d*") %>%
      str_trim() %>%
      as.numeric(),
    result_numeric = case_when(!is.na(lloq) ~ lloq/2,
                               result_numeric < -99 ~ NA,
                               T ~ result_numeric)) %>%
  select(-c(collection_date, collection_date_12, 
            collection_date_13, visit2)) %>%
  distinct(mrn, collection_date_combined, test_name, result_numeric, result_text,
           .keep_all = T)  %>%
  mutate(
    test_name = if_else(
      mrn == "1280776" & test_name == "MALBUR" & result_numeric == 210,
      "MALBUR_flag",
      test_name
    )
  )

# filter to UACR related vars
# unique(labs_df_clean$test_name)
urine_labs <- labs_df_clean %>%
  filter(test_name %in% c(urine_albumin_tests, "urine_cr"))

# create wide dataset
id_cols <- c(
  "source_file", "protocol_number", "mrn", "sex",
  "collection_month", "collection_year", "collection_day"
)

urine_wide <- urine_labs %>%
  pivot_wider(
    id_cols = all_of(id_cols),
    names_from = test_name,
    values_from = c(result_numeric, units, result_text),
    names_sep = "_"
  ) %>%
  group_by(mrn, collection_month, collection_year) %>%
  fill(result_numeric_MALBUT, result_numeric_MALBUR, result_numeric_urine_cr,
       .direction = "updown") %>%
  ungroup() %>%
  arrange(mrn, collection_year, collection_month, collection_day) %>%
  dplyr::mutate(
    uacr_mg_g = 100 * (result_numeric_MALBUR / result_numeric_urine_cr)
  ) %>%
  distinct(mrn, uacr_mg_g, .keep_all = T) %>%
  dplyr::group_by(mrn) %>%
  dplyr::mutate(visitn = min(collection_year) - collection_year,
         visit = case_when(visitn == 0 ~ "baseline",
                           visitn == -1 ~ "year_1",
                           visitn == -2 ~ "year_2",
                           visitn == -3 ~ "year_2"))


# pull in PANTHER from harm dataset to compare
harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na.strings = "")
dat <- harm_dat %>% filter((study == "PANTHER")) %>%
  mutate(visit = case_when(visit == "screening" ~ "baseline", T ~ visit),
         mrn = as.character(mrn)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(mrn, visit))

dat_uacr <- dat %>%
  dplyr::select(mrn, record_id, visit, group, acr_u,
                microalbumin_u, creatinine_u, u24_mab) %>%
  right_join(urine_wide %>% 
               dplyr::select(mrn, visit, acr_u_lab = uacr_mg_g, 
                             u24_mab_lab = result_numeric_MALBUT,
                             microalbumin_u_lab = result_numeric_MALBUR, 
                             creatinine_u_lab = result_numeric_urine_cr, 
                             collection_year, collection_month, collection_day)) %>%
  arrange(record_id, collection_year) %>%
  dplyr::mutate(uacr_diff = round(acr_u_lab - acr_u, 1),
                u24_mab_diff = round(u24_mab_lab - u24_mab, 1))

# write.csv(dat_uacr, file.path(root_path, "PANTHER/Data_Cleaned/panther_uacr_lab_comparison.csv"), na = "", row.names = F)

dat_mab <- dat_uacr %>%
  filter(u24_mab_diff != 0 | (is.na(u24_mab) & !is.na(u24_mab_lab))) %>%
  dplyr::select(mrn, collection_year, collection_month, collection_day, record_id, visit, u24_mab, u24_mab_lab, u24_mab_diff)
write.csv(dat_mab, file.path(root_path, "PANTHER/Data_Cleaned/panther_u24_mab_lab_comparison.csv"), na = "", row.names = F)
