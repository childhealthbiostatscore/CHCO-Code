library(jsonlite)
library(dplyr)
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

harm_dat <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/harmonized_dataset.csv"), na = "")
uuid_map <- read.csv(file.path(root_path, "Data Harmonization/Data Clean/mrn_uuid_map.csv"))

bx_dat <- harm_dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>%
  filter(!is.na(kit_id)) %>%
  dplyr::select(mrn, kit_id)

cardiac_mri <- harm_dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>%
  filter(!is.na(ultra_id) | !is.na(swht_id)) %>%
  dplyr::select(-rvr) %>%
  filter(
    if_any(
      matches("^(lv|rv)"),
      ~ !is.na(.x)
    )
  ) %>%
  dplyr::select(mrn, record_id, swht_id, ultra_id, rh_id, rh2_id, improve_id,
                starts_with("lv"), starts_with("rv")) 

cardiac_bx <- left_join(cardiac_mri, bx_dat, by = "mrn") %>%
  filter(!is.na(kit_id)) %>%
  distinct(kit_id, .keep_all = T) %>%
  left_join(uuid_map) %>%
  dplyr::select(-mrn)

write.csv(cardiac_bx, file.path(root_path, "ULTRA_SWEETHEART/Data Clean/cardiac_mri_bx_overlap.csv"), row.names = F, na = "")
