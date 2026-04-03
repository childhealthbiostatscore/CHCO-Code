library(Hmisc)
library(plyr)


label_harmonized <- function(data, dict) {
  # Dictionary file should have variable names as col names with corresponding description as one row.
  dict_list <- (t(dict))
  dict <- unlist(dict)
  names(dict) <- rownames(dict_list)
  Hmisc::label(data) = as.list(dict[match(names(data), names(dict))])
  data = as.data.frame(data)
  return(data)
}

label_harmonized_raw <- function(data, dict) {
  # Dictionary file should have two columns (variable name and description).
  dict <- setNames(data.frame(t(dict[ , - 1])), dict[ , 1]) # Transpose
  dict <- dict[intersect(names(data), names(dict))]
  dict[setdiff(names(data), names(dict))] <- ""
  dict_list <- (t(dict))
  dict <- unlist(dict)
  names(dict) <- rownames(dict_list)
  Hmisc::label(data) = as.list(dict[match(names(data), names(dict))])
  data = as.data.frame(data)
  return(data)
}

label_harmonized_dict <- function(data, dict) {
  dict <- dict %>% dplyr::select(-units, -notes)
  dict <- setNames(data.frame(t(dict[ , - 1])), dict[ , 1]) # Transpose
  dict <- dict[intersect(names(data), names(dict))]
  dict[setdiff(names(data), names(dict))] <- ""
  dict = as.data.frame(dict)
  return(dict)
}