#Collapsing function
collapse_by_id_visit <- function(df) {
  # Identify grouping and categorical variables
  group_vars <- c("record_id", "visit")
  categorical_vars <- c("study", "procedure", "group", "sex")
  
  # Get numeric variables (excluding grouping and categorical)
  all_vars <- names(df)
  numeric_vars <- setdiff(all_vars, c(group_vars, categorical_vars))
  
  df_collapsed <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      # For categorical variables, take first non-missing value
      across(all_of(categorical_vars[categorical_vars %in% all_vars]),
             ~dplyr::first(na.omit(.))),
      # For numeric variables, take mean of non-missing values
      across(all_of(numeric_vars[numeric_vars %in% all_vars]),
             ~if(all(is.na(.))) NA_real_ else mean(., na.rm = TRUE)),
      .groups = "drop"
    )
  
  return(df_collapsed)
}
# Define a custom rendering function without median
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)" = sprintf("%s (%s)", MEAN, SD)))
}