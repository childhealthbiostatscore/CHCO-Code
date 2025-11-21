# Zipcodes
library(readxl)
setwd(
  "/Users/tim/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Vigers/BDC/Janet Snell-Bergeon/CALICO"
)
zips = read_excel("./Data_Raw/Zipcodes_SVI.xlsx")
zcta = read_excel("./Data_Raw/ZIP Code to ZCTA Crosswalk.xlsx")
svi = read.csv("./Data_Raw/SVI_2022_US_ZCTA.csv")
colnames(svi) = svi[1, ]
svi = svi[-1, ]
svi$ZCTBA = gsub("ZCTA5 ", "", svi$LOCATION)
# Get SVI for each zipcode
zips$ZCTBA = zcta$zcta[match(zips$`Zip Code`, zcta$ZIP_CODE)]
zips$`SVI Percentile` = svi$RPL_THEMES[match(zips$ZCTBA, svi$ZCTBA)]
# Write
write.csv(
  zips,
  file = "./Data_Clean/zip_code_svi.csv",
  row.names = FALSE,
  na = ""
)
# LARC matching
library(tidyverse)
# Load data
load("./Data_Clean/analysis_data.RData")
# Get ever LARC/EC and age at first LARC/EC
df <- df %>%
  arrange(record_number, cv_monthssincepcosdx) %>%
  group_by(record_number) %>%
  mutate(
    larc_ever = any(larc == "Yes"),
    age_first_larc = first(cv_age[larc == "Yes"]),
    ec_ever = any(ec == "Yes"),
    age_first_ec = first(cv_age[ec == "Yes"])
  )
df$larc_ever <- factor(df$larc_ever, levels = c(F, T), labels = c("No", "Yes"))
df$ec_ever <- factor(df$ec_ever, levels = c(F, T), labels = c("No", "Yes"))
# Just get the columns we need for matching, remove duplicates
match_df = df |>
  select(
    record_number,
    site,
    larc_ever,
    age_first_larc,
    ec_ever,
    age_first_ec,
    ethnicity
  ) |>
  distinct()
# If they have a LARC start date, use that as age otherwise use EC age
match_df$age = coalesce(match_df$age_first_larc, match_df$age_first_ec)
# Some people have a missing value for site, so get the letters at the start of
# the ID
match_df = match_df |> drop_na(age)
t = matchit(
  larc_ever ~ age + site,
  data = match_df,
  method = "optimal"
)
matched = data.frame(
  "LARC" = match_df[as.numeric(rownames(t$match.matrix)), "record_number"],
  "EC" = match_df[as.numeric(t$match.matrix[, 1]), "record_number"]
)
