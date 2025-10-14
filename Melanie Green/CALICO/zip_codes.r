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
