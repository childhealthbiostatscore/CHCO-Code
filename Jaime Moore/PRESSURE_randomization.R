library(readxl)
library(blockrand)
library(dplyr)

set.seed(3654)

# obs
obs <- rep("OBS", 100)
obs <- as.data.frame(obs)
obs <- blockrand(n=nrow(obs), levels = c(0:1), block.sizes = 1)
obs$treatment_char <- NA
obs$treatment_char <- ifelse(obs$treatment==1, "SEMAGLUTIDE", "USUAL_CARE")
obs$intinclusion7 <- "OBS"

# existing patient pool
existing <- rep("EXISTING", 100)
existing <- as.data.frame(existing)
existing <- blockrand(n=nrow(existing), levels = c(0:1), block.sizes = 1)
existing$treatment_char <- NA
existing$treatment_char <- ifelse(existing$treatment==1, "SEMAGLUTIDE", "USUAL_CARE")
existing$intinclusion7 <- "EXISTING"

# combine
rand <- rbind(obs, existing)

# check
table(rand$treatment_char)
table(rand$intinclusion7, rand$treatment_char)

# output
write.csv(rand, 
          "/Volumes/dept/SOM/Peds/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Moore/PRESSURE/Randomization/PRESSURE_randomization.csv",
          row.names = F)




