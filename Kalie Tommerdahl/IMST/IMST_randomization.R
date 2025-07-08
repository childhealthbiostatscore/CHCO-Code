library(readxl)
library(blockrand)
library(dplyr)

set.seed(3654)

# male, hypertension
stratum <- rep("MALE, HYPERTENSION", 100)
male_hypertension <- as.data.frame(stratum)
temp <- blockrand(n=nrow(male_hypertension), levels = c(1:2), block.sizes = 1)
male_hypertension <- cbind(male_hypertension, temp)
male_hypertension$treatment_char <- NA
male_hypertension$treatment_char <- ifelse(male_hypertension$treatment==1, "IMST", "SHAM")
male_hypertension$scr_sex <- 1
male_hypertension$scr_hyper_130 <- 1

# male, no hypertension
stratum <- rep("MALE, NO HYPERTENSION", 100)
male_no_hypertension <- as.data.frame(stratum)
temp <- blockrand(n=nrow(male_no_hypertension), levels = c(1:2), block.sizes = 1)
male_no_hypertension <- cbind(male_no_hypertension, temp)
male_no_hypertension$treatment_char <- NA
male_no_hypertension$treatment_char <- ifelse(male_no_hypertension$treatment==1, "IMST", "SHAM")
male_no_hypertension$scr_sex <- 1
male_no_hypertension$scr_hyper_130 <- 2

# female, hypertension
stratum <- rep("FEMALE, HYPERTENSION", 100)
female_hypertension <- as.data.frame(stratum)
temp <- blockrand(n=nrow(female_hypertension), levels = c(1:2), block.sizes = 1)
female_hypertension <- cbind(female_hypertension, temp)
female_hypertension$treatment_char <- NA
female_hypertension$treatment_char <- ifelse(female_hypertension$treatment==1, "IMST", "SHAM")
female_hypertension$scr_sex <- 2
female_hypertension$scr_hyper_130 <- 1

# female, no hypertension
stratum <- rep("FEMALE, NO HYPERTENSION", 100)
female_no_hypertension <- as.data.frame(stratum)
temp <- blockrand(n=nrow(female_no_hypertension), levels = c(1:2), block.sizes = 1)
female_no_hypertension <- cbind(female_no_hypertension, temp)
female_no_hypertension$treatment_char <- NA
female_no_hypertension$treatment_char <- ifelse(female_no_hypertension$treatment==1, "IMST", "SHAM")
female_no_hypertension$scr_sex <- 2
female_no_hypertension$scr_hyper_130 <- 2

# other, hypertension
stratum <- rep("OTHER, HYPERTENSION", 100)
other_hypertension <- as.data.frame(stratum)
temp <- blockrand(n=nrow(other_hypertension), levels = c(1:2), block.sizes = 1)
other_hypertension <- cbind(other_hypertension, temp)
other_hypertension$treatment_char <- NA
other_hypertension$treatment_char <- ifelse(other_hypertension$treatment==1, "IMST", "SHAM")
other_hypertension$scr_sex <- 3
other_hypertension$scr_hyper_130 <- 1

# other, no hypertension
stratum <- rep("OTHER, NO HYPERTENSION", 100)
other_no_hypertension <- as.data.frame(stratum)
temp <- blockrand(n=nrow(other_no_hypertension), levels = c(1:2), block.sizes = 1)
other_no_hypertension <- cbind(other_no_hypertension, temp)
other_no_hypertension$treatment_char <- NA
other_no_hypertension$treatment_char <- ifelse(other_no_hypertension$treatment==1, "IMST", "SHAM")
other_no_hypertension$scr_sex <- 3
other_no_hypertension$scr_hyper_130 <- 2

# combine
rand <- rbind(male_hypertension, male_no_hypertension, female_hypertension, female_no_hypertension, other_hypertension, other_no_hypertension)

# check
table(rand$treatment_char)
table(rand$stratum, rand$treatment_char)

# format for redcap
rand$redcap_randomization_number <- ""
rand <- rand %>% select(redcap_randomization_number, treatment, scr_sex, scr_hyper_130)
colnames(rand) <- c("redcap_randomization_number", "redcap_randomization_group", 
                    "scr_sex", "scr_hyper_130")

# output
write.csv(rand, 
          "/Users/lpyle/Library/CloudStorage/OneDrive-UW/Tommerdahl/IMST/Randomization/IMST_randomization.csv",
          row.names = F)





