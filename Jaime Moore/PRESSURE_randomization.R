library(readxl)
library(blockrand)
library(dplyr)

set.seed(3654)

# notes: need to combine obs and existing into one variable
# make sure that variable and treatment assignment are named the same as in REDCap


# obs
obs <- rep("OBS", 100)
obs <- as.data.frame(obs)
obs <- blockrand(n=nrow(obs), levels = c(0:1), block.sizes = 1)
obs$treatment_char <- NA
obs$treatment_char <- ifelse(obs$treatment==1, "SEMAGLUTIDE", "USUAL_CARE")

# existing patient pool
existing <- rep("EXISTING", 100)
existing <- as.data.frame(existing)
existing <- blockrand(n=nrow(existing), levels = c(0:1), block.sizes = 1)
existing$treatment_char <- NA
existing$treatment_char <- ifelse(existing$treatment==1, "SEMAGLUTIDE", "USUAL_CARE")

# combine
rand <- rbind(obs, existing)

# check
table(rand$treatment_char)
table(rand$What.was.your.biological.sex.assigned.at.birth., rand$treatment_char)

# output
write.csv(rand, 
          "",
          row.names = F)





# read in participant info
participants <- read.csv("/Volumes/pylell/Hill/Randomization/23-2151 Cohort 1 Randomization.csv")

# stratify by biological sex and randomize to counseling or no counseling
males <- participants %>% filter(What.was.your.biological.sex.assigned.at.birth. == "Male")
females <- participants %>% filter(What.was.your.biological.sex.assigned.at.birth. == "Female")

# randomize males
males_rand <- blockrand(n=nrow(males), levels = c(0:1), block.sizes = 1)
males_rand <- males_rand[1:nrow(males),]
males <- cbind(males, males_rand$treatment)
males$treatment_char <- ifelse(males$`males_rand$treatment`==1, "Counseling", "No Counseling")
males$`males_rand$treatment` <- NULL

# randomize females
females_rand <- blockrand(n=nrow(females), levels = c(0:1), block.sizes = 1)
females_rand <- females_rand[1:nrow(females),]
females <- cbind(females, females_rand$treatment)
females$treatment_char <- ifelse(females$`females_rand$treatment`==1, "Counseling", "No Counseling")
females$`females_rand$treatment` <- NULL

# combine
rand <- rbind(males, females)

# check
table(rand$treatment_char)
table(rand$What.was.your.biological.sex.assigned.at.birth., rand$treatment_char)

# output
write.csv(rand, 
          "/Volumes/pylell/Hill/Randomization/23-2151 Cohort 1 Randomization with Assignments.csv",
          row.names = F)
