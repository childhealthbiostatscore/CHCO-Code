#1. Set Up ----
##a. Load Libraries & Directores----
#Libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Directories
# computer <- "mac studio"
computer <- "mac laptop"
if (computer == "mac studio") {
  user <- Sys.info()[["user"]]
  if (user == "hhampson") {
    dir.dat <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive")
    dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/PANTHER/Results")
    git_path <- "Users/hhampson/CHCO-Code/Petter Bjornstad"
  } else if (user == "shivaniramesh") {
    dir.dat <- c("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive")
    dir.results <- c("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver Pet IR/Results")
  } else if (user == "choiyej") {
    dir.dat <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
    dir.results <- paste0(dir.dat, "/PANTHER/Results")
    git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  } else if (user == "pylell") {
    dir.dat <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
    dir.results <- paste0(dir.dat, "/PANTHER/Results")
    git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  }
} else {
  user <- Sys.info()[["user"]]
  if (user == "hhampson") {
    dir.dat <- c("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive")
    dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/PANTHER/Results")
    git_path <- "Users/hhampson/CHCO-Code/Petter Bjornstad"
  } else {
    dir.dat <- c("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive")
    dir.results <- c("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver Pet IR/Results")
  }
}

# source(file.path(git_path, "PANTHER PFAS/Libraries.R"))
source("Libraries.R")

##b. Load Data----
#Load Metadata
dat <- read.csv(fs::path(dir.dat,"Data Harmonization","Data Clean","harmonized_dataset.csv"),na.strings = "")
lod <- read.csv(fs::path(dir.dat,"/PANTHER/Data_Raw/PFAS Emory 2025/8- PFAS Quantification/PFAS_LODs.csv"))


##c. Clean & Format Data----
pfas <- c("N.EtFOSAA","N.MeFOSAA","PFBA","PFDA","PFDoA","PFHpA","PFHxA","PFNA","PFBS","PFPeA","PFTeDA","PFTrDA","PFUnA","PFBS","PFDoS","PFHps","PFHxS","PFNS","PFOS","PFOSA","PFPeAS")
dat_collapsed1 <- dat %>%
  filter(study=="PANTHER") %>% 
  filter(visit=="baseline" | visit=="screening") %>% 
  # dplyr::select(record_id, date, rh2_id, visit,mrn, group, age, sex, bmi, hba1c,albuminuria_cat,microalbumin_u,starts_with("egfr"),cystatin_c_s,starts_with("creat"),all_of(pfas)) %>% 
  dplyr::select(record_id, dob, date, rh2_id, visit, mrn, group, age, sex, bmi, hba1c,albuminuria_cat,microalbumin_u,eGFR_CKiD_U25_avg,cystatin_c_s,creatinine_s,creatinine_u) %>% 
  # group_by(record_id, visit) %>%
  # tidyr::fill(date, .direction = "updown") %>% ungroup() %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>% 
  mutate(visit="baseline_screening")
dat_collapsed2 <- dat %>%
  filter(study=="PANTHER") %>% 
  filter(visit=="year_1" | visit=="year_2") %>% 
  # dplyr::select(record_id, date, rh2_id, visit,mrn, group, age, sex, bmi, hba1c,albuminuria_cat,microalbumin_u,starts_with("egfr"),cystatin_c_s,starts_with("creat"),all_of(pfas)) %>% 
  dplyr::select(record_id, dob, date, rh2_id, visit,mrn, group, age, sex, bmi, hba1c,albuminuria_cat,microalbumin_u,eGFR_CKiD_U25_avg,cystatin_c_s,creatinine_s,creatinine_u) %>% 
  # group_by(record_id, visit) %>%
  # tidyr::fill(date, .direction = "updown") %>% ungroup() %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) 
meta <- rbind(dat_collapsed1,dat_collapsed2)

pfas_dat <- dat %>% 
  filter(study=="PANTHER") %>% 
  filter(visit=="baseline") %>% 
  dplyr::select("record_id",all_of(pfas))%>% #127 panther baseline participants
  filter(!is.na(N.EtFOSAA))

#Impute PFAS Data
t_lod <- t(lod)
colnames(t_lod) <- t_lod[1,]
t_lod <- as.data.frame(t_lod)
t_lod <- t_lod[-1,]
lod <- t_lod
rm(t_lod)
lod <- lod %>% 
  mutate(across(everything(), ~as.numeric(.)))
colnames(lod) <- str_replace_all(colnames(lod),"-",".")

# pfas
N.EtFOSAA <- length(which(pfas_dat$N.EtFOSAA < lod$`N-EtFOSAA`))/length(unique(pfas_dat$record_id))
N.EtFOSAA < 2/3 #False, exclude
N.MeFOSAA <- length(which(pfas_dat$N.MeFOSAA < lod$`N-MeFOSAA`))/length(unique(pfas_dat$record_id))
N.MeFOSAA < 2/3 #True, impute 
PFBA <- length(which(pfas_dat$PFBA < lod$`PFBA`))/length(unique(pfas_dat$record_id))
PFBA < 2/3 #False, exclude
PFDA <- length(which(pfas_dat$PFDA < lod$`PFDA`))/length(unique(pfas_dat$record_id))
PFDA < 2/3 #True, impute
PFDoA <- length(which(pfas_dat$PFDoA < lod$`PFDoA`))/length(unique(pfas_dat$record_id))
PFDoA < 2/3 #False, exclude
PFHpA <- length(which(pfas_dat$PFHpA < lod$`PFHpA`))/length(unique(pfas_dat$record_id))
PFHpA < 2/3 #True, impute
PFHxA <- length(which(pfas_dat$PFHxA < lod$`PFHxA`))/length(unique(pfas_dat$record_id))
PFHxA < 2/3 #False, exclude
PFNA <- length(which(pfas_dat$PFNA < lod$`PFNA`))/length(unique(pfas_dat$record_id))
PFNA < 2/3 #True, impute
PFBS <- length(which(pfas_dat$PFBS < lod$`PFBS`))/length(unique(pfas_dat$record_id))
PFBS < 2/3 #True, impute
PFPeA <- length(which(pfas_dat$PFPeA < lod$`PFPeA`))/length(unique(pfas_dat$record_id))
PFPeA < 2/3 #False, exclude
PFTeDA <- length(which(pfas_dat$PFTeDA < lod$`PFTeDA`))/length(unique(pfas_dat$record_id))
PFTeDA < 2/3 #False, exclude
PFTrDA <- length(which(pfas_dat$PFTrDA < lod$`PFTrDA`))/length(unique(pfas_dat$record_id))
PFTrDA < 2/3 #False, exclude
PFUnA <- length(which(pfas_dat$PFUnA < lod$`PFUnA`))/length(unique(pfas_dat$record_id))
PFUnA < 2/3 #False, exclude
PFBS <- length(which(pfas_dat$PFBS < lod$`PFBS`))/length(unique(pfas_dat$record_id))
PFBS < 2/3 #True, impute
PFDoS <- length(which(pfas_dat$PFDoS < lod$`PFDoS`))/length(unique(pfas_dat$record_id))
PFDoS < 2/3 #False, exclude
PFHps <- length(which(pfas_dat$PFHps < lod$`PFHps`))/length(unique(pfas_dat$record_id))
PFHps < 2/3 #True, impute
PFHxS <- length(which(pfas_dat$PFHxS < lod$`PFHxS`))/length(unique(pfas_dat$record_id))
PFHxS < 2/3 #True, impute
PFNS <- length(which(pfas_dat$PFNS < lod$`PFNS`))/length(unique(pfas_dat$record_id))
PFNA < 2/3 #True, impute
PFOS <- length(which(pfas_dat$PFOS < lod$`PFOS`))/length(unique(pfas_dat$record_id))
PFOS < 2/3 #True, impute
PFOSA <- length(which(pfas_dat$PFOSA < lod$`PFOSA`))/length(unique(pfas_dat$record_id))
PFOSA < 2/3 #False, exclude
PFPeAS <- length(which(pfas_dat$PFPeAS < lod$`PFPeAS`))/length(unique(pfas_dat$record_id))
PFPeAS < 2/3 #True, impute

pfas_to_impute <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA", "PFBS", 
                    "PFHps", "PFHxS", "PFNS", "PFOS", "PFPeAS")
pfas_imputed <- pfas_dat %>% 
  dplyr::select(c("record_id",all_of(pfas_to_impute)))
for (pfas in pfas_to_impute){
  pfas_imputed[[pfas]] <- ifelse(pfas_imputed[[pfas]]<lod[[pfas]],(lod[[pfas]]/sqrt(2)),pfas_imputed[[pfas]])
}

#Save formatted and imputed pfas data
# saveRDS(pfas_imputed,fs::path(dir.dat,"PFAS_Data_Imputed.rds"))

# pfas_to_impute
table1(~N.MeFOSAA+PFDA+PFHpA+PFNA+PFBS+PFHps+PFHxS+PFNS+PFOS+PFPeAS,data=pfas_imputed,render.continuous = render.geometric)

#Format metadata
# meta <- dat_collapsed
dat <- tidylog::left_join(pfas_imputed,meta,by="record_id")
dat_baseline <- dat %>% 
  filter(visit=="baseline_screening")
# table1(~age+sex+bmi+hba1c+microalbumin_u+creatinine_s+creatinine_u+eGFR_CKiD_U25_avg | group,data=dat)



