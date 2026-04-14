library(dplyr)
library(stringr)

# specify user for paths
user <- Sys.info()[["user"]]
if (user == "laurapyle") {
  data_path <- "/Users/laurapyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward"
  github_path <- "/Users/laurapyle/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  root_path <- "/Users/laurapyle/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
} else if (user == "lpyle") {
  data_path <- "/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward"
  github_path <- "/Users/lpyle/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  root_path <- "/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
} else if (user == "pylell") {
  data_path <- "/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward"
  github_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  teen_labs_path <- "/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/Teen Labs/"
} else {
  stop("Unknown user: please specify root path for this user.")
}

setwd(data_path)

###############################################################
# Read in all raw data files                                  #
###############################################################

# TODAY forms
BASELINE  <- read.csv("./Clinical data/TODAY/BASELINE.csv")
CBL       <- read.csv("./Clinical data/TODAY/CBL.csv")
ADDCBL    <- read.csv("./Clinical data/TODAY/ADDCBL.csv")
PAT       <- read.csv("./Clinical data/TODAY/PAT.csv")
AGEBASE   <- read.csv("./Clinical data/TODAY/AGEBASE.csv")
PRIMOUT   <- read.csv("./Clinical data/TODAY/PRIMOUT.csv")
BW        <- read.csv("./Clinical data/TODAY/Birthweight.csv")
DXA       <- read.csv("./Clinical data/TODAY/DEXA.csv")
VISIT     <- read.csv("./Clinical data/TODAY/VISIT.csv")
comorb    <- read.csv("./Clinical data/COMORB.csv")

# TODAY2 forms
CBL_TODAY2    <- read.csv("./Clinical data/TODAY2/CBL.csv")
VISIT_TODAY2  <- read.csv("./Clinical data/TODAY2/VISIT.csv")
lipo          <- read.csv("./Clinical data/TODAY2/LIPO.csv")
tme           <- read.csv("./Clinical data/TODAY2/TME.csv")

# Adjudicated medical events (TODAY2)
AME           <- read.csv("./Clinical data/TODAY2/AME.csv")

# Quality of life (PedsQL) - TODAY
PEDSQLGC_TODAY <- read.csv("./Clinical data/TODAY/PEDSQLGC.csv")
PEDSQLGT_TODAY <- read.csv("./Clinical data/TODAY/PEDSQLGT.csv")

# Quality of life (PedsQL) - TODAY2
PEDSQLGA_TODAY2 <- read.csv("./Clinical data/TODAY2/PEDSQLGA.csv")
PEDSQLGY_TODAY2 <- read.csv("./Clinical data/TODAY2/PEDSQLGY.csv")
PEDSQLGT_TODAY2 <- read.csv("./Clinical data/TODAY2/PEDSQLGT.csv")

###############################################################
# COMORB dataset - outcomes with bariatric surgery censoring  #
###############################################################

comorb$MIC.OR.MAC <- ifelse(comorb$MAC==1 | comorb$MIC==1, 1,
                            ifelse(is.na(comorb$MAC) & is.na(comorb$MIC), NA, 0))
comorb <- comorb %>% mutate(DAYSTOMIC.OR.MAC = pmin(DAYSTOMIC, DAYSTOMAC))
comorb$releaseid <- comorb$RELEASEID
comorb$RELEASEID <- NULL

# Merge in bariatric surgery info
tme_bari <- tme %>% filter(TMETYPE==1) %>% select(RELEASEID, TMETYPE, DAYSTOTME)
colnames(tme_bari) <- c("releaseid", "TMETYPE", "DAYSTOTME")
comorb <- merge(comorb, tme_bari, by="releaseid", all.x=T, all.y=F)

# Censor all outcomes at time of bariatric surgery
# HTN
comorb <- comorb %>% mutate(DAYSTOHTN = case_when(
  HTN==0 ~ DAYSTOHTN,
  HTN==1 & is.na(DAYSTOTME) ~ DAYSTOHTN,
  HTN==1 & !is.na(DAYSTOTME) & DAYSTOHTN>DAYSTOTME ~ DAYSTOTME,
  HTN==1 & !is.na(DAYSTOTME) & DAYSTOHTN<=DAYSTOTME ~ DAYSTOHTN
))
comorb <- comorb %>% mutate(HTN = case_when(
  HTN==0 ~ HTN,
  HTN==1 & is.na(DAYSTOTME) ~ HTN,
  HTN==1 & !is.na(DAYSTOTME) & DAYSTOHTN>DAYSTOTME ~ as.integer(0),
  HTN==1 & !is.na(DAYSTOTME) & DAYSTOHTN<=DAYSTOTME ~ HTN
))
# LDL
comorb <- comorb %>% mutate(DAYSTOLDL = case_when(
  LDLDLP==0 ~ DAYSTOLDL,
  LDLDLP==1 & is.na(DAYSTOTME) ~ DAYSTOLDL,
  LDLDLP==1 & !is.na(DAYSTOTME) & DAYSTOLDL>DAYSTOTME ~ DAYSTOTME,
  LDLDLP==1 & !is.na(DAYSTOTME) & DAYSTOLDL<=DAYSTOTME ~ DAYSTOLDL
))
comorb <- comorb %>% mutate(LDLDLP = case_when(
  LDLDLP==0 ~ LDLDLP,
  LDLDLP==1 & is.na(DAYSTOTME) ~ LDLDLP,
  LDLDLP==1 & !is.na(DAYSTOTME) & DAYSTOLDL>DAYSTOTME ~ as.integer(0),
  LDLDLP==1 & !is.na(DAYSTOTME) & DAYSTOLDL<=DAYSTOTME ~ LDLDLP
))
# TG
comorb <- comorb %>% mutate(DAYSTOTG = case_when(
  TGDLP==0 ~ DAYSTOTG,
  TGDLP==1 & is.na(DAYSTOTME) ~ DAYSTOTG,
  TGDLP==1 & !is.na(DAYSTOTME) & DAYSTOTG>DAYSTOTME ~ DAYSTOTME,
  TGDLP==1 & !is.na(DAYSTOTME) & DAYSTOTG<=DAYSTOTME ~ DAYSTOTG
))
comorb <- comorb %>% mutate(TGDLP = case_when(
  TGDLP==0 ~ TGDLP,
  TGDLP==1 & is.na(DAYSTOTME) ~ TGDLP,
  TGDLP==1 & !is.na(DAYSTOTME) & DAYSTOTG>DAYSTOTME ~ as.integer(0),
  TGDLP==1 & !is.na(DAYSTOTME) & DAYSTOTG<=DAYSTOTME ~ TGDLP
))
# ANYDLP
comorb <- comorb %>% mutate(DAYSTOANYDLP = case_when(
  ANYDLP==0 ~ DAYSTOANYDLP,
  ANYDLP==1 & is.na(DAYSTOTME) ~ DAYSTOANYDLP,
  ANYDLP==1 & !is.na(DAYSTOTME) & DAYSTOANYDLP>DAYSTOTME ~ DAYSTOTME,
  ANYDLP==1 & !is.na(DAYSTOTME) & DAYSTOANYDLP<=DAYSTOTME ~ DAYSTOANYDLP
))
comorb <- comorb %>% mutate(ANYDLP = case_when(
  ANYDLP==0 ~ ANYDLP,
  ANYDLP==1 & is.na(DAYSTOTME) ~ ANYDLP,
  ANYDLP==1 & !is.na(DAYSTOTME) & DAYSTOANYDLP>DAYSTOTME ~ as.integer(0),
  ANYDLP==1 & !is.na(DAYSTOTME) & DAYSTOANYDLP<=DAYSTOTME ~ ANYDLP
))
# MIC
comorb <- comorb %>% mutate(DAYSTOMIC = case_when(
  MIC==0 ~ DAYSTOMIC,
  MIC==1 & is.na(DAYSTOTME) ~ DAYSTOMIC,
  MIC==1 & !is.na(DAYSTOTME) & DAYSTOMIC>DAYSTOTME ~ DAYSTOTME,
  MIC==1 & !is.na(DAYSTOTME) & DAYSTOMIC<=DAYSTOTME ~ DAYSTOMIC
))
comorb <- comorb %>% mutate(MIC = case_when(
  MIC==0 ~ MIC,
  MIC==1 & is.na(DAYSTOTME) ~ MIC,
  MIC==1 & !is.na(DAYSTOTME) & DAYSTOMIC>DAYSTOTME ~ as.integer(0),
  MIC==1 & !is.na(DAYSTOTME) & DAYSTOMIC<=DAYSTOTME ~ MIC
))
# MAC
comorb <- comorb %>% mutate(DAYSTOMAC = case_when(
  MAC==0 ~ DAYSTOMAC,
  MAC==1 & is.na(DAYSTOTME) ~ DAYSTOMAC,
  MAC==1 & !is.na(DAYSTOTME) & DAYSTOMAC>DAYSTOTME ~ DAYSTOTME,
  MAC==1 & !is.na(DAYSTOTME) & DAYSTOMAC<=DAYSTOTME ~ DAYSTOMAC
))
comorb <- comorb %>% mutate(MAC = case_when(
  MAC==0 ~ MAC,
  MAC==1 & is.na(DAYSTOTME) ~ MAC,
  MAC==1 & !is.na(DAYSTOTME) & DAYSTOMAC>DAYSTOTME ~ as.integer(0),
  MAC==1 & !is.na(DAYSTOTME) & DAYSTOMAC<=DAYSTOTME ~ MAC
))
# MIC.OR.MAC
comorb <- comorb %>% mutate(DAYSTOMIC.OR.MAC = case_when(
  MIC.OR.MAC==0 ~ DAYSTOMIC.OR.MAC,
  MIC.OR.MAC==1 & is.na(DAYSTOTME) ~ DAYSTOMIC.OR.MAC,
  MIC.OR.MAC==1 & !is.na(DAYSTOTME) & DAYSTOMIC.OR.MAC>DAYSTOTME ~ DAYSTOTME,
  MIC.OR.MAC==1 & !is.na(DAYSTOTME) & DAYSTOMIC.OR.MAC<=DAYSTOTME ~ DAYSTOMIC.OR.MAC
))
comorb <- comorb %>% mutate(MIC.OR.MAC = case_when(
  MIC.OR.MAC==0 ~ MIC.OR.MAC,
  MIC.OR.MAC==1 & is.na(DAYSTOTME) ~ MIC.OR.MAC,
  MIC.OR.MAC==1 & !is.na(DAYSTOTME) & DAYSTOMIC.OR.MAC>DAYSTOTME ~ as.integer(0),
  MIC.OR.MAC==1 & !is.na(DAYSTOTME) & DAYSTOMIC.OR.MAC<=DAYSTOTME ~ MIC.OR.MAC
))
# NEPHRO
comorb <- comorb %>% mutate(DAYSTONEPHRO = case_when(
  NEPHRO==0 ~ DAYSTONEPHRO,
  NEPHRO==1 & is.na(DAYSTOTME) ~ DAYSTONEPHRO,
  NEPHRO==1 & !is.na(DAYSTOTME) & DAYSTONEPHRO>DAYSTOTME ~ DAYSTOTME,
  NEPHRO==1 & !is.na(DAYSTOTME) & DAYSTONEPHRO<=DAYSTOTME ~ DAYSTONEPHRO
))
comorb <- comorb %>% mutate(NEPHRO = case_when(
  NEPHRO==0 ~ NEPHRO,
  NEPHRO==1 & is.na(DAYSTOTME) ~ NEPHRO,
  NEPHRO==1 & !is.na(DAYSTOTME) & DAYSTONEPHRO>DAYSTOTME ~ as.integer(0),
  NEPHRO==1 & !is.na(DAYSTOTME) & DAYSTONEPHRO<=DAYSTOTME ~ NEPHRO
))
# HYP
comorb <- comorb %>% mutate(DAYSTOHYP = case_when(
  HYP==0 ~ DAYSTOHYP,
  HYP==1 & is.na(DAYSTOTME) ~ DAYSTOHYP,
  HYP==1 & !is.na(DAYSTOTME) & DAYSTOHYP>DAYSTOTME ~ DAYSTOTME,
  HYP==1 & !is.na(DAYSTOTME) & DAYSTOHYP<=DAYSTOTME ~ DAYSTOHYP
))
comorb <- comorb %>% mutate(HYP = case_when(
  HYP==0 ~ HYP,
  HYP==1 & is.na(DAYSTOTME) ~ HYP,
  HYP==1 & !is.na(DAYSTOTME) & DAYSTOHYP>DAYSTOTME ~ as.integer(0),
  HYP==1 & !is.na(DAYSTOTME) & DAYSTOHYP<=DAYSTOTME ~ HYP
))
# FILAM
comorb <- comorb %>% mutate(DAYSTOFILAM = case_when(
  FILAM==0 ~ DAYSTOFILAM,
  FILAM==1 & is.na(DAYSTOTME) ~ DAYSTOFILAM,
  FILAM==1 & !is.na(DAYSTOTME) & DAYSTOFILAM>DAYSTOTME ~ DAYSTOTME,
  FILAM==1 & !is.na(DAYSTOTME) & DAYSTOFILAM<=DAYSTOTME ~ DAYSTOFILAM
))
comorb <- comorb %>% mutate(FILAM = case_when(
  FILAM==0 ~ FILAM,
  FILAM==1 & is.na(DAYSTOTME) ~ FILAM,
  FILAM==1 & !is.na(DAYSTOTME) & DAYSTOFILAM>DAYSTOTME ~ as.integer(0),
  FILAM==1 & !is.na(DAYSTOTME) & DAYSTOFILAM<=DAYSTOTME ~ FILAM
))
# NEURO
comorb <- comorb %>% mutate(DAYSTONEURO = case_when(
  NEURO==0 ~ DAYSTONEURO,
  NEURO==1 & is.na(DAYSTOTME) ~ DAYSTONEURO,
  NEURO==1 & !is.na(DAYSTOTME) & DAYSTONEURO>DAYSTOTME ~ DAYSTOTME,
  NEURO==1 & !is.na(DAYSTOTME) & DAYSTONEURO<=DAYSTOTME ~ DAYSTONEURO
))
comorb <- comorb %>% mutate(NEURO = case_when(
  NEURO==0 ~ NEURO,
  NEURO==1 & is.na(DAYSTOTME) ~ NEURO,
  NEURO==1 & !is.na(DAYSTOTME) & DAYSTONEURO>DAYSTOTME ~ as.integer(0),
  NEURO==1 & !is.na(DAYSTOTME) & DAYSTONEURO<=DAYSTOTME ~ NEURO
))
# RETINO
comorb <- comorb %>% mutate(DAYSTORETINO = case_when(
  RETINO==0 ~ DAYSTORETINO,
  RETINO==1 & is.na(DAYSTOTME) ~ DAYSTORETINO,
  RETINO==1 & !is.na(DAYSTOTME) & DAYSTORETINO>DAYSTOTME ~ DAYSTOTME,
  RETINO==1 & !is.na(DAYSTOTME) & DAYSTORETINO<=DAYSTOTME ~ DAYSTORETINO
))
comorb <- comorb %>% mutate(RETINO = case_when(
  RETINO==0 ~ RETINO,
  RETINO==1 & is.na(DAYSTOTME) ~ RETINO,
  RETINO==1 & !is.na(DAYSTOTME) & DAYSTORETINO>DAYSTOTME ~ as.integer(0),
  RETINO==1 & !is.na(DAYSTOTME) & DAYSTORETINO<=DAYSTOTME ~ RETINO
))
# MVD
comorb <- comorb %>% mutate(DAYSTOMVD = case_when(
  MVD==0 ~ DAYSTOMVD,
  MVD==1 & is.na(DAYSTOTME) ~ DAYSTOMVD,
  MVD==1 & !is.na(DAYSTOTME) & DAYSTOMVD>DAYSTOTME ~ DAYSTOTME,
  MVD==1 & !is.na(DAYSTOTME) & DAYSTOMVD<=DAYSTOTME ~ DAYSTOMVD
))
comorb <- comorb %>% mutate(MVD = case_when(
  MVD==0 ~ MVD,
  MVD==1 & is.na(DAYSTOTME) ~ MVD,
  MVD==1 & !is.na(DAYSTOTME) & DAYSTOMVD>DAYSTOTME ~ as.integer(0),
  MVD==1 & !is.na(DAYSTOTME) & DAYSTOMVD<=DAYSTOTME ~ MVD
))
# GLYC
comorb <- comorb %>% mutate(DAYSTOGLYC = case_when(
  GLYC==0 ~ DAYSTOGLYC,
  GLYC==1 & is.na(DAYSTOTME) ~ DAYSTOGLYC,
  GLYC==1 & !is.na(DAYSTOTME) & DAYSTOGLYC>DAYSTOTME ~ DAYSTOTME,
  GLYC==1 & !is.na(DAYSTOTME) & DAYSTOGLYC<=DAYSTOTME ~ DAYSTOGLYC
))
comorb <- comorb %>% mutate(GLYC = case_when(
  GLYC==0 ~ GLYC,
  GLYC==1 & is.na(DAYSTOTME) ~ GLYC,
  GLYC==1 & !is.na(DAYSTOTME) & DAYSTOGLYC>DAYSTOTME ~ as.integer(0),
  GLYC==1 & !is.na(DAYSTOTME) & DAYSTOGLYC<=DAYSTOTME ~ GLYC
))

# Drop bariatric surgery variables
comorb <- comorb %>% select(-c(TMETYPE, DAYSTOTME))

# Merge AME (adjudicated medical events) into comorb as wide per-person summary
# Keep first event of each type per person
ame_wide <- AME %>%
  filter(AMENAME == 1)  %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOARRHYTHMIA=DAYSTOAME) %>% mutate(ARRHYTHMIA=1) %>%
  full_join(AME %>% filter(AMENAME == 2)  %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOCAD=DAYSTOAME)       %>% mutate(CAD=1),       by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 3)  %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOCHF=DAYSTOAME)        %>% mutate(CHF=1),       by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 4)  %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOLVSD=DAYSTOAME)       %>% mutate(LVSD=1),      by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 5)  %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOMI=DAYSTOAME)         %>% mutate(MI=1),        by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 6)  %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOPAD=DAYSTOAME)        %>% mutate(PAD=1),       by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 8)  %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTODVT=DAYSTOAME)        %>% mutate(DVT=1),       by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 9)  %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOSTROKE=DAYSTOAME)     %>% mutate(STROKE=1),    by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 11) %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOTIA=DAYSTOAME)        %>% mutate(TIA=1),       by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 17) %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOCKD=DAYSTOAME)        %>% mutate(CKD=1),       by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 18) %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTOESKD=DAYSTOAME)       %>% mutate(ESKD=1),      by="releaseid") %>%
  full_join(AME %>% filter(AMENAME == 26) %>% arrange(releaseid, DAYSTOAME) %>% group_by(releaseid) %>% filter(row_number()==1) %>% ungroup() %>% select(releaseid, DAYSTOAME) %>% rename(DAYSTODEATH=DAYSTOAME)      %>% mutate(DEATH=1),     by="releaseid") %>%
  full_join(tme %>% filter(TMETYPE == 5)  %>% arrange(RELEASEID, DAYSTOTME) %>% group_by(RELEASEID) %>% filter(row_number()==1) %>% ungroup() %>% rename(releaseid=RELEASEID)  %>% select(releaseid, DAYSTOTME) %>% rename(DAYSTOULCER=DAYSTOTME) %>% mutate(ULCER=1), by="releaseid")

# Set event indicators to 0 for those without the event
ame_wide <- ame_wide %>%
  mutate(across(c(ARRHYTHMIA, CAD, CHF, LVSD, MI, PAD, DVT, STROKE, TIA, CKD, ESKD, DEATH, ULCER),
                ~ifelse(is.na(.), 0, .)))

# Merge into comorb
comorb <- merge(comorb, ame_wide, by="releaseid", all.x=TRUE, all.y=FALSE)

# Save
save(comorb, file="./Clinical data/comorb_testing_CKM.Rdata")

###############################################################
# Prepare shared reference datasets used across sections      #
###############################################################

# PAT - demographics
PAT$racedesc[PAT$race==1] <- "Non-Hispanic Black"
PAT$racedesc[PAT$race==2] <- "Hispanic"
PAT$racedesc[PAT$race==3] <- "Non-Hispanic White"
PAT$racedesc[PAT$race==4] <- "Other"
PAT$houseincdesc[PAT$houseinc==1] <- "<$24,999"
PAT$houseincdesc[PAT$houseinc==2] <- "$25,000-$49,999"
PAT$houseincdesc[PAT$houseinc==3] <- ">$50,000"
keepPAT <- PAT %>% select(releaseid, age, sex, dxtime, race, houseinc, racedesc, houseincdesc)
keepPAT$sex_char <- ifelse(keepPAT$sex == 1, "F", "M")

# AGEBASE - uncollapsed age at baseline
AGEBASE$releaseid <- AGEBASE$RELEASEID
AGEBASE$RELEASEID <- NULL

# PRIMOUT - treatment group assignment
PRIMOUT$txdesc[PRIMOUT$tx==1] <- "Metformin only"
PRIMOUT$txdesc[PRIMOUT$tx==2] <- "Metformin + rosiglitazone"
PRIMOUT$txdesc[PRIMOUT$tx==3] <- "Metformin + lifestyle"
keepPRIMOUT <- PRIMOUT %>% select(releaseid, tx, txdesc)

# BIRTH WEIGHT
BW[, 2:6] <- apply(BW[, 2:6], 2, as.numeric)

# DXA
DXA$visit <- DXA$mvisit
keepDXA <- DXA %>% select(releaseid, visit, WB_TOT_PFAT, WB_TOT_PFAT_P)
baseDXA <- keepDXA %>% filter(visit=="M00")

# LIPOSCIENCE
colnames(lipo)[5:ncol(lipo)] <- paste0("LIPO_", colnames(lipo)[5:ncol(lipo)])
lipo <- lipo %>% mutate(visit = case_when(
  LPMONTH == 0  ~ "M00",
  LPMONTH == 12 ~ "M12",
  LPMONTH == 24 ~ "M24",
  LPMONTH == 36 ~ "M36",
  LPMONTH == 48 ~ "M48",
  LPMONTH == 60 ~ "M60"
))
lipo$releaseid <- lipo$RELEASEID
lipo$RELEASEID <- NULL
base_lipo <- lipo %>% filter(visit == "M00")

###############################################################
# Create baseline dataset (baserisk)                          #
###############################################################

basecbl    <- CBL %>% filter(mvisit=="M00")
baseaddcbl <- ADDCBL %>% filter(mvisit=="M00")

baserisk <- merge(BASELINE, basecbl, by="releaseid", all.x=T, all.y=T)
baserisk <- merge(baserisk, baseaddcbl, by="releaseid", all.x=T, all.y=T)
baserisk$si_1_ins0 <- 1/baserisk$ins0min
baserisk$log_trig <- log(baserisk$Trig)
baserisk <- baserisk %>% select(releaseid, HbA1c, log_trig, sbp, dbp, uacid, si_1_ins0, UAlbCreat, bmi, bmipct, HDL, codi,
                                EstCreatClear, SerumCreat, serumcystc, glu0min, ins0min, ALT, AST, wastcirc, height)
baserisk$map <- baserisk$dbp + ((1/3)*(baserisk$sbp - baserisk$dbp))
baserisk <- merge(baserisk, keepPAT, by="releaseid", all.x=T, all.y=F)
baserisk$age <- NULL
baserisk <- merge(baserisk, AGEBASE, by="releaseid", all.x=T, all.y=F)
baserisk <- merge(baserisk, keepPRIMOUT, by="releaseid", all.x=T, all.y=F)
baserisk <- merge(baserisk, BW, by="releaseid", all.x=T, all.y=T)
baserisk <- merge(baserisk, baseDXA, by="releaseid", all.x=T, all.y=T)
baserisk <- baserisk %>% mutate(relative_fat_mass = case_when(
  sex_char == "F" ~  76 - ((20*height/wastcirc)),
  sex_char == "M" ~ 64 - ((20*height)/wastcirc)),
  TRUE = NA_real_
)
baserisk <- full_join(baserisk, base_lipo, by="releaseid")
baserisk$homa_ir    <- (baserisk$ins0min * baserisk$glu0min) / 405
baserisk$waist_height <- baserisk$wastcirc / baserisk$height

# Save
save(baserisk, file="./Clinical data/TODAY/baserisk_testing_CKM.Rdata")

###############################################################
# Create dataset of risk factors at 10-year visit (TODAY2)    #
###############################################################

CBL_TODAY2$visit <- CBL_TODAY2$pvisit

# 10-year visit labs
CBL_TODAY2_KEEP <- CBL_TODAY2 %>% filter(pvisit=="P120") %>%
  select(releaseid, hba1c, trig, estcreatclear, ualbcreat)

# Insulin was measured at 9-year visit, not 10-year
INS_TODAY2_KEEP <- CBL_TODAY2 %>% filter(pvisit=="P108") %>%
  select(releaseid, ins)

# BP and BMI at 10-year visit
VISIT_TODAY2$releaseid <- VISIT_TODAY2$RELEASEID
VISIT_TODAY2$sbp <- VISIT_TODAY2$SBP
VISIT_TODAY2$bmi <- VISIT_TODAY2$BMI
VISIT_TODAY2_KEEP_YR10 <- VISIT_TODAY2 %>% filter(PVISIT=="P120") %>%
  select(releaseid, sbp, bmi)

yr10risk <- merge(CBL_TODAY2_KEEP, INS_TODAY2_KEEP, by="releaseid", all.x=T, all.y=T)
yr10risk <- merge(yr10risk, VISIT_TODAY2_KEEP_YR10, by="releaseid", all.x=T, all.y=T)
yr10risk$si_1_ins0 <- 1/yr10risk$ins
yr10risk$log_trig  <- log(yr10risk$trig)

# Save
save(yr10risk, file="./Clinical data/TODAY/yr10risk_testing_CKM.Rdata")

###############################################################
# Create longitudinal dataset (long)                          #
###############################################################

# TODAY BASELINE - baseline BMI
BASELINE$visit      <- "M00"
BASELINE$visit_days <- BASELINE$days
BASELINE_keep <- BASELINE %>% select(releaseid, visit, bmi, visit_days, bmipct, wastcirc, height, weight)
BASELINE_keep$bsa_dubois <- 0.007184 * (BASELINE_keep$height^0.725) * (BASELINE_keep$weight^0.425)

# TODAY VISIT - BMI
VISIT$visit      <- VISIT$mvisit
VISIT$visit_days <- VISIT$days
VISIT_keep <- VISIT %>% select(releaseid, visit, bmi, height, weight, visit_days, bmipct, wastcirc)
VISIT_keep$bsa_dubois <- 0.007184 * (VISIT_keep$height^0.725) * (VISIT_keep$weight^0.425)

# TODAY CBL - eIS, labs, and inflammatory/vascular markers
CBL_keep  <- CBL %>%
  filter(!mvisit=="R") %>%
  rename(visit=mvisit, trig=Trig, hba1c=HbA1c, alt=ALT, ast=AST) %>%
  select(any_of(c(
    "releaseid", "visit", "ins0min", "trig", "hba1c", "alt", "ast", "glu0min",
    "UAlbCreat", "UAlb", "UCreat", "SerumCreat", "HDL", "LDL", "Chol",
    "LDLB", "LDLC", "LDLCB", "VLDL", "ApoB", "FFA", "FIB", "HOM", "PIN",
    "Rf", "VB12", "PAI1", "IL6", "il1", "hsCRP", "mcp1", "tnfa", "tnfr1",
    "tnfr2", "icam1", "vcam1", "eselectin", "vegf", "fgf23"
  )))

# TODAY ADDCBL - coDI and additional biomarkers
ADDCBL_keep  <- ADDCBL %>%
  rename(visit=mvisit) %>%
  select(any_of(c(
    "releaseid", "visit", "days", "codi", "bnp", "troponin", "adiponectin",
    "hmwa", "edol", "shbg", "testosterone", "ckd_gfr", "eGFR_FAS", "insinv"
  )))

# TODAY2 VISIT - BMI and waist circ
VISIT_TODAY2$bmi        <- VISIT_TODAY2$BMI
VISIT_TODAY2$visit      <- VISIT_TODAY2$PVISIT
VISIT_TODAY2$visit_days <- VISIT_TODAY2$DAYS
VISIT_TODAY2$bmipct     <- VISIT_TODAY2$BMIPCT
VISIT_TODAY2_KEEP <- VISIT_TODAY2 %>%
  select(releaseid, visit, bmi, HEIGHT, WEIGHT, visit_days, bmipct) %>%
  rename(height=HEIGHT, weight=WEIGHT)
VISIT_TODAY2_KEEP$wastcirc   <- NA
VISIT_TODAY2_KEEP$bsa_dubois <- 0.007184 * (VISIT_TODAY2_KEEP$height^0.725) * (VISIT_TODAY2_KEEP$weight^0.425)

# TODAY2 CBL - eIS, coDI, and additional markers (rename to match TODAY column names)
# Note: drop visit and ins0min if already present from yr10risk section above
CBL_TODAY2_KEEP <- CBL_TODAY2 %>%
  select(-any_of(c("visit", "ins0min"))) %>%
  rename(any_of(c(
    visit="pvisit", ins0min="ins",
    hba1c="hba1c", trig="trig", alt="alt", ast="ast",
    LDL="ldl", Chol="chol", hsCRP="hscrp", FFA="ffa", FIB="fib", HOM="hom",
    LDLB="ldlb", LDLC="ldlc", PIN="pin", Rf="rf", ApoB="apob",
    EstCreatClear="estcreatclear", HDL="hdl", LDLCB="ldlcb",
    SerumCreat="serumcreat", UAlb="ualb", UAlbCreat="ualbcreat",
    UCreat="ucreat", VB12="vb12", VLDL="vldl", IL6="il6", PAI1="pai1"
  ))) %>%
  select(any_of(c(
    "releaseid", "visit", "ins0min", "codi", "trig", "hba1c", "alt", "ast", "glu0min",
    "UAlbCreat", "UAlb", "UCreat", "SerumCreat", "HDL", "LDL", "Chol",
    "LDLB", "LDLC", "LDLCB", "VLDL", "ApoB", "FFA", "FIB", "HOM", "PIN",
    "Rf", "VB12", "PAI1", "IL6", "hsCRP"
  )))

# TODAY2 ADDCBL - same dataset as TODAY ADDCBL, already captured in ADDCBL_keep above

# Calculate length of follow-up for each person
ALL_VISITS_TODAY  <- VISIT %>% select(releaseid, visit)
ALL_VISITS_TODAY2 <- VISIT_TODAY2 %>% select(releaseid, visit)
ALL_VISITS <- rbind(ALL_VISITS_TODAY, ALL_VISITS_TODAY2)
agebase    <- baserisk %>% select(releaseid, AGEBASE)
ALL_VISITS <- merge(ALL_VISITS, agebase, by="releaseid", all.x=T, all.y=T)
ALL_VISITS$visitnum <- as.numeric(stringr::str_remove(stringr::str_remove(ALL_VISITS$visit, "P"), "M"))
ALL_VISITS <- arrange(ALL_VISITS, releaseid, desc(visitnum))
ALL_VISITS <- ALL_VISITS %>% group_by(releaseid) %>% filter(row_number()==1)
ALL_VISITS$fup_years    <- ALL_VISITS$visitnum / 12
ALL_VISITS <- ALL_VISITS %>% mutate(age_last_visit = sum(AGEBASE, fup_years))
ALL_VISITS$ge25yrs <- ifelse(ALL_VISITS$age_last_visit >= 25, 1, 0)
fup_length <- ALL_VISITS %>% select(releaseid, fup_years, age_last_visit, ge25yrs)

# Merge all longitudinal data
long <- bind_rows(BASELINE_keep, VISIT_keep, VISIT_TODAY2_KEEP)
labs <- merge(CBL_keep, ADDCBL_keep, by=c("releaseid","visit"), all.x=T, all.y=T)
labs <- bind_rows(labs, CBL_TODAY2_KEEP)
long <- merge(long, labs, by=c("releaseid","visit"), all.x=T, all.y=T)
long <- merge(long, keepDXA, by=c("releaseid","visit"), all.x=T, all.y=T)

# PedsQL - standardize visit column and column names, then combine and merge into long
pedsql_cols <- c("releaseid", "visit", "G01WALK", "G02RUN", "G03SPORT",
                 "G04LIFT", "G05BATH", "G06CHORE", "G07HURT", "G08ENERG")

# TODAY: uppercase column names, rename mvisit -> visit
names(PEDSQLGC_TODAY) <- toupper(names(PEDSQLGC_TODAY))
names(PEDSQLGT_TODAY) <- toupper(names(PEDSQLGT_TODAY))
PEDSQLGC_TODAY <- PEDSQLGC_TODAY %>% rename(releaseid=RELEASEID, visit=MVISIT)
PEDSQLGT_TODAY <- PEDSQLGT_TODAY %>% rename(releaseid=RELEASEID, visit=MVISIT)

# TODAY2: rename pvisit -> visit, releaseid standardization
PEDSQLGA_TODAY2 <- PEDSQLGA_TODAY2 %>% rename(releaseid=RELEASEID, visit=PVISIT)
PEDSQLGY_TODAY2 <- PEDSQLGY_TODAY2 %>% rename(releaseid=RELEASEID, visit=PVISIT)
PEDSQLGT_TODAY2 <- PEDSQLGT_TODAY2 %>% rename(releaseid=RELEASEID, visit=PVISIT)

pedsql_long <- bind_rows(
  PEDSQLGC_TODAY  %>% select(any_of(pedsql_cols)),
  PEDSQLGT_TODAY  %>% select(any_of(pedsql_cols)),
  PEDSQLGA_TODAY2 %>% select(any_of(pedsql_cols)),
  PEDSQLGY_TODAY2 %>% select(any_of(pedsql_cols)),
  PEDSQLGT_TODAY2 %>% select(any_of(pedsql_cols))
) %>% distinct(releaseid, visit, .keep_all=TRUE)

long <- merge(long, pedsql_long, by=c("releaseid","visit"), all.x=T, all.y=T)

# Add numeric visit variable
long$visit_num <- as.numeric(str_sub(long$visit, 2, length(long$visit)))

# Calculated variables
long$si_1_ins0  <- 1/long$ins0min
long <- merge(long, fup_length, by="releaseid", all.x=T, all.y=T)
long <- full_join(long, lipo, by=c("releaseid","visit"))
long$homa_ir     <- (long$ins0min * long$glu0min) / 405
long$waist_height <- long$wastcirc / long$height

# Diagnostics
long_unique <- long %>% select(releaseid, fup_years) %>% unique()
summary(long_unique$fup_years)

# Time to GLYC in those with the outcome
summary(comorb[comorb$GLYC==1,]$DAYSTOGLYC)

# Save
save(long, file="./Clinical data/TODAY/clinical_data_long_testing_CKM.Rdata")
