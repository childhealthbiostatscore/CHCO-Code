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

# COMORB dataset
comorb <- read.csv("./Clinical data/COMORB.csv")
comorb$MIC.OR.MAC <- ifelse(comorb$MAC==1 | comorb$MIC==1,1,
                            ifelse(is.na(comorb$MAC) & is.na(comorb$MIC), NA, 0))
comorb <- comorb %>% mutate(DAYSTOMIC.OR.MAC = pmin(DAYSTOMIC, DAYSTOMAC))
comorb$releaseid <- comorb$RELEASEID
comorb$RELEASEID <- NULL
# Merge in bariatric surgery info
# read in TME dataset to find those who had bariatric surgery
tme <- read.csv("./Clinical data/TODAY2/TME.csv")
tme <- tme %>% filter(TMETYPE==1)
tme <- tme %>% select(RELEASEID,TMETYPE,DAYSTOTME)
colnames(tme) <- c("releaseid","TMETYPE","DAYSTOTME")
comorb <- merge(comorb, tme, by="releaseid", all.x = T, all.y=F)
# now censor all outcomes at time of TME for those who had bariatric surgery
# need to reset the indicator variable as well as the time variable
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
# drop bariatric surgery variables
comorb <- comorb %>% select(-c(TMETYPE,DAYSTOTME))
# Save
save(comorb,file = "./Clinical data/comorb.Rdata")

# CBL
CBL <- read.csv("./Clinical data/TODAY/CBL.csv")

# ADDCBL
ADDCBL <- read.csv("./Clinical data/TODAY/ADDCBL.csv")

# BASELINE
BASELINE <- read.csv("./Clinical data/TODAY/BASELINE.csv")

# PAT
PAT <- read.csv("./Clinical data/TODAY/PAT.csv")
PAT$racedesc[PAT$race==1]<- "Non-Hispanic Black"
PAT$racedesc[PAT$race==2]<- "Hispanic"
PAT$racedesc[PAT$race==3]<- "Non-Hispanic White"
PAT$racedesc[PAT$race==4]<- "Other"
PAT$houseincdesc[PAT$houseinc==1]<- "<$24,999"
PAT$houseincdesc[PAT$houseinc==2]<- "$25,000-$49,999"
PAT$houseincdesc[PAT$houseinc==3]<- ">$50,000"
keepPAT <- PAT %>% select(releaseid,age,sex,dxtime,race,houseinc,racedesc,houseincdesc)
keepPAT$sex_char <- ifelse(keepPAT$sex == 1, "F", "M")

# AGEBASE - uncollapsed age at baseline
AGEBASE <- read.csv("./Clinical data/TODAY/AGEBASE.csv")
AGEBASE$releaseid <- AGEBASE$RELEASEID
AGEBASE$RELEASEID <- NULL

# PRIMOUT - treatment group assignment
PRIMOUT <- read.csv("./Clinical data/TODAY/PRIMOUT.csv")
PRIMOUT$txdesc[PRIMOUT$tx==1]<- "Metformin only"
PRIMOUT$txdesc[PRIMOUT$tx==2]<- "Metformin + rosiglitazone"
PRIMOUT$txdesc[PRIMOUT$tx==3]<- "Metformin + lifestyle"
keepPRIMOUT <- PRIMOUT %>% select(releaseid,tx,txdesc)

# BIRTH WEIGHT 
BW <- read.csv("./Clinical data/TODAY/Birthweight.csv")
BW[,2:6] <- apply(BW[,2:6], 2, as.numeric)

# DXA
DXA <- read.csv("./Clinical data/TODAY/DEXA.csv")
DXA$visit <- DXA$mvisit
keepDXA <- DXA %>% select(releaseid, visit, WB_TOT_PFAT, WB_TOT_PFAT_P)
baseDXA <- keepDXA %>% filter(visit=="M00")

# create new dataset of baseline risk factors
basecbl <- CBL %>% filter(mvisit=="M00")
baseaddcbl <- ADDCBL %>% filter(mvisit=="M00")
baserisk <- merge(BASELINE, basecbl, by="releaseid", all.x=T, ally=T)
baserisk <- merge(baserisk, baseaddcbl, by="releaseid", all.x=T, ally=T)
baserisk$si_1_ins0 <- 1/baserisk$ins0min
baserisk$log_trig <- log(baserisk$Trig)
baserisk <- baserisk %>% select(releaseid, HbA1c, log_trig, sbp, dbp, uacid, si_1_ins0, UAlbCreat, bmi, bmipct, HDL, codi,
                                EstCreatClear,SerumCreat,serumcystc,glu0min,ins0min,ALT,AST,wastcirc,height)
baserisk$map <- baserisk$dbp + ((1/3)*(baserisk$sbp - baserisk$dbp))
baserisk <- merge(baserisk,keepPAT,by="releaseid",all.x = T,all.y = F)
baserisk$age <- NULL
baserisk <- merge(baserisk,AGEBASE,by="releaseid",all.x = T,all.y = F)
baserisk <- merge(baserisk, keepPRIMOUT,by="releaseid",all.x = T,all.y = F)
baserisk <- merge(baserisk, BW, by="releaseid", all.x = T, all.y = T)
baserisk <- merge(baserisk, baseDXA, by="releaseid", all.x=T, all.y=T)
baserisk <- baserisk %>% mutate(relative_fat_mass = case_when(
  sex_char == "F" ~  76 - ((20*height/wastcirc)),
  sex_char == "M" ~ 64 - ((20*height)/wastcirc)),
  TRUE == NA_real_
)
  
# Save
save(baserisk,file = "./Clinical data/TODAY/baserisk.Rdata")

###############################################################
# create dataset of risk factors at 10 year visit from TODAY2 #
###############################################################

# CBL
CBL_TODAY2 <- read.csv("./Clinical data/TODAY2/CBL.csv")
CBL_TODAY2_KEEP <- CBL_TODAY2 %>% filter(pvisit=="P120") %>% select(releaseid, hba1c, trig, 
                                                                    estcreatclear, ualbcreat)
# INS 
# insulin was measured at 9 year visit not 10 year
INS_TODAY2 <- read.csv("./Clinical data/TODAY2/CBL.csv")
INS_TODAY2_KEEP <- INS_TODAY2 %>% filter(pvisit=="P108") %>% select(releaseid, ins)

# ADDCBL
ADDCBL_TODAY2 <- read.csv("./Clinical data/TODAY2/ADDCBL.csv")

# VISIT
VISIT_TODAY2 <- read.csv("./Clinical data/TODAY2/VISIT.csv")
VISIT_TODAY2$releaseid <- VISIT_TODAY2$RELEASEID
VISIT_TODAY2$sbp <- VISIT_TODAY2$SBP
VISIT_TODAY2$bmi <- VISIT_TODAY2$BMI
VISIT_TODAY2_KEEP <- VISIT_TODAY2 %>% filter(PVISIT=="P120") %>% select(releaseid, sbp, bmi)

yr10risk <- merge(CBL_TODAY2_KEEP, INS_TODAY2_KEEP, by="releaseid", all.x = T, all.y = T)
yr10risk <- merge(yr10risk, VISIT_TODAY2_KEEP, by="releaseid", all.x = T, all.y = T)
yr10risk$si_1_ins0 <- 1/yr10risk$ins
yr10risk$log_trig <- log(yr10risk$trig)

# Save
save(yr10risk,file = "./Clinical data/TODAY/yr10risk.Rdata")

###############################################################
# create longitudinal datasets from TODAY and TODAY2          #
###############################################################

# TODAY BASELINE - baseline BMI 
BASELINE <- read.csv("./Clinical data/TODAY/BASELINE.csv")
BASELINE$visit <- "M00"
BASELINE$visit_days <- BASELINE$days
BASELINE_keep <- BASELINE %>% select(releaseid,visit,bmi,visit_days,bmipct,wastcirc,height,weight)
#BASELINE_keep$height <- NA
#BASELINE_keep$weight <- NA

# TODAY VISIT - BMI 
VISIT <- read.csv("./Clinical data/TODAY/VISIT.csv")
VISIT$visit <- VISIT$mvisit
VISIT$visit_days <- VISIT$days
VISIT_keep <- VISIT %>% select(releaseid,visit,bmi,height,weight,visit_days,bmipct,wastcirc)

# TODAY CBL - eIS
CBL <- read.csv("./Clinical data/TODAY/CBL.csv")
CBL$visit <- CBL$mvisit
CBL_keep <- CBL %>% select(releaseid,visit,ins0min,Trig,HbA1c,ALT,AST)
CBL_keep <- CBL_keep %>% filter(!visit=="R")
CBL_keep$trig <- CBL_keep$Trig
CBL_keep$Trig <- NULL
CBL_keep$hba1c <- CBL_keep$HbA1c
CBL_keep$HbA1c <- NULL
CBL_keep$alt <- CBL_keep$ALT
CBL_keep$ALT <- NULL
CBL_keep$ast <- CBL_keep$AST
CBL_keep$AST <- NULL

# TODAY ADDCBL - coDI
ADDCBL <- read.csv("./Clinical data/TODAY/ADDCBL.csv")
ADDCBL$visit <- ADDCBL$mvisit
ADDCBL_keep <- ADDCBL %>% select(releaseid,visit,codi)

# TODAY2 VISIT - BMI and waist circ
VISIT_TODAY2 <- read.csv("./Clinical data/TODAY2/VISIT.csv")
VISIT_TODAY2$releaseid <- VISIT_TODAY2$RELEASEID
VISIT_TODAY2$bmi <- VISIT_TODAY2$BMI
VISIT_TODAY2$visit <- VISIT_TODAY2$PVISIT
VISIT_TODAY2$visit_days <- VISIT_TODAY2$DAYS
VISIT_TODAY2$bmipct <- VISIT_TODAY2$BMIPCT
VISIT_TODAY2_KEEP <- VISIT_TODAY2 %>% select(releaseid, visit, bmi, HEIGHT, WEIGHT, visit_days, bmipct)
VISIT_TODAY2_KEEP$height <- VISIT_TODAY2_KEEP$HEIGHT
VISIT_TODAY2_KEEP$weight <- VISIT_TODAY2_KEEP$WEIGHT
VISIT_TODAY2_KEEP$HEIGHT <- NULL
VISIT_TODAY2_KEEP$WEIGHT <- NULL
VISIT_TODAY2_KEEP$wastcirc <- NA

# TODAY2 CBL - eIS and coDI
CBL_TODAY2 <- read.csv("./Clinical data/TODAY2/CBL.csv")
CBL_TODAY2$visit <- CBL_TODAY2$pvisit
CBL_TODAY2$ins0min <- CBL_TODAY2$ins
CBL_TODAY2_KEEP <- CBL_TODAY2 %>% select(releaseid, visit, ins0min, codi,  trig, hba1c, alt, ast)

# TODAY2 DXA - no DXA in TODAY2

# calculate length of follow-up for each person
# we don't have visit dates, just visit numbers
ALL_VISITS_TODAY <- VISIT %>% select(releaseid, visit)
ALL_VISITS_TODAY2 <- VISIT_TODAY2 %>% select(releaseid, visit)
ALL_VISITS <- rbind(ALL_VISITS_TODAY, ALL_VISITS_TODAY2)
agebase <- baserisk %>% select(releaseid, AGEBASE)
ALL_VISITS <- merge(ALL_VISITS, agebase, by = "releaseid", all.x = T, all.y = T)
ALL_VISITS$visitnum <- stringr::str_remove(ALL_VISITS$visit, "P")
ALL_VISITS$visitnum <- stringr::str_remove(ALL_VISITS$visitnum, "M")
ALL_VISITS$visitnum <- as.numeric(ALL_VISITS$visitnum)
ALL_VISITS <- arrange(ALL_VISITS, releaseid, desc(visitnum))
ALL_VISITS <- ALL_VISITS %>% group_by(releaseid) %>% filter(row_number()==1)
ALL_VISITS$fup_years <- ALL_VISITS$visitnum/12
ALL_VISITS <- ALL_VISITS %>% mutate(age_last_visit = sum(AGEBASE, fup_years))
ALL_VISITS$ge25yrs <- ifelse(ALL_VISITS$age_last_visit >= 25, 1, 0)
fup_length <- ALL_VISITS %>% select(releaseid, fup_years, age_last_visit, ge25yrs)

# merge all data
long <- rbind(BASELINE_keep,VISIT_keep,VISIT_TODAY2_KEEP)
labs <- merge(CBL_keep,ADDCBL_keep,by=c("releaseid","visit"),all.x = T, all.y = T)
labs <- rbind(labs,CBL_TODAY2_KEEP)
long <- merge(long,labs,by=c("releaseid","visit"),all.x = T, all.y = T)
long <- merge(long,keepDXA,by=c("releaseid","visit"),all.x = T, all.y = T)
# add a numeric visit variable
long$visit_num <- as.numeric(str_sub(long$visit,2,length(long$visit)))
# calculated variables - eIS
long$si_1_ins0 <- 1/long$ins0min
long <- merge(long, fup_length, by = "releaseid", all.x = T, all.y = T)

# length of follow-up
long_unique <- long %>% select(releaseid, fup_years) %>% unique() 
summary(long_unique$fup_years)

# time to LOGC in those with the outcome
summary(comorb[comorb$GLYC==1,]$DAYSTOGLYC)

# write file
save(long,file = "./Clinical data/TODAY/clinical_data_long.Rdata")
