---
title: "CKM in TODAY"
author: "Laura Pyle / Jairo A Pinzon"
format: html
editor: visual
echo: FALSE
warning: FALSE
toc: TRUE
---

```{r}
library(dplyr)
library(arsenal)
library(tidyr)

getwd()

comorb <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/COMORB.csv")
load("/Users/jpcortes/Documents/Temp_files_TODAY/clinical_data_long.Rdata")

#ame adjudicated medical events
ame <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/AME.csv")

# TODAY ECHO
today_echo <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/ECHO_today.csv")
today_echo$TIMEPOINT <- "1"
# TODAY2 ECHO
today2_echo <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/ECHO_today2.csv")
today2_echo$TIMEPOINT <- "2"

## we need to decide which timepoint to use, for certain labs, or if AME would be enough so we don't redenie criteria for dx or certain conditions. 
cbl_1 <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/CBL_TODAY1.csv")
cbl_1$TIMEPOINT <- "1"

cbl_2 <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/CBL_TODAY2.csv")
cbl_2$TIMEPOINT <- "2"

addcbl_1 <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/ADDCBL_TODAY1.csv")
addcbl_1$TIMEPOINT <- "1"

addcbl_2 <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/ADDCBL_TODAY2.csv")
addcbl_2$TIMEPOINT <- "2"

# don't have GLS in TODAY, just TODAY2 - but there are two time points
speckle <- read.csv('/Users/jpcortes/Documents/Temp_files_TODAY/SPECKLE.csv')

# TODAY2 PWV
pwv <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/PWV.csv")

# baseline integration to estimate bsa, but also to keep abdominal circumference and other relevant variables. 

baseline <- read.csv(("/Users/jpcortes/Documents/Temp_files_TODAY/BASELINE.csv"))

## quality of life PEDSQL questionaire about health and activities (exercise tolerance and activities) for the 3 groups adults, young adults and teens ? from 0 to 4 in severity scale. 

pedsQL_A <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/PEDSQLGA.csv")
pedsQL_A <- read.csv("/Users/jarpincortes/Library/CloudStorage/GoogleDrive-jairo.pinzoncortes@monash.edu/My Drive/Post-Doctoral positions and Meetings/Temp_files_TODAY/PEDSQLGA.csv")
pedsQL_Y <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/PEDSQLGY.csv")
pedsQL_Y <- read.csv("/Users/jarpincortes/Library/CloudStorage/GoogleDrive-jairo.pinzoncortes@monash.edu/My Drive/Post-Doctoral positions and Meetings/Temp_files_TODAY/PEDSQLGY.csv")
pedsQL_T <- read.csv("/Users/jpcortes/Documents/Temp_files_TODAY/PEDSQLGT.csv")
pedsQL_T <- read.csv("/Users/jarpincortes/Library/CloudStorage/GoogleDrive-jairo.pinzoncortes@monash.edu/My Drive/Post-Doctoral positions and Meetings/Temp_files_TODAY/PEDSQLGT.csv")

```

```{r}
# read in COMORB adjudicated endpoints
comorb <- read.csv("/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/COMORB.csv")
# get fup time from visit dates - the problem is that these are not granular enough 
load("/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/TODAY/clinical_data_long.Rdata")

# in comorb, figure out length of follow-up by choosing number of days for an event the participant did not have
comorb$fup_time <- ifelse(comorb$HTN == 0, comorb$DAYSTOHTN, 
                          ifelse(comorb$LDLDLP == 0, comorb$DAYSTOLDL,
                                 ifelse(comorb$NEURO == 0, comorb$DAYSTONEURO, 
                                        ifelse(comorb$DNE == 0, comorb$DAYSTODNE, 
                                               ifelse(comorb$FILAM == 0, comorb$DAYSTOFILAM, 
                                                     ifelse(comorb$RETINO == 0, comorb$DAYSTORETINO, 
                                                            ifelse(comorb$TGDLP == 0, comorb$DAYSTOTG, 
                                                                   ifelse(comorb$NEPHRO == 0, comorb$DAYSTONEPHRO, 
                                                                          ifelse(comorb$HYP ==0, comorb$DAYSTOHYP, NA)))))))))
# one person () has every outcome, so fup_time is missing
comorb$fup_time <- ifelse(comorb$RELEASEID == "65-44824", 4425, comorb$fup_time)

# read in adjudicated medical event dataset
ame <- read.csv('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/TODAY2/AME.csv')

# reformat ame database so it's in the same format as COMORB (i.e., variable for event at baseline, variable for event during followup, and days to event)
# note that some people have multiple events of the same type (e.g., CAD) - for now, will retain all the events to allow the most flexibility later, but for most analyses, may need to use the first (to be parallel with COMORB)
# create new dataset ame_summary with one row per event
arrhythmia <- ame %>% filter(AMENAME == 1) %>% select(releaseid, DAYSTOAME)
colnames(arrhythmia) <- c("RELEASEID", "DAYSTOARRHYTHMIA")
arrhythmia$ARRHYTHMIA <- 1
arrhythmia <- arrhythmia %>% arrange(RELEASEID, DAYSTOARRHYTHMIA) %>% group_by(RELEASEID) %>% filter(row_number()==1)
cad <- ame %>% filter(AMENAME == 2) %>% select(releaseid, DAYSTOAME)
colnames(cad) <- c("RELEASEID", "DAYSTOCAD")
cad$CAD <- 1
cad <- cad %>% arrange(RELEASEID, DAYSTOCAD) %>% group_by(RELEASEID) %>% filter(row_number()==1)
chf <- ame %>% filter(AMENAME == 3) %>% select(releaseid, DAYSTOAME)
colnames(chf) <- c("RELEASEID", "DAYSTOCHF")
chf$CHF <- 1
chf <- chf %>% arrange(RELEASEID, DAYSTOCHF) %>% group_by(RELEASEID) %>% filter(row_number()==1)
LVSD <- ame %>% filter(AMENAME == 4) %>% select(releaseid, DAYSTOAME)
colnames(LVSD) <- c("RELEASEID", "DAYSTOLVSD")
LVSD$LVSD <- 1
LVSD <- LVSD %>% arrange(RELEASEID, DAYSTOLVSD) %>% group_by(RELEASEID) %>% filter(row_number()==1)
mi <- ame %>% filter(AMENAME == 5) %>% select(releaseid, DAYSTOAME)
colnames(mi) <- c("RELEASEID", "DAYSTOMI")
mi$MI <- 1
mi <- mi %>% arrange(RELEASEID, DAYSTOMI) %>% group_by(RELEASEID) %>% filter(row_number()==1)
pad <- ame %>% filter(AMENAME == 6) %>% select(releaseid, DAYSTOAME)
colnames(pad) <- c("RELEASEID", "DAYSTOPAD")
pad$PAD <- 1
pad <- pad %>% arrange(RELEASEID, DAYSTOPAD) %>% group_by(RELEASEID) %>% filter(row_number()==1)

#rad <- ame %>% filter(AMENAME == 7) %>% select(releaseid, DAYSTOAME)
#colnames(rad) <- c("RELEASEID", "DAYSTORAD")
#rad$RAD <- 1
dvt <- ame %>% filter(AMENAME == 8) %>% select(releaseid, DAYSTOAME)
colnames(dvt) <- c("RELEASEID", "DAYSTODVT")
dvt$DVT <- 1
dvt <- dvt %>% arrange(RELEASEID, DAYSTODVT) %>% group_by(RELEASEID) %>% filter(row_number()==1)
stroke <- ame %>% filter(AMENAME == 9) %>% select(releaseid, DAYSTOAME)
colnames(stroke) <- c("RELEASEID", "DAYSTOSTROKE")
stroke$STROKE <- 1
stroke <- stroke %>% arrange(RELEASEID, DAYSTOSTROKE) %>% group_by(RELEASEID) %>% filter(row_number()==1)
#cerebrovascular <- ame %>% filter(AMENAME == 10) %>% select(releaseid, DAYSTOAME)
#colnames(cerebrovascular) <- c("RELEASEID", "DAYSTOCEREBROVASCULAR")
#cerebrovascular$CEREBROVASCULAR <- 1
tia <- ame %>% filter(AMENAME == 11) %>% select(releaseid, DAYSTOAME)
colnames(tia) <- c("RELEASEID", "DAYSTOTIA")
tia$TIA <- 1
tia <- tia %>% arrange(RELEASEID, DAYSTOTIA) %>% group_by(RELEASEID) %>% filter(row_number()==1)

# not sure if we need those events since they might not represent outcomes we are interested for this current CKM analyses
pancreatitis <- ame %>% filter(AMENAME == 12) %>% select(releaseid, DAYSTOAME)
colnames(pancreatitis) <- c("RELEASEID", "DAYSTOPANCREATITIS")
pancreatitis$PANCREATITIS <- 1
pancreatitis <- pancreatitis %>% arrange(RELEASEID, DAYSTOPANCREATITIS) %>% group_by(RELEASEID) %>% filter(row_number()==1)
gallbladder <- ame %>% filter(AMENAME == 13) %>% select(releaseid, DAYSTOAME)
colnames(gallbladder) <- c("RELEASEID", "DAYSTOGALLBLADDER")
gallbladder$GALLBLADDER <- 1
gallbladder <- gallbladder %>% arrange(RELEASEID, DAYSTOGALLBLADDER) %>% group_by(RELEASEID) %>% filter(row_number()==1)
pneuro <- ame %>% filter(AMENAME == 14) %>% select(releaseid, DAYSTOAME)
colnames(pneuro) <- c("RELEASEID", "DAYSTOPNEURO")
pneuro$PNEURO <- 1
pneuro <- pneuro %>% arrange(RELEASEID, DAYSTOPNEURO) %>% group_by(RELEASEID) %>% filter(row_number()==1)
aneuro <- ame %>% filter(AMENAME == 15) %>% select(releaseid, DAYSTOAME)
colnames(aneuro) <- c("RELEASEID", "DAYSTOANEURO")
aneuro$ANEURO <- 1
aneuro <- aneuro %>% arrange(RELEASEID, DAYSTOANEURO) %>% group_by(RELEASEID) %>% filter(row_number()==1)
mononeuro <- ame %>% filter(AMENAME == 16) %>% select(releaseid, DAYSTOAME)
colnames(mononeuro) <- c("RELEASEID", "DAYSTOMONONEURO")
mononeuro$MONONEURO <- 1
mononeuro <- mononeuro %>% arrange(RELEASEID, DAYSTOMONONEURO) %>% group_by(RELEASEID) %>% filter(row_number()==1)

# restart here with CKD in case 
ckd <- ame %>% filter(AMENAME == 17) %>% select(releaseid, DAYSTOAME)
colnames(ckd) <- c("RELEASEID", "DAYSTOCKD")
ckd$CKD <- 1
ckd <- ckd %>% arrange(RELEASEID, DAYSTOCKD) %>% group_by(RELEASEID) %>% filter(row_number()==1)
eskd <- ame %>% filter(AMENAME == 18) %>% select(releaseid, DAYSTOAME)
colnames(eskd) <- c("RELEASEID", "DAYSTOESKD")
eskd$ESKD <- 1
eskd <- eskd %>% arrange(RELEASEID, DAYSTOESKD) %>% group_by(RELEASEID) %>% filter(row_number()==1)

npdr <- ame %>% filter(AMENAME == 19) %>% select(releaseid, DAYSTOAME)
colnames(npdr) <- c("RELEASEID", "DAYSTONPDR")
npdr$NPDR <- 1
npdr <- npdr %>% arrange(RELEASEID, DAYSTONPDR) %>% group_by(RELEASEID) %>% filter(row_number()==1)
pdr <- ame %>% filter(AMENAME == 20) %>% select(releaseid, DAYSTOAME)
colnames(pdr) <- c("RELEASEID", "DAYSTOPDR")
pdr$PDR <- 1
pdr <- pdr %>% arrange(RELEASEID, DAYSTOPDR) %>% group_by(RELEASEID) %>% filter(row_number()==1)
me <- ame %>% filter(AMENAME == 21) %>% select(releaseid, DAYSTOAME)
colnames(me) <- c("RELEASEID", "DAYSTOME")
me$ME <- 1
me <- me %>% arrange(RELEASEID, DAYSTOME) %>% group_by(RELEASEID) %>% filter(row_number()==1)
vh <- ame %>% filter(AMENAME == 22) %>% select(releaseid, DAYSTOAME)
colnames(vh) <- c("RELEASEID", "DAYSTOVH")
vh$VH <- 1
vh <- vh %>% arrange(RELEASEID, DAYSTOVH) %>% group_by(RELEASEID) %>% filter(row_number()==1)
blindness <- ame %>% filter(AMENAME == 23) %>% select(releaseid, DAYSTOAME)
colnames(blindness) <- c("RELEASEID", "DAYSTOBLINDNESS")
blindness$BLINDNESS <- 1
blindness <- blindness %>% arrange(RELEASEID, DAYSTOBLINDNESS) %>% group_by(RELEASEID) %>% filter(row_number()==1)
cataracts <- ame %>% filter(AMENAME == 24) %>% select(releaseid, DAYSTOAME)
colnames(cataracts) <- c("RELEASEID", "DAYSTOCATARACTS")
cataracts$CATARACTS <- 1
cataracts <- cataracts %>% arrange(RELEASEID, DAYSTOCATARACTS) %>% group_by(RELEASEID) %>% filter(row_number()==1)
glaucoma <- ame %>% filter(AMENAME == 25) %>% select(releaseid, DAYSTOAME)
colnames(glaucoma) <- c("RELEASEID", "DAYSTOGLAUCOMA")
glaucoma$GLAUCOMA <- 1
glaucoma <- glaucoma %>% arrange(RELEASEID, DAYSTOGLAUCOMA) %>% group_by(RELEASEID) %>% filter(row_number()==1)


death <- ame %>% filter(AMENAME == 26) %>% select(releaseid, DAYSTOAME)
colnames(death) <- c("RELEASEID", "DAYSTODEATH")
death$DEATH <- 1
death <- death %>% arrange(RELEASEID, DAYSTODEATH) %>% group_by(RELEASEID) %>% filter(row_number()==1)

ame_summary <- full_join(arrhythmia, cad, by = "RELEASEID")
ame_summary <- full_join(ame_summary, chf, by = "RELEASEID")
ame_summary <- full_join(ame_summary, LVSD, by = "RELEASEID")
ame_summary <- full_join(ame_summary, mi, by = "RELEASEID")
ame_summary <- full_join(ame_summary, pad, by = "RELEASEID")
#ame_summary <- full_join(ame_summary, rad, by = "RELEASEID")
ame_summary <- full_join(ame_summary, dvt, by = "RELEASEID")
ame_summary <- full_join(ame_summary, stroke, by = "RELEASEID")
#ame_summary <- full_join(ame_summary, cerebrovascular, by = "RELEASEID")
ame_summary <- full_join(ame_summary, tia, by = "RELEASEID")
ame_summary <- full_join(ame_summary, pancreatitis, by = "RELEASEID")
ame_summary <- full_join(ame_summary, gallbladder, by = "RELEASEID")
ame_summary <- full_join(ame_summary, pneuro, by = "RELEASEID")
ame_summary <- full_join(ame_summary, aneuro, by = "RELEASEID")
ame_summary <- full_join(ame_summary, mononeuro, by = "RELEASEID")
ame_summary <- full_join(ame_summary, ckd, by = "RELEASEID")
ame_summary <- full_join(ame_summary, eskd, by = "RELEASEID")
ame_summary <- full_join(ame_summary, npdr, by = "RELEASEID")
ame_summary <- full_join(ame_summary, pdr, by = "RELEASEID")
ame_summary <- full_join(ame_summary, me, by = "RELEASEID")
ame_summary <- full_join(ame_summary, vh, by = "RELEASEID")
ame_summary <- full_join(ame_summary, blindness, by = "RELEASEID")
ame_summary <- full_join(ame_summary, cataracts, by = "RELEASEID")
ame_summary <- full_join(ame_summary, glaucoma, by = "RELEASEID")
ame_summary <- full_join(ame_summary, death, by = "RELEASEID")

# record NAs to 0
ame_summary <- ame_summary %>% 
    mutate_at(c("ARRHYTHMIA", "CAD", "CHF", "LVSD", "MI", "PAD", "DVT", "STROKE", "TIA", "CKD", "ESKD", "DEATH"), ~replace_na(.,0))
ame_summary <- ame_summary %>% 
    mutate_at(c("ARRHYTHMIA", "CAD", "CHF", "LVSD", "MI", "PAD", "DVT", "STROKE", "TIA", "CKD", "ESKD", "DEATH"), ~as.factor(.))

dim(ame_summary)
ame_summary[,25]
ame_summary <- ame_summary[, -c(20, 25)]

ame_summary$clin_cvd <- 

# label AME types
ame$AMENAME <- factor(ame$AMENAME, labels = c("Arrhythmia", "Coronary artery disease", "Coronary heart failure", "Left ventricular systolic dysfunction", "Myocardial infarction", "PAD/vascular insufficiency", "Deep vein thrombosis", "Stroke", "TIA", "Chronic kidney disease", "ESKD", "Death"))

# categorising clinical CVD
ame_summary$clin_cvd <- ifelse(
  rowSums(ame_summary[, c("ARRHYTHMIA", "CAD", "CHF", "MI", "PAD", "DVT", "STROKE", "TIA", "DEATH")] == 1, na.rm = TRUE) > 0,
  1,
  0)

#joined dataframes
comorb_ckm <- merge(comorb, ame_summary, by = "RELEASEID", all.x = TRUE)

hist(comorb_ckm$clin_cvd)

# categorising subclinical cvd 
## creating a new data frame extracting bnp, troponins and egfr, I decided to keep the baseline levels, the highest value and the last value on each to keep for our further analyses.

# Define function to extract required rows per variable
extract_biomarker_values <- function(df, var) {
  df %>%
    arrange(releaseid, days) %>%
    group_by(releaseid) %>%
    summarise(
      baseline = first(.data[[var]]),
      highest  = max(.data[[var]], na.rm = TRUE),
      latest   = last(.data[[var]]),
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(var, "_", .), -releaseid)
}

# Apply to each biomarker
bnp_vals      <- extract_biomarker_values(addcbl_1, "bnp")
troponin_vals <- extract_biomarker_values(addcbl_1, "troponin")

addcbl_1$testosterone

adiponectin_vals <- extract_biomarker_values(addcbl_1, "adiponectin")
hmwa_vals <-extract_biomarker_values(addcbl_1, "hmwa")
edol_vals <-extract_biomarker_values(addcbl_1, "edol")
shbg_vals <- extract_biomarker_values(addcbl_1, "shbg")
test_vals <- extract_biomarker_values(addcbl_1, "testosterone")

adiponectin_vals <- adiponectin_vals %>% mutate(RELEASEID = toupper(releaseid))
hmwa_vals        <- hmwa_vals        %>% mutate(RELEASEID = toupper(releaseid))
edol_vals        <- edol_vals        %>% mutate(RELEASEID = toupper(releaseid))
shbg_vals        <- shbg_vals        %>% mutate(RELEASEID = toupper(releaseid))
test_vals <- test_vals %>% mutate(RELEASEID = toupper(releaseid))

### This is where I am at the moment with the coding.
comorb_ckm2 <- comorb_ckm2 %>%
  left_join(adiponectin_vals, by = "RELEASEID") %>%
  left_join(hmwa_vals,       by = "RELEASEID") %>%
  left_join(edol_vals,       by = "RELEASEID") %>%
  left_join(shbg_vals,       by = "RELEASEID") %>%
  left_join(test_vals,       by = "RELEASEID")

test_vals <- test_vals %>%
  select(-releaseid)

edol_vals <-edol_vals %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(test_vals, by = "RELEASEID")

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(edol_vals, by = "RELEASEID")

hist(comorb_ckm2$testosterone_baseline, if comorb_ckm2$sex==1)

par(mfrow = c(1,2))

sum(!is.na(comorb_ckm2$testosterone_baseline[comorb_ckm2$sex == 1]))

sum(comorb_ckm2$AGEBASE >= 17, na.rm = TRUE)

hist(
  comorb_ckm2$testosterone_baseline[comorb_ckm2$sex == 2],
  main = "Testosterone Baseline (Males)",
  xlab = "Testosterone",
  col = "lightblue"
)

hist(
  comorb_ckm2$testosterone_baseline[comorb_ckm2$sex == 1],
  main = "Testosterone Baseline (Females)",
  xlab = "Testosterone",
  col = "pink"
)

comorb_ckm2 %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    mean = mean(testosterone_baseline, na.rm = TRUE),
    sd = sd(testosterone_baseline, na.rm = TRUE),
    median = median(testosterone_baseline, na.rm = TRUE),
    p25 = quantile(testosterone_baseline, 0.25, na.rm = TRUE),
    p75 = quantile(testosterone_baseline, 0.75, na.rm = TRUE),
    IQR = IQR(testosterone_baseline, na.rm = TRUE),
    min = min(testosterone_baseline, na.rm = TRUE),
    max = max(testosterone_baseline, na.rm = TRUE)
  )

par(mfrow = c(1,1))  # reset layout


uacr_base <- extract_biomarker_values(baserisk, "UAlbCreat")

uacr_vals <- extract_biomarker_values(cbl_2_unique, "ualbcreat")
ualb_vals <- extract_biomarker_values(cbl_2_unique, "ualb")
ualb_vals <- ualb_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))
uacr_vals <- uacr_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

ualb_vals <- ualb_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)   # optional: drop original lowercase

uacr_vals <- uacr_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(ualb_vals, by = "RELEASEID") %>%
  left_join(uacr_vals, by = "RELEASEID")

il6_vals <- extract_biomarker_values(cbl_2_unique, "il6")

il6_vals <- il6_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

il6_vals <- il6_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(il6_vals, by = "RELEASEID")

il1_vals <- extract_biomarker_values(cbl_2_unique, "il1")

il1_vals <- il1_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

il1_vals <- il1_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(il1_vals, by = "RELEASEID")

hscrp_vals <- extract_biomarker_values(cbl_2_unique, "hscrp")

hscrp_vals <- hscrp_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

hscrp_vals <- hscrp_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(hscrp_vals, by = "RELEASEID")

apob_vals <- extract_biomarker_values(cbl_2_unique, "apob")

apob_vals <- apob_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

apob_vals <- apob_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(apob_vals, by = "RELEASEID")

mcp1_vals <- extract_biomarker_values(cbl_2_unique, "mcp1")

mcp1_vals <- mcp1_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

mcp1_vals <- mcp1_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(mcp1_vals, by = "RELEASEID")

tnfa_vals <- extract_biomarker_values(cbl_2_unique, "tnfa")

tnfa_vals <- tnfa_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

tnfa_vals <- tnfa_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(tnfa_vals, by = "RELEASEID")

tnfr1_vals <- extract_biomarker_values(cbl_2_unique, "tnfr1")

tnfr1_vals <- tnfr1_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

tnfr1_vals <- tnfr1_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(tnfr1_vals, by = "RELEASEID")

tnfr2_vals <- extract_biomarker_values(cbl_2_unique, "tnfr2")

tnfr2_vals <- tnfr2_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

tnfr2_vals <- tnfr2_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(tnfr2_vals, by = "RELEASEID")

hist(cbl_2_unique$vegf)

tnfr2_vals <- extract_biomarker_values(cbl_2_unique, "tnfr2")

tnfr2_vals <- tnfr2_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

tnfr2_vals <- tnfr2_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(tnfr2_vals, by = "RELEASEID")

icam1_vals <- extract_biomarker_values(cbl_2_unique, "icam1")

icam1_vals <- icam1_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

icam1_vals <- icam1_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(icam1_vals, by = "RELEASEID")

vcam1_vals <- extract_biomarker_values(cbl_2_unique, "vcam1")

vcam1_vals <- vcam1_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

vcam1_vals <- vcam1_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(vcam1_vals, by = "RELEASEID")

eselectin_vals <- extract_biomarker_values(cbl_2_unique, "eselectin")

eselectin_vals <- eselectin_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

eselectin_vals <- eselectin_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(eselectin_vals, by = "RELEASEID")

vegf_vals <- extract_biomarker_values(cbl_2_unique, "vegf")

vegf_vals <- vegf_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

vegf_vals <- vegf_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(vegf_vals, by = "RELEASEID")

fgf23_vals <- extract_biomarker_values(cbl_2_unique, "fgf23")

fgf23_vals <- fgf23_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

fgf23_vals <- fgf23_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(fgf23_vals, by = "RELEASEID")

bnp_vals <- bnp_vals %>%
  mutate(across(everything(), ~ifelse(. == -Inf, NA_real_, .)))

bnp_vals <- bnp_vals %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(bnp_vals, by = "RELEASEID")

baseline <- baseline %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

baserisk <- baserisk %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(-releaseid)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(pedsQL_summary, by = "RELEASEID")

# I haven't done this part just yet but I want to explore eGFR as a continous var to keep more of its variance and see how it fluctuates in the different CKM stages.

#for eGFR the extraction should be a little bit different, perhaps keeping the baseline, latest value, highest and lowest ever recoded ? 

min(addcbl_1$eGFR_FAS, na.rm = TRUE)
min(addcbl_1$ckd_gfr, na.rm = TRUE)

max(addcbl_1$eGFR_FAS, na.rm = TRUE)
max(addcbl_1$ckd_gfr, na.rm = TRUE)

hist(addcbl_1$eGFR_FAS)
hist(addcbl_1$ckd_gfr)

gfr_vals      <- extract_biomarker_values(addcbl_1, "ckd_gfr")
fas_vals      <- extract_biomarker_values(addcbl_1, "eGFR_FAS")

data.frame(
  metric = c("eGFR_FAS > 120", "eGFR_FAS > 130", "ckd_gfr > 120", "ckd_gfr > 130"),
  count = c(
    sum(addcbl_1$eGFR_FAS > 120, na.rm = TRUE),
    sum(addcbl_1$eGFR_FAS > 130, na.rm = TRUE),
    sum(addcbl_1$ckd_gfr > 120, na.rm = TRUE),
    sum(addcbl_1$ckd_gfr > 130, na.rm = TRUE)
  )
)

## The other important thing is to keep the hyperfiltration variable per stage as categorical, creating a sum

library(purrr)

# Merge all together on releaseid
biom_ckm <- list(bnp_vals, troponin_vals, gfr_vals, fas_vals) %>%
  reduce(full_join, by = "releaseid")

# After combining all summaries
biom_ckm[biom_ckm == -Inf] <- NA

# Identify columns to check (all except releaseid, mvist, days)
cols_to_check <- setdiff(names(biom_ckm), c("releaseid"))

# Filter out rows where all those columns are NA
biom_ckm <- biom_ckm %>%
  filter(!if_all(all_of(cols_to_check), is.na))

biom_ckm <- biom_ckm %>%
  left_join(baserisk %>% select(releaseid, sex, sex_char), by = "releaseid")

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(pwv_flags, by = "RELEASEID")

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(baserisk %>% select(RELEASEID, sex), by="RELEASEID")

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(baserisk %>% select(RELEASEID, AGEBASE), by="RELEASEID")

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(
    baserisk %>% 
      mutate(RELEASEID = toupper(releaseid)) %>%  # standardize to uppercase if needed
      select(RELEASEID, AGEBASE),
    by = "RELEASEID"
  )

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(baseline %>% select(RELEASEID, bmi, sbp, dbp, wastcirc), by= "RELEASEID")

range(biom_ckm$troponin_highest)

# lowest value =/ 0
min(biom_ckm$troponin_baseline[biom_ckm$troponin_baseline > 0], na.rm = TRUE)
min(biom_ckm$troponin_highest[biom_ckm$troponin_highest > 0], na.rm = TRUE)
min(biom_ckm$troponin_latest[biom_ckm$troponin_latest > 0], na.rm = TRUE)

# for the threshold of subclinical CVD troponin
biom_ckm <- biom_ckm %>%
  mutate(
    high_trop0 = ifelse(!is.na(troponin_baseline) & troponin_baseline > 0, 1, 0),
    high_trop1 = ifelse(
      (!is.na(troponin_highest) & troponin_highest > 0) |
      (!is.na(troponin_latest) & troponin_latest > 0),
      1, 0
    )
  )

# for the threshold of subclinical CVD troponin

names(comorb_ckm2)
names(biom_ckm)

biom_ckm <- biom_ckm %>%
  rename(RELEASEID = releaseid)

## Integration of CKM with biomarkers

comorb_ckm2 <- comorb_ckm %>%
  left_join(
    biom_ckm %>%
      select(RELEASEID, high_trop0, high_trop1, high_bnp0, high_bnp1),
    by = "RELEASEID"
  )

names(biom_ckm)

# create dfs with echoes, PWV, eGFR, BNP
# UACR will come from adjudicated events

# TODAY ECHO
today_echo <- read.csv("/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/TODAY/echo.csv")
today_echo$TIMEPOINT <- "1"
# TODAY2 ECHO
today2_echo <- read.csv("/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/TODAY2/echo.csv")
today2_echo$TIMEPOINT <- "2"
drop_cols <- c("dopplerqc","mmodeqc","overallqc","plaxqc","saxqc","apicalqc","mmodeqc")
today2_echo <- today2_echo %>% select(-one_of(drop_cols))
# combine echo datasets 
echo <- bind_rows(today2_echo, today_echo)
echo$TIMEPOINT <- as.factor(echo$TIMEPOINT)
echo$RELEASEID <- echo$releaseid
echo$releaseid <- NULL
# keep only needed vars, uppercase sensitive
echo_keep <- echo %>% select(RELEASEID, TIMEPOINT, lvsepem, lvem, lvsepratio, lvratio, peakvelo, laarea2d, lvmass, walthick, ivsdias, ivssyst, walldias, wallsyst)
echo_keep$average_E_Em <- (echo_keep$lvsepratio + echo_keep$lvratio) / 2

# SPECKLE TRACKING
# need to keep GLS_4CH
# don't have GLS in TODAY, just TODAY2 - but there are two time points
speckle <- read.csv('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/TODAY2/SPECKLE.csv')
## Adding the GLS 2ch and 3ch view 
speckle_keep <- speckle %>% select(RELEASEID, TIMEPOINT, GLS_4CH, GLS_2CH, GLS_3CH)

# Extract only the numeric part from "SPECKLE1" and "SPECKLE2"
speckle_keep$TIMEPOINT <- gsub("SPECKLE", "", speckle_keep$TIMEPOINT)

#rm(echo_keep)
echo_keep <- merge(echo_keep, speckle_keep,
                          by = c("RELEASEID", "TIMEPOINT"),
                          all.x = TRUE)  # Keep all rows from echo_keep

# Only merge columns that are not already present in echo_keep
speckle_merge <- speckle_keep %>%
  select(-any_of(c("GLS_4CH")))  # remove column if it might cause conflict

# Safe merge without overwriting GLS_4CH in echo_keep
echo_keep <- merge(echo_keep, speckle_merge,
                   by = c("RELEASEID", "TIMEPOINT"),
                   all.x = TRUE)

## What is the highest GLS number at both baseline and f/u.
# Ensure RELEASEID is in uppercase
speckle_keep <- speckle_keep %>%
  mutate(
    # Compute the highest GLS (i.e., least negative) value across the 3 columns
    GLS_highest = pmax(GLS_4CH, GLS_3CH, GLS_2CH, na.rm = TRUE)
  )

# Summarise GLS_highest by TIMEPOINT (one row per RELEASEID and TIMEPOINT)
speckle_summary <- speckle_keep %>%
  distinct(RELEASEID, TIMEPOINT, .keep_all = TRUE) %>%
  group_by(TIMEPOINT) %>%
  summarise(
    mean_GLS   = mean(GLS_highest, na.rm = TRUE),
    sd_GLS     = sd(GLS_highest, na.rm = TRUE),
    min_GLS    = min(GLS_highest, na.rm = TRUE),
    max_GLS    = max(GLS_highest, na.rm = TRUE),
    median_GLS = median(GLS_highest, na.rm = TRUE),
    IQR_GLS    = IQR(GLS_highest, na.rm = TRUE),
    p25_GLS    = quantile(GLS_highest, 0.25, na.rm = TRUE),
    p75_GLS    = quantile(GLS_highest, 0.75, na.rm = TRUE),
    n          = sum(!is.na(GLS_highest)),
    .groups = "drop"
  )

# TODAY2 PWV
pwv <- read.csv("/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward//Clinical data/TODAY2/pwv.csv")
pwv_keep <- pwv %>% select(RELEASEID, TIMEPOINT, PWVF)
pwv_keep$pwv_cf_gt10 <- ifelse(is.na(pwv_keep$PWVF), NA, 
                               ifelse(pwv_keep$PWVF >10, 1, 0))
pwv_keep <- pwv_keep %>% arrange(RELEASEID, desc(pwv_cf_gt10)) %>% group_by(RELEASEID) %>% filter(row_number()==1)

# Create high_pwv0 and high_pwv1 based on TIMEPOINT
pwv_keep$high_pwv0 <- ifelse(pwv_keep$TIMEPOINT == 1, pwv_keep$pwv_cf_gt10, NA)
pwv_keep$high_pwv1 <- ifelse(pwv_keep$TIMEPOINT == 2, pwv_keep$pwv_cf_gt10, NA)

pwv_flags <- pwv_keep %>%
  group_by(RELEASEID) %>%
  summarise(
    high_pwv0 = suppressWarnings(max(high_pwv0, na.rm = TRUE)),
    high_pwv1 = suppressWarnings(max(high_pwv1, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    high_pwv0 = ifelse(is.infinite(high_pwv0), NA, high_pwv0),
    high_pwv1 = ifelse(is.infinite(high_pwv1), NA, high_pwv1)
  )

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(pwv_flags, by = "RELEASEID")

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(
    pwv_keep %>%
      select(RELEASEID, pwv_cf_gt10),
    by = "RELEASEID"
  )

# Extract only the numeric part from "PWV1" and "PWV2"
# to later subcategorise the subclinical CVD staging stage 3

pwv_keep$TIMEPOINT <- gsub("PWV", "", pwv_keep$TIMEPOINT)

pwv_keep$TIMEPOINT <- gsub("PWV", "", pwv_keep$TIMEPOINT)

# Create HFpEF points scoring system
# Initialize the new variable with zeros

echo_keep$HFpEF_points <- NULL

echo_keep$HFpEF_points <- 0

any(echo_keep$peakvelo > 280, na.rm = TRUE)
sum(echo_keep$peakvelo > 280, na.rm = TRUE)

# Loop through each row in the dataset
# Functional segment of the HFA-PEFF score
# The scoring needs to be hierarchical and assign only max 2 points for any major or 2 point for any minor

for (i in 1:nrow(echo_keep)) {

  # Flag to check if any 2-point condition was scored
  scored_2pts <- FALSE

  # 2-point scoring block (OR logic)
  if (!is.na(echo_keep$lvsepem[i]) && echo_keep$lvsepem[i] < 7) {
    echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 2
    scored_2pts <- TRUE
  } else if (!is.na(echo_keep$lvem[i]) && echo_keep$lvem[i] < 10) {
    echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 2
    scored_2pts <- TRUE
  } else if (!is.na(echo_keep$average_E_Em[i]) && echo_keep$average_E_Em[i] >= 15) {
    echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 2
    scored_2pts <- TRUE
  } else if (!is.na(echo_keep$peakvelo[i]) && echo_keep$peakvelo[i] > 280) {
    echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 2
    scored_2pts <- TRUE
  }

  # 1-point scoring block (only if no 2-point condition was met)
  if (!scored_2pts) {
    if (
      (!is.na(echo_keep$GLS_4CH[i]) && echo_keep$GLS_4CH[i] > -16) ||
      (!is.na(echo_keep$GLS_3CH[i]) && echo_keep$GLS_3CH[i] > -16) ||
      (!is.na(echo_keep$GLS_2CH[i]) && echo_keep$GLS_2CH[i] > -16) ||
      (!is.na(echo_keep$average_E_Em[i]) && echo_keep$average_E_Em[i] >= 9 && echo_keep$average_E_Em[i] < 15)
    ) {
      echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 1
    }
  }
}

# Summarize the distribution of HFpEF_points
table(echo_keep$HFpEF_points)
summary(echo_keep$HFpEF_points)

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(echo_hfpef2, by = "RELEASEID")

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(echo_hfpef, by = "RELEASEID")

table(comorb_ckm2$hfpef1, useNA = "ifany")
table(comorb_ckm2$hfpef2, useNA = "ifany")

# Function to get the maximum BNP value, ignoring NAs
max_bnp <- function(x) {
  if(all(is.na(x))) return(NA)
  return(max(x, na.rm = TRUE))
}

# Get maximum BNP per RELEASEID
max_bnp_by_id <- aggregate(addcbl_1$bnp ~ addcbl_1$releaseid, data = addcbl_1, FUN = max_bnp)

# Explicitly set the column names to ensure they are correct
colnames(max_bnp_by_id) <- c("RELEASEID", "bnp")

# Step 2: Merge this information with echo_keep
echo_keep <- merge(echo_keep, max_bnp_by_id, 
                           by = "RELEASEID", 
                           all.x = TRUE)  # Keep all echo_keep records

echo_keep <- echo_keep %>%
  select(-bnp.y) %>%       # Remove duplicate
  rename(bnp = bnp.x)      # Rename bnp.x to bnp

sum(echo_keep$bnp >= 80)
sum(echo_keep$bnp >= 80, na.rm = TRUE)
sum(echo_keep$bnp > 80, na.rm = TRUE)

# Apply BNP scoring rules directly to echo_keep
for(i in 1:nrow(echo_keep)) {
  # Skip if BNP is NA
  if(!is.na(echo_keep$bnp[i])) {
    if(echo_keep$bnp[i] >= 80) {
      # Add 2 points if BNP > 80
      echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 2
    } else if(echo_keep$bnp[i] >= 35 && echo_keep$bnp[i] < 80) {
      # Add 1 point if BNP between 35-80
      echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 1
    }
    # No points added if BNP < 35
  }
}

hist(echo_keep$HFpEF_points)

# Apply scoring: 1 point if ANY of ivsdias, ivssyst, walldiast, or wallsyst > 1.2
for(i in 1:nrow(echo_keep)) {
  # Check if ANY of the four variables exceed 1.2
  any_exceed_threshold <- FALSE
  
  if(!is.na(echo_keep$ivsdias[i]) && echo_keep$ivsdias[i] > 1.2) {
    any_exceed_threshold <- TRUE
  } else if(!is.na(echo_keep$ivssyst[i]) && echo_keep$ivssyst[i] > 1.2) {
    any_exceed_threshold <- TRUE
  } else if(!is.na(echo_keep$walldias[i]) && echo_keep$walldias[i] > 1.2) {
    any_exceed_threshold <- TRUE
  } else if(!is.na(echo_keep$wallsyst[i]) && echo_keep$wallsyst[i] > 1.2) {
    any_exceed_threshold <- TRUE
  }
  
  # Add 1 point if any variable exceeds threshold
  if(any_exceed_threshold) {
    echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 1
  }
}

## Only if the scoring system qualifies for diastolic measures instead of systolic measurements
# Apply scoring: 1 point if ANY of ivsdias, walldiast > 1.2

# Complete morphological criteria

echo_keep <- echo_keep %>%
  mutate(
    morph_score = case_when(
      # 2-point conditions
      LAVI_est > 34 ~ 2,
      ((sex == 2 & lvmi >= 149) | (sex == 1 & lvmi >= 122)) & walthick > 0.42 ~ 2,

      # 1-point conditions (only if 2-point not met)
      between(LAVI_est, 29, 34) ~ 1,
      sex == 2 & lvmi >= 115 ~ 1,
      sex == 1 & lvmi >= 95 ~ 1,
      walthick > 0.42 ~ 1,
      ivsdias > 1.2 ~ 1,
      walldias > 1.2 ~ 1,

      TRUE ~ 0
    ),
    
    # Add to HFpEF points
    HFpEF_points = HFpEF_points + morph_score
  )


#for (i in 1:nrow(echo_keep)) {
  # Check if any wall-thickness criteria are met
 # if (
  #  (!is.na(echo_keep$ivsdias[i])   && echo_keep$ivsdias[i] > 1.2) ||
   # (!is.na(echo_keep$walldias[i])  && echo_keep$walldias[i] > 1.2) ||
#    (!is.na(echo_keep$walthick[i])  && echo_keep$walthick[i] > 0.42)
#  ) {
#    echo_keep$HFpEF_points[i] <- echo_keep$HFpEF_points[i] + 1
  }
}

# Create a variable for hfpef diagnosis, >= 5 points from the algorithm 
echo_keep$hfpef <-0
echo_keep$hfpef <- ifelse(echo_keep$HFpEF_points >= 5, 1, 0)
sum(echo_keep$hfpef == 1, na.rm = TRUE)

table(echo_keep$hfpef == 1, echo_keep$TIMEPOINT)

# Create hfpef1 and hfpef2 based on TIMEPOINT
echo_keep$hfpef1 <- ifelse(echo_keep$TIMEPOINT == 1, echo_keep$hfpef, NA)
echo_keep$hfpef2 <- ifelse(echo_keep$TIMEPOINT == 2, echo_keep$hfpef, NA)

length(unique(echo_keep$RELEASEID))

echo_hfpef <- echo_keep %>%
  filter(TIMEPOINT == 1) %>%                     # keep only timepoint 1 (for hfpef1)
  select(RELEASEID, hfpef1) %>%                  # keep only relevant columns
  filter(!is.na(hfpef1)) %>%                     # remove rows where hfpef1 is NA
  distinct(RELEASEID, .keep_all = TRUE)          # keep only one row per RELEASEID

echo_hfpef2 <- echo_keep %>%
filter(TIMEPOINT == 2) %>%                     # keep only timepoint 2
select(RELEASEID, hfpef2) %>%                  # keep only relevant columns
filter(!is.na(hfpef2)) %>%                     # remove rows where hfpef2 is NA
distinct(RELEASEID, .keep_all = TRUE)          # keep one row per RELEASEID

#CKM Staging:
# the stages include either baseline or outcome measurements and either HFpEF 1 or 2 
comorb_ckm2 <- comorb_ckm2 %>%
  mutate(
    CKM_syn = case_when(
      # Stage 4: most severe
      clin_cvd == 1 |
      CKD == 1 |
      hfpef1 == 1 | hfpef2 == 1 |
      MAC0 == 1 | MAC == 1 ~ "Stage4",

      # Stage 3: not Stage 4, but elevated biomarkers or moderate risk
      (high_bnp0 == 1 | high_bnp1 == 1 |
       high_trop0 == 1 | high_trop1 == 1 |
       high_pwv0 == 1 | high_pwv1 == 1 |
       CKD == 1 |
       MAC0 == 1 | MAC == 1) ~ "Stage3",

      # Stage 2 plus: not Stage 3 or 4, but classic risk factors
      (HTN0 == 1 | HTN == 1 |
       ANYDLP0 == 1 | ANYDLP == 1 |
       MAC0 == 1 | MAC == 1) ~ "Stage2_plus",

      # All else falls into Stage 2
      TRUE ~ "Stage2"
    )
  )

table(comorb_ckm2$CKM_syn, useNA = "ifany")

library(ggplot2)

ggplot(comorb_ckm2, aes(x = CKM_syn)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Distribution of CKM Stages",
    x = "CKM Stage",
    y = "Number of Patients"
  ) +
  theme_minimal(base_size = 14)

kruskal.test(HYP0 ~ CKM_syn, data = comorb_ckm2)
kruskal.test(HYP ~ CKM_syn, data = comorb_ckm2)

install.packages("FSA")
library(FSA)

dunn.test(HYP0 ~ CKM_syn, data = comorb_ckm2, method = "holm")
dunn.test(HYP ~ CKM_syn, data = comorb_ckm2, method = "holm")

ggplot(comorb_ckm2, aes(x = CKM_syn, y = HYP0)) +
  geom_boxplot(fill = "lightgray") +
  labs(title = "Distribution of HYP0 across CKM Stages", x = "CKM Stage", y = "HYP0") +
  theme_minimal()

comorb_ckm2$CKM_syn <- factor(comorb_ckm2$CKM_syn, levels = c("Stage2", "Stage2_plus", "Stage3", "Stage4"))

# Logistic regression for HYP0
model_hyp0 <- glm(HYP0 ~ CKM_syn, data = comorb_ckm2, family = binomial)

# Logistic regression for HYP
model_hyp <- glm(HYP ~ CKM_syn, data = comorb_ckm2, family = binomial)

summary(model_hyp0)
summary(model_hyp)

# How many participants presented HYP by each of the stages ?
# Count of HYP0 == 1 by CKM stage
comorb_ckm2 %>%
  group_by(CKM_syn) %>%
  summarise(
    count_HYP0_1 = sum(HYP0 == 1, na.rm = TRUE),
    total = n(),
    percent_HYP0_1 = round(100 * count_HYP0_1 / total, 1)
  )

# Count of HYP == 1 by CKM stage
comorb_ckm2 %>%
  group_by(CKM_syn) %>%
  summarise(
    count_HYP_1 = sum(HYP == 1, na.rm = TRUE),
    total = n(),
    percent_HYP_1 = round(100 * count_HYP_1 / total, 1)
  )

# Create contingency table
table_hyp <- table(comorb_ckm2$CKM_syn, comorb_ckm2$HYP)
print(table_hyp)
print(table_hyp0)
# Run chi-squared test
chisq.test(table_hyp)

fisher.test(table_hyp)

pairwise.prop.test(
  x = table_hyp[, 2],  # counts of HYP == 1
  n = rowSums(table_hyp),  # total per group
  p.adjust.method = "holm"
)

##joining stage3 and stage4, since both are CVD
library(dplyr)

comorb_ckm2 <- comorb_ckm2 %>%
  mutate(
    CKM_syn_collapsed = case_when(
      CKM_syn %in% c("Stage3", "Stage4") ~ "Adv_CVDCKM",
      CKM_syn == "Stage2_plus" ~ "Stage2_plus",
      CKM_syn == "Stage2" ~ "Stage2"
    )
  )

# Contingency tables
table_hyp0 <- table(comorb_ckm2$CKM_syn_collapsed, comorb_ckm2$HYP0)
table_hyp  <- table(comorb_ckm2$CKM_syn_collapsed, comorb_ckm2$HYP)

print(table_hyp)
print(table_hyp0)

# For HYP0
chisq.test(table_hyp0)

# For HYP
chisq.test(table_hyp)

fisher.test(table_hyp0)
fisher.test(table_hyp)


# Function to run McNemar's test by CKM stage
run_mcnemar <- function(stage) {
  tab <- comorb_ckm2 %>%
    filter(CKM_syn_collapsed == stage) %>%
    select(HYP0, HYP) %>%
    table(useNA = "no")
  
  cat("\nMcNemar's test for", stage, "\n")
  print(mcnemar.test(tab))
}

# Run for each stage
unique(comorb_ckm2$CKM_syn_collapsed) %>% lapply(run_mcnemar)

comorb_ckm2 %>%
  group_by(CKM_syn_collapsed) %>%
  summarise(
    prop_HYP0 = mean(HYP0 == 1, na.rm = TRUE),
    prop_HYP  = mean(HYP == 1, na.rm = TRUE),
    n = n()
  )

# Loop through each stage and print the 2x2 table
unique(comorb_ckm2$CKM_syn_collapsed) %>%
  lapply(function(stage) {
    cat("\n===== CKM Stage:", stage, "=====\n")
    
    # Filter for the current stage
    tab <- comorb_ckm2 %>%
      filter(CKM_syn_collapsed == stage) %>%
      select(HYP0, HYP) %>%
      table(useNA = "no")
    
    print(tab)
    
    # Optionally run McNemar's test
    if (all(dim(tab) == c(2, 2))) {
      cat("\nMcNemar's test p-value:\n")
      print(mcnemar.test(tab)$p.value)
    } else {
      cat("\n⚠️ Not a 2x2 table (may contain missing or invalid values)\n")
    }
  })

#CKM Staging: creating 2 stagings for the baseline point and for the prospective collection of endpoints

comorb_ckm2 <- comorb_ckm2 %>%
  mutate(

    # Baseline CKM staging (no clin_cvd or CKD, only baseline variables)
    CKM_syn_base = case_when(
      # Stage 4: severe features (baseline only)
      # hfpef1 == 1 |
      MAC0 == 1 ~ "Stage4",

      # Stage 3: elevated biomarkers (baseline only)
      high_bnp0 == 1 |
      high_trop0 == 1 |
      high_pwv0 == 1 |
      MAC0 == 1 ~ "Stage3",

      # Stage 2_plus: risk factors only
      HTN0 == 1 |
      ANYDLP0 == 1 |
      MAC0 == 1 ~ "Stage2_plus",

      # If none of the above, it's Stage 2
      TRUE ~ "Stage2"
    )
)
table(comorb_ckm2$CKM_syn_base, useNA = "ifany")

comorb_ckm2 <- comorb_ckm2 %>%
  mutate(
    CKM_syn_fu = case_when(
      # Stage 4: includes clinical disease and follow-up HFpEF or complications
      clin_cvd == 1 |
      CKD == 1 |
      hfpef2 == 1 |
      MAC == 1 ~ "Stage4",

      # Stage 3: elevated follow-up biomarkers or MAC
      high_bnp1 == 1 |
      high_trop1 == 1 |
      high_pwv1 == 1 |
      MAC == 1 ~ "Stage3",

      # Stage 2_plus: classic follow-up risk factors
      HTN == 1 |
      ANYDLP == 1 |
      MAC == 1 ~ "Stage2_plus",

      # Default: Stage 2
      TRUE ~ "Stage2"
    )
  )

table(comorb_ckm2$CKM_syn_fu, useNA = "ifany")

comorb_ckm2 <- comorb_ckm2 %>%
  mutate(

    # Baseline CKM staging (no clin_cvd or CKD, only baseline variables)
    CKM_syn_base = case_when(
      # Stage 4: severe features (baseline only)
      # hfpef1 == 1 |
      MAC0 == 1 ~ "Stage4",

      # Stage 3: elevated biomarkers (baseline only)
      high_bnp0 == 1 |
      high_trop0 == 1 |
      high_pwv0 == 1 |
      MAC0 == 1 ~ "Stage3",

      # Stage 2_plus: risk factors only
      HTN0 == 1 |
      ANYDLP0 == 1 |
      MAC0 == 1 ~ "Stage2_plus",

      # If none of the above, it's Stage 2
      TRUE ~ "Stage2"
    )
)

## Another CKM category, imputing patients with BMI separating the 2 groups. 

comorb_ckm2$bmi

comorb_ckm2 <- comorb_ckm2 %>%
  mutate(
    CKM_syn2 = case_when(
      # Stage 4: most severe
      clin_cvd == 1 |
      CKD == 1 |
      hfpef1 == 1 | hfpef2 == 1 |
      MAC0 == 1 | MAC == 1 ~ "Stage4",

      # Stage 3: not Stage 4, but elevated biomarkers or moderate risk
      (high_bnp0 == 1 | high_bnp1 == 1 |
       high_trop0 == 1 | high_trop1 == 1 |
       high_pwv0 == 1 | high_pwv1 == 1 |
       CKD == 1 |
       MAC0 == 1 | MAC == 1) ~ "Stage3",

      # Stage 2 plus: not Stage 3 or 4, but classic risk factors
      (HTN0 == 1 | HTN == 1 |
       ANYDLP0 == 1 | ANYDLP == 1 |
       MAC0 == 1 | MAC == 1 | bmi >= 30) ~ "Stage2_plus",

      # All else falls into Stage 2
      TRUE ~ "Stage2"
    )
  )

table(comorb_ckm2$CKM_syn2, useNA = "ifany")

comorb_ckm2 <- comorb_ckm2 %>%
  mutate(
    # Baseline CKM staging (no clin_cvd or CKD, only baseline variables)
    CKM_syn_base2 = case_when(
      # Stage 4: severe features
      MAC0 == 1 ~ "Stage4",

      # Stage 3: elevated biomarkers
      high_bnp0 == 1 |
      high_trop0 == 1 |
      high_pwv0 == 1 |
      MAC0 == 1 ~ "Stage3",

      # Stage 2_plus: risk factors OR BMI >=30
      HTN0 == 1 |
      ANYDLP0 == 1 |
      MAC0 == 1 |
      bmi >= 30 ~ "Stage2_plus",

      # If none of the above, Stage2
      TRUE ~ "Stage2"
    )
  )

table(comorb_ckm2$CKM_syn_base2, useNA = "ifany")

library(dplyr)

# Baseline: % by stage
baseline_pct <- comorb_ckm2 %>%
  group_by(CKM_syn_base) %>%
  summarise(
    count = n()
  ) %>%
  mutate(
    total = sum(count),
    percent = round(100 * count / total, 1),
    timepoint = "Baseline"
  )

# Follow-up: % by stage
followup_pct <- comorb_ckm2 %>%
  group_by(CKM_syn_fu) %>%
  summarise(
    count = n()
  ) %>%
  mutate(
    total = sum(count),
    percent = round(100 * count / total, 1),
    timepoint = "Follow-up"
  )

# Combine baseline and follow-up results
ckm_stage_pct <- bind_rows(
  baseline_pct %>% rename(stage = CKM_syn_base),
  followup_pct %>% rename(stage = CKM_syn_fu)
)

# View summary
ckm_stage_pct

# Combined percentage at baseline
combined_baseline <- comorb_ckm2 %>%
  filter(CKM_syn_base %in% c("Stage3", "Stage4")) %>%
  summarise(
    count = n(),
    total = nrow(comorb_ckm2),
    percent = round(100 * count / total, 1)
  ) %>%
  mutate(timepoint = "Baseline")

# Combined percentage at follow-up
combined_followup <- comorb_ckm2 %>%
  filter(CKM_syn_fu %in% c("Stage3", "Stage4")) %>%
  summarise(
    count = n(),
    total = nrow(comorb_ckm2),
    percent = round(100 * count / total, 1)
  ) %>%
  mutate(timepoint = "Follow-up")

# Combine both
combined_stages <- bind_rows(combined_baseline, combined_followup)

# View
combined_stages

## QoL questionnaire, keeping only the physical activity related variables. 
library(dplyr)

# Define the columns to keep
cols_to_keep <- c("RELEASEID", "DAYS", "G01WALK", "G02RUN", "G03SPORT",
                  "G04LIFT", "G05BATH", "G06CHORE", "G07HURT", "G08ENERG")

# Select and combine all three datasets
pedsQL_all <- bind_rows(
  pedsQL_Y %>% select(all_of(cols_to_keep)),
  pedsQL_T %>% select(all_of(cols_to_keep)),
  pedsQL_A %>% select(all_of(cols_to_keep))
)

# Define the correct PedsQL item columns
pedsQL_items <- c("G01WALK", "G02RUN", "G03SPORT", "G04LIFT", 
                  "G05BATH", "G06CHORE", "G07HURT", "G08ENERG")

## sum all the values from the physical activity
pedsQL_all <- pedsQL_all %>%
  rowwise() %>%
  mutate(physical_total = sum(c_across(all_of(pedsQL_items)), na.rm = TRUE)) %>%
  ungroup()

#Keep the first, maximum and latest scores for the pedsQL scoring system.

# Ensure DAYS is numeric (if it's not already)
pedsQL_all <- pedsQL_all %>%
  mutate(DAYS = as.numeric(DAYS))

# Get baseline (first), latest, and highest per RELEASEID
pedsQL_summary <- pedsQL_all %>%
  group_by(RELEASEID) %>%
  summarise(
    pedsql_baseline = physical_total[which.min(DAYS)],
    days_baseline = DAYS[which.min(DAYS)],
    
    pedsql_last = physical_total[which.max(DAYS)],
    days_last = DAYS[which.max(DAYS)],
    
    pedsql_highest = max(physical_total, na.rm = TRUE),
    days_highest = DAYS[which.max(physical_total)],
    
    .groups = "drop"
  )

hist(pedsQL_summary$pedsql_baseline)
hist(pedsQL_summary$pedsql_last)
hist(pedsQL_summary$pedsql_highest)

## Adding the height to estimate the BSA

library(dplyr)
library(stringr)

# Step 1: Clean and prepare visit data
long_clean <- long %>%
  mutate(
    RELEASEID = toupper(releaseid),
    visit_type = substr(visit, 1, 1),
    visit_num = as.numeric(str_extract(visit, "\\d+"))
  )

# Step 2: Get latest non-NA height for M visits (TIMEPOINT 1)
height_M <- long_clean %>%
  filter(visit_type == "M", !is.na(height)) %>%
  group_by(RELEASEID) %>%
  slice_max(order_by = visit_num, n = 1, with_ties = FALSE) %>%
  select(RELEASEID, height) %>%
  rename(height_M = height)

# Step 3: Get latest non-NA height for P visits (TIMEPOINT 2)
height_P <- long_clean %>%
  filter(visit_type == "P", !is.na(height)) %>%
  group_by(RELEASEID) %>%
  slice_max(order_by = visit_num, n = 1, with_ties = FALSE) %>%
  select(RELEASEID, height) %>%
  rename(height_P = height)

# Step 4: Join to echo_keep and assign height based on TIMEPOINT
echo_keep <- echo_keep %>%
  mutate(RELEASEID = toupper(RELEASEID)) %>%
  left_join(height_M, by = "RELEASEID") %>%
  left_join(height_P, by = "RELEASEID") %>%
  mutate(
    height_final = case_when(
      TIMEPOINT == 1 ~ height_M,
      TIMEPOINT == 2 ~ height_P,
      TRUE ~ NA_real_
    )
  ) %>%
  select(-height_M, -height_P)  # clean up

library(dplyr)

# Keep only one row per unique `days` in cbl_2 (e.g., the first)
cbl_2_unique <- cbl_2 %>%
  filter(!is.na(pvisit)) %>%
  group_by(days) %>%
  slice(1) %>%  # or use slice_max() if you want the last or largest visit
  ungroup()

# Join to echo
echo <- echo %>%
  left_join(
    cbl_2_unique %>% select(days, pvisit) %>% rename(visit = pvisit),
    by = "days"
  )

library(dplyr)

# Step 1: Deduplicate cbl_1 by days (keep one mvist per day)
cbl_1_unique <- cbl_1 %>%
  filter(!is.na(mvisit)) %>%
  group_by(days) %>%
  slice(1) %>%
  ungroup()

# Step 2: Join with echo, adding mvist as 'visit'
echo <- echo %>%
  left_join(
    cbl_1_unique %>% select(days, mvisit) %>% rename(visit = mvisit),
    by = "days"
  )

# Install if not already
install.packages("fuzzyjoin")

library(dplyr)
library(fuzzyjoin)

# Step 1: Filter only echo rows where both visit.x and visit.y are NA
echo_missing_visits <- echo %>%
  filter(is.na(visit.x) & is.na(visit.y)) %>%
  select(RELEASEID, days) %>%
  filter(!is.na(RELEASEID), !is.na(days))  # ensure valid keys

# Step 2: Prepare `long` with standardized uppercase RELEASEID and clean data
long_clean <- long %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  filter(!is.na(visit_days), !is.na(visit))

# Ensure numeric days in both dataframes
echo_missing_visits <- echo_missing_visits %>%
  mutate(days = as.numeric(days))

long_clean <- long_clean %>%
  mutate(visit_days = as.numeric(visit_days))

library(dplyr)

# Step 1: Prepare data
echo_missing_visits <- echo %>%
  filter(is.na(visit.x) & is.na(visit.y)) %>%
  select(RELEASEID, days) %>%
  mutate(days = as.numeric(days))

long_clean <- long %>%
  filter(!is.na(visit_days), !is.na(visit)) %>%
  mutate(
    RELEASEID = toupper(releaseid),
    visit_days = as.numeric(visit_days)
  ) %>% ()
  select(RELEASEID, visit_days, visit)

# Step 2: Perform manual closest join
closest_visit_fill <- echo_missing_visits %>%
  rowwise() %>%
  mutate(
    visit_closest = {
      sub_long <- long_clean %>% filter(RELEASEID == RELEASEID)
      if (nrow(sub_long) == 0) NA_character_
      else {
        closest_row <- sub_long[which.min(abs(sub_long$visit_days - days)), ]
        closest_row$visit
      }
    }
  ) %>%
  ungroup()

echo <- echo %>%
  left_join(closest_visit_fill, by = c("RELEASEID", "days")) %>%
  mutate(
    visit_final = coalesce(visit.x, visit.y, visit_closest)
  ) %>%
  select(-visit_closest)

# Ensure correct formats
echo <- echo %>%
  mutate(RELEASEID = toupper(RELEASEID),
         days = as.numeric(days))

long_clean <- long %>%
  mutate(
    releaseid = toupper(releaseid),  # match echo's case
    visit_days = as.numeric(visit_days)
  ) %>%
  filter(!is.na(visit_days), !is.na(visit))

# Match closest visit_days per RELEASEID and days
visit_match <- echo %>%
  rowwise() %>%
  mutate(
    visit2 = {
      sub_long <- long_clean %>% filter(releaseid == RELEASEID)
      if (nrow(sub_long) == 0) NA_character_
      else {
        closest_row <- sub_long[which.min(abs(sub_long$visit_days - days)), ]
        closest_row$visit
      }
    }
  ) %>%
  ungroup()

long_clean <- long %>%
  mutate(
    releaseid = toupper(releaseid),  # match echo's case
    visit_days = as.numeric(visit_days)
  ) %>%
  filter(!is.na(visit_days), !is.na(visit))

visit_match <- echo %>%
  rowwise() %>%
  mutate(
    visit2 = {
      sub_long <- long_clean %>% filter(releaseid == RELEASEID)
      if (nrow(sub_long) == 0) NA_character_
      else {
        closest_row <- sub_long[which.min(abs(sub_long$visit_days - days)), ]
        closest_row$visit
      }
    }
  ) %>%
  ungroup()

# Add visit2 to echo by matching on RELEASEID and TIMEPOINT
echo <- echo %>%
  left_join(
    visit_match %>% select(RELEASEID, TIMEPOINT, visit2),
    by = c("RELEASEID", "TIMEPOINT")
  )

echo_keep <- echo_keep %>%
  left_join(
    echo %>% select(RELEASEID, TIMEPOINT, visit2),
    by = c("RELEASEID", "TIMEPOINT")
  )

library(dplyr)

# Step 1: Prepare long with uppercase releaseid
long_clean <- long %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(RELEASEID, visit, height)

# Step 2: Join height into echo_keep by RELEASEID and visit2
echo_keep <- echo_keep %>%
  left_join(
    long_clean,
    by = c("RELEASEID", "visit2" = "visit")
  ) %>%
  rename(height2 = height)  # Optional: rename joined height

sum(is.na(echo_keep$height2))

echo_keep <- echo_keep %>%
  mutate(height2 = ifelse(is.na(height2), height_final, height2))

hist(echo_keep$height2)

# Step 1: Clean and prepare `long` with uppercase RELEASEID
long_clean <- long %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(RELEASEID, visit, visit_days, weight) %>%
  filter(!is.na(visit))

# Step 2: Join by RELEASEID and visit2 to get direct matches
echo_keep <- echo_keep %>%
  left_join(
    long_clean %>% select(RELEASEID, visit, weight),
    by = c("RELEASEID", "visit2" = "visit")
  ) %>%
  rename(weight = weight)

sum(is.na(echo_keep$weight))

echo_keep <- echo_keep %>%
  left_join(
    echo %>% select(RELEASEID, TIMEPOINT, days),
    by = c("RELEASEID", "TIMEPOINT")
  )

# Step 3: Fill NA weight using closest visit_days from `long`
# Prepare list of echo_keep with missing weights
echo_missing_weight <- echo_keep %>%
  filter(is.na(weight)) %>%
  select(RELEASEID, days)

# Match by closest visit_days for those with missing weight
long_weight_fill <- long_clean %>%
  filter(!is.na(weight), !is.na(visit_days)) %>%
  mutate(visit_days = as.numeric(visit_days))

echo_missing_weight <- echo_missing_weight %>%
  mutate(days = as.numeric(days))

# For each missing row, find closest weight from long
weight_closest <- echo_missing_weight %>%
  rowwise() %>%
  mutate(weight_filled = {
    sub_long <- long_weight_fill %>% filter(RELEASEID == RELEASEID)
    if (nrow(sub_long) == 0) NA_real_
    else {
      sub_long$weight[which.min(abs(sub_long$visit_days - days))]
    }
  }) %>%
  ungroup()

# Step 4: Fill in missing weights
echo_keep <- echo_keep %>%
  left_join(weight_closest, by = c("RELEASEID", "days")) %>%
  mutate(weight = coalesce(weight, weight_filled)) %>%
  select(-weight_filled)

# Du Bois & Du Bois: BSA = 0.007184 × Height(cm)^0.725 × Weight(kg)^0.425
echo_keep <- echo_keep %>%
  mutate(
    BSA_dubois = ifelse(!is.na(height2) & !is.na(weight),
                        0.007184 * (height2 ^ 0.725) * (weight ^ 0.425),
                        NA_real_)
  )

# Mosteller: BSA = sqrt((Height(cm) × Weight(kg)) / 3600)
echo_keep <- echo_keep %>%
  mutate(
    BSA_mosteller = ifelse(!is.na(height2) & !is.na(weight),
                           sqrt((height2 * weight) / 3600),
                           NA_real_)
  )

hist(echo_keep$BSA_dubois)
hist(echo_keep$BSA_mosteller)

echo_keep <- echo_keep %>%
  mutate(releaseid = tolower(RELEASEID)) %>%  # create lowercase match column
  left_join(
    baserisk %>% select(releaseid, sex),
    by = "releaseid"
  ) %>%
  select(-releaseid)  # optional: remove the helper column

echo_keep <- echo_keep %>%
  mutate(lvmi = ifelse(!is.na(lvmass) & !is.na(BSA_dubois),
                       lvmass / BSA_dubois,
                       NA_real_))

hist(echo_keep$lvmi)

echo_keep <- echo_keep %>%
  mutate(lvmi_score = case_when(
    sex == 2 & lvmi >= 149 ~ 2,
    sex == 2 & lvmi >= 122 ~ 1,
    sex == 2 ~ 0,
    sex == 1 & lvmi >= 115 ~ 2,
    sex == 1 & lvmi >= 95 ~ 1,
    sex == 1 ~ 0,
    TRUE ~ NA_real_
  ))

echo_keep <- echo_keep %>%
  mutate(
    lvmi_score = case_when(
      sex == 2 & lvmi >= 149 ~ 2,
      sex == 2 & lvmi >= 122 ~ 1,
      sex == 2 ~ 0,
      sex == 1 & lvmi >= 115 ~ 2,
      sex == 1 & lvmi >= 95  ~ 1,
      sex == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    HFpEF_points = HFpEF_points + coalesce(lvmi_score, 0)
  )

hist(echo_keep$lvmi_score)

hist(echo_keep$HFpEF_points)
sum(echo_keep$HFpEF_points >= 5) # 9 # 56
sum(echo_keep$HFpEF_points == 4) # 27 # 143
sum(echo_keep$HFpEF_points == 3) # 27 # 241
sum(echo_keep$HFpEF_points >= 2 & echo_keep$HFpEF_points <= 4, na.rm = TRUE) # 419 # 657

echo_keep <- echo_keep %>%
  left_join(
    echo %>% select(RELEASEID, ladimen),
    by = "RELEASEID"
  )

echo_keep <- echo_keep %>%
  distinct(RELEASEID, TIMEPOINT, .keep_all = TRUE)

## estimating the volumetric LAVI from the 2D measurements in the algorithm

echo_keep <- echo_keep %>%
  mutate(
    LA_volume_est = ifelse(!is.na(laarea2d) & !is.na(ladimen),
                           (8 / (3 * pi)) * (laarea2d^2 / ladimen),
                           NA_real_),
    LAVI_est = ifelse(!is.na(LA_volume_est) & !is.na(BSA_dubois),
                      LA_volume_est / BSA_dubois,
                      NA_real_)
  )

hist(echo_keep$LAVI_est)

library(dplyr)

comorb_ckm2 <- comorb_ckm2 %>%
  select(-CKM_syn_collapsed, -CKM_syn, -CKM_syn_base, -CKM_syn_fu, -hfpef1, -hfpef2)

comorb_ckm2
biom_ckm$eGFR_FAS_baseline biom_ckm$eGFR_FAS_highest biom_ckm$eGFR_FAS_latest biom_ckm$ckd_gfr_baseline biom_ckm$ckd_gfr_highest biom_ckm$ckd_gfr_latest

comorb_ckm2 <- comorb_ckm2 %>%
  left_join(
    biom_ckm %>%
      select(
        RELEASEID,
        eGFR_FAS_baseline,
        eGFR_FAS_highest,
        eGFR_FAS_latest,
        ckd_gfr_baseline,
        ckd_gfr_highest,
        ckd_gfr_latest
      ),
    by = "RELEASEID"
  )

library(ggplot2)

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn), y = eGFR_FAS_highest)) +
  geom_boxplot(fill = "skyblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, color = "darkblue") +
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red") +
  labs(
    x = "CKM Stage (CKM_syn)",
    y = "eGFR_FAS_highest (mL/min/1.73m²)",
    title = "Distribution of eGFR_FAS_highest by CKM Stage"
  ) +
  theme_minimal()

install.packages("ggpubr")  # Run only once
library(ggpubr)

ggboxplot(comorb_ckm2, x = "CKM_syn", y = "eGFR_FAS_highest",
          color = "CKM_syn", palette = "jco",
          add = "jitter", shape = "CKM_syn") +
  stat_compare_means(method = "kruskal.test", label.y = max(comorb_ckm2$eGFR_FAS_highest, na.rm = TRUE) + 10) +  # global p-value
  labs(
    title = "eGFR_FAS_highest Across CKM Stages",
    x = "CKM Stage",
    y = "eGFR_FAS_highest (mL/min/1.73m²)"
  )

comorb_ckm2$CKM_syn_base

kruskal.test(eGFR_FAS_baseline ~ as.factor(CKM_syn2), data = comorb_ckm2)
dunnTest(eGFR_FAS_baseline ~ as.factor(CKM_syn2), data = comorb_ckm2, method = "none")

kruskal.test(eGFR_FAS_highest ~ as.factor(CKM_syn), data = comorb_ckm2)

kruskal.test(eGFR_FAS_baseline ~ as.factor(CKM_syn_collapsed), data = comorb_ckm2)
dunnTest(eGFR_FAS_baseline ~ as.factor(CKM_syn_collapsed), data = comorb_ckm2, method = "none")

install.packages("ggpubr")
library(ggpubr)

ggboxplot(comorb_ckm2, x = "CKM_syn2", y = "eGFR_FAS_baseline",
          color = "CKM_syn2", palette = "jco",
          add = "jitter", shape = "CKM_syn2") +
  stat_compare_means(method = "kruskal.test", 
                     label.y = max(comorb_ckm2$eGFR_FAS_baseline, na.rm = TRUE) + 10) +
  labs(
    title = "Baseline eGFR (FAS) by CKM Stage",
    x = "CKM Stage (CKM_syn)",
    y = "eGFR_FAS_baseline (mL/min/1.73m²)"
  ) +
  theme_minimal()
install.packages("FSA")  # if not already installed
library(FSA)

## Another variable with eGFR <90 to see the differences by CKM Stages 
## Create the var and perform simple binary analyses 

dunnTest(eGFR_FAS_baseline ~ as.factor(CKM_syn), data = comorb_ckm2, method = "none")

comorb_ckm2$CKM_sy

kruskal.test(HYP ~ CKM_syn_base2, data = comorb_ckm2)

kruskal.test(HYP ~ CKM_syn2, data = comorb_ckm2)
dunnTest(HYP ~ CKM_syn2, data=comorb_ckm2, method = "holm")

kruskal.test(HYP0 ~ CKM_syn_base2, data = comorb_ckm2)
dunnTest(HYP0 ~ CKM_syn_base2, data=comorb_ckm2, method = "holm")


table(comorb_ckm2$HTN0, comorb_ckm2$CKM_syn_base)
table(comorb_ckm2$HTN, comorb_ckm2$CKM_syn)

kruskal.test(HTN0 ~ CKM_syn, data = comorb_ckm2)
dunnTest(HTN0 ~ CKM_syn2, data=comorb_ckm2, method = "holm")

kruskal.test(GLYC ~ CKM_syn, data = comorb_ckm2)

kruskal.test(MIC0 ~ CKM_syn2, data = comorb_ckm2)
dunnTest(MIC0 ~ CKM_syn2, data=comorb_ckm2, method = "holm")

kruskal.test(MIC ~ CKM_syn2, data = comorb_ckm2)
dunnTest(MIC ~ CKM_syn2, data=comorb_ckm2, method = "holm")

kruskal.test(DNE0 ~ CKM_syn2, data = comorb_ckm2) #NOT SIGN
kruskal.test(DNE ~ CKM_syn2, data = comorb_ckm2) #SIGNIFICANT
dunnTest(DNE ~ CKM_syn2, data=comorb_ckm2, method = "holm")

kruskal.test(FILAM ~ CKM_syn, data = comorb_ckm2) #SIGNIGICANT

kruskal.test(NEURO0 ~ CKM_syn2, data = comorb_ckm2)
dunnTest(NEURO0 ~ CKM_syn2, data=comorb_ckm2, method = "holm")
kruskal.test(NEURO ~ CKM_syn2, data = comorb_ckm2)
dunnTest(NEURO ~ CKM_syn, data=comorb_ckm2, method = "holm")

kruskal.test(RETINO ~ CKM_syn2, data = comorb_ckm2)
dunnTest(RETINO ~ CKM_syn2, data=comorb_ckm2, method = "holm")

kruskal.test(MVD0 ~ CKM_syn, data = comorb_ckm2)
dunnTest(MVD0 ~ CKM_syn, data=comorb_ckm2, method = "holm")

kruskal.test(MVD ~ CKM_syn, data = comorb_ckm2)
dunnTest(MVD ~ CKM_syn, data=comorb_ckm2, method = "holm")

install.packages("FSA")
library(FSA)

dunnTest(comorb_ckm2$HTN0 ~ comorb_ckm2$CKM_syn, data = comorb_ckm2, method = "holm") 
dunnTest(GLYC ~ CKM_syn, data=comorb_ckm2, method = "holm")
dunn.test(HYP ~ CKM_syn, data = comorb_ckm2, method = "holm")

comorb_ckm2$CKM_syn_collapsed
kruskal.test(ckd_gfr_highest ~ CKM_syn, data = comorb_ckm2)

kruskal.test(ckd_gfr_highest ~ CKM_syn, data = comorb_ckm2)

kruskal.test(ualb_baseline ~ CKM_syn2, data = comorb_ckm2)
dunnTest(ualb_baseline ~ CKM_syn2, data=comorb_ckm2, method = "holm")

kruskal.test(ualbcreat_baseline ~ CKM_syn2, data = comorb_ckm2)
dunnTest(ualbcreat_baseline ~ CKM_syn2, data=comorb_ckm2, method = "holm")

kruskal.test(il6_highest ~ CKM_syn2, data = comorb_ckm2)
kruskal.test(il6_highest ~ CKM_syn_collapsed, data = comorb_ckm2)
dunnTest(il6_highest ~ CKM_syn2, data=comorb_ckm2, method = "holm")
dunnTest(il6_highest ~ CKM_syn_collapsed, data=comorb_ckm2, method = "holm")

kruskal.test(il1_highest ~ CKM_syn2, data = comorb_ckm2)
kruskal.test(il1_highest ~ CKM_syn_collapsed, data = comorb_ckm2)
dunnTest(il1_highest ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(hscrp_highest ~ CKM_syn2, data = comorb_ckm2)
dunnTest(hscrp_highest ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(apob_highest ~ CKM_syn, data = comorb_ckm2)
dunnTest(apob_highest ~ CKM_syn, data = comorb_ckm2)

kruskal.test(mcp1_highest ~ CKM_syn2, data = comorb_ckm2)
dunnTest(mcp1_highest ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(tnfa_highest ~ CKM_syn2, data = comorb_ckm2)
dunnTest(tnfa_highest ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(tnfr1_highest ~ CKM_syn, data = comorb_ckm2)
dunnTest(tnfr1_highest ~ CKM_syn, data = comorb_ckm2)

kruskal.test(tnfr2_highest ~ CKM_syn, data = comorb_ckm2)
dunnTest(tnfr2_highest ~ CKM_syn, data = comorb_ckm2)

kruskal.test(icam1_highest ~ CKM_syn2, data = comorb_ckm2)
dunnTest(icam1_highest ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(vcam1_highest ~ CKM_syn2, data = comorb_ckm2)
dunnTest(vcam1_highest ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(eselectin_highest ~ CKM_syn2, data = comorb_ckm2)
dunnTest(eselectin_highest ~ CKM_syn, data = comorb_ckm2)

kruskal.test(vegf_highest ~ CKM_syn, data = comorb_ckm2)
dunnTest(vegf_highest ~ CKM_syn, data = comorb_ckm2)

kruskal.test(fgf23_highest ~ CKM_syn, data = comorb_ckm2)
dunnTest(fgf23_highest ~ CKM_syn, data = comorb_ckm2)

comorb_ckm2$bnp_highest
kruskal.test(bnp_baseline ~ CKM_syn2, data = comorb_ckm2)
dunnTest(bnp_baseline ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(bnp_highest ~ CKM_syn2, data = comorb_ckm2)
dunnTest(bnp_highest ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(bmi ~ CKM_syn2, data = comorb_ckm2)
dunnTest(bmi ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(wastcirc ~ CKM_syn2, data = comorb_ckm2)
dunnTest(wastcirc ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(sbp ~ CKM_syn, data = comorb_ckm2)
dunnTest(sbp ~ CKM_syn, data = comorb_ckm2)

kruskal.test(dbp ~ CKM_syn, data = comorb_ckm2)
dunnTest(dbp ~ CKM_syn, data = comorb_ckm2)

comorb_ckm2$hfpef
kruskal.test(pedsql_highest ~ CKM_syn2, data = comorb_ckm2)
dunnTest(pedsql_highest ~ CKM_syn2, data = comorb_ckm2)

kruskal.test(pedsql_highest ~ hfpef1, data = comorb_ckm2)
dunnTest(pedsql_highest ~ as.factor(hfpef2), data = comorb_ckm2)

comorb_ckm2$ckm

library(ggplot2)

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn2), y = ualb_baseline)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
  labs(
    title = "UALB by CKM Stage",
    x = "CKM Stage",
    y = "UALB (mg/dL or other units)"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn2), y = ualbcreat_baseline)) +
  geom_boxplot(fill = "lightgreen", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkgreen") +
  labs(
    title = "UACR by CKM Stage",
    x = "CKM Stage",
    y = "UACR (mg/g)"
  ) +
  theme_minimal()

comorb_ckm2$CKM_syn_collapsed

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn2), y = il6_highest)) +
  geom_boxplot(fill = "lightgreen", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkgreen") +
  labs(
    title = "IL-6 by CKM Stage",
    x = "CKM Stage",
    y = "IL-6 ()"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn2), y = hscrp_highest)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
  labs(
    title = "hs-CRP by CKM Stage",
    x = "CKM Stage",
    y = "hs-CRP ()"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn), y = apob_highest)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
  labs(
    title = "ApoB by CKM Stage",
    x = "CKM Stage",
    y = "ApoB (mg/dL)"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn), y = mcp1_highest)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
  labs(
    title = "MCP-1 by CKM Stage",
    x = "CKM Stage",
    y = "MCP-1 ()"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn2), y = bmi)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
  labs(
    title = "BMI by CKM Stage",
    x = "CKM Stage",
    y = "BMI (kg/m2)"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn2), y = sbp)) +
  geom_boxplot(fill = "pink", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkred") +
  labs(
    title = "SBP by CKM Stage",
    x = "CKM Stage",
    y = "SBP (mmHg)"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn), y = bnp_baseline)) +
  geom_boxplot(fill = "orange", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkorange3") +
  labs(
    title = "BNP-baseline by CKM Stage",
    x = "CKM Stage",
    y = "BNP (ng/mL)"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn), y = bnp_highest)) +
  geom_boxplot(fill = "orange", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkorange3") +
  labs(
    title = "BNP-highest by CKM Stage",
    x = "CKM Stage",
    y = "BNP (ng/mL)"
  ) +
  theme_minimal()

ggplot(comorb_ckm2, aes(x = as.factor(CKM_syn2), y = pedsql_highest)) +
  geom_boxplot(fill = "orange", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkorange3") +
  labs(
    title = "PedsQL by CKM Stage",
    x = "CKM Stage",
    y = "PedsQL"
  ) +
  theme_minimal()

hist(comorb_ckm2$bnp_highest)
hist(baserisk$c)

install.packages("ggalluvial")
library(ggalluvial)
library(ggplot2)

df_alluvial <- comorb_ckm2 %>%
  select(RELEASEID, CKM_syn_base, CKM_syn) %>%
  mutate(
    CKM_syn_base = as.factor(CKM_syn_base),
    CKM_syn = as.factor(CKM_syn)
  )

ggplot(df_alluvial,
       aes(axis1 = CKM_syn_base, axis2 = CKM_syn)) +
  geom_alluvium(aes(fill = CKM_syn_base), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Baseline CKM Stage", "Follow-up CKM Stage"), expand = c(.05, .05)) +
  labs(
    title = "Transitions Between CKM Stages",
    y = "Number of Participants"
  ) +
  theme_minimal()

## second Sankey if the patients are grouped according to the BMI

df_alluvial2 <- comorb_ckm2 %>%
  select(RELEASEID, CKM_syn_base2, CKM_syn2) %>%
  mutate(
    CKM_syn_base2 = as.factor(CKM_syn_base2),
    CKM_syn2 = as.factor(CKM_syn2)
  )

ggplot(df_alluvial2,
       aes(axis1 = CKM_syn_base2, axis2 = CKM_syn2)) +
  geom_alluvium(aes(fill = CKM_syn_base2), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Baseline CKM Stage", "Follow-up CKM Stage"), expand = c(.05, .05)) +
  labs(
    title = "Transitions Between CKM Stages",
    y = "Number of Participants"
  ) +
  theme_minimal()

comorb_ckm2$ckm

comorb_ckm2 %>%
  filter(CKM_syn == "Stage2") %>%
  summarise(
    BMI_ge_30 = sum(bmi >= 30, na.rm = TRUE),
    BMI_lt_30 = sum(bmi < 30, na.rm = TRUE)
  )

sarah <- comorb_ckm2 %>%
  select(RELEASEID, sex, AGEBASE, testosterone_baseline, testosterone_highest, testosterone_latest,
         edol_baseline, edol_highest, edol_latest)

glucose_ins_today2 <- cbl_2_unique %>%
  arrange(releaseid, days) %>%
  group_by(releaseid) %>%
  filter(!is.na(glucose) | !is.na(ins)) %>%  # keep rows with at least one value present
  slice(1) %>%
  ungroup() %>%
  select(releaseid, glucose, ins, pvisit, days)

glucose_ins_first <- cbl_1_unique %>%
  arrange(releaseid, days) %>%
  group_by(releaseid) %>%
  filter(!is.na(Glucose) | !is.na(ins)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(RELEASEID = toupper(releaseid)) %>%   # create uppercase RELEASEID
  select(RELEASEID, Glucose, ins, mvisit, days)

glucose_ins_today2 <- cbl_2_unique %>%
  arrange(releaseid, days) %>%
  group_by(releaseid) %>%
  filter(!is.na(glucose) | !is.na(ins)) %>% 
  slice(1) %>%
  ungroup() %>%
  mutate(RELEASEID = toupper(releaseid)) %>%
  select(RELEASEID, glucose2 = glucose, ins2 = ins, pvisit2 = pvisit, days2 = days)

sarah <- sarah %>%
  left_join(glucose_ins_today2, by = "RELEASEID")

sarah <- sarah %>%
  left_join(glucose_ins_first, by = "RELEASEID")

getwd()
write.csv(sarah, file = "sarah.csv", row.names = FALSE)

```

## Results

### Count of events (can be multiple events per person)

```{r results="asis"}
event_table <-tableby(data = ame, ~ as.factor(AMENAME))
print(summary(event_table))
```

### Count of people with each type of event

```{r results="asis"}
person_table <- tableby(data = ame_summary, ~ ARRHYTHMIA + CAD + CHF + LVSD + MI + PAD + DVT + STROKE + TIA + PANCREATITIS + GALLBLADDER + PNEURO + ANEURO + MONONEURO + CKD + ESKD + NPDR + PDR + ME + VH + BLINDNESS + CATARACTS + GLAUCOMA + DEATH)
print(summary(person_table))
```


