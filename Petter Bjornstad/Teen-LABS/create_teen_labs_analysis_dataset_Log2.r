library(tidyverse)
library(sas7bdat)
library(readxl)
library(SomaDataIO)
library(Hmisc)
library(sas7bdat)
library(dplyr)

#Directories
user <- "hhampson" #Replace with username
dir.dat <- c(paste0("/Users/",user,"/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Teen Labs"))

# Clinical Data
clinical <- read.sas7bdat(fs::path(dir.dat,"/Data_Cleaned/bjornstad_03_15_2023.sas7bdat"))

# Proteomics data
soma <- read_adat(fs::path(dir.dat,"/Data_Raw/WUS_22_007_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"))
soma <- soma %>%
  # Remove flagged rows, QC samples, and excluded samples
  filter(
    RowCheck == "PASS", SampleType == "Sample",
    !SampleDescription %in% c("FAF001634 0200", "FAF001855 0200")
  ) %>%
  # Remove unnecessary columns
  select(SampleDescription, contains("seq."))
# Log 2 transform
soma[,grep("seq",colnames(soma))] = lapply(soma[,grep("seq",colnames(soma))],log2)

# analytes
analytes <- getAnalyteInfo(soma)
analytes <- analytes %>% select(AptName,SeqId,SeqIdVersion,SomaId,TargetFullName,Target,UniProt,EntrezGeneID,EntrezGeneSymbol,Organism,Units,Type)
# remove fc mouse and no protein
drop <- analytes %>% filter(Target == "Fc_MOUSE" | Target == "No Protein" | !(Organism == "Human") | !(Type == "Protein"))
apt_drop <- drop$AptName
soma <- soma %>% select(!all_of(apt_drop))
analytes <- analytes %>% filter(!Target == "Fc_MOUSE")
analytes <- analytes %>% filter(!Target == "No Protein")
analytes <- analytes %>% filter(Organism == "Human")
analytes <- analytes %>% filter(Type == "Protein")

# Merge
df <- left_join(clinical, soma, by = c("SAMPLE_ID" = "SampleDescription"))
# Create new variables
df <- df %>%
  mutate(albuminuria = case_when(
    UACRATIO * 1000 < 30 ~ "A1",
    UACRATIO * 1000 >= 30 & UACRATIO * 1000 <= 300 ~ "A2",
    UACRATIO * 1000 > 300 ~ "A3"
  )) %>%
  dplyr::group_by(ID) %>%
  mutate(diab_resolved = case_when(
    sum(na.omit(diab)) > 0 & (last(na.omit(diab)) - first(na.omit(diab))) < 0 ~ 1,
    sum(na.omit(diab)) > 0 & (last(na.omit(diab)) - first(na.omit(diab))) == 0 ~ 0,
    sum(na.omit(diab)) == 0 ~ 2
  )) %>%
  mutate(race_ethnicity = case_when(
    RACE == 1 & ETHN == 2 ~ "Non-Hispanic White",
    RACE == 2 & ETHN == 2 ~ "Non-Hispanic Black",
    ETHN == 1 ~ "Hispanic",
    T ~ "Other"
  )) %>%
  group_by(ID) %>% mutate(diab_baseline = first(na.omit(diab))) %>% ungroup()

# Categorical data
df$SEX <- factor(df$SEX, levels = 1:2, labels = c("Male", "Female"))
df$ETHN <- factor(df$ETHN,
                  levels = 1:2,
                  labels = c("Hispanic", "Non-Hispanic")
)
df$RACE <- factor(df$RACE,
                  levels = 1:8,
                  labels = c(
                    "White or Caucasian",
                    "Black or African-American",
                    "Asian",
                    "American Indian or Alaska Native",
                    "Native Hawaiian or other Pacific Islander",
                    "Other",
                    "Unknown",
                    "More than one race"
                  )
)
df$SURG <- factor(df$SURG,
                  levels = c(1, 4, 5),
                  labels = c(
                    "Gastric bypass",
                    "Laparoscopic adjustable gastric band",
                    "Sleeve gastrectomy - initial stage"
                  )
)
df$POTENTIALPREG <- factor(df$POTENTIALPREG,
                           levels = 0:1,
                           labels = c("No", "Yes")
)
df$fasting <- factor(df$fasting,
                     levels = 0:1,
                     labels = c("Non-Fasting", "Fasting")
)
df$labid <- factor(df$labid, levels = 1:2, labels = c("NWRL", "Quest"))
df$diab <- factor(df$diab, levels = 0:1, labels = c("No", "Yes"))
df$diab_baseline <- factor(df$diab_baseline, levels = 0:1, labels = c("No", "Yes"))
df$months <- df$visit
df$visit <- factor(df$visit,
                   levels = c(1, 6, 12, 24, 36, 48, 60),
                   labels = c(
                     "Month 1", "Month 6", "Year 1", "Year 2",
                     "Year 3", "Year 4", "Year 5"
                   )
)
df$diab_resolved <- factor(df$diab_resolved,
                           levels = 0:2,
                           labels = c("No", "Yes", "Non-diabetic")
)

# calculate QCR values for FAS equations
df$age <- floor(df$age)
df$qcr[df$age==8] <- 0.46
df$qcr[df$age==9] <- 0.49
df$qcr[df$age==10] <- 0.51 
df$qcr[df$age==11] <- 0.53 
df$qcr[df$age==12] <- 0.57
df$qcr[df$age==13] <- 0.59
df$qcr[df$age==14] <- 0.61
# females
df$qcr[df$age==15 & df$SEX=="Female"] <- 0.64
# males
df$qcr[df$age==15 & df$SEX=="Male"] <- 0.72
# females
df$qcr[df$age==16 & df$SEX=="Female"] <- 0.67
# males
df$qcr[df$age==16 & df$SEX=="Male"] <- 0.78
# females
df$qcr[df$age==17 & df$SEX=="Female"] <- 0.69
# males
df$qcr[df$age==17 & df$SEX=="Male"] <- 0.82
# females
df$qcr[df$age==18 & df$SEX=="Female"] <- 0.69
# males
df$qcr[df$age==18 & df$SEX=="Male"] <- 0.85
# females
df$qcr[df$age==19 & df$SEX=="Female"] <- 0.70
# males
df$qcr[df$age==19 & df$SEX=="Male"] <- 0.88
# females
df$qcr[df$age>19 & df$SEX=="Female"] <- 0.70
# males
df$qcr[df$age>19 & df$SEX=="Male"] <- 0.90

# eGFR FAS creatinine
df$eGFR.fas_cr <-107.3/(df$CREAS/df$qcr)

# eGFR FAS combined creatinine and cystatin-C
df$f1 <- df$CREAS/df$qcr
df$f2 <- 1-0.5
df$f3 <- df$CYSC/0.82
df$eGFR.fas_cr_cysc <- 107.3 / ((0.5*df$f1) + (df$f2*df$f3))
label(df$eGFR.fas_cr_cysc) <- "eGFR FAS Cr Cys-C"
label(df$eGFR.fas_cr) <- "eGFR FAS Cr"

# read in NAFLD data
nafld <- read.csv(fs::path(dir.dat,"/Data_Cleaned/Pyle NAFLD 02_06_24.csv"))
nafld <- unique(nafld)
nafld <- nafld %>%
  mutate(visit = case_when(
    visit == 1 ~ "Month 1",
    visit == 6 ~ "Month 6",
    visit == 12 ~ "Year 1",
    visit == 36 ~ "Year 3"
  )) 
df <- left_join(df, nafld, by = c("ID", "visit"))

# read in diabetes data
diabetes <- read.csv(fs::path(dir.dat,"/Data_Cleaned/Teen-LABS_DM.csv"))
diabetes <- unique(diabetes)
diabetes$visit <- diabetes$Visit
diabetes$Visit <- NULL
diabetes <- diabetes %>%
  mutate(visit = case_when(
    visit == 1 ~ "Month 1",
    visit == 6 ~ "Month 6",
    visit == 12 ~ "Year 1",
    visit == 36 ~ "Year 3"
  )) 
diabetes <- diabetes %>%
  mutate(diabetes_duration = case_when(
    DBTOTAL == -2 ~ NA,
    DBTOTAL == -5 ~ NA,
    DBTOTAL == -10 ~ NA,
    .default = DBTOTAL
  )) 
# drop DIABETES variable as diab has missing values completed, also drop HBA1C
diabetes <- diabetes %>% select("ID", "visit", "DBTOTAL", "diabetes_duration")
df <- left_join(df, diabetes, by = c("ID", "visit"))

# read in ALT/AST data - "-5" means not done
alt_ast <- read.csv(fs::path(dir.dat,"/Data_Cleaned/Teen-LABS_AST_ALT.csv"), na.strings = c("-5"))
# "-7" means below lower limit of detection. Waiting to hear from Hailey what that is....for now, I will just put in a 1
alt_ast$AST <- ifelse(alt_ast$AST == "-7", 1, alt_ast$AST)
alt_ast$ALT <- ifelse(alt_ast$ALT == "-7", 1, alt_ast$ALT)
alt_ast$AST <- as.numeric(alt_ast$AST)
alt_ast$ALT <- as.numeric(alt_ast$ALT)
alt_ast$visit <- alt_ast$Visit
alt_ast$Visit <- NULL
alt_ast <- alt_ast %>%
  mutate(visit = case_when(
    visit == 1 ~ "Month 1",
    visit == 6 ~ "Month 6",
    visit == 12 ~ "Year 1",
    visit == 24 ~ "Year 2",
    visit == 36 ~ "Year 3",
    visit == 48 ~ "Year 4",
    visit == 60 ~ "Year 5",
    visit == 72 ~ "Year 6",
    visit == 96 ~ "Year 8"
  )) 
# getting warning about NAs introduced by coercion but I don't see where that is an issue
df <- left_join(df, alt_ast, by = c("ID", "visit"))

# read in Olink data
olink <- read.csv(fs::path(dir.dat,"/Olink data/TL_baseline_proteomics_processed_wide.csv"))
linking_file <- read.csv(fs::path(dir.dat,"/Olink data/Pyle 01_19_24.csv"))
linking_file <- linking_file %>% filter(visit == 1)
# drop NIH ID because we don't have that information in the Olink file
linking_file <- linking_file %>% select(key, visit, TLID)
linking_file <- unique(linking_file)
olink <- left_join(olink, linking_file, by = "key")
olink$ID <- olink$TLID
olink$visit <- ifelse(olink$visit == 1, "Month 1", olink$visit)
df <- left_join(df, olink, by = c("ID", "visit"))

# Labels
dict <- read_excel(fs::path(dir.dat,"/Data_Cleaned/DataDictionary.xlsx"))
analytes <- getAnalytes(soma)
analyte_info <- getAnalyteInfo(soma)
labels <- dict$Description[match(names(df), dict$Name)]
names(labels) <- colnames(df)
labels[analytes] <- paste0("log(",analyte_info$EntrezGeneSymbol[match(names(labels[analytes]), analyte_info$AptName)],")")
label(df) <- as.list(labels)
label(df$diab) <- "Diabetes"
label(df$diab_resolved) <- "Diabetes resolved?"
label(df$albuminuria) <- "Albuminuria level"
label(df$months) <- "Months"
label(df$diab_baseline) = "Diabetes at Baseline"

# As regular dataframe
df <- as.data.frame(df)
# Save
save(df,analyte_info, file = fs::path(dir.dat,"/Data_Cleaned/analysis_dataset_log2.RData"))
# write list of IDs
write.csv(unique(df$ID), file = fs::path(dir.dat,"/Data_Cleaned/id_list.csv"),
          row.names = F)
