library(tidyverse)
library(sas7bdat)
library(readxl)
library(SomaDataIO)
library(Hmisc)
library(sas7bdat)
library(dplyr)
# Clinical Data
clinical <- read.sas7bdat("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/Data_Cleaned/bjornstad_03_15_2023.sas7bdat")
# Proteomics data
soma <- read_adat("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/Data_Raw/WUS_22_007_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
soma <- soma %>%
  # Remove flagged rows, QC samples, and excluded samples
  filter(
    RowCheck == "PASS", SampleType == "Sample",
    !SampleDescription %in% c("FAF001634 0200", "FAF001855 0200")
  ) %>%
  # Remove unnecessary columns
  select(SampleDescription, contains("seq."))
# Log transform
soma[,grep("seq",colnames(soma))] = lapply(soma[,grep("seq",colnames(soma))],log)
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

# read in NAFLD data
nafld <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/Data_Cleaned/Pyle NAFLD 02_06_24.csv")
nafld <- unique(nafld)
nafld <- nafld %>%
  mutate(visit = case_when(
    visit == 1 ~ "Month 1",
    visit == 6 ~ "Month 6",
    visit == 12 ~ "Year 1",
    visit == 36 ~ "Year 3"
  )) 
df <- left_join(df, nafld, by = c("ID", "visit"))

# read in Olink data
olink <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/Olink data/TL_baseline_proteomics_processed_wide.csv")
linking_file <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/Olink data/Pyle 01_19_24.csv")
linking_file <- linking_file %>% filter(visit == 1)
# drop NIH ID because we don't have that information in the Olink file
linking_file <- linking_file %>% select(key, visit, TLID)
linking_file <- unique(linking_file)
olink <- left_join(olink, linking_file, by = "key")
olink$ID <- olink$TLID
olink$visit <- ifelse(olink$visit == 1, "Month 1", olink$visit)
df <- left_join(df, olink, by = c("ID", "visit"))

# Labels
dict <- read_excel("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/Data_Cleaned/DataDictionary.xlsx")
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
save(df,analyte_info, file = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/Data_Cleaned/analysis_dataset.RData")
