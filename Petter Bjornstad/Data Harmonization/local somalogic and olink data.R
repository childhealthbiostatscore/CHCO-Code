library(SomaDataIO)
library(stringr)
library(arsenal)
library(Hmisc)
library(dplyr)
library(tidyr)
library(purrr)

# if(Sys.info()["sysname"] == "Windows"){
#   home_dir = "E:/Petter Bjornstad"
# } else if (Sys.info()["sysname"] == "Linux"){
#   home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad"
# } else if (Sys.info()["sysname"] == "Darwin"){
#   home_dir = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/"
# }
# setwd(home_dir)

#setwd("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/")
setwd("/Users/laurapyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive")

###################
# ANML ADAT FILES #
###################

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("./Local cohort Somalogic data/WUS-22-002/WUS-22-002_v4.1_EDTAPlasma_hybNorm_medNormInt_plateScale_calibrate_anmlQC_qcCheck_anmlSMP.adat")
soma <- soma %>% select(-Optional2)
soma <- soma %>% filter(!is.na(SampleDescription))
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

# read in 2nd dataset
soma2 <- read_adat("./Local cohort Somalogic data/WUS-23-004/WUS_23_004_2023-11-15/WUS_23_004_v4.1_EDTAPlasma.hybNorm_medNormInt_plateScale_calibrate_anmlQC_qcCheck_anmlSMP.adat")
soma2 <- soma2 %>% filter(!is.na(SampleDescription))
# delete extraneous stuff in the sample description
soma2$SampleDescription <- str_sub(soma2$SampleDescription, 12)
soma2$SampleDescription <- str_remove(soma2$SampleDescription, " ")
soma2$SampleDescription <- str_remove_all(soma2$SampleDescription, "\\(\\S+")
# there is one duplicate sample - delete the second result
soma2 <- soma2 %>% filter(!SampleDescription=="RH-14-O")
# remove FC_Mouse and no protein
soma2 <- soma2 %>% select(!all_of(apt_drop))
analytes2 <- getAnalyteInfo(soma2)
# 2nd analytes file is esssentially the same as the first except for some batch specific information we don't need
# will keep the first file

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("./Local cohort Somalogic data/WUS-22-002/WUS-22-002_v4.1_EDTAPlasma_hybNorm_medNormInt_plateScale_calibrate_anmlQC_qcCheck_anmlSMP.adat")
soma <- soma %>% select(-Optional2)
soma <- soma %>% filter(!is.na(SampleDescription))
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

# read in PANTHER data
panther <- read_adat("./Local cohort Somalogic data/PANTHER/20240126_597_Bjornstad_SOMAscan7k_WUS-24-002_data_export/WUS-24-002_2024-01-26_Somalogic_standardized_files/WUS_24_002_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
panther <- panther %>% filter(!is.na(SampleDescription))
panther <- panther %>% select(-Optional1)
# remove fc mouse and no protein
panther <- panther %>% select(!all_of(apt_drop))
# one ID changed
panther$SampleDescription <- ifelse(panther$SampleDescription == "PAN-14-C", "PAN-14-O", panther$SampleDescription)

# read in ATTEMPT data
attempt <- read_adat("./Local cohort Somalogic data/ATTEMPT/20250313_SR006998_Bjornstad_SOMAscan_11K_WUS_25_003_data_export/WUS_25_003_somalogic_standardized_data_files/WUS_25_003_v5.0_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.20250312.adat")
attempt <- attempt %>% filter(!is.na(SampleDescription))
attempt <- attempt %>% select(-Optional1, -Optional2, -SampleNumber)
# remove fc mouse and no protein
attempt <- attempt %>% select(!all_of(apt_drop))

# ATTEMPT analytes file (9k)
analytes_attempt <- getAnalyteInfo(attempt)
analytes_attempt <- analytes_attempt %>% select(AptName,SeqId,SeqIdVersion,SomaId,TargetFullName,Target,UniProt,EntrezGeneID,EntrezGeneSymbol,Organism,Units,Type)
# remove fc mouse and no protein
drop <- analytes_attempt %>% filter(Target == "Fc_MOUSE" | Target == "No Protein" | !(Organism == "Human") | !(Type == "Protein"))
apt_drop <- drop$AptName
attempt <- attempt %>% select(!all_of(apt_drop))
analytes_attempt <- analytes_attempt %>% filter(!Target == "Fc_MOUSE")
analytes_attempt <- analytes_attempt %>% filter(!Target == "No Protein")
analytes_attempt <- analytes_attempt %>% filter(Organism == "Human")
analytes_attempt <- analytes_attempt %>% filter(Type == "Protein")
additional_proteins9k <- colnames(attempt)[colnames(attempt) %nin% colnames(soma)]
# add in ATTEMPT urine proteomics
attempt_urine <- read_adat("./Local cohort Somalogic data/ATTEMPT/20250527_Petter_SOMAscan7k_Urine_WUS_25_015_SR007698_data_export/WUS_25_015_v4.1_Urine.hybNorm.medNormInt.plateScale.anmlQC.qcCheck.AddLOD.anmlSMP.adat")
attempt_urine <- attempt_urine %>% filter(!is.na(SampleDescription))
attempt_urine <- attempt_urine %>% select(-SampleNumber)
urine_analytes <- getAnalyteInfo(attempt_urine)
# remove fc mouse and no protein
attempt_urine <- attempt_urine %>% select(!all_of(apt_drop))
#attempt_urine[additional_proteins9k] <- NA
attempt_urine[,c("NormScale_20","NormScale_0_005","NormScale_0_5","ANMLFractionUsed_20","ANMLFractionUsed_0_005","ANMLFractionUsed_0_5")] <- NA
drop <- urine_analytes %>% filter(Target == "Fc_MOUSE" | Target == "No Protein" | !(Organism == "Human") | !(Type == "Protein"))
apt_drop <- drop$AptName
attempt_urine <- attempt_urine %>% select(!all_of(apt_drop))
urine_analytes <- urine_analytes %>% filter(!Target == "Fc_MOUSE")
urine_analytes <- urine_analytes %>% filter(!Target == "No Protein")
urine_analytes <- urine_analytes %>% filter(Organism == "Human")
urine_analytes <- urine_analytes %>% filter(Type == "Protein")
# add _urine to colnames in the urine proteomics dataset
new_names <- paste0(colnames(attempt_urine), "_urine")
colnames(attempt_urine) <- new_names
# recode one ID
attempt_urine$SampleDescription_urine <- ifelse(attempt_urine$SampleDescription_urine == "a10397", "10397", attempt_urine$SampleDescription_urine)
# merge ATTEMPT blood and urine proteomics by SampleDescription and TimePoint (need to remove _urine from urine variable names)
attempt_urine$SampleDescription <- as.numeric(attempt_urine$SampleDescription_urine)
attempt_urine$TimePoint <- attempt_urine$TimePoint_urine
rownames(attempt) <- NULL
rownames(attempt_urine) <- NULL
additional_proteins7k <- colnames(attempt_urine)[colnames(attempt_urine) %nin% colnames(soma)]
attempt <- full_join(attempt, attempt_urine, by = c("SampleDescription", "TimePoint"))
attempt$SampleDescription <- as.character(attempt$SampleDescription)

#soma[additional_proteins9k] <- NA
#soma2[additional_proteins9k] <- NA
#panther[additional_proteins9k] <- NA
# combine blood proteomics dfs
#soma <- rbind(soma,soma2,panther,attempt)
rownames(soma) <- NULL
rownames(soma2) <- NULL
soma <- rbind(soma, soma2)
rownames(panther) <- NULL
soma <- rbind(soma, panther)
#soma$VolumeSubmitted <- as.numeric(soma$VolumeSubmitted)
# BELOW IS NOT WORKING
soma <- bind_rows(soma, attempt)
# for some reason, colnames with same names are not getting merged
# delete Pima data
soma <- soma %>% filter(!str_detect(SampleDescription,"CKDS"))
# fix sample IDs on a few RH2 participants who changed groups
soma <- soma %>%  mutate(
  SampleDescription = case_when(
    SampleDescription == "RH2-13-T" ~ "RH2-13-O",
    SampleDescription == "RH2-06-T" ~ "RH2-06-O",
    SampleDescription == "RH2-07-T" ~ "RH2-07-O",
    SampleDescription == "RH2-38-O" ~ "RH2-38-T",
    .default = SampleDescription
  )
)
soma[,c("NormScale_75_S1","NormScale_75_S3","NormScale_75_S2","ANMLFractionUsed_75_S1","ANMLFractionUsed_75_S3","ANMLFractionUsed_75_S2")] <- NA

# read in harmonized dataset to get group information
# df <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", 
#                na.strings = c(" ", "", "-9999",-9999))
# coenroll_id <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Renal HERITAGE/Data_Cleaned/coenrolled_ids.csv") %>%
#   pivot_longer(cols = 'improve_id':'crc_id',
#                values_to = "record_id") %>% 
#   dplyr::select(merged_id, record_id) %>%
#   filter(record_id != "")
# df <- df %>%
#   dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
#                    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
#                    .by = c(record_id, visit)) %>%
#   filter(participation_status!="Removed"|is.na(participation_status)) %>%
#   left_join(coenroll_id)
# df <- df %>% 
#   dplyr::select(record_id, co_enroll_id, merged_id, study, visit, date, group)
# # Merge
# soma_harmonized$record_id <- soma_harmonized$SampleDescription
# soma_harmonized <- soma_harmonized %>% mutate(
#   record_id = case_when(
#     str_detect(record_id, "IT2D") ~ str_replace(record_id, "IT2D-", "IT_"),
#     !str_detect(record_id, "IT2D") ~ record_id
#   )
# )
# soma_harmonized <- soma_harmonized %>%  mutate(
#   visit = case_when(
#     TimePoint == "BL" ~ "baseline",
#     TimePoint == "Baseline" ~ "baseline",
#     TimePoint == "3M" ~ "3_months_post_surgery",
#     TimePoint == "12M" ~ "12_months_post_surgery"
#   )
# )
# soma_harmonized <- left_join(soma_harmonized, df, by = c("record_id", "visit"))

# some extra spaces in the sample IDs
soma_combined$SampleDescription <- str_trim(soma_combined$SampleDescription)
# add labels
labels <- paste0("log(",analytes_attempt$EntrezGeneSymbol[match(colnames(soma_combined), analytes_attempt$AptName)],")")
labels[labels == "log(NA)"] <- ""
label(soma_combined) <- as.list(labels)

# write copy of the entire local SomaScan dataframe to the data harmonization folder
save(soma_combined, file = "./Data Harmonization/Combined SomaScan/soma_combined_anml_2.RData")
save(analytes_attempt, file = "./Data Harmonization/Combined SomaScan/analytes_2.RData")
