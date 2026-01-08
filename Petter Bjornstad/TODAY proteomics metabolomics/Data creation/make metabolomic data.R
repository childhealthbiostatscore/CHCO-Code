library(dplyr)
library(berryFunctions)
library(stringr)
library(openxlsx)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward"
}

home_dir = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/TODAY subaward"

setwd(home_dir)

####################
# normalized urine #
####################
nih_urine <- openxlsx::read.xlsx("./Metabolomic data/NIDDK_AA_20220427_20220620_nM_AF.xlsx", sheet = "Urine",
                                 startRow = 2,colNames = TRUE)

# add OA data
nih_urine_oa <- openxlsx::read.xlsx("./Metabolomic data/LEAD_NIDDK_OA_08292022.xlsx", sheet = "NIDDK Urine",
                                        startRow = 2,colNames = TRUE)
xtra <- openxlsx::read.xlsx("./Metabolomic data/6217497 Urine.xlsx",startRow = 2,colNames = TRUE)
nih_urine_oa <- rbind(nih_urine_oa,xtra)
# ID scheme (sample name and freezerworks ID) appears to be switched from the previous file
# sent email to Anthony and Kumar to confirm
nih_urine_oa$Freezerworks.ID <- nih_urine_oa$Sample.Name
nih_urine_oa$Sample.Name <- NULL
nih_urine <- merge(nih_urine,nih_urine_oa,by="Freezerworks.ID",all.x = T,all.y = T)
nih_urine$site <- "NIH"

# read in LEAD samples
lead_urine <- openxlsx::read.xlsx("./Metabolomic data/Lead_AA_20220322_20220620_nM_AF.xlsx", sheet = "Urine",
                                 startRow = 2,colNames = TRUE)
# add OA data
lead_urine_oa <- openxlsx::read.xlsx("./Metabolomic data/LEAD_NIDDK_OA_08292022.xlsx", sheet = "LEAD Urine",
                                    startRow = 2,colNames = TRUE)
lead_urine_oa$Sample.Name <- NULL
# one sample shows up in both the LEAD urine and plasma tabs of the OA, remove from urine
lead_urine_oa <- lead_urine_oa %>% filter(!Freezerworks.ID==165196)
lead_urine <- merge(lead_urine, lead_urine_oa, by="Freezerworks.ID",all.x = T,all.y = T)
lead_urine$site <- "LEAD"

####################
# plasma           #
####################

# read in NIH samples
nih_plasma <- openxlsx::read.xlsx("./Metabolomic data/NIDDK_AA_20220427_20220620_nM_AF.xlsx", sheet = "Plasma",
                                  startRow = 2,colNames = TRUE)
# add OA data
nih_plasma_oa <- openxlsx::read.xlsx("./Metabolomic data/LEAD_NIDDK_OA_08292022.xlsx", sheet = "NIDDK Plasma",
                                     startRow = 2,colNames = TRUE)
nih_plasma_oa$Freezerworks.ID <- nih_plasma_oa$Sample.Name
nih_plasma_oa$Sample.Name <- NULL
nih_plasma_oa$Creatinine.in.uM <- NULL
nih_plasma <- merge(nih_plasma,nih_plasma_oa,by="Freezerworks.ID",all.x = T,all.y = T)
nih_plasma$site <- "NIH"

# read in LEAD samples
lead_plasma <- openxlsx::read.xlsx("./Metabolomic data/Lead_AA_20220322_20220620_nM_AF.xlsx", sheet = "Plasma",
                                  startRow = 2,colNames = TRUE)
# add OA data
lead_plamsa_oa <- openxlsx::read.xlsx("./Metabolomic data/LEAD_NIDDK_OA_08292022.xlsx", sheet = "LEAD Plasma",
                                      startRow = 2,colNames = TRUE)
lead_plamsa_oa$Sample.Name <- NULL
lead_plamsa_oa$Creatinine.in.uM <- NULL
lead_plasma <- merge(lead_plasma,lead_plamsa_oa,by="Freezerworks.ID")
lead_plasma$site <- "LEAD"

# read in LEAD OA samples (2025)
lead_plasma_oa_2025 <- openxlsx::read.xlsx("./Metabolomic data/200251008_TODAY_LEAD_Organic_Acid Panel Metabolites_Final.xlsx", sheet = "TODAY_LEAD_OA",
                                        startRow = 2,colNames = TRUE)
lead_plasma_oa_2025$Sample.Name <- gsub("_Plasma", "", lead_plasma_oa_2025$File.Name)
lead_plasma <- merge(lead_plasma,lead_plasma_oa_2025,by="Sample.Name")


######################
# link IDs and merge #
######################

# read in the files that will link repository ID to sample ID
ids_lead <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at Colorado LEAD Center - UT Health San Antonio.csv")
ids_lead$bsi_id <- NA
ids_niddk_today <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY - UT Health San Antonio.csv")
ids_niddk_today$MASK.ID <- NA
ids_niddk_today2 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY2 - UT Health San Antonio.csv")
ids_niddk_today2$MASK.ID <- NA
ids_niddk <- rbind(ids_niddk_today, ids_niddk_today2)

# merge to urine
# first NIH - merge Freezerworks ID in results file to current label in NIDDK file
nih_urine$current_label <- nih_urine$Freezerworks.ID
nih_urine <- merge(nih_urine,ids_niddk,by="current_label",all.x = T, all.y = F)
nih_urine$SAMPLE_ID <- NA
# for LEAD, need to merge MASK.ID in ID file to sample name in results file
ids_lead$Sample.Name <- str_sub(ids_lead$MASK.ID, 1, 8)
ids_lead$Sample.Name <- str_replace(ids_lead$Sample.Name, "-", "_")
ids_lead$t <- ifelse(str_trim(ids_lead$material_type)=="Urine","_U","")
ids_lead$Sample.Name <- paste0(ids_lead$Sample.Name,ids_lead$t)
ids_lead$t <- NULL
lead_urine <- merge(lead_urine,ids_lead,by="Sample.Name",all.x=T, all.y = F)
lead_urine$current_label <- NA
# combine NIH and LEAD
urine <- rbind(nih_urine,lead_urine)

# merge to plasma
# first NIH
nih_plasma$current_label <- nih_plasma$Freezerworks.ID
nih_plasma <- merge(nih_plasma,ids_niddk,by="current_label",all.x = T, all.y = F)
nih_plasma$SAMPLE_ID <- NA
# LEAD
lead_plasma <- merge(lead_plasma,ids_lead,by="Sample.Name",all.x=T, all.y = F)
lead_plasma$current_label <- NA
# combine NIH and LEAD
# plasma <- rbind(nih_plasma,lead_plasma)
nih_plasma$Sample.Name <- as.character(nih_plasma$Sample.Name)
plasma <- bind_rows(nih_plasma,lead_plasma)

# Save
# save(urine,file = "./Metabolomic data/urine.Rdata")
# save(plasma,file = "./Metabolomic data/plasma.Rdata")
write.csv(plasma, file = "./Metabolomic data/today_plasma_metabolomics.csv", row.names = F, na = "")
