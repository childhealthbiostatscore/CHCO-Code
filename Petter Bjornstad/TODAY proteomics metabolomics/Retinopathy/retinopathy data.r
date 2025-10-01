library(SomaDataIO)
library(limma)
library(dplyr)
library(caret)
library(purrr)
library(multtest)
library(openxlsx)
library(tableone)
library(EnhancedVolcano)
library(knitr)
library(survival)
library(broom)
library(emmeans)
library(ggvenn)

knitr::opts_chunk$set(echo = FALSE,warning = FALSE)

knitr::opts_chunk$set(echo = FALSE)
if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward"
}
knitr::opts_knit$set(root.dir = home_dir)

setwd(home_dir)
# load somalogic data, with QC samples already excluded
load("./Somalogic data raw/soma.Rdata")

# load analyte info
load("./Somalogic data raw/analytes.Rdata")

# load comorbidity data
load("./Clinical data/comorb.Rdata")

# load baseline risk factors
load("./Clinical data/TODAY/baserisk.Rdata")

# take only the baseline soma samples
# can we just take the earliest or is it possible someone would have a follow-up sample but not baseline?
# probably can take the first and then check the years to make sure compatible with TODAY
base <- soma %>% arrange(releaseid,Date.Drawn) %>% dplyr::group_by(releaseid) %>% dplyr::filter(row_number()==1)
# these 3 release IDs have first sample after end of recruitment, so they must be missing the baseline visit
base <- base %>% filter(!releaseid %in% c("65-85903","65-47984","65-25901"))

# the above section of code used to work but seems to have broken....this should accomplish the same thing?
#base <- soma %>% filter(visit==1)

# merge in complication data
base <- left_join(base, comorb, by="releaseid")
# this was previously:
# base <- merge(base, comorb, by="releaseid",all.x=T, all.y=F)

# merge in baseline risk factors
base <- left_join(base, baserisk, by="releaseid")
# this was previously:
#base <- merge(base, baserisk, by="releaseid",all.x=T, all.y=F)

# log transform UAlbCreat
base$log_UAlbCreat <- log(base$UAlbCreat + 0.0000001)

# identify columns corresponding to proteins
#is_seq <- function(.x) grepl("^seq\\.[0-9]{4}", .x) # regex for analytes
is_seq <- function(.x) grepl("seq", .x)
seq <- is_seq(names(base))

# convert to numeric
base[,seq] <- apply(base[,seq],2,as.numeric)

# are there proteins with low variability?
no_var = caret::nearZeroVar(base[,seq])
# none

# log transform
base_log <- base %>% modify_if(is_seq(names(.)), log)

# scale by SD
base_log_scale <- base_log
predictors <- colnames(base_log[seq])
for (i in 1:length(predictors)) {
  base_log_scale[,paste0(predictors[i])] <- base_log_scale[,paste0(predictors[i])]/sd(unlist(base_log[,paste0(predictors[i])]))
}

# read in fundus data to cross check progression variable
fundus <- read.csv('/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/TODAY/FUNDUS.csv')
fundus2 <- read.csv('/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/TODAY2/FUNDUS.csv')

fundus <- fundus %>% rename(RELEASEID = releaseid, DAYS = days, EYE = Eye)
all_fundus <- full_join(fundus, fundus2, by = c("RELEASEID", "DAYS", "EYE"))

# summarize TODAY2 fundus progression variable - collapse to 1 row/ppt across eyes
fundus2_summary <- fundus2 %>% select(RELEASEID, DAYS, EYE, DR_PROG3S) 
fundus2_summary <- fundus2_summary %>% group_by(RELEASEID, DAYS) %>% 
  summarise(max_prog = max(DR_PROG3S)) %>% ungroup()
