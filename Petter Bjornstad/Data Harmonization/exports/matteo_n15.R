##### data pull for Matteo: 5/26/2026 #####
library(dplyr)
library(purrr)

dat<-read.csv("/Users/kristenmiller/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv")

# for IT_19, only want baseline:
dat<-subset(dat,!(dat$record_id=="IT_19" & dat$visit=="12_months_post_surgery"))
dat<-subset(dat,!(dat$record_id=="IT_19" & dat$visit=="3_months_post_surgery"))

# for each patient, condense to one row, as these are all 1 visit studies (or baseline for IT2D)
#if not numeric column, take value in last row
#if numeric column, average across all rows
dat.harm.matteo <- dat %>% filter(record_id %in%c("RH2-21-T", #
                                          "RH2-14-T", #
                                          "RH-91-T", ## RH-91 -> is RH-91-T?
                                          "RH-93-T", ## RH-93 -> is RH-93-T?
                                          "RH2-07-O", #
                                          #"IT2D_19_BL", #IT2D_19_BL -> is IT_19?
                                          "IT_19", 
                                          "RH-67-T",  # 
                                          "CRC-11", #
                                          "RH2-19-T", #
                                          "RH2-38-T", #
                                          "RH2-23-T", #
                                          "RH2-22-T", #
                                          "RH2-53-T", #
                                          "RH2-51-T", #
                                          "RH2-59-T")) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>%
  arrange(record_id) %>%
  select(-mrn)

write.csv(dat.harm.matteo,"/Users/kristenmiller/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Exports/matteo_n15_05282026.csv")

#RH-67-T is the same as RH2-19-T: 
# I would only include RH2-19-T unless there's aspects of the data that's missing in RH but is available in RH2 
# (for example, if RH didn't collect MRI but RH2 did).. but I think he's mostly interested in the clinical stuff 
#   that would all be available as part of RH. If not, I would start by looking at the date of visit for RH-67-T 
#   and RH2-19-T to see if they are the same visit, 
# and if they're not the same visit, RH2-19-T is probably a follow up visit from RH-67-T, so their data would 
# not be the same and would probably be ~ 3 years apart

#RH-67-T date: 2019-12-16
#RH2-19-T date: 2023-03-13

