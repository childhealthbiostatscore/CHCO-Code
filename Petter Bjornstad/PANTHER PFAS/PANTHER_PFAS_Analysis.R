#1. Set Up ----
##a. Load Libraries & Directores----
#Libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Directories
computer <- "mac studio"
# computer <- "mac laptop"
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
source("Functions.R")

##b. Load Data----
#Load Metadata
dat <- read.csv(fs::path(dir.dat,"Data Harmonization","Data Clean","harmonized_dataset.csv"),na.strings = "")
lod <- read.csv(fs::path(dir.dat,"/PANTHER/Data_Raw/PFAS Emory 2025/8- PFAS Quantification/PFAS_LODs.csv"))

##c. Clean & Format Data----
# dat <- dat %>% 
#   mutate(group=ifelse(record_id=="PAN-110-O","Obese Control",
#                       ifelse(record_id=="PAN-61-C","Lean Control",
#                              ifelse(record_id=="PAN-71-C","Lean Control",group))))

pfas <- c("N.EtFOSAA","N.MeFOSAA","PFBA","PFDA","PFDoA","PFHpA","PFHxA","PFNA","PFOA","PFBS","PFPeA","PFTeDA","PFTrDA","PFUnA","PFBS","PFDoS","PFHps","PFHxS","PFNS","PFOS","PFOSA","PFPeAS")
dat_collapsed1 <- dat %>%
  filter(study=="PANTHER") %>% 
  filter(visit=="baseline" | visit=="screening") %>% 
  # dplyr::select(record_id, date, rh2_id, visit,mrn, group, age, sex, bmi, hba1c,albuminuria_cat,microalbumin_u,starts_with("egfr"),cystatin_c_s,starts_with("creat"),all_of(pfas)) %>% 
  dplyr::select(record_id, dob, date, rh2_id, visit, mrn, group, age, sex, bmi,race,ethnicity,race_ethnicity, hba1c,albuminuria_cat,acr_u,microalbumin_u,eGFR_CKiD_U25_avg,cystatin_c_s,creatinine_s,creatinine_u) %>% 
  # group_by(record_id, visit) %>%
  # tidyr::fill(date, .direction = "updown") %>% ungroup() %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>% 
  mutate(visit="baseline_screening") %>% 
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other"))
dat_collapsed2 <- dat %>%
  filter(study=="PANTHER") %>% 
  filter(visit=="year_1" | visit=="year_2") %>% 
  # dplyr::select(record_id, date, rh2_id, visit,mrn, group, age, sex, bmi, hba1c,albuminuria_cat,microalbumin_u,starts_with("egfr"),cystatin_c_s,starts_with("creat"),all_of(pfas)) %>% 
  dplyr::select(record_id, dob, date, rh2_id, visit,mrn, group, age, sex, bmi, race,ethnicity, race_ethnicity, hba1c,acr_u,albuminuria_cat,microalbumin_u,eGFR_CKiD_U25_avg,cystatin_c_s,creatinine_s,creatinine_u) %>% 
  # group_by(record_id, visit) %>%
  # tidyr::fill(date, .direction = "updown") %>% ungroup() %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>% 
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" &
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other"))
meta <- rbind(dat_collapsed1,dat_collapsed2)

pfas_dat <- dat %>% 
  filter(study=="PANTHER") %>% 
  filter(visit=="baseline") %>% 
  dplyr::select("record_id",all_of(pfas)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>% 
  filter(!is.na(N.EtFOSAA)) #91 participants with baseline pfas data

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
PFOA <- length(which(pfas_dat$PFOA < lod$`PFOA`))/length(unique(pfas_dat$record_id))
PFOA < 2/3 #True, impute
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
PFNS < 2/3 #FALSE, exclude
PFOS <- length(which(pfas_dat$PFOS < lod$`PFOS`))/length(unique(pfas_dat$record_id))
PFOS < 2/3 #True, impute
PFOSA <- length(which(pfas_dat$PFOSA < lod$`PFOSA`))/length(unique(pfas_dat$record_id))
PFOSA < 2/3 #False, exclude
PFPeAS <- length(which(pfas_dat$PFPeAS < lod$`PFPeAS`))/length(unique(pfas_dat$record_id))
PFPeAS < 2/3 #True, impute

pfas_to_impute <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA","PFOA" ,"PFBS", 
                    "PFHps", "PFHxS", "PFOS", "PFPeAS")
pfas_imputed <- pfas_dat %>% 
  dplyr::select(c("record_id",all_of(pfas_to_impute)))
for (pfas in pfas_to_impute){
  pfas_imputed[[pfas]] <- ifelse(pfas_imputed[[pfas]]<lod[[pfas]],(lod[[pfas]]/sqrt(2)),pfas_imputed[[pfas]])
}

#Save formatted and imputed pfas data
# saveRDS(pfas_imputed,fs::path(dir.dat,"PFAS_Data_Imputed_11_03.rds"))

# pfas_to_impute
table1(~N.MeFOSAA+PFDA+PFHpA+PFNA+PFOA+PFBS+PFHps+PFHxS+PFOS+PFPeAS,data=pfas_imputed,render.continuous = render.geometric)

#Format metadata
# meta <- dat_collapsed
dat <- tidylog::left_join(pfas_imputed,meta,by="record_id")

#Particpants 
dat_baseline <- dat %>% 
  filter(visit=="baseline_screening")
table1(~age+sex+bmi+race_ethnicity_condensed+hba1c+microalbumin_u+acr_u+creatinine_s+creatinine_u+eGFR_CKiD_U25_avg|group,data=dat_baseline)

#Folowoup Year 1
dat_followup1 <- dat %>% 
  filter(visit=="year_1") %>% 
  dplyr::select(c("record_id","microalbumin_u","acr_u","eGFR_CKiD_U25_avg")) %>% 
  dplyr::rename(microalbumin_u_yr1=microalbumin_u,
                acr_u_yr1=acr_u,
                eGFR_CKiD_U25_avg_yr1=eGFR_CKiD_U25_avg)

#Followup Year 2
dat_followup2 <- dat %>% 
  filter(visit=="year_2") %>% 
  dplyr::select(c("record_id","microalbumin_u","acr_u","eGFR_CKiD_U25_avg")) %>% 
  dplyr::rename(microalbumin_u_yr2=microalbumin_u,
                acr_u_yr2=acr_u,
                eGFR_CKiD_U25_avg_yr2=eGFR_CKiD_U25_avg)

dat_all <- tidylog::left_join(dat_baseline,dat_followup1,by="record_id")
dat_all <- tidylog::left_join(dat_all,dat_followup2,by="record_id")

table1(~N.MeFOSAA+PFDA+PFHpA+PFNA+PFBS+PFHps+PFHxS+PFNS+PFOS+PFPeAS,data=pfas_imputed,render.continuous = render.geometric)
table1(~age+sex+bmi+race_ethnicity_condensed+hba1c+microalbumin_u+microalbumin_u_yr1+microalbumin_u_yr2+acr_u+
         acr_u_yr1+acr_u_yr2+creatinine_s+creatinine_u+eGFR_CKiD_U25_avg+eGFR_CKiD_U25_avg_yr1+eGFR_CKiD_U25_avg_yr2|group,data=dat_all2)
table1(~N.MeFOSAA+PFDA+PFHpA+PFNA+PFBS+PFHps+PFHxS+PFOA+PFOS+PFPeAS|group,data=dat_all2,render.continuous = render.geometric)
# #Filter out albumin and acru outlier
# dat_baseline2 <- dat_baseline %>% 
#   mutate(microalbumin_u=ifelse(microalbumin_u>3000,NA,microalbumin_u)) %>% 
#   mutate(acr_u=ifelse(acr_u>2367,NA,acr_u))
# #Filter out albumin outlier and acru outlier
# dat2 <- dat %>% 
#   mutate(microalbumin_u=ifelse(microalbumin_u>3000,NA,microalbumin_u)) %>% 
#   mutate(acr_u=ifelse(acr_u>2367,NA,acr_u))
dat_all2 <- dat_all %>% 
  mutate(microalbumin_u=ifelse(microalbumin_u>3000,NA,microalbumin_u)) %>% 
  mutate(acr_u=ifelse(acr_u>2367,NA,acr_u))

#Exposure assessment: 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/PANTHER/Results/PFAS Results/pfas_acru_plots_no_outlier_withPFOA.pdf", width = 8, height = 6)
#Baseline Analysis
pfas_vars <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA", "PFBS", 
               "PFHps", "PFHxS", "PFOA","PFOS", "PFPeAS")
for (pfas in pfas_vars) {
  m0 <- as.formula(paste0("acr_u~",pfas,"+age+sex"))
  m1 <- lm(m0,dat=dat_baseline2)
  beta <- summary(m1)$coef[2,1]
  pval <- summary(m1)$coef[2,4]
  
  p <- ggplot(data=dat_baseline2,aes(x=.data[[pfas]],y=acr_u))+
    geom_point()+
    geom_smooth(method="lm",se=F)+
    theme_bw()+
    labs(title = c("Association between PFAS and ACRu, adj. Age & Sex"),
         subtitle = paste0("Beta = ", round(beta,3), ", pval = ", round(pval,3)),
         x = pfas,  # This will use the variable name from pfas
         y = "ACRu")
  print(p)
  
}
dev.off()

pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/PANTHER/Results/PFAS Results/pfas_microalbumin_plots_no_outlier.pdf", width = 8, height = 6)
#Baseline Analysis
pfas_vars <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA", "PFBS", 
                    "PFHps", "PFHxS", "PFOS", "PFPeAS")
for (pfas in pfas_vars) {
  m0 <- as.formula(paste0("microalbumin_u~",pfas,"+age+sex"))
  m1 <- lm(m0,dat=dat_baseline2)
  beta <- summary(m1)$coef[2,1]
  pval <- summary(m1)$coef[2,4]
  
  p <- ggplot(data=dat_baseline2,aes(x=.data[[pfas]],y=microalbumin_u))+
    geom_point()+
    geom_smooth(method="lm",se=F)+
    theme_bw()+
    labs(title = c("Association between PFAS and Microalbuminuria, adj. Age & Sex"),
         subtitle = paste0("Beta = ", round(beta,3), ", pval = ", round(pval,3)),
         x = pfas,  # This will use the variable name from pfas
         y = "Microalbumin (U)")
  print(p)
  
}
dev.off()

pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/PANTHER/Results/PFAS Results/pfas_eGFR_plots.pdf", width = 8, height = 6)
#Baseline Analysis
pfas_vars <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA", "PFBS", 
               "PFHps", "PFHxS", "PFOS", "PFPeAS")
for (pfas in pfas_vars) {
  m0 <- as.formula(paste0("eGFR_CKiD_U25_avg~",pfas,"+age+sex"))
  m1 <- lm(m0,dat=dat_baseline2)
  beta <- summary(m1)$coef[2,1]
  pval <- summary(m1)$coef[2,4]
  
  p <- ggplot(data=dat_baseline2,aes(x=.data[[pfas]],y=eGFR_CKiD_U25_avg))+
    geom_point()+
    geom_smooth(method="lm",se=F)+
    theme_bw()+
    labs(title = c("Association between PFAS and eGFR, adj. Age & Sex"),
         subtitle = paste0("Beta = ", round(beta,3), ", pval = ", round(pval,3)),
         x = pfas,  # This will use the variable name from pfas
         y = "eGFR")
  print(p)
  
}
dev.off()

#Mixed Model Analysis with Time
dat2$time_numeric <- as.numeric(factor(dat2$visit)) - 1  # Start at 0 for baseline
m0 <- as.formula(paste0("acr_u ~ time + (1|record_id)"))
model_simple <- lmer(m0, data = dat2)
summary(model_simple)$coef

# Scatter plot with individual trajectories
p3 <- ggplot(dat2, aes(x = time, y = acr_u)) +
  geom_point(aes(color = factor(visit)), alpha = 0.6) +
  geom_path(aes(group = record_id), alpha = 0.2) +
  geom_smooth(method = "lm", color = "red4", linewidth= 0.5,se=F) +
  scale_color_viridis_d(name = "Time") +
  theme_minimal() +
  labs(title = "ACRu Levels over Time",
       subtitle = "Beta = -3.15, P-value = 0.418",
       x = "Time",
       y = "acr_u Level") +
  theme(plot.title = element_text(face = "bold"))
# legend.position="none")

print(p3)

#Baseline PFAS with Follow-Up ACRu
#Exposure assessment: 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/PANTHER/Results/PFAS Results/baseline_pfas_followup_acru_yr1.pdf", width = 8, height = 6)
#Baseline Analysis
pfas_vars <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA", "PFBS", 
               "PFHps", "PFHxS", "PFOA","PFOS", "PFPeAS")
for (pfas in pfas_vars) {
  m0 <- as.formula(paste0("acr_u_yr1~",pfas,"+age+sex+acr_u"))
  m1 <- lm(m0,dat=dat_all2)
  beta <- summary(m1)$coef[2,1]
  pval <- summary(m1)$coef[2,4]
  
  p <- ggplot(data=dat_all2,aes(x=.data[[pfas]],y=acr_u_yr1))+
    geom_point()+
    geom_smooth(method="lm",se=F)+
    theme_bw()+
    labs(title = c("Association between Baseline PFAS and ACRu 1 Year Follow-Up, adj. Age, Sex and Baseline ACRu"),
         subtitle = paste0("Beta = ", round(beta,3), ", pval = ", round(pval,3)),
         x = pfas,  # This will use the variable name from pfas
         y = "ACRu (1 Yr Follow-Up)")
  print(p)
  
}
dev.off()

pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/PANTHER/Results/PFAS Results/baseline_pfas_followup_acru_yr2.pdf", width = 8, height = 6)
#Baseline Analysis
pfas_vars <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA", "PFBS", 
               "PFHps", "PFHxS", "PFOA","PFOS", "PFPeAS")
for (pfas in pfas_vars) {
  m0 <- as.formula(paste0("acr_u_yr1~",pfas,"+age+sex+acr_u"))
  m1 <- lm(m0,dat=dat_all2)
  beta <- summary(m1)$coef[2,1]
  pval <- summary(m1)$coef[2,4]
  
  p <- ggplot(data=dat_all2,aes(x=.data[[pfas]],y=acr_u_yr2))+
    geom_point()+
    geom_smooth(method="lm",se=F)+
    theme_bw()+
    labs(title = c("Association between Baseline PFAS and ACRu 2 Year Follow-Up, adj. Age, Sex and Baseline ACRu"),
         subtitle = paste0("Beta = ", round(beta,3), ", pval = ", round(pval,3)),
         x = pfas,  # This will use the variable name from pfas
         y = "ACRu (2 Yr Follow-Up)")
  print(p)
  
}
dev.off()