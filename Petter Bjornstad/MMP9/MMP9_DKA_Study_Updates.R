#### DKA MMP9 Study 

library(dplyr)
library(stringr)
library(ggplot2)
# Mixed-Effects Regression Analysis for MMP9 Claims
library(dplyr)
library(lme4)      # For mixed-effects models
library(lmerTest)  # For p-values in mixed models
library(broom.mixed) # For tidy output of mixed models
library(emmeans)   # For estimated marginal means

library(tidyverse)
library(nlme)
library(emmeans)
library(ggpubr)
library(Hmisc)
library(nephro)
library(knitr)
library(readxl)
library(ggsignif)
library(arsenal)
library(readxl)


raw_data <- data.table::fread('C:/Users/netio/Documents/UofW/Projects/MMP9/Full 16-1403 DKA RedCap Raw Data Set_04-24edited.csv')


names_rawdata <- names(raw_data)
mmp9_index <- str_which(names_rawdata, 'mmp9')

data_dictionary <- readxl::read_xlsx('C:/Users/netio/Documents/UofW/Projects/MMP9/Full 16-1403 DKA RedCap Data Key_04-24edited.xlsx')


final_data <- data.table::fread('C:/Users/netio/Documents/UofW/Projects/MMP9/Full_DKA_Data Set_06-26_final.csv')
kidney_data <- data.table::fread('C:/Users/netio/Documents/UofW/Projects/MMP9/Kidney_DATA_2020-12-03_1614.csv')

qcr <- data.table::fread('C:/Users/netio/Documents/UofW/Projects/MMP9/qcr.csv')

mmp_conc <- readxl::read_xlsx("C:/Users/netio/Documents/UofW/Projects/MMP9/16-1403 MMP9 results 2020.xlsx", skip = 1)
names(mmp_conc) <- c('record_id', 'delete', 'delete1', 'mrn', 'collection_dat', 'time_point', 'dilution', 'raw_result', 
                     'result', 'units')
mmp_conc <- mmp_conc %>% dplyr::select(-delete, -delete1)

mmp_conc$time_point[which(mmp_conc$time_point == 'Follow-up')] <- '3 months'

ggplot(mmp_conc, aes(x=time_point, y = result))+geom_boxplot()+theme_classic()+
  labs(x = 'Time Point', y = 'MMP9 Concentration')



#find variable

# Apply table() to each column and display results
#lapply(final_data, table)

# Or to see them more clearly one by one:
#for(col in names(final_data)) {
#  cat("\nColumn:", col, "\n")
#  print(table(final_data[[col]]))
#}

# To find columns with specific counts (16, 15, 9):
##target_counts <- c(16, 15, 9)
#for(col in names(final_data)) {
#  tab <- table(final_data[[col]])
#  if(any(tab %in% target_counts)) {
#    cat("\nColumn", col, "has target counts:\n")
#    print(tab)
#  }
#}


final_data$severity <- ""
final_data$severity[((final_data$vbg_1 >= 7.3) & (final_data$bmp_1 >= 15))] <- "Mild"
final_data$severity[((final_data$vbg_1 < 7.3) | (final_data$bmp_1 < 15))] <- "Mild"
final_data$severity[((final_data$vbg_1 < 7.2) | (final_data$bmp_1 < 10))] <- "Moderate"
final_data$severity[((final_data$vbg_1 < 7.1) | (final_data$bmp_1 < 5))] <- "Severe"
final_data$severity[final_data$severity == ""] <- "Unknown"


#demographics
demographics <- final_data %>% 
  dplyr::select(sex, age, hba1c = a1c_er, height = ht, weight = wt, #bmi, 
                systolic_bp = sbp, diastolic_bp = dbp, heartrate = hr, a1c_er, a1c_3mo,
                mmp9_0_8hr, mmp9_12_24hr, mmp9_3mo, t1d_status, severity #DKA severity,
                ) %>%
  mutate(bmi = (weight / (height/100)^2),
         age_years = age/12,
         sex_corrected = ifelse(sex == 1, 'Male',
                                ifelse(sex == 2, 'Female', NA)),
         new_t1d = ifelse(t1d_status == 1, 'Yes', 'No')) %>% 
  filter(!is.na(age_years))

demographics$hba1c <- demographics$hba1c %>% str_replace(pattern = '>', replacement = '') %>%
  as.numeric()




desc_table1_fixed <- demographics %>%
  select(age_years, sex_corrected, height, weight, bmi, hba1c, new_t1d, systolic_bp, diastolic_bp, heartrate, severity) %>%
  tbl_summary(
    type = list(
      age_years ~ "continuous",
      bmi ~ "continuous", 
      weight ~ 'continuous',
      height ~ 'continuous',
      hba1c ~ "continuous",
      sex_corrected ~ "categorical",
      new_t1d ~ 'categorical', 
      systolic_bp ~ 'continuous', 
      diastolic_bp ~ 'continuous',
      heartrate ~ 'continuous', 
      severity ~ 'categorical'
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age_years ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age_years ~ "Age (years)",
      sex_corrected ~ "Sex", 
      height ~ "Height (cm)",
      weight ~ "Weight (kg)",
      bmi ~ "BMI (kg/m²)",
      new_t1d ~ 'New Type 1 Diabetes Diagnosis',
      systolic_bp ~ 'Systolic Blood Pressure (mmHg)',
      diastolic_bp ~ 'Diastolic Blood Pressure (mmHg)',
      heartrate ~ 'Heart Rate (bpm)',
      severity ~ 'DKA Severity',
      hba1c ~ "HbA1c (%)"
    ),
    missing_text = "Missing"
  ) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save the table
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Projects/MMP9/MMP9_demographics.png", 
         vwidth = 1200, vheight = 800)















data_set$scr_q <- 107.3/data_set$eGFR_crea
data_set$q <- data_set$screatinine / data_set$scr_q
data_set$aki <- NA
data_set$aki <- data_set$scr_q > 1.5











#table1::table1( ~ sex_corrected + race + age_years + height + weight + bmi + t1d_status + systolic_bp + diastolic_bp + heartrate + severity + a1c_er + mmp9_0_8hr,
#               data = demographics)

library(tidyverse)
library(nlme)
library(emmeans)
library(ggpubr)
library(Hmisc)
library(nephro)
library(knitr)
library(readxl)
library(ggsignif)
library(arsenal)
library(readxl)

#Analysis & Plotting 

dat <- read.csv("C:/Users/netio/Documents/UofW/Projects/MMP9/Full_DKA_Data Set_06-26_final.csv")
dat$male <- 1 - dat$female
dat$age_y <- dat$age/12
dat$height_m <- dat$ht/100
dat$black <- 0

dat$ivinsulinstop <- strptime(dat$ivinsulinstop, format = "%m/%d/%Y %H:%M")
dat$ivinsulin_start <- strptime(dat$ivinsulin_start, format = "%m/%d/%Y %H:%M")


dat$ivinsulin_duration <- as.numeric(dat$ivinsulinstop - dat$ivinsulin_start)
dat$ivinsdurtile <- cut2(dat$ivinsulin_duration, g=3)

dat$severity <- ""

dat$prev_dka <- F
dat$prev_dka <- ifelse(dat$prev_dka_1 == "Yes" , T, dat$prev_dka)
dat$prev_dka <- ifelse(dat$prev_dka_mult == "Yes" , T, dat$prev_dka)

dat$severity[((dat$vbg_1 >= 7.3) & (dat$bmp_1 >= 15))] <- "Mild"
dat$severity[((dat$vbg_1 < 7.3) | (dat$bmp_1 < 15))] <- "Mild"
dat$severity[((dat$vbg_1 < 7.2) | (dat$bmp_1 < 10))] <- "Moderate"
dat$severity[((dat$vbg_1 < 7.1) | (dat$bmp_1 < 5))] <- "Severe"
dat$severity[dat$severity == ""] <- "Unknown"

dat$t1d_status[dat$t1d_status == 1] <- "New"
dat$t1d_status[dat$t1d_status == 2] <- "Known"

dat$a1c_er <- as.numeric(gsub(">","",dat$a1c_er))

dat_time_vars <- dat %>% select(record_id, 
                                t1d_status,
                                prev_dka,
                                severity, 
                                ivinsulin_duration,
                                ivinsdurtile,
                                male,
                                age_y,
                                black,
                                vbg_1,
                                bmp_1,
                                a1c_er,
                                contains("3mo"),
                                contains("0_8"),
                                contains("08"), 
                                contains("12_24"),
                                contains("1224"))


names(dat_time_vars) <- gsub("autoab_3mo_fu___0","autoab_gada_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___1","autoab_gad65_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___2","autoab_miaa_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___3","autoab_ia2_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___4","autoab_znt8_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___5","autoab_ica512_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("12_24","1224", names(dat_time_vars))
names(dat_time_vars) <- gsub("1224hr","1224", names(dat_time_vars))
names(dat_time_vars) <- gsub("0_8","08", names(dat_time_vars))
names(dat_time_vars) <- gsub("08hr","08", names(dat_time_vars))
names(dat_time_vars) <- gsub("08_hr","08", names(dat_time_vars))
names(dat_time_vars) <- gsub("3mo_fu","3mo", names(dat_time_vars))


threemo <- dat_time_vars %>% select(record_id, 
                                    t1d_status,
                                    prev_dka,
                                    severity, 
                                    ivinsulin_duration,
                                    ivinsdurtile,
                                    male,
                                    age_y,
                                    black,
                                    vbg_1,
                                    bmp_1,
                                    a1c_er,
                                    contains("3mo"))
zero8 <- dat_time_vars %>% select(record_id, 
                                  t1d_status,
                                  prev_dka,
                                  severity, 
                                  ivinsulin_duration,
                                  ivinsdurtile,
                                  male,
                                  age_y,
                                  black,
                                  vbg_1,
                                  bmp_1,
                                  a1c_er,
                                  contains("08"))
twelve24 <- dat_time_vars %>% select(record_id, 
                                     t1d_status,
                                     prev_dka,
                                     severity, 
                                     ivinsulin_duration,
                                     ivinsdurtile,
                                     male,
                                     age_y,
                                     black,
                                     vbg_1,
                                     bmp_1,
                                     a1c_er,
                                     contains("1224"))

threemo$time <- "3 months"
zero8$time <- "0-8 hours"
twelve24$time <- "12-24 hours"

names(threemo) <- gsub("_3mo","", names(threemo))
names(zero8) <- gsub("_08","", names(zero8))
names(twelve24) <- gsub("_1224","", names(twelve24))

merge_set <- merge(merge(zero8,twelve24, all = T),threemo, all = T)

# new variable derivations

merge_set$frac_ua <- (merge_set$uua * merge_set$screatinine) / (merge_set$sua * merge_set$ucreatinine)


qcr_df <- read.csv("C:/Users/netio/Documents/UofW/Projects/MMP9/qcr.csv")

get_eGFR <- function(age,sex,creatinine,cystatinc){
  age_year <- floor(age)
  qcr <- ifelse(sex == 1, qcr_df$m_qcr[qcr_df$age == age_year],qcr_df$f_qcr[qcr_df$age == age_year])
  f1 <- creatinine/qcr
  f2 <- 0.5
  f3 <- cystatinc/0.82
  fas <- 107.3 / ((0.5*f1) + (f2*f3))
  fas
}

get_eGFR_crea <- function(age,sex,creatinine){
  age_year <- floor(age)
  qcr <- ifelse(sex == 1, qcr_df$m_qcr[qcr_df$age == age_year],qcr_df$f_qcr[qcr_df$age == age_year])
  f1 <- creatinine/qcr
  fas <- 107.3 / f1
  fas
}


for (i in 1:nrow(merge_set)){
  merge_set$eGFR[i] <- get_eGFR(merge_set$age_y[i],merge_set$male[i],merge_set$screatinine[i],merge_set$cystatinc[i])
}

for (i in 1:nrow(merge_set)){
  merge_set$eGFR_crea[i] <- get_eGFR_crea(merge_set$age_y[i],merge_set$male[i],merge_set$screatinine[i])
}

dat <- read.csv("C:/Users/netio/Documents/UofW/Projects/MMP9/Full_DKA_Data Set_06-26_final.csv")
dat$male <- 1 - dat$female
dat$age_y <- dat$age/12
dat$height_m <- dat$ht/100
dat$black <- 0

dat$ivinsulinstop <- strptime(dat$ivinsulinstop, format = "%m/%d/%Y %H:%M")
dat$ivinsulin_start <- strptime(dat$ivinsulin_start, format = "%m/%d/%Y %H:%M")


dat$ivinsulin_duration <- as.numeric(dat$ivinsulinstop - dat$ivinsulin_start)
dat$ivinsdurtile <- cut2(dat$ivinsulin_duration, g=3)

dat$severity <- ""

dat$prev_dka <- F
dat$prev_dka <- ifelse(dat$prev_dka_1 == "Yes" , T, dat$prev_dka)
dat$prev_dka <- ifelse(dat$prev_dka_mult == "Yes" , T, dat$prev_dka)

dat$severity[((dat$vbg_1 >= 7.3) & (dat$bmp_1 >= 15))] <- "Mild"
dat$severity[((dat$vbg_1 < 7.3) | (dat$bmp_1 < 15))] <- "Mild"
dat$severity[((dat$vbg_1 < 7.2) | (dat$bmp_1 < 10))] <- "Moderate"
dat$severity[((dat$vbg_1 < 7.1) | (dat$bmp_1 < 5))] <- "Severe"
dat$severity[dat$severity == ""] <- "Unknown"

dat$t1d_status[dat$t1d_status == 1] <- "New"
dat$t1d_status[dat$t1d_status == 2] <- "Known"

dat$a1c_er <- as.numeric(gsub(">","",dat$a1c_er))

dat_time_vars <- dat %>% select(record_id, 
                                t1d_status,
                                prev_dka,
                                severity, 
                                ivinsulin_duration,
                                ivinsdurtile,
                                male,
                                age_y,
                                black,
                                vbg_1,
                                bmp_1,
                                a1c_er,
                                contains("3mo"),
                                contains("0_8"),
                                contains("08"), 
                                contains("12_24"),
                                contains("1224"))


names(dat_time_vars) <- gsub("autoab_3mo_fu___0","autoab_gada_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___1","autoab_gad65_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___2","autoab_miaa_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___3","autoab_ia2_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___4","autoab_znt8_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("autoab_3mo_fu___5","autoab_ica512_3mo", names(dat_time_vars))
names(dat_time_vars) <- gsub("12_24","1224", names(dat_time_vars))
names(dat_time_vars) <- gsub("1224hr","1224", names(dat_time_vars))
names(dat_time_vars) <- gsub("0_8","08", names(dat_time_vars))
names(dat_time_vars) <- gsub("08hr","08", names(dat_time_vars))
names(dat_time_vars) <- gsub("08_hr","08", names(dat_time_vars))
names(dat_time_vars) <- gsub("3mo_fu","3mo", names(dat_time_vars))


threemo <- dat_time_vars %>% select(record_id, 
                                    t1d_status,
                                    prev_dka,
                                    severity, 
                                    ivinsulin_duration,
                                    ivinsdurtile,
                                    male,
                                    age_y,
                                    black,
                                    vbg_1,
                                    bmp_1,
                                    a1c_er,
                                    contains("3mo"))
zero8 <- dat_time_vars %>% select(record_id, 
                                  t1d_status,
                                  prev_dka,
                                  severity, 
                                  ivinsulin_duration,
                                  ivinsdurtile,
                                  male,
                                  age_y,
                                  black,
                                  vbg_1,
                                  bmp_1,
                                  a1c_er,
                                  contains("08"))
twelve24 <- dat_time_vars %>% select(record_id, 
                                     t1d_status,
                                     prev_dka,
                                     severity, 
                                     ivinsulin_duration,
                                     ivinsdurtile,
                                     male,
                                     age_y,
                                     black,
                                     vbg_1,
                                     bmp_1,
                                     a1c_er,
                                     contains("1224"))

threemo$time <- "3 months"
zero8$time <- "0-8 hours"
twelve24$time <- "12-24 hours"

names(threemo) <- gsub("_3mo","", names(threemo))
names(zero8) <- gsub("_08","", names(zero8))
names(twelve24) <- gsub("_1224","", names(twelve24))

merge_set <- merge(merge(zero8,twelve24, all = T),threemo, all = T)

# new variable derivations

merge_set$frac_ua <- (merge_set$uua * merge_set$screatinine) / (merge_set$sua * merge_set$ucreatinine)


qcr_df <- read.csv("C:/Users/netio/Documents/UofW/Projects/MMP9/qcr.csv")

get_eGFR <- function(age,sex,creatinine,cystatinc){
  age_year <- floor(age)
  qcr <- ifelse(sex == 1, qcr_df$m_qcr[qcr_df$age == age_year],qcr_df$f_qcr[qcr_df$age == age_year])
  f1 <- creatinine/qcr
  f2 <- 0.5
  f3 <- cystatinc/0.82
  fas <- 107.3 / ((0.5*f1) + (f2*f3))
  fas
}

get_eGFR_crea <- function(age,sex,creatinine){
  age_year <- floor(age)
  qcr <- ifelse(sex == 1, qcr_df$m_qcr[qcr_df$age == age_year],qcr_df$f_qcr[qcr_df$age == age_year])
  f1 <- creatinine/qcr
  fas <- 107.3 / f1
  fas
}


for (i in 1:nrow(merge_set)){
  merge_set$eGFR[i] <- get_eGFR(merge_set$age_y[i],merge_set$male[i],merge_set$screatinine[i],merge_set$cystatinc[i])
}

for (i in 1:nrow(merge_set)){
  merge_set$eGFR_crea[i] <- get_eGFR_crea(merge_set$age_y[i],merge_set$male[i],merge_set$screatinine[i])
}

data_set <- merge_set

data_set$scr_q <- 107.3/data_set$eGFR_crea
data_set$q <- data_set$screatinine / data_set$scr_q
data_set$aki <- NA
data_set$aki <- data_set$scr_q > 1.5


data_set$sngal <- data_set$sngal/1000
data_set$sykl40 <- data_set$sykl40/1000

aki_frame <- NULL
for (i in unique(data_set$record_id)){
  id_row <- data_set %>% filter(record_id == i) %>% filter(time == '0-8 hours' | time == '12-24 hours') %>% select(record_id, time, aki)
  aki <- sum(id_row$aki, na.rm=T) > 0
  aki_frame <- rbind(aki_frame, c(i,aki))
}
aki_frame <- as.data.frame(aki_frame)
names(aki_frame) <- c("record_id","aki_0_24")

data_set <- merge(data_set,aki_frame)
data_set$aki_0_24 <- as.logical(data_set$aki_0_24)

data_set$aki_0_24 <- ifelse(data_set$aki_0_24, "Yes","No")

mmp9_dat <- read_excel("C:/Users/netio/Documents/UofW/Projects/MMP9/16-1403 MMP9 results 2020.xlsx", sheet = 1, skip = 1)
mmp9_dat <- mmp9_dat %>% select(SID, `Time Point`, `Actual Result`)
names(mmp9_dat) <- c("record_id", "time", "mmp9_actual")
mmp9_dat$time[mmp9_dat$time == "Follow-up"] <- "3 months"

data_set <- merge(data_set,mmp9_dat, all = T)
data_set$mmp9_egfr <- (data_set$mmp9_actual*100) / data_set$eGFR
data_set$fru[data_set$fru > 100] <- NA

cd93_dat <- read_excel("C:/Users/netio/Documents/UofW/Projects/MMP9/DKA ELISA cd93 stat.xlsx", sheet = 1)
names(cd93_dat) <- c("record_id", "time", "cd93")

data_set <- merge(data_set,cd93_dat, all = T)
data_set$cd93_egfr <- (data_set$cd93*100) / data_set$eGFR
data_set$fru_crea <- data_set$fru*100/data_set$crea




graph1 <- ggplot(data_set %>% filter(time %in% c('0-8 hours', '3 months')), aes(x=time, y = mmp9_actual, fill = severity))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = c("#fff9ec", "#fcb1a6", "#fb6376"))+
  labs(x = '', y = 'MMP9', fill = 'Severity')+
  theme(text = element_text(size = 16))

  
graph2 <- ggplot(data_set %>% filter(time %in% c('0-8 hours', '3 months')), aes(x=time, y = mmp9_actual, fill = ivinsdurtile))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = c("#fff9ec", "#fcb1a6", "#fb6376"))+
  labs(x = '', y = 'MMP9', fill = 'Days IV Insulin')+
  theme(text = element_text(size = 16))

graph3 <- ggplot(data_set %>% filter(time %in% c('0-8 hours', '3 months')), aes(x=time, y = mmp9_actual, fill = aki_0_24))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = c('#4CAF50', '#F44336'))+
  labs(x = '', y = 'MMP9', fill = 'AKI')+
  theme(text = element_text(size = 16))

graph4 <- ggplot(data_set %>% filter(time %in% c('0-8 hours', '3 months')), aes(x=time, y = mmp9_actual, fill = t1d_status))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = c('#4CAF50', '#F44336'))+
  labs(x = '', y = 'MMP9', fill = 'T1D Status')+
  theme(text = element_text(size = 16))



combined_boxplots <- (graph1 | graph2) / 
  (graph3 | graph4)

# Add overall title
combined_boxplots <- combined_boxplots + 
  plot_annotation(title = "",
                  tag_levels = 'A', 
                  theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))

# Display the combined plot
print(combined_boxplots)

# Save the combined plot
ggsave("C:/Users/netio/Documents/UofW/Projects/MMP9/MMP9_boxplots_by_groups.png", combined_boxplots,
       width = 12, height = 10, dpi = 300, units = "in")




##################################################################################### Analysis

# Mixed-Effects Regression Analysis for MMP9 Claims
library(dplyr)
library(lme4)      # For mixed-effects models
library(lmerTest)  # For p-values in mixed models
library(broom.mixed) # For tidy output of mixed models
library(emmeans)   # For estimated marginal means

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Create analysis dataset
analysis_data <- data_set %>%
  filter(time %in% c('0-8 hours', '3 months')) %>%
  filter(!is.na(mmp9_actual) & !is.na(record_id)) %>%
  mutate(
    time_binary = ifelse(time == '0-8 hours', 0, 1),  # 0 = DKA, 1 = 3 months
    time_factor = factor(time, levels = c('0-8 hours', '3 months'))
  )

# Check data structure
cat("Data structure:\n")
cat("Total observations:", nrow(analysis_data), "\n")
cat("Unique participants:", length(unique(analysis_data$record_id)), "\n")
cat("Observations per time point:\n")
print(table(analysis_data$time))

# Check how many participants have data at both time points
participant_timepoints <- analysis_data %>%
  group_by(record_id) %>%
  summarise(
    n_timepoints = n(),
    has_dka = any(time == '0-8 hours'),
    has_followup = any(time == '3 months'),
    .groups = 'drop'
  )

cat("\nParticipant time point coverage:\n")
cat("Participants with both time points:", sum(participant_timepoints$has_dka & participant_timepoints$has_followup), "\n")
cat("Participants with DKA only:", sum(participant_timepoints$has_dka & !participant_timepoints$has_followup), "\n")
cat("Participants with follow-up only:", sum(!participant_timepoints$has_dka & participant_timepoints$has_followup), "\n")

# ============================================================================
# CLAIM 1: MMP9 concentrations higher during DKA than 3 months follow-up
# Expected: mean±SE: 1504.6±137 vs. 668.7±159 ng/mL, p=0.0003
# ============================================================================

# Mixed-effects model for time comparison
model1 <- lmer(mmp9_egfr ~ time_factor + (1|record_id), data = analysis_data)

# Model summary
print("Mixed-effects model summary:")
print(summary(model1))

# Get estimated marginal means (least squares means)
emm1 <- emmeans(model1, ~ time_factor)
print("\nEstimated marginal means:")
print(emm1)

# Pairwise comparisons
pairs1 <- pairs(emm1)
print("\nPairwise comparison:")
print(pairs1)

# Extract coefficients with confidence intervals
coef1 <- tidy(model1, effects = "fixed", conf.int = TRUE)
print("\nFixed effects coefficients:")
print(coef1)

cat("\nExpected from paper: 0-8 hours: 1504.6±137, 3 months: 668.7±159, p=0.0003\n\n")

# ============================================================================
# CLAIM 2: At 0-8 hours, participants with AKI had higher MMP9
# Expected: 2256.9±310.1 vs. 1344.7±143.5 ng/mL, p=0.01
# ============================================================================

cat(rep("=", 70) + "\n")
cat("CLAIM 2: MMP9 by AKI Status at 0-8 hours\n")
cat(rep("=", 70) + "\n")

# Filter for DKA time point only and check AKI variable
dka_data <- analysis_data %>% 
  filter(time == '0-8 hours', !is.na(aki_0_24))

cat("AKI variable levels at 0-8 hours:", unique(dka_data$aki_0_24), "\n")
cat("Sample sizes by AKI status:", table(dka_data$aki_0_24), "\n")

if(length(unique(dka_data$aki_0_24)) >= 2) {
  # Mixed-effects model for AKI comparison (though with single time point, random effect may not be needed)
  # But keeping it for consistency and in case some participants have multiple DKA measurements
  model2 <- aov(mmp9_egfr ~ aki_0_24, data = dka_data)
  
  print("Mixed-effects model summary for AKI:")
  print(summary(model2))
  
  # Estimated marginal means
  emm2 <- emmeans(model2, ~ aki_0_24)
  print("\nEstimated marginal means by AKI status:")
  print(emm2)
  
  # Pairwise comparisons
  pairs2 <- pairs(emm2)
  print("\nPairwise comparison:")
  print(pairs2)
  
  # Coefficients
  coef2 <- tidy(model2, effects = "fixed", conf.int = TRUE)
  print("\nFixed effects coefficients:")
  print(coef2)
  
} else {
  cat("Insufficient AKI groups for comparison\n")
}

#3 months
dka_data <- analysis_data %>% 
  filter(time == '3 months', !is.na(aki_0_24))

cat("AKI variable levels at 0-8 hours:", unique(dka_data$aki_0_24), "\n")
cat("Sample sizes by AKI status:", table(dka_data$aki_0_24), "\n")


# Mixed-effects model for AKI comparison (though with single time point, random effect may not be needed)
# But keeping it for consistency and in case some participants have multiple DKA measurements
model2 <- aov(mmp9_egfr ~ aki_0_24, data = dka_data)

print("Mixed-effects model summary for AKI:")
print(summary(model2))

# Estimated marginal means
emm2 <- emmeans(model2, ~ aki_0_24)
print("\nEstimated marginal means by AKI status:")
print(emm2)

# Pairwise comparisons
pairs2 <- pairs(emm2)
print("\nPairwise comparison:")
print(pairs2)

# Coefficients
coef2 <- tidy(model2, effects = "fixed", conf.int = TRUE)
print("\nFixed effects coefficients:")
print(coef2)





##################################### Grouping analysis 


ngal_model <- lme(mmp9_egfr ~ time*scopeptin, data = data_set,
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ time*sua, data = data_set,
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)



#aki yes
ngal_model <- lme(mmp9_egfr ~ time*scopeptin, data = data_set %>% filter(aki_0_24 == 'Yes'),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ time*sua, data = data_set %>% filter(aki_0_24 == 'Yes'),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


#aki no

ngal_model <- lme(mmp9_egfr ~ time*scopeptin, data = data_set %>% filter(aki_0_24 == 'No'),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ time*sua, data = data_set %>% filter(aki_0_24 == 'No'),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)
















# Filter for DKA time point only and check AKI variable
dka_data <- analysis_data %>% 
  filter(time == '0-8 hours') %>%
  filter(aki_0_24 == 'Yes')

  model_aki_baseline <- lm(mmp9_egfr ~ scopeptin, data = dka_data)
  summary(model_aki_baseline)
  
  model_aki_baseline <- lm(mmp9_egfr ~ sua, data = dka_data)
  summary(model_aki_baseline)
  
  # Estimated marginal means
  emm2 <- emmeans(model_aki_baseline, ~ scopeptin)

  # Pairwise comparisons
  pairs2 <- pairs(emm2)
 
  
  # Coefficients
  coef2 <- tidy(model, effects = "fixed", conf.int = TRUE)
  print("\nFixed effects coefficients:")
  print(coef2)
  

#3 months
dka_data <- analysis_data %>% 
  filter(time == '3 months', !is.na(aki_0_24))

cat("AKI variable levels at 0-8 hours:", unique(dka_data$aki_0_24), "\n")
cat("Sample sizes by AKI status:", table(dka_data$aki_0_24), "\n")


# Mixed-effects model for AKI comparison (though with single time point, random effect may not be needed)
# But keeping it for consistency and in case some participants have multiple DKA measurements
model2 <- aov(mmp9_egfr ~ aki_0_24, data = dka_data)

print("Mixed-effects model summary for AKI:")
print(summary(model2))

# Estimated marginal means
emm2 <- emmeans(model2, ~ aki_0_24)
print("\nEstimated marginal means by AKI status:")
print(emm2)

# Pairwise comparisons
pairs2 <- pairs(emm2)
print("\nPairwise comparison:")
print(pairs2)

# Coefficients
coef2 <- tidy(model2, effects = "fixed", conf.int = TRUE)
print("\nFixed effects coefficients:")
print(coef2)






























  
  
  
  
  
  
  
  
  
  



library(dplyr)
library(tidyr)
library(ggplot2)

# Method 1: Keep grouping variables in the analysis
mmp9_change <- data_set %>%
  # Keep only the timepoints we need
  filter(time %in% c("0-8 hours", "3 months")) %>%
  # Select relevant columns INCLUDING grouping variables
  select(record_id, time, mmp9_actual, ivinsdurtile, aki_0_24, severity, t1d_status) %>%
  # Remove rows with missing MMP9 data
  filter(!is.na(mmp9_actual)) %>%
  # Group by participant and grouping variables to preserve them
  group_by(record_id, ivinsdurtile, aki_0_24, severity, t1d_status) %>%
  # Check that we have both timepoints for this participant
  filter(n() == 2) %>%
  # Calculate the change within each group
  summarise(
    mmp9_baseline = mmp9_actual[time == "0-8 hours"],
    mmp9_3months = mmp9_actual[time == "3 months"],
    mmp9_change = mmp9_3months - mmp9_baseline,
    .groups = 'drop'
  )

# View the results
print("MMP9 change with grouping variables:")
print(head(mmp9_change))

# Overall density plot
overall_density <- ggplot(mmp9_change, aes(x = mmp9_change)) + 
  geom_density() + 
  theme_classic() +
  labs(title = "Overall Distribution of MMP9 Changes",
       x = "MMP9 Change (3 months - baseline)",
       y = "Density")

print(overall_density)

# Density plot by IV insulin duration
iv_density <- ggplot(mmp9_change, aes(x = mmp9_change, fill = ivinsdurtile)) + 
  geom_density(alpha = 0.6) + 
  theme_classic() +
  labs(title = "MMP9 Changes by IV Insulin Duration",
       x = "MMP9 Change (3 months - baseline)",
       y = "Density",
       fill = "IV Insulin Duration") +
  theme(legend.position = "bottom")

print(iv_density)

# Density plot by AKI status
aki_density <- ggplot(mmp9_change, aes(x = mmp9_change, fill = aki_0_24)) + 
  geom_density(alpha = 0.6) + 
  theme_classic() +
  labs(title = "MMP9 Changes by AKI Status",
       x = "MMP9 Change (3 months - baseline)",
       y = "Density",
       fill = "AKI Status") +
  theme(legend.position = "bottom")

print(aki_density)

# Density plot by T1D status
t1d_density <- ggplot(mmp9_change, aes(x = mmp9_change, fill = t1d_status)) + 
  geom_density(alpha = 0.6) + 
  theme_classic() +
  labs(title = "MMP9 Changes by T1D Status",
       x = "MMP9 Change (3 months - baseline)",
       y = "Density",
       fill = "T1D Status") +
  theme(legend.position = "bottom")

print(t1d_density)

# Density plot by severity
severity_density <- ggplot(mmp9_change, aes(x = mmp9_change, fill = severity)) + 
  geom_density(alpha = 0.6) + 
  theme_classic() +
  labs(title = "MMP9 Changes by Severity",
       x = "MMP9 Change (3 months - baseline)",
       y = "Density",
       fill = "Severity") +
  theme(legend.position = "bottom")

print(severity_density)


combined_density_plots <- (iv_density | aki_density) / 
  (severity_density | t1d_density)

# Add overall title
combined_density_plots <- combined_density_plots + 
  plot_annotation(title = "MMP9 Changes (3 months - baseline) by Clinical Variables",
                  theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))

# Display the combined plot
print(combined_density_plots)

# Save the combined plot
ggsave("C:/Users/netio/Documents/UofW/Projects/MMP9/MMP9_density_plots_by_groups.png", combined_density_plots, 
       width = 12, height = 10, dpi = 300, units = "in")



























plot_df <- mmp_conc %>% left_join(final_data %>% dplyr::select(record_id, severity, ivinsulin_start, ivinsulinstop), 
                                  by = 'record_id')

ggplot(plot_df, aes(x=time_point, y = result, fill = severity))+geom_boxplot()+theme_classic()+
  labs(x = 'Time Point', y = 'MMP9 Concentration', )




#analysis dataframe

analysis_df <- final_data %>% 
  dplyr::select(record_id, scopeptin_0_8hr, scopeptin_3mo, 
                serum_uric_0_8hr = sua_0_8hr, serum_uric_3mo = sua_3mo,  
                urine_uric_0_8hr = uua_0_8_hr, urine_uric_3mo = uua_3mo, 
                ivinsulin_start, ivinsulinstop, severity)

# Reshape the data
library(dplyr)
library(tidyr)

long_data <- analysis_df %>%
  # Convert all measurement columns to character to handle mixed types
  mutate(across(ends_with(c("_0_8hr", "_3mo")), as.character)) %>%
  pivot_longer(
    cols = -c(record_id, ivinsulin_start, ivinsulinstop, severity),  # Add severity here
    names_to = c("variable", "timepoint"), 
    names_pattern = "(.+)_(0_8hr|3mo)",
    values_to = "value"
  ) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("0_8hr", "3mo"), 
                       labels = c("0-8 hours", "3 months")),
    # Convert back to numeric where possible, keeping character for ">1600" etc.
    value_numeric = as.numeric(value)
  )

# Handle the >1600 values
long_data$value_numeric[which(long_data$value == '>1600')] <- 1600

# View the result
head(long_data)




plot1 <- ggplot(long_data %>% filter(variable == 'mmp9'), 
                aes(x=timepoint, y = as.numeric(value_numeric), fill = severity)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  labs(x='Time Point', y = 'MMP9') +
  theme_classic()

#Why are there so few data? 


ggplot(long_data %>% filter(variable == 'mmp9'), 
       aes(x=timepoint, y = as.numeric(value_numeric), fill = severity)) +
  geom_violin()+geom_boxplot(width = 0.3, fill = NA)+
  labs(x='Time Point', y = 'MMP9') +
  theme_classic()


















