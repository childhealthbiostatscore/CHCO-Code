#### DKA MMP9 Study 

library(dplyr)
library(stringr)
library(ggplot2)
library(lme4)     
library(lmerTest) 
library(broom.mixed) 
library(emmeans) 

library(tidyverse)
library(nlme)
library(ggpubr)
library(Hmisc)
library(nephro)
library(knitr)
library(readxl)
library(ggsignif)
library(arsenal)
library(readxl)
library(gtsummary)
library(gt)
library(tidyr)
library(patchwork)
library(ggsignif)





working_folder <- '/Users/netio/Documents/UofW/Projects/MMP9/'

setwd(working_folder)

raw_data <- data.table::fread('Full 16-1403 DKA RedCap Raw Data Set_04-24edited.csv')


names_rawdata <- names(raw_data)
mmp9_index <- str_which(names_rawdata, 'mmp9')

data_dictionary <- readxl::read_xlsx('Full 16-1403 DKA RedCap Data Key_04-24edited.xlsx')


final_data <- data.table::fread('Full_DKA_Data Set_06-26_final.csv')
kidney_data <- data.table::fread('Kidney_DATA_2020-12-03_1614.csv')

qcr <- data.table::fread('qcr.csv')

#mmp_conc <- readxl::read_xlsx("16-1403 MMP9 results 2020.xlsx", skip = 1)
#names(mmp_conc) <- c('record_id', 'delete', 'delete1', 'mrn', 'collection_dat', 'time_point', 'dilution', 'raw_result', 
 #                    'result', 'units')
#mmp_conc <- mmp_conc %>% dplyr::select(-delete, -delete1)


mmp_conc <- data.table::fread('16-1403 MMP9 results 2020.csv')


mmp_conc$time_point[which(mmp_conc$time_point == 'Follow-up')] <- '3 months'



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
  gtsave("MMP9_demographics.png", 
         vwidth = 1200, vheight = 800)



##Analysis & Plotting 

dat <- read.csv("Full_DKA_Data Set_06-26_final.csv")
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


qcr_df <- read.csv("qcr.csv")

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

dat <- read.csv("Full_DKA_Data Set_06-26_final.csv")
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


qcr_df <- read.csv("qcr.csv")

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

mmp9_dat <- read_excel("16-1403 MMP9 results 2020.xlsx", sheet = 1, skip = 1)
mmp9_dat <- mmp9_dat %>% select(SID, `Time Point`, `Actual Result`)
names(mmp9_dat) <- c("record_id", "time", "mmp9_actual")
mmp9_dat$time[mmp9_dat$time == "Follow-up"] <- "3 months"

data_set <- merge(data_set,mmp9_dat, all = T)
data_set$mmp9_egfr <- (data_set$mmp9_actual*100) / data_set$eGFR
data_set$fru[data_set$fru > 100] <- NA

cd93_dat <- read_excel("DKA ELISA cd93 stat.xlsx", sheet = 1)
names(cd93_dat) <- c("record_id", "time", "cd93")

data_set <- merge(data_set,cd93_dat, all = T)
data_set$cd93_egfr <- (data_set$cd93*100) / data_set$eGFR
data_set$fru_crea <- data_set$fru*100/data_set$crea


###### Plotting

library(ggplot2)
library(patchwork)
library(nlme)
library(emmeans)
library(dplyr)
library(ggsignif)

# Function to run mixed effects model and get pairwise comparisons
get_mixed_model_pvalues <- function(data, grouping_var) {
  formula_str <- paste0("mmp9_egfr ~ time * ", grouping_var)
  
  model <- lme(as.formula(formula_str), 
               random = ~1|record_id,  # Change to your subject ID variable
               data = data,
               na.action = na.omit)
  
  # Within-timepoint comparisons (compare groups at each time)
  emm_within <- emmeans(model, specs = as.formula(paste0("~ ", grouping_var, " | time")))
  within_time <- pairs(emm_within, adjust = "none")  # Changed from "tukey"
  
  # Across-time comparisons (compare times within each group)
  emm_across <- emmeans(model, specs = as.formula(paste0("~ time | ", grouping_var)))
  across_time <- pairs(emm_across, adjust = "none")  # Changed from "tukey"
  
  return(list(
    within_time = within_time,
    across_time = across_time,
    model = model
  ))
}

# Function to get significance stars
get_sig_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "ns")))
}

# Function to create plot with significance brackets
create_plot_with_brackets <- function(data, grouping_var, results, 
                                      fill_colors, legend_title) {
  
  # Get groups and create position mapping
  groups <- sort(unique(data[[grouping_var]]))
  n_groups <- length(groups)
  times <- c("0-8 hours", "3 months")
  
  # Base plot
  p <- ggplot(data, aes(x = time, y = mmp9_egfr, fill = .data[[grouping_var]])) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_classic() +
    scale_fill_manual(values = fill_colors) +
    labs(x = '', y = 'MMP9 (eGFR-adj)', fill = legend_title) +
    theme(text = element_text(size = 16))
  
  # Calculate y positions
  y_max <- max(data$mmp9_egfr, na.rm = TRUE)
  y_range <- diff(range(data$mmp9_egfr, na.rm = TRUE))
  y_step <- y_range * 0.12
  current_y <- y_max + y_step * 0.5
  
  # Convert to data frames and get proper structure
  within_df <- as.data.frame(results$within_time)
  across_df <- as.data.frame(results$across_time)
  
  # Debug output
  cat("\n--- Processing", legend_title, "---\n")
  cat("Groups:", paste(groups, collapse = ", "), "\n")
  cat("Row names of within comparisons:\n")
  print(head(rownames(within_df)))
  
  # For within-time comparisons, the emmeans object with | time creates sections
  # We need to manually track which rows belong to which timepoint
  # Typically the first n_comparisons rows are for first timepoint, next n for second
  n_groups_actual <- length(groups)
  n_comparisons_per_time <- n_groups_actual * (n_groups_actual - 1) / 2
  
  # Add within-time comparisons
  for (t_idx in 1:length(times)) {
    time_point <- times[t_idx]
    
    # Calculate which rows correspond to this timepoint
    start_row <- (t_idx - 1) * n_comparisons_per_time + 1
    end_row <- t_idx * n_comparisons_per_time
    
    if (end_row <= nrow(within_df)) {
      cat("\nTimepoint:", time_point, "- Rows", start_row, "to", end_row, "\n")
      
      for (row_idx in start_row:end_row) {
        p_val <- within_df$p.value[row_idx]
        sig <- get_sig_stars(p_val)
        contrast_text <- as.character(within_df$contrast[row_idx])
        
        cat("  Comparison:", contrast_text, "- p =", p_val, "- sig:", sig, "\n")
        
        if (sig != "ns") {
          # Parse contrast to get group names
          parts <- strsplit(contrast_text, " - ")[[1]]
          group1 <- trimws(parts[1])
          group2 <- trimws(parts[2])
          
          # Find group indices
          g1_idx <- which(groups == group1)
          g2_idx <- which(groups == group2)
          
          if (length(g1_idx) > 0 && length(g2_idx) > 0) {
            # Calculate x positions
            dodge_width <- 0.8
            x_base <- t_idx
            x1 <- x_base + (g1_idx - (n_groups + 1)/2) * (dodge_width/n_groups)
            x2 <- x_base + (g2_idx - (n_groups + 1)/2) * (dodge_width/n_groups)
            
            cat("    Adding bracket: x1 =", x1, ", x2 =", x2, ", y =", current_y, "\n")
            
            # Add bracket
            p <- p + 
              annotate("segment", x = x1, xend = x1, 
                       y = current_y, yend = current_y + y_step*0.3,
                       color = "black", linewidth = 0.8) +
              annotate("segment", x = x1, xend = x2, 
                       y = current_y + y_step*0.3, yend = current_y + y_step*0.3,
                       color = "black", linewidth = 0.8) +
              annotate("segment", x = x2, xend = x2, 
                       y = current_y + y_step*0.3, yend = current_y,
                       color = "black", linewidth = 0.8) +
              annotate("text", x = (x1 + x2)/2, y = current_y + y_step*0.55,
                       label = sig, size = 5, fontface = "bold")
            
            current_y <- current_y + y_step
          }
        }
      }
    }
  }
  
  # Add across-time comparisons
  # Each row in across_df corresponds to one group
  cat("\nAcross-time comparisons:\n")
  for (row_idx in 1:min(nrow(across_df), length(groups))) {
    group_name <- groups[row_idx]
    p_val <- across_df$p.value[row_idx]
    sig <- get_sig_stars(p_val)
    
    cat("Group:", group_name, "- p =", p_val, "- sig:", sig, "\n")
    
    if (sig != "ns") {
      # Calculate x positions for this group across both timepoints
      dodge_width <- 0.8
      x1 <- 1 + (row_idx - (n_groups + 1)/2) * (dodge_width/n_groups)
      x2 <- 2 + (row_idx - (n_groups + 1)/2) * (dodge_width/n_groups)
      
      cat("  Adding bracket: x1 =", x1, ", x2 =", x2, ", y =", current_y, "\n")
      
      # Add bracket
      p <- p + 
        annotate("segment", x = x1, xend = x1, 
                 y = current_y, yend = current_y + y_step*0.3,
                 color = "black", linewidth = 0.8) +
        annotate("segment", x = x1, xend = x2, 
                 y = current_y + y_step*0.3, yend = current_y + y_step*0.3,
                 color = "black", linewidth = 0.8) +
        annotate("segment", x = x2, xend = x2, 
                 y = current_y + y_step*0.3, yend = current_y,
                 color = "black", linewidth = 0.8) +
        annotate("text", x = (x1 + x2)/2, y = current_y + y_step*0.55,
                 label = sig, size = 5, fontface = "bold")
      
      current_y <- current_y + y_step
    }
  }
  
  # Expand y-axis to show brackets
  p <- p + coord_cartesian(ylim = c(0, current_y + y_step * 0.5))
  
  return(p)
}

# Filter data
data_filtered <- data_set %>% filter(time %in% c('0-8 hours', '3 months'))

# Run models with specific filtering for each variable
# Severity - remove Unknown and NAs
data_severity <- data_filtered %>% 
  filter(!is.na(mmp9_egfr), !is.na(severity), severity != "Unknown")
results_severity <- get_mixed_model_pvalues(data_severity, "severity")

# IV Insulin - remove NAs only
data_insulin <- data_filtered %>% 
  filter(!is.na(mmp9_egfr), !is.na(ivinsdurtile))
results_insulin <- get_mixed_model_pvalues(data_insulin, "ivinsdurtile")

# AKI - remove NAs only
data_aki <- data_filtered %>% 
  filter(!is.na(mmp9_egfr), !is.na(aki_0_24))
results_aki <- get_mixed_model_pvalues(data_aki, "aki_0_24")

# T1D Status - remove NAs only
data_t1d <- data_filtered %>% 
  filter(!is.na(mmp9_egfr), !is.na(t1d_status))
results_t1d <- get_mixed_model_pvalues(data_t1d, "t1d_status")

# Print results
cat("\n=== SEVERITY MODEL ===\n")
cat("\nWithin-Timepoint Comparisons:\n")
print(results_severity$within_time)
cat("\nAcross-Time Comparisons:\n")
print(results_severity$across_time)

cat("\n=== IV INSULIN MODEL ===\n")
cat("\nWithin-Timepoint Comparisons:\n")
print(results_insulin$within_time)
cat("\nAcross-Time Comparisons:\n")
print(results_insulin$across_time)

cat("\n=== AKI MODEL ===\n")
cat("\nWithin-Timepoint Comparisons:\n")
print(results_aki$within_time)
cat("\nAcross-Time Comparisons:\n")
print(results_aki$across_time)

cat("\n=== T1D STATUS MODEL ===\n")
cat("\nWithin-Timepoint Comparisons:\n")
print(results_t1d$within_time)
cat("\nAcross-Time Comparisons:\n")
print(results_t1d$across_time)

# Create plots with brackets
graph1 <- create_plot_with_brackets(data_filtered, "severity", results_severity,
                                    c("#fff9ec", "#fcb1a6", "#fb6376"), "Severity")

graph2 <- create_plot_with_brackets(data_filtered, "ivinsdurtile", results_insulin,
                                    c("#fff9ec", "#fcb1a6", "#fb6376"), "Days IV Insulin")

graph3 <- create_plot_with_brackets(data_filtered, "aki_0_24", results_aki,
                                    c('#4CAF50', '#F44336'), "AKI")

graph4 <- create_plot_with_brackets(data_filtered, "t1d_status", results_t1d,
                                    c('#4CAF50', '#F44336'), "T1D Status")

# Create scatter plots with regression lines
graph5 <- ggplot(data_filtered %>% filter(!is.na(scopeptin)), 
                 aes(x = scopeptin, y = mmp9_egfr, color = time)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = 'lm', se = TRUE, linewidth = 1) +
  theme_classic() +
  scale_color_manual(values = c("0-8 hours" = "#2E86AB", "3 months" = "#A23B72"),
                     name = "Time") +
  labs(x = 'Serum Copeptin', y = 'MMP9 (eGFR adj)') +
  theme(text = element_text(size = 16),
        legend.position = "bottom")

graph6 <- ggplot(data_filtered %>% filter(!is.na(sua)), 
                 aes(x = sua, y = mmp9_egfr, color = time)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = 'lm', se = TRUE, linewidth = 1) +
  theme_classic() +
  scale_color_manual(values = c("0-8 hours" = "#2E86AB", "3 months" = "#A23B72"),
                     name = "Time") +
  labs(x = 'Serum Uric Acid', y = 'MMP9 (eGFR adj)') +
  theme(text = element_text(size = 16),
        legend.position = "bottom")

# Combine all 6 plots (3 rows x 2 columns)
combined_all_plots <- (graph1 | graph2) / 
  (graph3 | graph4) /
  (graph5 | graph6)

combined_all_plots <- combined_all_plots + 
  plot_annotation(tag_levels = 'A', 
                  theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))

print(combined_all_plots)

ggsave("MMP9_all_plots.png", combined_all_plots,
       width = 14, height = 18, dpi = 300, units = "in")




##################################################################################### Analysis

# Create analysis dataset
analysis_data <- data_set %>%
  filter(time %in% c('0-8 hours', '3 months')) %>%
  filter(!is.na(mmp9_actual) & !is.na(record_id)) %>%
  mutate(
    time_binary = ifelse(time == '0-8 hours', 0, 1),  # 0 = DKA, 1 = 3 months
    time_factor = factor(time, levels = c('0-8 hours', '3 months'))
  )



# Mixed-effects model for time comparison
model1 <- lmer(mmp9_egfr ~ time_factor + (1|record_id), data = analysis_data)

# Model summary
print("Mixed-effects model summary:")
print(summary(model1))

# Get estimated marginal means (least squares means)
emm1 <- emmeans(model1, ~ time_factor)
print("\nEstimated marginal means:")
print(emm1)

as.data.frame(emm1)


# Pairwise comparisons
pairs1 <- pairs(emm1)
print("\nPairwise comparison:")
print(pairs1)

# Extract coefficients with confidence intervals
coef1 <- tidy(model1, effects = "fixed", conf.int = TRUE)
print("\nFixed effects coefficients:")
print(coef1)



# Filter for DKA time point only and check AKI variable
dka_data <- analysis_data %>% 
  filter(time == '0-8 hours', !is.na(aki_0_24))





model2 <- aov(mmp9_egfr ~ aki_0_24, data = dka_data)

print("Mixed-effects model summary for AKI:")
print(summary(model2))

# Estimated marginal means
emm2 <- emmeans(model2, ~ aki_0_24)
print("\nEstimated marginal means by AKI status:")
print(emm2)
as.data.frame(emm2)

# Pairwise comparisons
pairs2 <- pairs(emm2)
print("\nPairwise comparison:")
print(pairs2)

# Coefficients
coef2 <- tidy(model2, effects = "fixed", conf.int = TRUE)
print("\nFixed effects coefficients:")
print(coef2)


##### 3 months
dka_data <- analysis_data %>% 
  filter(time == '3 months', !is.na(aki_0_24))

cat("AKI variable levels at 0-8 hours:", unique(dka_data$aki_0_24), "\n")
cat("Sample sizes by AKI status:", table(dka_data$aki_0_24), "\n")


model2 <- aov(mmp9_egfr ~ aki_0_24, data = dka_data)

print("Mixed-effects model summary for AKI:")
print(summary(model2))

# Estimated marginal means
emm2 <- emmeans(model2, ~ aki_0_24)
print("\nEstimated marginal means by AKI status:")
print(emm2)
as.data.frame(emm2)


# Pairwise comparisons
pairs2 <- pairs(emm2)
print("\nPairwise comparison:")
print(pairs2)

# Coefficients
coef2 <- tidy(model2, effects = "fixed", conf.int = TRUE)
print("\nFixed effects coefficients:")
print(coef2)







##################################### Grouping analysis 


ngal_model <- lme(mmp9_egfr ~ time*scopeptin, data = data_set %>% filter(time %in% c('0-8 hours', '3 months')),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ time*sua, data = data_set %>% filter(time %in% c('0-8 hours', '3 months')),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)



#aki yes
ngal_model <- lme(mmp9_egfr ~ time*scopeptin, data = data_set %>% filter(aki_0_24 == 'Yes') %>% filter(time %in% c('0-8 hours', '3 months')),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ time*sua, data = data_set %>% filter(aki_0_24 == 'Yes') %>% filter(time %in% c('0-8 hours', '3 months')),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


#aki no

ngal_model <- lme(mmp9_egfr ~ time*scopeptin, data = data_set %>% filter(aki_0_24 == 'No') %>% filter(time %in% c('0-8 hours', '3 months')),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ time*sua, data = data_set %>% filter(aki_0_24 == 'No') %>% filter(time %in% c('0-8 hours', '3 months')),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


#three-way interaction

ngal_model <- lme(mmp9_egfr ~ time*scopeptin*aki_0_24, data = data_set,
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ time*sua*aki_0_24, data = data_set,
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


#three way but not interaction 
ngal_model <- lme(mmp9_egfr ~ time + scopeptin*aki_0_24, data = data_set,
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ time + sua*aki_0_24, data = data_set,
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)



ggplot(data_set %>% filter(time %in% c('0-8 hours', '3 months')), aes(x=mmp9_egfr, y= sua))+
  geom_point()+geom_smooth(method='lm') +facet_wrap(time ~ aki_0_24)


ggplot(data_set %>% filter(time %in% c('0-8 hours', '3 months')), aes(x=mmp9_egfr, y= scopeptin))+
  geom_point()+geom_smooth(method='lm') +facet_wrap(time ~ aki_0_24)







#split by time 

ngal_model <- lme(mmp9_egfr ~ scopeptin*aki_0_24, data = data_set %>% filter(time == '0-8 hours'),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ sua*aki_0_24, data = data_set %>% filter(time == '0-8 hours') ,
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


#3 months
ngal_model <- lme(mmp9_egfr ~ scopeptin*aki_0_24, data = data_set %>% filter(time == '3 months'),
                  random = ~1|record_id,
                  na.action = na.omit)

kable(anova(ngal_model))
kable(summary(ngal_model)$tTable)


ngal_model <- lme(mmp9_egfr ~ sua*aki_0_24, data = data_set %>% filter(time == '3 months') ,
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












### Analyses from plots 

get_mixed_model_pvalues <- function(data, grouping_var) {
  formula_str <- paste0("mmp9_egfr ~ time * ", grouping_var)
  
  model <- lme(as.formula(formula_str), 
               random = ~1|record_id,  # Change to your subject ID variable
               data = data,
               na.action = na.omit)
  
  # Within-timepoint comparisons (compare groups at each time)
  emm_within <- emmeans(model, specs = as.formula(paste0("~ ", grouping_var, " | time")))
  within_time <- pairs(emm_within, adjust = "none")  # Changed from "tukey"
  
  # Across-time comparisons (compare times within each group)
  emm_across <- emmeans(model, specs = as.formula(paste0("~ time | ", grouping_var)))
  across_time <- pairs(emm_across, adjust = "none")  # Changed from "tukey"
  
  return(list(
    within_time = within_time,
    across_time = across_time,
    model = model
  ))
}


data_filtered <- data_set %>% filter(time %in% c('0-8 hours', '3 months'))

# Run models with specific filtering for each variable
# Severity - remove Unknown and NAs
data_severity <- data_filtered %>% 
  filter(!is.na(mmp9_egfr), !is.na(severity), severity != "Unknown")
results_severity <- get_mixed_model_pvalues(data_severity, "severity")

# IV Insulin - remove NAs only
data_insulin <- data_filtered %>% 
  filter(!is.na(mmp9_egfr), !is.na(ivinsdurtile))
results_insulin <- get_mixed_model_pvalues(data_insulin, "ivinsdurtile")

# AKI - remove NAs only
data_aki <- data_filtered %>% 
  filter(!is.na(mmp9_egfr), !is.na(aki_0_24))
results_aki <- get_mixed_model_pvalues(data_aki, "aki_0_24")

# T1D Status - remove NAs only
data_t1d <- data_filtered %>% 
  filter(!is.na(mmp9_egfr), !is.na(t1d_status))
results_t1d <- get_mixed_model_pvalues(data_t1d, "t1d_status")
