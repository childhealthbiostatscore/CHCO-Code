#### DKA MMP9 Study 

library(dplyr)
library(stringr)
library(ggplot2)



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








#table1::table1( ~ sex_corrected + race + age_years + height + weight + bmi + t1d_status + systolic_bp + diastolic_bp + heartrate + severity + a1c_er + mmp9_0_8hr,
#               data = demographics)


#Plotting 

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


















