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


#find variable

# Apply table() to each column and display results
lapply(final_data, table)

# Or to see them more clearly one by one:
for(col in names(final_data)) {
  cat("\nColumn:", col, "\n")
  print(table(final_data[[col]]))
}

# To find columns with specific counts (16, 15, 9):
target_counts <- c(16, 15, 9)
for(col in names(final_data)) {
  tab <- table(final_data[[col]])
  if(any(tab %in% target_counts)) {
    cat("\nColumn", col, "has target counts:\n")
    print(tab)
  }
}


final_data$severity <- ""
final_data$severity[((final_data$vbg_1 >= 7.3) & (final_data$bmp_1 >= 15))] <- "Mild"
final_data$severity[((final_data$vbg_1 < 7.3) | (final_data$bmp_1 < 15))] <- "Mild"
final_data$severity[((final_data$vbg_1 < 7.2) | (final_data$bmp_1 < 10))] <- "Moderate"
final_data$severity[((final_data$vbg_1 < 7.1) | (final_data$bmp_1 < 5))] <- "Severe"
final_data$severity[final_data$severity == ""] <- "Unknown"


#demographics
demographics <- final_data %>% 
  dplyr::select(sex, race, age, height = ht, weight = wt, #bmi, 
                systolic_bp = sbp, diastolic_bp = dbp, heartrate = hr, a1c_er, a1c_3mo,
                mmp9_0_8hr, mmp9_12_24hr, mmp9_3mo, t1d_status #DKA severity,
                ) %>%
  mutate(bmi = (weight / (height/100)^2),
         age_years = age/12)






#analysis dataframe

analysis_df <- final_data %>% 
  dplyr::select(record_id, mmp9_0_8hr, mmp9_3mo, scopeptin_0_8hr, scopeptin_3mo, 
                serum_uric_0_8hr = sua_0_8hr, serum_uric_3mo = sua_3mo,  
                urine_uric_0_8hr = uua_0_8_hr, urine_uric_3mo = uua_3mo, 
                ivinsulin_start, ivinsulinstop)


























