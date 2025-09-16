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


demographics <- final_data %>% 
  dplyr::select(sex, race, age, height = ht, weight = wt, #bmi, 
                systolic_bp = sbp, diastolic_bp = dbp, heartrate = hr, )



























