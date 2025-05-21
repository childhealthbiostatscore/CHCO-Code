library(fastmap)
library(dplyr)
library(stringr)
library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(tidyr)
library(purrr)




setwd('C:/Users/netio/Documents/Harmonized_data/')

#Function to calculate most often answer
calculate_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

determine_any_meds <- function(x){
  str_detect(x, pattern='Yes') %>% sum() > 0 
}




dir.dat <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/"
harmonized_data <- read.csv(fs::path(dir.dat,"Data Harmonization","Data Clean","harmonized_dataset.csv"),na="")
dat <- harmonized_data %>%
  arrange(screen_date) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(mrn, visit))





#medication dataset

medications <- readxl::read_xlsx("Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>% select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), pattern = '_1', replacement = '')
medications$Any <- apply(medications, 1, determine_any_meds)
medications$Any <- ifelse(medications$Any == TRUE, 'Yes', 'No')

med_list <- names(medications)[2:16] %>% str_replace(pattern='_1', replacement = '')
med_list <- c('Any', med_list, 'Ignore')

med_desc <- readxl::read_xlsx('Biopsies_w_mrn_Oct3.xlsx', sheet = 2, col_names=F)
med_desc <- med_desc[seq(from = 1, to = 57, by =4),]
names(med_desc) <- c('Medication', 'Description')
med_desc$Medication <- str_replace(med_desc$Medication, pattern='_1', replacement = '')

med_desc_2 <- data.frame(Medication = c('Any', 'Ignore'), 
                         Description = c('If participant is taking any relevant meds during the 1st kidney biopsy.',
                                         'Medication use is not being shown in data'))

med_desc <- rbind(med_desc_2, med_desc)







#Obtain data and study names 
harmonized_data <- data.table::fread("harmonized_dataset.csv")
harmonized_data$sex[which(harmonized_data$sex == '')] <- 'Other'
study_names <- unique(harmonized_data$study)
study_names <- c(study_names, 'All')

study_descriptions <- readxl::read_xlsx('study_descriptions.xlsx')






#create ui

ui <- fluidPage(
  shinythemes::shinytheme('journal'),
  h1('Outlier Identification App'),
  sidebarPanel(
    sliderInput('integer', 'Minimum SD Away from Median:', min = 0, max = 10, step = 1, value = 2),
    selectInput('x', 'Select study: ', choices = study_names, 'All', multiple = TRUE),
    selectInput('MRNs', 'Which specific MRNs are you interested in?', choices = unique(harmonized_data$mrn), 
                multiple=TRUE)
    ),
  mainPanel(
    tabsetPanel(
      
    ),
    tabsetPanel(
      
    )
  )
  
)


server <- function(input, output){
  
  
  
  
  
}



shinyApp(ui = ui, server = server)








