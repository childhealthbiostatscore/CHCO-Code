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


#Overall outlier evaluations



#create ui

ui <- fluidPage(
  shinythemes::shinytheme('journal'),
  h1('Outlier Identification App'),
  sidebarPanel(header = 'Outlier Investigations',
    sliderInput('integer', 'Minimum SD Away from Median:', min = 0, max = 10, step = 1, value = 2),
    selectInput('x', 'Select study: ', choices = study_names, 'All', multiple = TRUE),
    selectInput('MRNs', 'Which specific MRNs are you interested in?', choices = unique(harmonized_data$mrn), 
                multiple=TRUE)
    ),
  mainPanel(
    tabsetPanel(
      tabPanel(dataTableOutput('Outliers'), 'Outlier IDs')
      
    ),
    tabsetPanel(header = 'MRN Investigations',
      tabPanel(dataTableOutput('IDs'), 'ID/Study Summary')
    )
  )
  
)


server <- function(input, output){
  
  #Output to find outliers
  output$Outliers <- renderDataTable({
    
    find_outliers <- function(){
      
      
      
      
      
      
      
      
      
    }
    
    

    find_outliers()
    
  })
  
  
  
  #Output for IDs
  output$IDs <- renderDataTable({
    #Selected Individual analyses 
    
    get_ind_table <- function(){
      tmp_ids <- input$MRNs
      tmp_meds <- medications %>% 
        filter(mrn %in% tmp_ids)
      tmp_harm <- harmonized_data %>% 
        select(1:13, 16:19, 24, 1180, age, 24) %>% 
        filter(mrn %in% tmp_ids)
      tmp_results_df <- data.frame(MRN = tmp_ids, Study_IDs = NA, 
                                   Studies = NA, Age =NA, Sex = NA, Group = NA,
                                   Medications = NA, Race_Ethnicity = NA)
      
      for(i in c(1:nrow(tmp_results_df))){
        iter_id <- tmp_ids[i]
        #study IDs 
        tmp_harm_iter <- tmp_harm %>% filter(mrn == iter_id)
        
        tmp_results_df$Study_IDs[i] <- tmp_harm_iter %>% select(1:12) %>% 
          as.matrix() %>% as.vector() %>% unique() %>% paste(collapse = ', ')
        
        tmp_results_df$Studies[i] <- paste(unique(tmp_harm_iter$study), collapse=', ')
        
        age_sd <- sd(tmp_harm_iter$age, na.rm=T) %>% round(digits = 2)
        age_mean <- mean(tmp_harm_iter$age, na.rm=T) %>% round(digits = 2)
        tmp_results_df$Age[i] <- paste0('Mean: ', age_mean, ', SD:', age_sd,
                                        ', Values: ', paste(unique(tmp_harm_iter$age), collapse = ', '))
        
        tmp_results_df$Sex[i] <- paste0(names(table(tmp_harm_iter$sex)), ': ', table(tmp_harm_iter$sex),
                                        collapse = ', ')
        
        tmp_results_df$Group[i] <- paste0(names(table(tmp_harm_iter$group)), ': ', table(tmp_harm_iter$group),
                                          collapse = ', ')
        tmp_results_df$Race_Ethnicity[i] <- paste0(names(table(tmp_harm_iter$race_ethnicity)), 
                                                   ': ', table(tmp_harm_iter$race_ethnicity), collapse = ', ')
        
        
        
        #medication
        iter_med_vec <- tmp_meds %>% filter(mrn == iter_id) %>% str_detect('Yes')
        tmp_results_df$Medications[i] <- paste(names(tmp_meds)[iter_med_vec], collapse=', ')
        
        
      }
      #push the data.frame
      tmp_results_df$Studies <-  str_replace(tmp_results_df$Studies, pattern=', $', replacement ='')
      tmp_results_df$Studies <-  str_replace(tmp_results_df$Studies, pattern='NA,', replacement ='')
      tmp_results_df$Study_IDs <- str_replace(tmp_results_df$Study_IDs, pattern = 'NA, ', replacement = '')
      tmp_results_df$Study_IDs <- str_replace(tmp_results_df$Study_IDs, pattern=', $', replacement = '')
      tmp_results_df$Age <- str_replace(tmp_results_df$Age, pattern= 'NA, ', replacement ='')
      
      
      tmp_results_df
      
    }
    
    
    
    
    get_ind_table()
  })
  
}



shinyApp(ui = ui, server = server)








