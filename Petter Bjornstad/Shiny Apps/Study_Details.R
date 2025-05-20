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


#dir.dat <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/"

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
  h1('Demographics of Studies App'),
  #add to sidebar: description of studies
  sidebarPanel(
    helpText('Choose a study of interest from the menu'),
    selectInput('x', 'Select study: ', choices = study_names, 'All', multiple = TRUE),
    htmlOutput('Study Description title'),
    htmlOutput('Study Descriptions'),
    selectInput('medication', 'Select medication: ', choices = med_list, 'Ignore'),
    textOutput('Medication Description')
  ),
  #might split these into demographics and more data driven. Or add blood/proteomics data
  mainPanel(tabsetPanel(header = 'Study Demographics',
                        tabPanel('Sex', plotly::plotlyOutput('Sex_Med_Plot')),
                        tabPanel('Age', plotOutput('Age_Med_Plot')), 
                        tabPanel("Height", plotly::plotlyOutput('Height_Med_Plot')),
                        tabPanel("Weight", plotly::plotlyOutput('Weight_Med_Plot')),
                        tabPanel('BMI', plotOutput('BMI_Med_Plot')),
                        tabPanel('Diabetes Duration', plotOutput('Diabetes_Med_Plot')), 
                        tabPanel('Race and Ethnicity', plotly::plotlyOutput('RaceEthnicity_Plot'))
  ),
  tabsetPanel(header = 'Research Data',
              tabPanel('Procedures', plotly::plotlyOutput("Procedure_Plot")),
              tabPanel('Participant Groups', plotly::plotlyOutput('Study_Groups')),
              tabPanel('albumin_creatinine ratio', plotly::plotlyOutput('Acr_u_Plot')),
              tabPanel('creatine', plotly::plotlyOutput('Crea_Plot')),
              tabPanel('cystatin', plotly::plotlyOutput('Cys_Plot')),
              tabPanel('HBA1C', plotly::plotlyOutput('HPA1C')),
              tabPanel('GFR', plotly::plotlyOutput('Gfr_Plot')),
              tabPanel('ERPF', plotly::plotlyOutput('Erpf_Plot')),
              tabPanel("Trigylcerides", plotly::plotlyOutput('Trigly_Plot')),
              tabPanel('HDL', plotly::plotlyOutput('HDL_Plot')),
              tabPanel('LDL', plotly::plotlyOutput('LDL_Plot'))
              
  )
  )
)


#server side functions
server <- function(input, output){
  
  #Creating the html title 
  output$`Study Description title` <- renderUI({
    HTML(paste0('You have selected ','<b>', input$x, '</b>'))
  })
  
  #Creating the text descriptions
  
  output$`Study Descriptions` <- renderUI({
    HTML(paste0( '<b>', study_descriptions$Study[which(study_descriptions$Study %in% input$x)], '</b>: ',
                 study_descriptions$Description[which(study_descriptions$Study %in% input$x)],
                 collapse = '<br/><br/>'))
  })
  
  
  
  
  #Plotting race/ethnicity. Need to combine the multiple repeat options (like for Hispanic/latino)
  plot_race <- function(){
    if('All' %in% input$x){
      race_ethn <- harmonized_data %>% group_by(record_id) %>% 
        summarize(Race = calculate_mode(race_ethnicity))
      race_ethn_df <- table(race_ethn$Race) %>% t() %>% as.data.frame()
      race_ethn_df$Var2 <- as.character(race_ethn_df$Var2)
      plot_ly(race_ethn_df, labels = ~Var2, values = race_ethn_df$Freq, type = 'pie') %>% 
        layout(title = paste0('Race/Ethnicity Data in All'), 
               xaxis = list(showgrid = F, zeroline = F, showticklabels = F),
               yaxis = list(showgrid = F, zeroline = F, showticklabels = F))
    }else{
      tmp_data <- harmonized_data %>% 
        filter(study %in% input$x)
      tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
      
      race_ethn <- tmp_data %>% group_by(record_id) %>% 
        summarize(Race = calculate_mode(race_ethnicity))
      race_ethn_df <- table(race_ethn$Race) %>% t() %>% as.data.frame()
      race_ethn_df$Var2 <- as.character(race_ethn_df$Var2)
      plot_ly(race_ethn_df, labels = ~Var2, values = race_ethn_df$Freq, type = 'pie') %>% 
        layout(title = paste0('Race/Ethnicity in ', paste(input$x, collapse= ', ')), 
               xaxis = list(showgrid = F, zeroline = F, showticklabels = F),
               yaxis = list(showgrid = F, zeroline = F, showticklabels = F))
      
    } 
  }
  
  output$RaceEthnicity_Plot <- plotly::renderPlotly({
    plot_race()
  })     
  
  
  #Plotting research groups 
  plot_group <- function(){
    if('All' %in% input$x){
      group_df <- harmonized_data %>% group_by(record_id) %>% 
        summarize(Group = calculate_mode(group))
      group_df <- table(group_df$Group) %>% t() %>% as.data.frame()
      group_df$Var2 <- as.character(group_df$Var2)
      plot_ly(group_df, labels = ~Var2, values = group_df$Freq, type = 'pie') %>% 
        layout(title = paste0('Research Groups in All'), 
               xaxis = list(showgrid = F, zeroline = F, showticklabels = F),
               yaxis = list(showgrid = F, zeroline = F, showticklabels = F))
    }else{
      tmp_data <- harmonized_data %>% 
        filter(study %in% input$x)
      tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
      
      group_df <- tmp_data %>% group_by(record_id) %>% 
        summarize(Group = calculate_mode(group))
      group_df <- table(group_df$Group) %>% t() %>% as.data.frame()
      group_df$Var2 <- as.character(group_df$Var2)
      plot_ly(group_df, labels = ~Var2, values = group_df$Freq, type = 'pie') %>% 
        layout(title = paste0('Research Groups in ', paste(input$x, collapse= ', ')), 
               xaxis = list(showgrid = F, zeroline = F, showticklabels = F),
               yaxis = list(showgrid = F, zeroline = F, showticklabels = F))
      
    } 
  }
  
  output$Study_Groups <- plotly::renderPlotly({
    plot_group()
  })     
  
  #plotting list of procedures. Need to make sure this is correct, and the information we want to share. 
  plot_procedure <- function(){
    if('All' %in% input$x){
      procedure_df <- table(harmonized_data$procedure) %>% t() %>% as.data.frame()
      plotly::plot_ly(x=procedure_df$Var2, y=procedure_df$Freq,
                      names='Data Available by Procedures')%>% 
        layout(title = paste0('Procedure Data Available in All'))
    }else{
      tmp_data <- harmonized_data %>% 
        filter(study %in% input$x)
      tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
      
      procedure_df <- table(tmp_data$procedure) %>% t() %>% as.data.frame()
      plotly::plot_ly(x=procedure_df$Var2, y=procedure_df$Freq,
                      names='Data Available by Procedures') %>% 
        layout(title = paste0('Data Available by Procedures in ', paste(input$x, collapse= ', ')))
    }
  }
  
  output$Procedure_Plot <- plotly::renderPlotly({
    plot_procedure()
  })
  
  
  
  
  
  
  #####MEDICATION SPECIFICATION PLOTS ####
  
  #sex
  
  plot_sex_med <- function(){
    if(input$medication == 'Ignore'){
      if('All' %in% input$x){
        age_sex <- harmonized_data %>% group_by(record_id) %>% 
          summarize(Age = mean(age, na.rm=T),
                    Sex = calculate_mode(sex))
        sex_df <- table(age_sex$Sex) %>% t() %>% as.data.frame()
        sex_df$Var2 <- as.character(sex_df$Var2)
        plot_ly(sex_df, labels = ~Var2, values = sex_df$Freq, type = 'pie') %>% 
          layout(title = paste0('Sex Distribution in All'), 
                 xaxis = list(showgrid = F, zeroline = F, showticklabels = F),
                 yaxis = list(showgrid = F, zeroline = F, showticklabels = F))
      }else{
        tmp_data <- harmonized_data %>% 
          filter(study %in% input$x)
        tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
        
        age_sex <- tmp_data %>% group_by(record_id) %>% 
          summarize(Age = mean(age, na.rm=T),
                    Sex = calculate_mode(sex))
        sex_df <- table(age_sex$Sex) %>% t() %>% as.data.frame()
        sex_df$Var2 <- as.character(sex_df$Var2)
        plot_ly(sex_df, labels = ~Var2, values = sex_df$Freq, type = 'pie') %>% 
          layout(title =  paste0('Sex Distribution in ', paste(input$x, collapse= ', ')), 
                 xaxis = list(showgrid = F, zeroline = F, showticklabels = F),
                 yaxis = list(showgrid = F, zeroline = F, showticklabels = F))
        
      } 
    } else if(input$medication != 'Ignore' & 'All' %in% input$x){
      age_sex <- harmonized_data %>% group_by(mrn) %>% 
        summarize(Age = mean(age, na.rm=T),
                  Sex = calculate_mode(sex))
      
      age_sex$mrn <- as.character(age_sex$mrn)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      age_sex <- age_sex %>% left_join(meds_df, by='mrn')
      
      age_sex <- age_sex %>% group_by(medication, Sex) %>% summarize(Count = n())
      age_sex <- age_sex %>% spread(key = Sex, value = Count, fill = 0)
      
      
      plot_ly(age_sex, x=~medication, y=~Female,
              type='bar', name = 'Female') %>% 
        add_trace(y = ~Male, name = 'Male') %>% 
        add_trace(y = ~Other, name = 'Other') %>% 
        layout(yaxis = list(title = paste0('Sex Counts in All')), 
               xaxis = list(title = paste0('Use of ', input$medication)),
               barmode='stack')
      
    }else if(input$medication != 'Ignore'){
      
      tmp_data <- harmonized_data %>% 
        filter(study %in% input$x)
      tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
      
      age_sex <- tmp_data %>% group_by(mrn) %>% 
        summarize(Age = mean(age, na.rm=T),
                  Sex = calculate_mode(sex))
      
      age_sex$mrn <- as.character(age_sex$mrn)
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      age_sex <- age_sex %>% left_join(meds_df, by='mrn')
      
      age_sex <- age_sex %>% group_by(medication, Sex) %>% summarize(Count = n())
      age_sex <- age_sex %>% spread(key = Sex, value = Count, fill = 0)
      
      
      figure <- plot_ly(age_sex, x=~medication, y=~Female,
                        type='bar', name = 'Female') %>% 
        add_trace(y = ~Male, name = 'Male')
      
      if('Other' %in% names(age_sex)){
        figure <- figure %>% add_trace(y = ~Other, name = 'Other') 
      }
      figure %>%
        layout(yaxis = list(title = paste0('Sex Counts in ', paste(input$x, collapse= ', '))), 
               xaxis = list(title = paste0('Use of ', input$medication)),
               barmode='stack')
      
      
      
    }  
  }
  
  output$Sex_Med_Plot <- renderPlotly({
    plot_sex_med()
  })
  
  
  
  
  
  
  #age 
  
  plot_age_med <- function(){
    if(input$medication == 'Ignore'){
      if('All' %in% input$x){
        age_sex <- harmonized_data %>% group_by(record_id) %>% 
          summarize(Age = mean(age, na.rm=T),
                    Sex = calculate_mode(sex))
        ggplot(age_sex, aes(x=Age))+geom_density(alpha=0.2)+
          scale_x_continuous(limits = c(8, 80))+
          theme_classic()+labs(title = paste0('Age Distribution in All'))
      }else{
        tmp_data <- harmonized_data %>% 
          filter(study %in% input$x)
        tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
        
        age_sex <- tmp_data %>% group_by(record_id) %>% 
          summarize(Age = mean(age, na.rm=T),
                    Sex = calculate_mode(sex))
        ggplot(age_sex, aes(x=Age))+geom_density(alpha=0.2)+
          scale_x_continuous(limits = c(8, 80))+
          theme_classic()+labs(title = paste0('Age Distribution in ', paste(input$x, collapse= ', ')))
      }
    } else if(input$medication != 'Ignore' & 'All' %in% input$x){
      age_sex <- harmonized_data %>% group_by(mrn) %>% 
        summarize(Age = mean(age, na.rm=T),
                  Sex = calculate_mode(sex))
      
      age_sex$mrn <- as.character(age_sex$mrn)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      age_sex <- age_sex %>% left_join(meds_df, by='mrn')
      
      ggplot(age_sex, aes(x=Age, color=medication))+geom_density(alpha=0.2)+
        scale_x_continuous(limits = c(8, 80))+
        theme_classic()+labs(title = paste0('Age Distribution in All'))+
        labs(color = paste0(input$medication, ' use'))
      
    }else if(input$medication != 'Ignore'){
      
      tmp_data <- harmonized_data %>% 
        filter(study %in% input$x)
      tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
      
      age_sex <- tmp_data %>% group_by(mrn) %>% 
        summarize(Age = mean(age, na.rm=T),
                  Sex = calculate_mode(sex))
      
      age_sex$mrn <- as.character(age_sex$mrn)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      age_sex <- age_sex %>% left_join(meds_df, by='mrn')
      
      ggplot(age_sex, aes(x=Age, color=medication))+geom_density(alpha=0.2)+
        scale_x_continuous(limits = c(8, 80))+
        theme_classic()+labs(title = paste0('Age Distributions in ', paste(input$x, collapse= ', ')))+
        labs(color = paste0(input$medication, ' use')) 
      
      
      
      
    }  
  }
  
  output$Age_Med_Plot <- renderPlot({
    plot_age_med()
  })
  
  #BMI
  
  plot_bmi_med <- function(){
    if(input$medication == 'Ignore'){
      if('All' %in% input$x){
        bmi <- harmonized_data %>% group_by(record_id) %>% 
          summarize(BMI = mean(bmi, na.rm=T))
        ggplot(bmi, aes(x=BMI))+geom_density(alpha=0.2)+
          theme_classic()+labs(title = paste0('BMI Distribution in All'))
      }else{
        tmp_data <- harmonized_data %>% 
          filter(study %in% input$x)
        tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
        
        bmi <- tmp_data %>% group_by(record_id) %>% 
          summarize(BMI = mean(bmi, na.rm=T))
        ggplot(bmi, aes(x=BMI))+geom_density(alpha=0.2)+
          theme_classic()+labs(title =  paste0('BMI Distribution in ', paste(input$x, collapse= ', ')))
      }
    } else if(input$medication != 'Ignore' & 'All' %in% input$x){
      
      bmi <- harmonized_data %>% group_by(mrn) %>% 
        summarize(BMI = mean(bmi, na.rm=T))
      bmi$mrn <- as.character(bmi$mrn)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      
      bmi <- bmi %>% left_join(meds_df, by='mrn')
      
      ggplot(bmi, aes(x=medication, y=BMI))+geom_violin(aes(fill=medication))+
        geom_boxplot(width=0.2, color='grey', alpha=0.2)+
        theme_classic()+
        labs(title =  paste0('BMI Distribution in ', paste(input$x, collapse= ', ')),
             x=paste('Use of ', input$medication))
      
    }else if(input$medication != 'Ignore'){
      
      tmp_data <- harmonized_data %>% 
        filter(study %in% input$x)
      tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      bmi <- tmp_data %>% group_by(mrn) %>% 
        summarize(BMI = mean(bmi, na.rm=T))
      bmi$mrn <- as.character(bmi$mrn)
      
      bmi <- bmi %>% left_join(meds_df, by='mrn')
      
      ggplot(bmi, aes(x=medication, y=BMI))+geom_violin(aes(fill=medication))+
        geom_boxplot(width=0.2, color='grey', alpha=0.2)+
        theme_classic()+
        labs(title =  paste0('BMI Distribution in ', paste(input$x, collapse= ', ')),
             x=paste('Use of ', input$medication))
      
      
      
    }  
  }
  
  output$BMI_Med_Plot <- renderPlot({
    plot_bmi_med()
  })
  
  
  
  
  
  #Diabetes duration 
  plot_diab_med <- function(){
    if(input$medication == 'Ignore'){
      if('All' %in% input$x){
        diabetes <- harmonized_data %>% group_by(record_id) %>% 
          summarize(Dia_Dur = mean(diabetes_duration, na.rm=T))
        ggplot(diabetes, aes(x=Dia_Dur))+geom_density(alpha=0.2)+
          theme_classic()+labs(title = paste0('Length of Time of Diabetes Diagnosis in All'),
                               x='Duration of Diabetes Diagnosis')
      }else{
        tmp_data <- harmonized_data %>% 
          filter(study %in% input$x)
        tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
        
        diabetes <- tmp_data %>% group_by(record_id) %>% 
          summarize(Dia_Dur = mean(diabetes_duration, na.rm=T))
        ggplot(diabetes, aes(x=Dia_Dur))+geom_density(alpha=0.2)+
          theme_classic()+labs(title =  paste0('Length of Time of Diabetes Diagnosis in ', paste(input$x, collapse= ', ')),
                               x='Duration of Diabetes Diagnosis')
      }
    } else if(input$medication != 'Ignore' & 'All' %in% input$x){
      
      diab <- harmonized_data %>% group_by(mrn) %>% 
        summarize(Diab_Dur = mean(diabetes_duration, na.rm=T))
      diab$mrn <- as.character(diab$mrn)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      
      diab <- diab %>% left_join(meds_df, by='mrn')
      
      ggplot(diab, aes(x=medication, y=Diab_Dur))+geom_violin(aes(fill=medication))+
        geom_boxplot(width=0.2, color='grey', alpha=0.2)+
        theme_classic()+
        labs(title =  paste0('Diabetes Distribution in ', paste(input$x, collapse= ', ')),
             x=paste('Use of ', input$medication), y='Diabetes Duration')
      
    }else if(input$medication != 'Ignore'){
      
      tmp_data <- harmonized_data %>% 
        filter(study %in% input$x)
      tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      diab <- tmp_data %>% group_by(mrn) %>% 
        summarize(Diab_Dur = mean(diabetes_duration, na.rm=T))
      diab$mrn <- as.character(diab$mrn)
      
      diab <- diab %>% left_join(meds_df, by='mrn')
      
      ggplot(diab, aes(x=medication, y=Diab_Dur))+geom_violin(aes(fill=medication))+
        geom_boxplot(width=0.2, color='grey', alpha=0.2)+
        theme_classic()+
        labs(title =  paste0('Diabetes Distribution in ', paste(input$x, collapse= ', ')),
             x=paste('Use of ', input$medication), y='Diabetes Duration')
      
      
      
    }  
  }
  
  output$Diabetes_Med_Plot <- renderPlot({
    plot_diab_med()
  })
  
  
  
  
  #Make a generalized script to accept names of variables, and output boxplots, etc. 
  
  generalized_test_plot <- function(variable){
    if(input$medication == 'Ignore'){
      if('All' %in% input$x){
        
        index <- which(names(harmonized_data) == variable)
        tmp_data <- harmonized_data %>% select(record_id, mrn, index)
        names(tmp_data)[3] <- 'Variable'
        
        
        tmp_data <- tmp_data %>% group_by(record_id) %>% 
          summarize(Variable = mean(Variable, na.rm=T))
        
        ggplot(tmp_data, aes(x=Variable))+geom_density(alpha=0.2)+
          theme_classic()+labs(title = paste0(variable, ' Distribution in All'),
                               x=variable)
      }else{
        tmp_data <- harmonized_data %>% 
          filter(study %in% input$x)
        tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
        
        
        index <- which(names(tmp_data) == variable)
        tmp_data <- tmp_data %>% select(record_id, mrn, index)
        names(tmp_data)[3] <- 'Variable'
        tmp_data <- tmp_data %>% group_by(record_id) %>% 
          summarize(Variable = mean(Variable, na.rm=T))
        
        ggplot(tmp_data, aes(x=Variable))+geom_density(alpha=0.2)+
          theme_classic()+labs(title =  paste0(variable, ' Distribution in ', 
                                               paste(input$x, collapse= ', ')),
                               x = variable) %>%
          ggplotly()
      }
    } else if(input$medication != 'Ignore' & 'All' %in% input$x){
      
      index <- which(names(harmonized_data) == variable)
      tmp_data <- harmonized_data %>% select(record_id, mrn, index)
      names(tmp_data)[3] <- 'Variable'
      
      
      tmp_data <- tmp_data %>% group_by(mrn) %>% 
        summarize(Variable = mean(Variable, na.rm=T))
      tmp_data$mrn <- as.character(tmp_data$mrn)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      
      tmp_data <- tmp_data %>% left_join(meds_df, by='mrn')
      
      
      plot_ly(tmp_data, x = ~medication, y = ~Variable, split = ~medication, 
              type = 'violin', box = list(visible = TRUE),
              meanline = list(visible = TRUE)) %>% 
        layout(xaxis = list(title = paste('Use of ', input$medication)), 
               yaxis = list(title = variable),
               title =  paste0(variable, ' Distribution in ', paste(input$x, collapse= ', ')))
      
    }else if(input$medication != 'Ignore'){
      
      tmp_data <- harmonized_data %>% 
        filter(study %in% input$x)
      tmp_data <- harmonized_data %>% filter(record_id %in% tmp_data$record_id)
      
      meds_df <- medications %>% select(mrn, which(names(medications) == input$medication))
      names(meds_df) <- c('mrn', 'medication')
      
      index <- which(names(tmp_data) == variable)
      tmp_data <- tmp_data %>% select(record_id, mrn, index)
      names(tmp_data)[3] <- 'Variable'
      
      tmp_data <- tmp_data %>% group_by(mrn) %>% 
        summarize(Variable = mean(Variable, na.rm=T))
      
      
      tmp_data$mrn <- as.character(tmp_data$mrn)
      tmp_data <- tmp_data %>% left_join(meds_df, by='mrn')
      
      
      plot_ly(tmp_data, x = ~medication, y = ~Variable, split = ~medication, 
              type = 'violin', box = list(visible = TRUE),
              meanline = list(visible = TRUE)) %>% 
        layout(xaxis = list(title = paste('Use of ', input$medication)), 
               yaxis = list(title = variable),
               title =  paste0(variable, ' Distribution in ', paste(input$x, collapse= ', ')))
      
      
    }  
    
    
    
    
    
    
    
    
  }
  
  
  
  output$Height_Med_Plot <- renderPlotly({
    generalized_test_plot('height')
  })
  output$Weight_Med_Plot <- renderPlotly({
    generalized_test_plot('weight')
  })
  output$Acr_u_Plot <- renderPlotly({
    generalized_test_plot('acr_u')
  })
  output$Crea_Plot <- renderPlotly({
    generalized_test_plot('creatinine_s')
  })
  output$Cys_Plot <- renderPlotly({
    generalized_test_plot('cystatin_c_s')
  })
  output$HPA1C <- renderPlotly({
    generalized_test_plot('hba1c')
  })
  output$Gfr_Plot <- renderPlotly({
    generalized_test_plot('gfr_raw_plasma')
  })
  output$Erpf_Plot <- renderPlotly({
    generalized_test_plot('erpf_raw_plasma')
  })
  output$Trigly_Plot <- renderPlotly({
    generalized_test_plot('triglycerides')
  })
  output$HDL_Plot <-renderPlotly({
    generalized_test_plot('hdl')
  })
  output$LDL_Plot <-renderPlotly({
    generalized_test_plot('ldl')
  })
  
  
  
  
  
  
  
  
  
  
}



#Running the app 
shinyApp(ui = ui, server = server)

