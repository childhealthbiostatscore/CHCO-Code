# Load required libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(reactable))
suppressPackageStartupMessages(library(sparkline))

# Load data
progress             = readRDS(here("data/progress.rds"))
top_snps_description = fread(here("data", "top_snps_description.txt"), sep = "\t", header = TRUE, data.table = FALSE)

# Load Manhattan plot function
source(here("R/manhattan_plot.R"))
source(here("R/build_summary_table_shiny.R"))


# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Dependent Dropdown Menus"),
  
  tags$head(
    # Add sparkline dependencies
    sparklineOutput("dummy")  # This ensures sparkline is loaded
  ),
  
  sidebarLayout(
    sidebarPanel(
      # First dropdown: Select vector name
      selectInput("qtl_analysis", 
                  "Select a QTL analysis:",
                  choices = names(progress),
                  selected = names(progress)[1]),
      
      # Second dropdown: Select element from chosen vector
      selectInput("phenotype", 
                  "Select a phenotype:",
                  choices = NULL),  # Will be populated dynamically
      
      br(),
      h4("Instructions:"),
      p("1. First, choose a QTL analysis from the dropdown above"),
      p("2. Then select a specific phenotype"),
      p("3. Your selections will be displayed in the main panel")
    ),
    
    mainPanel(
      # Add tabs to organize your outputs
        tabsetPanel(
        tabPanel("Summary statistics report", 
                 verbatimTextOutput("selection_output")),
      
        tabPanel(
          "GWAS Results",
          br(),
          div(
            style = "margin: 20px;",
            h4("Top GWAS Association Results"),
            p("Summary statistics showing the most significant associations from metabolomics and proteomics GWAS analyses."),
            reactableOutput("gwas_table", height = "600px")
          )
        ),
        
        tabPanel("Manhattan Plot", 
                 plotlyOutput("manhattan_plot", height = "600px")),
        
        tabPanel("Summary statistics table", 
                 tableOutput("summary_table"))
        )
      )
    )
  )


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Dummy sparkline to load dependencies
  output$dummy <- renderSparkline({
    sparkline(c(1,2,3))
  })
  
  # Update the second dropdown based on first dropdown selection
  observeEvent(input$qtl_analysis, {
    # Get the vector corresponding to the selected name
    phenotypes <- progress[[input$qtl_analysis]][["summary_stats"]]
    
    # Update the choices in the second dropdown
    updateSelectInput(session, 
                      "phenotype",
                      choices = phenotypes,
                      selected = phenotypes[1])  # Select first element by default
  })
  
  # Display the selected values
  output$selection_output <- renderText({
    if (!is.null(input$phenotype)) {
      paste("Selected Category:", input$qtl_analysis,
            "\nSelected Element:", input$phenotype,
            "\nIndex in Vector:", which(progress[[input$qtl_analysis]][["summary_stats"]] == input$phenotype))
    }
  })
  
  output$gwas_table <- renderReactable({
    create_gwas_table(top_snps_description)
    })

  sumstats      = reactive({
    req(input$qtl_analysis, input$phenotype)
    sumstats_file = paste(progress[[input$qtl_analysis]][["folder_path"]], "summary_statistics", paste(input$phenotype, "summary_stats.rds", sep = "_"), sep = "/")
    
    readRDS(sumstats_file)
  })
  
  output$summary_table <- renderTable({
    data <- sumstats() %>% 
      arrange(P) %>% 
      select(CHROM, POS, REF, ALT, BETA, SE, P) %>%
      mutate(P = format(P, scientific = TRUE, digits = 3))
    
    head(data, 100)  # Show only first 100 rows
  }, 
  striped  = TRUE,           # Alternating row colors
  hover    = TRUE,             # Hover effects
  bordered = TRUE,          # Table borders
  spacing  = 'xs',           # Compact spacing
  width    = '100%',           # Full width
  align    = 'c')              # Center align
  
  output$manhattan_plot <- renderPlotly({
    req(input$qtl_analysis, input$phenotype)
    
    data <- sumstats()
    
    create_manhattan_plot(
      data = data,
      significance_threshold = 5e-8,
      suggestive_threshold = 1e-5,
      point_size = 3,
      title = paste(input$phenotype)
    )
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
