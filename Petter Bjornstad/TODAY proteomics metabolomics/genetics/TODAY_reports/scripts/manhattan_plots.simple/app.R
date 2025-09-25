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

################################################################################
# Load data
progress             = readRDS(here("data/progress.rds"))
top_snps_description = fread(here("data", "top_snps_description.txt"), sep = "\t", header = TRUE, data.table = FALSE) %>%
  mutate(url_snp   = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr", CHROM, "%3A", POS)) %>%
  mutate(url_locus = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr", CHROM, "%3A", START, "-", END))

# Load Manhattan plot function
source(here("R/build_summary_table_shiny.R"))
source(here("R/create_manhattan_plot_summary.R"))


################################################################################
# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("QTL analysis"),
  
  tags$head(
    # Add sparkline dependencies
    #sparklineOutput("dummy"),  # This ensures sparkline is loaded
    tags$style(HTML("
      .col-sm-4 { width: 25% !important; }  /* Sidebar width */
      .col-sm-8 { width: 75% !important; }  /* Main content width */
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      h4("QTL summary statistics:"),
      p("Lead variant of all the QTL analyses"),
      
      # Genomic region input
      div(
        style = "margin-bottom: 15px;",
        h5("Genomic Region"),
        textInput(
          "genomic_region", 
          label = NULL,
          placeholder = "e.g., 1:1000-50000 or 22:1-1000000",
          value = ""
        ),
        helpText("Format: CHR:START-END (leave blank for all)")
      ),
      
      # P-value threshold slider
      div(
        style = "margin-bottom: 15px;",
        h5("P-value Threshold"),
        sliderInput(
          "p_threshold",
          label = NULL,
          min = -100, max = -1,
          value = log10(5e-8),
          step = 1,
          post = " (10^x)"
        ),
        div(
          style = "font-size: 12px; color: #666; margin-top: -10px;",
          textOutput("p_threshold_display")
        )
      ),
      # Apply filters button
      div(
        style = "margin-top: 15px;",
        actionButton("apply_filters", "Apply Filters", 
                     class = "btn-primary btn-block")
      ),
      
      # Clear filters button
      div(
        style = "margin-top: 10px;",
        actionButton("clear_filters", "Clear All Filters", 
                     class = "btn-secondary btn-block")
      )      
    ),
    
    mainPanel(
      # Add tabs to organize your outputs
        tabsetPanel(
          tabPanel(
            "Summary results",
            br(),
            div(
              style = "background: #f8f9fa; padding: 10px; margin-bottom: 15px; border-radius: 5px;",
              fluidRow(
                column(3, div(h6("Total signals:"    ), textOutput("total_snps"         , inline = TRUE))),
                column(3, div(h6("Region:"           ), textOutput("genomic_region"     , inline = TRUE))),
                column(3, div(h6("P-value threshold:"), textOutput("p_threshold_display", inline = TRUE)))
                )
              )
            ),      
          
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
          
          tabPanel(
            "Manhattan plot",
            br(),
            div(
              style = "margin: 20px;",
              h4("Manhattan plot showing all associations"),
              
              plotlyOutput("manhattan_plot", height = "600px")
            )
          ),      
          
          tabPanel(
            "Volcano plot",
            br(),
            div(
              style = "margin: 20px;",
              h4("Volcano plot showing all associations"),
              
              plotlyOutput("volcano_plot", height = "600px")
            )
          )
        )
      )
    )
  )


################################################################################
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Dummy sparkline to load dependencies
  output$dummy <- renderSparkline({
    sparkline(c(1,2,3))
  })
  
  # Reactive expression for P-value threshold display
  output$p_threshold_display <- renderText({
    threshold <- 10^input$p_threshold
    paste("P <", formatC(threshold, format = "e", digits = 1))
  })
  
  # Parse genomic region input
  parse_genomic_region <- function(region_string) {
    if (is.null(region_string) || region_string == "") return(NULL)
    
    # Remove spaces
    region_string <- gsub("\\s+", "", region_string)
    
    # Pattern: CHR:START-END
    pattern <- "^(\\d+|X|Y):([0-9,]+)-([0-9,]+)$"
    
    if (grepl(pattern, region_string)) {
      matches <- regmatches(region_string, regexec(pattern, region_string))[[1]]
      
      return(list(
        chr = matches[2],
        start = as.numeric(gsub(",", "", matches[3])),
        end = as.numeric(gsub(",", "", matches[4]))
      ))
    }
    
    return(NULL)
  }
  
  # Reactive data filtering
  filtered_data <- eventReactive(input$apply_filters, {
    data <- top_snps_description
    
    # Apply genomic region filter
    region <- parse_genomic_region(input$genomic_region)
    
    if (!is.null(region)) {
      data <- data[data$CHROM == region$chr & 
                     data$POS >= region$start & 
                     data$POS <= region$end, ]
    } 
    
    # Apply p-value threshold
    p_threshold <- 10^input$p_threshold
    data <- data[data$P <= p_threshold, ]
    
    return(data)
  }, ignoreNULL = FALSE)
  
  # Initialize with full data
  observeEvent(input$apply_filters, {}, ignoreInit = FALSE, once = TRUE)
  
  # Summary statistics
  output$total_snps <- renderText({
    format(nrow(top_snps_description), big.mark = ",")
  })
  
  output$current_region <- renderText({
    region <- parse_genomic_region(input$genomic_region)
    if (!is.null(region)) {
      paste0("Chr", region$chr, ":", 
             format(region$start, big.mark = ","), "-",
             format(region$end, big.mark = ","))
    } else if (input$chromosome != "") {
      paste0("Chr", input$chromosome)
    } else {
      "All chromosomes"
    }
  })
  
  # Quick threshold buttons
  observeEvent(input$genome_wide, {
    updateSliderInput(session, "p_threshold", value = log10(5e-8))
  })
  
  observeEvent(input$suggestive, {
    updateSliderInput(session, "p_threshold", value = log10(1e-5))
  })
  
  observeEvent(input$nominal, {
    updateSliderInput(session, "p_threshold", value = log10(0.05))
  })
  
  # Clear filters
  observeEvent(input$clear_filters, {
    updateTextInput(session, "genomic_region", value = "")
    updateSliderInput(session, "p_threshold", value = -2)
  })
  
  # Main table
  output$gwas_table <- renderReactable({
    if (input$apply_filters == 0) {
      data <- top_snps_description  
    } else {
      data <- filtered_data()
    }  
      create_gwas_table(data)
    })
  
  # Manhattan plot
  output$manhattan_plot <- renderPlotly({
    if (input$apply_filters == 0) {
      data <- top_snps_description  
    } else {
      data <- filtered_data()
    }  
    create_manhattan_plot(data)
  })
  
  # Volcano plot
  output$volcano_plot <- renderPlotly({
    if (input$apply_filters == 0) {
      data <- top_snps_description  
    } else {
      data <- filtered_data()
    }  
    create_volcano_plot(data)
  })
  
}

################################################################################
# Run the application 
shinyApp(ui = ui, server = server)
