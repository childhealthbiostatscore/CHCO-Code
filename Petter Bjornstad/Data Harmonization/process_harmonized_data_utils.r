#
# This Shiny web application provides a dashboard for exploring the
# 'harmonized_dataset.csv' file with multiple interactive filters.
#
# To run this app, place this file and your `harmonized_dataset.csv` file
# in the correct directory, then run the following in your R console:
# shiny::runApp()
#

# --- 1. Load Required Libraries ---
# If these packages are not installed, uncomment and run the following lines:
# install.packages(c("shiny", "dplyr", "ggplot2", "DT", "readr", "tidyr", "shinyjs", "shinyauthr"))
library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(readr)
library(tidyr) # Required for pivot functions in descriptive table
library(shinyjs) # Required for UI toggling and password-only login
library(shinyauthr) # Required for authentication
library(RColorBrewer) # Added for plot coloring palettes

# --- 1.5 Define Authentication Credentials ---
# IMPORTANT: CHANGE THIS PASSWORD TO A STRONG ONE FOR PRODUCTION!
user_base <- tibble::tibble(
  user = c("harmonized", "temp"), # Fixed dummy username for password-only access
  password = c("k!dn3y", "k1dn3y_t3mp"), # <-- UPDATED: Use this password to log in
  permissions = c("standard", "standard"),
  name = c("Data Access", "Temporary Access")
)
# --- 2. Load and Prepare Data ---

# Detect if running on shinyapps.io or locally
if (file.exists("harmonized_dataset.csv")) {
  # Running on shinyapps.io or files are in the same directory as app
  df_orig <- read_csv("harmonized_dataset.csv", show_col_types = FALSE)
} else if (file.exists("data/harmonized_dataset.csv")) {
  # Common local setup with a data folder
  df_orig <- read_csv("data/harmonized_dataset.csv", show_col_types = FALSE)
} else {
  # Fallback for testing/dummy data if file is missing
  df_orig <- data.frame(
    record_id = paste0("R", 101:120),
    age = sample(20:70, 20, replace = TRUE),
    acr_u = runif(20, 10, 200),
    group = factor(sample(c("CKD", "Healthy", "Diabetes"), 20, replace = TRUE)),
    sex = factor(sample(c("Male", "Female"), 20, replace = TRUE)),
    measurement_date = seq(as.Date("2023/01/01"), by = "month", length.out = 20)
  )
  warning("Using dummy data: 'harmonized_dataset.csv' not found in expected locations.")
}

# Ensure all column names are valid R variables (replace spaces/special chars with '_')
names(df_orig) <- make.names(names(df_orig))


# --- 3. Define UI ---
ui <- fluidPage(
  useShinyjs(), # Enable shinyjs functions
  tags$head(
    # Custom CSS for a cleaner, modern look
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
      body {
        font-family: 'Inter', sans-serif;
        background-color: #f3f4f6; /* Light gray background */
      }
      .container-fluid {
        max-width: 1400px;
        margin-left: auto;
        margin-right: auto;
      }
      .sidebar {
        background-color: #ffffff;
        padding: 20px;
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
      }
      .data-output {
        background-color: #ffffff;
        padding: 20px;
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        margin-bottom: 20px;
      }
      .nav-tabs {
        margin-top: 10px;
      }
      .shiny-input-container {
        margin-bottom: 15px;
      }
      .btn {
        background-color: #10b981; /* Emerald green */
        color: white;
        border: none;
        border-radius: 5px;
        padding: 8px 15px;
        transition: background-color 0.3s;
      }
      .btn:hover {
        background-color: #059669;
      }
    "))
  ),

  # Authentication UI
  shinyauthr::loginUI(id = "login", title = "Data Explorer Login"),

  # Main application UI, hidden until authenticated
  hidden(
    div(id = "app_content",
      navbarPage("Data Explorer Dashboard",
        # --- Data Filtering Tab ---
        tabPanel("Data Explorer",
          sidebarLayout(
            sidebarPanel(
              class = "sidebar",
              h3("Data Filtering"),
              selectInput("group_filter", "Filter by Group", choices = c("All", unique(df_orig$group)), selected = "All"),
              sliderInput("age_range", "Filter by Age",
                          min = floor(min(df_orig$age, na.rm = TRUE)),
                          max = ceiling(max(df_orig$age, na.rm = TRUE)),
                          value = c(floor(min(df_orig$age, na.rm = TRUE)), ceiling(max(df_orig$age, na.rm = TRUE)))
              ),
              selectInput("sex_filter", "Filter by Sex", choices = c("All", unique(df_orig$sex)), selected = "All"),
              dateRangeInput("date_range", "Filter by Date Range",
                             start = min(df_orig$measurement_date, na.rm = TRUE),
                             end = max(df_orig$measurement_date, na.rm = TRUE),
                             min = min(df_orig$measurement_date, na.rm = TRUE),
                             max = max(df_orig$measurement_date, na.rm = TRUE)
              ),
              actionButton("reset_filters", "Reset All Filters", class = "btn")
            ),
            mainPanel(
              class = "data-output",
              h3("Filtered Data (Showing First 25 Rows by default)"),
              verbatimTextOutput("num_records"),
              DTOutput("filtered_data_table")
            )
          )
        ),

        # --- Descriptive Statistics Tab ---
        tabPanel("Descriptive Stats",
          fluidRow(
            column(12,
              h3("Descriptive Statistics for Filtered Numeric Variables"),
              DTOutput("descriptive_table")
            )
          )
        ),

        # --- Plots Tab ---
        tabPanel("Plots",
          sidebarLayout(
            sidebarPanel(
              class = "sidebar",
              h3("Plot Settings"),
              selectInput("plot_type", "Select Plot Type",
                          choices = c("Histogram", "Scatter Plot", "Box Plot"),
                          selected = "Histogram"
              ),
              uiOutput("plot_ui"), # UI for variable selection
              selectInput("group_variable", "Categorical Grouping/Coloring Variable", choices = NULL), # New input for grouping/coloring
              actionButton("generate_plot", "Generate Plot", class = "btn")
            ),
            mainPanel(
              class = "data-output",
              h3("Data Visualization"),
              textOutput("plot_message"),
              plotOutput("main_plot")
            )
          )
        ),

        # --- Unique Record IDs Tab (NEW FEATURE) ---
        tabPanel("Unique Record IDs",
          fluidRow(
            column(12,
              class = "data-output",
              h3("List of Unique Record IDs"),
              p("This list shows all unique record_id values after applying any filters."),
              div(style = "max-height: 70vh; overflow-y: auto; padding: 10px; border: 1px solid #e5e7eb; border-radius: 5px;",
                  verbatimTextOutput("unique_record_ids_output")
              )
            )
          )
        )
      )
    )
  )
)


# --- 4. Define Server Logic ---
server <- function(input, output, session) {

  # --- 4.1. Authentication ---
  credentials <- shinyauthr::loginServer(
    id = "login",
    data = user_base,
    user_col = user,
    pwd_col = password,
    log_out = reactive(FALSE) # Keep login visible after success
  )

  # Check authentication status and toggle UI visibility
  user_auth <- reactive({
    req(credentials$user_auth)
    if (credentials$user_auth) {
      show("app_content")
      hide("login_box") # Hide login form after success
    } else {
      hide("app_content")
      show("login_box")
    }
    credentials
  })

  # --- 4.2. Data Filtering ---
  filtered_data <- reactive({
    req(user_auth$authenticated)
    df <- df_orig

    # Group Filter
    if (input$group_filter != "All") {
      df <- df %>% filter(group == input$group_filter)
    }

    # Age Filter
    df <- df %>% filter(age >= input$age_range[1], age <= input$age_range[2])

    # Sex Filter
    if (input$sex_filter != "All") {
      df <- df %>% filter(sex == input$sex_filter)
    }

    # Date Filter
    df <- df %>% filter(measurement_date >= input$date_range[1], measurement_date <= input$date_range[2])

    df
  })

  # Reset filters
  observeEvent(input$reset_filters, {
    updateSelectInput(session, "group_filter", selected = "All")
    updateSliderInput(session, "age_range", value = c(min(df_orig$age, na.rm = TRUE), max(df_orig$age, na.rm = TRUE)))
    updateSelectInput(session, "sex_filter", selected = "All")
    updateDateRangeInput(session, "date_range",
                         start = min(df_orig$measurement_date, na.rm = TRUE),
                         end = max(df_orig$measurement_date, na.rm = TRUE)
    )
  })


  # --- 4.3. Data Outputs ---

  # Number of records
  output$num_records <- renderText({
    paste("Number of records:", nrow(filtered_data()))
  })

  # Main Data Table (Improvement 4: Increase pageLength to 25)
  output$filtered_data_table <- renderDT({
    req(user_auth$authenticated)
    df <- filtered_data()
    datatable(df, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
  })

  # Descriptive Statistics Table (Improvement 4: Increase pageLength to 25)
  output$descriptive_table <- renderDT({
    req(user_auth$authenticated)
    df <- filtered_data() %>%
      select(where(is.numeric)) # Only select numeric columns

    if (ncol(df) == 0) {
      return(NULL)
    }

    descriptive_stats <- df %>%
      summarise(
        across(everything(), list(
          N = ~ sum(!is.na(.)),
          Mean = ~ mean(., na.rm = TRUE),
          SD = ~ sd(., na.rm = TRUE),
          Min = ~ min(., na.rm = TRUE),
          Max = ~ max(., na.rm = TRUE)
        ), .names = "{.fn}_{.col}")
      ) %>%
      pivot_longer(
        cols = everything(),
        names_to = c(".value", "Variable"),
        names_sep = "_"
      ) %>%
      select(Variable, N, Mean, SD, Min, Max)

    datatable(descriptive_stats,
              options = list(pageLength = 25, scrollX = TRUE), # Increased pageLength
              rownames = FALSE
    ) %>%
      formatRound(columns = c("Mean", "SD", "Min", "Max"), digits = 2)
  })


  # --- 4.4. Plots ---

  # Dynamically update plot variable choices based on filtered data
  observeEvent(filtered_data(), {
    df <- filtered_data()
    all_vars <- names(df)
    numeric_vars <- names(df)[sapply(df, is.numeric)]
    # Non-numeric variables are candidates for categorical grouping/coloring
    categorical_vars <- names(df)[!sapply(df, is.numeric)]

    updateSelectInput(session, "plot_variable", choices = c("", numeric_vars))
    updateSelectInput(session, "y_variable", choices = c("", numeric_vars))
    # Update the new categorical grouping input
    updateSelectInput(session, "group_variable", choices = c("None", categorical_vars))
  })


  # Conditional UI for y_variable
  output$plot_ui <- renderUI({
    req(input$plot_type)
    if (input$plot_type == "Scatter Plot") {
      list(
        selectInput("plot_variable", "X-Axis Variable", choices = NULL),
        selectInput("y_variable", "Y-Axis Variable", choices = NULL)
      )
    } else if (input$plot_type %in% c("Histogram", "Box Plot")) {
      selectInput("plot_variable", "Numeric Variable", choices = NULL)
    }
  })

  # Plot generation
  output$main_plot <- renderPlot({
    req(user_auth$authenticated)
    df <- filtered_data()
    req(input$plot_type)
    req(input$plot_variable)

    # Check for empty data frame
    if (nrow(df) == 0) {
      output$plot_message <- renderText({"No data to plot after filtering."})
      return(NULL)
    }

    # Reset message
    output$plot_message <- renderText({""})

    # Prepare grouping variable if selected
    group_var <- input$group_variable
    use_grouping <- !is.null(group_var) && group_var != "None" && group_var != "" && group_var %in% names(df)
    
    # Plot generation switch
    p <- switch(input$plot_type,
                "Histogram" = {
                  if (!is.numeric(df[[input$plot_variable]])) {
                    output$plot_message <- renderText({"Please select a numeric variable for a Histogram."})
                    return(NULL)
                  }
                  ggplot(df, aes(x = !!sym(input$plot_variable))) +
                    geom_histogram(fill = "#10b981", color = "#1f2937", bins = 20) +
                    labs(title = paste("Histogram of", input$plot_variable), x = input$plot_variable, y = "Frequency") +
                    theme_minimal()
                },
                
                "Scatter Plot" = {
                  # --- IMPROVEMENT 2: Color Scatter Plots ---
                  req(input$y_variable)
                  if (!(input$y_variable %in% names(df))) {
                    return(NULL)
                  }
                  if (!is.numeric(df[[input$plot_variable]]) || !is.numeric(df[[input$y_variable]])) {
                    output$plot_message <- renderText({"Please select two numeric variables for a Scatter Plot."})
                    return(NULL)
                  }
                  
                  # Base aesthetic with optional color mapping
                  base_aes <- aes(x = !!sym(input$plot_variable), y = !!sym(input$y_variable))
                  color_map <- NULL
                  
                  if (use_grouping && !is.numeric(df[[group_var]])) {
                    # Only apply color if group_var is categorical
                    base_aes <- aes(x = !!sym(input$plot_variable), y = !!sym(input$y_variable), color = !!sym(group_var))
                    color_map <- scale_color_brewer(palette = "Set1")
                  }
                  
                  ggplot(df, base_aes) +
                    geom_point(alpha = 0.7) +
                    labs(title = paste("Scatter Plot of", input$plot_variable, "vs", input$y_variable, 
                                       if(use_grouping && !is.numeric(df[[group_var]])) paste("colored by", group_var)), 
                         x = input$plot_variable, y = input$y_variable) +
                    theme_minimal() +
                    color_map
                },
                
                "Box Plot" = {
                  # --- IMPROVEMENT 1: Box Plot ---
                  req(group_var)
                  if (!use_grouping || is.numeric(df[[group_var]])) {
                    output$plot_message <- renderText({"Please select a categorical grouping variable (X-axis) for a Box Plot."})
                    return(NULL)
                  }
                  if (!is.numeric(df[[input$plot_variable]])) {
                    output$plot_message <- renderText({"Please select a numeric variable (Y-axis) for a Box Plot."})
                    return(NULL)
                  }
                  
                  ggplot(df, aes(x = !!sym(group_var), y = !!sym(input$plot_variable), fill = !!sym(group_var))) +
                    geom_boxplot(outlier.shape = 16, outlier.alpha = 0.5) +
                    geom_jitter(color = "#4b5563", alpha = 0.3, width = 0.2) + # Add jittered points
                    labs(title = paste("Box Plot of", input$plot_variable, "by", group_var),
                         x = group_var, y = input$plot_variable, fill = group_var) +
                    theme_minimal() +
                    scale_fill_brewer(palette = "Pastel1")
                },
                {
                  return(NULL)
                }
    )

    print(p)
  }, height = 500)
  
  # --- 4.5 Unique Record IDs (NEW FEATURE) ---
  output$unique_record_ids_output <- renderPrint({
    req(user_auth$authenticated)
    df <- filtered_data()
    
    # Check if 'record_id' column exists
    if ("record_id" %in% names(df)) {
        unique_ids <- sort(unique(df$record_id))
        
        if (length(unique_ids) > 0) {
            # Print the IDs in a clear, column-friendly format
            cat(paste("Total Unique Records:", length(unique_ids), "\n\n"))
            # Print as a formatted list, 5 IDs per line
            ids_per_line <- 5
            matrix(unique_ids, ncol = ids_per_line, byrow = TRUE) %>%
                apply(1, function(row) {
                    cat(paste(sprintf("%-15s", row), collapse = "  "), "\n")
                })
        } else {
            cat("No unique record IDs found after filtering.")
        }
    } else {
        cat("The 'record_id' column is not present in the dataset.")
    }
  })
}

# --- 5. Run the App ---
shinyApp(ui = ui, server = server)
