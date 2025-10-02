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
  # Running on shinyapps.io or files are in the same directory as app.R
  cat("Loading data from app directory\n")
  data_path <- "harmonized_dataset.csv"
  dictionary_path <- "data_dictionary_master.csv"
} else {
  # Running locally - use user-specific paths
  os_user <- Sys.info()[["user"]]
  
  if (os_user == "choiyej") {
    base_data_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
    data_path <- file.path(base_data_path, "Data Harmonization/Data Clean/harmonized_dataset.csv")
    dictionary_path <- file.path(base_data_path, "Data Harmonization/data_dictionary_master.csv")
  } else if (os_user == "pylell") {
    base_data_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
    data_path <- file.path(base_data_path, "Data Harmonization/Data Clean/harmonized_dataset.csv")
    dictionary_path <- file.path(base_data_path, "Data Harmonization/data_dictionary_master.csv")
  } else if (os_user == "shivaniramesh") {
    base_data_path <- "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean"
    data_path <- file.path(base_data_path, "harmonized_dataset.csv")
    dictionary_path <- file.path(dirname(base_data_path), "Data Clean/data_dictionary_master.csv")
  } else {
    stop("Unknown user: please specify root path for this user.")
  }
}

# Read the data files
cat("Attempting to load data from:", data_path, "\n")
data <- read_csv(data_path)
cat("Successfully loaded data with", nrow(data), "rows\n")

cat("Attempting to load dictionary from:", dictionary_path, "\n")
dictionary_data <- read_csv(dictionary_path)
cat("Successfully loaded dictionary with", nrow(dictionary_data), "rows\n")


# The file reading commands are now moved outside the if/else block
# since data_path and dictionary_path are guaranteed to be defined for all users.
data <- read_csv(data_path)
dictionary_data <- read_csv(dictionary_path)

# --- DE-IDENTIFICATION: Replace MRN with UUID ---
if ("mrn" %in% names(data)) {
  cat("De-identifying data: replacing MRN with UUID\n")
  
  # Check if we're on shinyapps.io (map file in same directory) or local
  if (file.exists("mrn_uuid_map.csv")) {
    mrn_map_path <- "mrn_uuid_map.csv"
  } else {
    # For local users, assuming map is in same directory as harmonized data
    mrn_map_path <- file.path(dirname(data_path), "mrn_uuid_map.csv")
  }
  
  # Load the mapping - read all columns as character
  mrn_uuid_map <- read_csv(mrn_map_path, col_types = cols(.default = col_character()))
  
  cat("Loaded", nrow(mrn_uuid_map), "MRN-UUID mappings\n")
  
  # Standardize both MRN columns: convert to character and trim whitespace
  data <- data %>%
    mutate(mrn = trimws(as.character(mrn)))
  
  # FIX: Remove .0 from mapping file MRNs and trim whitespace
  mrn_uuid_map <- mrn_uuid_map %>%
    mutate(mrn = trimws(gsub("\\.0$", "", mrn)))
  
  cat("Sample data MRNs:", paste(head(unique(data$mrn), 3), collapse=", "), "\n")
  cat("Sample map MRNs (after cleanup):", paste(head(mrn_uuid_map$mrn, 3), collapse=", "), "\n")
  
  # Check for matches
  matches <- sum(data$mrn %in% mrn_uuid_map$mrn)
  cat("Matching MRNs found:", matches, "out of", length(unique(data$mrn)), "unique MRNs in data\n")
  
  # Perform the join
  data <- data %>%
    left_join(mrn_uuid_map, by = "mrn")
  
  # Replace MRN with UUID
  if ("uuid" %in% names(data)) {
    data <- data %>%
      select(-mrn) %>%
      rename(mrn = uuid)
    
    cat("Successfully replaced MRN with UUID\n")
    cat("Sample UUID values:", paste(head(na.omit(data$mrn), 3), collapse=", "), "\n")
    cat("Rows with missing UUID mapping:", sum(is.na(data$mrn)), "\n")
  } else {
    cat("ERROR: 'uuid' column not found after join!\n")
    stop("De-identification failed")
  }
  
  # Reorder columns
  key_cols <- c("record_id", "study", "group", "sex", "visit", "procedure", "mrn")
  existing_key_cols <- intersect(key_cols, names(data))
  other_cols <- setdiff(names(data), existing_key_cols)
  
  data <- data %>%
    select(all_of(existing_key_cols), all_of(other_cols))
  
} else {
  cat("No MRN column found in data\n")
}
# --- END DE-IDENTIFICATION ---

# Trim whitespace from the study column to ensure correct filtering
data <- data %>%
  mutate(study = trimws(as.character(study)))

# Get unique values and ranges for filters
# The `study` column is now used for filtering
unique_studies <- sort(unique(data$study))

# Get unique form names from the dictionary
unique_form_names <- sort(unique(dictionary_data$form_name))

# Get all data variable names for manual selection (excluding key IDs that are always shown)
vars_to_exclude_from_selector <- c("record_id", "visit", "procedure", "study")
all_data_variables_for_selection <- sort(setdiff(names(data), vars_to_exclude_from_selector))

# Calculate overall statistics
total_studies <- length(unique_studies)
total_rows <- nrow(data)
avg_rows_per_study <- round(total_rows / total_studies, 0)

# --- 3. Define the User Interface (UI) ---
ui <- fluidPage(
  
  # 1. Initialize shinyjs and define the login UI
  shinyjs::useShinyjs(),
  shinyauthr::loginUI(id = "login", title = "Harmonized Data Login"),
  
  # 2. Wrap all main content in a div that is initially hidden
  div(
    id = "main_content",
    class = "hidden", # This class will be removed by shinyjs on successful login
    
    # Set up a title and page layout
    tags$head(
      tags$style(HTML("
        /* ADDED: Define hidden class for initial state */
        .hidden { display: none; } 

        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
        body { font-family: 'Inter', sans-serif; background-color: #f3f4e6; color = #374151; }
        /* MODIFIED: Added flex for layout of title and logout button */
        .main-header { background-color: #4f46e5; color: white; padding: 1rem; box-shadow: 0 4px 6px rgba(0,0,0,0.1); display: flex; justify-content: space-between; align-items: center; }
        .card { background-color: #ffffff; border-radius: 12px; box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1); padding: 24px; }
        .sidebar { background-color: #ffffff; border-right: 1px solid #e5e7eb; border-radius: 12px; padding: 24px; }
        .stat-box { border-radius: 8px; padding: 1.5rem; text-align: center; }
        .dt-container { overflow-x: auto; }
        .tab-content { padding: 1.5rem; }
        .btn-danger { background-color: #dc2626; border-color: #b91c1c; color: white; border-radius: 6px; }
        .btn-danger:hover { background-color: #b91c1c; border-color: #991b1b; }
        /* Utility classes for text (since the original used Tailwind classes) */
        .text-3xl { font-size: 1.875rem; line-height: 2.25rem; }
        .text-xl { font-size: 1.25rem; line-height: 1.75rem; }
        .text-lg { font-size: 1.125rem; line-height: 1.75rem; }
        .text-sm { font-size: 0.875rem; line-height: 1.25rem; }
        .font-bold { font-weight: 700; }
        .font-semibold { font-weight: 600; }
        .mb-4 { margin-bottom: 1rem; }
        .mb-2 { margin-bottom: 0.5rem; }
        .mb-6 { margin-bottom: 1.5rem; }
        .mt-4 { margin-top: 1rem; }
        .text-indigo-700 { color: #4338ca; }
      "))
    ),
    
    # App Title and Logout Button
    tags$div(class = "main-header", 
             h1("Harmonized Dataset Dashboard", class = "text-3xl font-bold"),
             shinyauthr::logoutUI(id = "logout", label = "Logout", class = "btn btn-danger") # Added Logout
    ),
    
    # Main layout with sidebar and main panel
    sidebarLayout(
      
      # Sidebar panel for user inputs
      sidebarPanel(
        class = "sidebar",
        h2("Filter Controls", class = "text-xl font-bold mb-4"),
        
        # Checkboxes to select one or more studies
        checkboxGroupInput(
          "selected_study",
          "Select Study:",
          choices = unique_studies,
          selected = character(0)
        ),
        
        # Filter for record ID
        textInput(
          "selected_record_id",
          "Search by Record ID:"
        ),
        
        # Dynamic filters for group and sex
        uiOutput("dynamic_filters"),
        
        # New dynamic filter for procedure
        uiOutput("procedure_ui"),
        
        # --- Variable Categories Section (Moved Back to Sidebar) ---
        tags$hr(),
        h2("Select Variable Categories", class = "text-xl font-bold mb-2"),
        checkboxGroupInput(
          "selected_form_names",
          NULL,
          choices = unique_form_names,
          selected = character(0)
        ),
        # --- End Variable Categories Section ---
        
        # --- Visuals Section ---
        tags$hr(),
        h2("Visualizations", class = "text-xl font-bold mb-4"),
        
        # Select plot type
        selectInput(
          "plot_type",
          "Chart Type:",
          choices = c("Bar Chart", "Histogram", "Scatter Plot")
        ),
        
        # Dynamic UI for plot variables
        uiOutput("plot_variable_ui"),
        uiOutput("y_variable_ui"),
        
        # This will be used to hold PET variables for backwards compatibility
        uiOutput("pet_variables_ui")
      ),
      
      # Main panel for displaying outputs
      mainPanel(
        class = "main-content",
        
        # Data Overview section with statistics
        tags$div(class = "card mb-6",
                 h2("Data Overview", class = "text-2xl font-bold mb-4"),
                 fluidRow(
                   column(4, tags$div(class = "stat-box bg-blue-50",
                                      p(total_rows, class = "text-xl font-bold text-blue-600"),
                                      p("Total Rows", class = "text-sm font-medium text-gray-500")
                   )),
                   column(4, tags$div(class = "stat-box bg-green-50",
                                      p(total_studies, class = "text-xl font-bold text-green-600"),
                                      p("Total Studies", class = "text-sm font-medium text-gray-500")
                   )),
                   column(4, tags$div(class = "stat-box bg-yellow-50",
                                      p(avg_rows_per_study, class = "text-xl font-bold text-yellow-600"),
                                      p("Average Rows per Study", class = "text-sm font-medium text-gray-500")
                   ))
                 )
        ),
        
        # Tabs for different views
        tags$div(class = "card",
                 tabsetPanel(
                   id = "main_tabs",
                   tabPanel("Descriptive Statistics", 
                            DTOutput("descriptive_table"),
                            downloadButton("download_descriptive", "Download Table")
                   ),
                   tabPanel("Data Availability",
                            h4("Procedure Availability"),
                            DTOutput("data_availability_table"),
                            tags$hr(),
                            h4("Variable Availability Search"),
                            textInput("variable_search", "Search for a variable:"),
                            DTOutput("searched_variable_availability_table"),
                            tags$hr(),
                            h4("Variable Availability by Study"),
                            DTOutput("searched_variable_availability_by_study_table"),
                            downloadButton("download_availability", "Download Table")
                   ),
                   tabPanel("Raw Data Table",
                            
                            # --- START: Manual Variable Search Kept Here ---
                            tags$div(class = "card p-4 mb-4 bg-gray-50 border border-gray-200",
                                     h3("Manual Variable Selection", class = "text-xl font-semibold mb-2 text-indigo-700"),
                                     p("Use the search bar below to add specific variables to the table, in addition to the categories selected in the sidebar.", class="text-sm text-gray-600 mb-4"),
                                     selectizeInput(
                                       "manual_variables_select",
                                       "Type or Select Specific Variables:",
                                       choices = all_data_variables_for_selection,
                                       multiple = TRUE,
                                       selected = character(0),
                                       options = list(placeholder = "Start typing variable names...")
                                     )
                            ),
                            # --- END: Manual Variable Search Kept Here ---
                            
                            DTOutput("raw_data_table"),
                            downloadButton("download_raw", "Download Raw Data")
                   ),
                   tabPanel("Data Visuals",
                            textOutput("plot_message"),
                            plotOutput("data_plot"),
                            downloadButton("download_plot", "Download Plot")
                   ),
                   # New tab for the data dictionary table
                   tabPanel("Data Dictionary",
                            DTOutput("dictionary_table")
                   )
                 )
        )
      )
    )
  ) # End of div#main_content
)

# --- 4. Define the Server Logic ---
server <- function(input, output, session) {
  
  # --- Authentication Logic Start ---
  
  # 1. Define the logout reactive value signal
  logout_init <- shinyauthr::logoutServer(
    id = "logout",
    active = reactive(credentials()$user_auth) 
  )
  
  # 2. Set up the authentication mechanism
  credentials <- shinyauthr::loginServer(
    id = "login",
    data = user_base,
    user_col = "user",
    pwd_col = "password",
    log_out = reactive(logout_init())
  )
  
  # 3. Password-Only Access Setup: Hide username field and auto-populate
  fixed_user <- user_base$user[1]
  js_injected <- reactiveVal(FALSE) 
  
  observe({
    # Only run this if the JS has not been injected yet AND if not authenticated 
    if (!credentials()$user_auth && !js_injected()) {
      shinyjs::runjs(sprintf("
        // Set the fixed username value
        var user_input = document.getElementById('login-user');
        if (user_input) {
          user_input.value = '%s';
        }

        // Hide the entire form group containing the username input and label
        // This is necessary to achieve a password-only login.
        $('#login-user').closest('.form-group').hide();
      ", fixed_user))
      js_injected(TRUE)
    }
  })
  
  # 4. Observe authentication status and toggle main content visibility
  observe({
    if (credentials()$user_auth) {
      shinyjs::removeClass(selector = "#main_content", class = "hidden")
      shinyjs::hide(id = "login")
      cat("User logged in:", credentials()$info$name, "\n")
    } else {
      shinyjs::addClass(selector = "#main_content", class = "hidden")
      shinyjs::show(id = "login")
      js_injected(FALSE) # Reset the flag on logout
    }
  })
  
  # --- Authentication Logic End ---
  
  
  # Observe changes in selected study for debugging
  observe({
    req(credentials()$user_auth) # Require authentication
    print("--- User Input Debug ---")
    print(paste("Selected studies:", paste(input$selected_study, collapse = ", ")))
  })
  
  # Reactive expression to dynamically create filters based on data columns
  output$dynamic_filters <- renderUI({
    req(credentials()$user_auth) # Require authentication
    tagList(
      if ("group" %in% names(data)) {
        checkboxGroupInput(
          "selected_group",
          "Select Group(s):",
          choices = sort(unique(data$group)),
          selected = character(0)
        )
      },
      if ("sex" %in% names(data)) {
        checkboxGroupInput(
          "selected_sex",
          "Select Sex:",
          choices = sort(unique(data$sex)),
          selected = character(0)
        )
      }
    )
  })
  
  # Reactive expression to dynamically create procedure filters
  output$procedure_ui <- renderUI({
    req(credentials()$user_auth) # Require authentication
    if ("procedure" %in% names(data)) {
      tagList(
        h3("Select Procedures", class = "text-lg font-bold mt-4 mb-2"),
        checkboxGroupInput(
          "selected_procedure",
          NULL,
          choices = sort(unique(data$procedure)),
          selected = character(0)
        )
      )
    }
  })
  
  # Reactive expression to get variable names from the data dictionary based on selected form names
  selected_form_variables <- reactive({
    req(credentials()$user_auth) # Require authentication
    if (is.null(input$selected_form_names)) {
      return(character(0))
    }
    
    dictionary_data %>%
      filter(form_name %in% input$selected_form_names) %>%
      pull(variable_name) %>%
      unique()
  })
  
  # Reactive expression to filter data based on user inputs
  filtered_data <- reactive({
    req(credentials()$user_auth) # Require authentication
    
    # --- DEBUGGING START ---
    print("--- Filtered Data Reactive Expression Started ---")
    data_filtered <- data
    print(paste("Initial number of rows:", nrow(data_filtered)))
    # --- DEBUGGING END ---
    
    # Filter by selected record ID
    if (!is.null(input$selected_record_id) && input$selected_record_id != "") {
      data_filtered <- data_filtered %>%
        filter(record_id == input$selected_record_id)
      # --- DEBUGGING START ---
      print(paste("Number of rows after record ID filter:", nrow(data_filtered)))
      # --- DEBUGGING END ---
    }
    
    # Filter by selected studies
    # This filter is only applied if a study is selected
    if (length(input$selected_study) > 0) {
      data_filtered <- data_filtered %>%
        filter(study %in% input$selected_study)
      # --- DEBUGGING START ---
      print(paste("Number of rows after study filter:", nrow(data_filtered)))
      # --- DEBUGGING END ---
    }
    
    # Filter by selected groups (if filter exists)
    if (length(input$selected_group) > 0) {
      data_filtered <- data_filtered %>%
        filter(group %in% input$selected_group)
      # --- DEBUGGING START ---
      print(paste("Number of rows after group filter:", nrow(data_filtered)))
      # --- DEBUGGING END ---
    }
    
    # Filter by selected sex (if filter exists)
    if (length(input$selected_sex) > 0) {
      data_filtered <- data_filtered %>%
        filter(sex %in% input$selected_sex)
      # --- DEBUGGING START ---
      print(paste("Number of rows after sex filter:", nrow(data_filtered)))
      # --- DEBUGGING END ---
    }
    
    # Filter by selected procedures (if filter exists)
    if (length(input$selected_procedure) > 0) {
      data_filtered <- data_filtered %>%
        filter(procedure %in% input$selected_procedure)
      # --- DEBUGGING START ---
      print(paste("Number of rows after procedure filter:", nrow(data_filtered)))
      # --- DEBUGGING END ---
    }
    
    # --- DEBUGGING START ---
    print("--- Filtered Data Reactive Expression Finished ---")
    print(paste("Final number of rows:", nrow(data_filtered)))
    # --- DEBUGGING END ---
    
    data_filtered
  })
  
  # Reactive expression to create the descriptive statistics table
  output$descriptive_table <- renderDT({
    req(filtered_data())
    
    # Create a local copy to avoid repeated function calls
    df_temp <- filtered_data()
    
    # NEW FIX: Handle NA group values
    df_temp$group[is.na(df_temp$group)] <- "NA Group"
    
    # Check if 'group' column exists before grouping
    if (!"group" %in% names(df_temp)) {
      descriptive_summary <- data.frame("Note" = "Group column not found in filtered data.")
      return(datatable(descriptive_summary, options = list(dom = 't')))
    }
    
    # Calculate group-level summaries and grand total in one step
    descriptive_summary <- df_temp %>%
      bind_rows(mutate(., group = "Grand Total")) %>%
      group_by(group) %>%
      summarise(
        `Record ID (N)` = length(unique(record_id)),
        # Use a check for each column before summarizing to prevent errors
        `Age (avg)` = if("age" %in% names(df_temp)) mean(age, na.rm = TRUE) else NA_real_,
        `BMI (avg)` = if("bmi" %in% names(df_temp)) mean(bmi, na.rm = TRUE) else NA_real_,
        `HbA1c (avg)` = if("hba1c" %in% names(df_temp)) mean(hba1c, na.rm = TRUE) else NA_real_,
        `ACR-U (med)` = if("acr_u" %in% names(df_temp)) median(acr_u, na.rm = TRUE) else NA_real_,
        `Female (%)` = if("sex" %in% names(df_temp)) mean(sex == "Female", na.rm = TRUE) * 100 else NA_real_,
        `Hispanic (%)` = if("ethnicity" %in% names(df_temp)) mean(ethnicity == "Hispanic or Latino", na.rm = TRUE) * 100 else NA_real_,
        `eGFR_CKD_epi (avg)` = if("eGFR_CKD_epi" %in% names(df_temp)) mean(eGFR_CKD_epi, na.rm = TRUE) else NA_real_,
        `Total Kidney Volume (avg)` = if("total_kidney_volume_ml" %in% names(df_temp)) mean(total_kidney_volume_ml, na.rm = TRUE) else NA_real_,
        `Manual Total Kidney Volume (avg)` = if("total_kidney_volume_ml_manual" %in% names(df_temp)) mean(total_kidney_volume_ml_manual, na.rm = TRUE) else NA_real_,
        `ERPF Raw Plasma (avg)` = if("erpf raw plasma" %in% names(df_temp)) mean(`erpf raw plasma`, na.rm = TRUE) else NA_real_,
        `ERPF BSA Plasma (avg)` = if("erpf bsa plasma" %in% names(df_temp)) mean(`erpf bsa plasma`, na.rm = TRUE) else NA_real_,
        `PGLO (avg)` = if("pglo" %in% names(df_temp)) mean(pglo, na.rm = TRUE) else NA_real_,
        `RA (avg)` = if("ra" %in% names(df_temp)) mean(ra, na.rm = TRUE) else NA_real_,
        `RE (avg)` = if("re" %in% names(df_temp)) mean(re, na.rm = TRUE) else NA_real_,
        .groups = 'drop'
      ) %>%
      # Round all numeric columns to 2 decimal places
      mutate(across(where(is.numeric), ~ round(., 2)))
    
    # Transpose the table to have groups as columns and metrics as rows
    descriptive_table_display <- descriptive_summary %>%
      rename(Category = group) %>%
      tidyr::pivot_longer(cols = -Category, names_to = "Metric", values_to = "Count") %>%
      tidyr::pivot_wider(names_from = Category, values_from = Count) %>%
      select(Metric, everything())
    
    datatable(descriptive_table_display, options = list(dom = 't'))
  })
  
  # Reactive expression to create the main data availability table (procedures)
  main_availability_table <- reactive({
    req(filtered_data())
    
    df <- filtered_data()
    
    # NEW FIX: Handle NA group values
    df$group[is.na(df$group)] <- "NA Group"
    
    # Define the rows and their corresponding calculation logic
    availability_rows <- c(
      "Participants",
      "Brain_Biomarkers", "Clamp", "Dxa", "MRI", "Ivgtt", "MINMOD",
      "Kidney Biopsy", "Kidney Morphometris", "PET Scan",
      "Metabolomics - Blood", "Metabolomics - Tissue", "Renal Clearance",
      "Hemodynamic Outcomes", "Proteomics"
    )
    
    # Define the calculation logic for each category
    calculate_availability <- function(group_data) {
      if(nrow(group_data) == 0) return(rep(0, length(availability_rows)))
      
      counts <- c(
        Participants = length(unique(group_data$record_id)),
        Brain_Biomarkers = length(unique(group_data$record_id[group_data$procedure == "brain_biomarkers"])),
        Clamp = length(unique(group_data$record_id[group_data$procedure == "clamp"])),
        Dxa = length(unique(group_data$record_id[group_data$procedure == "dxa"])),
        MRI = length(unique(group_data$record_id[group_data$procedure %in% c("bold_mri", "cardio_abdominal_mri", "mri", "pc_mri_2d")])),
        Ivgtt = length(unique(group_data$record_id[group_data$procedure == "ivgtt"])),
        MINMOD = {
          vars <- dictionary_data %>% filter(form_name == "ivgtt_minmod") %>% pull(variable_name)
          if(length(vars) > 0) {
            matched_cols <- intersect(vars, names(group_data))
            if(length(matched_cols) > 0) length(unique(group_data$record_id[rowSums(is.na(group_data[, matched_cols, drop = FALSE])) < length(matched_cols)])) else 0
          } else 0
        },
        `Kidney Biopsy` = length(unique(group_data$record_id[group_data$procedure == "kidney_biopsy"])),
        `Kidney Morphometris` = length(unique(group_data$record_id[grepl("kidney", group_data$procedure, ignore.case = TRUE)])),
        `PET Scan` = length(unique(group_data$record_id[group_data$procedure %in% c("liver_pet_scan", "pet_scan")])),
        `Metabolomics - Blood` = length(unique(group_data$record_id[group_data$procedure == "metabolomics_blood"])),
        `Metabolomics - Tissue` = length(unique(group_data$record_id[group_data$procedure == "metabolomics_tissue"])),
        `Renal Clearance` = length(unique(group_data$record_id[group_data$procedure == "renal_clearance_testing"])),
        `Hemodynamic Outcomes` = length(unique(group_data$record_id[grepl("hemodynamic", group_data$procedure, ignore.case = TRUE)])),
        Proteomics = {
          vars <- dictionary_data %>% filter(form_name == "proteomics") %>% pull(variable_name)
          if(length(vars) > 0) {
            matched_cols <- intersect(vars, names(group_data))
            if(length(matched_cols) > 0) length(unique(group_data$record_id[rowSums(is.na(group_data[, matched_cols, drop = FALSE])) < length(matched_cols)])) else 0
          } else 0
        }
      )
      
      return(counts)
    }
    
    # Calculate group-level summaries
    group_summary <- df %>%
      group_by(group) %>%
      summarise(counts = list(calculate_availability(cur_data()))) %>%
      tidyr::unnest_wider(counts)
    
    # Calculate grand total
    grand_total_counts <- as.data.frame(t(calculate_availability(df)))
    grand_total_counts$group <- "Grand Total"
    
    # Bind rows and arrange the final table
    final_table <- bind_rows(group_summary, grand_total_counts)
    
    # Make a final table for display
    final_table_display <- final_table %>%
      rename(Category = group) %>%
      tidyr::pivot_longer(cols = -Category, names_to = "Metric", values_to = "Count") %>%
      tidyr::pivot_wider(names_from = Category, values_from = Count) %>%
      select(Metric, everything())
    
    final_table_display
  })
  
  # Reactive expression to create the variable-specific data availability table
  searched_variable_availability <- reactive({
    req(credentials()$user_auth) # Require authentication
    req(input$variable_search)
    df <- filtered_data()
    
    if (nchar(input$variable_search) > 0) {
      # Find variables that match the search term
      matched_vars <- names(df)[grepl(input$variable_search, names(df), ignore.case = TRUE)]
      
      if (length(matched_vars) > 0) {
        # NEW FIX: Handle NA group values
        df$group[is.na(df$group)] <- "NA Group"
        
        # Calculate availability for each matched variable, grouped by group
        availability_list <- lapply(matched_vars, function(var) {
          df %>%
            group_by(group) %>%
            summarise(
              Variable = var,
              `Participants with Data (N)` = length(unique(record_id[!is.na(get(var))])),
              .groups = "drop"
            )
        })
        
        # Combine results into a single data frame
        result_df <- bind_rows(availability_list) %>%
          tidyr::pivot_wider(names_from = group, values_from = `Participants with Data (N)`) %>%
          select(Variable, everything())
        
        # Add a grand total column
        grand_total <- df %>%
          summarise(
            Variable = matched_vars,
            `Grand Total` = sapply(matched_vars, function(var) length(unique(record_id[!is.na(get(var))])))
          )
        
        # Join the grand total to the pivoted table
        final_table <- left_join(result_df, grand_total, by = "Variable")
        
        return(final_table)
      }
    }
    NULL
  })
  
  # New reactive expression for the study-specific data availability table
  searched_variable_availability_by_study <- reactive({
    req(credentials()$user_auth) # Require authentication
    req(input$variable_search)
    df <- filtered_data()
    
    if (nchar(input$variable_search) > 0) {
      df$group[is.na(df$group)] <- "NA Group"
      matched_vars <- names(df)[grepl(input$variable_search, names(df), ignore.case = TRUE)]
      if (length(matched_vars) > 0) {
        availability_list <- lapply(matched_vars, function(var) {
          df %>%
            group_by(study) %>%
            summarise(
              Variable = var,
              `Participants with Data (N)` = length(unique(record_id[!is.na(get(var))])),
              .groups = "drop"
            )
        })
        
        result_df <- bind_rows(availability_list) %>%
          tidyr::pivot_wider(names_from = study, values_from = `Participants with Data (N)`) %>%
          select(Variable, everything())
        
        grand_total <- df %>%
          summarise(
            Variable = matched_vars,
            `Total Participants` = sapply(matched_vars, function(var) length(unique(record_id[!is.na(get(var))])))
          )
        
        final_table <- left_join(result_df, grand_total, by = "Variable")
        
        # Filter to keep only columns where at least one value is greater than 0
        cols_to_keep <- c("Variable", "Total Participants")
        numeric_cols <- names(final_table)[sapply(final_table, is.numeric)]
        numeric_cols_to_keep <- numeric_cols[sapply(final_table[, numeric_cols], function(col) any(col > 0, na.rm = TRUE))]
        
        final_table <- final_table[, c(cols_to_keep, numeric_cols_to_keep)]
        
        return(final_table)
      }
    }
    NULL
  })
  
  # Render the main availability table (procedures)
  output$data_availability_table <- renderDT({
    req(credentials()$user_auth) # Require authentication
    datatable(main_availability_table(), options = list(dom = 't'))
  })
  
  # Render the variable search availability table
  output$searched_variable_availability_table <- renderDT({
    req(credentials()$user_auth) # Require authentication
    # Only render this table if a search term is present
    if (!is.null(searched_variable_availability())) {
      datatable(searched_variable_availability(), options = list(dom = 't'))
    } else {
      # Return an empty table if there's no search term
      return(datatable(data.frame(), options = list(dom = 't')))
    }
  })
  
  # Render the new variable availability by study table
  output$searched_variable_availability_by_study_table <- renderDT({
    req(credentials()$user_auth) # Require authentication
    if (!is.null(searched_variable_availability_by_study())) {
      datatable(searched_variable_availability_by_study(), options = list(dom = 't'))
    } else {
      datatable(data.frame(), options = list(dom = 't'))
    }
  })
  
  # Reactive expression to filter columns based on selected form names AND manual selection
  raw_data_display_columns <- reactive({
    req(credentials()$user_auth) # Require authentication
    df_filtered <- filtered_data()
    
    # 1. Get variables from selected forms (categories)
    # This comes from the sidebar
    vars_from_forms <- selected_form_variables()
    
    # 2. Get variables from manual search/selection
    # This comes from the Raw Data Table tab
    vars_from_manual <- input$manual_variables_select
    if (is.null(vars_from_manual)) {
      vars_from_manual <- character(0)
    }
    
    # Combine the two sets of selected variables
    selected_variables <- unique(c(vars_from_forms, vars_from_manual))
    
    # Define mandatory key columns to always be included
    mandatory_cols <- c("record_id", "study", "group", "sex", "visit", "procedure")
    
    if (length(selected_variables) > 0) {
      # Include mandatory columns and selected variables
      columns_to_keep <- unique(c(mandatory_cols, selected_variables))
      
      # Select only the columns that exist in the filtered dataframe
      existing_columns <- intersect(columns_to_keep, names(df_filtered))
      
      # Return the dataframe with the selected columns
      return(df_filtered %>% select(all_of(existing_columns)))
    } else {
      # If no categories or manual variables are selected, return a default set of key columns
      default_cols <- intersect(mandatory_cols, names(df_filtered))
      return(df_filtered %>% select(all_of(default_cols)))
    }
  })
  
  # Reactive expression to create the raw data table with integrated column selection
  output$raw_data_table <- renderDT({
    req(credentials()$user_auth) # Require authentication
    req(raw_data_display_columns())
    
    datatable(
      raw_data_display_columns(),
      extensions = 'Buttons',
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = list(
          'colvis' # Column visibility button
        )
      )
    )
  })
  
  # Reactive expression to create the data dictionary table with the new 'location' column
  dictionary_with_location <- reactive({
    req(credentials()$user_auth) # Require authentication
    # Get the column names from the harmonized data
    harmonized_vars <- names(data)
    
    dictionary_data %>%
      mutate(location = if_else(variable_name %in% harmonized_vars, "harmonized", "REDCap"))
  })
  
  # Render the new data dictionary table
  output$dictionary_table <- renderDT({
    req(credentials()$user_auth) # Require authentication
    datatable(
      dictionary_with_location(),
      extensions = 'Buttons',
      options = list(
        pageLength = 10,
        dom = 'Bfrtip',
        buttons = list(
          list(extend = 'csv', text = 'Download Dictionary')
        )
      )
    )
  })
  
  # Server-side download handler for the descriptive table
  output$download_descriptive <- downloadHandler(
    filename = function() {
      paste("descriptive-stats-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      df_temp <- filtered_data()
      if (!"group" %in% names(df_temp)) {
        descriptive_summary <- data.frame("Note" = "Group column not found in filtered data.")
      } else {
        # NEW FIX: Handle NA group values
        df_temp$group[is.na(df_temp$group)] <- "NA Group"
        descriptive_summary <- df_temp %>%
          bind_rows(mutate(., group = "Grand Total")) %>%
          group_by(group) %>%
          summarise(
            `Record ID (N)` = length(unique(record_id)),
            `Age (avg)` = if("age" %in% names(df_temp)) mean(age, na.rm = TRUE) else NA_real_,
            `BMI (avg)` = if("bmi" %in% names(df_temp)) mean(bmi, na.rm = TRUE) else NA_real_,
            `HbA1c (avg)` = if("hba1c" %in% names(df_temp)) mean(hba1c, na.rm = TRUE) else NA_real_,
            `ACR-U (med)` = if("acr_u" %in% names(df_temp)) median(acr_u, na.rm = TRUE) else NA_real_,
            `Female (%)` = if("sex" %in% names(df_temp)) mean(sex == "Female", na.rm = TRUE) * 100 else NA_real_,
            `Hispanic (%)` = if("ethnicity" %in% names(df_temp)) mean(ethnicity == "Hispanic or Latino", na.rm = TRUE) * 100 else NA_real_,
            `eGFR_CKD_epi (avg)` = if("eGFR_CKD_epi" %in% names(df_temp)) mean(eGFR_CKD_epi, na.rm = TRUE) else NA_real_,
            `Total Kidney Volume (avg)` = if("total_kidney_volume_ml" %in% names(df_temp)) mean(total_kidney_volume_ml, na.rm = TRUE) else NA_real_,
            `Manual Total Kidney Volume (avg)` = if("total_kidney_volume_ml_manual" %in% names(df_temp)) mean(total_kidney_volume_ml_manual, na.rm = TRUE) else NA_real_,
            `ERPF Raw Plasma (avg)` = if("erpf raw plasma" %in% names(df_temp)) mean(`erpf raw plasma`, na.rm = TRUE) else NA_real_,
            `ERPF BSA Plasma (avg)` = if("erpf bsa plasma" %in% names(df_temp)) mean(`erpf bsa plasma`, na.rm = TRUE) else NA_real_,
            `PGLO (avg)` = if("pglo" %in% names(df_temp)) mean(pglo, na.rm = TRUE) else NA_real_,
            `RA (avg)` = if("ra" %in% names(df_temp)) mean(ra, na.rm = TRUE) else NA_real_,
            `RE (avg)` = if("re" %in% names(df_temp)) mean(re, na.rm = TRUE) else NA_real_,
            .groups = 'drop'
          ) %>%
          mutate(across(where(is.numeric), ~ round(., 2))) %>%
          rename(Category = group) %>%
          tidyr::pivot_longer(cols = -Category, names_to = "Metric", values_to = "Count") %>%
          tidyr::pivot_wider(names_from = Category, values_from = Count) %>%
          select(Metric, everything())
      }
      write_csv(descriptive_summary, file)
    }
  )
  
  # Server-side download handler for the availability table
  output$download_availability <- downloadHandler(
    filename = function() {
      if (!is.null(searched_variable_availability_by_study())) {
        paste("variable-availability-by-study-", Sys.Date(), ".csv", sep = "")
      } else if (!is.null(searched_variable_availability())) {
        paste("variable-availability-by-group-", Sys.Date(), ".csv", sep = "")
      } else {
        paste("data-availability-", Sys.Date(), ".csv", sep = "")
      }
    },
    content = function(file) {
      if (!is.null(searched_variable_availability_by_study())) {
        write_csv(searched_variable_availability_by_study(), file)
      } else if (!is.null(searched_variable_availability())) {
        write_csv(searched_variable_availability(), file)
      } else {
        write_csv(main_availability_table(), file)
      }
    }
  )
  
  # New download handler for the raw data table
  output$download_raw <- downloadHandler(
    filename = function() {
      paste("raw-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write_csv(raw_data_display_columns(), file)
    }
  )
  
  # Server-side download handler for the plot
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      # This requires the plot to be assigned to a variable within the renderPlot scope 
      # or by passing the plot object explicitly.
      # Since renderPlot returns the last expression (the plot object), we rely on that.
      # For a true ggsave, we'd need to re-generate the plot logic here, but for simplicity 
      # and Shiny's internal handling, we'll use the default plot device.
      # A simple workaround to ensure the plot is available for ggsave:
      plot_obj <- output$data_plot
      
      # Re-generate the plot object inside the handler to ensure it exists
      df <- filtered_data()
      p <- switch(input$plot_type,
                  "Bar Chart" = {
                    if (is.numeric(df[[input$plot_variable]])) return(NULL)
                    ggplot(df, aes(x = !!sym(input$plot_variable))) + geom_bar(fill = "#4f46e5") + theme_minimal() + labs(title = paste("Distribution of", input$plot_variable), x = input$plot_variable, y = "Count") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                  },
                  "Histogram" = {
                    if (!is.numeric(df[[input$plot_variable]])) return(NULL)
                    ggplot(df, aes(x = !!sym(input$plot_variable))) + geom_histogram(fill = "#10b981", color = "#1f2937", bins = 20) + labs(title = paste("Histogram of", input$plot_variable), x = input$plot_variable, y = "Frequency") + theme_minimal()
                  },
                  "Scatter Plot" = {
                    if (is.null(input$y_variable) || !is.numeric(df[[input$plot_variable]]) || !is.numeric(df[[input$y_variable]])) return(NULL)
                    ggplot(df, aes(x = !!sym(input$plot_variable), y = !!sym(input$y_variable))) + geom_point(color = "#f59e0b", alpha = 0.7) + labs(title = paste("Scatter Plot of", input$plot_variable, "vs", input$y_variable), x = input$plot_variable, y = input$y_variable) + theme_minimal()
                  },
                  NULL
      )
      
      if (!is.null(p)) {
        # Save the plot object to the file path
        ggsave(file, plot = p, device = "png", width = 10, height = 7)
      }
    }
  )
  
  # Reactive list of columns suitable for plotting
  plot_variables <- reactive({
    req(credentials()$user_auth) # Require authentication
    names(filtered_data()) %>%
      setdiff(c("record_id", "visit"))
  })
  
  # Dynamic UI for the x-axis variable
  output$plot_variable_ui <- renderUI({
    req(credentials()$user_auth) # Require authentication
    selectizeInput(
      "plot_variable",
      "Variable (X-Axis):",
      choices = plot_variables(),
      selected = NULL,
      options = list(placeholder = "Select a variable...", onInitialize = I('function() { this.setValue(""); }'))
    )
  })
  
  # Dynamic UI for the y-axis variable (for scatter plots)
  output$y_variable_ui <- renderUI({
    req(credentials()$user_auth) # Require authentication
    if (input$plot_type == "Scatter Plot") {
      selectizeInput(
        "y_variable",
        "Variable (Y-Axis):",
        choices = plot_variables(),
        selected = NULL,
        options = list(placeholder = "Select a variable...", onInitialize = I('function() { this.setValue(""); }'))
      )
    }
  })
  
  # Render the data visual based on user selections
  output$data_plot <- renderPlot({
    req(credentials()$user_auth) # Require authentication
    req(input$plot_variable, filtered_data())
    
    # Check if the selected variable exists in the filtered data
    if (!(input$plot_variable %in% names(filtered_data()))) {
      return(NULL)
    }
    
    # Check if the plot type and variable type are compatible
    df <- filtered_data()
    x_var_type <- class(df[[input$plot_variable]])
    
    output$plot_message <- renderText({""})
    
    # Handle plotting logic based on the selected type
    p <- switch(input$plot_type,
                "Bar Chart" = {
                  if (is.numeric(df[[input$plot_variable]])) {
                    output$plot_message <- renderText({"Please select a non-numeric variable for a Bar Chart."})
                    return(NULL)
                  }
                  ggplot(df, aes(x = !!sym(input$plot_variable))) +
                    geom_bar(fill = "#4f46e5") +
                    labs(title = paste("Distribution of", input$plot_variable), x = input$plot_variable, y = "Count") +
                    theme_minimal() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                },
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
                  req(input$y_variable)
                  if (!(input$y_variable %in% names(df))) {
                    return(NULL)
                  }
                  if (!is.numeric(df[[input$plot_variable]]) || !is.numeric(df[[input$y_variable]])) {
                    output$plot_message <- renderText({"Please select two numeric variables for a Scatter Plot."})
                    return(NULL)
                  }
                  ggplot(df, aes(x = !!sym(input$plot_variable), y = !!sym(input$y_variable))) +
                    geom_point(color = "#f59e0b", alpha = 0.7) +
                    labs(title = paste("Scatter Plot of", input$plot_variable, "vs", input$y_variable), x = input$plot_variable, y = input$y_variable) +
                    theme_minimal()
                },
                {
                  return(NULL)
                }
    )
    
    print(p)
  })
}

# --- 5. Run the App ---
shinyApp(ui = ui, server = server)
