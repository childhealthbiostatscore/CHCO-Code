# Load required libraries
library(shiny)
library(reactable)
library(htmltools)
library(dplyr)

# Create the GWAS summary statistics table for Shiny app
create_gwas_table <- function(data) {
  
  # Helper function to format p-values
  format_pvalue <- function(p) {
    case_when(
      p < 1e-10 ~ sprintf("%.2e", p),
      p < 0.001 ~ sprintf("%.1e", p),
      p < 0.01 ~ sprintf("%.4f", p),
      TRUE ~ sprintf("%.3f", p)
    )
  }
  
  # Helper function to get p-value color
  get_pvalue_color <- function(p) {
    case_when(
      p < 5e-8 ~ "#d73027",   # Genome-wide significant - red
      p < 1e-5 ~ "#fc8d59",   # Suggestive - orange  
      p < 0.05 ~ "#fee08b",   # Nominally significant - yellow
      TRUE ~ "#e0e0e0"        # Non-significant - gray
    )
  }
  
  # Create the reactable
  reactable(
    data,
    columns = list(
      phenotype = colDef(
        name = "Phenotype",
        minWidth = 200,
        cell = function(value, index) {
          if(grepl("^proteomics", data$analysis[index]) == TRUE)
          {
            div(title = value, 
                #style = "overflow: hidden; text-overflow: ellipsis; white-space: font-family: monospace; font-size: 11px;",
                tagList(
                  data$TargetFullName[index],
                  " (",
                  value,
                  "; ",
                  tags$a(href = paste("https://www.uniprot.org/uniprotkb", data$UniProt[index], "entry", sep = "/"), target = "_blank", data$UniProt[index]),
                  "; ",
                  tags$a(href = paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", data$EntrezGeneID[index]), target = "_blank", data$EntrezGeneSymbol[index]),
                  ")"
                ))
            
          }else
          {
            clean_name <- gsub("\\.", " ", value)
            clean_name <- gsub("_", " ", clean_name)
            
            div(title = clean_name, 
                style = "overflow: hidden; text-overflow: ellipsis; white-space: nowrap;",
                clean_name)
          }
        }
      ),
      
      analysis = colDef(
        name = "Analysis Type",
        maxWidth = 200,
        cell = function(value) {
          type <- case_when(
            grepl("metabolomics_plasma", value) ~ "Metabolomics (Plasma)",
            grepl("metabolomics_urine", value) ~ "Metabolomics (Urine)",
            grepl("proteomics", value) ~ "Proteomics",
            TRUE ~ "Other"
          )
          span(type, 
               style = paste0("padding: 2px 6px; border-radius: 3px; font-size: 11px; font-weight: bold; ",
                              if (type == "Metabolomics (Plasma)") "background-color: #AA4488; color: #ffffff;" 
                              else if (type == "Metabolomics (Urine)") "background-color: #F0027F; color: #ffffff;"
                              else if (type == "Proteomics") "background-color: #377EB8; color: #ffffff;"
                              else "background-color: #f5f5f5; color: #616161;"))
        }
      ),
      
      ID = colDef(
        maxWidth = 200,
        cell = function(value, index) {
          profile_url <- data$url_locus[index]
          tags$a(href   = profile_url, 
                 target = "_blank", 
                 style  = "font-family: monospace; ",
                 paste0(data$CHROM[index], ":", value, "-", data$END[index]))
        }
      ),
      
      BETA = colDef(
        name = "Effect Size (±SE)",
        maxWidth = 200,
        cell = function(value, index) {
          beta <- data$BETA[index]
          se <- data$SE[index]
          span(paste(beta, se, sep = "±"),
               style  = "font-family: monospace;",)
        }
      ),
      
      P = colDef(
        name = "P-value",
        maxWidth = 200,
        align = "right",
        cell = function(value) {
          color <- get_pvalue_color(value)
          formatted_p <- format_pvalue(value)
          
          span(formatted_p,
               style = sprintf("color: %s; font-weight: bold; padding: 2px 4px; border-radius: 2px; background-color: %s;",
                               if (value < 0.05) "white" else "#333",
                               if (value < 0.05) color else "transparent"))
        },
        sortNALast = TRUE
      ),
      
      START = colDef(
        name = "Locus",
        maxWidth = 200,
        cell = function(value, index) {
          profile_url <- data$url_snp[index]
          tags$a(href   = profile_url, 
                 target = "_blank", 
                 style  = "font-family: monospace; ",
                 value)
        }
      ),

      # Hide technical columns that aren't needed for display
      TargetFullName   = colDef(show = FALSE),
      Target           = colDef(show = FALSE),
      sumstats_file    = colDef(show = FALSE),
      UniProt          = colDef(show = FALSE),
      EntrezGeneID     = colDef(show = FALSE),
      EntrezGeneSymbol = colDef(show = FALSE),
      CHROM            = colDef(show = FALSE),
      POS              = colDef(show = FALSE),
      REF              = colDef(show = FALSE),
      ALT              = colDef(show = FALSE),
      FDR              = colDef(show = FALSE),
      url_snp          = colDef(show = FALSE),
      url_locus        = colDef(show = FALSE),
      END              = colDef(show = FALSE),
      N_VARIANTS       = colDef(show = FALSE),
      SIGNAL_ID        = colDef(show = FALSE),
      SE               = colDef(show = FALSE)  # Hide SE column since it's shown in the bar plot
    ),
    
    # Table options
    defaultSorted = list(P = "asc"),
    defaultPageSize = 20,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(10, 20, 50, 100),
    searchable = TRUE,
    filterable = TRUE,
    resizable = TRUE,
    showSortable = TRUE,
    highlight = TRUE,
    striped = TRUE,
    
    # Theme and styling
    theme = reactableTheme(
      headerStyle = list(
        background = "#f8f9fa",
        borderColor = "#dee2e6",
        fontSize = "12px",
        fontWeight = "bold"
      ),
      cellStyle = list(fontSize = "11px"),
      searchInputStyle = list(width = "300px"),
      pageButtonStyle = list(fontSize = "12px")
    )
    
  )
}
