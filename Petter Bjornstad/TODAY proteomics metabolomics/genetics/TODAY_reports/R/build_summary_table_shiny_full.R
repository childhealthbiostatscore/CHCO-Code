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
  
  # Helper function to create effect size bar plot
  create_effect_bar <- function(beta, se) {
    # Normalize bar width based on absolute beta values
    max_beta <- max(abs(beta), na.rm = TRUE)
    
    bars <- purrr::map2(beta, se, function(b, s) {
      if (is.na(b) | is.na(s)) return("")
      
      # Calculate bar width as percentage
      width_pct <- abs(b) / max_beta * 100
      
      # Color based on effect direction
      bar_color <- if (b > 0) "#2166ac" else "#d73027"  # Blue for positive, red for negative
      
      # Create bar with error indication
      error_width <- (s / max_beta) * 100
      
      div(
        style = sprintf(
          "display: flex; align-items: center; height: 20px; width: 100px;"
        ),
        div(
          style = sprintf(
            "background-color: %s; height: 12px; width: %.1f%%; margin-right: 2px; position: relative;",
            bar_color, width_pct
          ),
          # Add error bar indication
          div(
            style = sprintf(
              "position: absolute; right: -%.1fpx; top: 0; height: 12px; width: 2px; background-color: black; opacity: 0.6;",
              error_width * (width_pct / 100)
            )
          )
        ),
        span(sprintf("%.3f±%.3f", b, s), style = "font-size: 10px; margin-left: 4px;")
      )
    })
    
    return(bars)
  }
  
  # Create the reactable
  reactable(
    data,
    columns = list(
      phenotype = colDef(
        name = "Phenotype",
        width = 200,
        cell = function(value) {
          # Clean up phenotype names
          clean_name <- gsub("\\.", " ", value)
          clean_name <- gsub("_", " ", clean_name)
          div(title = clean_name, 
              style = "overflow: hidden; text-overflow: ellipsis; white-space: nowrap;",
              clean_name)
        }
      ),
      
      TargetFullName = colDef(
        name = "Target Name",
        width = 180,
        cell = function(value) {
          if (is.na(value)) return("—")
          div(title = value,
              style = "overflow: hidden; text-overflow: ellipsis; white-space: nowrap;",
              value)
        }
      ),
      
      Target = colDef(
        name = "Symbol",
        width = 80,
        cell = function(value) if (is.na(value)) "—" else value
      ),
      
      analysis = colDef(
        name = "Analysis Type",
        width = 120,
        cell = function(value) {
          type <- case_when(
            grepl("metabolomics_plasma", value) ~ "Metabolomics (Plasma)",
            grepl("metabolomics_urine", value) ~ "Metabolomics (Urine)",
            grepl("proteomics", value) ~ "Proteomics",
            TRUE ~ "Other"
          )
          span(type, 
               style = paste0("padding: 2px 6px; border-radius: 3px; font-size: 11px; font-weight: bold; ",
                              if (type == "Metabolomics") "background-color: #e1f5fe; color: #0277bd;" 
                              else if (type == "Proteomics") "background-color: #f3e5f5; color: #7b1fa2;"
                              else "background-color: #f5f5f5; color: #616161;"))
        }
      ),
      
      CHROM = colDef(
        name = "Chr",
        width = 50,
        align = "center"
      ),
      
      POS = colDef(
        name = "Position",
        width = 100,
        align = "right",
        cell = function(value) format(value, big.mark = ",")
      ),
      
      ID = colDef(
        name = "Variant ID",
        width = 140,
        cell = function(value) {
          span(value, style = "font-family: monospace; font-size: 11px;")
        }
      ),
      
      REF = colDef(name = "Ref", width = 50, align = "center"),
      ALT = colDef(name = "Alt", width = 50, align = "center"),
      
      BETA = colDef(
        name = "Effect Size (±SE)",
        width = 150,
        cell = function(value, index) {
          beta <- data$BETA[index]
          se <- data$SE[index]
          create_effect_bar(beta, se)[[1]]
        }
      ),
      
      SE = colDef(show = FALSE),  # Hide SE column since it's shown in the bar plot
      
      P = colDef(
        name = "P-value",
        width = 100,
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
      
      FDR = colDef(
        name = "FDR",
        width = 100,
        align = "right",
        cell = function(value) {
          if (is.na(value)) return("—")
          formatted_fdr <- format_pvalue(value)
          
          color <- case_when(
            value < 0.001 ~ "#d73027",
            value < 0.01 ~ "#fc8d59", 
            value < 0.05 ~ "#fee08b",
            TRUE ~ "#666"
          )
          
          span(formatted_fdr, style = sprintf("color: %s; font-weight: %s;",
                                              color, if (value < 0.05) "bold" else "normal"))
        }
      ),
      
      # Hide technical columns that aren't needed for display
      sumstats_file = colDef(show = FALSE),
      UniProt = colDef(show = FALSE),
      EntrezGeneID = colDef(show = FALSE),
      EntrezGeneSymbol = colDef(show = FALSE)
    ),
    
    # Table options
    defaultSorted = list(P = "asc"),
    defaultPageSize = 20,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(10, 20, 50, 100),
    searchable = TRUE,
    filterable = FALSE,
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
