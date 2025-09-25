get_pvalue_color_gwas <- function(p) {
  case_when(
    p < 5e-8 ~ "#d73027",   # Genome-wide significant - red
    p < 1e-5 ~ "#fc8d59",   # Suggestive - orange  
    p < 0.05 ~ "#fee08b",   # Nominally significant - yellow
    TRUE ~ "#e0e0e0"        # Non-significant - gray
  )
}

get_pvalue_color <- function(p, min_p = NULL, n = 100,
                             palette = "Reds") {
  # Clip values to (0,1]
  p <- pmax(p, .Machine$double.xmin)
  p <- pmin(p, 1)
  
  # Transform to -log10 scale
  lp <- -log10(p)
  
  # Define the range: p=1 (lp=0) → strongest signal
  if (is.null(min_p)) {
    min_p <- min(p[p < 0.05], 1e-6, na.rm = TRUE)
  }
  max_lp <- -log10(min_p)
  
  # Scale to [0,1]
  scaled <- lp / max_lp
  scaled[scaled > 1] <- 1
  
  # Get a sequential palette from colorspace
  pal <- sequential_hcl(n, palette, rev = TRUE)
  
  # Map scaled values to colors
  pal[ceiling(scaled * (n - 1)) + 1]
}

convert_phenotype_to_full_name = function(value, index, data)
{
  if(grepl("^proteomics", data$analysis[index]) == TRUE)
  {
    div(title = value, 
        tagList(
          data$TargetFullName[index],
          " (",
          value,
          "; ",
          data$Target[index],
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

color_analysis = function(value)
{
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

convert_id_to_url = function(value, index, data)
{
  profile_url <- data$url_locus[index]
  tags$a(href   = profile_url, 
         target = "_blank", 
         style  = "font-family: monospace; ",
         value)
}

combine_beta_se = function(value, index, data)
{
  beta <- signif(data$BETA[index], 3)
  se   <- signif(data$SE[index], 3)
  span(paste(beta, se, sep = "±"),
       style  = "font-family: monospace;")
}

annotate_pval = function(value, ...)
{
  color <- get_pvalue_color(value, ...)
  formatted_p <- format(value, scientific = TRUE, digits = 3)
  
  span(formatted_p,
       style = sprintf("color: %s; font-weight: bold; padding: 2px 4px; border-radius: 2px; background-color: %s;",
                       if (value < 0.05) "white" else "#333",
                       if (value < 0.05) color else "transparent"))
  
}

annotate_model = function(value)
{
  my_color = case_when(
    value == "m1" ~ "#1B263B",
    value == "m2" ~ "#2D6A4F",
    value == "m3" ~ "#7B2CBF",
    value == "m4" ~ "#8B0000",
    TRUE ~ "#000000"
  )
  span(value,
       style = sprintf("color: %s; font-weight: bold; padding: 2px 4px; border-radius: 2px; background-color: %s;",
                       "white", 
                       my_color))
  
}

annotate_term = function(value)
{
  value = sub("gt", "genotype", value)
  my_color = case_when(
    value == "genotype" ~ "#AEC6CF",
    value == "value" ~ "#B5EAD7",
    value == "genotype:value" ~ "#E6A8D7",
    TRUE ~ "#ffffff"
  )
  span(value,
       style = sprintf("color: %s; font-weight: bold; padding: 2px 4px; border-radius: 2px; background-color: %s;",
                       "black", 
                       my_color))
  
}

get_auc_color <- function(x, n = 100, palette = "Reds") 
{
  x = ifelse(x < 0.5, yes = 0, no = (x - 0.5) * 2)
  

  # Get a sequential palette from colorspace
  pal <- sequential_hcl(n, palette, rev = TRUE)
  
  # Map scaled values to colors
  pal[ceiling(x * (n - 1)) + 1]
}


annotate_auc_color = function(value, ...)
{
  color <- get_auc_color(value, ...)
  formatted_auc <- signif(value, digits = 3)
  
  span(formatted_auc,
       style = sprintf("color: %s; font-weight: bold; padding: 2px 4px; border-radius: 2px; background-color: %s;",
                       "#333",
                       color))
}


################################################################################
# Create full summary statistics table
create_full_sumstats_table = function(data)
{
  #data = data[1:10,]
  
  reactable(
    data    = data,
    columns = list(
      phenotype = colDef(
        name     = "Phenotype",
        minWidth = 200,
        cell     = function(value, index) {convert_phenotype_to_full_name(value, index, data)}
      ),
      
      analysis = colDef(
        name = "Analysis Type",
        maxWidth = 200,
        cell = function(value) {color_analysis(value)}
      ),
      
      ID = colDef(
        name     = "ID",
        minWidth = 200,
        cell     = function(value, index) {convert_id_to_url(value, index, data)}
      ),
      
      BETA = colDef(
        name     = "Effect Size (±SE)",
        minWidth = 200,
        cell     = function(value, index) {combine_beta_se(value, index, data)}
      ),
      
      P = colDef(
        name = "P-value",
        maxWidth = 200,
        align = "right",
        cell = function(value) {annotate_pval(value, palette = "Rocket", min_p = 1e-100)},
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
                 paste0(data$CHROM[index], ":", value, "-", data$END[index]))
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
      file_model       = colDef(show = FALSE),
      complete         = colDef(show = FALSE),
      pheno_name       = colDef(show = FALSE),
      SE               = colDef(show = FALSE)  # Hide SE column since it's shown in the bar plot
    ),
    
    # Table options
    defaultSorted       = list(P = "asc"),
    defaultPageSize     = 20,
    showPageSizeOptions = TRUE,
    pageSizeOptions     = c(10, 20, 50, 100),
    searchable          = TRUE,
    filterable          = TRUE,
    resizable           = TRUE,
    showSortable        = TRUE,
    highlight           = TRUE,
    striped             = TRUE,
    
    # Theme and styling
    theme = reactableTheme(
      headerStyle = list(
        background = "#f8f9fa",
        borderColor = "#dee2e6",
        fontSize = "12px",
        fontWeight = "bold"),
      cellStyle = list(fontSize = "11px"),
      searchInputStyle = list(width = "300px"),
      pageButtonStyle = list(fontSize = "12px")
    )
  )
}

################################################################################
# Create full models summary statistics table
create_full_models_table = function(data)
{
  #data = data[1:10,]
  data = data %>% rename(BETA = estimate, 
                         SE   = std_error,
                         P    = p_value)
  
  reactable(
    data = data,
    columns = list(
      phenotype = colDef(
        name     = "Phenotype",
        minWidth = 200,
        cell     = function(value, index) {convert_phenotype_to_full_name(value, index, data)}
      ),
      
      analysis = colDef(
        name = "Analysis Type",
        maxWidth = 200,
        cell = function(value) {color_analysis(value)}
      ),
      
      ID = colDef(
        name     = "ID",
        minWidth = 200
      ),
      
      outcome   = colDef(
        name = "Outcome",
      ),
      
      model = colDef(
        name = "Model",
        maxWidth = 200,
        cell = function(value) {annotate_model(value)}
      ),
      
      term = colDef(
        name = "Term",
        maxWidth = 200,
        cell = function(value) {annotate_term(value)}
      ),
      
      BETA = colDef(
        name     = "Effect Size (±SE)",
        minWidth = 200,
        cell     = function(value, index) {combine_beta_se(value, index, data)}
      ),
      
      P = colDef(
        name = "P-value (GLM)",
        maxWidth = 200,
        align = "right",
        cell = function(value) {annotate_pval(value, palette = "Purp", min_p = 1e-10)},
        sortNALast = TRUE
      ),
      
      auc = colDef(
        name = "ROC AUC",
        maxWidth = 200,
        cell = function(value) {annotate_auc_color(value, palette = "Peach")}
      ),
      
      aic = colDef(
        name = "AIC",
        maxWidth = 200,
        cell = function(value) {trunc(value)}
      ),
      
      bic = colDef(
        name = "BIC",
        maxWidth = 200,
        cell = function(value) {trunc(value)}
      ),
      
      lr_pvalue = colDef(
        name = "P-value (LRT)",
        maxWidth = 200,
        align = "right",
        cell = function(value) {annotate_pval(value, palette = "RdPu", min_p = 1e-10)},
        sortNALast = TRUE
      ),
      
      # Hide technical columns that aren't needed for display
      TargetFullName   = colDef(show = FALSE),
      Target           = colDef(show = FALSE),
      UniProt          = colDef(show = FALSE),
      EntrezGeneID     = colDef(show = FALSE),
      EntrezGeneSymbol = colDef(show = FALSE),
      file_model       = colDef(show = FALSE),
      SE               = colDef(show = FALSE)
    ),
    
    # Table options
    defaultSorted       = list(auc = "desc"),
    defaultPageSize     = 20,
    showPageSizeOptions = TRUE,
    pageSizeOptions     = c(10, 20, 50, 100),
    searchable          = TRUE,
    filterable          = TRUE,
    resizable           = TRUE,
    showSortable        = TRUE,
    highlight           = TRUE,
    striped             = TRUE,
    
    # Theme and styling
    theme = reactableTheme(
      headerStyle = list(
        background = "#f8f9fa",
        borderColor = "#dee2e6",
        fontSize = "12px",
        fontWeight = "bold"),
      cellStyle = list(fontSize = "11px"),
      searchInputStyle = list(width = "300px"),
      pageButtonStyle = list(fontSize = "12px")
    )
  )
}


################################################################################
# Create summary table filtered by phenotype
create_summary_table = function(data)
{
  data = data[1:10,]
  
  reactable(
    data = data,
    columns = list(
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
    defaultSorted       = list(P = "asc"),
    defaultPageSize     = 20,
    showPageSizeOptions = TRUE,
    pageSizeOptions     = c(10, 20, 50, 100),
    searchable          = TRUE,
    filterable          = TRUE,
    resizable           = TRUE,
    showSortable        = TRUE,
    highlight           = TRUE,
    striped             = TRUE,
    
    # Theme and styling
    theme = reactableTheme(
      headerStyle = list(
        background = "#f8f9fa",
        borderColor = "#dee2e6",
        fontSize = "12px",
        fontWeight = "bold"),
      cellStyle = list(fontSize = "11px"),
      searchInputStyle = list(width = "300px"),
      pageButtonStyle = list(fontSize = "12px")
    )
  )
}


################################################################################
# Create model table per phenotype
create_variant_table = function(data)
{
  data = data[1:10,]
  
  reactable(
    data = data,
    columns = list(
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
    defaultSorted       = list(P = "asc"),
    defaultPageSize     = 20,
    showPageSizeOptions = TRUE,
    pageSizeOptions     = c(10, 20, 50, 100),
    searchable          = TRUE,
    filterable          = TRUE,
    resizable           = TRUE,
    showSortable        = TRUE,
    highlight           = TRUE,
    striped             = TRUE,
    
    # Theme and styling
    theme = reactableTheme(
      headerStyle = list(
        background = "#f8f9fa",
        borderColor = "#dee2e6",
        fontSize = "12px",
        fontWeight = "bold"),
      cellStyle = list(fontSize = "11px"),
      searchInputStyle = list(width = "300px"),
      pageButtonStyle = list(fontSize = "12px")
    )
  )
}






