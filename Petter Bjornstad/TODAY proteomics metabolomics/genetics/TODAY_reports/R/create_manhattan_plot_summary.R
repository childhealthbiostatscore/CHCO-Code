create_manhattan_plot = function(data)
{
  point_size = 8
  
  data = data %>% 
    mutate(logp  = -log10(P),
           analysis = case_when(
             grepl("metabolomics_plasma", analysis) ~ "Metabolomics (Plasma)",
             grepl("metabolomics_urine" , analysis) ~ "Metabolomics (Urine)",
             grepl("proteomics"         , analysis) ~ "Proteomics",
             TRUE ~ "Other")) %>%
    mutate(my_color = case_when(
             analysis == "Metabolomics (Plasma)" ~ "#AA4488",
             analysis == "Metabolomics (Urine)"  ~ "#F0027F",
             analysis == "Proteomics"            ~ "#377EB8",
             TRUE ~ "#e0e0e0"),
           text  = ifelse(analysis == "Proteomics", 
                          yes = paste("<br>Full name:"     , TargetFullName, 
                                      "<br>UniProt:"       , UniProt, 
                                      "<br>Entrez Gene ID:", EntrezGeneID, 
                                      "<br>Gene Symbol:"   , EntrezGeneSymbol),
           no  = paste("<br>Full name:", gsub("_", " ", gsub("\\.", " ", phenotype))))
           )

  # Create cumulative base pair positions for plotting
  data <- data %>%
    arrange(CHROM, POS) %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(POS), .groups = 'drop') %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    left_join(data, ., by = "CHROM") %>%
    arrange(CHROM, POS) %>%
    mutate(BPcum = POS + tot)
  # Get chromosome center positions for x-axis labels
  axisdf <- data %>%
    group_by(CHROM) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = 'drop')

  
  p <- plot_ly(data, 
               x = ~BPcum, 
               y = ~logp,
               color = ~I(my_color),
               text = ~paste("SNP:", ID,
                             text,
                             "<br>P-value:", format(P, scientific = TRUE, digits = 3)),
               hovertemplate = "%{text}<extra></extra>",
               type = "scatter",
               mode = "markers",
               marker = list(size = point_size, opacity = 0.7)) %>%

  layout(
    xaxis = list(
      title = "Chromosome",
      tickmode = "array",
      tickvals = axisdf$center,
      ticktext = axisdf$CHROM,
      showgrid = FALSE,
      zeroline = FALSE
    ),
    yaxis = list(
      title = "-log₁₀(P-value)",
      showgrid = TRUE,
      zeroline = TRUE
    ),
    
    showlegend = FALSE,
    hovermode = "closest",
    plot_bgcolor = "white",
    paper_bgcolor = "white"
  )  

  return(p)
}

################################################################################
# Volcano plot
create_volcano_plot = function(data)
{
  point_size = 8
  
  data = data %>% 
    mutate(logp  = -log10(P),
           analysis = case_when(
             grepl("metabolomics_plasma", analysis) ~ "Metabolomics (Plasma)",
             grepl("metabolomics_urine" , analysis) ~ "Metabolomics (Urine)",
             grepl("proteomics"         , analysis) ~ "Proteomics",
             TRUE ~ "Other")) %>%
    mutate(my_color = case_when(
      analysis == "Metabolomics (Plasma)" ~ "#AA4488",
      analysis == "Metabolomics (Urine)"  ~ "#F0027F",
      analysis == "Proteomics"            ~ "#377EB8",
      TRUE ~ "#e0e0e0"),
      text  = ifelse(analysis == "Proteomics", 
                     yes = paste("<br>Full name:"     , TargetFullName, 
                                 "<br>UniProt:"       , UniProt, 
                                 "<br>Entrez Gene ID:", EntrezGeneID, 
                                 "<br>Gene Symbol:"   , EntrezGeneSymbol),
                     no  = paste("<br>Full name:", gsub("_", " ", gsub("\\.", " ", phenotype))))
    )
  
  p <- plot_ly(data, 
               x = ~BETA, 
               y = ~logp,
               color = ~I(my_color),
               text = ~paste("SNP:", ID,
                             text,
                             "<br>P-value:", format(P, scientific = TRUE, digits = 3)),
               hovertemplate = "%{text}<extra></extra>",
               type = "scatter",
               mode = "markers",
               marker = list(size = point_size, opacity = 0.7)) %>%
    
    layout(
      xaxis = list(
        title = "Effect size",
        showgrid = FALSE,
        zeroline = TRUE
      ),
      yaxis = list(
        title = "-log₁₀(P-value)",
        showgrid = TRUE,
        zeroline = TRUE
      ),
      
      showlegend = FALSE,
      hovermode = "closest",
      plot_bgcolor = "white",
      paper_bgcolor = "white"
    )  
  
  return(p)
}