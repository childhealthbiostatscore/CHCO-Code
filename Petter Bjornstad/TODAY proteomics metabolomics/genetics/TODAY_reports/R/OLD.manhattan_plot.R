create_manhattan_plot <- function(data,
                                  significance_threshold = 5e-8,
                                  suggestive_threshold = 1e-5,
                                  title = "Manhattan Plot",
                                  point_size = 3) {
  data = data %>% 
    filter(P < 1e-4 & (CHROM %in% 1:22 | CHROM %in% as.character(1:22))) %>%
    mutate(LOG10P = -log10(P))
  
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
  
  # Create alternating colors for chromosomes
  data$color <- ifelse(data$CHROM %% 2 == 0, "Even", "Odd")
  
  # Define significance lines
  sig_line <- -log10(significance_threshold)
  sug_line <- -log10(suggestive_threshold)
  
  p <- plot_ly(data, 
               x = ~BPcum, 
               y = ~LOG10P,
               color = ~color,
               colors = c("Even" = "#2E86AB", "Odd" = "#A23B72"),
               text = ~paste("SNP:", ID,
                             "<br>CHROM:", CHROM,
                             "<br>POS:", format(POS, big.mark = ","),
                             "<br>P-value:", format(P, scientific = TRUE, digits = 3)),
               hovertemplate = "%{text}<extra></extra>",
               type = "scatter",
               mode = "markers",
               marker = list(size = point_size, opacity = 0.7)) %>%
    
    # Add significance threshold line
    # Add significance threshold line
    add_segments(x = min(data$BPcum), xend = max(data$BPcum), 
                 y = sig_line, yend = sig_line,
              line = list(color = "red", dash = "dash", width = 2),
              showlegend = FALSE,
              hoverinfo = "skip") %>%
    
    # Add suggestive threshold line
    add_segments(x = min(data$BPcum), xend = max(data$BPcum), 
              y = sug_line, yend = sug_line,
              line = list(color = "blue", dash = "dot", width = 2),
              showlegend = FALSE,
              hoverinfo = "skip") %>%
    
    # Customize layout
    layout(
      title = list(text = title, font = list(size = 16)),
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
        zeroline = FALSE
      ),
      showlegend = FALSE,
      hovermode = "closest",
      plot_bgcolor = "white",
      paper_bgcolor = "white"
    )  
  
  return(p)
}
  