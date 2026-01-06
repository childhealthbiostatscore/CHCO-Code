# Install ggsignif if needed
# install.packages("ggsignif")
library(ggsignif)
library(patchwork)

# FSOC by Group with ALL pairwise p-values on plot
plot_fsoc_by_group <- function(dat) {
  
  fsoc_long <- dat %>%
    filter(!is.na(group)) %>%
    dplyr::select(record_id, group, sex,
                  whole_kidney_fsoc_abs, medullary_fsoc_abs) %>%
    pivot_longer(cols = contains("fsoc"),
                 names_to = "fsoc_type",
                 values_to = "fsoc_value") %>%
    filter(!is.na(fsoc_value)) %>%
    filter(fsoc_value >= 0) %>%
    mutate(fsoc_type = case_match(fsoc_type,
                                  "whole_kidney_fsoc_abs" ~ "Whole Kidney FSOC",
                                  "medullary_fsoc_abs" ~ "Medullary FSOC"))
  
  # Get unique groups (excluding PKD if no data)
  groups <- fsoc_long %>% 
    filter(fsoc_type == "Whole Kidney FSOC") %>%
    pull(group) %>% 
    unique() %>%
    as.character()
  
  # Create all pairwise comparisons
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  # Calculate pairwise p-values
  calc_pairwise_p <- function(data, fsoc_name) {
    subset_data <- data %>% filter(fsoc_type == fsoc_name)
    p_values <- sapply(comparisons, function(pair) {
      g1 <- subset_data %>% filter(group == pair[1]) %>% pull(fsoc_value)
      g2 <- subset_data %>% filter(group == pair[2]) %>% pull(fsoc_value)
      if (length(g1) > 1 & length(g2) > 1) {
        wilcox.test(g1, g2)$p.value
      } else {
        NA
      }
    })
    return(p_values)
  }
  
  wk_pvals <- calc_pairwise_p(fsoc_long, "Whole Kidney FSOC")
  med_pvals <- calc_pairwise_p(fsoc_long, "Medullary FSOC")
  
  # Format p-values for plot
  format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("p<0.001")
    if (p < 0.01) return(sprintf("p=%.3f", p))
    return(sprintf("p=%.2f", p))
  }
  
  # Print to console
  cat("\n--- Pairwise Wilcoxon Tests (Raw p-values) ---\n")
  cat("WHOLE KIDNEY:\n")
  for (i in seq_along(comparisons)) {
    cat(sprintf("  %s vs %s: %s\n", comparisons[[i]][1], comparisons[[i]][2], format_p(wk_pvals[i])))
  }
  cat("\nMEDULLARY:\n")
  for (i in seq_along(comparisons)) {
    cat(sprintf("  %s vs %s: %s\n", comparisons[[i]][1], comparisons[[i]][2], format_p(med_pvals[i])))
  }
  
  # Whole Kidney plot
  wk_data <- fsoc_long %>% filter(fsoc_type == "Whole Kidney FSOC")
  
  p_wk <- ggplot(wk_data, aes(x = group, y = fsoc_value, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    geom_signif(
      comparisons = comparisons,
      annotations = sapply(wk_pvals, format_p),
      step_increase = 0.12,
      tip_length = 0.02,
      textsize = 2.8,
      vjust = 0.3,
      map_signif_level = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    ) +
    labs(x = "Disease Group", 
         y = "FSOC Value (s⁻¹)",
         title = "Whole Kidney FSOC") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  # Medullary plot
  med_data <- fsoc_long %>% filter(fsoc_type == "Medullary FSOC")
  
  p_med <- ggplot(med_data, aes(x = group, y = fsoc_value, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    geom_signif(
      comparisons = comparisons,
      annotations = sapply(med_pvals, format_p),
      step_increase = 0.12,
      tip_length = 0.02,
      textsize = 2.8,
      vjust = 0.3,
      map_signif_level = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    ) +
    labs(x = "Disease Group", 
         y = "FSOC Value (s⁻¹)",
         title = "Medullary FSOC") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  # Combine
  combined <- p_wk + p_med +
    plot_annotation(
      title = "FSOC Endpoints by Disease Group",
      subtitle = "Raw pairwise Wilcoxon p-values; negative FSOC values excluded",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, face = "italic", hjust = 0.5)
      )
    )
  
  return(combined)
}

# Run it
p_group <- plot_fsoc_by_group(dat)
print(p_group)

# Save
ggsave(file.path(OUTPUT_DIR, "fsoc_by_group_pairwise.pdf"), p_group, width = 14, height = 8)
ggsave(file.path(OUTPUT_DIR, "fsoc_by_group_pairwise.tiff"), p_group, width = 14, height = 8, dpi = 300)