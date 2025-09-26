#Libraries
library(ggrepel)
# library(qpcR)
library(ggpubr)
library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)

#Directories
user <- "hhampson" 
dir.dat <- c(paste0("/Users/",user,"/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Liver project/IPA Results"))
dir.results <- c(paste0("/Users/",user,"/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Liver project/NEBULA Results"))

color_table <- tibble(
  up_down = c(-1, 0, 1),
  color = c("indianred", "gray", "steelblue")
)

ipa_plot <- function(data){
  data <- data %>% mutate(up_down = case_when(
    abs(`z-score`) < 0.0001 | is.na(`z-score`) ~ 0,
    `z-score` <= -0.0001 ~ -1,
    `z-score` >= 0.0001 ~ 1,
    .default = 0
  ))
  data$up_down <- factor(data$up_down, levels = unique(data$up_down))
  p <- ggplot(data, aes(reorder(`Ingenuity Canonical Pathways`, `-log(p-value)`), `-log(p-value)`, fill = up_down))+
    geom_col() +
    geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "black") +
    theme_classic()+
    # scale_y_continuous(limits = c(0, 6)) +
    ylab("-log(p-value)") +
    xlab("Pathway") + 
    # theme(legend.position="none") +
    scale_fill_manual(values = c("-1" = "steelblue", "0" = "grey", "1" = "indianred")) 
  p <- p + coord_flip()
}

#TODAY
##AST
ast <- readxl::read_xls(fs::path(dir.dat,"Today_AST_log2_ipa.xls"), skip = 1)

ast <- ast %>% arrange(desc("-log(p-value)")) 
ast_keep <- ast[1:40,]

ast_plot <- ipa_plot(ast_keep)
ast_plot 
png(fs::path(dir.results,"AST_pathways_TeenLabs.png"),width=3000,height=2000,res=300)
plot(ast_plot)
dev.off()

##ALT
alt <- readxl::read_xls(fs::path(dir.dat,"Today_ALT_log2_ipa.xls"), skip = 1)

alt <- alt %>% arrange(desc("-log(p-value)")) 
alt_keep <- alt[1:40,]

alt_plot <- ipa_plot(alt_keep)
alt_plot 
png(fs::path(dir.results,"ALT_pathways_TeenLabs.png"),width=3000,height=2000,res=300)
plot(alt_plot)
dev.off()

##Ratio
ratio <- readxl::read_xls(fs::path(dir.dat,"Today_AST_ALT_ratio_log2_ipa.xls"), skip = 1)

ratio <- ratio %>% arrange(desc("-log(p-value)")) 
ratio_keep <- ratio[1:40,]

ratio_plot <- ipa_plot(ratio_keep)
ratio_plot 
png(fs::path(dir.results,"AST_ALT_Ratio_pathways_TeenLabs.png"),width=3000,height=2000,res=300)
plot(ratio_plot)
dev.off()

#Try alternative formatting
library(tidyverse)
library(ggplot2)
library(ggtext)
library(scales)

# Enhanced IPA plotting function
ipa_plot <- function(data, 
                     title = NULL,
                     subtitle = NULL,
                     n_pathways = 40,
                     pval_cutoff = 0.05,
                     zscore_cutoff = 2) {
  
  # Prepare data
  data <- data %>%
    # Calculate direction based on z-score
    mutate(
      up_down = case_when(
        is.na(`z-score`) | abs(`z-score`) < 0.01 ~ "Not Significant",
        `z-score` > 0 ~ "Activated",
        `z-score` < 0 ~ "Inhibited"
      ),
      # Create significance indicator
      sig = ifelse(`-log(p-value)` > -log(pval_cutoff), "*", ""),
      # Truncate long pathway names for display
      pathway_label = str_wrap(`Ingenuity Canonical Pathways`, width = 40),
      pathway_short = ifelse(
        nchar(`Ingenuity Canonical Pathways`) > 45,
        paste0(substr(`Ingenuity Canonical Pathways`, 1, 42), "..."),
        `Ingenuity Canonical Pathways`
      )
    ) %>%
    # Keep top pathways by p-value
    arrange(desc(`-log(p-value)`)) %>%
    slice_head(n = n_pathways) %>%
    # Reorder by z-score for plotting
    mutate(pathway_label = fct_reorder(pathway_label, `z-score`, .na_rm = TRUE))
  
  # Set factor levels for consistent coloring
  data$up_down <- factor(data$up_down, 
                         levels = c("Inhibited", "Not Significant", "Activated"))
  
  # Create the plot
  p <- ggplot(data, aes(x = pathway_label, y = `z-score`, fill = up_down)) +
    # Add bars
    geom_col(width = 0.7, alpha = 0.9) +
    
    # Add significance stars
    geom_text(aes(label = sig, 
                  y = `z-score` + sign(`z-score`) * 0.1),
              size = 5, color = "black") +
    
    # Add reference lines
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-zscore_cutoff, zscore_cutoff), 
               linetype = "dashed", color = "gray40", alpha = 0.5) +
    
    # Styling
    scale_fill_manual(
      values = c("Inhibited" = "#3B82C4",      # Nice blue
                 "Not Significant" = "#9CA3AF", # Gray
                 "Activated" = "#DC2626"),      # Nice red
      name = "Pathway Status"
    ) +
    
    scale_y_continuous(
      breaks = pretty_breaks(n = 6),
      limits = c(min(c(-3, min(data$`z-score`, na.rm = TRUE) * 1.1)),
                 max(c(3, max(data$`z-score`, na.rm = TRUE) * 1.1)))
    ) +
    
    # Labels
    labs(
      x = NULL,
      y = "Activation z-score",
      title = title %||% "IPA Canonical Pathway Analysis",
      subtitle = subtitle %||% paste0("Top ", n_pathways, " pathways by significance"),
      caption = paste0("* p < ", pval_cutoff, " | Dashed lines at z = ±", zscore_cutoff)
    ) +
    
    # Theme
    theme_minimal(base_size = 11) +
    theme(
      # Title and labels
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 1),
      
      # Axes
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11, face = "bold"),
      
      # Legend
      legend.position = "top",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.7, "cm"),
      
      # Grid
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      
      # Background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "gray98", color = NA)
    ) +
    
    # Flip coordinates for horizontal bars
    coord_flip()
  
  return(p)
}

# Alternative version with pathway names on bars
ipa_plot_labeled <- function(data, 
                             title = NULL,
                             subtitle = NULL,
                             n_pathways = 20,  # Fewer for cleaner look
                             pval_cutoff = 0.05,
                             zscore_cutoff = 2) {
  
  # Prepare data
  data <- data %>%
    mutate(
      up_down = case_when(
        is.na(`z-score`) | abs(`z-score`) < 0.01 ~ "Not Significant",
        `z-score` > 0 ~ "Activated",
        `z-score` < 0 ~ "Inhibited"
      ),
      sig = ifelse(`-log(p-value)` > -log(pval_cutoff), "*", ""),
      # Truncate pathway names
      pathway_display = ifelse(
        nchar(`Ingenuity Canonical Pathways`) > 35,
        paste0(substr(`Ingenuity Canonical Pathways`, 1, 32), "..."),
        `Ingenuity Canonical Pathways`
      )
    ) %>%
    arrange(desc(`-log(p-value)`)) %>%
    slice_head(n = n_pathways) %>%
    arrange(`z-score`)  # Order by z-score for vertical plot
  
  # Add row numbers for positioning
  data$position <- seq_len(nrow(data))
  
  # Set factor levels
  data$up_down <- factor(data$up_down, 
                         levels = c("Inhibited", "Not Significant", "Activated"))
  
  # Create the plot
  p <- ggplot(data, aes(x = position, y = `z-score`, fill = up_down)) +
    # Add bars
    geom_col(width = 0.8, alpha = 0.9) +
    
    # Add pathway names on bars
    geom_text(aes(label = pathway_display,
                  y = 0.05 * sign(`z-score`)),
              hjust = ifelse(data$`z-score` > 0, 0, 1),
              size = 3.5, fontface = "bold", color = "white") +
    
    # Add z-score values at the end of bars
    geom_text(aes(label = round(`z-score`, 2),
                  y = `z-score` + sign(`z-score`) * 0.15),
              size = 3, fontface = "bold") +
    
    # Add significance indicators
    geom_text(aes(label = sig,
                  y = `z-score` + sign(`z-score`) * 0.4),
              size = 5, color = "red") +
    
    # Reference lines
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    geom_hline(yintercept = c(-zscore_cutoff, zscore_cutoff), 
               linetype = "dotted", color = "gray40", alpha = 0.7) +
    
    # Styling
    scale_fill_manual(
      values = c("Inhibited" = "#2563EB",      # Bright blue
                 "Not Significant" = "#94A3B8", # Slate gray
                 "Activated" = "#EF4444"),      # Bright red
      name = "Status"
    ) +
    
    scale_y_continuous(
      breaks = pretty_breaks(n = 6),
      limits = c(min(c(-4, min(data$`z-score`, na.rm = TRUE) * 1.2)),
                 max(c(4, max(data$`z-score`, na.rm = TRUE) * 1.2)))
    ) +
    
    scale_x_continuous(
      breaks = NULL,
      expand = c(0.02, 0.02)
    ) +
    
    # Labels
    labs(
      x = NULL,
      y = "Pathway Activation Score (z-score)",
      title = title %||% "IPA Canonical Pathway Analysis",
      subtitle = subtitle %||% paste0("Top ", n_pathways, " significant pathways"),
      caption = paste0("* p < ", pval_cutoff, " | Reference lines at z = ±", zscore_cutoff)
    ) +
    
    # Theme
    theme_minimal(base_size = 12) +
    theme(
      # Title
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, color = "gray40", hjust = 0.5),
      plot.caption = element_text(size = 9, color = "gray50"),
      
      # Axes
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      
      # Legend
      legend.position = c(0.85, 0.92),
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm"),
      
      # Grid and background
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#FAFAFA", color = NA),
      
      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

# Usage with your data
# Read data
ast <- readxl::read_xls(fs::path(dir.dat, "Today_AST_log2_ipa.xls"), skip = 1)

# Version 1: Horizontal bars with names on y-axis
ast_plot_v1 <- ipa_plot(ast, 
                        title = "AST-Associated Canonical Pathways",
                        subtitle = "Pathway activation analysis",
                        n_pathways = 40)

# Version 2: Vertical bars with names on bars (cleaner for fewer pathways)
ast_plot_v2 <- ipa_plot_labeled(ast,
                                title = "AST-Associated Canonical Pathways",
                                subtitle = "Top activated and inhibited pathways",
                                n_pathways = 20)

# Display plots
print(ast_plot_v1)
print(ast_plot_v2)

# Save plots
ggsave("ast_pathways_horizontal.pdf", ast_plot_v1, width = 12, height = 10, dpi = 300)
ggsave("ast_pathways_vertical.pdf", ast_plot_v2, width = 14, height = 8, dpi = 300)

# Combine multiple plots if you have ALT data too
# library(patchwork)
# combined_plot <- ast_plot_v1 / alt_plot_v1 + 
#   plot_annotation(title = "Liver Enzyme Associated Pathways",
#                   theme = theme(plot.title = element_text(size = 16, face = "bold")))