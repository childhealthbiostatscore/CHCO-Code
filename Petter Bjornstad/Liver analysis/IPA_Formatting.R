#Libraries
library(ggrepel)
# library(qpcR)
library(ggpubr)
library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)
library(tidyverse)
library(ggtext)
library(scales)

#Directories
user <- "hhampson" 
dir.dat <- c(paste0("/Users/",user,"/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Liver project/IPA Results"))
dir.results <- c(paste0("/Users/",user,"/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Liver project/NEBULA Results"))

color_table <- tibble(
  up_down = c(-1, 0, 1),
  color = c("#990000", "gray", "#003366")
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
    scale_fill_manual(values = c("-1" = "#003366", "0" = "grey", "1" = "indianred")) 
  p <- p + coord_flip()
}

# Enhanced IPA plotting function
ipa_plot_enhanced <- function(data, 
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
    # geom_hline(yintercept = c(-3, 3), 
    #            linetype = "dashed", color = "gray40", alpha = 0.5) +
    # geom_hline(yintercept = c(-2, 2), 
    #            linetype = "dashed", color = "gray40", alpha = 0.5) +
    # geom_hline(yintercept = c(-1, 1), 
    #            linetype = "dashed", color = "gray40", alpha = 0.5) +
    
    # Styling
    scale_fill_manual(
      values = c("Inhibited" = "#003366",      # Nice blue
                 "Not Significant" = "#9CA3AF", # Gray
                 "Activated" = "#990000"),      # Nice red
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

#TODAY
##AST
ast <- readxl::read_xls(fs::path(dir.dat,"Today_AST_log2_ipa.xls"), skip = 1)

ast <- ast %>% arrange(desc("-log(p-value)")) 
ast_keep <- ast[1:30,]
ast_keep <- ast_keep %>%
  mutate(`z-score`=ifelse(is.na(`z-score`),0,`z-score`))

ast_plot <- ipa_plot(ast_keep)
ast_plot 
# png(fs::path(dir.results,"AST_pathways_TODAY.png"),width=3000,height=2000,res=300)
# plot(ast_plot)
# dev.off()

ast2 <- ipa_plot_enhanced(ast_keep, 
                        title = "AST-Associated Canonical Pathways",
                        subtitle = "Pathway activation analysis",
                        n_pathways = 30)
png(fs::path(dir.results,"AST_pathways_TODAY.png"),width=3500,height=3000,res=300)
plot(ast2)
dev.off()

##ALT
alt <- readxl::read_xls(fs::path(dir.dat,"Today_ALT_log2_ipa.xls"), skip = 1)

alt <- alt %>% arrange(desc("-log(p-value)")) 
alt_keep <- alt[1:30,]
alt_keep <- alt_keep %>%
  mutate(`z-score`=ifelse(is.na(`z-score`),0,`z-score`))

alt_plot <- ipa_plot(alt_keep)
alt_plot 
# png(fs::path(dir.results,"ALT_pathways_TODAY.png"),width=3000,height=2000,res=300)
# plot(alt_plot)
# dev.off()

alt2 <- ipa_plot_enhanced(alt_keep, 
                 title = "ALT-Associated Canonical Pathways",
                 subtitle = "Pathway activation analysis",
                 n_pathways = 30)
png(fs::path(dir.results,"ALT_pathways_TODAY.png"),width=3500,height=3000,res=300)
plot(alt2)
dev.off()

##Ratio
ratio <- readxl::read_xls(fs::path(dir.dat,"Today_AST_ALT_ratio_log2_ipa.xls"), skip = 1)

ratio <- ratio %>% arrange(desc("-log(p-value)")) 
ratio_keep <- ratio[1:30,]

ratio_plot <- ipa_plot(ratio_keep)
ratio_plot 
# png(fs::path(dir.results,"AST_ALT_Ratio_pathways_TeenLabs.png"),width=3000,height=2000,res=300)
# plot(ratio_plot)
# dev.off()

ratio2 <- ipa_plot_enhanced(ratio, 
                 title = "AST/ALT-Associated Canonical Pathways",
                 subtitle = "Pathway activation analysis",
                 n_pathways = 30)
png(fs::path(dir.results,"AST_ALT_ratio_pathways_TODAY.png"),width=3500,height=3000,res=300)
plot(ratio2)
dev.off()

#TeenLabs
##AST
ast <- readxl::read_xls(fs::path(dir.dat,"TeenLabs_AST_log2_ipa.xls"), skip = 1)

ast <- ast %>% arrange(desc("-log(p-value)")) 
# ast <- ast %>% #Filter out NA zscore pathways
ast_keep <- ast[1:30,]
ast_keep <- ast_keep %>% 
  mutate(`z-score`=ifelse(is.na(`z-score`),0,`z-score`))
ast_plot <- ipa_plot(ast_keep)
ast_plot 
# png(fs::path(dir.results,"AST_pathways_TeenLabs.png"),width=3000,height=2000,res=300)
# plot(ast_plot)
# dev.off()

ast2 <- ipa_plot_enhanced(ast_keep, 
                 title = "AST-Associated Canonical Pathways",
                 subtitle = "Pathway activation analysis",
                 n_pathways = 30)
png(fs::path(dir.results,"AST_pathways_TeenLabs.png"),width=3500,height=3000,res=300)
plot(ast2)
dev.off()

##ALT
alt <- readxl::read_xls(fs::path(dir.dat,"TeenLabs_ALT_log2_ipa.xls"), skip = 1)

alt <- alt %>% arrange(desc("-log(p-value)")) 
# alt <- alt %>% #Filter out NA zscore pathways
alt_keep <- alt[1:30,]
alt_keep <- alt_keep %>%
  mutate(`z-score`=ifelse(is.na(`z-score`),0,`z-score`))
alt_plot <- ipa_plot(alt_keep)
alt_plot 
# png(fs::path(dir.results,"alt_pathways_TeenLabs.png"),width=3000,height=2000,res=300)
# plot(alt_plot)
# dev.off()

alt2 <- ipa_plot_enhanced(alt_keep, 
                          title = "ALT-Associated Canonical Pathways",
                          subtitle = "Pathway activation analysis",
                          n_pathways = 30)
png(fs::path(dir.results,"ALT_pathways_TeenLabs.png"),width=3500,height=3000,res=300)
plot(alt2)
dev.off()
