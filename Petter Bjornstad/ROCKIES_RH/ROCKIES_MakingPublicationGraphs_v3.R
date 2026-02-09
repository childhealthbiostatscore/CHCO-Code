########################################################################
# ROCKIES Figure 1 - SGLT2 Inhibition Reduces Kidney Oxidative Metabolism
########################################################################

library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(lme4)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

########################################################################
# DATA LOADING - European-formatted Excel file
########################################################################

# Read the file - try xlsx first, fall back to csv
tryCatch({
  raw <- read_excel("C:/Users/netio/Downloads/PET data ROCKIES variables(Blad1).csv")
}, error = function(e) {
  message("Trying as CSV with European formatting...")
})

# If it's actually a CSV with European formatting (semicolon-delimited, comma decimals):
raw <- read.csv2("C:/Users/netio/Downloads/PET data ROCKIES variables(Blad1).csv",
                 stringsAsFactors = FALSE)

# If read.csv2 doesn't work, try this:
# raw <- read.csv("C:/Users/netio/Downloads/PET data ROCKIES variables(Blad1).csv",
#                  sep = ";", dec = ",", stringsAsFactors = FALSE)

# Inspect the data
str(raw)
head(raw)
names(raw)

########################################################################
# DATA CLEANING
# Convert any remaining European decimals and ensure numeric columns
########################################################################

# Helper: convert European decimal strings to numeric
euro_to_numeric <- function(x) {
  if (is.character(x)) {
    x <- gsub(",", ".", x)
    return(as.numeric(x))
  }
  return(as.numeric(x))
}

# Apply to all columns that should be numeric (adjust column names after inspection)
dat <- raw %>%
  mutate(across(where(is.character), ~ {
    # Try converting to numeric; if mostly NA, keep as character
    converted <- suppressWarnings(euro_to_numeric(.x))
    if (sum(!is.na(converted)) > sum(!is.na(.x)) * 0.5) converted else .x
  }))

# Print column names to identify variables
cat("Column names:\n")
print(names(dat))
cat("\nFirst few rows:\n")
print(head(dat))

########################################################################
# VARIABLE IDENTIFICATION
# Adjust these mappings based on actual column names in your data
# Run the code above first, then update these variable names
########################################################################

# Expected variables (update after seeing actual column names):
# - Subject/patient ID
# - Treatment period (ertugliflozin vs placebo)
# - Cortical k2 (two-compartment model)
# - Medullary k2
# - Cortical K_mono (monoexponential model)
# - Medullary K_mono
# - Whole-kidney K_mono
# - Cortical F (perfusion)
# - Insulin sensitivity measures (fasting + postprandial)
# - Tubular sodium reabsorption
# - mGFR / eGFR
# - ERPF

# ====================================================================
# PLACEHOLDER VARIABLE MAPPING - UPDATE THESE AFTER INSPECTING DATA
# ====================================================================
# Uncomment and modify after you see the actual column names:

# id_var        <- "subject_id"        # participant ID
# treatment_var <- "treatment"         # "ertugliflozin" vs "placebo"
# cortical_k2   <- "cortical_k2"      # cortical oxidative rate (2-compartment)
# medullary_k2  <- "medullary_k2"     # medullary oxidative rate
# cortical_kmono <- "cortical_kmono"  # cortical K_mono
# medullary_kmono <- "medullary_kmono"
# whole_kmono   <- "whole_kidney_kmono"
# cortical_F    <- "cortical_F"       # cortical perfusion
# insulin_sens  <- "insulin_sensitivity" # e.g., Matsuda index or M-value
# sodium_reabs  <- "sodium_reabsorption"
# mgfr_var      <- "mGFR"

########################################################################
# THEME FOR PUBLICATION FIGURES
########################################################################

theme_rockies <- theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 10, 5, 10)
  )

# Color palette
col_placebo <- "#4A90D9"     # blue for placebo
col_ertu    <- "#E74C3C"     # red for ertugliflozin
col_line    <- "gray60"      # paired lines

########################################################################
# FIGURE 1G: Cortical k2 - Ertugliflozin vs Placebo (paired)
########################################################################

# Once you identify the correct column names, use this template:
# Reshape data to have one row per subject per treatment

make_paired_plot <- function(data, id, treatment, value, 
                             ylab, title, tag,
                             placebo_label = "Placebo", 
                             drug_label = "Ertugliflozin") {
  
  df <- data %>%
    dplyr::select(id = all_of(id), 
                  trt = all_of(treatment), 
                  val = all_of(value)) %>%
    filter(!is.na(val)) %>%
    mutate(trt = factor(trt, levels = c(placebo_label, drug_label)))
  
  # Summary stats
  summ <- df %>% group_by(trt) %>% 
    summarise(m = mean(val), se = sd(val)/sqrt(n()), .groups = "drop")
  
  # Paired test
  wide <- df %>% pivot_wider(names_from = trt, values_from = val)
  pval <- wilcox.test(wide[[placebo_label]], wide[[drug_label]], paired = TRUE)$p.value
  pval_label <- ifelse(pval < 0.001, "p<0.001", 
                       ifelse(pval < 0.01, paste0("p=", formatC(pval, format="f", digits=3)),
                              paste0("p=", formatC(pval, format="f", digits=3))))
  
  ggplot(df, aes(x = trt, y = val)) +
    # Paired lines
    geom_line(aes(group = id), color = col_line, alpha = 0.4, linewidth = 0.4) +
    # Individual points
    geom_point(aes(color = trt), size = 2.5, alpha = 0.7) +
    # Mean ± SEM bars
    geom_errorbar(data = summ, aes(x = trt, y = m, ymin = m - se, ymax = m + se),
                  width = 0.15, linewidth = 0.8, inherit.aes = FALSE) +
    geom_point(data = summ, aes(x = trt, y = m), size = 4, shape = 18,
               inherit.aes = FALSE) +
    # P-value annotation
    annotate("text", x = 1.5, y = max(df$val) * 1.05, 
             label = pval_label, size = 3.5, fontface = "italic") +
    annotate("segment", x = 1, xend = 2, 
             y = max(df$val) * 1.03, yend = max(df$val) * 1.03,
             linewidth = 0.4) +
    scale_color_manual(values = c(col_placebo, col_ertu)) +
    labs(x = NULL, y = ylab, title = title, tag = tag) +
    theme_rockies
}

########################################################################
# FIGURE 1C-E: Correlation plots (placebo conditions)
########################################################################

make_corr_plot <- function(data, xvar, yvar, xlab, ylab, title, tag,
                           method = "spearman") {
  
  df <- data %>%
    dplyr::select(x = all_of(xvar), y = all_of(yvar)) %>%
    filter(!is.na(x) & !is.na(y))
  
  cor_test <- cor.test(df$x, df$y, method = method)
  rho <- cor_test$estimate
  pval <- cor_test$p.value
  
  label <- paste0("rho = ", formatC(rho, format = "f", digits = 2),
                  "\np = ", ifelse(pval < 0.001, "<0.001", 
                                   formatC(pval, format = "f", digits = 3)))
  
  ggplot(df, aes(x = x, y = y)) +
    geom_point(size = 2.5, alpha = 0.7, color = col_placebo) +
    geom_smooth(method = "lm", se = TRUE, color = "black", 
                linewidth = 0.8, fill = "gray85") +
    annotate("text", x = min(df$x) + diff(range(df$x)) * 0.05,
             y = max(df$y) - diff(range(df$y)) * 0.05,
             label = label, hjust = 0, size = 3.2, fontface = "italic") +
    labs(x = xlab, y = ylab, title = title, tag = tag) +
    theme_rockies
}

########################################################################
# FIGURE 1I-J: Delta correlations (change with treatment)
########################################################################

make_delta_corr_plot <- function(data, id, treatment, xvar, yvar, 
                                 xlab, ylab, title, tag,
                                 placebo_label = "Placebo",
                                 drug_label = "Ertugliflozin") {
  
  df <- data %>%
    dplyr::select(id = all_of(id), trt = all_of(treatment),
                  x = all_of(xvar), y = all_of(yvar)) %>%
    filter(!is.na(x) & !is.na(y))
  
  # Calculate deltas (drug - placebo)
  wide_x <- df %>% dplyr::select(id, trt, x) %>%
    pivot_wider(names_from = trt, values_from = x)
  wide_y <- df %>% dplyr::select(id, trt, y) %>%
    pivot_wider(names_from = trt, values_from = y)
  
  deltas <- tibble(
    id = wide_x$id,
    dx = wide_x[[drug_label]] - wide_x[[placebo_label]],
    dy = wide_y[[drug_label]] - wide_y[[placebo_label]]
  ) %>% filter(!is.na(dx) & !is.na(dy))
  
  cor_test <- cor.test(deltas$dx, deltas$dy, method = "spearman")
  rho <- cor_test$estimate
  pval <- cor_test$p.value
  
  label <- paste0("rho = ", formatC(rho, format = "f", digits = 2),
                  "\np = ", ifelse(pval < 0.001, "<0.001",
                                   formatC(pval, format = "f", digits = 3)))
  
  ggplot(deltas, aes(x = dx, y = dy)) +
    geom_point(size = 2.5, alpha = 0.7, color = "#8E44AD") +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linewidth = 0.8, fill = "gray85") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    annotate("text", x = min(deltas$dx) + diff(range(deltas$dx)) * 0.05,
             y = max(deltas$dy) - diff(range(deltas$dy)) * 0.05,
             label = label, hjust = 0, size = 3.2, fontface = "italic") +
    labs(x = paste0("\u0394 ", xlab), y = paste0("\u0394 ", ylab),
         title = title, tag = tag) +
    theme_rockies
}

########################################################################
# ASSEMBLE FIGURE 1
# ---------------------------------------------------------------
# After inspecting column names, fill in the variable mappings and
# uncomment the plotting calls below.
########################################################################

# Example assembly (update variable names!):

# --- Panels C-D: K_mono vs insulin resistance (placebo only) ---
# dat_placebo <- dat %>% filter(!!sym(treatment_var) == "Placebo")
#
# fig1c <- make_corr_plot(dat_placebo, insulin_sens, cortical_kmono,
#                          "Insulin Resistance", 
#                          expression(bold("Cortical K"["mono"]*" (min"^{-1}*")")),
#                          "Cortex", "C")
#
# fig1d <- make_corr_plot(dat_placebo, insulin_sens, whole_kmono,
#                          "Insulin Resistance",
#                          expression(bold("Whole-Kidney K"["mono"]*" (min"^{-1}*")")),
#                          "Whole Kidney", "D")

# --- Panel E: Sodium reabsorption vs medullary K_mono (placebo) ---
# fig1e <- make_corr_plot(dat_placebo, sodium_reabs, medullary_kmono,
#                          "Tubular Na Reabsorption",
#                          expression(bold("Medullary K"["mono"]*" (min"^{-1}*")")),
#                          "Medulla", "E")

# --- Panel G: Cortical k2 paired plot ---
# fig1g <- make_paired_plot(dat, id_var, treatment_var, cortical_k2,
#                            expression(bold("Cortical k"[2]*" (min"^{-1}*")")),
#                            "Cortical Oxidative Rate", "G")

# --- Panel H: Medullary k2 paired plot ---
# fig1h <- make_paired_plot(dat, id_var, treatment_var, medullary_k2,
#                            expression(bold("Medullary k"[2]*" (min"^{-1}*")")),
#                            "Medullary Oxidative Rate", "H")

# --- Panel I: Delta insulin sensitivity vs delta cortical k2 ---
# fig1i <- make_delta_corr_plot(dat, id_var, treatment_var,
#                                insulin_sens, cortical_k2,
#                                "Insulin Sensitivity",
#                                expression("Cortical k"[2]*" (min"^{-1}*")"),
#                                "Cortex", "I")

# --- Panel J: Delta insulin sensitivity vs delta medullary k2 ---
# fig1j <- make_delta_corr_plot(dat, id_var, treatment_var,
#                                insulin_sens, medullary_k2,
#                                "Insulin Sensitivity",
#                                expression("Medullary k"[2]*" (min"^{-1}*")"),
#                                "Medulla", "J")

# ====================================================================
# COMPOSITE FIGURE
# ====================================================================
# Panels A (study design) and B (PET schematic) are typically made in
# PowerPoint/Illustrator, then combined with data panels.
# Panel F (metabolomics) would come from separate analysis.
#
# For the data-driven panels:

# fig1_data <- (fig1c | fig1d | fig1e) /
#              (fig1g | fig1h | plot_spacer()) /
#              (fig1i | fig1j | plot_spacer()) +
#   plot_annotation(
#     title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
#     theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
#   )
#
# ggsave(paste0(base_path, "Figure1_ROCKIES_data_panels.pdf"),
#        fig1_data, width = 14, height = 12, dpi = 300)
# ggsave(paste0(base_path, "Figure1_ROCKIES_data_panels.png"),
#        fig1_data, width = 14, height = 12, dpi = 300)

cat("\n============================================================\n")
cat("NEXT STEPS:\n")
cat("1. Run the data loading section above\n")
cat("2. Check the printed column names\n")
cat("3. Update the variable mappings in the PLACEHOLDER section\n")
cat("4. Uncomment and run the plotting calls\n")
cat("============================================================\n")