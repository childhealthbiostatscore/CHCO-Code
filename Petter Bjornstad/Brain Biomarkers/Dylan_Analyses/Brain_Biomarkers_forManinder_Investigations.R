library(tidyverse)
library(broom)
library(patchwork)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(knitr)
library(kableExtra)
library(moments)  # For skewness calculation

# Define Quanterix brain biomarkers
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Define clinical predictors of interest
clinical_predictors <- list(
  "Glycemic Control" = c("hba1c", "fbg"),
  "Insulin Sensitivity" = c("avg_m_fsoc", "homa_ir", "adipose_ir", "search_eis"),
  "CGM Metrics" = c("cgm_mean_glucose", "cgm_sd", "cgm_cv", "time_in_range", 
                    "time_above_range", "time_below_range"),
  "Blood Pressure" = c("sbp", "dbp", "map") 
)

# Flatten the list for easier use
all_predictors <- unlist(clinical_predictors, use.names = FALSE)

# ============================================================
# DATA LOADING AND PREPARATION
# ============================================================

cat("\n##########################################################")
cat("\n### LOADING AND PREPARING DATA")
cat("\n##########################################################\n\n")

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

# Summarize data by participant and visit
dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

# ============================================================
# DETERMINE SAMPLE TYPE (SERUM/PLASMA) BY STUDY
# ============================================================

cat("Determining sample type by study/record_id prefix...\n")

dat <- dat %>%
  mutate(
    sample_type = case_when(
      # CROCODILE participants (CRC IDs) - plasma
      grepl("^CRC", record_id, ignore.case = TRUE) ~ "plasma",
      # PENGUIN participants (PEN IDs) - plasma
      grepl("^PEN", record_id, ignore.case = TRUE) ~ "plasma",
      # RH2 participants - serum
      study == "RENAL-HEIRitage" ~ "serum",
      study == "RENAL-HEIR" ~ "serum",
      # Default to study-based assignment
      study == "CROCODILE" ~ "plasma",
      study == "PENGUIN" ~ "plasma",
      TRUE ~ NA_character_
    )
  )



# ============================================================
# FILTER DATA - BASELINE WITH QUANTERIX DATA
# ============================================================

dat_baseline <- dat %>% 
  filter(visit == 'baseline') %>%
  filter(!is.na(ab40_avg_conc)) %>%  # Has Quanterix data
  filter(!is.na(sample_type)) %>%     # Known sample type
  filter(!is.na(age) & !is.na(sex) & !is.na(bmi))  # Has covariates


dup_ids <- dat_baseline$record_id[which(dat_baseline$mrn %in% dat_baseline$mrn[which(duplicated(dat_baseline$mrn))])]


dat_baseline_dup <- dat_baseline %>% filter(record_id %in% dup_ids) %>% 
  dplyr::select(record_id, mrn, group, study, ab40_avg_conc, ab42_avg_conc, adipose_ir, search_eis) %>%
  arrange(by = mrn)



library(ggplot2)
library(tidyr)

# Define your variables
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Reshape data to long format
dat_long <- dat_baseline %>%
  pivot_longer(cols = all_of(qx_var),
               names_to = "biomarker",
               values_to = "concentration")

# Create the plot
p <- ggplot(dat_long, aes(x = group, y = concentration, fill = sample_type)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("plasma" = "#56B4E9", "serum" = "#E69F00")) +
  theme_minimal() +
  labs(x = "Group", 
       y = "Concentration",
       fill = "Sample Type") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10))

# Save as JPEG with high resolution to specified path
ggsave("C:/Users/netio/Downloads/violin_plot_biomarkers.jpeg", 
       plot = p,
       width = 14, 
       height = 10, 
       dpi = 300,
       units = "in")



library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# Define your variables
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Filter for lean controls and reshape data to long format
dat_long <- dat_baseline %>%
  filter(group == "Lean Control") %>%
  pivot_longer(cols = all_of(qx_var),
               names_to = "biomarker",
               values_to = "concentration")

# Create the plot
p <- ggplot(dat_long, aes(x = sample_type, y = concentration, fill = sample_type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  facet_wrap(~ biomarker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("plasma" = "#56B4E9", "serum" = "#E69F00")) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("plasma", "serum")),
                     label = "p.signif",
                     tip.length = 0.02) +
  theme_minimal() +
  labs(x = "Sample Type", 
       y = "Concentration",
       fill = "Sample Type",
       title = "Lean Controls - Biomarker Concentrations by Sample Type") +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10))

# Save as JPEG with high resolution to specified path
ggsave("C:/Users/netio/Downloads/violin_plot_lean_controls_biomarkers.jpeg", 
       plot = p,
       width = 14, 
       height = 8, 
       dpi = 300,
       units = "in")












