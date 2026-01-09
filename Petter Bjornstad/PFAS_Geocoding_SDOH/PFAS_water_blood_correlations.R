###PFAS Correlations



library(dplyr)
library(stringr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(tidygeocoder)



imputed_pfas <- readRDS("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/PANTHER/Data_Cleaned/PFAS_Data_Imputed_11_03.rds")

base.dir <- 'C:/Users/netio/Documents/UofW/Projects/PFAS_Water/'

water_data <- data.table::fread(paste0(base.dir, 'participants_with_study_classification.csv'))

harmonized_path <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"
harmonized_data <- read.csv(harmonized_path, na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))




####Analysis 

### PFAS Correlations
library(dplyr)
library(stringr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(corrplot)

# Load data
imputed_pfas <- readRDS("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/PANTHER/Data_Cleaned/PFAS_Data_Imputed_11_03.rds")
base.dir <- 'C:/Users/netio/Documents/UofW/Projects/PFAS_Water/'
water_data <- data.table::fread(paste0(base.dir, 'participants_with_study_classification.csv'))
harmonized_path <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"
harmonized_data <- read.csv(harmonized_path, na = '')

# Summarize harmonized data
dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

# Convert to plain data.frame to avoid any tibble/data.table issues
dat <- as.data.frame(dat)
imputed_pfas <- as.data.frame(imputed_pfas)
water_data <- as.data.frame(water_data)

# Define column names explicitly
serum_pfas_cols <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA", "PFOA", 
                     "PFBS", "PFHps", "PFHxS", "PFOS", "PFPeAS")

water_pfas_cols <- c(
  "hexafluoropropylene_oxide_dimer_acid",
  "perfluorooctanesulfonic_acid_pfos_",
  "perfluoroundecanoic_acid_pfuna_",
  "n_methyl_perfluorooctanesulfonamidoacetic_acid_mefosaa_",
  "n_ethyl_perfluorooctanesulfonamidoacetic_acid_etfosaa_",
  "perfluorohexanoic_acid_pfhxa_",
  "perfluorododecanoic_acid_pfdoa_",
  "perfluorooctanoic_acid_pfoa_",
  "perfluorodecanoic_acid_pfda_",
  "perfluorohexanesulfonic_acid_pfhxs_",
  "perfluorobutanesulfonic_acid_pfbs_",
  "perfluoroheptanoic_acid_pfhpa_",
  "perfluorononanoic_acid_pfna_",
  "perfluorotetradecanoic_pftea_",
  "perfluorotridecanoic_acid_pftria_",
  "9_chlorohexadecafluoro_3_oxanone_1_sulfonic_acid_9cl_pf3ons_",
  "11_chloroeicosafluoro_3_oxaundecane_1_sulfonic_acid_11cl_pf3ouds_",
  "4_8_dioxa_3h_perfluorononanoic_acid_adona_",
  "pfos_and_pfoa_total_combined",
  "total_pfas"
)

outcome_vars <- c("eGFR_CKD_epi", "acr_u")

# --- MERGE USING BASE R merge() ---
# This is more reliable than dplyr joins for mixed data sources
cat("=== MERGING DATA ===\n")

# Step 1: Merge dat with imputed_pfas
merged_data <- merge(dat, imputed_pfas, by = "record_id", all.x = TRUE)
cat("After merging with serum PFAS:\n")
cat("  Rows:", nrow(merged_data), "\n")
cat("  PFOS non-NA:", sum(!is.na(merged_data$PFOS)), "\n")

# Step 2: Merge with water_data
merged_data <- merge(merged_data, water_data, by = "record_id", all.x = TRUE)
cat("After merging with water PFAS:\n")
cat("  Rows:", nrow(merged_data), "\n")
cat("  total_pfas non-NA:", sum(!is.na(merged_data$total_pfas)), "\n\n")

# Create matching pairs (serum PFAS to water PFAS)
pfas_pairs <- tribble(
  ~serum,       ~water,
  "PFOS",       "perfluorooctanesulfonic_acid_pfos_",
  "PFOA",       "perfluorooctanoic_acid_pfoa_",
  "PFDA",       "perfluorodecanoic_acid_pfda_",
  "PFHxS",      "perfluorohexanesulfonic_acid_pfhxs_",
  "PFBS",       "perfluorobutanesulfonic_acid_pfbs_",
  "PFHpA",      "perfluoroheptanoic_acid_pfhpa_",
  "PFNA",       "perfluorononanoic_acid_pfna_",
  "N.MeFOSAA",  "n_methyl_perfluorooctanesulfonamidoacetic_acid_mefosaa_"
)

# --- Minimum sample size threshold ---
MIN_N <- 10

# --- Correlation function with error handling ---
safe_cor_test <- function(x, y, method = "spearman", min_n = MIN_N) {
  complete_idx <- complete.cases(x, y)
  n <- sum(complete_idx)
  
  if (n < min_n) {
    return(tibble(rho = NA_real_, p_value = NA_real_, n = n, 
                  skipped = TRUE, skip_reason = paste0("n < ", min_n)))
  }
  
  x_complete <- x[complete_idx]
  y_complete <- y[complete_idx]
  
  if (sd(x_complete, na.rm = TRUE) == 0 || sd(y_complete, na.rm = TRUE) == 0) {
    return(tibble(rho = NA_real_, p_value = NA_real_, n = n,
                  skipped = TRUE, skip_reason = "Zero variance"))
  }
  
  tryCatch({
    test <- cor.test(x_complete, y_complete, method = method)
    tibble(rho = test$estimate, p_value = test$p.value, n = n,
           skipped = FALSE, skip_reason = NA_character_)
  }, error = function(e) {
    tibble(rho = NA_real_, p_value = NA_real_, n = n,
           skipped = TRUE, skip_reason = as.character(e$message))
  })
}

# --- Function to summarize data availability ---
summarize_data_availability <- function(data, vars) {
  map_dfr(vars, function(v) {
    if (!v %in% names(data)) {
      return(tibble(variable = v, n_total = NA_integer_, n_nonmissing = 0L, 
                    pct_available = 0, status = "Column not found"))
    }
    x <- data[[v]]
    n_total <- length(x)
    n_nonmissing <- sum(!is.na(x))
    pct <- round(100 * n_nonmissing / n_total, 1)
    status <- case_when(
      n_nonmissing < MIN_N ~ "Insufficient data",
      pct < 10 ~ "Very sparse (<10%)",
      pct < 50 ~ "Sparse (<50%)",
      TRUE ~ "OK"
    )
    tibble(variable = v, n_total = n_total, n_nonmissing = n_nonmissing,
           pct_available = pct, status = status)
  })
}

# --- Check data availability ---
cat("=== DATA AVAILABILITY SUMMARY ===\n\n")

cat("Serum PFAS:\n")
serum_availability <- summarize_data_availability(merged_data, serum_pfas_cols)
print(serum_availability, n = Inf)

cat("\nDrinking Water PFAS:\n")
water_availability <- summarize_data_availability(merged_data, water_pfas_cols)
print(water_availability, n = Inf)

cat("\nKidney Outcomes:\n")
outcome_availability <- summarize_data_availability(merged_data, outcome_vars)
print(outcome_availability, n = Inf)

# Save availability summary
all_availability <- bind_rows(
  serum_availability %>% mutate(source = "Serum PFAS"),
  water_availability %>% mutate(source = "Water PFAS"),
  outcome_availability %>% mutate(source = "Outcome")
)
write.csv(all_availability, paste0(base.dir, "data_availability_summary.csv"), row.names = FALSE)

# Filter to usable variables
usable_serum <- serum_availability %>% filter(n_nonmissing >= MIN_N) %>% pull(variable)
usable_water <- water_availability %>% filter(n_nonmissing >= MIN_N) %>% pull(variable)
usable_outcomes <- outcome_availability %>% filter(n_nonmissing >= MIN_N) %>% pull(variable)

cat("\n=== USABLE VARIABLES ===\n")
cat("Serum PFAS (n >=", MIN_N, "):", length(usable_serum), "of", length(serum_pfas_cols), "\n")
cat("Water PFAS (n >=", MIN_N, "):", length(usable_water), "of", length(water_pfas_cols), "\n")
cat("Outcomes (n >=", MIN_N, "):", length(usable_outcomes), "of", length(outcome_vars), "\n\n")

# Initialize empty result tibbles
serum_water_cors <- tibble()
serum_outcome_cors <- tibble()
water_outcome_cors <- tibble()
skipped_sw <- tibble()
skipped_so <- tibble()
skipped_wo <- tibble()

# --- 1. Serum PFAS ~ Drinking Water PFAS correlations ---
usable_pairs <- pfas_pairs %>%
  filter(serum %in% usable_serum, water %in% usable_water)

cat("Running serum vs water correlations for", nrow(usable_pairs), "pairs...\n")

if (nrow(usable_pairs) > 0) {
  serum_water_cors <- usable_pairs %>%
    rowwise() %>%
    mutate(cor_result = list(safe_cor_test(merged_data[[serum]], merged_data[[water]]))) %>%
    unnest(cor_result) %>%
    ungroup() %>%
    mutate(comparison = "Serum vs Water PFAS")
  
  skipped_sw <- serum_water_cors %>% filter(skipped)
  if (nrow(skipped_sw) > 0) {
    cat("\nSkipped serum-water pairs:\n")
    print(skipped_sw %>% dplyr::select(serum, water, n, skip_reason))
  }
  
  cat("\nSerum PFAS vs Drinking Water PFAS:\n")
  print(serum_water_cors %>% filter(!skipped) %>% dplyr::select(-skipped, -skip_reason))
} else {
  cat("No serum-water pairs with sufficient data\n")
}

# --- 2. Serum PFAS ~ eGFR and Albuminuria ---
cat("\nRunning serum PFAS vs outcome correlations...\n")

if (length(usable_serum) > 0 && length(usable_outcomes) > 0) {
  serum_outcome_cors <- expand_grid(pfas = usable_serum, outcome = usable_outcomes) %>%
    rowwise() %>%
    mutate(cor_result = list(safe_cor_test(merged_data[[pfas]], merged_data[[outcome]]))) %>%
    unnest(cor_result) %>%
    ungroup() %>%
    mutate(p_adjusted = p.adjust(p_value, method = "BH"), significant = p_adjusted < 0.05) %>%
    arrange(p_value)
  
  skipped_so <- serum_outcome_cors %>% filter(skipped)
  if (nrow(skipped_so) > 0) {
    cat("\nSkipped serum-outcome pairs:\n")
    print(skipped_so %>% dplyr::select(pfas, outcome, n, skip_reason))
  }
  
  cat("\nSerum PFAS vs Kidney Outcomes:\n")
  print(serum_outcome_cors %>% filter(!skipped) %>% dplyr::select(-skipped, -skip_reason))
} else {
  cat("Insufficient data for serum-outcome correlations\n")
}

# --- 3. Drinking Water PFAS ~ eGFR and Albuminuria ---
cat("\nRunning water PFAS vs outcome correlations...\n")

if (length(usable_water) > 0 && length(usable_outcomes) > 0) {
  water_outcome_cors <- expand_grid(pfas = usable_water, outcome = usable_outcomes) %>%
    rowwise() %>%
    mutate(cor_result = list(safe_cor_test(merged_data[[pfas]], merged_data[[outcome]]))) %>%
    unnest(cor_result) %>%
    ungroup() %>%
    mutate(p_adjusted = p.adjust(p_value, method = "BH"), significant = p_adjusted < 0.05) %>%
    arrange(p_value)
  
  skipped_wo <- water_outcome_cors %>% filter(skipped)
  if (nrow(skipped_wo) > 0) {
    cat("\nSkipped water-outcome pairs:\n")
    print(skipped_wo %>% dplyr::select(pfas, outcome, n, skip_reason))
  }
  
  cat("\nDrinking Water PFAS vs Kidney Outcomes:\n")
  print(water_outcome_cors %>% filter(!skipped) %>% dplyr::select(-skipped, -skip_reason))
} else {
  cat("Insufficient data for water-outcome correlations\n")
}

# --- Combine all results ---
all_correlations <- bind_rows(
  if (nrow(serum_water_cors) > 0) serum_water_cors %>% mutate(analysis = "Serum-Water") else NULL,
  if (nrow(serum_outcome_cors) > 0) serum_outcome_cors %>% rename(serum = pfas) %>% mutate(analysis = "Serum-Outcome") else NULL,
  if (nrow(water_outcome_cors) > 0) water_outcome_cors %>% rename(water = pfas) %>% mutate(analysis = "Water-Outcome") else NULL
)

# --- Save results ---
write.csv(all_correlations, paste0(base.dir, "pfas_correlation_results.csv"), row.names = FALSE)

# --- PLOTS ---
plot_dir <- paste0(base.dir, "plots/")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Log-transform UACR for visualization
merged_data$log_acr_u <- log10(merged_data$acr_u + 1)

cat("\n=== GENERATING PLOTS ===\n")
cat("Saving to:", plot_dir, "\n\n")

# --- Plot 1: Serum PFAS vs Water PFAS ---
if (nrow(serum_water_cors %>% filter(!skipped)) > 0) {
  plottable_pairs <- serum_water_cors %>% filter(!skipped)
  
  pdf(paste0(plot_dir, "serum_vs_water_pfas.pdf"), width = 12, height = 10)
  plots <- plottable_pairs %>%
    rowwise() %>%
    mutate(plot = list({
      df <- data.frame(x = merged_data[[water]], y = merged_data[[serum]]) %>% filter(!is.na(x) & !is.na(y))
      ggplot(df, aes(x = x, y = y)) +
        geom_point(alpha = 0.5, color = "steelblue") +
        geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
        labs(x = paste("Water:", serum), y = paste("Serum:", serum),
             title = paste0(serum, ": Serum vs Water"),
             subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", rho, p_value, n)) +
        theme_bw() + theme(plot.title = element_text(face = "bold"))
    }))
  gridExtra::grid.arrange(grobs = plots$plot, ncol = 3)
  dev.off()
  cat("  - serum_vs_water_pfas.pdf\n")
} else {
  cat("  - serum_vs_water_pfas.pdf SKIPPED\n")
}

# --- Plot 2: Serum PFAS vs eGFR ---
if (nrow(serum_outcome_cors %>% filter(!skipped, outcome == "eGFR_CKD_epi")) > 0) {
  plot_data <- serum_outcome_cors %>% filter(!skipped, outcome == "eGFR_CKD_epi")
  
  pdf(paste0(plot_dir, "serum_pfas_vs_eGFR.pdf"), width = 12, height = 10)
  plots <- map(plot_data$pfas, function(p) {
    row <- plot_data %>% filter(pfas == p)
    df <- data.frame(x = merged_data[[p]], y = merged_data$eGFR_CKD_epi) %>% filter(!is.na(x) & !is.na(y))
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
      labs(x = p, y = "eGFR (CKD-EPI)", title = paste0(p, " vs eGFR"),
           subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", row$rho, row$p_value, row$n)) +
      theme_bw() + theme(plot.title = element_text(face = "bold"))
  })
  gridExtra::grid.arrange(grobs = plots, ncol = 3)
  dev.off()
  cat("  - serum_pfas_vs_eGFR.pdf\n")
} else {
  cat("  - serum_pfas_vs_eGFR.pdf SKIPPED\n")
}

# --- Plot 3: Serum PFAS vs UACR ---
if (nrow(serum_outcome_cors %>% filter(!skipped, outcome == "acr_u")) > 0) {
  plot_data <- serum_outcome_cors %>% filter(!skipped, outcome == "acr_u")
  
  pdf(paste0(plot_dir, "serum_pfas_vs_UACR.pdf"), width = 12, height = 10)
  plots <- map(plot_data$pfas, function(p) {
    row <- plot_data %>% filter(pfas == p)
    df <- data.frame(x = merged_data[[p]], y = merged_data$log_acr_u) %>% filter(!is.na(x) & !is.na(y) & is.finite(y))
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "darkorange") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
      labs(x = p, y = "log10(UACR + 1)", title = paste0(p, " vs UACR"),
           subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", row$rho, row$p_value, row$n)) +
      theme_bw() + theme(plot.title = element_text(face = "bold"))
  })
  gridExtra::grid.arrange(grobs = plots, ncol = 3)
  dev.off()
  cat("  - serum_pfas_vs_UACR.pdf\n")
} else {
  cat("  - serum_pfas_vs_UACR.pdf SKIPPED\n")
}

# --- Plot 4: Water PFAS vs eGFR ---
water_main <- c("total_pfas", "pfos_and_pfoa_total_combined", "perfluorooctanesulfonic_acid_pfos_", 
                "perfluorooctanoic_acid_pfoa_", "perfluorohexanesulfonic_acid_pfhxs_", "perfluorononanoic_acid_pfna_")
if (nrow(water_outcome_cors %>% filter(!skipped, outcome == "eGFR_CKD_epi", pfas %in% water_main)) > 0) {
  plot_data <- water_outcome_cors %>% filter(!skipped, outcome == "eGFR_CKD_epi", pfas %in% water_main)
  
  pdf(paste0(plot_dir, "water_pfas_vs_eGFR.pdf"), width = 12, height = 8)
  plots <- map(plot_data$pfas, function(p) {
    row <- plot_data %>% filter(pfas == p)
    df <- data.frame(x = merged_data[[p]], y = merged_data$eGFR_CKD_epi) %>% filter(!is.na(x) & !is.na(y))
    short <- str_extract(p, "(?<=_)[a-z]+(?=_$)|total_pfas|pfos_and_pfoa") %>% toupper() %>% coalesce(p)
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "forestgreen") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
      labs(x = paste("Water:", short), y = "eGFR (CKD-EPI)", title = paste0("Water ", short, " vs eGFR"),
           subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", row$rho, row$p_value, row$n)) +
      theme_bw() + theme(plot.title = element_text(face = "bold"))
  })
  gridExtra::grid.arrange(grobs = plots, ncol = 3)
  dev.off()
  cat("  - water_pfas_vs_eGFR.pdf\n")
} else {
  cat("  - water_pfas_vs_eGFR.pdf SKIPPED\n")
}

# --- Plot 5: Water PFAS vs UACR ---
if (nrow(water_outcome_cors %>% filter(!skipped, outcome == "acr_u", pfas %in% water_main)) > 0) {
  plot_data <- water_outcome_cors %>% filter(!skipped, outcome == "acr_u", pfas %in% water_main)
  
  pdf(paste0(plot_dir, "water_pfas_vs_UACR.pdf"), width = 12, height = 8)
  plots <- map(plot_data$pfas, function(p) {
    row <- plot_data %>% filter(pfas == p)
    df <- data.frame(x = merged_data[[p]], y = merged_data$log_acr_u) %>% filter(!is.na(x) & !is.na(y) & is.finite(y))
    short <- str_extract(p, "(?<=_)[a-z]+(?=_$)|total_pfas|pfos_and_pfoa") %>% toupper() %>% coalesce(p)
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "purple") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
      labs(x = paste("Water:", short), y = "log10(UACR + 1)", title = paste0("Water ", short, " vs UACR"),
           subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", row$rho, row$p_value, row$n)) +
      theme_bw() + theme(plot.title = element_text(face = "bold"))
  })
  gridExtra::grid.arrange(grobs = plots, ncol = 3)
  dev.off()
  cat("  - water_pfas_vs_UACR.pdf\n")
} else {
  cat("  - water_pfas_vs_UACR.pdf SKIPPED\n")
}

# --- Plot 6: Correlation Heatmap ---
if (nrow(serum_outcome_cors %>% filter(!skipped)) > 0) {
  pdf(paste0(plot_dir, "correlation_heatmap.pdf"), width = 10, height = 8)
  
  heatmap_data <- serum_outcome_cors %>%
    filter(!skipped) %>%
    mutate(outcome = case_when(outcome == "eGFR_CKD_epi" ~ "eGFR", outcome == "acr_u" ~ "UACR", TRUE ~ outcome)) %>%
    dplyr::select(pfas, outcome, rho) %>%
    pivot_wider(names_from = outcome, values_from = rho) %>%
    column_to_rownames("pfas") %>%
    as.matrix()
  
  corrplot(heatmap_data, method = "color", is.corr = FALSE,
           col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           tl.col = "black", tl.cex = 1, cl.cex = 0.8,
           addCoef.col = "black", number.cex = 0.8,
           title = "Serum PFAS Correlations with Kidney Outcomes\n(Spearman's rho)",
           mar = c(0, 0, 3, 0))
  dev.off()
  cat("  - correlation_heatmap.pdf\n")
} else {
  cat("  - correlation_heatmap.pdf SKIPPED\n")
}

cat("\n=== COMPLETE ===\n")
cat("Results saved to:", base.dir, "\n")
cat("Plots saved to:", plot_dir, "\n")























