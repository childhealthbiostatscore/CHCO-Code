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

# Serum PFAS (from imputed_pfas)
serum_pfas_cols <- c("N.MeFOSAA", "PFDA", "PFHpA", "PFNA", "PFOA", 
                     "PFBS", "PFHps", "PFHxS", "PFOS", "PFPeAS")

# Drinking water PFAS (from water_data)
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

# Kidney outcomes: eGFR_CKD_epi and acr_u

# Merge all datasets
merged_data <- dat %>%
  left_join(imputed_pfas, by = "record_id") %>%
  left_join(water_data, by = "record_id")

# Kidney outcome variables
outcome_vars <- c("eGFR_CKD_epi", "acr_u")

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
MIN_N <- 10  # Minimum observations required for correlation

# --- Correlation function with error handling ---
safe_cor_test <- function(x, y, method = "spearman", min_n = MIN_N) {
  complete_idx <- complete.cases(x, y)
  n <- sum(complete_idx)
  
  # Check for sufficient data
  if (n < min_n) {
    return(tibble(rho = NA_real_, p_value = NA_real_, n = n, 
                  skipped = TRUE, skip_reason = paste0("n < ", min_n)))
  }
  
  # Check for zero variance
  x_complete <- x[complete_idx]
  y_complete <- y[complete_idx]
  
  if (sd(x_complete, na.rm = TRUE) == 0 || sd(y_complete, na.rm = TRUE) == 0) {
    return(tibble(rho = NA_real_, p_value = NA_real_, n = n,
                  skipped = TRUE, skip_reason = "Zero variance"))
  }
  
  # Run correlation
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
      return(tibble(variable = v, n_total = NA, n_nonmissing = 0, 
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

# --- Check data availability before running correlations ---
cat("=== DATA AVAILABILITY SUMMARY ===\n\n")

cat("Serum PFAS:\n")
serum_availability <- summarize_data_availability(merged_data, serum_pfas_cols)
print(serum_availability, n = Inf)

cat("\nDrinking Water PFAS:\n
")
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

# --- 1. Serum PFAS ~ Drinking Water PFAS correlations ---
# Filter pfas_pairs to usable variables
usable_pairs <- pfas_pairs %>%
  filter(serum %in% usable_serum, water %in% usable_water)

cat("Running serum vs water correlations for", nrow(usable_pairs), "pairs...\n")

serum_water_cors <- usable_pairs %>%
  rowwise() %>%
  mutate(
    cor_result = list(safe_cor_test(
      merged_data[[serum]], 
      merged_data[[water]]
    ))
  ) %>%
  unnest(cor_result) %>%
  ungroup() %>%
  mutate(comparison = "Serum vs Water PFAS")

# Report skipped
skipped_sw <- serum_water_cors %>% filter(skipped)
if (nrow(skipped_sw) > 0) {
  cat("\nSkipped serum-water pairs:\n")
  print(skipped_sw %>% dplyr::select(serum, water, n, skip_reason))
}

cat("\nSerum PFAS vs Drinking Water PFAS:\n")
print(serum_water_cors %>% filter(!skipped) %>% dplyr::select(-skipped, -skip_reason))

# --- 2. Serum PFAS ~ eGFR and Albuminuria ---
cat("\nRunning serum PFAS vs outcome correlations...\n")

serum_outcome_cors <- expand_grid(
  pfas = usable_serum,
  outcome = usable_outcomes
) %>%
  rowwise() %>%
  mutate(
    cor_result = list(safe_cor_test(
      merged_data[[pfas]], 
      merged_data[[outcome]]
    ))
  ) %>%
  unnest(cor_result) %>%
  ungroup() %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05
  ) %>%
  arrange(p_value)

skipped_so <- serum_outcome_cors %>% filter(skipped)
if (nrow(skipped_so) > 0) {
  cat("\nSkipped serum-outcome pairs:\n")
  print(skipped_so %>% dplyr::select(pfas, outcome, n, skip_reason))
}

cat("\nSerum PFAS vs Kidney Outcomes:\n")
print(serum_outcome_cors %>% filter(!skipped) %>% dplyr::select(-skipped, -skip_reason))

# --- 3. Drinking Water PFAS ~ eGFR and Albuminuria ---
cat("\nRunning water PFAS vs outcome correlations...\n")

water_outcome_cors <- expand_grid(
  pfas = usable_water,
  outcome = usable_outcomes
) %>%
  rowwise() %>%
  mutate(
    cor_result = list(safe_cor_test(
      merged_data[[pfas]], 
      merged_data[[outcome]]
    ))
  ) %>%
  unnest(cor_result) %>%
  ungroup() %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05
  ) %>%
  arrange(p_value)

skipped_wo <- water_outcome_cors %>% filter(skipped)
if (nrow(skipped_wo) > 0) {
  cat("\nSkipped water-outcome pairs:\n")
  print(skipped_wo %>% dplyr::select(pfas, outcome, n, skip_reason))
}

cat("\nDrinking Water PFAS vs Kidney Outcomes:\n")
print(water_outcome_cors %>% filter(!skipped) %>% dplyr::select(-skipped, -skip_reason))

# --- Combine all results ---
all_correlations <- bind_rows(
  serum_water_cors %>% mutate(analysis = "Serum-Water"),
  serum_outcome_cors %>% rename(serum = pfas) %>% mutate(analysis = "Serum-Outcome"),
  water_outcome_cors %>% rename(water = pfas) %>% mutate(analysis = "Water-Outcome")
)

# --- Visualization: Heatmap of serum PFAS vs outcomes ---
serum_outcome_matrix <- serum_outcome_cors %>%
  dplyr::select(pfas, outcome, rho) %>%
  pivot_wider(names_from = outcome, values_from = rho) %>%
  column_to_rownames("pfas") %>%
  as.matrix()

corrplot(serum_outcome_matrix, method = "color", is.corr = FALSE,
         col = colorRampPalette(c("blue", "white", "red"))(100),
         tl.col = "black", tl.cex = 0.8,
         title = "Serum PFAS Correlations with Kidney Outcomes",
         mar = c(0, 0, 2, 0))

# --- Save results ---
write.csv(all_correlations, paste0(base.dir, "pfas_correlation_results.csv"), row.names = FALSE)

# --- PLOTS ---

# Create plots directory
plot_dir <- paste0(base.dir, "plots/")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Log-transform UACR for better visualization
merged_data <- merged_data %>%
  mutate(log_acr_u = log10(acr_u + 1))

# --- 1. Serum PFAS vs Drinking Water PFAS scatterplots ---
# Only plot pairs that weren't skipped
plottable_pairs <- usable_pairs %>%
  left_join(serum_water_cors %>% dplyr::select(serum, water, skipped), by = c("serum", "water")) %>%
  filter(!skipped)

if (nrow(plottable_pairs) > 0) {
  pdf(paste0(plot_dir, "serum_vs_water_pfas.pdf"), width = 12, height = 10)
  plots_serum_water <- plottable_pairs %>%
    rowwise() %>%
    mutate(
      plot = list({
        df <- merged_data %>% 
          dplyr::select(x = !!sym(water), y = !!sym(serum)) %>%
          filter(!is.na(x) & !is.na(y))
        
        cor_res <- safe_cor_test(df$x, df$y)
        
        ggplot(df, aes(x = x, y = y)) +
          geom_point(alpha = 0.5, color = "steelblue") +
          geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
          labs(
            x = paste("Water:", serum),
            y = paste("Serum:", serum),
            title = paste0(serum, ": Serum vs Drinking Water"),
            subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", cor_res$rho, cor_res$p_value, cor_res$n)
          ) +
          theme_bw() +
          theme(plot.title = element_text(face = "bold"))
      })
    )
  
  gridExtra::grid.arrange(grobs = plots_serum_water$plot, ncol = 3)
  dev.off()
  cat("  - serum_vs_water_pfas.pdf\n")
} else {
  cat("  - serum_vs_water_pfas.pdf SKIPPED (no plottable pairs)\n")
}

# --- 2. Serum PFAS vs eGFR scatterplots ---
plottable_serum_egfr <- serum_outcome_cors %>%
  filter(outcome == "eGFR_CKD_epi", !skipped) %>%
  pull(pfas)

if (length(plottable_serum_egfr) > 0) {
  pdf(paste0(plot_dir, "serum_pfas_vs_eGFR.pdf"), width = 12, height = 10)
  plots_serum_egfr <- map(plottable_serum_egfr, function(pfas) {
    df <- merged_data %>% 
      dplyr::select(x = !!sym(pfas), y = eGFR_CKD_epi) %>%
      filter(!is.na(x) & !is.na(y))
    
    cor_res <- safe_cor_test(df$x, df$y)
    
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
      labs(
        x = pfas,
        y = "eGFR (CKD-EPI)",
        title = paste0(pfas, " vs eGFR"),
        subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", cor_res$rho, cor_res$p_value, cor_res$n)
      ) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"))
  })
  
  gridExtra::grid.arrange(grobs = plots_serum_egfr, ncol = 3)
  dev.off()
  cat("  - serum_pfas_vs_eGFR.pdf\n")
} else {
  cat("  - serum_pfas_vs_eGFR.pdf SKIPPED (insufficient data)\n")
}

# --- 3. Serum PFAS vs UACR (log-transformed) scatterplots ---
plottable_serum_uacr <- serum_outcome_cors %>%
  filter(outcome == "acr_u", !skipped) %>%
  pull(pfas)

if (length(plottable_serum_uacr) > 0) {
  pdf(paste0(plot_dir, "serum_pfas_vs_UACR.pdf"), width = 12, height = 10)
  plots_serum_uacr <- map(plottable_serum_uacr, function(pfas) {
    df <- merged_data %>% 
      dplyr::select(x = !!sym(pfas), y = log_acr_u) %>%
      filter(!is.na(x) & !is.na(y) & is.finite(y))
    
    cor_res <- safe_cor_test(df$x, df$y)
    
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "darkorange") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
      labs(
        x = pfas,
        y = "log10(UACR + 1)",
        title = paste0(pfas, " vs UACR"),
        subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", cor_res$rho, cor_res$p_value, cor_res$n)
      ) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"))
  })
  
  gridExtra::grid.arrange(grobs = plots_serum_uacr, ncol = 3)
  dev.off()
  cat("  - serum_pfas_vs_UACR.pdf\n")
} else {
  cat("  - serum_pfas_vs_UACR.pdf SKIPPED (insufficient data)\n")
}

# --- 4. Water PFAS vs eGFR (total and major compounds) ---
water_pfas_main <- c("total_pfas", "pfos_and_pfoa_total_combined",
                     "perfluorooctanesulfonic_acid_pfos_", 
                     "perfluorooctanoic_acid_pfoa_",
                     "perfluorohexanesulfonic_acid_pfhxs_",
                     "perfluorononanoic_acid_pfna_")

plottable_water_egfr <- water_outcome_cors %>%
  filter(outcome == "eGFR_CKD_epi", !skipped, pfas %in% water_pfas_main) %>%
  pull(pfas)

if (length(plottable_water_egfr) > 0) {
  pdf(paste0(plot_dir, "water_pfas_vs_eGFR.pdf"), width = 12, height = 8)
  plots_water_egfr <- map(plottable_water_egfr, function(pfas) {
    df <- merged_data %>% 
      dplyr::select(x = !!sym(pfas), y = eGFR_CKD_epi) %>%
      filter(!is.na(x) & !is.na(y))
    
    cor_res <- safe_cor_test(df$x, df$y)
    short_name <- str_extract(pfas, "(?<=_)[a-z]+(?=_$)|total_pfas|pfos_and_pfoa_total_combined") %>%
      toupper() %>% coalesce(pfas)
    
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "forestgreen") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
      labs(
        x = paste("Water:", short_name),
        y = "eGFR (CKD-EPI)",
        title = paste0("Water ", short_name, " vs eGFR"),
        subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", cor_res$rho, cor_res$p_value, cor_res$n)
      ) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"))
  })
  
  gridExtra::grid.arrange(grobs = plots_water_egfr, ncol = 3)
  dev.off()
  cat("  - water_pfas_vs_eGFR.pdf\n")
} else {
  cat("  - water_pfas_vs_eGFR.pdf SKIPPED (insufficient data)\n")
}

# --- 5. Water PFAS vs UACR (log-transformed) ---
plottable_water_uacr <- water_outcome_cors %>%
  filter(outcome == "acr_u", !skipped, pfas %in% water_pfas_main) %>%
  pull(pfas)

if (length(plottable_water_uacr) > 0) {
  pdf(paste0(plot_dir, "water_pfas_vs_UACR.pdf"), width = 12, height = 8)
  plots_water_uacr <- map(plottable_water_uacr, function(pfas) {
    df <- merged_data %>% 
      dplyr::select(x = !!sym(pfas), y = log_acr_u) %>%
      filter(!is.na(x) & !is.na(y) & is.finite(y))
    
    cor_res <- safe_cor_test(df$x, df$y)
    short_name <- str_extract(pfas, "(?<=_)[a-z]+(?=_$)|total_pfas|pfos_and_pfoa_total_combined") %>%
      toupper() %>% coalesce(pfas)
    
    ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "purple") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
      labs(
        x = paste("Water:", short_name),
        y = "log10(UACR + 1)",
        title = paste0("Water ", short_name, " vs UACR"),
        subtitle = sprintf("rho = %.3f, p = %.3g, n = %d", cor_res$rho, cor_res$p_value, cor_res$n)
      ) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"))
  })
  
  gridExtra::grid.arrange(grobs = plots_water_uacr, ncol = 3)
  dev.off()
  cat("  - water_pfas_vs_UACR.pdf\n")
} else {
  cat("  - water_pfas_vs_UACR.pdf SKIPPED (insufficient data)\n")
}

# --- 6. Correlation heatmap (only non-skipped) ---
heatmap_cors <- serum_outcome_cors %>% filter(!skipped)

if (nrow(heatmap_cors) > 0) {
  pdf(paste0(plot_dir, "correlation_heatmap.pdf"), width = 10, height = 8)
  
  # Build matrix for heatmap
  heatmap_data <- heatmap_cors %>%
    mutate(outcome = case_when(
      outcome == "eGFR_CKD_epi" ~ "eGFR",
      outcome == "acr_u" ~ "UACR",
      TRUE ~ outcome
    )) %>%
    dplyr::select(pfas, outcome, rho) %>%
    pivot_wider(names_from = outcome, values_from = rho) %>%
    column_to_rownames("pfas") %>%
    as.matrix()
  
  corrplot(heatmap_data, method = "color", is.corr = FALSE,
           col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           tl.col = "black", tl.cex = 1,
           cl.cex = 0.8,
           addCoef.col = "black", number.cex = 0.8,
           title = "Serum PFAS Correlations with Kidney Outcomes\n(Spearman's rho)",
           mar = c(0, 0, 3, 0))
  
  dev.off()
  cat("  - correlation_heatmap.pdf\n")
} else {
  cat("  - correlation_heatmap.pdf SKIPPED (insufficient data)\n")
}

# --- Summary of skipped analyses ---
all_skipped <- bind_rows(
  skipped_sw %>% mutate(analysis = "Serum-Water"),
  skipped_so %>% mutate(analysis = "Serum-Outcome"),
  skipped_wo %>% mutate(analysis = "Water-Outcome")
)

if (nrow(all_skipped) > 0) {
  cat("\n=== SKIPPED ANALYSES SUMMARY ===\n")
  cat("Total skipped:", nrow(all_skipped), "\n")
  print(all_skipped %>% dplyr::select(analysis, everything(), -comparison))
  write.csv(all_skipped, paste0(base.dir, "skipped_analyses.csv"), row.names = FALSE)
  cat("\nSkipped analyses saved to: skipped_analyses.csv\n")
}

cat("\nPlots saved to:", plot_dir, "\n")
cat("Files created:\n")
cat("  - serum_vs_water_pfas.pdf\n")
cat("  - serum_pfas_vs_eGFR.pdf\n")
cat("  - serum_pfas_vs_UACR.pdf\n")
cat("  - water_pfas_vs_eGFR.pdf\n")
cat("  - water_pfas_vs_UACR.pdf\n")
cat("  - correlation_heatmap.pdf\n")
























