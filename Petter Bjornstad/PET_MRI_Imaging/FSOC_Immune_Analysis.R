# ============================================================================
# FSOC ANALYSIS WITH SOMASCAN PROTEOMICS - FINAL CORRECT VERSION
# Uses seq.NUMBER.NUMBER format (periods, not underscores!)
# ============================================================================

library(tidyverse)
library(readxl)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(ggrepel)

# ============================================================================
# SETUP
# ============================================================================

OUTPUT_DIR <- "C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani/FSOC_Proteomics"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("=== LOADING SOMA HARMONIZED DATASET ===\n")

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", na = '')

data_dictionary <- readxl::read_xlsx('C:/Users/netio/Downloads/data_dictionary_master.xlsx')

cat("Data loaded:", nrow(harmonized_data), "rows x", ncol(harmonized_data), "columns\n")

# Find seq. columns (with PERIODS!)
all_cols <- names(harmonized_data)

# Get plasma proteomics (seq.NUMBER.NUMBER with no suffix)
plasma_seq_cols <- all_cols[str_detect(all_cols, "^seq\\.[0-9]+\\.[0-9]+$")]

# Get urine proteomics
urine_seq_cols <- all_cols[str_detect(all_cols, "^seq\\.[0-9]+\\.[0-9]+_urine$")]
urine_cradj_cols <- all_cols[str_detect(all_cols, "^seq\\.[0-9]+\\.[0-9]+_urine_cradj$")]

cat("\nFound proteomics columns:\n")
cat("  Plasma:", length(plasma_seq_cols), "\n")
cat("  Urine:", length(urine_seq_cols), "\n")
cat("  Urine (Cr-adjusted):", length(urine_cradj_cols), "\n")

# Use plasma proteomics for analysis
seq_cols <- plasma_seq_cols

cat("\nUsing", length(seq_cols), "plasma proteomics columns\n")
cat("Sample columns:\n")
print(head(seq_cols, 20))

# ============================================================================
# 2. MAP SEQ IDs TO PROTEIN NAMES
# ============================================================================

cat("\n=== MAPPING SEQ IDs TO PROTEIN NAMES ===\n")

# Get proteomics dictionary
proteomics_dict <- data_dictionary %>%
  filter(form_name == 'proteomics')

cat("Proteomics entries in dictionary:", nrow(proteomics_dict), "\n")

# Convert data column names to dictionary format
# seq.10000.28 -> seq_10000_28
seq_mapping <- data.frame(
  column_name = seq_cols,
  dict_name = str_replace_all(seq_cols, "\\.", "_"),
  stringsAsFactors = FALSE
) %>%
  left_join(
    proteomics_dict %>% select(variable_name, label),
    by = c("dict_name" = "variable_name")
  )

cat("Successfully mapped:", sum(!is.na(seq_mapping$label)), "of", nrow(seq_mapping), "proteins\n\n")

cat("Sample mappings:\n")
print(head(seq_mapping %>% filter(!is.na(label)) %>% select(column_name, label), 20))

write.csv(seq_mapping, file.path(OUTPUT_DIR, "seq_to_protein_mapping.csv"), 
          row.names = FALSE)

# ============================================================================
# 3. IDENTIFY INFLAMMATORY/IMMUNE PROTEINS
# ============================================================================

cat("\n=== IDENTIFYING INFLAMMATORY/IMMUNE PROTEINS ===\n")

inflammatory_keywords <- c(
  "interleukin", "IL-", "IL6", "IL1", "IL8", "IL10", "IL12", "IL17", "IL18",
  "TNF", "tumor necrosis",
  "chemokine", "CCL", "CXCL", "MCP", "MIP", "RANTES",
  "C-reactive", "CRP",
  "complement", " C3", " C5", " C4",
  "adhesion", "ICAM", "VCAM", "selectin", "SELE", "SELL",
  "myeloperoxidase", "MPO",
  "kidney injury", "KIM", "NGAL", "lipocalin",
  "cystatin", "clusterin",
  "uromodulin", "osteopontin", "pentraxin",
  "haptoglobin", "serum amyloid", "SAA",
  "fibrinogen", "interferon", "IFN",
  "CD40", "CD14", "CD163",
  "matrix metalloproteinase", "MMP"
)

inflammatory_proteins <- seq_mapping %>%
  filter(!is.na(label)) %>%
  filter(str_detect(tolower(label), 
                    paste(tolower(inflammatory_keywords), collapse = "|")))

cat("Found", nrow(inflammatory_proteins), "inflammatory/immune proteins\n\n")

if (nrow(inflammatory_proteins) > 0) {
  cat("Inflammatory proteins found:\n")
  print(inflammatory_proteins %>% select(column_name, label))
  
  write.csv(inflammatory_proteins, 
            file.path(OUTPUT_DIR, "inflammatory_proteins_list.csv"),
            row.names = FALSE)
} else {
  cat("WARNING: No inflammatory proteins found with keywords!\n")
  cat("Will use all proteomics for analysis.\n")
}

# ============================================================================
# 4. PREPARE DATA FOR FSOC ANALYSIS
# ============================================================================

cat("\n=== PREPARING DATA ===\n")

# Aggregate by subject
dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

# Check FSOC columns
cat("FSOC columns available:\n")
fsoc_cols_found <- names(dat)[str_detect(tolower(names(dat)), "fsoc")]
print(fsoc_cols_found)

# Calculate medullary FSOC
dat <- dat %>%
  mutate(
    medullary_fsoc = case_when(
      !is.na(avg_m_fsoc) ~ avg_m_fsoc,
      !is.na(fsoc_l_medulla) & !is.na(fsoc_r_medulla) ~ (fsoc_l_medulla + fsoc_r_medulla) / 2,
      !is.na(fsoc_l_medulla) ~ fsoc_l_medulla,
      !is.na(fsoc_r_medulla) ~ fsoc_r_medulla,
      TRUE ~ NA_real_
    ),
    whole_kidney_fsoc = case_when(
      !is.na(avg_k_fsoc) ~ avg_k_fsoc,
      !is.na(fsoc_l_kidney) & !is.na(fsoc_r_kidney) ~ (fsoc_l_kidney + fsoc_r_kidney) / 2,
      !is.na(fsoc_l_kidney) ~ fsoc_l_kidney,
      !is.na(fsoc_r_kidney) ~ fsoc_r_kidney,
      TRUE ~ NA_real_
    )
  )

# Filter to valid FSOC
dat_fsoc <- dat %>%
  filter(visit == 'baseline') %>%
  filter(!is.na(medullary_fsoc)) %>%
  filter(medullary_fsoc >= 0) %>%
  filter(whole_kidney_fsoc < 15 | is.na(whole_kidney_fsoc))

cat("Subjects with valid FSOC:", nrow(dat_fsoc), "\n")

# Decide which proteins to test
proteomics_cols <- if (nrow(inflammatory_proteins) >= 10) {
  inflammatory_proteins$column_name
} else {
  seq_cols
}

# Filter to subjects with proteomics
dat_fsoc_proteomics <- dat_fsoc %>%
  filter(rowSums(!is.na(select(., any_of(proteomics_cols)))) > 10)

cat("Subjects with FSOC and proteomics:", nrow(dat_fsoc_proteomics), "\n")

# Classify FSOC
dat_fsoc_proteomics <- dat_fsoc_proteomics %>%
  mutate(
    fsoc_binary = ifelse(medullary_fsoc < median(medullary_fsoc, na.rm = TRUE), 
                         "Impaired", "Normal"),
    fsoc_binary = factor(fsoc_binary, levels = c("Normal", "Impaired"))
  )

cat("\nFSOC Classification:\n")
print(table(dat_fsoc_proteomics$fsoc_binary))

# ============================================================================
# 5. DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

results <- data.frame()

for(var in proteomics_cols) {
  
  if (!var %in% names(dat_fsoc_proteomics)) next
  
  # Get protein name
  protein_name <- seq_mapping$label[seq_mapping$column_name == var]
  if (length(protein_name) == 0 || is.na(protein_name)) protein_name <- var
  
  temp_data <- dat_fsoc_proteomics %>%
    select(fsoc_binary, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value), !is.na(fsoc_binary))
  
  if(nrow(temp_data) < 10) next
  
  impaired <- temp_data$value[temp_data$fsoc_binary == "Impaired"]
  normal <- temp_data$value[temp_data$fsoc_binary == "Normal"]
  
  if(length(impaired) < 3 || length(normal) < 3) next
  
  test_result <- wilcox.test(value ~ fsoc_binary, data = temp_data)
  
  results <- rbind(results, data.frame(
    variable = var,
    protein_name = protein_name,
    p_value = test_result$p.value,
    mean_impaired = mean(impaired, na.rm = TRUE),
    mean_normal = mean(normal, na.rm = TRUE),
    median_impaired = median(impaired, na.rm = TRUE),
    median_normal = median(normal, na.rm = TRUE),
    n_impaired = length(impaired),
    n_normal = length(normal),
    stringsAsFactors = FALSE
  ))
}

# Adjust for multiple testing
results <- results %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "fdr"),
    fold_change = mean_impaired / mean_normal,
    log2_fc = log2(fold_change),
    difference = mean_impaired - mean_normal
  ) %>%
  arrange(p_value)

cat("\nProteins tested:", nrow(results), "\n")
cat("Significant (p < 0.05):", sum(results$p_value < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR < 0.05):", sum(results$p_fdr < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR < 0.10):", sum(results$p_fdr < 0.10, na.rm = TRUE), "\n\n")

cat("Top 30 proteins:\n")
print(head(results %>% select(protein_name, p_value, p_fdr, fold_change), 30))

write.csv(results, 
          file.path(OUTPUT_DIR, "fsoc_proteomics_results_FINAL.csv"),
          row.names = FALSE)

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================

if (nrow(results) > 0) {
  
  # Volcano plot
  sig_proteins <- results %>% filter(p_fdr < 0.10)
  
  p_volcano <- ggplot(results, 
                      aes(x = log2_fc, y = -log10(p_value), 
                          color = p_fdr < 0.10)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
               color = "blue", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"),
                       name = "FDR < 0.10") +
    theme_bw(base_size = 12) +
    labs(x = "Log2 Fold Change (Impaired / Normal FSOC)",
         y = "-log10(P-value)",
         title = "Proteomics: Impaired vs Normal Medullary FSOC",
         subtitle = paste0("n=", nrow(dat_fsoc_proteomics), " subjects, ", 
                           nrow(results), " proteins tested")) +
    theme(legend.position = "top")
  
  if (nrow(sig_proteins) > 0 && nrow(sig_proteins) <= 40) {
    p_volcano <- p_volcano +
      geom_text_repel(data = sig_proteins,
                      aes(label = protein_name), 
                      size = 2, max.overlaps = 40, segment.size = 0.2)
  }
  
  ggsave(file.path(OUTPUT_DIR, "volcano_plot_FINAL.pdf"), 
         p_volcano, width = 14, height = 10)
  ggsave(file.path(OUTPUT_DIR, "volcano_plot_FINAL.png"), 
         p_volcano, width = 14, height = 10, dpi = 300)
  
  cat("\nVolcano plot saved.\n")
  
  # Boxplots
  if (sum(results$p_value < 0.05, na.rm = TRUE) > 0) {
    top_n <- min(12, sum(results$p_value < 0.05, na.rm = TRUE))
    top_vars <- head(results$variable, top_n)
    
    plot_list <- list()
    
    for(i in 1:length(top_vars)) {
      var <- top_vars[i]
      var_info <- results %>% filter(variable == var)
      
      plot_data <- dat_fsoc_proteomics %>%
        select(fsoc_binary, all_of(var)) %>%
        rename(value = all_of(var)) %>%
        filter(!is.na(value), !is.na(fsoc_binary))
      
      display_name <- if(nchar(var_info$protein_name) > 45) {
        paste0(substr(var_info$protein_name, 1, 42), "...")
      } else {
        var_info$protein_name
      }
      
      p_label <- sprintf("p=%.4f, FDR=%.3f, FC=%.2f", 
                         var_info$p_value, var_info$p_fdr, var_info$fold_change)
      
      p <- ggplot(plot_data, aes(x = fsoc_binary, y = value, fill = fsoc_binary)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
        scale_fill_manual(values = c("Impaired" = "#E6550D", "Normal" = "#3182BD")) +
        theme_bw(base_size = 9) +
        labs(title = display_name,
             subtitle = p_label,
             x = "FSOC Status",
             y = "RFU") +
        theme(legend.position = "none",
              plot.title = element_text(size = 8, face = "bold"),
              plot.subtitle = element_text(size = 6.5))
      
      plot_list[[i]] <- p
    }
    
    png(file.path(OUTPUT_DIR, 'top_proteins_boxplots_FINAL.png'), 
        height = 15, width = 14, units = 'in', res = 300)
    grid.arrange(grobs = plot_list, ncol = 3)
    dev.off()
    
    cat("Boxplots saved.\n")
  }
}

# ============================================================================
# 7. SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", OUTPUT_DIR, "\n\n")

cat("=== SUMMARY ===\n")
cat("Total subjects analyzed:", nrow(dat_fsoc_proteomics), "\n")
cat("  Impaired FSOC:", sum(dat_fsoc_proteomics$fsoc_binary == "Impaired"), "\n")
cat("  Normal FSOC:", sum(dat_fsoc_proteomics$fsoc_binary == "Normal"), "\n\n")

cat("Proteins tested:", nrow(results), "\n")
cat("Significant (p < 0.05):", sum(results$p_value < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR < 0.10):", sum(results$p_fdr < 0.10, na.rm = TRUE), "\n\n")

if (sum(results$p_fdr < 0.10, na.rm = TRUE) > 0) {
  cat("=== TOP SIGNIFICANT PROTEINS (FDR < 0.10) ===\n")
  sig_results <- results %>% filter(p_fdr < 0.10)
  print(sig_results %>% 
          select(protein_name, fold_change, p_value, p_fdr) %>%
          head(20))
} else {
  cat("=== TOP PROTEINS BY P-VALUE ===\n")
  print(head(results %>% select(protein_name, fold_change, p_value, p_fdr), 20))
}

cat("\n\nFILES CREATED:\n")
cat("1. seq_to_protein_mapping.csv - Protein name mappings\n")
if (nrow(inflammatory_proteins) > 0) {
  cat("2. inflammatory_proteins_list.csv - Inflammatory markers\n")
}
cat("3. fsoc_proteomics_results_FINAL.csv - Full results\n")
cat("4. volcano_plot_FINAL.pdf/.png - Volcano plot\n")
if (sum(results$p_value < 0.05, na.rm = TRUE) > 0) {
  cat("5. top_proteins_boxplots_FINAL.png - Boxplots\n")
}

cat("\n=== DONE! ===\n")



































