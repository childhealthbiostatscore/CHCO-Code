# ============================================================================
# COMPLETE FSOC PROTEOMICS ANALYSIS PIPELINE
# From data loading to immune-focused visualizations with real protein names
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
# STEP 1: LOAD SOMASCAN ANNOTATION FILE (THE KEY!)
# ============================================================================

cat("========================================\n")
cat("STEP 1: LOADING SOMASCAN ANNOTATION\n")
cat("========================================\n\n")

# UPDATE THIS PATH to your SomaScan annotation file
# It should have columns: SeqId, TargetFullName, UniProt, EntrezGeneSymbol
ANNOTATION_FILE <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Combined SomaScan/analytes_2.RData"

# Try loading as CSV first
if (file.exists(ANNOTATION_FILE)) {
  if (str_detect(ANNOTATION_FILE, "\\.xlsx$")) {
    somascan_annotation <- readxl::read_excel(ANNOTATION_FILE)
  } else {
    somascan_annotation <- load(ANNOTATION_FILE)
  }
  
  cat("✓ Loaded annotation file:", nrow(somascan_annotation), "entries\n")
  cat("Columns:", paste(names(somascan_annotation), collapse = ", "), "\n\n")
  
  # Show sample
  cat("Sample entries:\n")
  print(head(somascan_annotation, 10))
} else {
  stop("\n❌ ERROR: Annotation file not found!\n",
       "Please update ANNOTATION_FILE path at line 30\n",
       "Expected file with columns: SeqId, TargetFullName, UniProt, EntrezGeneSymbol\n")
}

# ============================================================================
# STEP 2: LOAD HARMONIZED DATA
# ============================================================================

cat("\n========================================\n")
cat("STEP 2: LOADING HARMONIZED DATA\n")
cat("========================================\n\n")

harmonized_data <- read.csv(
  "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", 
  na = ''
)

cat("Data loaded:", nrow(harmonized_data), "rows x", ncol(harmonized_data), "columns\n")

# Find plasma proteomics columns (seq.NUMBER.NUMBER)
all_cols <- names(harmonized_data)
plasma_seq_cols <- all_cols[str_detect(all_cols, "^seq\\.[0-9]+\\.[0-9]+$")]

cat("Found", length(plasma_seq_cols), "plasma proteomics columns\n\n")

# ============================================================================
# STEP 3: PREPARE DATA AND CALCULATE FSOC
# ============================================================================

cat("========================================\n")
cat("STEP 3: PREPARING DATA\n")
cat("========================================\n\n")

# Aggregate by subject
dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

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

# Filter to subjects with proteomics
dat_fsoc_proteomics <- dat_fsoc %>%
  filter(rowSums(!is.na(select(., any_of(plasma_seq_cols)))) > 10)

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
# STEP 4: DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

cat("\n========================================\n")
cat("STEP 4: DIFFERENTIAL EXPRESSION\n")
cat("========================================\n\n")

results <- data.frame()

for(var in plasma_seq_cols) {
  
  if (!var %in% names(dat_fsoc_proteomics)) next
  
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

cat("Proteins tested:", nrow(results), "\n")
cat("Significant (p < 0.05):", sum(results$p_value < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR < 0.10):", sum(results$p_fdr < 0.10, na.rm = TRUE), "\n\n")

# ============================================================================
# STEP 5: ANNOTATE WITH REAL PROTEIN NAMES
# ============================================================================

cat("========================================\n")
cat("STEP 5: ANNOTATING WITH REAL NAMES\n")
cat("========================================\n\n")

# Extract SeqId from column names
results <- results %>%
  mutate(
    seq_id = str_remove(variable, "^seq\\."),
    seq_id_dash = str_replace(seq_id, "\\.", "-")
  )

# Try matching (adjust based on your annotation file format)
results_annotated <- results %>%
  left_join(
    somascan_annotation,
    by = c("seq_id_dash" = "SeqId")
  )

# If that didn't work, try alternative formats
if (sum(!is.na(results_annotated$TargetFullName)) < 100) {
  cat("Trying alternative SeqId format...\n")
  results_annotated <- results %>%
    left_join(
      somascan_annotation,
      by = c("seq_id" = "SeqId")
    )
}

n_matched <- sum(!is.na(results_annotated$TargetFullName))
cat("Matched", n_matched, "of", nrow(results), "proteins to annotation\n\n")

if (n_matched < 100) {
  cat("⚠️ WARNING: Few matches. Check SeqId format in annotation file.\n")
  cat("Column names in annotation:", paste(names(somascan_annotation), collapse = ", "), "\n")
}

# Save full results
write.csv(results_annotated, 
          file.path(OUTPUT_DIR, "all_proteins_results_annotated.csv"),
          row.names = FALSE)

# ============================================================================
# STEP 6: IDENTIFY IMMUNE/INFLAMMATORY PROTEINS
# ============================================================================

cat("========================================\n")
cat("STEP 6: IDENTIFYING IMMUNE PROTEINS\n")
cat("========================================\n\n")

# COMPREHENSIVE immune/inflammatory keywords
immune_keywords <- c(
  # Interleukins
  "interleukin", "IL-", "IL1", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL8", 
  "IL9", "IL10", "IL11", "IL12", "IL13", "IL15", "IL17", "IL18", "IL21", "IL23", "IL33",
  
  # TNF family
  "TNF", "tumor necrosis", "TNFR", "TNFSF", "TRAIL", "TWEAK", "lymphotoxin", "FAS",
  
  # Chemokines
  "chemokine", "CCL", "CXCL", "CX3CL", "MCP", "MIP", "RANTES", "eotaxin", "IP-10",
  "fractalkine", "SDF", "GRO", "ENA",
  
  # Inflammatory markers
  "C-reactive", "CRP", "serum amyloid", "SAA", "pentraxin", "PTX",
  "calprotectin", "S100", "HMGB1", "ferritin",
  
  # Complement
  "complement", " C3", " C4", " C5", "factor B", "factor D", "factor H", "properdin",
  
  # Adhesion molecules
  "ICAM", "VCAM", "selectin", "SELE", "SELL", "SELP", "integrin", "PECAM", "cadherin",
  
  # Growth factors
  "VEGF", "PDGF", "FGF", "TGF", "EGF", "HGF", "IGF", "BDNF", "NGF",
  
  # Matrix & proteases
  "metalloproteinase", "MMP", "TIMP", "ADAM", "cathepsin", "elastase",
  
  # Cell markers
  "CD40", "CD14", "CD163", "CD28", "CD80", "CD86", "PD-1", "PD-L1",
  "CD3", "CD4", "CD8", "CD11", "CD19", "CD20", "CD25",
  
  # Interferons
  "interferon", "IFN",
  
  # Colony stimulating
  "CSF", "G-CSF", "GM-CSF", "M-CSF",
  
  # Kidney injury
  "kidney injury", "KIM", "NGAL", "lipocalin", "cystatin", "clusterin", 
  "uromodulin", "nephrin",
  
  # Oxidative stress
  "myeloperoxidase", "MPO", "NADPH oxidase", "superoxide",
  
  # Adipokines
  "osteopontin", "resistin", "adiponectin", "leptin", "visfatin"
)

# Search in TargetFullName AND EntrezGeneSymbol
immune_results <- results_annotated %>%
  filter(!is.na(TargetFullName)) %>%
  filter(
    str_detect(tolower(TargetFullName), paste(tolower(immune_keywords), collapse = "|")) |
      str_detect(tolower(EntrezGeneSymbol), paste(tolower(immune_keywords), collapse = "|"))
  )

cat("Immune proteins found:", nrow(immune_results), "\n")
cat("Significant (p < 0.05):", sum(immune_results$p_value < 0.05), "\n")
cat("Significant (FDR < 0.10):", sum(immune_results$p_fdr < 0.10), "\n\n")

cat("Top 20 immune proteins by p-value:\n")
print(immune_results %>% 
        select(EntrezGeneSymbol, TargetFullName, p_value, p_fdr, fold_change) %>%
        arrange(p_value) %>%
        head(20))

write.csv(immune_results, 
          file.path(OUTPUT_DIR, "immune_proteins_results.csv"),
          row.names = FALSE)

# ============================================================================
# STEP 7: VOLCANO PLOT - IMMUNE PROTEINS
# ============================================================================

cat("\n========================================\n")
cat("STEP 7: CREATING VOLCANO PLOT\n")
cat("========================================\n\n")

sig_immune <- immune_results %>% filter(p_fdr < 0.10)
top_immune <- immune_results %>% filter(p_value < 0.05) %>% arrange(p_value)

p_volcano <- ggplot(immune_results, 
                    aes(x = log2_fc, y = -log10(p_value))) +
  geom_point(aes(color = p_fdr < 0.10, size = -log10(p_value)), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "blue", alpha = 0.5, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#E41A1C"),
                     name = "FDR < 0.10") +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  theme_bw(base_size = 12) +
  labs(x = "Log2 Fold Change (Impaired / Normal FSOC)",
       y = "-log10(P-value)",
       title = "Immune-Related Proteins: Impaired vs Normal Medullary FSOC",
       subtitle = paste0(nrow(immune_results), " immune proteins, ",
                         sum(immune_results$p_value < 0.05), " p<0.05, ",
                         sum(immune_results$p_fdr < 0.10), " FDR<0.10")) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

# Add labels
if (nrow(sig_immune) > 0) {
  p_volcano <- p_volcano +
    geom_text_repel(data = sig_immune,
                    aes(label = EntrezGeneSymbol), 
                    size = 3, fontface = "bold",
                    max.overlaps = 30, segment.size = 0.3)
} else if (nrow(top_immune) > 0) {
  label_data <- head(top_immune, 20)
  p_volcano <- p_volcano +
    geom_text_repel(data = label_data,
                    aes(label = EntrezGeneSymbol), 
                    size = 2.5, max.overlaps = 20)
}

ggsave(file.path(OUTPUT_DIR, "volcano_immune_proteins.pdf"), 
       p_volcano, width = 12, height = 10)
ggsave(file.path(OUTPUT_DIR, "volcano_immune_proteins.png"), 
       p_volcano, width = 12, height = 10, dpi = 300)

cat("Volcano plot saved.\n")

# ============================================================================
# STEP 8: BOXPLOTS OF TOP IMMUNE PROTEINS
# ============================================================================

cat("\n========================================\n")
cat("STEP 8: CREATING BOXPLOTS\n")
cat("========================================\n\n")

# Get top proteins
if (sum(immune_results$p_fdr < 0.10) >= 12) {
  top_proteins <- immune_results %>% 
    filter(p_fdr < 0.10) %>%
    arrange(p_value) %>%
    head(12)
} else {
  top_proteins <- immune_results %>% 
    arrange(p_value) %>%
    head(min(12, nrow(immune_results)))
}

if (nrow(top_proteins) > 0) {
  plot_list <- list()
  
  for(i in 1:nrow(top_proteins)) {
    var <- top_proteins$variable[i]
    gene <- top_proteins$EntrezGeneSymbol[i]
    full_name <- top_proteins$TargetFullName[i]
    
    plot_data <- dat_fsoc_proteomics %>%
      select(fsoc_binary, all_of(var)) %>%
      rename(value = all_of(var)) %>%
      filter(!is.na(value), !is.na(fsoc_binary))
    
    # Create title with gene symbol
    if (!is.na(gene) && gene != "" && gene != "NA") {
      title <- paste0(gene, ": ", str_trunc(full_name, 35))
    } else {
      title <- str_trunc(full_name, 45)
    }
    
    p_label <- sprintf("p=%.2e, FDR=%.3f, FC=%.2f", 
                       top_proteins$p_value[i], 
                       top_proteins$p_fdr[i], 
                       top_proteins$fold_change[i])
    
    p <- ggplot(plot_data, aes(x = fsoc_binary, y = value, fill = fsoc_binary)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.6) +
      geom_jitter(width = 0.25, alpha = 0.6, size = 2) +
      scale_fill_manual(values = c("Normal" = "#3182BD", "Impaired" = "#E6550D")) +
      theme_bw(base_size = 10) +
      labs(title = title,
           subtitle = p_label,
           x = "FSOC Status",
           y = "Expression (RFU)") +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold"),
            plot.subtitle = element_text(size = 7.5),
            panel.grid.major.x = element_blank())
    
    plot_list[[i]] <- p
  }
  
  png(file.path(OUTPUT_DIR, 'top_immune_proteins_boxplots.png'), 
      height = 15, width = 16, units = 'in', res = 300)
  grid.arrange(grobs = plot_list, ncol = 3,
               top = grid::textGrob("Top Immune Proteins by FSOC Status", 
                                    gp = grid::gpar(fontsize = 16, fontface = "bold")))
  dev.off()
  
  cat("Boxplots saved.\n")
}

# ============================================================================
# STEP 9: SUMMARY TABLE
# ============================================================================

cat("\n========================================\n")
cat("STEP 9: CREATING SUMMARY TABLE\n")
cat("========================================\n\n")

summary_table <- immune_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  mutate(
    Direction = ifelse(log2_fc > 0, "↑ Impaired", "↓ Impaired"),
    `P-value` = formatC(p_value, format = "e", digits = 2),
    `FDR` = formatC(p_fdr, format = "f", digits = 3),
    `Fold Change` = formatC(fold_change, format = "f", digits = 2)
  ) %>%
  select(EntrezGeneSymbol, TargetFullName, Direction, `P-value`, FDR, `Fold Change`) %>%
  head(30)

write.csv(summary_table, 
          file.path(OUTPUT_DIR, "immune_proteins_summary_table.csv"),
          row.names = FALSE)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("========================================\n\n")

cat("COHORT:\n")
cat("  Total subjects:", nrow(dat_fsoc_proteomics), "\n")
cat("  Impaired FSOC:", sum(dat_fsoc_proteomics$fsoc_binary == "Impaired"), "\n")
cat("  Normal FSOC:", sum(dat_fsoc_proteomics$fsoc_binary == "Normal"), "\n\n")

cat("PROTEINS:\n")
cat("  Total tested:", nrow(results_annotated), "\n")
cat("  With annotations:", sum(!is.na(results_annotated$TargetFullName)), "\n")
cat("  Immune-related:", nrow(immune_results), "\n\n")

cat("SIGNIFICANCE:\n")
cat("  All proteins p<0.05:", sum(results_annotated$p_value < 0.05, na.rm = TRUE), "\n")
cat("  All proteins FDR<0.10:", sum(results_annotated$p_fdr < 0.10, na.rm = TRUE), "\n")
cat("  Immune proteins p<0.05:", sum(immune_results$p_value < 0.05), "\n")
cat("  Immune proteins FDR<0.10:", sum(immune_results$p_fdr < 0.10), "\n\n")

cat("FILES CREATED in", OUTPUT_DIR, ":\n")
cat("  1. all_proteins_results_annotated.csv\n")
cat("  2. immune_proteins_results.csv\n")
cat("  3. volcano_immune_proteins.pdf/.png\n")
if (nrow(top_proteins) > 0) {
  cat("  4. top_immune_proteins_boxplots.png\n")
}
cat("  5. immune_proteins_summary_table.csv\n")

cat("\n========================================\n")



