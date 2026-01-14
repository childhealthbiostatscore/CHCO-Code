# ============================================================================
# FSOC IMMUNE-FOCUSED PROTEOMICS ANALYSIS
# Enhanced visualizations with real protein names and immune enrichment
# ============================================================================

library(tidyverse)
library(readxl)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(ggrepel)
library(enrichR)

OUTPUT_DIR <- "C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani/FSOC_Proteomics"

# ============================================================================
# 1. LOAD PREVIOUS RESULTS
# ============================================================================

cat("=== LOADING RESULTS ===\n")

# Load the results from previous analysis
results <- read.csv(file.path(OUTPUT_DIR, "fsoc_proteomics_results_FINAL.csv"))
seq_mapping <- read.csv(file.path(OUTPUT_DIR, "seq_to_protein_mapping.csv"))

cat("Total proteins analyzed:", nrow(results), "\n")
cat("Proteins with names:", sum(!is.na(results$protein_name)), "\n\n")

# Load harmonized data for plotting
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", na = '')

data_dictionary <- readxl::read_xlsx('C:/Users/netio/Downloads/data_dictionary_master.xlsx')

# ============================================================================
# 2. PREPARE DATA (same as before)
# ============================================================================

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  ) %>%
  mutate(
    medullary_fsoc = case_when(
      !is.na(avg_m_fsoc) ~ avg_m_fsoc,
      !is.na(fsoc_l_medulla) & !is.na(fsoc_r_medulla) ~ (fsoc_l_medulla + fsoc_r_medulla) / 2,
      !is.na(fsoc_l_medulla) ~ fsoc_l_medulla,
      !is.na(fsoc_r_medulla) ~ fsoc_r_medulla,
      TRUE ~ NA_real_
    )
  )

dat_fsoc_proteomics <- dat %>%
  filter(visit == 'baseline') %>%
  filter(!is.na(medullary_fsoc), medullary_fsoc >= 0) %>%
  mutate(
    fsoc_binary = factor(
      ifelse(medullary_fsoc < median(medullary_fsoc, na.rm = TRUE), "Impaired", "Normal"),
      levels = c("Normal", "Impaired")
    )
  )

# ============================================================================
# 3. IDENTIFY IMMUNE-RELATED PROTEINS
# ============================================================================

cat("=== IDENTIFYING IMMUNE-RELATED PROTEINS ===\n")

# COMPREHENSIVE immune/inflammatory keywords - EXPANDED
immune_keywords <- c(
  # Interleukins (all)
  "interleukin", "IL-", "IL1", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL8", 
  "IL9", "IL10", "IL11", "IL12", "IL13", "IL14", "IL15", "IL16", "IL17", "IL18",
  "IL19", "IL20", "IL21", "IL22", "IL23", "IL24", "IL25", "IL26", "IL27", "IL28",
  "IL29", "IL30", "IL31", "IL32", "IL33", "IL34", "IL35", "IL36", "IL37",
  
  # TNF superfamily
  "TNF", "tumor necrosis", "TNFR", "TNFSF", "TRAIL", "TWEAK", "lymphotoxin",
  "FAS", "FASL", "CD95", "OX40", "4-1BB", "GITR", "LIGHT", "APRIL", "BAFF",
  
  # Chemokines (comprehensive)
  "chemokine", "CCL", "CXCL", "CX3CL", "XCL", 
  "MCP", "MIP", "RANTES", "eotaxin", "IP-10", "MIG", "fractalkine", "SDF",
  "GRO", "ENA", "NAP", "I-TAC", "Mig", "lymphotactin",
  
  # Inflammatory cytokines & markers
  "C-reactive", "CRP", "serum amyloid", "SAA", "pentraxin", "PTX",
  "calprotectin", "S100A", "HMGB1", "ferritin", "procalcitonin",
  
  # Complement system
  "complement", " C1", " C2", " C3", " C4", " C5", " C6", " C7", " C8", " C9",
  "factor B", "factor D", "factor H", "factor I", "properdin", "mannose-binding",
  "MBL", "ficolin", "C1q", "C3a", "C3b", "C5a", "MAC",
  
  # Adhesion molecules
  "adhesion", "ICAM", "VCAM", "selectin", "SELE", "SELL", "SELP",
  "integrin", "PECAM", "CD31", "PSGL-1", "LFA-1", "VLA-4", "cadherin",
  
  # Acute phase proteins
  "haptoglobin", "fibrinogen", "alpha-1-antitrypsin", "ceruloplasmin",
  "alpha-2-macroglobulin", "transferrin", "albumin",
  
  # Growth factors & angiogenesis
  "growth factor", "VEGF", "PDGF", "FGF", "TGF", "EGF", "HGF", "IGF",
  "NGF", "BDNF", "GDNF", "angiopoietin", "angiogenin", "endoglin",
  "PlGF", "bFGF", "SCF", "thrombopoietin",
  
  # Matrix & proteases
  "matrix metalloproteinase", "MMP", "TIMP", "ADAM", "ADAMTS",
  "cathepsin", "elastase", "collagenase", "gelatinase", "stromelysin",
  "kallikrein", "urokinase", "plasminogen", "tissue factor",
  
  # Cell surface markers & receptors
  "CD40", "CD14", "CD163", "CD80", "CD86", "CD28", "CTLA", "PD-1", "PD-L1", "PD-L2",
  "CD3", "CD4", "CD8", "CD11", "CD16", "CD19", "CD20", "CD25", "CD56", "CD68",
  "CD138", "HLA-DR", "Toll-like", "TLR", "NOD", "RAGE", "scavenger receptor",
  
  # Interferons
  "interferon", "IFN", "IFN-gamma", "IFN-alpha", "IFN-beta", "IFN-lambda",
  
  # Colony stimulating & hematopoietic factors
  "CSF", "G-CSF", "GM-CSF", "M-CSF", "erythropoietin", "EPO", "thrombopoietin",
  
  # Immune cell markers
  "macrophage", "monocyte", "neutrophil", "lymphocyte", "eosinophil", 
  "basophil", "dendritic", "NK cell", "mast cell", "T cell", "B cell",
  
  # Adipokines & metabolic inflammation
  "osteopontin", "resistin", "adiponectin", "leptin", "visfatin", "omentin",
  "apelin", "chemerin", "vaspin", "RBP4", "fetuin", "PAI-1",
  
  # Kidney injury & fibrosis
  "kidney injury", "KIM-1", "NGAL", "lipocalin", "cystatin", "clusterin", 
  "uromodulin", "nephrin", "podocin", "synaptopodin", "CTGF", "periostin",
  
  # Oxidative stress & ROS
  "myeloperoxidase", "MPO", "NADPH oxidase", "NOX", "superoxide dismutase",
  "catalase", "glutathione", "thioredoxin", "peroxiredoxin",
  
  # Coagulation & hemostasis
  "von Willebrand", "vWF", "thrombin", "fibrin", "D-dimer", "plasmin",
  "protein C", "protein S", "antithrombin", "tissue factor",
  
  # Apoptosis & cell death
  "caspase", "BCL-2", "BAX", "cytochrome c", "PARP", "annexin",
  "survivin", "XIAP", "Smac", "DIABLO",
  
  # Signaling molecules
  "JAK", "STAT", "MAPK", "NF-kB", "PI3K", "AKT", "mTOR", "AMPK",
  
  # Specialized immune
  "perforin", "granzyme", "defensin", "cathelicidin", "lactoferrin",
  "lysozyme", "mucin", "surfactant protein", "collectin",
  
  # Autoimmune & rheumatologic
  "rheumatoid factor", "anti-CCP", "ANA", "dsDNA", "Sm antigen",
  "Jo-1", "Scl-70", "centromere", "histone",
  
  # Other inflammation-related
  "prostaglandin", "leukotriene", "lipoxin", "thromboxane",
  "bradykinin", "histamine", "serotonin", "substance P",
  "neopterin", "tryptase", "chymase", "nitric oxide", "iNOS"
)

# Find immune proteins
immune_results <- results %>%
  filter(!is.na(protein_name)) %>%
  filter(str_detect(tolower(protein_name), 
                    paste(tolower(immune_keywords), collapse = "|")))

cat("Immune-related proteins found:", nrow(immune_results), "\n")
cat("Significant (p < 0.05):", sum(immune_results$p_value < 0.05), "\n")
cat("Significant (FDR < 0.10):", sum(immune_results$p_fdr < 0.10), "\n\n")

# Save immune results
write.csv(immune_results, 
          file.path(OUTPUT_DIR, "immune_proteins_results.csv"),
          row.names = FALSE)

# ============================================================================
# 4. VOLCANO PLOT - IMMUNE PROTEINS ONLY
# ============================================================================

cat("=== CREATING IMMUNE-FOCUSED VOLCANO PLOT ===\n")

sig_immune <- immune_results %>% filter(p_fdr < 0.10)
top_immune <- immune_results %>% filter(p_value < 0.05)

# Extract gene symbols from protein names (first word or in parentheses)
immune_results <- immune_results %>%
  mutate(
    gene_symbol = case_when(
      str_detect(protein_name, "\\(([A-Z0-9-]+)\\)") ~ 
        str_extract(protein_name, "(?<=\\()[A-Z0-9-]+(?=\\))"),
      TRUE ~ word(protein_name, 1)
    )
  )

p_volcano_immune <- ggplot(immune_results, 
                           aes(x = log2_fc, y = -log10(p_value))) +
  geom_point(aes(color = p_fdr < 0.10, size = -log10(p_value)), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "blue", alpha = 0.5, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 0.5) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", 
             color = "gray50", alpha = 0.3) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#E41A1C"),
                     name = "FDR < 0.10",
                     labels = c("No", "Yes")) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  theme_bw(base_size = 12) +
  labs(x = "Log2 Fold Change (Impaired / Normal FSOC)",
       y = "-log10(P-value)",
       title = "Immune-Related Proteins: Impaired vs Normal Medullary FSOC",
       subtitle = paste0(nrow(immune_results), " immune proteins tested, ",
                         sum(immune_results$p_value < 0.05), " nominal p<0.05, ",
                         sum(immune_results$p_fdr < 0.10), " FDR<0.10")) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

# Label significant proteins
if (nrow(sig_immune) > 0) {
  # Use proteins with FDR < 0.10
  label_data <- sig_immune %>%
    mutate(
      gene_symbol = case_when(
        str_detect(protein_name, "\\(([A-Z0-9-]+)\\)") ~ 
          str_extract(protein_name, "(?<=\\()[A-Z0-9-]+(?=\\))"),
        TRUE ~ word(protein_name, 1)
      )
    )
  
  p_volcano_immune <- p_volcano_immune +
    geom_text_repel(data = label_data,
                    aes(label = gene_symbol), 
                    size = 3, fontface = "bold",
                    max.overlaps = 30, 
                    segment.size = 0.3,
                    segment.color = "gray40",
                    min.segment.length = 0,
                    box.padding = 0.5)
} else if (nrow(top_immune) > 0) {
  # If no FDR<0.10, label top 20 by p-value
  label_data <- head(top_immune, 20) %>%
    mutate(
      gene_symbol = case_when(
        str_detect(protein_name, "\\(([A-Z0-9-]+)\\)") ~ 
          str_extract(protein_name, "(?<=\\()[A-Z0-9-]+(?=\\))"),
        TRUE ~ word(protein_name, 1)
      )
    )
  
  p_volcano_immune <- p_volcano_immune +
    geom_text_repel(data = label_data,
                    aes(label = gene_symbol), 
                    size = 2.5,
                    max.overlaps = 20, 
                    segment.size = 0.2)
}

ggsave(file.path(OUTPUT_DIR, "volcano_immune_proteins.pdf"), 
       p_volcano_immune, width = 12, height = 10)
ggsave(file.path(OUTPUT_DIR, "volcano_immune_proteins.png"), 
       p_volcano_immune, width = 12, height = 10, dpi = 300)

print(p_volcano_immune)

# ============================================================================
# 5. ENHANCED BOXPLOTS WITH REAL PROTEIN NAMES
# ============================================================================

cat("\n=== CREATING ENHANCED BOXPLOTS ===\n")

# Get top proteins (prioritize FDR, then p-value)
if (sum(immune_results$p_fdr < 0.10) >= 12) {
  top_proteins <- immune_results %>% 
    filter(p_fdr < 0.10) %>%
    arrange(p_value) %>%
    head(12)
} else {
  top_proteins <- immune_results %>% 
    arrange(p_value) %>%
    head(12)
}

plot_list <- list()

for(i in 1:nrow(top_proteins)) {
  var <- top_proteins$variable[i]
  
  plot_data <- dat_fsoc_proteomics %>%
    select(fsoc_binary, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value), !is.na(fsoc_binary))
  
  # Extract clean protein name
  full_name <- top_proteins$protein_name[i]
  
  # Try to extract gene symbol in parentheses
  if (str_detect(full_name, "\\([A-Z0-9-]+\\)")) {
    gene_symbol <- str_extract(full_name, "(?<=\\()[A-Z0-9-]+(?=\\))")
    # Get protein name before parentheses
    protein_part <- str_extract(full_name, "^[^\\(]+") %>% str_trim()
    # Truncate if too long
    if (nchar(protein_part) > 30) {
      protein_part <- paste0(substr(protein_part, 1, 27), "...")
    }
    display_name <- paste0(gene_symbol, ": ", protein_part)
  } else {
    # Just use first 40 characters
    display_name <- if(nchar(full_name) > 40) {
      paste0(substr(full_name, 1, 37), "...")
    } else {
      full_name
    }
  }
  
  # Statistical annotation
  p_val <- top_proteins$p_value[i]
  p_fdr <- top_proteins$p_fdr[i]
  fc <- top_proteins$fold_change[i]
  
  p_label <- sprintf("p=%.2e, FDR=%.3f, FC=%.2f", p_val, p_fdr, fc)
  
  # Add significance stars
  stars <- if(p_fdr < 0.001) "***" else if(p_fdr < 0.01) "**" else if(p_fdr < 0.05) "*" else ""
  
  p <- ggplot(plot_data, aes(x = fsoc_binary, y = value, fill = fsoc_binary)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.6) +
    geom_jitter(width = 0.25, alpha = 0.6, size = 2) +
    scale_fill_manual(values = c("Normal" = "#3182BD", "Impaired" = "#E6550D")) +
    theme_bw(base_size = 10) +
    labs(title = display_name,
         subtitle = paste0(p_label, " ", stars),
         x = "FSOC Status",
         y = "Expression (RFU)") +
    theme(legend.position = "none",
          plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7.5),
          panel.grid.major.x = element_blank())
  
  # Add p-value annotation on plot
  y_max <- max(plot_data$value, na.rm = TRUE)
  y_min <- min(plot_data$value, na.rm = TRUE)
  y_range <- y_max - y_min
  
  if (p_val < 0.05) {
    p <- p + annotate("text", x = 1.5, y = y_max + 0.05 * y_range,
                      label = stars, size = 6, fontface = "bold")
  }
  
  plot_list[[i]] <- p
}

png(file.path(OUTPUT_DIR, 'top_immune_proteins_boxplots.png'), 
    height = 15, width = 16, units = 'in', res = 300)
grid.arrange(grobs = plot_list, ncol = 3,
             top = grid::textGrob("Top Immune Proteins by FSOC Status", 
                                  gp = grid::gpar(fontsize = 16, fontface = "bold")))
dev.off()

cat("Enhanced boxplots saved.\n")

# ============================================================================
# 6. PATHWAY ENRICHMENT ANALYSIS
# ============================================================================

cat("\n=== PATHWAY ENRICHMENT ANALYSIS ===\n")

# Extract gene symbols for enrichment
sig_genes <- immune_results %>%
  filter(p_value < 0.05) %>%
  mutate(
    gene_symbol = case_when(
      str_detect(protein_name, "\\(([A-Z0-9-]+)\\)") ~ 
        str_extract(protein_name, "(?<=\\()[A-Z0-9-]+(?=\\))"),
      TRUE ~ word(protein_name, 1)
    )
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

cat("Genes for enrichment:", length(sig_genes), "\n")

if (length(sig_genes) >= 5) {
  # Set enrichR databases
  dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", 
           "Reactome_2022", "WikiPathways_2019_Human")
  
  tryCatch({
    # Run enrichment
    enriched <- enrichr(sig_genes, dbs)
    
    # Process results for each database
    for (db_name in names(enriched)) {
      result_df <- enriched[[db_name]]
      
      if (nrow(result_df) > 0) {
        # Filter significant
        sig_pathways <- result_df %>%
          filter(Adjusted.P.value < 0.05) %>%
          arrange(Adjusted.P.value) %>%
          head(20)
        
        if (nrow(sig_pathways) > 0) {
          cat("\n=== ", db_name, " ===\n", sep = "")
          print(sig_pathways %>% select(Term, Overlap, Adjusted.P.value))
          
          # Save
          write.csv(sig_pathways, 
                    file.path(OUTPUT_DIR, paste0("enrichment_", db_name, ".csv")),
                    row.names = FALSE)
        }
      }
    }
    
    # Create enrichment plot for top pathways
    if (nrow(enriched[[1]]) > 0) {
      top_pathways <- enriched[[1]] %>%
        filter(Adjusted.P.value < 0.05) %>%
        arrange(Adjusted.P.value) %>%
        head(15) %>%
        mutate(
          Term = str_trunc(Term, 60),
          logP = -log10(Adjusted.P.value)
        )
      
      if (nrow(top_pathways) > 0) {
        p_enrich <- ggplot(top_pathways, 
                           aes(x = reorder(Term, logP), y = logP)) +
          geom_col(fill = "#3182BD", alpha = 0.8) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                     color = "red") +
          coord_flip() +
          theme_bw(base_size = 11) +
          labs(x = "", y = "-log10(Adjusted P-value)",
               title = "Pathway Enrichment: Immune Proteins (p<0.05)",
               subtitle = paste0(length(sig_genes), " genes analyzed")) +
          theme(panel.grid.major.y = element_blank())
        
        ggsave(file.path(OUTPUT_DIR, "pathway_enrichment.pdf"), 
               p_enrich, width = 10, height = 8)
        ggsave(file.path(OUTPUT_DIR, "pathway_enrichment.png"), 
               p_enrich, width = 10, height = 8, dpi = 300)
        
        cat("\nEnrichment plot saved.\n")
      }
    }
    
  }, error = function(e) {
    cat("Enrichment analysis failed:", e$message, "\n")
    cat("This is optional - main analysis is complete.\n")
  })
} else {
  cat("Not enough genes (n=", length(sig_genes), ") for enrichment analysis\n")
}

# ============================================================================
# 7. SUMMARY TABLE
# ============================================================================

cat("\n=== CREATING SUMMARY TABLE ===\n")

summary_table <- immune_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  mutate(
    gene_symbol = case_when(
      str_detect(protein_name, "\\(([A-Z0-9-]+)\\)") ~ 
        str_extract(protein_name, "(?<=\\()[A-Z0-9-]+(?=\\))"),
      TRUE ~ word(protein_name, 1)
    ),
    Direction = ifelse(log2_fc > 0, "↑ Impaired", "↓ Impaired"),
    `P-value` = formatC(p_value, format = "e", digits = 2),
    `FDR` = formatC(p_fdr, format = "f", digits = 3),
    `Fold Change` = formatC(fold_change, format = "f", digits = 2)
  ) %>%
  select(gene_symbol, protein_name, Direction, `P-value`, FDR, `Fold Change`) %>%
  head(30)

write.csv(summary_table, 
          file.path(OUTPUT_DIR, "immune_proteins_summary_table.csv"),
          row.names = FALSE)

cat("\n=== FINAL SUMMARY ===\n")
cat("Immune proteins analyzed:", nrow(immune_results), "\n")
cat("Nominally significant (p<0.05):", sum(immune_results$p_value < 0.05), "\n")
cat("FDR significant (q<0.10):", sum(immune_results$p_fdr < 0.10), "\n")
cat("Upregulated in Impaired:", sum(immune_results$log2_fc > 0 & immune_results$p_value < 0.05), "\n")
cat("Downregulated in Impaired:", sum(immune_results$log2_fc < 0 & immune_results$p_value < 0.05), "\n")

cat("\n\n=== FILES CREATED ===\n")
cat("1. immune_proteins_results.csv - Full immune protein results\n")
cat("2. volcano_immune_proteins.pdf/.png - Immune-focused volcano plot\n")
cat("3. top_immune_proteins_boxplots.png - Boxplots with real names\n")
cat("4. immune_proteins_summary_table.csv - Summary table\n")
if (length(sig_genes) >= 5) {
  cat("5. enrichment_*.csv - Pathway enrichment results\n")
  cat("6. pathway_enrichment.pdf/.png - Enrichment visualization\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")



