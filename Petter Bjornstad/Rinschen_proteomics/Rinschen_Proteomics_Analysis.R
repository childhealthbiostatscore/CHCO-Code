##############################################################################
# Proteomics + Phosphoproteomics Pipeline
# Human Kidney Biopsies | T2D ± SGLT2 inhibitor
# Data: Rinschen lab (woB1 = batch 1 excluded, already log2-transformed)
# Files: hKidneyBiopsies_PROT_Expression_woB1.xlsx
#        hKidneyBiopsies_PHOS_Expression_woB1.xlsx
#        Samples Shipped.xlsx
##############################################################################

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────
pkgs <- c(
  "tidyverse", "readxl", "limma", "vsn",
  "pheatmap", "ggrepel", "patchwork",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "corrplot", "broom", "janitor", "openxlsx"
)
# BiocManager::install(setdiff(pkgs, rownames(installed.packages())))
invisible(lapply(pkgs, library, character.only = TRUE))

# ── Resolve common Bioconductor masking conflicts ─────────────────────────────
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::slice_min)
conflicts_prefer(base::intersect)

# ── 1. PATHS ─────────────────────────────────────────────────────────────────
base_dir <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Rinschen/"
dir.create(file.path(base_dir, "results"), showWarnings = FALSE)
dir.create(file.path(base_dir, "plots"),   showWarnings = FALSE)

# ── 2. LOAD RAW FILES ────────────────────────────────────────────────────────
prot_raw <- read_xlsx(file.path(base_dir, "hKidneyBiopsies_PROT_Expression_woB1.xlsx"))
phos_raw <- read_xlsx(file.path(base_dir, "hKidneyBiopsies_PHOS_Expression_woB1.xlsx"))
meta_raw <- read_xlsx(file.path(base_dir, "Samples Shipped.xlsx")) %>% clean_names()

##############################################################################
# SECTION 1: METADATA WRANGLING
##############################################################################

meta <- meta_raw %>%
  # Drop the one sample (CRC-04) with near-complete NA metadata
  filter(!is.na(treatment)) %>%
  mutate(
    # Strip dashes so S-1907-005614 → S1907005614 for column matching
    aliquot_nodash = str_remove_all(rn_alater_id, "-"),
    # Cohort from study_id prefix
    cohort = str_extract(study_id, "^[A-Z]+"),
    # Binary treatment group
    group  = if_else(treatment == "SGLT2 inhibitor", "SGLT2i", "Control"),
    sex    = factor(sex),
    age    = as.numeric(age)
  )

cat("Samples after QC:", nrow(meta), "\n")
cat("Group breakdown:\n"); print(table(meta$group, meta$cohort))

# ── Helper: extract aliquot ID from long DIA-NN column names ─────────────────
# Column pattern: "...hKidneyBiopsy_S1907005614_F_B2.raw.PG.Quantity"
extract_aliquot <- function(col_names) {
  str_extract(col_names, "(?<=hKidneyBiopsy_)S\\d+")
}

##############################################################################
# SECTION 2: BUILD ABUNDANCE MATRICES
##############################################################################

build_matrix <- function(df, id_col) {
  # Identify sample columns and their aliquot IDs
  quant_cols  <- names(df)[str_detect(names(df), "Quantity")]
  aliquot_ids <- extract_aliquot(quant_cols)
  
  # Keep only columns matching our metadata
  keep        <- aliquot_ids %in% meta$aliquot_nodash
  quant_cols  <- quant_cols[keep]
  aliquot_ids <- aliquot_ids[keep]
  
  # Subset, drop NA/duplicate IDs, convert to matrix
  mat_df <- df[!is.na(df[[id_col]]), c(id_col, quant_cols)]
  mat_df <- mat_df[!duplicated(mat_df[[id_col]]), ]
  
  mat             <- as.matrix(mat_df[, quant_cols])
  rownames(mat)   <- mat_df[[id_col]]
  
  # Rename columns: aliquot ID → study_id
  col_to_study    <- setNames(meta$study_id, meta$aliquot_nodash)
  colnames(mat)   <- col_to_study[aliquot_ids]
  
  # Only keep columns that matched — some study_ids may not be in matrix
  matched <- base::intersect(meta$study_id, colnames(mat))
  cat(id_col, "— matched", length(matched), "of", nrow(meta), "samples\n")
  mat[, matched, drop = FALSE]
}

prot_mat <- build_matrix(prot_raw, "GENE")    # rows = gene symbols
phos_mat <- build_matrix(phos_raw, "PSITE")   # rows = phosphosite IDs (e.g. TNS2_S762)

# Lookup table: phosphosite → gene symbol
psite_gene <- phos_raw %>%
  select(PSITE, GENE) %>%
  distinct() %>%
  rename(psite_id = PSITE, gene_symbol = GENE)

cat("Proteomics matrix:      ", nrow(prot_mat), "proteins ×", ncol(prot_mat), "samples\n")
cat("Phosphoproteomics matrix:", nrow(phos_mat), "sites ×",   ncol(phos_mat), "samples\n")

# Note: data are already log2-transformed by DIA-NN — do NOT log-transform again

##############################################################################
# SECTION 3: QC AND NORMALIZATION
##############################################################################

run_qc_norm <- function(mat, meta_df, label = "PROT") {
  
  # -- Missingness --
  miss_pct <- apply(mat, 1, function(x) mean(is.na(x)) * 100)
  
  p_miss <- ggplot(data.frame(pct = miss_pct), aes(pct)) +
    geom_histogram(bins = 40, fill = "#4E79A7") +
    labs(title = paste(label, "— Feature missingness"),
         x = "% Missing", y = "Count") +
    theme_bw()
  ggsave(file.path(base_dir, paste0("plots/", label, "_missingness.pdf")),
         p_miss, width = 6, height = 4)
  
  # Keep features present in ≥ 70% of samples
  mat_filt <- mat[miss_pct <= 30, ]
  cat(label, "— features after missingness filter:", nrow(mat_filt), "\n")
  
  # -- Median centering (data already log2; VSN not needed) --
  col_meds  <- apply(mat_filt, 2, median, na.rm = TRUE)
  grand_med <- median(col_meds)
  mat_norm  <- sweep(mat_filt, 2, col_meds - grand_med, "-")
  
  # Distribution boxplot
  box_df <- mat_norm %>% as.data.frame() %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "study_id", values_to = "intensity") %>%
    left_join(meta_df %>% select(study_id, group, cohort), by = "study_id")
  
  p_box <- ggplot(box_df, aes(study_id, intensity, fill = group)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual(values = c(SGLT2i = "#E15759", Control = "#4E79A7")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    labs(title = paste(label, "— Normalized intensities"), x = NULL)
  ggsave(file.path(base_dir, paste0("plots/", label, "_boxplot.pdf")),
         p_box, width = 10, height = 4)
  
  # -- PCA --
  mat_pca <- mat_norm
  row_meds <- apply(mat_pca, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(mat_pca)))
    mat_pca[i, is.na(mat_pca[i, ])] <- row_meds[i]
  
  pca     <- prcomp(t(mat_pca), scale. = TRUE)
  var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  pca_df  <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("study_id") %>%
    left_join(meta_df, by = "study_id")
  
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = cohort,
                              label = study_id)) +
    geom_point(size = 3) +
    geom_text_repel(size = 2.5, max.overlaps = 15) +
    scale_color_manual(values = c(SGLT2i = "#E15759", Control = "#4E79A7")) +
    labs(title = paste(label, "— PCA"),
         x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw()
  ggsave(file.path(base_dir, paste0("plots/", label, "_PCA.pdf")),
         p_pca, width = 7, height = 5)
  
  # -- MinProb imputation --
  mat_imp <- mat_norm
  for (j in seq_len(ncol(mat_imp))) {
    col <- mat_imp[, j]; idx <- is.na(col)
    if (any(idx)) {
      lv <- quantile(col, 0.01, na.rm = TRUE)
      mat_imp[idx, j] <- rnorm(sum(idx), lv, sd(col, na.rm = TRUE) * 0.3)
    }
  }
  
  mat_imp
}

prot_imp <- run_qc_norm(prot_mat, meta, "PROT")
phos_imp <- run_qc_norm(phos_mat, meta, "PHOS")

##############################################################################
# SECTION 4: DIFFERENTIAL ABUNDANCE — SGLT2i vs. Control
##############################################################################

run_limma <- function(mat, meta_df, label = "PROT") {
  
  md <- meta_df %>%
    filter(study_id %in% colnames(mat)) %>%
    arrange(match(study_id, colnames(mat))) %>%
    mutate(group = factor(group, levels = c("Control", "SGLT2i")),
           sex   = factor(sex))
  
  # Covariates: age and sex (cohort subsumed by small N; add if model converges)
  design  <- model.matrix(~ 0 + group + age + sex, data = md)
  colnames(design) <- make.names(colnames(design))
  
  cont_mat <- makeContrasts(
    SGLT2i_vs_Control = groupSGLT2i - groupControl,
    levels = design
  )
  
  fit  <- lmFit(mat[, md$study_id], design)
  fit2 <- contrasts.fit(fit, cont_mat)
  fit2 <- eBayes(fit2)
  
  res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P") %>%
    rownames_to_column("feature_id") %>%
    mutate(data_type = label)
  
  n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  cat(label, "— significant features (FDR < 0.05):", n_sig, "\n")
  res
}

res_prot <- run_limma(prot_imp, meta, "PROT")
res_phos <- run_limma(phos_imp, meta, "PHOS") %>%
  left_join(psite_gene, by = c("feature_id" = "psite_id"))

write.xlsx(list(Proteomics        = res_prot,
                Phosphoproteomics = res_phos),
           file.path(base_dir, "results/differential_abundance_PROT_PHOS.xlsx"))

# ── Volcano plots ─────────────────────────────────────────────────────────────
plot_volcano <- function(res, title, fc_cut = 0.5, fdr_cut = 0.05) {
  res %>%
    mutate(
      dir   = case_when(
        adj.P.Val < fdr_cut & logFC >  fc_cut ~ "Up (SGLT2i)",
        adj.P.Val < fdr_cut & logFC < -fc_cut ~ "Down (SGLT2i)",
        TRUE ~ "NS"
      ),
      label = if_else(adj.P.Val < fdr_cut & abs(logFC) > fc_cut,
                      feature_id, NA_character_)
    ) %>%
    ggplot(aes(logFC, -log10(P.Value), color = dir)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20) +
    scale_color_manual(values = c("Up (SGLT2i)"   = "#E15759",
                                  "Down (SGLT2i)" = "#4E79A7",
                                  NS             = "grey70")) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05),        linetype = "dashed", color = "grey50") +
    labs(title = title, x = "log2 FC (SGLT2i vs Control)",
         y = "-log10(p-value)", color = NULL) +
    theme_bw()
}

p_volc <- plot_volcano(res_prot, "Proteomics — SGLT2i vs Control") /
  plot_volcano(res_phos, "Phosphoproteomics — SGLT2i vs Control")
ggsave(file.path(base_dir, "plots/volcano_PROT_PHOS.pdf"), p_volc, width = 9, height = 12)

##############################################################################
# SECTION 5: PHOSPHO / PROTEIN NORMALIZATION
# Separates phospho changes driven by total protein abundance from
# those reflecting true changes in site-specific phosphorylation
##############################################################################

phos_norm <- res_phos %>%
  left_join(res_prot %>% select(feature_id, logFC) %>%
              rename(gene_symbol = feature_id, prot_logFC = logFC),
            by = "gene_symbol") %>%
  mutate(
    logFC_adjusted = logFC - replace_na(prot_logFC, 0),
    site_type = case_when(
      is.na(prot_logFC)                                  ~ "No protein data",
      abs(prot_logFC) > 0.5 &
        sign(logFC) == sign(prot_logFC)                  ~ "Protein-driven",
      TRUE                                               ~ "Phospho-specific"
    )
  )

write.csv(phos_norm,
          file.path(base_dir, "results/phospho_normalized_to_protein.csv"),
          row.names = FALSE)

ggplot(phos_norm %>% filter(!is.na(prot_logFC)),
       aes(prot_logFC, logFC, color = site_type)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Protein-driven"   = "#E15759",
                                "Phospho-specific" = "#4E79A7",
                                "No protein data"  = "grey70")) +
  labs(title = "Phospho FC vs. Protein FC (SGLT2i vs Control)",
       x = "log2FC (Protein)", y = "log2FC (Phosphosite)", color = NULL) +
  theme_bw()
ggsave(file.path(base_dir, "plots/phospho_vs_protein_FC.pdf"), width = 7, height = 6)

##############################################################################
# SECTION 6: PATHWAY ENRICHMENT — PROT and PHOS
##############################################################################

run_go_kegg <- function(sig_genes, universe_genes, label = "") {
  
  sig_genes  <- intersect(sig_genes, universe_genes)
  if (length(sig_genes) < 5) {
    cat("Too few significant genes for enrichment:", label, "\n"); return(NULL)
  }
  
  ego <- tryCatch(
    enrichGO(gene = sig_genes, universe = universe_genes,
             OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) { message("GO failed: ", e$message); NULL }
  )
  
  if (!is.null(ego) && nrow(ego) > 0) {
    p <- dotplot(ego, showCategory = 20,
                 title = paste("GO BP —", label, "(SGLT2i vs Control)"))
    ggsave(file.path(base_dir, paste0("plots/GO_", label, ".pdf")),
           p, width = 9, height = 8)
  }
  
  # KEGG
  ez      <- bitr(sig_genes,      "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
  univ_ez <- bitr(universe_genes, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
  ekegg   <- tryCatch(
    enrichKEGG(gene = ez$ENTREZID, universe = univ_ez$ENTREZID,
               organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05),
    error = function(e) { message("KEGG failed: ", e$message); NULL }
  )
  
  list(go = ego, kegg = ekegg)
}

# With n=12 samples, FDR < 0.05 is very stringent — use nominal p < 0.05
# PROT
sig_up_prot   <- res_prot %>% filter(P.Value < 0.05, logFC >  0.5) %>% pull(feature_id)
sig_down_prot <- res_prot %>% filter(P.Value < 0.05, logFC < -0.5) %>% pull(feature_id)
run_go_kegg(sig_up_prot,   rownames(prot_imp), "PROT_up")
run_go_kegg(sig_down_prot, rownames(prot_imp), "PROT_down")

# PHOS (use gene symbols of significant phosphosites)
sig_up_phos   <- phos_norm %>%
  filter(P.Value < 0.05, logFC_adjusted >  0.5) %>%
  pull(gene_symbol) %>% unique()
sig_down_phos <- phos_norm %>%
  filter(P.Value < 0.05, logFC_adjusted < -0.5) %>%
  pull(gene_symbol) %>% unique()
run_go_kegg(sig_up_phos,   unique(psite_gene$gene_symbol), "PHOS_up")
run_go_kegg(sig_down_phos, unique(psite_gene$gene_symbol), "PHOS_down")

##############################################################################
# SECTION 7: CROSS-COHORT CONSISTENCY CHECK
# With only 12 samples, formal per-cohort limma is underpowered,
# but we can check whether effect directions replicate across cohorts
##############################################################################

# Get top 30 PROT hits by nominal p-value (FDR too stringent at n=12)
top_prot <- res_prot %>%
  slice_min(P.Value, n = 30) %>%
  pull(feature_id)

cohort_fc <- map_dfr(unique(meta$cohort), function(coh) {
  md_c <- meta %>% filter(cohort == coh, study_id %in% colnames(prot_imp))
  grp_counts <- table(md_c$group)
  
  # Need both groups present AND at least 2 samples per group for eBayes
  if (length(grp_counts) < 2 || any(grp_counts < 2)) {
    cat("Skipping cohort", coh, "— insufficient group sizes:", 
        paste(names(grp_counts), grp_counts, sep = "=", collapse = ", "), "\n")
    return(NULL)
  }
  mat_c  <- prot_imp[top_prot, md_c$study_id, drop = FALSE]
  design <- model.matrix(~ group, data = md_c %>%
                           mutate(group = factor(group, levels = c("Control","SGLT2i"))))
  fit    <- eBayes(lmFit(mat_c, design))
  topTable(fit, coef = 2, number = Inf) %>%
    rownames_to_column("feature_id") %>%
    mutate(cohort = coh)
}) %>%
  bind_rows()

if (nrow(cohort_fc) > 0) {
  fc_wide <- cohort_fc %>%
    select(feature_id, cohort, logFC) %>%
    pivot_wider(names_from = cohort, values_from = logFC)
  
  forest_df <- cohort_fc %>%
    mutate(lower = logFC - 1.96 * sqrt(1 / AveExpr),
           upper = logFC + 1.96 * sqrt(1 / AveExpr))
  
  ggplot(forest_df, aes(logFC, feature_id, color = cohort)) +
    geom_point(position = position_dodge(0.5), size = 2.5) +
    geom_linerange(aes(xmin = lower, xmax = upper),
                   position = position_dodge(0.5)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(RH = "#4E79A7", CRC = "#E15759", IT = "#76B7B2")) +
    labs(title = "Effect sizes by cohort — top PROT hits",
         x = "log2 FC (SGLT2i vs Control)", y = NULL) +
    theme_bw()
  ggsave(file.path(base_dir, "plots/forest_by_cohort.pdf"), width = 9, height = 8)
}

##############################################################################
# SECTION 8: scRNA-SEQ INTEGRATION (uncomment when data available)
##############################################################################

# scrna_hits <- read_csv(file.path(base_dir, "scrna_top_hits.csv"))
# Expected cols: gene_symbol, cell_type, log2FC_rna, fdr_rna

# concordance <- res_prot %>%
#   filter(feature_id %in% scrna_hits$gene_symbol) %>%
#   left_join(scrna_hits, by = c("feature_id" = "gene_symbol")) %>%
#   mutate(
#     concordant = sign(logFC) == sign(log2FC_rna),
#     sig_both   = adj.P.Val < 0.05 & fdr_rna < 0.05
#   )
#
# ggplot(concordance, aes(log2FC_rna, logFC, color = sig_both)) +
#   geom_point() +
#   geom_text_repel(data = filter(concordance, sig_both),
#                   aes(label = feature_id), size = 2.5) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_color_manual(values = c("TRUE" = "#E15759", "FALSE" = "grey60")) +
#   labs(title = "mRNA–protein concordance (SGLT2i effect)",
#        x = "log2FC (scRNA-seq)", y = "log2FC (Proteomics)") +
#   theme_bw()

##############################################################################
# SECTION 9: SUMMARY TABLE
##############################################################################

summary_table <- res_prot %>%
  filter(adj.P.Val < 0.05) %>%
  left_join(
    phos_norm %>%
      filter(adj.P.Val < 0.05) %>%
      count(gene_symbol, name = "n_sig_phosphosites"),
    by = c("feature_id" = "gene_symbol")
  ) %>%
  arrange(P.Value) %>%
  select(feature_id, logFC, AveExpr, t, P.Value, adj.P.Val,
         n_sig_phosphosites, data_type)

write.xlsx(list(Summary = summary_table),
           file.path(base_dir, "results/summary_table_ASN_abstract.xlsx"))

cat("\n=== Pipeline complete ===\n",
    "results/differential_abundance_PROT_PHOS.xlsx\n",
    "results/phospho_normalized_to_protein.csv\n",
    "results/summary_table_ASN_abstract.xlsx\n",
    "plots/: QC, volcano, pathway, forest plots\n")








##############################################################################
# SECTION 10: PATHO PLEX COMPARISON
##############################################################################

library(ggnewscale)

# ── Paths to PathoPlex supplementary tables ───────────────────────────────────
base_patho <- "C:/Users/netio/Downloads/41586_2025_9225_MOESM3_ESM/2023-11-20558D-s3/"
s3_path    <- file.path(base_patho, "R3_Table_S3.xlsx")
s5_path    <- file.path(base_patho, "R3_Table_S5.xlsx")

# ── Load and clean S3 antibody panel ─────────────────────────────────────────
s3_raw <- read_excel(s3_path, skip = 1)
colnames(s3_raw) <- c("protein_label", "target")
s3 <- s3_raw %>%
  filter(!is.na(protein_label), protein_label != "Protein") %>%
  mutate(
    gene_symbol = case_when(
      str_detect(protein_label, "\\(") ~
        str_extract(protein_label, "(?<=\\()([^)]+)(?=\\))"),
      TRUE ~ protein_label
    ),
    gene_symbol = str_trim(str_to_upper(gene_symbol))
  )

# ── Alias map: protein labels → gene symbols ──────────────────────────────────
alias_map <- c(
  "Β-CATENIN"                = "CTNNB1",
  "C-FOS"                    = "FOS",
  "CALPAIN SMALL SUBUNIT  1" = "CAPNS1",
  "CALPASTATIN"              = "CAST",
  "CALRETICULIN"             = "CALR",
  "CD3"                      = "CD247",
  "CD41"                     = "ITGA2B",
  "CD42B"                    = "GP1BA",
  "CLAUDIN-1"                = "CLDN1",
  "COLLAGEN IV"              = "COL4A1",
  "COLLAGEN TYPE III"        = "COL3A1",
  "COLLAGEN V"               = "COL5A1",
  "CYCLIN B1"                = "CCNB1",
  "CYTOKERATIN 19"           = "KRT19",
  "CYTOKERATIN 8"            = "KRT8",
  "FIBRONECTIN"              = "FN1",
  "GLUCOCORTICOID RECEPTOR"  = "NR3C1",
  "GLYCOPHORIN A"            = "GYPA",
  "GRP78"                    = "HSPA5",
  "IBA1"                     = "AIF1",
  "IL-1RA"                   = "IL1RN",
  "INTEGRIN-Β1"              = "ITGB1",
  "KIM-1"                    = "HAVCR1",
  "LC3B"                     = "MAP1LC3B",
  "NEPHRIN"                  = "NPHS1",
  "P-C-JUN"                  = "JUN",
  "P-EZRIN"                  = "EZR",
  "P-RIBOSOMAL PROTEIN S6"   = "RPS6",
  "P-STAT3"                  = "STAT3",
  "P62/SQSTM1"               = "SQSTM1",
  "PDI"                      = "P4HB",
  "PHOSPHO-HISTONE H3"       = "H3-3A",
  "PHOSPHO-ERK1/2"           = "MAPK3",
  "PROTEASOME 20S LMP7"      = "PSMB9",
  "RAB7"                     = "RAB7A",
  "SR-B1"                    = "SCARB1",
  "TALIN1"                   = "TLN1",
  "UBIQUITYL-HISTONE H2B"    = "H2BC11",
  "VIMENTIN"                 = "VIM"
)

s3 <- s3 %>%
  mutate(
    gene_symbol_mapped = coalesce(alias_map[gene_symbol], gene_symbol),
    gene_symbol_mapped = na_if(gene_symbol_mapped, "NA")
  )

# ── Load S5 scRNA-seq results ─────────────────────────────────────────────────
s5_raw <- read_excel(s5_path, skip = 1)
colnames(s5_raw)[1:9] <- c("drop", "cell_type", "group", "comparison_group",
                           "gene_symbol", "mean_group", "mean_comparison",
                           "p_value", "log2fc")
s5 <- s5_raw %>%
  filter(!is.na(gene_symbol), gene_symbol != "Gene symbol") %>%
  select(-drop) %>%
  mutate(across(c(p_value, log2fc, mean_group, mean_comparison), as.numeric),
         gene_symbol = str_trim(str_to_upper(gene_symbol)))

s5_hits <- s5 %>%
  filter(gene_symbol %in% na.omit(s3$gene_symbol_mapped))

# ── Overlap: S3 vs proteomics (using limma results) ───────────────────────────
prot_genes <- rownames(prot_imp)

overlap <- s3 %>%
  mutate(
    in_proteomics = gene_symbol_mapped %in% prot_genes,
    in_phospho    = gene_symbol_mapped %in% rownames(phos_imp)
  )

cat("\n── PathoPlex S3 Overlap ─────────────────────────────────────────────────\n")
cat(sprintf("In proteomics:        %d / %d (%.1f%%)\n",
            sum(overlap$in_proteomics, na.rm=TRUE), nrow(s3),
            100*mean(overlap$in_proteomics, na.rm=TRUE)))
cat(sprintf("In phosphoproteomics: %d / %d (%.1f%%)\n",
            sum(overlap$in_phospho, na.rm=TRUE), nrow(s3),
            100*mean(overlap$in_phospho, na.rm=TRUE)))

# ── Merge limma FC with PathoPlex proteins ────────────────────────────────────
# Use res_prot (limma, age+sex adjusted) instead of simple mean FC
patho_limma <- overlap %>%
  filter(in_proteomics) %>%
  inner_join(
    res_prot %>% select(feature_id, logFC, P.Value, adj.P.Val),
    by = c("gene_symbol_mapped" = "feature_id")
  )

cat("\nPathoPlex proteins in proteomics with limma results:\n")
print(patho_limma %>%
        select(protein_label, gene_symbol_mapped, logFC, P.Value, adj.P.Val) %>%
        arrange(P.Value), n = Inf)

# ── Reversal analysis ─────────────────────────────────────────────────────────
s5_dkd_avg <- s5_hits %>%
  filter(p_value < 0.05, group == "DKD") %>%
  group_by(gene_symbol) %>%
  summarise(log2fc_dkd_scrna = mean(log2fc, na.rm = TRUE), .groups = "drop")

reversal <- patho_limma %>%
  inner_join(s5_dkd_avg, by = c("gene_symbol_mapped" = "gene_symbol")) %>%
  mutate(
    reversal_score  = logFC * log2fc_dkd_scrna,
    reversed        = reversal_score < 0,
    direction_label = case_when(
      reversed & log2fc_dkd_scrna > 0 ~ "Up in DKD → Down with SGLT2i",
      reversed & log2fc_dkd_scrna < 0 ~ "Down in DKD → Up with SGLT2i",
      !reversed & log2fc_dkd_scrna > 0 ~ "Up in DKD → Up with SGLT2i",
      !reversed & log2fc_dkd_scrna < 0 ~ "Down in DKD → Down with SGLT2i"
    )
  )

cat("\n── Reversal Summary ─────────────────────────────────────────────────────\n")
print(reversal %>% count(direction_label, reversed) %>% arrange(desc(n)))
cat(sprintf("\n%d / %d proteins (%.1f%%) show reversal of DKD direction with SGLT2i\n",
            sum(reversal$reversed), nrow(reversal),
            100 * mean(reversal$reversed)))

cat("\nTop reversed proteins (strongest reversal score):\n")
print(reversal %>% filter(reversed) %>%
        arrange(reversal_score) %>%
        select(protein_label, gene_symbol_mapped, log2fc_dkd_scrna,
               logFC, P.Value, reversal_score, direction_label), n = Inf)

# ── Dotplot: all overlapping proteins, z-scored, using limma FC ───────────────
s5_sglt2i_avg <- s5_hits %>%
  filter(p_value < 0.05, group == "SGLT2i") %>%
  group_by(gene_symbol) %>%
  summarise(log2fc_scrna = mean(log2fc, na.rm = TRUE), .groups = "drop")

# Z-score within modality
prot_scaled <- patho_limma %>%
  select(GENE = gene_symbol_mapped, protein_label, logFC) %>%
  mutate(log2fc_prot_z = as.numeric(scale(logFC)))

scrna_scaled <- bind_rows(
  s5_sglt2i_avg %>% mutate(comparison = "SGLT2i vs Control"),
  s5_dkd_avg    %>% rename(log2fc_scrna = log2fc_dkd_scrna) %>%
    mutate(comparison = "DKD vs Control")
) %>%
  group_by(comparison) %>%
  mutate(log2fc_scrna_z = as.numeric(scale(log2fc_scrna))) %>%
  ungroup()

combined <- prot_scaled %>%
  left_join(scrna_scaled, by = c("GENE" = "gene_symbol")) %>%
  mutate(
    comparison   = replace_na(comparison, "Proteomics only"),
    has_scrna    = !is.na(log2fc_scrna_z),
    concordant_label = case_when(
      !has_scrna ~ NA_character_,
      sign(log2fc_prot_z) == sign(log2fc_scrna_z) ~ "Concordant",
      TRUE ~ "Discordant"
    )
  ) %>%
  mutate(GENE = factor(GENE, levels = prot_scaled %>%
                         arrange(log2fc_prot_z) %>% pull(GENE)))

pt_all <- bind_rows(
  combined %>% select(GENE, comparison, log2fc_prot_z) %>%
    mutate(modality = "Proteomics (Rinschen)", log2fc_z = log2fc_prot_z) %>%
    distinct(),
  combined %>% filter(has_scrna) %>%
    select(GENE, comparison, log2fc_scrna_z) %>%
    mutate(modality = "scRNA-seq (S5)", log2fc_z = log2fc_scrna_z)
)

seg_all <- combined %>%
  filter(has_scrna) %>%
  select(GENE, comparison, concordant_label, log2fc_prot_z, log2fc_scrna_z)

p_patho_final <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_segment(data = seg_all,
               aes(x = log2fc_scrna_z, xend = log2fc_prot_z,
                   y = GENE, yend = GENE, color = concordant_label),
               linewidth = 0.8, alpha = 0.6) +
  scale_color_manual(name = "Direction",
                     values = c("Concordant" = "#00A087", "Discordant" = "#E64B35"),
                     na.value = NA) +
  new_scale_color() +
  geom_point(data = pt_all,
             aes(x = log2fc_z, y = GENE, color = modality, shape = modality),
             size = 2.5, alpha = 0.9, position = position_dodge(width = 0.5)) +
  scale_color_manual(name = "Modality",
                     values = c("Proteomics (Rinschen)" = "#F39B7F",
                                "scRNA-seq (S5)"        = "#3C5488")) +
  scale_shape_manual(name = "Modality",
                     values = c("Proteomics (Rinschen)" = 17, "scRNA-seq (S5)" = 16)) +
  facet_wrap(~ comparison, nrow = 1) +
  labs(title    = "PathoPlex Proteins: Rinschen Proteomics vs S5 scRNA-seq",
       subtitle = "Limma log2FC (age+sex adjusted), z-scored within modality",
       x = "Scaled log2FC (z-score)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8)) +
  guides(color = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 1))

print(p_patho_final)
ggsave(file.path(base_dir, "plots/patho_dotplot_limma.pdf"),
       p_patho_final, width = 14, height = 10)

# ── Save reversal results ─────────────────────────────────────────────────────
write.xlsx(list(
  Overlap    = overlap %>% select(protein_label, gene_symbol_mapped,
                                  target, in_proteomics, in_phospho),
  Limma_FC   = patho_limma %>% select(protein_label, gene_symbol_mapped,
                                      target, logFC, P.Value, adj.P.Val),
  Reversal   = reversal %>% select(protein_label, gene_symbol_mapped,
                                   log2fc_dkd_scrna, logFC, P.Value,
                                   reversal_score, reversed, direction_label)
), file.path(base_dir, "results/PathoPlex_comparison.xlsx"), overwrite = TRUE)

cat("\nSection 10 complete — saved plots/patho_dotplot_limma.pdf\n")
cat("and results/PathoPlex_comparison.xlsx\n")



