##############################################################################
# Proteomics + Phosphoproteomics Pipeline
# Human Kidney Biopsies | Diabetes / SGLT2i
# Data: Rinschen lab (woB1 = batch 1 excluded)
# Files: hKidneyBiopsies_PROT_Expression_woB1.xlsx
#        hKidneyBiopsies_PHOS_Expression_woB1.xlsx
#        Samples Shipped.xlsx
##############################################################################

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────
pkgs <- c(
  "tidyverse", "readxl", "limma", "vsn",
  "pheatmap", "ggrepel", "patchwork",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "corrplot", "broom", "janitor", "writexl"
)
# BiocManager::install(setdiff(pkgs, rownames(installed.packages())))
invisible(lapply(pkgs, library, character.only = TRUE))

# ── 1. PATHS ─────────────────────────────────────────────────────────────────
base_dir <- "/Users/pylell/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/Rinschen"
dir.create(file.path(base_dir, "results"), showWarnings = FALSE)
dir.create(file.path(base_dir, "plots"),   showWarnings = FALSE)

# ── 2. LOAD DATA ──────────────────────────────────────────────────────────────
# -- Proteomics (protein-level abundance)
prot_raw <- read_xlsx(file.path(base_dir, "hKidneyBiopsies_PROT_Expression_woB1.xlsx"))

# -- Phosphoproteomics (phosphosite-level abundance)
phos_raw <- read_xlsx(file.path(base_dir, "hKidneyBiopsies_PHOS_Expression_woB1.xlsx"))

# -- Sample metadata
meta_raw <- read_xlsx(file.path(base_dir, "Samples Shipped.xlsx")) %>%
  clean_names()

# ── INSPECT BEFORE PROCEEDING ─────────────────────────────────────────────────
# Run these lines first and adjust column names below to match your actual data
glimpse(prot_raw)
glimpse(phos_raw)
glimpse(meta_raw)

##############################################################################
# SECTION 0: DATA WRANGLING
# !!! Adjust these column name assignments after inspecting the files !!!
##############################################################################

# -- Expected: first column(s) = protein/site identifiers, rest = sample columns
# Proteomics identifier column (e.g., "Gene.names", "Protein.IDs", "gene_symbol")
PROT_ID_COL  <- "gene_symbol"     # <-- change to match actual column name

# Phosphoproteomics identifier columns
PHOS_ID_COL  <- "phosphosite_id"  # <-- e.g. "Gene_pSite", "Sequence.window"
PHOS_GENE_COL <- "gene_symbol"    # <-- gene name column in phospho file

# Metadata columns
META_SAMPLE_COL  <- "sample_id"   # <-- sample ID column matching matrix column names
META_GROUP_COL   <- "group"       # <-- treatment/disease group (e.g., "SGLT2i", "DM", "Control")
META_COHORT_COL  <- "cohort"      # <-- cohort label (RH, CRC, IT)

# Clinical variables available in metadata (adjust to what's actually there)
CLINICAL_VARS <- c("gfr", "uacr", "hba1c", "egfr_cric")

# ── Build abundance matrices ──────────────────────────────────────────────────
make_matrix <- function(df, id_col) {
  df %>%
    select(-any_of(setdiff(names(df)[!names(df) %in% meta_raw[[META_SAMPLE_COL]]], id_col))) %>%
    column_to_rownames(id_col) %>%
    as.matrix()
}

# Separate identifier columns from numeric sample columns
prot_ids  <- prot_raw %>% select(all_of(PROT_ID_COL))
prot_mat  <- prot_raw %>%
  select(all_of(PROT_ID_COL),
         any_of(meta_raw[[META_SAMPLE_COL]])) %>%
  column_to_rownames(PROT_ID_COL) %>%
  as.matrix()

phos_ids  <- phos_raw %>% select(all_of(c(PHOS_ID_COL, PHOS_GENE_COL)))
phos_mat  <- phos_raw %>%
  select(all_of(PHOS_ID_COL),
         any_of(meta_raw[[META_SAMPLE_COL]])) %>%
  column_to_rownames(PHOS_ID_COL) %>%
  as.matrix()

meta <- meta_raw %>% rename(
  sample_id = all_of(META_SAMPLE_COL),
  group     = all_of(META_GROUP_COL),
  cohort    = all_of(META_COHORT_COL)
) %>%
  filter(sample_id %in% colnames(prot_mat))

cat("Proteomics:     ", nrow(prot_mat), "proteins x", ncol(prot_mat), "samples\n")
cat("Phosphoproteomics:", nrow(phos_mat), "sites x", ncol(phos_mat), "samples\n")
cat("Groups:", table(meta$group), "\n")

##############################################################################
# SECTION 1: QC AND PREPROCESSING (applied to both PROT and PHOS)
##############################################################################

run_qc_and_norm <- function(mat, meta_df, label = "PROT") {
  
  cat("\n--- QC:", label, "---\n")
  
  # Missingness
  miss <- apply(mat, 1, function(x) mean(is.na(x))) * 100
  cat("Median missingness per feature:", round(median(miss), 1), "%\n")
  
  hist_p <- ggplot(data.frame(pct_missing = miss), aes(pct_missing)) +
    geom_histogram(bins = 40, fill = "#4E79A7") +
    labs(title = paste(label, "— Missingness distribution"),
         x = "% Missing per feature", y = "Count") +
    theme_bw()
  
  ggsave(file.path(base_dir, paste0("plots/", label, "_missingness.pdf")),
         hist_p, width = 6, height = 4)
  
  # Filter: keep features present in ≥ 70% of samples
  mat_filt <- mat[miss <= 30, ]
  cat("Features after missingness filter:", nrow(mat_filt), "\n")
  
  # Log2 transform if values look like raw intensities (max > 100)
  if (max(mat_filt, na.rm = TRUE) > 100) {
    mat_filt <- log2(mat_filt + 1)
    cat("Log2 transformation applied\n")
  }
  
  # VSN normalization
  mat_vsn <- justvsn(mat_filt)
  
  # PCA
  mat_imp_pca <- mat_vsn
  row_meds    <- apply(mat_imp_pca, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(mat_imp_pca)))
    mat_imp_pca[i, is.na(mat_imp_pca[i, ])] <- row_meds[i]
  
  pca     <- prcomp(t(mat_imp_pca), scale. = TRUE)
  var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  pca_df  <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("sample_id") %>%
    left_join(meta_df, by = "sample_id")
  
  pca_p <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = cohort)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 12) +
    labs(title = paste(label, "— PCA (VSN normalized)"),
         x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw()
  
  ggsave(file.path(base_dir, paste0("plots/", label, "_PCA.pdf")),
         pca_p, width = 7, height = 5)
  
  # Distribution boxplot
  box_df <- mat_vsn %>% as.data.frame() %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample_id", values_to = "intensity") %>%
    left_join(meta_df %>% select(sample_id, group, cohort), by = "sample_id")
  
  box_p <- ggplot(box_df, aes(sample_id, intensity, fill = group)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    labs(title = paste(label, "— Normalized intensities"), x = NULL)
  
  ggsave(file.path(base_dir, paste0("plots/", label, "_boxplot_norm.pdf")),
         box_p, width = 10, height = 4)
  
  # MinProb imputation for downstream analysis
  mat_final <- mat_vsn
  for (j in seq_len(ncol(mat_final))) {
    col      <- mat_final[, j]
    miss_idx <- is.na(col)
    if (any(miss_idx)) {
      lv <- quantile(col, 0.01, na.rm = TRUE)
      sd_v <- sd(col, na.rm = TRUE) * 0.3
      mat_final[miss_idx, j] <- rnorm(sum(miss_idx), lv, sd_v)
    }
  }
  
  mat_final
}

prot_imp <- run_qc_and_norm(prot_mat, meta, label = "PROT")
phos_imp <- run_qc_and_norm(phos_mat, meta, label = "PHOS")

##############################################################################
# SECTION 2: DIFFERENTIAL ABUNDANCE — PROT and PHOS
# Primary contrast: SGLT2i vs. Diabetic Control
# Adjust group levels to match your actual labels in the metadata
##############################################################################

run_limma <- function(mat, meta_df,
                      case_label,        # e.g. "SGLT2i"
                      ctrl_label,        # e.g. "DM_control"
                      covariates = c("age", "sex", "bmi"),
                      label = "PROT") {
  
  samps   <- intersect(colnames(mat), meta_df$sample_id)
  mat_s   <- mat[, samps]
  md      <- meta_df %>%
    filter(sample_id %in% samps) %>%
    arrange(match(sample_id, colnames(mat_s))) %>%
    mutate(group = factor(group, levels = c(ctrl_label, case_label)))
  
  # Drop covariates not in metadata
  covariates <- intersect(covariates, names(md))
  cov_str    <- if (length(covariates)) paste("+", paste(covariates, collapse = " + ")) else ""
  design     <- model.matrix(as.formula(paste("~ 0 + group", cov_str)), data = md)
  colnames(design) <- make.names(colnames(design))
  
  grp_case <- make.names(paste0("group", case_label))
  grp_ctrl <- make.names(paste0("group", ctrl_label))
  cont_mat  <- makeContrasts(
    contrasts = paste0(grp_case, " - ", grp_ctrl),
    levels    = design
  )
  
  fit  <- lmFit(mat_s, design)
  fit2 <- contrasts.fit(fit, cont_mat)
  fit2 <- eBayes(fit2)
  
  res  <- topTable(fit2, coef = 1, number = Inf, sort.by = "P") %>%
    rownames_to_column("feature_id") %>%
    mutate(data_type = label)
  
  cat(label, "—", sum(res$adj.P.Val < 0.05), "significant features (FDR < 0.05)\n")
  res
}

# !!! Set your group labels here !!!
CASE_LABEL <- "SGLT2i"       # treatment / post group
CTRL_LABEL <- "DM_control"  # comparator group

res_prot <- run_limma(prot_imp, meta,
                      case_label = CASE_LABEL, ctrl_label = CTRL_LABEL,
                      label = "PROT")

res_phos <- run_limma(phos_imp, meta,
                      case_label = CASE_LABEL, ctrl_label = CTRL_LABEL,
                      label = "PHOS")

write_xlsx(list(Proteomics      = res_prot,
                Phosphoproteomics = res_phos),
           file.path(base_dir, "results/differential_abundance_PROT_PHOS.xlsx"))

# ── 2A. Volcano plots ─────────────────────────────────────────────────────────
plot_volcano <- function(res, title = "Volcano", fc_cut = 0.5, fdr_cut = 0.05) {
  res %>%
    mutate(
      sig   = adj.P.Val < fdr_cut & abs(logFC) > fc_cut,
      dir   = case_when(
        sig & logFC >  fc_cut ~ "Up (SGLT2i)",
        sig & logFC < -fc_cut ~ "Down (SGLT2i)",
        TRUE                  ~ "NS"
      ),
      label = ifelse(sig, feature_id, NA_character_)
    ) %>%
    ggplot(aes(logFC, -log10(P.Value), color = dir)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20) +
    scale_color_manual(values = c("Up (SGLT2i)"   = "#E15759",
                                  "Down (SGLT2i)" = "#4E79A7",
                                  NS              = "grey70")) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(title = title, x = "log2 Fold Change (SGLT2i vs DM)",
         y = "-log10(p-value)", color = NULL) +
    theme_bw()
}

p_volc <- plot_volcano(res_prot, "Proteomics — SGLT2i vs DM") |
  plot_volcano(res_phos, "Phosphoproteomics — SGLT2i vs DM")
ggsave(file.path(base_dir, "plots/volcano_PROT_PHOS.pdf"), p_volc, width = 14, height = 6)

##############################################################################
# SECTION 3: PHOS/PROT NORMALIZATION
# Adjusts phosphosite changes for corresponding total protein changes
# (i.e., is the phosphorylation change driven by protein abundance?)
##############################################################################

# Join phospho gene to protein-level logFC
phos_norm <- res_phos %>%
  left_join(phos_ids, by = c("feature_id" = PHOS_ID_COL)) %>%
  left_join(res_prot %>% select(feature_id, logFC) %>%
              rename(prot_logFC = logFC, gene_symbol = feature_id),
            by = c(PHOS_GENE_COL := "gene_symbol")) %>%
  mutate(
    logFC_adjusted = logFC - prot_logFC,   # phospho signal independent of protein level
    site_driven    = abs(prot_logFC) > 0.5 & sign(logFC) == sign(prot_logFC)
  )

write_csv(phos_norm, file.path(base_dir, "results/phospho_normalized_to_protein.csv"))

# Scatter: raw phospho FC vs protein FC
ggplot(phos_norm %>% filter(!is.na(prot_logFC)),
       aes(prot_logFC, logFC, color = site_driven)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("TRUE" = "#E15759", "FALSE" = "#4E79A7"),
                     labels = c("TRUE" = "Protein-driven", "FALSE" = "Phospho-specific")) +
  labs(title = "Phospho FC vs. Protein FC (SGLT2i vs DM)",
       x = "log2FC (Protein)", y = "log2FC (Phosphosite)", color = NULL) +
  theme_bw()
ggsave(file.path(base_dir, "plots/phospho_vs_protein_FC.pdf"), width = 7, height = 6)

##############################################################################
# SECTION 4: KINASE ENRICHMENT ANALYSIS (PhosphoSitePlus / KSEA)
# Identifies which kinases are differentially active based on substrate changes
##############################################################################

# Install if needed: remotes::install_github("casecpb/KSEA")
# library(KSEA)
#
# ksea_input <- phos_norm %>%
#   select(feature_id, logFC_adjusted, adj.P.Val) %>%
#   rename(SUB_GENE_PTM = feature_id, log2FC = logFC_adjusted, p = adj.P.Val)
#
# ksea_res <- KSEA.Scores(KSEAapp::PhosphoSitePlus,
#                          ksea_input,
#                          NetworKIN     = TRUE,
#                          NetworKIN.cutoff = 5)
# write_csv(ksea_res, file.path(base_dir, "results/kinase_enrichment_KSEA.csv"))

# Manual alternative: enrichment via known kinase-substrate databases
# Useful if KSEA package is unavailable
kinase_substrate_db <- read_csv("https://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz")
# (requires free PSP account; save locally first)

##############################################################################
# SECTION 5: PATHWAY ENRICHMENT — PROT and PHOS
##############################################################################

run_enrichment <- function(sig_features, universe, label = "") {
  ego <- tryCatch(
    enrichGO(gene          = sig_features,
             universe      = universe,
             OrgDb         = org.Hs.eg.db,
             keyType       = "SYMBOL",
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             readable      = TRUE),
    error = function(e) { message("GO failed: ", e$message); NULL }
  )
  
  entrez  <- bitr(sig_features, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
  univ_ez <- bitr(universe, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
  
  ekegg <- tryCatch(
    enrichKEGG(gene         = entrez$ENTREZID,
               universe      = univ_ez$ENTREZID,
               organism      = "hsa",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05),
    error = function(e) { message("KEGG failed: ", e$message); NULL }
  )
  
  if (!is.null(ego) && nrow(ego) > 0) {
    ggsave(file.path(base_dir, paste0("plots/GO_BP_", label, ".pdf")),
           dotplot(ego, showCategory = 20,
                   title = paste("GO BP —", label)),
           width = 9, height = 8)
  }
  
  list(go = ego, kegg = ekegg)
}

# PROT enrichment
sig_prot_up   <- res_prot %>% filter(adj.P.Val < 0.05, logFC >  0.5) %>% pull(feature_id)
sig_prot_down <- res_prot %>% filter(adj.P.Val < 0.05, logFC < -0.5) %>% pull(feature_id)
run_enrichment(sig_prot_up,   rownames(prot_imp), "PROT_up_SGLT2i")
run_enrichment(sig_prot_down, rownames(prot_imp), "PROT_down_SGLT2i")

# PHOS enrichment (use gene symbols from phosphosite IDs)
sig_phos_genes_up   <- phos_norm %>%
  filter(adj.P.Val < 0.05, logFC_adjusted >  0.5) %>%
  pull(all_of(PHOS_GENE_COL)) %>% unique()
sig_phos_genes_down <- phos_norm %>%
  filter(adj.P.Val < 0.05, logFC_adjusted < -0.5) %>%
  pull(all_of(PHOS_GENE_COL)) %>% unique()
run_enrichment(sig_phos_genes_up,   unique(phos_ids[[PHOS_GENE_COL]]), "PHOS_up_SGLT2i")
run_enrichment(sig_phos_genes_down, unique(phos_ids[[PHOS_GENE_COL]]), "PHOS_down_SGLT2i")

##############################################################################
# SECTION 6: CLINICAL CORRELATIONS — GFR, uACR, etc.
##############################################################################

correlate_with_clinical <- function(mat, meta_df, clin_vars, label = "PROT") {
  sig_feats <- switch(label,
                      "PROT" = res_prot %>% filter(adj.P.Val < 0.05) %>% pull(feature_id),
                      "PHOS" = phos_norm %>% filter(adj.P.Val < 0.05) %>% pull(feature_id)
  )
  sig_feats <- intersect(sig_feats, rownames(mat))
  
  expand_grid(feature = sig_feats, clinical = clin_vars) %>%
    rowwise() %>%
    mutate(
      test = list({
        x <- mat[feature, meta_df$sample_id]
        y <- meta_df[[clinical]]
        n <- sum(!is.na(x) & !is.na(y))
        if (n >= 8) {
          ct <- cor.test(x, y, method = "spearman")
          tibble(r = ct$estimate, p = ct$p.value, n = n)
        } else tibble(r = NA_real_, p = NA_real_, n = n)
      })
    ) %>%
    unnest(test) %>%
    ungroup() %>%
    group_by(clinical) %>%
    mutate(fdr = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    mutate(data_type = label)
}

avail_clin <- intersect(CLINICAL_VARS, names(meta))
corr_prot  <- correlate_with_clinical(prot_imp, meta, avail_clin, "PROT")
corr_phos  <- correlate_with_clinical(phos_imp, meta, avail_clin, "PHOS")

write_xlsx(list(PROT = corr_prot, PHOS = corr_phos),
           file.path(base_dir, "results/clinical_correlations_PROT_PHOS.xlsx"))

# Heatmap for top protein-clinical correlations
top_corr <- corr_prot %>%
  filter(fdr < 0.05, !is.na(r)) %>%
  select(feature, clinical, r) %>%
  pivot_wider(names_from = clinical, values_from = r) %>%
  column_to_rownames("feature") %>%
  as.matrix()

if (nrow(top_corr) >= 2) {
  pheatmap(top_corr, na_col = "grey90",
           main = "Protein–clinical correlations (FDR < 0.05)",
           filename = file.path(base_dir, "plots/clinical_corr_heatmap_PROT.pdf"),
           width = 8, height = max(4, nrow(top_corr) * 0.15 + 2))
}

##############################################################################
# SECTION 7: scRNA-SEQ INTEGRATION
##############################################################################

# Load your top scRNA/snRNA-seq hits
# scrna_hits <- read_csv(file.path(base_dir, "scrna_top_hits.csv"))
# Expected cols: gene_symbol, cell_type, log2FC_rna, fdr_rna

# Concordance between protein and mRNA fold changes
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
#        x = "log2FC scRNA-seq", y = "log2FC Proteomics") +
#   theme_bw()

##############################################################################
# SECTION 8: SUMMARY FOR ASN ABSTRACT
##############################################################################

summary_table <- res_prot %>%
  filter(adj.P.Val < 0.05) %>%
  left_join(
    corr_prot %>%
      filter(clinical %in% c("gfr", "uacr"), fdr < 0.05) %>%
      select(feature, clinical, r) %>%
      pivot_wider(names_from = clinical, values_from = r,
                  names_prefix = "r_"),
    by = c("feature_id" = "feature")
  ) %>%
  left_join(
    phos_norm %>%
      filter(adj.P.Val < 0.05) %>%
      count(!!sym(PHOS_GENE_COL), name = "n_sig_phosphosites") %>%
      rename(feature_id = all_of(PHOS_GENE_COL)),
    by = "feature_id"
  ) %>%
  arrange(adj.P.Val)

write_xlsx(list(Summary = summary_table),
           file.path(base_dir, "results/summary_table_ASN_abstract.xlsx"))

cat("\n=== Pipeline complete ===\n",
    "results/differential_abundance_PROT_PHOS.xlsx\n",
    "results/phospho_normalized_to_protein.csv\n",
    "results/clinical_correlations_PROT_PHOS.xlsx\n",
    "results/summary_table_ASN_abstract.xlsx\n",
    "plots/ — QC, volcano, pathway, correlation plots\n")