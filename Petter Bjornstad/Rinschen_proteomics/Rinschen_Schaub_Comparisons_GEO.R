################################################################################
# PROTEOMICS vs. scRNA-seq COMPARISON - Mac version
# Three proteomics datasets vs Schaub GEO NEBULA results:
#   A) Rinschen kidney proteomics - re-run limma from raw
#   B) Rinschen kidney proteomics - pre-computed OneDrive results
#   C) SomaScan blood proteomics - from harmonized dataset
# Output: /Users/jdw/Documents/UofW/Projects/Rinschen_Comparisons/
################################################################################

# ── 0. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("tidyverse", "readxl", "limma", "openxlsx",
          "ggrepel", "patchwork", "janitor")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))

# ── 1. Paths ──────────────────────────────────────────────────────────────────
rinschen_dir  <- path.expand(
  "~/Documents/OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/Rinschen/")
harmonized_dir <- path.expand(
  "~/Documents/OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/")
output_dir    <- "/Users/jdw/Documents/UofW/Projects/Rinschen_Comparisons/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# NEBULA results from GEO replication script
nebula_file <- file.path(output_dir, "GEO_replication_results.xlsx")
if (!file.exists(nebula_file))
  stop("Cannot find GEO_replication_results.xlsx in:\n  ", output_dir,
       "\nRun Schaub_GEO_NEBULA_comparison.R first.")

# ── 2. Pathway gene sets ──────────────────────────────────────────────────────
pathway_genes <- list(
  Glycolysis       = c("PKLR","PFKFB3","PFKL","ALDOC","HK2",
                       "ENO2","PGK1","PGAM1","TPI1","GAPDH"),
  Gluconeogenesis  = c("SLC25A10","GOT2","GOT1","FBP1","SLC25A11","PCK1","MDH1"),
  TCA_cycle        = c("SDHB","SUCLG1","PDK2","ACO2","IDH3G","SUCLA2",
                       "HAGH","PDHB","LDHA"),
  Glutathione      = c("CNDP2","GSTM4","GSTT2B","GSTO1","GGCT","GSTM3","AKR1A1"),
  Metallothioneins = c("MT1G","MT1X","MT1H","MT2A")
)
all_pathway_genes <- unlist(pathway_genes)

gene_to_pathway <- function(g) case_when(
  g %in% pathway_genes$Glycolysis       ~ "Glycolysis",
  g %in% pathway_genes$Gluconeogenesis  ~ "Gluconeogenesis",
  g %in% pathway_genes$TCA_cycle        ~ "TCA cycle",
  g %in% pathway_genes$Glutathione      ~ "Glutathione",
  g %in% pathway_genes$Metallothioneins ~ "Metallothioneins",
  TRUE ~ "Other"
)

# ── 3. Load NEBULA results ────────────────────────────────────────────────────
nebula_res <- read_xlsx(nebula_file, sheet = "NEBULA_GEO") %>%
  dplyr::select(gene_symbol, log2FC_nebula, pval_nebula, padj_nebula, pathway)
cat("NEBULA results loaded:", nrow(nebula_res), "genes\n")

# ── 4. Shared helpers ─────────────────────────────────────────────────────────
pathway_colors <- c("Glycolysis"       = "#E76F51",
                    "Gluconeogenesis"  = "#F4A261",
                    "TCA cycle"        = "#2A9D8F",
                    "Glutathione"      = "#457B9D",
                    "Metallothioneins" = "#6A0572")

make_combined <- function(prot_res, nebula_df, label) {
  prot_res %>%
    dplyr::filter(pathway != "Other") %>%
    left_join(nebula_df %>% dplyr::select(gene_symbol, log2FC_nebula,
                                          padj_nebula),
              by = "gene_symbol") %>%
    mutate(concordant = sign(log2FC_prot) == sign(log2FC_nebula),
           pathway    = factor(pathway, levels = c("Glycolysis","Gluconeogenesis",
                                                   "TCA cycle","Glutathione",
                                                   "Metallothioneins")),
           analysis   = label)
}

print_concordance <- function(combined_df, label) {
  cat("\n---", label, "---\n")
  cat("Overall:", sum(combined_df$concordant, na.rm=TRUE), "/",
      sum(!is.na(combined_df$concordant)), "genes concordant\n")
  combined_df %>%
    group_by(pathway) %>%
    summarise(concordant = sum(concordant, na.rm=TRUE),
              n          = n(),
              pct        = round(100*concordant/n, 1),
              .groups    = "drop") %>%
    print()
}

make_scatter <- function(combined_df, title_label, y_label = "Proteomics limma") {
  n_conc <- sum(combined_df$concordant, na.rm = TRUE)
  n_tot  <- sum(!is.na(combined_df$concordant))
  combined_df %>%
    dplyr::filter(!is.na(log2FC_nebula), !is.na(log2FC_prot)) %>%
    ggplot(aes(x = log2FC_nebula, y = log2FC_prot,
               color = pathway, label = gene_symbol)) +
    geom_point(size = 3) +
    geom_text_repel(size = 2.8, max.overlaps = 25) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30") +
    scale_color_manual(values = pathway_colors) +
    labs(title    = title_label,
         subtitle = paste0("Concordance: ", n_conc, "/", n_tot,
                           " genes (", round(100*n_conc/n_tot,1), "%)"),
         x     = "log2FC - scRNA NEBULA (GEO, SGLT2i+ vs SGLT2i-)",
         y     = paste0("log2FC - ", y_label, " (SGLT2i vs Control)"),
         color = "Pathway") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

make_bar <- function(combined_df, title_label, prot_label = "Proteomics limma") {
  method_colors <- c("scRNA NEBULA (GEO)" = "#2A9D8F")
  method_colors[prot_label] <- "#E07B39"
  combined_df %>%
    dplyr::select(gene_symbol, pathway,
                  `scRNA NEBULA (GEO)` = log2FC_nebula) %>%
    dplyr::mutate(!!prot_label := combined_df$log2FC_prot) %>%
    pivot_longer(-c(gene_symbol, pathway),
                 names_to = "method", values_to = "log2FC") %>%
    mutate(method      = factor(method, levels = names(method_colors)),
           gene_symbol = factor(gene_symbol,
                                levels = combined_df$gene_symbol[
                                  order(combined_df$pathway)])) %>%
    ggplot(aes(x = gene_symbol, y = log2FC, fill = method)) +
    geom_col(position = position_dodge(0.75), width = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    facet_wrap(~ pathway, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = method_colors, name = NULL) +
    labs(title = title_label,
         subtitle = "SGLT2i+ vs SGLT2i- | Positive = higher in SGLT2i",
         x = NULL, y = "log2 Fold Change") +
    theme_bw(base_size = 11) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
          legend.position  = "bottom",
          strip.background = element_rect(fill = "grey92"),
          strip.text       = element_text(face = "bold", size = 10),
          plot.title       = element_text(face = "bold"))
}

run_limma_on_mat <- function(mat, meta_sub, genes_use) {
  genes_use <- genes_use[genes_use %in% rownames(mat)]
  if (length(genes_use) == 0) {
    cat("  No pathway genes found in matrix\n"); return(NULL)
  }
  mat_sub <- mat[genes_use, meta_sub$study_id, drop = FALSE]
  if (max(mat_sub, na.rm = TRUE) > 50) {
    cat("  Log2 transforming\n"); mat_sub <- log2(mat_sub + 1)
  }
  # Remove proteins > 50% missing, impute remainder
  mat_sub <- mat_sub[rowMeans(is.na(mat_sub)) <= 0.5, , drop = FALSE]
  for (i in seq_len(nrow(mat_sub))) {
    na_i <- is.na(mat_sub[i, ])
    if (any(na_i)) mat_sub[i, na_i] <- min(mat_sub[i, !na_i], na.rm = TRUE)
  }
  design <- model.matrix(~ group, data = meta_sub)
  fit    <- limma::lmFit(mat_sub, design)
  fit    <- limma::eBayes(fit)
  limma::topTable(fit, coef = "groupSGLT2i", number = Inf,
                  adjust.method = "BH", sort.by = "none") %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    dplyr::rename(log2FC_prot = logFC, padj_prot = adj.P.Val) %>%
    dplyr::select(gene_symbol, log2FC_prot, padj_prot) %>%
    mutate(pathway = gene_to_pathway(gene_symbol))
}

results_list <- list()

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: KIDNEY PROTEOMICS (Rinschen hKidneyBiopsies)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n════ PART 1: Rinschen Kidney Proteomics ════\n")

prot_file <- file.path(rinschen_dir, "hKidneyBiopsies_PROT_Expression_woB1.xlsx")
meta_file <- file.path(rinschen_dir, "Samples Shipped.xlsx")

if (file.exists(prot_file) && file.exists(meta_file)) {
  
  prot_raw <- read_xlsx(prot_file)
  meta_raw <- read_xlsx(meta_file) %>% clean_names()
  
  # Detect gene column
  gene_col <- names(prot_raw)[grepl("^GENE$|^gene$|^Gene$",
                                    names(prot_raw))][1]
  cat("Gene column:", gene_col, "\n")
  
  # Build matrix
  quant_cols  <- names(prot_raw)[grepl("Quantity", names(prot_raw))]
  aliquot_ids <- stringr::str_extract(quant_cols,
                                      "(?<=hKidneyBiopsy_)S\\d+")
  meta <- meta_raw %>%
    filter(!is.na(treatment)) %>%
    mutate(aliquot_nodash = str_remove_all(rn_alater_id, "-"),
           group = if_else(str_detect(treatment,
                                      regex("sglt2|inhibit", ignore_case=TRUE)),
                           "SGLT2i", "Control"))
  cat("Kidney proteomics - SGLT2i:", sum(meta$group=="SGLT2i"),
      "| Control:", sum(meta$group=="Control"), "\n")
  
  keep      <- !is.na(aliquot_ids) & aliquot_ids %in% meta$aliquot_nodash
  col_map   <- meta %>% dplyr::select(aliquot_nodash, study_id) %>% deframe()
  
  # Deduplicate gene names - keep row with highest mean abundance
  prot_dedup <- prot_raw %>%
    dplyr::select(all_of(gene_col), all_of(quant_cols[keep])) %>%
    dplyr::rename(gene = all_of(gene_col)) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    mutate(row_mean = rowMeans(dplyr::across(where(is.numeric)), na.rm = TRUE)) %>%
    group_by(gene) %>%
    dplyr::slice_max(row_mean, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(-row_mean)
  
  n_total <- nrow(prot_raw)
  n_dedup <- nrow(prot_dedup)
  if (n_total != n_dedup)
    cat("  Deduplicated gene names:", n_total, "->", n_dedup, "rows\n")
  
  prot_mat  <- prot_dedup %>% column_to_rownames("gene") %>% as.matrix()
  colnames(prot_mat) <- col_map[aliquot_ids[keep]]
  matched_ids <- intersect(meta$study_id, colnames(prot_mat))
  prot_mat    <- prot_mat[, matched_ids, drop = FALSE]
  meta_k      <- meta %>% filter(study_id %in% matched_ids) %>%
    arrange(match(study_id, matched_ids))
  
  
  in_prot <- all_pathway_genes[all_pathway_genes %in% rownames(prot_mat)]
  cat("Pathway genes in kidney proteomics:", length(in_prot), "/",
      length(all_pathway_genes), "\n")
  
  # ── Analysis A: Re-run limma from raw ──────────────────────────────────────
  cat("\n─ Analysis A: Re-run limma from raw ─\n")
  prot_a <- run_limma_on_mat(prot_mat, meta_k, in_prot)
  if (!is.null(prot_a)) {
    combined_a <- make_combined(prot_a, nebula_res, "Kidney prot - re-run")
    results_list[["kidney_rerun"]] <- combined_a
    print_concordance(combined_a, "Kidney proteomics (re-run)")
    ggsave(file.path(output_dir, "kidney_prot_scatter_rerun.pdf"),
           make_scatter(combined_a, "Kidney Proteomics (Re-run) vs NEBULA"), width=9, height=7)
    ggsave(file.path(output_dir, "kidney_prot_bar_rerun.pdf"),
           make_bar(combined_a, "Kidney Proteomics (Re-run) vs NEBULA"), width=18, height=7)
  }
  
  # ── Analysis B: Pre-computed OneDrive results ──────────────────────────────
  cat("\n─ Analysis B: Pre-computed OneDrive results ─\n")
  precomp_file <- file.path(rinschen_dir,
                            "results/differential_abundance_PROT_PHOS.xlsx")
  if (file.exists(precomp_file)) {
    precomp_raw <- read_xlsx(precomp_file, sheet = "Proteomics") %>%
      clean_names()
    id_col  <- names(precomp_raw)[grepl("feature|gene|symbol",
                                        names(precomp_raw), ignore.case=TRUE)][1]
    fc_col  <- names(precomp_raw)[grepl("^logfc$|^log2fc$",
                                        names(precomp_raw), ignore.case=TRUE)][1]
    fdr_col <- names(precomp_raw)[grepl("adj|fdr|padj",
                                        names(precomp_raw), ignore.case=TRUE)][1]
    cat("Detected columns - ID:", id_col, "| FC:", fc_col, "| FDR:", fdr_col, "\n")
    
    prot_b <- precomp_raw %>%
      dplyr::rename(gene_symbol  = all_of(id_col),
                    log2FC_prot  = all_of(fc_col),
                    padj_prot    = all_of(fdr_col)) %>%
      dplyr::select(gene_symbol, log2FC_prot, padj_prot) %>%
      dplyr::filter(gene_symbol %in% all_pathway_genes) %>%
      mutate(pathway = gene_to_pathway(gene_symbol))
    
    combined_b <- make_combined(prot_b, nebula_res, "Kidney prot - pre-computed")
    results_list[["kidney_precomp"]] <- combined_b
    print_concordance(combined_b, "Kidney proteomics (pre-computed)")
    ggsave(file.path(output_dir, "kidney_prot_scatter_precomp.pdf"),
           make_scatter(combined_b, "Kidney Proteomics (Pre-computed) vs NEBULA"), width=9, height=7)
    ggsave(file.path(output_dir, "kidney_prot_bar_precomp.pdf"),
           make_bar(combined_b, "Kidney Proteomics (Pre-computed) vs NEBULA"), width=18, height=7)
  } else {
    cat("Pre-computed file not found:", precomp_file, "\n")
  }
  
} else {
  cat("Kidney proteomics files not found - skipping Part 1\n")
  cat("Expected:\n  ", prot_file, "\n  ", meta_file, "\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: BLOOD PROTEOMICS (SomaScan from harmonized dataset)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n════ PART 2: SomaScan Blood Proteomics ════\n")

# SomaScan data lives in the harmonized dataset as seq_XXXXX_XX columns.
# We need the data dictionary to map seq IDs → gene symbols.
soma_file <- file.path(harmonized_dir, "soma_harmonized_dataset.csv")
if (!file.exists(soma_file))
  soma_file <- file.path(harmonized_dir, "harmonized_dataset.csv")

dict_file <- path.expand(
  "~/Downloads/data_dictionary_master.xlsx")
# Update dict_file if your data dictionary is elsewhere

if (file.exists(soma_file)) {
  cat("Loading harmonized data:", soma_file, "\n")
  soma_raw <- read.csv(soma_file, na.strings = c("", "NA"))
  cat("Dimensions:", nrow(soma_raw), "rows x", ncol(soma_raw), "cols\n")
  
  # Find seq_ columns
  seq_cols <- names(soma_raw)[grepl("^seq[_\\.]\\d+", names(soma_raw))]
  cat("SomaScan seq columns found:", length(seq_cols), "\n")
  
  if (length(seq_cols) == 0) {
    cat("No seq_ columns found. Available column prefixes:\n")
    print(sort(table(sub("_.*", "", names(soma_raw))), decreasing=TRUE)[1:20])
  }
  
  # ── Load SomaScan annotation ──────────────────────────────────────────────
  annotation <- NULL
  if (file.exists(dict_file)) {
    cat("Loading data dictionary:", dict_file, "\n")
    dict <- read_xlsx(dict_file)
    # Look for the proteomics section with SeqId → gene symbol mapping
    cat("Dictionary sheets/columns:\n"); print(names(dict)[1:10])
    
    # Try to find gene symbol and SeqId columns
    seq_id_col  <- names(dict)[grepl("seq.?id|seqid",
                                     names(dict), ignore.case=TRUE)][1]
    gene_id_col <- names(dict)[grepl("gene.?symbol|entrez.*gene|target",
                                     names(dict), ignore.case=TRUE)][1]
    if (!is.na(seq_id_col) && !is.na(gene_id_col)) {
      annotation <- dict %>%
        dplyr::select(seq_id   = all_of(seq_id_col),
                      gene_sym = all_of(gene_id_col)) %>%
        dplyr::filter(!is.na(gene_sym), gene_sym != "") %>%
        # Normalise SeqId to match column name format (seq_XXXXX_XX)
        mutate(col_name = paste0("seq_", str_replace_all(seq_id, "-", "_")))
      cat("Annotation loaded:", nrow(annotation), "entries\n")
    } else {
      cat("Could not auto-detect SeqId/GeneSymbol columns in dictionary.\n")
      cat("Dictionary column names:\n"); print(names(dict))
    }
  } else {
    cat("Data dictionary not found at:", dict_file, "\n")
    cat("Update dict_file path above.\n")
  }
  
  # ── Find pathway gene seq IDs ─────────────────────────────────────────────
  if (!is.null(annotation)) {
    pathway_seqs <- annotation %>%
      dplyr::filter(gene_sym %in% all_pathway_genes) %>%
      dplyr::filter(col_name %in% seq_cols)
    
    cat("\nPathway genes with SomaScan seq IDs:", nrow(pathway_seqs), "\n")
    if (nrow(pathway_seqs) > 0) print(pathway_seqs)
    
    if (nrow(pathway_seqs) > 0) {
      # ── Check SGLT2i variable ───────────────────────────────────────────────
      # Common names for SGLT2i in harmonized data
      sglt_col <- names(soma_raw)[grepl("sglt2|sglt_2|sglt2i",
                                        names(soma_raw), ignore.case=TRUE)][1]
      id_col   <- names(soma_raw)[grepl("^mrn$|^record_id$|^study_id$|^id$",
                                        names(soma_raw), ignore.case=TRUE)][1]
      cat("\nDetected ID column:", id_col, "| SGLT2i column:", sglt_col, "\n")
      cat("SGLT2i values:\n")
      print(table(soma_raw[[sglt_col]], useNA="always"))
      
      # Build metadata
      meta_soma <- soma_raw %>%
        dplyr::select(study_id = all_of(id_col),
                      sglt2i   = all_of(sglt_col)) %>%
        dplyr::filter(!is.na(sglt2i)) %>%
        mutate(group = case_when(
          sglt2i %in% c(1, "1", "Yes", "YES", "yes", "SGLT2i") ~ "SGLT2i",
          sglt2i %in% c(0, "0", "No",  "NO",  "no",  "Control", "None") ~ "Control",
          TRUE ~ NA_character_
        )) %>%
        dplyr::filter(!is.na(group))
      
      cat("\nBlood proteomics - SGLT2i:", sum(meta_soma$group=="SGLT2i"),
          "| Control:", sum(meta_soma$group=="Control"), "\n")
      
      # Build matrix: rows = gene symbols, cols = study_id
      soma_mat <- soma_raw %>%
        dplyr::filter(!is.na(.data[[id_col]])) %>%
        dplyr::select(study_id = all_of(id_col),
                      all_of(pathway_seqs$col_name)) %>%
        column_to_rownames("study_id") %>%
        t() %>% as.matrix()
      
      # Map rownames from seq_IDs to gene symbols
      seq_to_gene <- pathway_seqs %>%
        dplyr::select(col_name, gene_sym) %>% deframe()
      rownames(soma_mat) <- seq_to_gene[rownames(soma_mat)]
      
      # Align to metadata
      shared <- intersect(meta_soma$study_id, colnames(soma_mat))
      soma_mat   <- soma_mat[, shared, drop = FALSE]
      meta_soma  <- meta_soma %>% dplyr::filter(study_id %in% shared) %>%
        arrange(match(study_id, shared))
      
      cat("Final SomaScan matrix:", nrow(soma_mat), "genes x",
          ncol(soma_mat), "samples\n")
      
      # Run limma
      prot_c <- run_limma_on_mat(soma_mat, meta_soma, rownames(soma_mat))
      
      if (!is.null(prot_c) && nrow(prot_c) > 0) {
        combined_c <- make_combined(prot_c, nebula_res, "Blood SomaScan")
        results_list[["blood_soma"]] <- combined_c
        print_concordance(combined_c, "Blood SomaScan proteomics")
        ggsave(file.path(output_dir, "blood_soma_scatter.pdf"),
               make_scatter(combined_c, "Blood SomaScan vs NEBULA",
                            "SomaScan blood proteomics"), width=9, height=7)
        ggsave(file.path(output_dir, "blood_soma_bar.pdf"),
               make_bar(combined_c, "Blood SomaScan vs NEBULA",
                        "SomaScan blood proteomics"), width=18, height=7)
      }
    } else {
      cat("None of the pathway genes found in SomaScan annotation.\n")
      cat("These genes may not be measured in the SomaScan panel:\n")
      print(all_pathway_genes)
    }
  }
} else {
  cat("Harmonized SomaScan file not found.\n")
  cat("Tried:\n  ", soma_file, "\n")
  cat("Update soma_file path in PART 2 above.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 3: SUMMARY - All analyses side by side
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n════ SUMMARY: All concordance results ════\n")

if (length(results_list) > 0) {
  summary_df <- bind_rows(results_list) %>%
    group_by(analysis, pathway) %>%
    summarise(n_genes    = n(),
              concordant = sum(concordant, na.rm=TRUE),
              pct        = round(100*concordant/n_genes, 1),
              .groups    = "drop")
  print(summary_df)
  
  # Multi-panel scatter: one panel per analysis
  scatter_plots <- lapply(names(results_list), function(nm) {
    df    <- results_list[[nm]]
    label <- unique(df$analysis)
    make_scatter(df, label)
  })
  
  if (length(scatter_plots) >= 2) {
    p_all <- wrap_plots(scatter_plots, ncol = min(length(scatter_plots), 3)) +
      plot_annotation(
        title    = "Multi-Omics Concordance with scRNA NEBULA",
        subtitle = "Rinschen kidney proteomics + SomaScan blood proteomics vs Schaub GEO NEBULA",
        theme    = theme(plot.title    = element_text(face="bold", size=14),
                         plot.subtitle = element_text(size=10))
      )
    ggsave(file.path(output_dir, "multiomics_concordance_panel.pdf"),
           p_all, width = 9 * min(length(scatter_plots), 3), height = 7)
    cat("Saved: multiomics_concordance_panel.pdf\n")
  }
  
  # Excel output
  wb <- createWorkbook()
  for (nm in names(results_list)) {
    sheet <- substr(nm, 1, 31)  # Excel sheet name limit
    addWorksheet(wb, sheet); writeData(wb, sheet, results_list[[nm]])
  }
  addWorksheet(wb, "Concordance_summary"); writeData(wb, "Concordance_summary", summary_df)
  saveWorkbook(wb, file.path(output_dir, "multiomics_scrna_comparison.xlsx"), overwrite=TRUE)
  cat("Saved: multiomics_scrna_comparison.xlsx\n")
}

cat("\n=== Done ===\n")
cat("Output directory:", output_dir, "\n")