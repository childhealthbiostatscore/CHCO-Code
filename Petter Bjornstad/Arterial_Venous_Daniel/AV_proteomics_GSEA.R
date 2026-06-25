# ============================================================
# ROCKIES — Pathway Enrichment for Correlation Results
# 1. enrichR (split by direction: higher / lower)
# 2. fgsea (ranked by rho — directional NES plots)
# ============================================================

library(tidyverse)
library(enrichR)
library(biomaRt)
library(fgsea)
library(msigdbr)
library(ggplot2)

path_results <- "Projects/Art_Ven_Daniel/"
path_out     <- "Projects/Art_Ven_Daniel/Pathways/"
dir.create(path_out, showWarnings = FALSE, recursive = TRUE)

# enrichR databases
dbs <- c(
  "GO_Biological_Process_2025",
  "KEGG_2026",
  "Reactome_Pathways_2024",
  "WikiPathways_2024_Human"
)


# ============================================================
# 1. UNIPROT → HGNC MAP
# ============================================================

cat("Building UniProt → HGNC map via biomaRt...\n")

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

all_files <- list.files(path_results, pattern = "^correlation_.*\\.csv$", full.names = TRUE)

uniprot_vec <- map_dfr(all_files, read_csv, show_col_types = FALSE) %>%
  distinct(UniProt) %>%
  filter(!is.na(UniProt), UniProt != "") %>%
  pull(UniProt) %>%
  str_split("\\|") %>%
  unlist() %>%
  unique() %>%
  .[. != ""]

up_map <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters    = "uniprotswissprot",
  values     = uniprot_vec,
  mart       = mart
) %>%
  filter(hgnc_symbol != "") %>%
  rename(UniProt_single = uniprotswissprot, gene_symbol = hgnc_symbol)

cat("Mapped", nrow(up_map), "UniProt → HGNC entries\n\n")


# ============================================================
# 2. GENE SETS FOR FGSEA
# ============================================================

cat("Loading gene sets for fgsea...\n")

make_pathways <- function(collection, subcollection = NULL) {
  args <- list(species = "Homo sapiens", collection = collection)
  if (!is.null(subcollection)) args$subcollection <- subcollection
  db <- do.call(msigdbr, args)
  db %>%
    select(gs_name, gene_symbol) %>%
    group_by(gs_name) %>%
    summarise(genes = list(gene_symbol), .groups = "drop") %>%
    deframe()
}

pathway_collections <- list(
  GO_BP        = make_pathways("C5", "GO:BP"),
  KEGG         = make_pathways("C2", "CP:KEGG_LEGACY"),
  Reactome     = make_pathways("C2", "CP:REACTOME"),
  WikiPathways = make_pathways("C2", "CP:WIKIPATHWAYS")
)

cat("Gene set sizes:\n")
walk2(pathway_collections, names(pathway_collections),
      ~cat(" ", .y, ":", length(.x), "pathways\n"))
cat("\n")


# ============================================================
# 3. HELPERS
# ============================================================

# HGNC symbols for a directional subset
get_symbols <- function(res_df) {
  res_df %>%
    mutate(UniProt_single = str_split(UniProt, "\\|")) %>%
    unnest(UniProt_single) %>%
    left_join(up_map, by = "UniProt_single") %>%
    filter(!is.na(gene_symbol)) %>%
    distinct(Feature_ID, gene_symbol) %>%
    pull(gene_symbol) %>%
    unique()
}

# Ranked named vector (rho) for fgsea — all proteins
get_ranked <- function(res_df) {
  res_df %>%
    mutate(UniProt_single = str_split(UniProt, "\\|")) %>%
    unnest(UniProt_single) %>%
    left_join(up_map, by = "UniProt_single") %>%
    filter(!is.na(gene_symbol)) %>%
    group_by(gene_symbol) %>%
    slice_max(abs(rho), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(rho)) %>%
    select(gene_symbol, rho) %>%
    deframe()
}


# ── enrichR (one direction) ──────────────────────────────────
run_enrichr_dir <- function(genes, label, filename_prefix) {
  cat("    enrichR:", label, "| N =", length(genes), "\n")
  if (length(genes) < 5) { cat("      Skipping\n"); return(NULL) }
  
  enr <- enrichr(genes, dbs)
  
  results_all <- imap_dfr(enr, function(df, db_name) {
    if (nrow(df) == 0) return(NULL)
    df %>%
      mutate(Database = db_name) %>%
      filter(Adjusted.P.value < 0.05) %>%
      arrange(Adjusted.P.value) %>%
      slice_head(n = 15)
  })
  
  write_csv(results_all, paste0(path_out, filename_prefix, "_enrichr.csv"))
  if (nrow(results_all) == 0) { cat("      No significant terms\n"); return(NULL) }
  
  plot_df <- results_all %>%
    mutate(
      GeneRatio = map_dbl(Overlap, function(x) {
        p <- str_split(x, "/")[[1]]
        as.numeric(p[1]) / as.numeric(p[2])
      }),
      Term = str_trunc(Term, 55),
      Term = fct_reorder(Term, -log10(Adjusted.P.value))
    ) %>%
    group_by(Database) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  p <- ggplot(plot_df, aes(x = GeneRatio, y = Term,
                           color = -log10(Adjusted.P.value),
                           size = Combined.Score)) +
    geom_point() +
    scale_color_gradient(low = "steelblue", high = "firebrick",
                         name = "-log10\n(adj. p)") +
    scale_size_continuous(name = "Combined\nScore", range = c(2, 8)) +
    facet_wrap(~ Database, scales = "free_y", ncol = 1) +
    labs(title = label, x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey90"),
          axis.text.y = element_text(size = 8))
  
  plot_h <- max(6, 1.5 + nrow(plot_df) * 0.28 + n_distinct(plot_df$Database) * 0.5)
  ggsave(paste0(path_out, filename_prefix, "_dotplot.pdf"),
         p, width = 10, height = plot_h, limitsize = FALSE)
  
  cat("      Saved:", nrow(results_all), "terms\n")
  return(results_all)
}


# ── fgsea ────────────────────────────────────────────────────
run_fgsea <- function(ranked_vec, label, filename_prefix) {
  cat("    fgsea:", label, "| N proteins =", length(ranked_vec), "\n")
  
  all_gsea <- imap_dfr(pathway_collections, function(pathways, db_name) {
    fgsea(
      pathways    = pathways,
      stats       = ranked_vec,
      minSize     = 10,
      maxSize     = 500,
      nPermSimple = 10000
    ) %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      mutate(Database = db_name)
  })
  
  # Save CSV
  all_gsea %>%
    mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) %>%
    write_csv(paste0(path_out, filename_prefix, "_fgsea.csv"))
  
  if (nrow(all_gsea) == 0) {
    cat("      No significant fgsea terms\n")
    return(NULL)
  }
  
  # NES plot
  plot_df <- all_gsea %>%
    group_by(Database) %>%
    slice_max(abs(NES), n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      pathway   = str_trunc(pathway, 55),
      pathway   = fct_reorder(pathway, NES),
      Direction = ifelse(NES > 0,
                         "Positive (↑ with outcome)",
                         "Negative (↓ with outcome)")
    )
  
  p <- ggplot(plot_df, aes(x = NES, y = pathway,
                           fill = Direction,
                           size = -log10(padj))) +
    geom_point(shape = 21, color = "grey30") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(
      values = c("Positive (↑ with outcome)" = "#d73027",
                 "Negative (↓ with outcome)" = "#4575b4"),
      name = NULL
    ) +
    scale_size_continuous(name = "-log10\n(adj. p)", range = c(2, 9)) +
    facet_wrap(~ Database, scales = "free_y", ncol = 1) +
    labs(title = paste("fgsea —", label),
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey90"),
          axis.text.y = element_text(size = 8),
          legend.position = "bottom")
  
  plot_h <- max(6, 1.5 + nrow(plot_df) * 0.28 + n_distinct(plot_df$Database) * 0.5)
  ggsave(paste0(path_out, filename_prefix, "_fgsea_NES.pdf"),
         p, width = 10, height = plot_h, limitsize = FALSE)
  
  cat("      fgsea saved:", nrow(all_gsea), "significant pathways\n")
  return(all_gsea)
}


# ============================================================
# 4. MAIN LOOP
# ============================================================

fdr_file <- "correlation_Arterial_FEgluc.csv"

for (f in all_files) {
  
  base_name  <- tools::file_path_sans_ext(basename(f))
  label      <- str_replace_all(str_remove(base_name, "^correlation_"), "_", " ")
  out_prefix <- str_remove(base_name, "^correlation_")
  res        <- read_csv(f, show_col_types = FALSE)
  
  cat("\n══════════════════════════════════════════\n")
  cat("", label, "\n")
  cat("══════════════════════════════════════════\n")
  
  # ── enrichR: split by direction ───────────────────────────
  if (basename(f) == fdr_file) {
    higher <- res %>% filter(p_adj < 0.05, direction == "Positive")
    lower  <- res %>% filter(p_adj < 0.05, direction == "Negative")
  } else {
    higher <- res %>% filter(p < 0.05, direction == "Positive")
    lower  <- res %>% filter(p < 0.05, direction == "Negative")
  }
  
  run_enrichr_dir(get_symbols(higher),
                  paste0("Higher with ", label),
                  paste0(out_prefix, "_higher"))
  
  run_enrichr_dir(get_symbols(lower),
                  paste0("Lower with ", label),
                  paste0(out_prefix, "_lower"))
  
  # ── fgsea: ranked by rho ──────────────────────────────────
  run_fgsea(get_ranked(res), label, out_prefix)
}

cat("\n✓ All done. Outputs saved to:", path_out, "\n")