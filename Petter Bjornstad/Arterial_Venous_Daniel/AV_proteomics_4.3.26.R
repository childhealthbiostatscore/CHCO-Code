###AV Analysis

library(tidyverse)
library(readxl)
library(haven)
library(ggrepel)
library(pheatmap)
library(limma)

# ============================================================
# PATHS
# ============================================================

path_prot  <- "../../Downloads/T2D.v.HC/T2D.v.HC_Differential_Expression.xlsx"
path_clin  <- "../../Downloads/ROCKIES_Hoofddatabase__jan2026.sav"
path_out   <- "Projects/Art_Ven_Daniel/"
dir.create(path_out, showWarnings = FALSE)


# ============================================================
# 1. LOAD PROTEOMICS + EXTRACT PARTICIPANT IDs
# ============================================================

de_raw   <- read_xlsx(path_prot)
get_id   <- function(x) sub("sample\\.(ROCK\\d+)_.*", "\\1", x)
art_cols <- grep("_Arterial$", colnames(de_raw), value = TRUE)
ven_cols <- grep("_Venous$",   colnames(de_raw), value = TRUE)
shared   <- intersect(get_id(art_cols), get_id(ven_cols))


# ============================================================
# 2. LOAD + ALIGN CLINICAL DATA
# ============================================================

clin_raw <- read_sav(path_clin)

clin <- clin_raw %>%
  transmute(
    participant_id = as.character(Participant),
    group          = ifelse(as.numeric(Group) == 1, "T2D", "NC")
  )

common_ids <- intersect(shared, clin$participant_id)
cat("Proteomics N:", length(shared),
    "| Clinical N:", nrow(clin),
    "| Matched N:", length(common_ids), "\n")
cat("Missing from clinical:", paste(setdiff(shared, clin$participant_id), collapse = ", "), "\n")

shared <- common_ids
clin   <- clin %>%
  filter(participant_id %in% common_ids) %>%
  arrange(match(participant_id, common_ids))

stopifnot(all(clin$participant_id == shared))
cat("Group breakdown:\n"); print(table(clin$group))


# ============================================================
# 3. BUILD MATRICES (subset to matched participants)
# ============================================================

annot <- de_raw %>%
  select(Feature_ID, Target, TargetFullName, UniProt)

art_mat <- as.matrix(de_raw[, paste0("sample.", shared, "_Arterial")])
ven_mat <- as.matrix(de_raw[, paste0("sample.", shared, "_Venous")])
av_mat  <- art_mat - ven_mat

rownames(art_mat) <- de_raw$Feature_ID
rownames(ven_mat) <- de_raw$Feature_ID
rownames(av_mat)  <- de_raw$Feature_ID
colnames(art_mat) <- colnames(ven_mat) <- colnames(av_mat) <- shared


# ============================================================
# 4. LIMMA FUNCTION (reusable for all three analyses)
# ============================================================

# Significance threshold — FDR < 0.05 is ideal but may be too strict
# with N=9/group. Adjust here if needed: try 0.10 or 0.20
# Also use nominal p < 0.05 as fallback for exploratory results
fdr_thresh <- 0.05
p_thresh   <- 0.05   # nominal, unadjusted fallback

run_limma <- function(mat, group_vec, annot_df, label) {
  design <- model.matrix(~ factor(group_vec, levels = c("NC", "T2D")))
  fit    <- eBayes(lmFit(mat, design))
  res    <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("Feature_ID") %>%
    left_join(annot_df, by = "Feature_ID") %>%
    mutate(
      analysis  = label,
      sig       = adj.P.Val < fdr_thresh,
      sig_nom   = P.Value   < p_thresh,    # nominal p<0.05 fallback
      direction = case_when(
        sig_nom & logFC > 0 ~ "Higher in T2D",
        sig_nom & logFC < 0 ~ "Lower in T2D",
        TRUE                ~ "ns"
      )
    )
  cat("\n---", label, "---\n")
  cat("Significant FDR<", fdr_thresh, ":", sum(res$sig), "\n")
  cat("Nominal p<0.05:    ", sum(res$sig_nom), "\n")
  res
}

res_art <- run_limma(art_mat, clin$group, annot, "Arterial")
res_ven <- run_limma(ven_mat, clin$group, annot, "Venous")
res_av  <- run_limma(av_mat,  clin$group, annot, "AV Difference")

# Diagnostics — check p-value distributions
cat("\n--- P-value diagnostics ---\n")
for (res in list(res_art, res_ven, res_av)) {
  cat(res$analysis[1], "\n")
  cat("  Min raw p-value:     ", min(res$P.Value, na.rm = TRUE), "\n")
  cat("  Min adj p-value:     ", min(res$adj.P.Val, na.rm = TRUE), "\n")
  cat("  Hits at p<0.05:      ", sum(res$P.Value < 0.05, na.rm = TRUE), "\n")
  cat("  Hits at FDR<0.20:    ", sum(res$adj.P.Val < 0.20, na.rm = TRUE), "\n")
  cat("  Hits at FDR<0.05:    ", sum(res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
}


# ============================================================
# 5. VOLCANO PLOT FUNCTION (reusable)
# ============================================================

make_volcano <- function(res, title_label) {
  # Use FDR if anything passes, otherwise fall back to nominal p<0.05
  use_sig <- if (any(res$sig)) "sig" else "sig_nom"
  thresh_label <- if (use_sig == "sig") paste("FDR <", fdr_thresh) else "nominal p < 0.05"
  
  top_ids <- res %>% filter(.data[[use_sig]]) %>% slice_min(P.Value, n = 20) %>% pull(Feature_ID)
  res <- res %>% mutate(label = ifelse(Feature_ID %in% top_ids, Target, NA))
  
  ggplot(res, aes(x = logFC, y = -log10(P.Value), color = direction)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 20,
                    na.rm = TRUE, box.padding = 0.4, segment.color = "grey60") +
    scale_color_manual(values = c(
      "Higher in T2D" = "firebrick",
      "Lower in T2D"  = "steelblue",
      "ns"            = "grey70"
    )) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(
      title    = paste("T2D vs NC —", title_label, "(T=90)"),
      subtitle = paste0(sum(res[[use_sig]]), " proteins at ", thresh_label),
      x        = "logFC (T2D − NC)", y = "-log10(p-value)", color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

v_art <- make_volcano(res_art, "Arterial")
v_ven <- make_volcano(res_ven, "Venous")
v_av  <- make_volcano(res_av,  "AV Difference")

ggsave(paste0(path_out, "volcano_Arterial.pdf"),     v_art, width = 8, height = 7)
ggsave(paste0(path_out, "volcano_Venous.pdf"),       v_ven, width = 8, height = 7)
ggsave(paste0(path_out, "volcano_AV_Difference.pdf"),v_av,  width = 8, height = 7)


# ============================================================
# 6. HEATMAP FUNCTION (reusable)
# ============================================================

make_heatmap <- function(mat, res, clin_df, title_label, filename) {
  use_sig  <- if (any(res$sig)) "sig" else "sig_nom"
  top_prots <- res %>% filter(.data[[use_sig]]) %>% slice_min(P.Value, n = 50) %>% pull(Feature_ID)
  if (length(top_prots) == 0) { cat("No significant proteins for heatmap:", title_label, "\n"); return(NULL) }
  
  prot_labels    <- res %>% select(Feature_ID, Target) %>% deframe()
  m              <- mat[top_prots, ]
  rownames(m)    <- prot_labels[rownames(m)]
  
  ann_col        <- data.frame(Group = clin_df$group, row.names = clin_df$participant_id)
  ann_colors     <- list(Group = c(NC = "steelblue", T2D = "firebrick"))
  
  pheatmap(m,
           annotation_col    = ann_col,
           annotation_colors = ann_colors,
           scale             = "row",
           show_colnames     = FALSE,
           fontsize_row      = 7,
           clustering_method = "ward.D2",
           main              = paste("Top 50 Significant Proteins —", title_label),
           filename          = filename,
           width = 8, height = 10
  )
}

make_heatmap(art_mat, res_art, clin, "Arterial",     paste0(path_out, "heatmap_Arterial.pdf"))
make_heatmap(ven_mat, res_ven, clin, "Venous",       paste0(path_out, "heatmap_Venous.pdf"))
make_heatmap(av_mat,  res_av,  clin, "AV Difference",paste0(path_out, "heatmap_AV_Difference.pdf"))


# ============================================================
# 7. SAVE RESULTS TABLES
# ============================================================

save_results <- function(res, filename) {
  use_sig <- if (any(res$sig)) "sig" else "sig_nom"
  res %>%
    filter(.data[[use_sig]]) %>%
    arrange(adj.P.Val) %>%
    select(Feature_ID, Target, TargetFullName, UniProt,
           logFC, P.Value, adj.P.Val, direction) %>%
    write.csv(filename, row.names = FALSE)
}

save_results(res_art, paste0(path_out, "significant_Arterial.csv"))
save_results(res_ven, paste0(path_out, "significant_Venous.csv"))
save_results(res_av,  paste0(path_out, "significant_AV_Difference.csv"))

cat("\nAll outputs saved to:", path_out, "\n")

# ============================================================
# PATHS (add to your existing PATHS block at the top)
# ============================================================

path_omics <- "../../Downloads/Omics_subgroup_ROCKIES_OGTT_Database_final_june24_DVR.sav"


# ============================================================
# 8. LOAD T=90 CLINICAL VARIABLES (FEgluc + Feins)
# ============================================================

omics_raw <- read_sav(path_omics)

# ── Check available variable names ──────────────────────────
# Run this once to confirm exact column names in your .sav:
# names(omics_raw)

omics_t90 <- omics_raw %>%
  transmute(
    participant_id = as.character(Participant),   # adjust if ID column differs
    FEgluc_T90     = FEgluc,                      # adjust to exact column name
    Feins_T90      = Feins                         # adjust to exact column name
  ) %>%
  filter(!is.na(FEgluc_T90) | !is.na(Feins_T90))

# Merge with matched participants from main analysis
clin_t90 <- clin %>%
  left_join(omics_t90, by = "participant_id")

cat("\nT=90 variable availability:\n")
cat("  FEgluc available N:", sum(!is.na(clin_t90$FEgluc_T90)), "\n")
cat("  Feins  available N:", sum(!is.na(clin_t90$Feins_T90)),  "\n")


# ============================================================
# 9. SPEARMAN CORRELATION FUNCTION (proteins ~ FEgluc / Feins)
# ============================================================

run_correlation <- function(mat, outcome_vec, outcome_name, annot_df, label) {
  
  # Only use participants with non-missing outcome
  keep      <- !is.na(outcome_vec)
  mat_sub   <- mat[, keep]
  out_sub   <- outcome_vec[keep]
  cat("\n---", label, "~", outcome_name, "| N =", sum(keep), "---\n")
  
  res_cor <- apply(mat_sub, 1, function(vals) {
    tt <- cor.test(vals, out_sub, method = "spearman", exact = FALSE)
    c(rho = unname(tt$estimate), p = tt$p.value)
  }) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Feature_ID") %>%
    mutate(
      p_adj     = p.adjust(p, method = "BH"),
      sig       = p_adj < fdr_thresh,
      sig_nom   = p     < p_thresh,
      direction = case_when(
        sig_nom & rho > 0 ~ "Positive",
        sig_nom & rho < 0 ~ "Negative",
        TRUE              ~ "ns"
      ),
      outcome   = outcome_name,
      analysis  = label
    ) %>%
    left_join(annot_df, by = "Feature_ID") %>%
    arrange(p_adj)
  
  cat("  Significant FDR <", fdr_thresh, ":", sum(res_cor$sig),    "\n")
  cat("  Nominal p < 0.05: ", sum(res_cor$sig_nom), "\n")
  res_cor
}

# Run for all matrix × outcome combinations
cor_results <- list(
  run_correlation(art_mat, clin_t90$FEgluc_T90, "FEgluc", annot, "Arterial"),
  run_correlation(ven_mat, clin_t90$FEgluc_T90, "FEgluc", annot, "Venous"),
  run_correlation(av_mat,  clin_t90$FEgluc_T90, "FEgluc", annot, "AV Difference"),
  run_correlation(art_mat, clin_t90$Feins_T90,  "Feins",  annot, "Arterial"),
  run_correlation(ven_mat, clin_t90$Feins_T90,  "Feins",  annot, "Venous"),
  run_correlation(av_mat,  clin_t90$Feins_T90,  "Feins",  annot, "AV Difference")
)


# ============================================================
# 10. CORRELATION VOLCANO PLOTS (rho vs -log10 p)
# ============================================================

make_cor_volcano <- function(res, mat_label, outcome_label) {
  use_sig     <- if (any(res$sig)) "sig" else "sig_nom"
  thresh_label <- if (use_sig == "sig") paste("FDR <", fdr_thresh) else "nominal p < 0.05"
  
  top_ids <- res %>% filter(.data[[use_sig]]) %>% slice_min(p, n = 20) %>% pull(Feature_ID)
  res     <- res %>% mutate(label = ifelse(Feature_ID %in% top_ids, Target, NA))
  
  ggplot(res, aes(x = rho, y = -log10(p), color = direction)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 20,
                    na.rm = TRUE, box.padding = 0.4, segment.color = "grey60") +
    scale_color_manual(values = c(
      "Positive" = "firebrick",
      "Negative" = "steelblue",
      "ns"       = "grey70"
    )) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(
      title    = paste(mat_label, "proteins ~", outcome_label, "(T=90)"),
      subtitle = paste0(sum(res[[use_sig]]), " proteins at ", thresh_label),
      x        = "Spearman ρ", y = "-log10(p-value)", color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

# Generate and save all 6 volcano plots
cor_plot_params <- list(
  list("Arterial",     "FEgluc"), list("Venous",       "FEgluc"), list("AV Difference", "FEgluc"),
  list("Arterial",     "Feins"),  list("Venous",       "Feins"),  list("AV Difference", "Feins")
)

for (i in seq_along(cor_results)) {
  p     <- cor_plot_params[[i]]
  fname <- paste0(path_out, "volcano_cor_", gsub(" ", "_", p[[1]]), "_", p[[2]], ".pdf")
  ggsave(fname, make_cor_volcano(cor_results[[i]], p[[1]], p[[2]]), width = 8, height = 7)
}


# ============================================================
# 11. SAVE CORRELATION RESULTS TABLES
# ============================================================

for (i in seq_along(cor_results)) {
  res     <- cor_results[[i]]
  use_sig <- if (any(res$sig)) "sig" else "sig_nom"
  p       <- cor_plot_params[[i]]
  fname   <- paste0(path_out, "correlation_", gsub(" ", "_", p[[1]]), "_", p[[2]], ".csv")
  
  res %>%
    filter(.data[[use_sig]]) %>%
    arrange(p_adj) %>%
    select(Feature_ID, Target, TargetFullName, UniProt,
           rho, p, p_adj, direction, outcome, analysis) %>%
    write.csv(fname, row.names = FALSE)
}


# ============================================================
# 12. PATHWAY ENRICHMENT (enrichR) ON SIGNIFICANT PROTEINS
# ============================================================

# install.packages("enrichR")   # run once if needed
library(enrichR)

# Databases to query — covers GO, KEGG, Reactome, WikiPathways
enrich_dbs <- c(
  "GO_Biological_Process_2023",
  "KEGG_2021_Human",
  "Reactome_2022",
  "WikiPathway_2023_Human"
)

run_enrichment <- function(res, mat_label, outcome_label) {
  use_sig  <- if (any(res$sig)) "sig" else "sig_nom"
  sig_genes <- res %>% filter(.data[[use_sig]]) %>% pull(Target)
  
  if (length(sig_genes) < 3) {
    cat("Too few significant proteins for enrichment:", mat_label, outcome_label, "\n")
    return(NULL)
  }
  cat("\nRunning enrichment for", mat_label, "~", outcome_label,
      "| N genes =", length(sig_genes), "\n")
  
  enrich_res <- enrichr(sig_genes, enrich_dbs)
  
  # Combine databases, keep significant terms
  combined <- bind_rows(enrich_res, .id = "database") %>%
    filter(Adjusted.P.value < 0.05) %>%
    arrange(Adjusted.P.value) %>%
    mutate(mat = mat_label, outcome = outcome_label)
  
  cat("  Enriched terms (adj.p < 0.05):", nrow(combined), "\n")
  
  # Save table
  fname <- paste0(path_out, "enrichment_", gsub(" ", "_", mat_label), "_", outcome_label, ".csv")
  write.csv(combined, fname, row.names = FALSE)
  
  # Dot plot of top 15 terms per database
  if (nrow(combined) > 0) {
    plot_data <- combined %>%
      group_by(database) %>%
      slice_min(Adjusted.P.value, n = 15) %>%
      ungroup() %>%
      mutate(
        Term       = str_trunc(Term, 50),
        Gene_count = as.integer(str_extract(Overlap, "^\\d+"))
      )
    
    p <- ggplot(plot_data, aes(x = -log10(Adjusted.P.value), y = reorder(Term, -Adjusted.P.value),
                               size = Gene_count, color = database)) +
      geom_point(alpha = 0.8) +
      labs(
        title    = paste("Pathway Enrichment —", mat_label, "~", outcome_label),
        subtitle = "Top 15 terms per database, adj.p < 0.05",
        x        = "-log10(adjusted p-value)", y = NULL,
        size     = "Gene count", color = "Database"
      ) +
      theme_bw(base_size = 10) +
      theme(axis.text.y = element_text(size = 7), legend.position = "right")
    
    ggsave(paste0(path_out, "enrichment_dotplot_", gsub(" ", "_", mat_label), "_", outcome_label, ".pdf"),
           p, width = 10, height = 8)
  }
  combined
}

# Run enrichment for all 6 result sets
enrich_results <- mapply(
  function(res, p) run_enrichment(res, p[[1]], p[[2]]),
  cor_results, cor_plot_params,
  SIMPLIFY = FALSE
)

cat("\nAll correlation and pathway outputs saved to:", path_out, "\n")


