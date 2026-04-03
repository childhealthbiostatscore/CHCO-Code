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
# 8. GLUCOSE CLEARANCE CORRELATIONS (uncomment when available)
# ============================================================

# clin$glucose_clearance <- c(...)  # add per-participant values from Daniel
#
# run_correlation <- function(mat, gc, annot_df, filename) {
#   apply(mat, 1, function(vals) {
#     t <- cor.test(vals, gc, method = "spearman")
#     c(rho = unname(t$estimate), p = t$p.value)
#   }) %>%
#     t() %>% as.data.frame() %>%
#     rownames_to_column("Feature_ID") %>%
#     mutate(p_adj = p.adjust(p, method = "BH")) %>%
#     left_join(annot_df, by = "Feature_ID") %>%
#     arrange(p_adj) %>%
#     write.csv(filename, row.names = FALSE)
# }
#
# run_correlation(art_mat, clin$glucose_clearance, annot,
#                 paste0(path_out, "correlation_gc_Arterial.csv"))
# run_correlation(ven_mat, clin$glucose_clearance, annot,
#                 paste0(path_out, "correlation_gc_Venous.csv"))
# run_correlation(av_mat,  clin$glucose_clearance, annot,
#                 paste0(path_out, "correlation_gc_AV_Difference.csv"))