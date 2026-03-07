library(Seurat)
library(nebula)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tidyr)
library(purrr)
library(fgsea)
library(msigdbr)
library(stringr)
library(gtsummary)

# ══════════════════════════════════════════════════════════════════════════════
# 0. LOAD & PREPARE SEURAT OBJECT
# ══════════════════════════════════════════════════════════════════════════════

load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')
so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

# ── Load & collapse harmonized data ─────────────────────────────────────────
harmonized_data <- read.csv(
  "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv",
  na = ''
)

dat <- harmonized_data %>%
  dplyr::select(-dob) %>%
  arrange(date_of_screen) %>%
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

# ── Merge clinical variables into Seurat metadata via match() ───────────────
vars_to_merge <- c("baseline_ffa", "ffa_suppression", "age", "bmi",
                   "acr_u", "hba1c",
                   "epic_insulin_1", "epic_sglti2_1", "epic_glp1ra_1",
                   "epic_tzd_1", "epic_mfm_1", "epic_raasi_1",
                   "epic_statin_1", "epic_fibrate_1")

merge_df <- dat %>%
  dplyr::select(any_of(c("record_id", vars_to_merge))) %>%
  distinct(record_id, .keep_all = TRUE)

for (v in vars_to_merge) {
  if (v %in% colnames(merge_df)) {
    so_subset@meta.data[[v]] <- merge_df[[v]][
      match(so_subset@meta.data$record_id, merge_df$record_id)
    ]
  }
}

cat("baseline_ffa non-NA cells:", sum(!is.na(so_subset@meta.data$baseline_ffa)), "\n")
cat("ffa_suppression non-NA cells:", sum(!is.na(so_subset@meta.data$ffa_suppression)), "\n")

# ── Cell type annotations ────────────────────────────────────────────────────
so_subset$celltype1 <- case_when(
  grepl("PT-",  so_subset$celltype_rpca) ~ "PT",
  grepl("TAL-", so_subset$celltype_rpca) ~ "TAL",
  grepl("EC-",  so_subset$celltype_rpca) ~ "EC",
  grepl("POD",  so_subset$celltype_rpca) ~ "POD",
  grepl("MAC",  so_subset$celltype_rpca) ~ "MAC",
  grepl("MON",  so_subset$celltype_rpca) ~ "MON",
  grepl("PC-",  so_subset$celltype_rpca) ~ "PC",
  grepl("FIB",  so_subset$celltype_rpca) ~ "FIB_MC_VSMC",
  grepl("DTL",  so_subset$celltype_rpca) ~ "DTL",
  so_subset$celltype_rpca == "DCT"        ~ "DCT",
  so_subset$celltype_rpca == "ATL"        ~ "ATL",
  so_subset$celltype_rpca == "B"          ~ "B",
  so_subset$celltype_rpca == "T"          ~ "T"
)
so_subset$celltype1      <- as.character(so_subset$celltype1)
so_subset$KPMP_celltype2 <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2      <- ifelse(
  so_subset$KPMP_celltype %in% c("aPT", "PT-S1/S2", "PT-S3"), "PT",
  ifelse(grepl("TAL", so_subset$KPMP_celltype), "TAL",
         ifelse(grepl("EC-", so_subset$KPMP_celltype), "EC",
                so_subset$KPMP_celltype2)))
so_subset$DCT_celltype <- ifelse(
  so_subset$KPMP_celltype %in% c("DCT", "dDCT"), "DCT", "Non-DCT")

# ── Derive microalbuminuria from acr_u ───────────────────────────────────────
# Standard thresholds: ≥30 mg/g = microalbuminuria, ≥300 mg/g = macroalbuminuria
so_subset@meta.data <- so_subset@meta.data %>%
  mutate(microalbuminuria = case_when(
    is.na(acr_u)    ~ NA_character_,
    acr_u >= 300    ~ "Macroalbuminuria",
    acr_u >= 30     ~ "Microalbuminuria",
    TRUE            ~ "Normoalbuminuria"
  ))

# ── Subset to T2D, exclude CRC-55 ───────────────────────────────────────────
so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type_2_Diabetes')

dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/NEFA_NEBULA/'
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. DEMOGRAPHICS TABLE
# ══════════════════════════════════════════════════════════════════════════════

demo_df <- so_subset@meta.data %>%
  distinct(record_id, .keep_all = TRUE) %>%
  filter(!is.na(baseline_ffa)) %>%
  mutate(
    across(c(age, bmi, hba1c, baseline_ffa, ffa_suppression, acr_u),
           ~ as.numeric(as.character(.x))),
    across(c(epic_insulin_1, epic_sglti2_1, epic_glp1ra_1, epic_tzd_1,
             epic_mfm_1, epic_raasi_1, epic_statin_1, epic_fibrate_1,
             microalbuminuria), as.factor)
  )

demo_table <- demo_df %>%
  dplyr::select(age, sex, bmi, hba1c,
                baseline_ffa, ffa_suppression,
                acr_u, microalbuminuria,
                epic_insulin_1, epic_sglti2_1, epic_glp1ra_1,
                epic_tzd_1, epic_mfm_1, epic_raasi_1,
                epic_statin_1, epic_fibrate_1) %>%
  tbl_summary(
    type = list(
      age              ~ "continuous",
      bmi              ~ "continuous",
      hba1c            ~ "continuous",
      baseline_ffa     ~ "continuous",
      ffa_suppression  ~ "continuous",
      acr_u            ~ "continuous2",   # continuous2 allows custom statistic per var
      microalbuminuria ~ "categorical",
      epic_insulin_1   ~ "categorical",
      epic_sglti2_1    ~ "categorical",
      epic_glp1ra_1    ~ "categorical",
      epic_tzd_1       ~ "categorical",
      epic_mfm_1       ~ "categorical",
      epic_raasi_1     ~ "categorical",
      epic_statin_1    ~ "categorical",
      epic_fibrate_1   ~ "categorical"
    ),
    statistic = list(
      age             ~ "{mean} ({sd})",
      bmi             ~ "{mean} ({sd})",
      hba1c           ~ "{mean} ({sd})",
      baseline_ffa    ~ "{mean} ({sd})",
      ffa_suppression ~ "{mean} ({sd})",
      acr_u           ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age             ~ 1,
      bmi             ~ 1,
      hba1c           ~ 1,
      baseline_ffa    ~ 1,
      ffa_suppression ~ 1,
      acr_u           ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age              ~ "Age, years",
      sex              ~ "Sex",
      bmi              ~ "BMI, kg/m²",
      hba1c            ~ "HbA1c, %",
      baseline_ffa     ~ "Baseline FFA, µmol/L",
      ffa_suppression  ~ "FFA suppression, %",
      acr_u            ~ "Urine albumin, mg/L",
      microalbuminuria ~ "Albuminuria category",
      epic_insulin_1   ~ "Insulin",
      epic_sglti2_1    ~ "SGLT2i",
      epic_glp1ra_1    ~ "GLP-1RA",
      epic_tzd_1       ~ "PPARγ agonist (TZD)",
      epic_mfm_1       ~ "Metformin",
      epic_raasi_1     ~ "ACEi/ARB",
      epic_statin_1    ~ "Statin",
      epic_fibrate_1   ~ "Fenofibrate"
    ),
    missing_text = "Missing"
  ) %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) unless otherwise noted; urine albumin shown as median (IQR)")

print(demo_table)

# Save as Word document
library(flextable)
library(officer)

demo_table %>%
  as_flex_table() %>%
  flextable::save_as_docx(
    path = file.path(dir.results, "demographics_table.docx")
  )

# Also save as CSV as backup
demo_table %>%
  as_tibble() %>%
  write.csv(file.path(dir.results, "demographics_table.csv"), row.names = FALSE)

cat("Demographics table saved to:", file.path(dir.results, "demographics_table.docx"), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. COMPUTE OFFSET (before any cell-type subsetting)
# ══════════════════════════════════════════════════════════════════════════════

sce_full <- SingleCellExperiment(
  assays = list(counts = round(GetAssayData(so_subset, layer = "counts")))
)
sce_full <- computeSumFactors(sce_full)
so_subset@meta.data$pooled_offset <- sizeFactors(sce_full)
cat("Offset computed:", ncol(so_subset), "cells\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. HELPERS
# ══════════════════════════════════════════════════════════════════════════════

get_celltype_column <- function(ct) {
  if (ct %in% c("DCT", "dDCT"))              return("DCT_celltype")
  if (ct %in% c("IC-A", "IC-B", "CNT"))      return("KPMP_celltype2")
  return("celltype2")
}

cell_types_of_interest <- c("EC")
covariates <- NULL   # unadjusted; add c("sex","age","bmi") when available

# ══════════════════════════════════════════════════════════════════════════════
# 4. NEBULA FUNCTION (generic — works for any continuous FFA variable)
# ══════════════════════════════════════════════════════════════════════════════

run_nebula_ffa <- function(seurat_obj, cell_type, ffa_col, covariates = NULL,
                           label = NULL) {
  if (is.null(label)) label <- ffa_col
  
  cat("\n──────────────────────────────────────────\n")
  cat("Variable:", label, "| Cell type:", cell_type, "\n")
  
  ct_col <- get_celltype_column(cell_type)
  so_ct  <- subset(seurat_obj,
                   cells = which(seurat_obj@meta.data[[ct_col]] == cell_type))
  cat("Cells:", ncol(so_ct), "\n")
  
  if (ncol(so_ct) < 20) { cat("  Too few cells — skipping.\n"); return(NULL) }
  
  count_mat <- round(GetAssayData(so_ct, layer = "counts"))
  meta      <- so_ct@meta.data
  
  if (!ffa_col %in% colnames(meta)) {
    cat("  Column", ffa_col, "not found — skipping.\n"); return(NULL)
  }
  
  pred_vars <- c(ffa_col, covariates)
  keep      <- complete.cases(meta[, pred_vars, drop = FALSE])
  cat("  Complete cases:", sum(keep), "/", nrow(meta), "\n")
  
  if (sum(keep) < 20) { cat("  Too few complete cases — skipping.\n"); return(NULL) }
  
  count_mat <- count_mat[, keep, drop = FALSE]
  meta      <- meta[keep, ]
  
  meta$ffa_scaled <- scale(meta[[ffa_col]])[, 1]
  
  formula_str <- if (!is.null(covariates)) {
    paste("~ffa_scaled +", paste(covariates, collapse = " + "))
  } else "~ffa_scaled"
  
  pred_mat <- model.matrix(as.formula(formula_str), data = meta)
  cat("  Formula:", formula_str, "\n")
  cat("  Subjects:", length(unique(meta$record_id)), "\n")
  
  # Filter lowly expressed genes (expressed in ≥10 cells)
  expr_cells <- rowSums(count_mat > 0)
  count_mat  <- count_mat[expr_cells >= 10, , drop = FALSE]
  cat("  Genes after filter:", nrow(count_mat), "\n")
  
  library_sizes <- meta$pooled_offset
  
  data_g <- group_cell(count  = count_mat,
                       id     = meta$record_id,
                       pred   = pred_mat,
                       offset = library_sizes)
  if (is.null(data_g)) {
    data_g <- list(count   = count_mat, id = meta$record_id,
                   pred    = pred_mat,  library = library_sizes)
  }
  
  offset_use <- if ("library" %in% names(data_g)) data_g$library else data_g$offset
  
  res <- nebula(count     = data_g$count, id = data_g$id,
                pred      = data_g$pred,  offset = offset_use,
                ncore     = 1, reml = TRUE, model = "NBLMM",
                output_re = FALSE, covariance = TRUE)
  
  df_res          <- as.data.frame(res$summary)
  df_res$cell_type <- cell_type
  df_res$ffa_var   <- label
  df_res
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. RUN NEBULA — BASELINE FFA
# ══════════════════════════════════════════════════════════════════════════════

results_baseline <- lapply(cell_types_of_interest, function(ct) {
  run_nebula_ffa(so_subset, ct, ffa_col = "baseline_ffa",
                 covariates = covariates, label = "baseline_ffa")
})
results_baseline_df <- do.call(rbind, results_baseline)

nefa_results <- results_baseline_df %>%
  group_by(cell_type) %>%
  mutate(FDR = p.adjust(p_ffa_scaled, method = "BH")) %>%
  ungroup() %>%
  arrange(cell_type, FDR)

write.csv(nefa_results,
          file.path(dir.results, "NEFA_NEBULA_baseline_ffa_FDR.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 6. RUN NEBULA — FFA SUPPRESSION
# ══════════════════════════════════════════════════════════════════════════════

results_suppression <- lapply(cell_types_of_interest, function(ct) {
  run_nebula_ffa(so_subset, ct, ffa_col = "ffa_suppression",
                 covariates = covariates, label = "ffa_suppression")
})
results_suppression_df <- do.call(rbind, results_suppression)

suppression_results <- results_suppression_df %>%
  group_by(cell_type) %>%
  mutate(FDR = p.adjust(p_ffa_scaled, method = "BH")) %>%
  ungroup() %>%
  arrange(cell_type, FDR)

write.csv(suppression_results,
          file.path(dir.results, "NEFA_NEBULA_ffa_suppression_FDR.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 7. VOLCANO PLOTS (for both analyses)
# ══════════════════════════════════════════════════════════════════════════════

plot_volcano_nefa <- function(df, ct, title_var = "FFA", fdr_thresh = 0.05,
                              logfc_thresh = 0.1) {
  df_ct <- df %>%
    filter(cell_type == ct) %>%
    mutate(
      sig       = !is.na(FDR) & !is.na(logFC_ffa_scaled) & FDR < fdr_thresh & abs(logFC_ffa_scaled) > logfc_thresh,
      direction = case_when(
        sig & logFC_ffa_scaled > 0 ~ paste("Up with", title_var),
        sig & logFC_ffa_scaled < 0 ~ paste("Down with", title_var),
        TRUE                       ~ "NS"
      )
    )
  
  top_genes <- if (sum(df_ct$sig) >= 3) {
    df_ct %>% filter(sig) %>% slice_min(FDR, n = 20)
  } else {
    df_ct %>% slice_min(p_ffa_scaled, n = 20)
  }
  
  ggplot(df_ct, aes(x = logFC_ffa_scaled, y = -log10(p_ffa_scaled),
                    color = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 3, max.overlaps = 30) +
    scale_color_manual(values = setNames(
      c("#e63946", "#457b9d", "grey70"),
      c(paste("Up with", title_var), paste("Down with", title_var), "NS")
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh),
               linetype = "dashed", color = "grey40") +
    labs(title    = paste(title_var, "vs. Transcriptome —", ct, "(T2D)"),
         subtitle = paste0("FDR < ", fdr_thresh, " | |logFC| > ", logfc_thresh,
                           "\nUnadjusted"),
         x        = paste("Log2 FC per 1 SD", title_var),
         y        = "-log10(p-value)", color = NULL) +
    theme_classic(base_size = 13) +
    theme(legend.position = "top")
}

for (ct in cell_types_of_interest) {
  p1 <- plot_volcano_nefa(nefa_results, ct, title_var = "Baseline FFA")
  ggsave(file.path(dir.results, paste0("volcano_baseline_ffa_", ct, ".pdf")),
         p1, width = 7, height = 6)
  
  p2 <- plot_volcano_nefa(suppression_results, ct, title_var = "FFA Suppression")
  ggsave(file.path(dir.results, paste0("volcano_ffa_suppression_", ct, ".pdf")),
         p2, width = 7, height = 6)
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. GSEA FUNCTION (shared)
# ══════════════════════════════════════════════════════════════════════════════

# Build gene sets once
msig_h  <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
msig_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)
gene_sets <- bind_rows(msig_h, msig_bp) %>%
  split(x = .$gene_symbol, f = .$gs_name)

run_gsea_ffa <- function(df, ct, label = "FFA") {
  cat("\nRunning GSEA:", label, "—", ct, "\n")
  
  df_ct <- df %>% filter(cell_type == ct)
  if (nrow(df_ct) == 0) { cat("  No results — skipping.\n"); return(NULL) }
  
  ranked <- df_ct %>%
    filter(!is.na(p_ffa_scaled), !is.na(logFC_ffa_scaled)) %>%
    mutate(stat = -log10(p_ffa_scaled) * sign(logFC_ffa_scaled)) %>%
    arrange(desc(stat)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    { setNames(.$stat, .$gene) }
  
  set.seed(123)
  res <- fgsea(pathways = gene_sets, stats = ranked,
               minSize = 15, maxSize = 500)
  res$cell_type <- ct
  res$ffa_var   <- label
  res %>% arrange(pval)
}

# ── GSEA: baseline FFA ───────────────────────────────────────────────────────
gsea_baseline <- lapply(cell_types_of_interest, function(ct) {
  run_gsea_ffa(nefa_results, ct, label = "baseline_ffa")
})
gsea_baseline_df <- do.call(rbind, gsea_baseline) %>%
  group_by(cell_type) %>%
  mutate(FDR_gsea = p.adjust(pval, method = "BH")) %>%
  ungroup() %>% arrange(cell_type, pval)

write.csv(gsea_baseline_df %>% dplyr::select(-leadingEdge),
          file.path(dir.results, "GSEA_baseline_ffa.csv"), row.names = FALSE)

# ── GSEA: FFA suppression ────────────────────────────────────────────────────
gsea_suppression <- lapply(cell_types_of_interest, function(ct) {
  run_gsea_ffa(suppression_results, ct, label = "ffa_suppression")
})
gsea_suppression_df <- do.call(rbind, gsea_suppression) %>%
  group_by(cell_type) %>%
  mutate(FDR_gsea = p.adjust(pval, method = "BH")) %>%
  ungroup() %>% arrange(cell_type, pval)

write.csv(gsea_suppression_df %>% dplyr::select(-leadingEdge),
          file.path(dir.results, "GSEA_ffa_suppression.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 9. GSEA DOT PLOTS
# ══════════════════════════════════════════════════════════════════════════════

plot_gsea_dot <- function(df, ct, title_var, n_top = 20) {
  df_plot <- df %>%
    filter(cell_type == ct) %>%
    slice_min(pval, n = n_top) %>%
    mutate(pathway = str_wrap(gsub("_", " ", pathway), 40))
  
  ggplot(df_plot, aes(x = NES, y = reorder(pathway, NES),
                      size = size, color = pval)) +
    geom_point() +
    scale_color_gradient(low = "#e63946", high = "grey80", name = "p-value") +
    scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(title    = paste("GSEA —", title_var, "—", ct, "(T2D)"),
         subtitle = paste0("Top ", n_top, " pathways | Hallmark + GO:BP"),
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
}

for (ct in cell_types_of_interest) {
  p1 <- plot_gsea_dot(gsea_baseline_df, ct, "Baseline FFA")
  ggsave(file.path(dir.results, paste0("GSEA_baseline_ffa_", ct, ".pdf")),
         p1, width = 9, height = 7)
  
  p2 <- plot_gsea_dot(gsea_suppression_df, ct, "FFA Suppression")
  ggsave(file.path(dir.results, paste0("GSEA_ffa_suppression_", ct, ".pdf")),
         p2, width = 9, height = 7)
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. COMPARE BASELINE FFA vs. FFA SUPPRESSION GSEA
# ══════════════════════════════════════════════════════════════════════════════

# Merge both GSEA results by pathway and cell type
gsea_compare <- gsea_baseline_df %>%
  dplyr::select(pathway, cell_type, NES_baseline = NES,
                pval_baseline = pval, FDR_baseline = FDR_gsea) %>%
  full_join(
    gsea_suppression_df %>%
      dplyr::select(pathway, cell_type, NES_suppression = NES,
                    pval_suppression = pval, FDR_suppression = FDR_gsea),
    by = c("pathway", "cell_type")
  ) %>%
  mutate(
    sig_baseline    = pval_baseline < 0.05,
    sig_suppression = pval_suppression < 0.05,
    concordant      = sign(NES_baseline) == sign(NES_suppression),
    overlap_group   = case_when(
      sig_baseline & sig_suppression & concordant  ~ "Both (concordant)",
      sig_baseline & sig_suppression & !concordant ~ "Both (discordant)",
      sig_baseline & !sig_suppression              ~ "Baseline FFA only",
      !sig_baseline & sig_suppression              ~ "FFA suppression only",
      TRUE                                         ~ "Neither"
    )
  ) %>%
  arrange(cell_type, pval_baseline)

write.csv(gsea_compare,
          file.path(dir.results, "GSEA_comparison_baseline_vs_suppression.csv"),
          row.names = FALSE)

# Summary of overlap
cat("\n── GSEA overlap summary ──\n")
print(table(gsea_compare$overlap_group))

# NES scatter plot: baseline FFA vs. FFA suppression
for (ct in cell_types_of_interest) {
  df_plot <- gsea_compare %>%
    filter(cell_type == ct, !is.na(NES_baseline), !is.na(NES_suppression)) %>%
    mutate(label = ifelse(sig_baseline | sig_suppression, 
                          str_wrap(gsub("_", " ", pathway), 30), NA))
  
  p <- ggplot(df_plot, aes(x = NES_baseline, y = NES_suppression,
                           color = overlap_group)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                    na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    scale_color_manual(values = c(
      "Both (concordant)"    = "#2a9d8f",
      "Both (discordant)"    = "#e76f51",
      "Baseline FFA only"    = "#457b9d",
      "FFA suppression only" = "#e63946",
      "Neither"              = "grey80"
    )) +
    labs(title    = paste("GSEA: Baseline FFA vs. FFA Suppression —", ct),
         subtitle = "Each point = one gene set | Diagonal = perfect concordance",
         x        = "NES (Baseline FFA)",
         y        = "NES (FFA Suppression)",
         color    = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(dir.results, paste0("GSEA_scatter_comparison_", ct, ".pdf")),
         p, width = 8, height = 7)
  print(p)
}

# ══════════════════════════════════════════════════════════════════════════════
# 11. CORRELATION PLOTS — FFA VARIABLES
# ══════════════════════════════════════════════════════════════════════════════

library(ggpubr)

# Get one row per participant with all FFA variables
ffa_cor_df <- so_subset@meta.data %>%
  distinct(record_id, .keep_all = TRUE) %>%
  filter(!is.na(baseline_ffa)) %>%
  dplyr::select(record_id, baseline_ffa, ffa_suppression, ffa_suppression_combined) %>%
  mutate(across(c(baseline_ffa, ffa_suppression, ffa_suppression_combined),
                ~ as.numeric(as.character(.x))))

# Helper: scatter with regression line, r, and p-value
scatter_cor <- function(df, x_var, y_var, x_lab, y_lab) {
  df_clean <- df %>% filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))
  ct <- cor.test(df_clean[[x_var]], df_clean[[y_var]])
  r  <- round(ct$estimate, 2)
  p  <- round(ct$p.value, 3)
  p_label <- ifelse(p < 0.001, "p < 0.001", paste0("p = ", p))
  
  ggplot(df_clean, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(size = 2.5, alpha = 0.8, color = "#457b9d") +
    geom_smooth(method = "lm", se = TRUE, color = "#e63946", fill = "#e63946",
                alpha = 0.15, linewidth = 1) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("r = ", r, "\n", p_label),
             hjust = 1.1, vjust = 1.5, size = 4, color = "grey30") +
    labs(x = x_lab, y = y_lab) +
    theme_classic(base_size = 12)
}

p_base_supp <- scatter_cor(
  ffa_cor_df, "baseline_ffa", "ffa_suppression",
  "Baseline FFA (µmol/L)", "FFA Suppression (%)"
)

p_base_supp_combined <- scatter_cor(
  ffa_cor_df, "baseline_ffa", "ffa_suppression_combined",
  "Baseline FFA (µmol/L)", "FFA Suppression Combined (%)"
)

p_supp_supp_combined <- scatter_cor(
  ffa_cor_df, "ffa_suppression", "ffa_suppression_combined",
  "FFA Suppression (%)", "FFA Suppression Combined (%)"
)

# Combine into one figure
library(patchwork)
p_cor_combined <- (p_base_supp | p_base_supp_combined | p_supp_supp_combined) +
  plot_annotation(
    title    = "Correlations between FFA variables — T2D biopsy cohort",
    subtitle = paste0("n = ", nrow(ffa_cor_df %>% filter(!is.na(baseline_ffa))))
  )

ggsave(file.path(dir.results, "FFA_variable_correlations.pdf"),
       p_cor_combined, width = 12, height = 4)
print(p_cor_combined)

cat("\nAll analyses complete. Results saved to:", dir.results, "\n")







####Plotting

library(patchwork)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(fgsea)
library(stringr)

# ══════════════════════════════════════════════════════════════════════════════
# LOAD PRE-COMPUTED RESULTS (saved by your NEBULA script)
# ══════════════════════════════════════════════════════════════════════════════

dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/NEFA_NEBULA/'

nefa_results        <- read.csv(file.path(dir.results, "NEFA_NEBULA_baseline_ffa_FDR.csv"))
suppression_results <- read.csv(file.path(dir.results, "NEFA_NEBULA_ffa_suppression_FDR.csv"))
gsea_baseline_df    <- read.csv(file.path(dir.results, "GSEA_baseline_ffa.csv"))
gsea_suppression_df <- read.csv(file.path(dir.results, "GSEA_ffa_suppression.csv"))
gsea_compare        <- read.csv(file.path(dir.results, "GSEA_comparison_baseline_vs_suppression.csv"))

ct <- "EC"

# ══════════════════════════════════════════════════════════════════════════════
# SHARED THEME
# ══════════════════════════════════════════════════════════════════════════════

theme_fig <- theme_classic(base_size = 10) +
  theme(
    legend.position  = "top",
    legend.key.size  = unit(0.4, "cm"),
    legend.text      = element_text(size = 8),
    plot.title       = element_text(size = 10, face = "bold"),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    axis.text        = element_text(size = 8),
    axis.title       = element_text(size = 9)
  )

fdr_thresh  <- 0.05
lfc_thresh  <- 0.1
n_top_gsea  <- 20

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — Volcano: Baseline FFA
# ══════════════════════════════════════════════════════════════════════════════

make_volcano <- function(df, ct, title_var, fdr_t = fdr_thresh, lfc_t = lfc_thresh) {
  up_lab   <- paste("Up with", title_var)
  down_lab <- paste("Down with", title_var)
  
  df_ct <- df %>%
    filter(cell_type == ct) %>%
    filter(!(abs(logFC_ffa_scaled) > 10 & (is.na(FDR) | FDR >= fdr_t))) %>%
    mutate(
      sig = !is.na(FDR) & !is.na(logFC_ffa_scaled) &
        FDR < fdr_t & abs(logFC_ffa_scaled) > lfc_t,
      direction = case_when(
        sig & logFC_ffa_scaled > 0 ~ up_lab,
        sig & logFC_ffa_scaled < 0 ~ down_lab,
        TRUE ~ "NS"
      )
    )
  
  top_genes <- if (sum(df_ct$sig) >= 1) {
    df_ct %>% filter(sig) %>% slice_min(FDR, n = 20)
  } else {
    df_ct %>% slice_min(p_ffa_scaled, n = 5)
  }
  
  ggplot(df_ct, aes(logFC_ffa_scaled, -log10(p_ffa_scaled), color = direction)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 25, show.legend = FALSE) +
    scale_color_manual(values = setNames(
      c("#e63946", "#457b9d", "grey75"),
      c(up_lab, down_lab, "NS")
    ), name = NULL) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = c(-lfc_t, lfc_t), linetype = "dashed", linewidth = 0.4, color = "grey50") +
    labs(
      title    = paste0(title_var, " vs. Transcriptome — ", ct, " (T2D)"),
      subtitle = paste0("FDR < ", fdr_t, " | |logFC| > ", lfc_t, "  ·  Unadjusted"),
      x        = paste("Log2 FC per 1 SD", title_var),
      y        = expression(-log[10](p))
    ) +
    theme_fig
}

pA <- make_volcano(nefa_results,        ct, "Baseline FFA")
pB <- make_volcano(suppression_results, ct, "FFA Suppression")

# ══════════════════════════════════════════════════════════════════════════════
# PANELS C & D — GSEA dot plots
# ══════════════════════════════════════════════════════════════════════════════

make_gsea_dot <- function(df, ct, title_var, n_top = n_top_gsea) {
  df_plot <- df %>%
    filter(cell_type == ct) %>%
    slice_min(pval, n = n_top) %>%
    mutate(
      pathway  = str_wrap(gsub("_", " ", gsub("GOBP_|HALLMARK_", "", pathway)), 38),
      neg_logp = -log10(pval)
    )
  
  ggplot(df_plot, aes(NES, reorder(pathway, NES), size = size, color = pval)) +
    geom_point() +
    scale_color_gradient(low = "#e63946", high = "grey80",
                         name = "p-value", guide = guide_colorbar(barwidth = 4, barheight = 0.5)) +
    scale_size_continuous(name = "Set size", range = c(2, 7),
                          guide = guide_legend(override.aes = list(color = "grey50"))) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
    labs(
      title    = paste0("GSEA — ", title_var, " — ", ct, " (T2D)"),
      subtitle = paste0("Top ", n_top, " pathways | Hallmark + GO:BP"),
      x = "Normalized Enrichment Score (NES)", y = NULL
    ) +
    theme_fig +
    theme(axis.text.y = element_text(size = 7))
}

pC <- make_gsea_dot(gsea_baseline_df,    ct, "Baseline FFA")
pD <- make_gsea_dot(gsea_suppression_df, ct, "FFA Suppression")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL E — NES concordance scatter
# ══════════════════════════════════════════════════════════════════════════════

overlap_colors <- c(
  "Both (concordant)"    = "#2a9d8f",
  "Both (discordant)"    = "#e76f51",
  "Baseline FFA only"    = "#457b9d",
  "FFA suppression only" = "#e63946",
  "Neither"              = "grey82"
)

df_e <- gsea_compare %>%
  filter(cell_type == ct, !is.na(NES_baseline), !is.na(NES_suppression)) %>%
  mutate(
    label = ifelse(
      (sig_baseline | sig_suppression) & overlap_group != "Neither",
      str_wrap(gsub("_", " ", gsub("GOBP_|HALLMARK_", "", pathway)), 28),
      NA_character_
    )
  )

pE <- ggplot(df_e, aes(NES_baseline, NES_suppression, color = overlap_group)) +
  geom_point(alpha = 0.55, size = 1.5) +
  geom_text_repel(aes(label = label), size = 2.3, max.overlaps = 22,
                  na.rm = TRUE, lineheight = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.4, color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey60") +
  scale_color_manual(values = overlap_colors, name = NULL) +
  labs(
    title    = paste0("GSEA: Baseline FFA vs. FFA Suppression — ", ct),
    subtitle = "Each point = one gene set  |  Diagonal = perfect concordance",
    x = "NES (Baseline FFA)", y = "NES (FFA Suppression)"
  ) +
  theme_fig +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 7))

# ══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE COMBINED FIGURE
# ══════════════════════════════════════════════════════════════════════════════
# Layout:
#   Row 1: A (volcano baseline) | B (volcano suppression)
#   Row 2: C (GSEA dot baseline) | D (GSEA dot suppression)
#   Row 3: E (concordance scatter) — centered, narrower

layout <- "
AABB
CCDD
EEEE
"

combined <- pA + pB + pC + pD + pE +
  plot_layout(design = layout, heights = c(1, 1.4, 1.2)) +
  plot_annotation(
    title   = "Free Fatty Acids and the Kidney Endothelial Cell Transcriptome in Youth-Onset T2D",
    caption = "NEBULA mixed-effects model; unadjusted; FDR < 0.05 | |log2FC| > 0.1",
    theme   = theme(
      plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 8, color = "grey50", hjust = 1)
    )
  ) &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

# ══════════════════════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════════════════════

ggsave(
  filename = file.path(dir.results, "Figure_FFA_EC_Transcriptome_combined.pdf"),
  plot     = combined,
  width    = 16, height = 20,
  device   = cairo_pdf
)

ggsave(
  filename = file.path(dir.results, "Figure_FFA_EC_Transcriptome_combined.png"),
  plot     = combined,
  width    = 16, height = 20, dpi = 300
)

cat("Combined figure saved to:", dir.results, "\n")






















#### FFA suppression combined 


library(Seurat)
library(nebula)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(purrr)
library(fgsea)
library(msigdbr)
library(stringr)
library(patchwork)

# ══════════════════════════════════════════════════════════════════════════════
# 0. LOAD & PREPARE SEURAT OBJECT
# ══════════════════════════════════════════════════════════════════════════════

load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')
so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

# ── Load & collapse harmonized data ─────────────────────────────────────────
harmonized_data <- read.csv(
  "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv",
  na = ''
)

dat <- harmonized_data %>%
  dplyr::select(-dob) %>%
  arrange(date_of_screen) %>%
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric),         ~ ifelse(all(is.na(.x)), NA_real_,      mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

# ── Merge clinical variables ─────────────────────────────────────────────────
vars_to_merge <- c("baseline_ffa", "ffa_suppression", "ffa_suppression_combined",
                   "age", "bmi", "acr_u", "hba1c",
                   "epic_insulin_1", "epic_sglti2_1", "epic_glp1ra_1",
                   "epic_tzd_1", "epic_mfm_1", "epic_raasi_1",
                   "epic_statin_1", "epic_fibrate_1")

merge_df <- dat %>%
  dplyr::select(any_of(c("record_id", vars_to_merge))) %>%
  distinct(record_id, .keep_all = TRUE)

for (v in vars_to_merge) {
  if (v %in% colnames(merge_df)) {
    so_subset@meta.data[[v]] <- merge_df[[v]][
      match(so_subset@meta.data$record_id, merge_df$record_id)
    ]
  }
}

cat("baseline_ffa non-NA cells:",            sum(!is.na(so_subset@meta.data$baseline_ffa)), "\n")
cat("ffa_suppression non-NA cells:",         sum(!is.na(so_subset@meta.data$ffa_suppression)), "\n")
cat("ffa_suppression_combined non-NA cells:", sum(!is.na(so_subset@meta.data$ffa_suppression_combined)), "\n")

# ── Cell type annotations ────────────────────────────────────────────────────
so_subset$celltype1 <- case_when(
  grepl("PT-",  so_subset$celltype_rpca) ~ "PT",
  grepl("TAL-", so_subset$celltype_rpca) ~ "TAL",
  grepl("EC-",  so_subset$celltype_rpca) ~ "EC",
  grepl("POD",  so_subset$celltype_rpca) ~ "POD",
  grepl("MAC",  so_subset$celltype_rpca) ~ "MAC",
  grepl("MON",  so_subset$celltype_rpca) ~ "MON",
  grepl("PC-",  so_subset$celltype_rpca) ~ "PC",
  grepl("FIB",  so_subset$celltype_rpca) ~ "FIB_MC_VSMC",
  grepl("DTL",  so_subset$celltype_rpca) ~ "DTL",
  so_subset$celltype_rpca == "DCT"        ~ "DCT",
  so_subset$celltype_rpca == "ATL"        ~ "ATL",
  so_subset$celltype_rpca == "B"          ~ "B",
  so_subset$celltype_rpca == "T"          ~ "T"
)
so_subset$celltype1      <- as.character(so_subset$celltype1)
so_subset$KPMP_celltype2 <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2      <- ifelse(
  so_subset$KPMP_celltype %in% c("aPT", "PT-S1/S2", "PT-S3"), "PT",
  ifelse(grepl("TAL", so_subset$KPMP_celltype), "TAL",
         ifelse(grepl("EC-", so_subset$KPMP_celltype), "EC",
                so_subset$KPMP_celltype2)))
so_subset$DCT_celltype <- ifelse(
  so_subset$KPMP_celltype %in% c("DCT", "dDCT"), "DCT", "Non-DCT")

# ── Subset to T2D, exclude CRC-55 ───────────────────────────────────────────
so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type_2_Diabetes')

dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/NEFA_NEBULA_combined/'
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. COMPUTE OFFSET
# ══════════════════════════════════════════════════════════════════════════════

sce_full <- SingleCellExperiment(
  assays = list(counts = round(GetAssayData(so_subset, layer = "counts")))
)
sce_full <- computeSumFactors(sce_full)
so_subset@meta.data$pooled_offset <- sizeFactors(sce_full)
cat("Offset computed:", ncol(so_subset), "cells\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. HELPERS
# ══════════════════════════════════════════════════════════════════════════════

get_celltype_column <- function(ct) {
  if (ct %in% c("DCT", "dDCT"))         return("DCT_celltype")
  if (ct %in% c("IC-A", "IC-B", "CNT")) return("KPMP_celltype2")
  return("celltype2")
}

cell_types_of_interest <- c("EC")
covariates <- NULL

# ══════════════════════════════════════════════════════════════════════════════
# 3. NEBULA FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

run_nebula_ffa <- function(seurat_obj, cell_type, ffa_col,
                           covariates = NULL, label = NULL) {
  if (is.null(label)) label <- ffa_col
  
  cat("\n──────────────────────────────────────────\n")
  cat("Variable:", label, "| Cell type:", cell_type, "\n")
  
  ct_col <- get_celltype_column(cell_type)
  so_ct  <- subset(seurat_obj,
                   cells = which(seurat_obj@meta.data[[ct_col]] == cell_type))
  cat("Cells:", ncol(so_ct), "\n")
  
  if (ncol(so_ct) < 20) { cat("  Too few cells — skipping.\n"); return(NULL) }
  
  count_mat <- round(GetAssayData(so_ct, layer = "counts"))
  meta      <- so_ct@meta.data
  
  if (!ffa_col %in% colnames(meta)) {
    cat("  Column", ffa_col, "not found — skipping.\n"); return(NULL)
  }
  
  pred_vars <- c(ffa_col, covariates)
  keep      <- complete.cases(meta[, pred_vars, drop = FALSE])
  cat("  Complete cases:", sum(keep), "/", nrow(meta), "\n")
  
  if (sum(keep) < 20) { cat("  Too few complete cases — skipping.\n"); return(NULL) }
  
  count_mat <- count_mat[, keep, drop = FALSE]
  meta      <- meta[keep, ]
  
  meta$ffa_scaled <- scale(meta[[ffa_col]])[, 1]
  
  formula_str <- if (!is.null(covariates)) {
    paste("~ffa_scaled +", paste(covariates, collapse = " + "))
  } else "~ffa_scaled"
  
  pred_mat <- model.matrix(as.formula(formula_str), data = meta)
  cat("  Formula:", formula_str, "\n")
  cat("  Subjects:", length(unique(meta$record_id)), "\n")
  
  expr_cells <- rowSums(count_mat > 0)
  count_mat  <- count_mat[expr_cells >= 10, , drop = FALSE]
  cat("  Genes after filter:", nrow(count_mat), "\n")
  
  data_g <- group_cell(count  = count_mat,
                       id     = meta$record_id,
                       pred   = pred_mat,
                       offset = meta$pooled_offset)
  if (is.null(data_g)) {
    data_g <- list(count = count_mat, id = meta$record_id,
                   pred  = pred_mat,  library = meta$pooled_offset)
  }
  
  offset_use <- if ("library" %in% names(data_g)) data_g$library else data_g$offset
  
  res <- nebula(count     = data_g$count, id = data_g$id,
                pred      = data_g$pred,  offset = offset_use,
                ncore     = 1, reml = TRUE, model = "NBLMM",
                output_re = FALSE, covariance = TRUE)
  
  df_res           <- as.data.frame(res$summary)
  df_res$cell_type <- cell_type
  df_res$ffa_var   <- label
  df_res
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. RUN NEBULA — FFA SUPPRESSION COMBINED
# ══════════════════════════════════════════════════════════════════════════════

results_combined <- lapply(cell_types_of_interest, function(ct) {
  run_nebula_ffa(so_subset, ct,
                 ffa_col    = "ffa_suppression_combined",
                 covariates = covariates,
                 label      = "ffa_suppression_combined")
})
results_combined_df <- do.call(rbind, results_combined) %>%
  group_by(cell_type) %>%
  mutate(FDR = p.adjust(p_ffa_scaled, method = "BH")) %>%
  ungroup() %>%
  arrange(cell_type, FDR)

write.csv(results_combined_df,
          file.path(dir.results, "NEFA_NEBULA_ffa_suppression_combined_FDR.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 5. VOLCANO PLOT
# ══════════════════════════════════════════════════════════════════════════════

theme_fig <- theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.4, "cm"),
    legend.text     = element_text(size = 8),
    plot.title      = element_text(size = 10, face = "bold"),
    plot.subtitle   = element_text(size = 8, color = "grey40"),
    axis.text       = element_text(size = 8),
    axis.title      = element_text(size = 9)
  )

fdr_thresh <- 0.05
lfc_thresh <- 0.1

make_volcano <- function(df, ct, title_var,
                         fdr_t = fdr_thresh, lfc_t = lfc_thresh) {
  up_lab   <- paste("Up with", title_var)
  down_lab <- paste("Down with", title_var)
  
  df_ct <- df %>%
    filter(cell_type == ct) %>%
    filter(!(abs(logFC_ffa_scaled) > 10 & (is.na(FDR) | FDR >= fdr_t))) %>%
    mutate(
      sig = !is.na(FDR) & !is.na(logFC_ffa_scaled) &
        FDR < fdr_t & abs(logFC_ffa_scaled) > lfc_t,
      direction = case_when(
        sig & logFC_ffa_scaled > 0 ~ up_lab,
        sig & logFC_ffa_scaled < 0 ~ down_lab,
        TRUE                       ~ "NS"
      )
    )
  
  top_genes <- if (sum(df_ct$sig) >= 1) {
    df_ct %>% filter(sig) %>% slice_min(FDR, n = 20)
  } else {
    df_ct %>% slice_min(p_ffa_scaled, n = 5)
  }
  
  ggplot(df_ct, aes(logFC_ffa_scaled, -log10(p_ffa_scaled), color = direction)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 25, show.legend = FALSE) +
    scale_color_manual(values = setNames(
      c("#e63946", "#457b9d", "grey75"),
      c(up_lab, down_lab, "NS")
    ), name = NULL) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = c(-lfc_t, lfc_t), linetype = "dashed",
               linewidth = 0.4, color = "grey50") +
    labs(
      title    = paste0(title_var, " vs. Transcriptome — ", ct, " (T2D)"),
      subtitle = paste0("FDR < ", fdr_t, " | |logFC| > ", lfc_t, "  ·  Unadjusted"),
      x        = paste("Log2 FC per 1 SD", title_var),
      y        = expression(-log[10](p))
    ) +
    theme_fig
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. GSEA
# ══════════════════════════════════════════════════════════════════════════════

msig_h  <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
msig_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)
gene_sets <- bind_rows(msig_h, msig_bp) %>%
  split(x = .$gene_symbol, f = .$gs_name)

run_gsea_ffa <- function(df, ct, label = "FFA") {
  cat("\nRunning GSEA:", label, "—", ct, "\n")
  df_ct <- df %>% filter(cell_type == ct)
  if (nrow(df_ct) == 0) { cat("  No results — skipping.\n"); return(NULL) }
  
  ranked <- df_ct %>%
    filter(!is.na(p_ffa_scaled), !is.na(logFC_ffa_scaled)) %>%
    mutate(stat = -log10(p_ffa_scaled) * sign(logFC_ffa_scaled)) %>%
    arrange(desc(stat)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    { setNames(.$stat, .$gene) }
  
  set.seed(123)
  res <- fgsea(pathways = gene_sets, stats = ranked,
               minSize = 15, maxSize = 500)
  res$cell_type <- ct
  res$ffa_var   <- label
  res %>% arrange(pval)
}

gsea_combined_df <- lapply(cell_types_of_interest, function(ct) {
  run_gsea_ffa(results_combined_df, ct, label = "ffa_suppression_combined")
}) %>%
  do.call(rbind, .) %>%
  group_by(cell_type) %>%
  mutate(FDR_gsea = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  arrange(cell_type, pval)

write.csv(gsea_combined_df %>% dplyr::select(-leadingEdge),
          file.path(dir.results, "GSEA_ffa_suppression_combined.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 7. GSEA DOT PLOT
# ══════════════════════════════════════════════════════════════════════════════

make_gsea_dot <- function(df, ct, title_var, n_top = 20) {
  df_plot <- df %>%
    filter(cell_type == ct) %>%
    slice_min(pval, n = n_top) %>%
    mutate(
      collection = case_when(
        str_starts(pathway, "HALLMARK_") ~ "Hallmark",
        str_starts(pathway, "GOBP_")     ~ "GO:BP",
        TRUE                             ~ "Other"
      ),
      pathway = str_wrap(gsub("_", " ", gsub("GOBP_|HALLMARK_", "", pathway)), 38)
    )
  
  ggplot(df_plot, aes(NES, reorder(pathway, NES),
                      size = size, color = pval, shape = collection)) +
    geom_point() +
    scale_color_gradient(low = "#e63946", high = "grey80", name = "p-value",
                         guide = guide_colorbar(barwidth = 4, barheight = 0.5)) +
    scale_size_continuous(name = "Set size", range = c(2, 7),
                          guide = guide_legend(override.aes = list(color = "grey50"))) +
    scale_shape_manual(values = c("GO:BP" = 16, "Hallmark" = 17, "Other" = 15),
                       name = "Collection") +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.4, color = "grey50") +
    labs(
      title    = paste0("GSEA — ", title_var, " — ", ct, " (T2D)"),
      subtitle = "Top 20 pathways | Hallmark + GO:BP",
      x = "Normalized Enrichment Score (NES)", y = NULL
    ) +
    theme_fig +
    theme(axis.text.y = element_text(size = 7))
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. ASSEMBLE & SAVE COMBINED FIGURE
# ══════════════════════════════════════════════════════════════════════════════

for (ct in cell_types_of_interest) {
  pA <- make_volcano(results_combined_df, ct, "FFA Suppression Combined")
  pB <- make_gsea_dot(gsea_combined_df,   ct, "FFA Suppression Combined")
  
  combined <- pA + pB +
    plot_layout(ncol = 2, widths = c(1, 1.4)) +
    plot_annotation(
      title   = "FFA Suppression (Combined) and the Kidney Endothelial Cell Transcriptome in Youth-Onset T2D",
      caption = "NEBULA mixed-effects model; unadjusted; FDR < 0.05 | |log2FC| > 0.1",
      theme   = theme(
        plot.title   = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.caption = element_text(size = 8, color = "grey50", hjust = 1)
      )
    ) &
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 11, face = "bold"))
  
  ggsave(
    filename = file.path(dir.results, paste0("Figure_FFA_suppression_combined_", ct, ".pdf")),
    plot     = combined,
    width = 16, height = 9, device = cairo_pdf
  )
  ggsave(
    filename = file.path(dir.results, paste0("Figure_FFA_suppression_combined_", ct, ".png")),
    plot     = combined,
    width = 16, height = 9, dpi = 300
  )
  print(combined)
}

cat("\nAll analyses complete. Results saved to:", dir.results, "\n")