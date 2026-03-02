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


# ══════════════════════════════════════════════════════════════════════════════
# 0. LOAD & PREPARE SEURAT OBJECT (T2D only, matching your established setup)
# ══════════════════════════════════════════════════════════════════════════════

load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')
so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

# ── Merge baseline_ffa from harmonized dataset ──────────────────────────────
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

# Check how many T2D participants have baseline_ffa
cat("Participants with baseline_ffa (T2D):",
    sum(!is.na(dat$baseline_ffa[dat$group == "Type_2_Diabetes"])), "\n")

# Merge into Seurat metadata by record_id (one value per participant → all cells get same value)
ffa_merge <- dat %>%
  dplyr::select(record_id, baseline_ffa) %>%
  distinct()

so_subset@meta.data <- so_subset@meta.data %>%
  dplyr::select(-any_of("baseline_ffa")) %>%   # drop old all-NA column
  left_join(ffa_merge, by = "record_id")

cat("Cells with baseline_ffa after merge:",
    sum(!is.na(so_subset@meta.data$baseline_ffa)), "\n")
# ────────────────────────────────────────────────────────────────────────────

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
so_subset$celltype1       <- as.character(so_subset$celltype1)
so_subset$KPMP_celltype2  <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2       <- ifelse(
  so_subset$KPMP_celltype == "aPT" |
    so_subset$KPMP_celltype == "PT-S1/S2" |
    so_subset$KPMP_celltype == "PT-S3", "PT",
  ifelse(grepl("TAL", so_subset$KPMP_celltype), "TAL",
         ifelse(grepl("EC-", so_subset$KPMP_celltype), "EC",
                so_subset$KPMP_celltype2)))
so_subset$DCT_celltype <- ifelse(
  so_subset$KPMP_celltype == "DCT" | so_subset$KPMP_celltype == "dDCT",
  "DCT", "Non-DCT")

# Remove excluded participant and restrict to T2D
so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type_2_Diabetes')

dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/NEFA_NEBULA/'
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. USER SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

nefa_col   <- "baseline_ffa"          # baseline (fasting) free fatty acids
covariates <- c("sex", "age", "bmi")  # set to NULL to run unadjusted

# Primary focus: EC. Add others if desired e.g. c("EC","PT","TAL","POD")
cell_types_of_interest <- c("EC")

# ══════════════════════════════════════════════════════════════════════════════
# 2. COMPUTE OFFSET ON FULL T2D OBJECT (before any cell-type subsetting)
# ══════════════════════════════════════════════════════════════════════════════

sce_full <- SingleCellExperiment(
  assays = list(counts = round(GetAssayData(so_subset, layer = "counts")))
)
sce_full <- computeSumFactors(sce_full)
so_subset@meta.data$pooled_offset <- sizeFactors(sce_full)
cat("Offset computed on T2D object:", ncol(so_subset), "cells\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. HELPER: cell-type column logic (matches your established pattern)
#    - EC, PT, TAL, POD etc. → celltype2
#    - DCT subtypes          → DCT_celltype
#    - Fine-grained KPMP     → KPMP_celltype2
# ══════════════════════════════════════════════════════════════════════════════

get_celltype_column <- function(ct) {
  dct_types  <- c("DCT", "dDCT")
  kpmp_types <- c("IC-A", "IC-B", "DCT1", "DCT2", "CNT")
  if (ct %in% dct_types)  return("DCT_celltype")
  if (ct %in% kpmp_types) return("KPMP_celltype2")
  return("celltype2")
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. NEBULA FUNCTION — baseline FFA as continuous predictor
# ══════════════════════════════════════════════════════════════════════════════

run_nebula_nefa <- function(seurat_obj, cell_type, nefa_col, covariates = NULL) {
  
  cat("\n──────────────────────────────────────────\n")
  cat("Cell type:", cell_type, "\n")
  
  # Subset to cell type
  ct_col <- get_celltype_column(cell_type)
  so_ct  <- subset(seurat_obj,
                   cells = which(seurat_obj@meta.data[[ct_col]] == cell_type))
  cat("Cells:", ncol(so_ct), "\n")
  
  if (ncol(so_ct) < 20) {
    cat("  Too few cells — skipping.\n")
    return(NULL)
  }
  
  count_mat <- round(GetAssayData(so_ct, layer = "counts"))
  meta      <- so_ct@meta.data
  
  # Check baseline_ffa is present
  if (!nefa_col %in% colnames(meta)) {
    cat("  Column", nefa_col, "not found in metadata — skipping.\n")
    return(NULL)
  }
  
  # Drop rows with NA in FFA or covariates
  pred_vars <- c(nefa_col, covariates)
  keep      <- complete.cases(meta[, pred_vars, drop = FALSE])
  cat("  Complete cases:", sum(keep), "/", nrow(meta), "\n")
  
  if (sum(keep) < 20) {
    cat("  Too few complete cases — skipping.\n")
    return(NULL)
  }
  
  count_mat <- count_mat[, keep, drop = FALSE]
  meta      <- meta[keep, ]
  
  # Scale FFA (mean 0, SD 1) for model stability
  # logFC interpreted as transcriptional change per 1 SD increase in baseline FFA
  meta$ffa_scaled <- scale(meta[[nefa_col]])[, 1]
  
  formula_str <- if (!is.null(covariates)) {
    paste("~ffa_scaled +", paste(covariates, collapse = " + "))
  } else {
    "~ffa_scaled"
  }
  pred_mat <- model.matrix(as.formula(formula_str), data = meta)
  cat("  Formula:", formula_str, "\n")
  cat("  Subjects (record_id):", length(unique(meta$record_id)), "\n")
  
  library_sizes <- meta$pooled_offset
  gene_names    <- rownames(count_mat)
  results_list  <- vector("list", length(gene_names))
  
  for (i in seq_along(gene_names)) {
    g          <- gene_names[i]
    count_gene <- count_mat[g, , drop = FALSE]
    
    data_g <- group_cell(count  = count_gene,
                         id     = meta$record_id,
                         pred   = pred_mat,
                         offset = library_sizes)
    
    if (is.null(data_g)) {
      data_g <- list(count   = count_gene,
                     id      = meta$record_id,
                     pred    = pred_mat,
                     library = library_sizes)
    }
    
    offset_use <- if ("library" %in% names(data_g)) data_g$library else data_g$offset
    
    tryCatch({
      res    <- nebula(count      = data_g$count,
                       id         = data_g$id,
                       pred       = data_g$pred,
                       offset     = offset_use,
                       ncore      = 1,
                       reml       = TRUE,
                       model      = "NBLMM",
                       output_re  = FALSE,
                       covariance = TRUE)
      df_res            <- as.data.frame(res$summary)
      df_res$gene       <- g
      df_res$cell_type  <- cell_type
      results_list[[i]] <- df_res
    }, error = function(e) { })
  }
  
  do.call(rbind, results_list)
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. RUN
# ══════════════════════════════════════════════════════════════════════════════

all_results <- lapply(cell_types_of_interest, function(ct) {
  run_nebula_nefa(seurat_obj = so_subset,
                  cell_type  = ct,
                  nefa_col   = nefa_col,
                  covariates = covariates)
})

all_results_df <- do.call(rbind, all_results)

write.csv(all_results_df,
          file = file.path(dir.results, "NEFA_NEBULA_all_results.csv"),
          row.names = FALSE)

# After running, check column names and adjust below:
# colnames(all_results_df)

# ══════════════════════════════════════════════════════════════════════════════
# 6. EXTRACT FFA COEFFICIENT & FDR CORRECTION
#    *** Run colnames(all_results_df) first and adjust column names below ***
# ══════════════════════════════════════════════════════════════════════════════

# Typical NEBULA output columns look like:
#   logFC_ffa_scaled   p_ffa_scaled   (confirm with colnames)

nefa_results <- all_results_df %>%
  group_by(cell_type) %>%
  mutate(FDR = p.adjust(p_ffa_scaled, method = "BH")) %>%  # adjust if col name differs
  ungroup() %>%
  arrange(cell_type, FDR)

write.csv(nefa_results,
          file = file.path(dir.results, "NEFA_NEBULA_ffa_coeff_FDR.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 7. VOLCANO PLOT
# ══════════════════════════════════════════════════════════════════════════════

plot_volcano_nefa <- function(df, ct, fdr_thresh = 0.05, logfc_thresh = 0.1) {
  df_ct <- df %>%
    filter(cell_type == ct) %>%
    mutate(
      sig       = FDR < fdr_thresh & abs(logFC_ffa_scaled) > logfc_thresh,
      direction = case_when(
        sig & logFC_ffa_scaled > 0 ~ "Up with FFA",
        sig & logFC_ffa_scaled < 0 ~ "Down with FFA",
        TRUE                       ~ "NS"
      )
    )
  
  top_genes <- df_ct %>% filter(sig) %>% slice_min(FDR, n = 20)
  
  ggplot(df_ct, aes(x = logFC_ffa_scaled, y = -log10(p_ffa_scaled), color = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("Up with FFA"   = "#e63946",
                                  "Down with FFA" = "#457b9d",
                                  "NS"            = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh),
               linetype = "dashed", color = "grey40") +
    labs(title    = paste("Baseline FFA vs. Transcriptome —", ct, "(T2D)"),
         subtitle = paste0("FDR < ", fdr_thresh, " | |logFC| > ", logfc_thresh,
                           "\nAdjusted for sex, age, BMI"),
         x        = "Log2 FC per 1 SD baseline FFA",
         y        = "-log10(p-value)",
         color    = NULL) +
    theme_classic(base_size = 13) +
    theme(legend.position = "top")
}

for (ct in cell_types_of_interest) {
  p <- plot_volcano_nefa(nefa_results, ct)
  ggsave(filename = file.path(dir.results, paste0("volcano_FFA_", ct, ".pdf")),
         plot = p, width = 7, height = 6)
  print(p)
}


# ══════════════════════════════════════════════════════════════════════════════
# 8. GSEA — MSigDB Hallmark + GO:BP gene sets
# ══════════════════════════════════════════════════════════════════════════════

library(fgsea)
library(msigdbr)

# Build gene sets: Hallmark + GO Biological Process
msig_h  <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
msig_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets <- bind_rows(msig_h, msig_bp) %>%
  split(x = .$gene_symbol, f = .$gs_name)

run_gsea_nefa <- function(df, ct) {
  cat("\nRunning GSEA for:", ct, "\n")
  
  df_ct <- df %>% filter(cell_type == ct)
  
  if (nrow(df_ct) == 0) {
    cat("  No results for", ct, "— skipping.\n")
    return(NULL)
  }
  
  # Ranking statistic: -log10(p) * sign(logFC) — prioritises significance + direction
  ranked <- df_ct %>%
    filter(!is.na(p_ffa_scaled), !is.na(logFC_ffa_scaled)) %>%
    mutate(stat = -log10(p_ffa_scaled) * sign(logFC_ffa_scaled)) %>%
    arrange(desc(stat)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    { setNames(.$stat, .$gene) }
  
  set.seed(123)
  res <- fgsea(pathways  = gene_sets,
               stats     = ranked,
               minSize   = 15,
               maxSize   = 500)
  
  res$cell_type <- ct
  res %>% arrange(pval)
}

gsea_results <- lapply(cell_types_of_interest, run_gsea_nefa)
gsea_results_df <- do.call(rbind, gsea_results) %>%
  group_by(cell_type) %>%
  mutate(FDR_gsea = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  arrange(cell_type, pval)

write.csv(gsea_results_df %>% dplyr::select(-leadingEdge),
          file = file.path(dir.results, "NEFA_GSEA_results.csv"),
          row.names = FALSE)

# Dot plot: top 20 pathways by nominal p-value per cell type
for (ct in cell_types_of_interest) {
  df_plot <- gsea_results_df %>%
    filter(cell_type == ct) %>%
    slice_min(pval, n = 20) %>%
    mutate(pathway = gsub("_", " ", pathway),
           pathway = stringr::str_wrap(pathway, 40),
           direction = ifelse(NES > 0, "Positive NES", "Negative NES"))
  
  p <- ggplot(df_plot, aes(x = NES, y = reorder(pathway, NES),
                           size = size, color = pval)) +
    geom_point() +
    scale_color_gradient(low = "#e63946", high = "grey80", name = "p-value") +
    scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(title    = paste("GSEA — Baseline FFA —", ct, "(T2D)"),
         subtitle = "Top 20 pathways by p-value | Hallmark + GO:BP",
         x        = "Normalized Enrichment Score (NES)",
         y        = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(filename = file.path(dir.results, paste0("GSEA_dotplot_FFA_", ct, ".pdf")),
         plot = p, width = 9, height = 7)
  print(p)
}

cat("\nGSEA complete. Results saved to:", dir.results, "\n")



