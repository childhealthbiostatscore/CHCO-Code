#### Pioglitazone scRNAseq - Covariate Sensitivity Analysis
#### Tests BMI, sex, metformin, and GLP-1RA as individual covariates
#### Compares results against base model (pioglitazone only)

library(scran)
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggpubr)
library(rstatix)
library(data.table)
library(pheatmap)
library(readxl)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(dplyr)
library(nebula)
library(Matrix)
library(ggplot2)
library(gt)
library(gtsummary)
library(fgsea)
library(msigdbr)
library(stringr)
library(gridExtra)
library(grid)
library(ggrepel)

# ============================================================================
# SETUP
# ============================================================================

setwd('C:/Users/netio/Documents/Harmonized_data/')
dir.results <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/Sensitivity/'

dir.create(paste0(dir.results, 'NEBULA/'),      recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'Comparison/'),  recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'VolcanoPlots/'),recursive = TRUE, showWarnings = FALSE)

COVARIATES <- list(
  base      = NULL,
  bmi       = "bmi",
  sex       = "sex",
  metformin = "epic_mfm",
  glp1ra    = "epic_glp1ra"
)

celltypes_vec <- c('All', 'PT', 'TAL', 'EC', 'POD', 'DCT', 'IC')

# ============================================================================
# DATA LOADING
# ============================================================================

load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')
so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

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
  so_subset$celltype_rpca == "DCT"       ~ "DCT",
  so_subset$celltype_rpca == "ATL"       ~ "ATL",
  so_subset$celltype_rpca == "B"         ~ "B",
  so_subset$celltype_rpca == "T"         ~ "T"
)
so_subset$celltype1 <- as.character(so_subset$celltype1)

so_subset$KPMP_celltype2 <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2 <- ifelse(
  so_subset$KPMP_celltype %in% c("aPT", "PT-S1/S2", "PT-S3"), "PT",
  ifelse(grepl("TAL", so_subset$KPMP_celltype), "TAL",
         ifelse(grepl("EC-", so_subset$KPMP_celltype), "EC", so_subset$KPMP_celltype2))
)

so_subset$DCT_celltype <- ifelse(
  so_subset$KPMP_celltype %in% c("DCT", "dDCT"), "DCT", "Non-DCT"
)

so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type_2_Diabetes')

dir.results <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/Sensitivity/'

# ============================================================================
# MEDICATION & COVARIATE DATA
# ============================================================================


harmonized_data <- data.table::fread("harmonized_dataset.csv")
harmonized_data <- harmonized_data %>% filter(group == 'Type 2 Diabetes')

medications <- readxl::read_xlsx("Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>%
  dplyr::select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), '_1', '')

meta.data <- so_subset@meta.data

harmonized_data <- harmonized_data %>%
  filter(rh_id %in% meta.data$record_id | croc_id %in% meta.data$record_id |
           improve_id %in% meta.data$record_id | penguin_id %in% meta.data$record_id |
           rh2_id %in% meta.data$record_id) %>%
  dplyr::select(record_id, mrn, group, bmi, sex)

medications <- medications %>%
  dplyr::select(mrn, tzd, epic_mfm_1=mfm, epic_glp1ra_1=glp1ra) %>%
  filter(mrn %in% harmonized_data$mrn)

harmonized_data$mrn <- as.character(harmonized_data$mrn)

final_df <- medications %>% left_join(harmonized_data, by = 'mrn')
final_df$combined_id <- paste0(final_df$mrn, '_', final_df$record_id)
final_df <- final_df %>% filter(!duplicated(combined_id))

final_df <- medications %>% left_join(harmonized_data, by = 'mrn')
final_df$combined_id <- paste0(final_df$mrn, '_', final_df$record_id)
final_df <- final_df %>% filter(!duplicated(combined_id))


scrna_small <- subset(so_subset, record_id %in% final_df$record_id)
scrna_small <- subset(scrna_small,
                      subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

meta.data      <- scrna_small@meta.data
covar_lookup   <- final_df %>%
  dplyr::select(record_id, tzd, bmi, sex, epic_mfm_1, epic_glp1ra_1) %>%
  distinct(record_id, .keep_all = TRUE)

meta_joined <- meta.data %>%
  left_join(covar_lookup, by = 'record_id')

scrna_small$group_labels <- meta_joined$tzd
scrna_small$bmi          <- meta_joined$bmi.x
scrna_small$sex          <- meta_joined$sex.y
scrna_small$epic_mfm     <- meta_joined$epic_mfm.x
scrna_small$epic_glp1ra  <- meta_joined$epic_glp1ra.x

# ============================================================================
# POOLED OFFSET
# ============================================================================

counts_full <- round(GetAssayData(scrna_small, layer = 'counts'))
scrna_small$library_size <- Matrix::colSums(counts_full)
sce_full <- SingleCellExperiment(assays = list(counts = counts_full))
sce_full <- computeSumFactors(sce_full)
scrna_small$pooled_offset <- sizeFactors(sce_full)
remove(sce_full, counts_full)

# ============================================================================
# NEBULA SENSITIVITY FUNCTION
# ============================================================================

run_nebula_sensitivity <- function(so_obj, dir.results, celltype,
                                   covariate_name = NULL, model_label = "base") {
  if (celltype == 'All') {
    so_celltype <- so_obj
  } else if (celltype %in% c('TAL', 'EC', 'PT')) {
    so_celltype <- subset(so_obj, celltype2 == celltype)
  } else if (celltype == 'IC') {
    so_celltype <- subset(so_obj, celltype1 %in% c("B", "T", "MON", "MAC"))
  } else if (celltype == 'DCT') {
    so_celltype <- subset(so_obj, DCT_celltype == "DCT")
  } else if (celltype == 'POD') {
    so_celltype <- subset(so_obj, celltype1 == "POD")
  } else {
    if (celltype %in% unique(so_obj$celltype2)) {
      so_celltype <- subset(so_obj, celltype2 == celltype)
    } else if (celltype %in% unique(so_obj$celltype1)) {
      so_celltype <- subset(so_obj, celltype1 == celltype)
    } else {
      cat(paste0("Cell type '", celltype, "' not found. Skipping.\n"))
      return(invisible(NULL))
    }
  }
  
  DefaultAssay(so_celltype) <- "RNA"
  ct_clean <- str_replace_all(str_replace_all(celltype, "/", "_"), "-", "_")
  cat(paste0("\n=== ", ct_clean, " | model: ", model_label, " ===\n"))
  
  counts_mat <- round(GetAssayData(so_celltype, layer = "counts"))
  meta_gene  <- so_celltype@meta.data
  
  required_cols <- if (is.null(covariate_name)) "group_labels" else c("group_labels", covariate_name)
  complete_idx  <- complete.cases(meta_gene[, required_cols, drop = FALSE])
  cat("Cells with complete data:", sum(complete_idx), "\n")
  
  meta_gene  <- meta_gene[complete_idx, ]
  counts_mat <- counts_mat[, complete_idx]
  
  tmp_df  <- meta_gene %>% dplyr::select(record_id, group_labels) %>% distinct(record_id, .keep_all = TRUE)
  num_yes <- sum(tmp_df$group_labels == 'Yes', na.rm = TRUE)
  num_no  <- sum(tmp_df$group_labels == 'No',  na.rm = TRUE)
  cat("Participants - Yes:", num_yes, "No:", num_no, "\n")
  
  if (is.null(covariate_name)) {
    formula_str <- "~group_labels"
  } else {
    formula_str <- paste0("~group_labels + ", covariate_name)
  }
  pred_gene <- model.matrix(as.formula(formula_str), data = meta_gene)
  
  lib    <- meta_gene$pooled_offset
  data_g <- group_cell(count = counts_mat, id = meta_gene$record_id,
                       pred = pred_gene, offset = lib)
  if (is.null(data_g)) {
    data_g <- list(count = counts_mat, id = meta_gene$record_id,
                   pred = pred_gene, offset = lib)
  }
  
  result <- nebula(count = data_g$count, id = data_g$id,
                   pred = data_g$pred, ncore = 1, reml = TRUE,
                   model = "NBLMM", output_re = TRUE, covariance = TRUE,
                   offset = data_g$offset)
  
  full_results <- as.data.frame(result)
  full_results$num_cells   <- nrow(meta_gene)
  full_results$num_pio_yes <- num_yes
  full_results$num_pio_no  <- num_no
  full_results$model_label <- model_label
  full_results$covariate   <- ifelse(is.null(covariate_name), "none", covariate_name)
  
  out_file <- paste0(dir.results, 'NEBULA/NEBULA_', ct_clean, '_', model_label, '.csv')
  write.csv(full_results, out_file, row.names = FALSE)
  cat("Saved:", out_file, "\n")
  
  return(invisible(full_results))
}

# ============================================================================
# RUN ALL MODELS x ALL CELL TYPES
# ============================================================================

for (ct in celltypes_vec) {
  for (model_label in names(COVARIATES)) {
    run_nebula_sensitivity(
      so_obj         = scrna_small,
      dir.results    = dir.results,
      celltype       = ct,
      covariate_name = COVARIATES[[model_label]],
      model_label    = model_label
    )
  }
}

# ============================================================================
# HELPER: LOAD PIO EFFECT FROM NEBULA OUTPUT
# ============================================================================

load_pio_effect <- function(dir.results, ct_clean, model_label) {
  fp <- paste0(dir.results, 'NEBULA/NEBULA_', ct_clean, '_', model_label, '.csv')
  if (!file.exists(fp)) return(NULL)
  df <- read.csv(fp)
  df %>%
    dplyr::select(
      gene   = summary.gene,
      logFC  = summary.logFC_group_labelsYes,
      pvalue = summary.p_group_labelsYes
    ) %>%
    mutate(model = model_label, fdr = p.adjust(pvalue, method = "BH"))
}

# ============================================================================
# 1. LogFC CORRELATION MATRIX
# ============================================================================

cat("\n=== Building LogFC correlation matrices ===\n")

for (ct in celltypes_vec) {
  ct_clean  <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  model_dfs <- lapply(names(COVARIATES), function(m) load_pio_effect(dir.results, ct_clean, m))
  names(model_dfs) <- names(COVARIATES)
  model_dfs <- Filter(Negate(is.null), model_dfs)
  if (length(model_dfs) < 2) next
  
  wide <- Reduce(function(a, b) full_join(a, b, by = "gene"),
                 lapply(names(model_dfs), function(m) {
                   model_dfs[[m]] %>% dplyr::select(gene, logFC) %>%
                     rename_with(~ m, "logFC")
                 }))
  
  cor_mat <- cor(wide[, names(model_dfs)], use = "pairwise.complete.obs")
  
  png(paste0(dir.results, 'Comparison/Corr_LogFC_', ct_clean, '.png'),
      width = 900, height = 800, res = 150)
  pheatmap::pheatmap(
    cor_mat,
    display_numbers = TRUE, number_format = "%.3f",
    color  = colorRampPalette(c("white", "#4393c3"))(50),
    breaks = seq(0.8, 1, length.out = 51),
    main   = paste0(ct, ": PIO LogFC Correlation Across Models"),
    fontsize = 12
  )
  dev.off()
  cat("Saved correlation heatmap for", ct, "\n")
}

# ============================================================================
# 2. EFFECT STABILITY SCATTER PLOTS
# ============================================================================

cat("\n=== Building effect-stability scatter plots ===\n")

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  base_df  <- load_pio_effect(dir.results, ct_clean, "base")
  if (is.null(base_df)) next
  
  plot_list <- list()
  for (m in setdiff(names(COVARIATES), "base")) {
    cov_df <- load_pio_effect(dir.results, ct_clean, m)
    if (is.null(cov_df)) next
    
    merged <- inner_join(
      base_df %>% dplyr::select(gene, logFC_base = logFC, pval_base = pvalue),
      cov_df  %>% dplyr::select(gene, logFC_cov  = logFC, pval_cov  = pvalue),
      by = "gene"
    ) %>%
      mutate(
        status = case_when(
          pval_base < 0.05 & pval_cov < 0.05  ~ "Sig in both",
          pval_base < 0.05 & pval_cov >= 0.05 ~ "Lost after covariate",
          pval_base >= 0.05 & pval_cov < 0.05 ~ "Gained after covariate",
          TRUE                                 ~ "NS in both"
        ),
        delta = abs(logFC_cov - logFC_base),
        label = ifelse(rank(-delta) <= 8 & status != "NS in both", gene, NA)
      )
    
    r_val  <- round(cor(merged$logFC_base, merged$logFC_cov, use = "complete.obs"), 3)
    n_lost <- sum(merged$status == "Lost after covariate")
    n_gain <- sum(merged$status == "Gained after covariate")
    
    plot_list[[m]] <- ggplot(merged, aes(x = logFC_base, y = logFC_cov, color = status)) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
      ggrepel::geom_text_repel(aes(label = label), size = 2.5, color = "black", max.overlaps = 15) +
      scale_color_manual(values = c(
        "Sig in both"            = "#d73027",
        "Lost after covariate"   = "#fc8d59",
        "Gained after covariate" = "#4575b4",
        "NS in both"             = "grey80"
      )) +
      theme_bw() +
      labs(x = "LogFC (base model)", y = paste0("LogFC (+", m, ")"), color = NULL,
           title = paste0(ct, ": base vs +", m),
           subtitle = paste0("r = ", r_val, "  |  Lost: ", n_lost, "  |  Gained: ", n_gain)) +
      theme(plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 9),
            legend.position = "bottom")
  }
  
  if (length(plot_list) > 0) {
    combined <- wrap_plots(plot_list, ncol = 2) +
      plot_annotation(title = paste0(ct, " — PIO Effect: Base vs Covariate-Adjusted Models"),
                      theme = theme(plot.title = element_text(face = "bold", size = 14)))
    ggsave(paste0(dir.results, 'Comparison/Scatter_BaseVsCov_', ct_clean, '.pdf'), combined, width = 14, height = 10)
    ggsave(paste0(dir.results, 'Comparison/Scatter_BaseVsCov_', ct_clean, '.png'), combined, width = 14, height = 10, dpi = 200)
    cat("Saved scatter plots for", ct, "\n")
  }
}

# ============================================================================
# 3. SIGNIFICANCE CHANGE SUMMARY TABLE
# ============================================================================

cat("\n=== Building significance summary table ===\n")

summary_rows <- list()
for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  base_df  <- load_pio_effect(dir.results, ct_clean, "base")
  if (is.null(base_df)) next
  
  for (m in setdiff(names(COVARIATES), "base")) {
    cov_df <- load_pio_effect(dir.results, ct_clean, m)
    if (is.null(cov_df)) next
    
    merged <- inner_join(
      base_df %>% dplyr::select(gene, logFC_base = logFC, pval_base = pvalue),
      cov_df  %>% dplyr::select(gene, logFC_cov  = logFC, pval_cov  = pvalue),
      by = "gene"
    )
    
    summary_rows[[paste0(ct, "_", m)]] <- data.frame(
      celltype       = ct,
      covariate      = m,
      n_genes_tested = nrow(merged),
      n_sig_base     = sum(merged$pval_base < 0.05, na.rm = TRUE),
      n_sig_adjusted = sum(merged$pval_cov  < 0.05, na.rm = TRUE),
      n_sig_both     = sum(merged$pval_base < 0.05 & merged$pval_cov < 0.05, na.rm = TRUE),
      n_lost         = sum(merged$pval_base < 0.05 & merged$pval_cov >= 0.05, na.rm = TRUE),
      n_gained       = sum(merged$pval_base >= 0.05 & merged$pval_cov < 0.05, na.rm = TRUE),
      r_logFC        = round(cor(merged$logFC_base, merged$logFC_cov, use = "complete.obs"), 4)
    )
  }
}

summary_df <- bind_rows(summary_rows)
write.csv(summary_df, paste0(dir.results, 'Comparison/Significance_Change_Summary.csv'), row.names = FALSE)
print(summary_df)

# ============================================================================
# 4. VOLCANO SIDE-BY-SIDE
# ============================================================================

cat("\n=== Building side-by-side volcano plots ===\n")

for (ct in celltypes_vec) {
  ct_clean  <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  vol_list  <- list()
  
  for (m in names(COVARIATES)) {
    df <- load_pio_effect(dir.results, ct_clean, m)
    if (is.null(df)) next
    
    df <- df %>%
      mutate(diffexp = case_when(
        pvalue < 0.05 & logFC > 0 ~ "Up",
        pvalue < 0.05 & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )) %>%
      arrange(pvalue) %>%
      mutate(label = ifelse(row_number() <= 8, gene, NA)) %>%
      filter(abs(logFC) < 10)
    
    title_str <- if (m == "base") "Base model" else paste0("+", m)
    
    vol_list[[m]] <- ggplot(df, aes(x = logFC, y = -log10(pvalue), color = diffexp, label = label)) +
      geom_point(size = 1.2, alpha = 0.7) +
      ggrepel::geom_text_repel(size = 2.2, color = "black", max.overlaps = 12) +
      scale_color_manual(values = c(Down = "orange", NS = "grey70", Up = "purple"),
                         labels = c(Down = "Down in Pio", NS = "NS", Up = "Up in Pio")) +
      geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed", linewidth = 0.4) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
      theme_classic() +
      labs(x = "LogFC", y = "-log10(p)", color = NULL, title = paste0(ct, ": ", title_str)) +
      theme(plot.title = element_text(face = "bold", size = 10),
            legend.position = "bottom", aspect.ratio = 1)
  }
  
  if (length(vol_list) > 1) {
    combined_vol <- wrap_plots(vol_list, ncol = 3) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
    ggsave(paste0(dir.results, 'VolcanoPlots/Volcano_AllModels_', ct_clean, '.pdf'), combined_vol, width = 18, height = 6)
    ggsave(paste0(dir.results, 'VolcanoPlots/Volcano_AllModels_', ct_clean, '.png'), combined_vol, width = 18, height = 6, dpi = 200)
    cat("Saved volcano grid for", ct, "\n")
  }
}

# ============================================================================
# 5. ROBUSTLY SIGNIFICANT GENES
# ============================================================================

cat("\n=== Identifying robustly significant genes ===\n")

robust_results <- list()
for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  all_dfs  <- lapply(names(COVARIATES), function(m) load_pio_effect(dir.results, ct_clean, m))
  names(all_dfs) <- names(COVARIATES)
  all_dfs <- Filter(Negate(is.null), all_dfs)
  if (length(all_dfs) < 2) next
  
  gene_sig_counts <- Reduce(function(a, b) full_join(a, b, by = "gene"),
                            lapply(names(all_dfs), function(m) {
                              all_dfs[[m]] %>%
                                mutate(sig = as.integer(pvalue < 0.05)) %>%
                                dplyr::select(gene, sig) %>%
                                rename_with(~ paste0("sig_", m), "sig")
                            }))
  
  sig_cols <- paste0("sig_", names(all_dfs))
  gene_sig_counts$n_models_sig <- rowSums(gene_sig_counts[, sig_cols], na.rm = TRUE)
  
  robust <- gene_sig_counts %>%
    left_join(all_dfs[["base"]] %>% dplyr::select(gene, logFC_base = logFC, pval_base = pvalue), by = "gene") %>%
    arrange(desc(n_models_sig), pval_base) %>%
    mutate(celltype = ct)
  
  robust_results[[ct]] <- robust
  write.csv(robust, paste0(dir.results, 'Comparison/RobustGenes_', ct_clean, '.csv'), row.names = FALSE)
}

robust_combined <- bind_rows(robust_results)
write.csv(robust_combined, paste0(dir.results, 'Comparison/RobustGenes_AllCelltypes.csv'), row.names = FALSE)

cat("\n=== Sensitivity analysis complete! ===\n")
cat("Outputs saved to:", dir.results, "\n")














###GSEA

#### Pioglitazone scRNAseq - GSEA Sensitivity Analysis
#### Runs fgsea for each covariate model and compares pathway results
#### Assumes NEBULA sensitivity CSVs already exist from sensitivity script

library(tidyverse)
library(data.table)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(stringr)
library(ggrepel)
library(RColorBrewer)

# ============================================================================
# SETUP
# ============================================================================

dir.nebula  <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/Sensitivity/NEBULA/'
dir.results <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/Sensitivity/GSEA/'

dir.create(paste0(dir.results, 'DotPlots/'),   recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'Heatmaps/'),   recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'Comparison/'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'Raw/'),        recursive = TRUE, showWarnings = FALSE)

COVARIATES    <- c("base", "bmi", "sex", "metformin", "glp1ra")
celltypes_vec <- c('All', 'PT', 'TAL', 'EC', 'POD', 'DCT', 'IC')

GSEA_PARAMS <- list(minSize = 15, maxSize = 500, nPermSimple = 10000)
P_THRESH    <- 0.05   # nominal p-value threshold for "significant"
NES_THRESH  <- 0      # filter for directional plots (keep all if 0)

# ============================================================================
# GENE SETS
# ============================================================================

cat("Loading gene sets...\n")
hallmark_list <- split(msigdbr(species = "Homo sapiens", category = "H")$gene_symbol,
                       msigdbr(species = "Homo sapiens", category = "H")$gs_name)
go_bp_list    <- split(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")$gene_symbol,
                       msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")$gs_name)

geneset_types <- list(Hallmark = hallmark_list, GO_BP = go_bp_list)
# Add GO_CC / GO_MF if desired — kept to Hallmark + GO_BP for manageability

# ============================================================================
# HELPER: LOAD NEBULA RESULTS & BUILD RANKED GENE LIST
# ============================================================================

load_ranked_genes <- function(dir.nebula, ct_clean, model_label) {
  fp <- paste0(dir.nebula, 'NEBULA_', ct_clean, '_', model_label, '.csv')
  if (!file.exists(fp)) return(NULL)
  df <- read.csv(fp)
  ranked <- df %>%
    dplyr::select(gene = summary.gene, logFC = summary.logFC_group_labelsYes) %>%
    filter(!is.na(logFC)) %>%
    arrange(desc(logFC)) %>%
    { setNames(.$logFC, .$gene) }
  return(ranked)
}

# ============================================================================
# HELPER: CLEAN PATHWAY NAMES
# ============================================================================

clean_pathway <- function(x) {
  x %>%
    gsub("HALLMARK_|GOBP_|GOCC_|GOMF_", "", .) %>%
    gsub("_", " ", .) %>%
    tools::toTitleCase() %>%
    tolower() %>%
    tools::toTitleCase()
}

# ============================================================================
# RUN GSEA FOR ALL CELL TYPES x MODELS x GENE SETS
# ============================================================================

cat("\n=== Running GSEA for all models ===\n")
all_gsea <- list()  # named: celltype_model_geneset

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  for (model in COVARIATES) {
    ranked_genes <- load_ranked_genes(dir.nebula, ct_clean, model)
    if (is.null(ranked_genes)) {
      cat("  Missing:", ct_clean, model, "— skipping\n")
      next
    }
    
    for (gs_name in names(geneset_types)) {
      key <- paste(ct, model, gs_name, sep = "__")
      cat("  Running:", key, "\n")
      
      set.seed(42)
      res <- fgsea(
        pathways  = geneset_types[[gs_name]],
        stats     = ranked_genes,
        minSize   = GSEA_PARAMS$minSize,
        maxSize   = GSEA_PARAMS$maxSize,
        nPermSimple = GSEA_PARAMS$nPermSimple
      )
      
      res$celltype     <- ct
      res$model        <- model
      res$geneset_type <- gs_name
      res$pathway_clean <- clean_pathway(res$pathway)
      
      all_gsea[[key]] <- res
    }
  }
}

# Save raw results
combined_raw <- bind_rows(all_gsea) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(combined_raw, paste0(dir.results, 'Raw/all_gsea_sensitivity.csv'), row.names = FALSE)
cat("Raw GSEA saved.\n")

# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================

# Pull NES / pval for a given celltype + geneset, wide across models
make_wide_table <- function(all_gsea, ct, gs_name, value_col = "NES") {
  keys <- paste(ct, COVARIATES, gs_name, sep = "__")
  dfs  <- lapply(COVARIATES, function(m) {
    key <- paste(ct, m, gs_name, sep = "__")
    df  <- all_gsea[[key]]
    if (is.null(df)) return(NULL)
    df %>% dplyr::select(pathway, pathway_clean, !!value_col) %>%
      rename_with(~ m, all_of(value_col))
  })
  Reduce(function(a, b) full_join(a, b, by = c("pathway", "pathway_clean")),
         Filter(Negate(is.null), dfs))
}

# ============================================================================
# 1. NES HEATMAP: base vs covariate models, for pathways sig in base
# ============================================================================

cat("\n=== Building NES heatmaps ===\n")

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  for (gs_name in names(geneset_types)) {
    
    # Get pathways significant in the base model
    base_key <- paste(ct, "base", gs_name, sep = "__")
    base_df  <- all_gsea[[base_key]]
    if (is.null(base_df)) next
    
    sig_paths <- base_df %>% filter(pval < P_THRESH) %>% pull(pathway)
    if (length(sig_paths) == 0) {
      cat("  No significant pathways in base for", ct, gs_name, "— skipping heatmap\n")
      next
    }
    
    # Wide NES table
    wide_nes  <- make_wide_table(all_gsea, ct, gs_name, "NES")
    wide_pval <- make_wide_table(all_gsea, ct, gs_name, "pval")
    
    mat <- wide_nes %>%
      filter(pathway %in% sig_paths) %>%
      column_to_rownames("pathway_clean") %>%
      dplyr::select(all_of(COVARIATES)) %>%
      as.matrix()
    
    if (nrow(mat) == 0) next
    
    # Significance annotation (asterisks) from pval matrix
    pval_mat <- wide_pval %>%
      filter(pathway %in% sig_paths) %>%
      column_to_rownames("pathway_clean") %>%
      dplyr::select(all_of(COVARIATES)) %>%
      as.matrix()
    
    # Cell labels: NES rounded, asterisk if p<0.05
    cell_labels <- matrix(
      ifelse(pval_mat < 0.05,
             paste0(round(mat, 2), "*"),
             round(mat, 2)),
      nrow = nrow(mat), ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
    
    col_lim <- max(abs(mat), na.rm = TRUE)
    
    png(paste0(dir.results, 'Heatmaps/NES_Heatmap_', ct_clean, '_', gs_name, '.png'),
        width = 1000, height = max(400, nrow(mat) * 22 + 200), res = 120)
    pheatmap::pheatmap(
      mat,
      display_numbers = cell_labels,
      number_fontsize = 7,
      color  = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
      breaks = seq(-col_lim, col_lim, length.out = 101),
      cluster_cols = FALSE,
      cluster_rows = nrow(mat) > 1,
      border_color = "white",
      main = paste0(ct, " — ", gs_name, "\nNES across models (* = p<0.05)"),
      fontsize_row = 8,
      fontsize_col = 10
    )
    dev.off()
    cat("  Saved NES heatmap:", ct, gs_name, "\n")
  }
}

# ============================================================================
# 2. NES SCATTER: base NES vs covariate-adjusted NES (per pathway)
# ============================================================================

cat("\n=== Building NES scatter plots ===\n")

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  for (gs_name in names(geneset_types)) {
    
    base_key <- paste(ct, "base", gs_name, sep = "__")
    base_df  <- all_gsea[[base_key]]
    if (is.null(base_df)) next
    
    cov_models <- setdiff(COVARIATES, "base")
    plot_list  <- list()
    
    for (m in cov_models) {
      cov_key <- paste(ct, m, gs_name, sep = "__")
      cov_df  <- all_gsea[[cov_key]]
      if (is.null(cov_df)) next
      
      merged <- inner_join(
        base_df %>% dplyr::select(pathway, pathway_clean, NES_base = NES, pval_base = pval),
        cov_df  %>% dplyr::select(pathway, NES_cov = NES, pval_cov = pval),
        by = "pathway"
      ) %>%
        mutate(
          status = case_when(
            pval_base < P_THRESH & pval_cov < P_THRESH ~ "Sig in both",
            pval_base < P_THRESH & pval_cov >= P_THRESH ~ "Lost",
            pval_base >= P_THRESH & pval_cov < P_THRESH ~ "Gained",
            TRUE ~ "NS in both"
          ),
          delta  = abs(NES_cov - NES_base),
          label  = ifelse(
            (status %in% c("Sig in both", "Lost", "Gained")) & rank(-delta) <= 10,
            pathway_clean, NA
          )
        )
      
      r_val  <- round(cor(merged$NES_base, merged$NES_cov, use = "complete.obs"), 3)
      n_lost  <- sum(merged$status == "Lost")
      n_gained <- sum(merged$status == "Gained")
      
      p <- ggplot(merged, aes(x = NES_base, y = NES_cov, color = status, label = label)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
        ggrepel::geom_text_repel(size = 2.3, color = "black", max.overlaps = 12) +
        scale_color_manual(values = c(
          "Sig in both" = "#d73027",
          "Lost"        = "#fc8d59",
          "Gained"      = "#4575b4",
          "NS in both"  = "grey80"
        )) +
        theme_bw() +
        labs(
          x     = "NES (base model)",
          y     = paste0("NES (+", m, ")"),
          color = NULL,
          title = paste0("+", m),
          subtitle = paste0("r = ", r_val, "  |  Lost: ", n_lost, "  |  Gained: ", n_gained)
        ) +
        theme(plot.title = element_text(face = "bold", size = 10),
              plot.subtitle = element_text(size = 8),
              legend.position = "bottom",
              aspect.ratio = 1)
      
      plot_list[[m]] <- p
    }
    
    if (length(plot_list) > 0) {
      combined <- wrap_plots(plot_list, ncol = 2) +
        plot_annotation(
          title = paste0(ct, " — ", gs_name, ": NES Stability Across Covariate Models"),
          theme = theme(plot.title = element_text(face = "bold", size = 13))
        ) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
      
      ggsave(paste0(dir.results, 'Comparison/NES_Scatter_', ct_clean, '_', gs_name, '.pdf'),
             combined, width = 12, height = 10)
      ggsave(paste0(dir.results, 'Comparison/NES_Scatter_', ct_clean, '_', gs_name, '.png'),
             combined, width = 12, height = 10, dpi = 200)
      cat("  Saved NES scatter:", ct, gs_name, "\n")
    }
  }
}

# ============================================================================
# 3. ROBUST PATHWAY TABLE
#    Pathways significant in base + how many covariate models they survive
# ============================================================================

cat("\n=== Identifying robust pathways ===\n")

robust_pathway_rows <- list()

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  for (gs_name in names(geneset_types)) {
    
    # Collect sig flags across all models
    sig_wide <- lapply(COVARIATES, function(m) {
      key <- paste(ct, m, gs_name, sep = "__")
      df  <- all_gsea[[key]]
      if (is.null(df)) return(NULL)
      df %>%
        mutate(sig = as.integer(pval < P_THRESH)) %>%
        dplyr::select(pathway, pathway_clean, NES, pval, sig) %>%
        rename_with(~ paste0(c("NES", "pval", "sig"), "_", m),
                    c("NES", "pval", "sig"))
    })
    sig_wide <- Filter(Negate(is.null), sig_wide)
    if (length(sig_wide) == 0) next
    
    merged_wide <- Reduce(function(a, b) full_join(a, b, by = c("pathway", "pathway_clean")),
                          sig_wide)
    
    sig_cols <- paste0("sig_", COVARIATES)
    sig_cols_present <- intersect(sig_cols, names(merged_wide))
    merged_wide$n_models_sig <- rowSums(merged_wide[, sig_cols_present], na.rm = TRUE)
    
    # NES direction consistency: flag if any covariate model flips sign
    nes_cols <- paste0("NES_", COVARIATES)
    nes_cols_present <- intersect(nes_cols, names(merged_wide))
    if (length(nes_cols_present) > 1) {
      nes_matrix <- merged_wide[, nes_cols_present]
      merged_wide$direction_consistent <- apply(nes_matrix, 1, function(x) {
        x <- na.omit(x)
        if (length(x) == 0) return(NA)
        all(x > 0) | all(x < 0)
      })
    }
    
    robust <- merged_wide %>%
      filter(!is.na(NES_base)) %>%  # must have base model result
      arrange(desc(n_models_sig), pval_base) %>%
      mutate(celltype = ct, geneset_type = gs_name)
    
    key_out <- paste0(ct, "__", gs_name)
    robust_pathway_rows[[key_out]] <- robust
  }
}

robust_all <- bind_rows(robust_pathway_rows)
write.csv(robust_all, paste0(dir.results, 'Comparison/RobustPathways_All.csv'), row.names = FALSE)

# Focus table: sig in base AND all covariate models
fully_robust <- robust_all %>%
  filter(n_models_sig == length(COVARIATES)) %>%
  dplyr::select(celltype, geneset_type, pathway_clean,
                NES_base, pval_base, n_models_sig,
                any_of("direction_consistent")) %>%
  arrange(celltype, geneset_type, pval_base)

write.csv(fully_robust, paste0(dir.results, 'Comparison/FullyRobustPathways.csv'), row.names = FALSE)
cat("Fully robust pathways (sig in all models):", nrow(fully_robust), "\n")
print(fully_robust %>% dplyr::select(celltype, geneset_type, pathway_clean, NES_base, pval_base))

# ============================================================================
# 4. SIGNIFICANCE SUMMARY TABLE (n sig pathways per model)
# ============================================================================

cat("\n=== Building pathway count summary ===\n")

pathway_summary <- lapply(celltypes_vec, function(ct) {
  lapply(names(geneset_types), function(gs_name) {
    lapply(COVARIATES, function(m) {
      key <- paste(ct, m, gs_name, sep = "__")
      df  <- all_gsea[[key]]
      if (is.null(df)) return(NULL)
      data.frame(
        celltype     = ct,
        geneset_type = gs_name,
        model        = m,
        n_sig        = sum(df$pval < P_THRESH, na.rm = TRUE),
        n_sig_up     = sum(df$pval < P_THRESH & df$NES > 0, na.rm = TRUE),
        n_sig_down   = sum(df$pval < P_THRESH & df$NES < 0, na.rm = TRUE)
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()

write.csv(pathway_summary, paste0(dir.results, 'Comparison/PathwayCount_Summary.csv'), row.names = FALSE)

# Print as wide pivot
pathway_summary_wide <- pathway_summary %>%
  pivot_wider(names_from = model, values_from = c(n_sig, n_sig_up, n_sig_down))
print(pathway_summary_wide)

# ============================================================================
# 5. DOT PLOT: top pathways across models for a given cell type
#    One panel per covariate model; rows = top 20 pathways from base
# ============================================================================

cat("\n=== Building dot plots ===\n")

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  for (gs_name in names(geneset_types)) {
    
    base_key <- paste(ct, "base", gs_name, sep = "__")
    base_df  <- all_gsea[[base_key]]
    if (is.null(base_df)) next
    
    top_paths <- base_df %>%
      filter(pval < P_THRESH) %>%
      arrange(pval) %>%
      slice_head(n = 20) %>%
      pull(pathway)
    
    if (length(top_paths) == 0) {
      # Fall back to top 20 by p-value even if none pass threshold
      top_paths <- base_df %>% arrange(pval) %>% slice_head(n = 20) %>% pull(pathway)
      if (length(top_paths) == 0) next
    }
    
    # Collect NES + pval for all models, filtered to top_paths
    plot_df <- lapply(COVARIATES, function(m) {
      key <- paste(ct, m, gs_name, sep = "__")
      df  <- all_gsea[[key]]
      if (is.null(df)) return(NULL)
      df %>%
        filter(pathway %in% top_paths) %>%
        dplyr::select(pathway, pathway_clean, NES, pval) %>%
        mutate(model = m)
    }) %>%
      bind_rows() %>%
      mutate(
        model = factor(model, levels = COVARIATES),
        pathway_clean = factor(pathway_clean,
                               levels = base_df %>%
                                 filter(pathway %in% top_paths) %>%
                                 arrange(NES) %>%
                                 pull(pathway_clean))
      )
    
    p <- ggplot(plot_df, aes(x = model, y = pathway_clean,
                             color = NES, size = -log10(pval))) +
      geom_point() +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                            name = "NES") +
      scale_size_continuous(range = c(1, 8), name = "-log10(p)") +
      geom_point(data = plot_df %>% filter(pval < P_THRESH),
                 aes(x = model, y = pathway_clean),
                 shape = 8, color = "black", size = 1.5) +
      theme_bw() +
      theme(
        axis.text.x  = element_text(angle = 30, hjust = 1, size = 10),
        axis.text.y  = element_text(size = 8),
        plot.title   = element_text(face = "bold", size = 11),
        panel.grid   = element_line(color = "grey90")
      ) +
      labs(
        x     = "Model",
        y     = NULL,
        title = paste0(ct, " — ", gs_name,
                       "\nTop base-model pathways across covariate adjustments"),
        caption = "★ = p < 0.05"
      )
    
    h <- max(5, length(top_paths) * 0.35 + 2)
    ggsave(paste0(dir.results, 'DotPlots/DotPlot_', ct_clean, '_', gs_name, '.pdf'),
           p, width = 9, height = h)
    ggsave(paste0(dir.results, 'DotPlots/DotPlot_', ct_clean, '_', gs_name, '.png'),
           p, width = 9, height = h, dpi = 200)
    cat("  Saved dot plot:", ct, gs_name, "\n")
  }
}

cat("\n=== GSEA sensitivity analysis complete! ===\n")
cat("Output directory:", dir.results, "\n\n")
cat("Key outputs:\n")
cat("  Raw/all_gsea_sensitivity.csv         — all fgsea results\n")
cat("  Heatmaps/NES_Heatmap_*              — NES across models for base-sig pathways\n")
cat("  Comparison/NES_Scatter_*            — NES stability scatter (base vs each covariate)\n")
cat("  DotPlots/DotPlot_*                  — dot plot of top pathways across all models\n")
cat("  Comparison/PathwayCount_Summary.csv — n sig pathways per model\n")
cat("  Comparison/RobustPathways_All.csv   — all pathways with n_models_sig count\n")
cat("  Comparison/FullyRobustPathways.csv  — pathways sig in ALL models\n")

