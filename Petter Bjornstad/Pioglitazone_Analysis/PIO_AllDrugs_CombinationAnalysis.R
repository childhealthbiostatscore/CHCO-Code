#### ============================================================
#### Multi-Drug & Combination GSEA – Endothelial Cells
#### Drugs: Pioglitazone (TZD), GLP-1RA, SGLT2i, Metformin
#### Compares: each drug alone + all pairwise combinations
#### ============================================================

library(Seurat); library(nebula); library(Matrix)
library(SingleCellExperiment); library(scran)
library(dplyr); library(tidyr); library(stringr)
library(ggplot2); library(ggrepel); library(patchwork)
library(fgsea); library(msigdbr)
library(readxl); library(pheatmap)
library(RColorBrewer); library(viridis)
library(ComplexHeatmap); library(circlize)

# ── Paths ────────────────────────────────────────────────────────────────────
setwd("C:/Users/netio/Documents/Harmonized_data/")
dir.base <- "C:/Users/netio/Documents/UofW/Projects/drug_combinations_EC/"
for (d in c("NEBULA/", "GSEA/raw/", "Figures/Heatmaps/", "Figures/Bubbles/",
            "Figures/NES_Matrices/", "Figures/Synergy/", "Figures/Volcanos/")) {
  dir.create(paste0(dir.base, d), recursive = TRUE, showWarnings = FALSE)
}

# ── Load Seurat object ────────────────────────────────────────────────────────
load("C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData")
seu <- so_kpmp_sc
rm(so_kpmp_sc)

# ── Cell type definitions ─────────────────────────────────────────────────────
seu$celltype1 <- case_when(
  grepl("PT-",  seu$celltype_rpca) ~ "PT",
  grepl("TAL-", seu$celltype_rpca) ~ "TAL",
  grepl("EC-",  seu$celltype_rpca) ~ "EC",
  grepl("POD",  seu$celltype_rpca) ~ "POD",
  grepl("MAC",  seu$celltype_rpca) ~ "MAC",
  grepl("MON",  seu$celltype_rpca) ~ "MON",
  grepl("PC-",  seu$celltype_rpca) ~ "PC",
  grepl("FIB",  seu$celltype_rpca) ~ "FIB_MC_VSMC",
  grepl("DTL",  seu$celltype_rpca) ~ "DTL",
  seu$celltype_rpca == "DCT"       ~ "DCT",
  seu$celltype_rpca == "ATL"       ~ "ATL",
  seu$celltype_rpca == "B"         ~ "B",
  seu$celltype_rpca == "T"         ~ "T"
)
seu$celltype1 <- as.character(seu$celltype1)
seu$KPMP_celltype2 <- as.character(seu$KPMP_celltype)
seu$celltype2 <- ifelse(
  grepl("PT",  seu$KPMP_celltype), "PT",
  ifelse(grepl("TAL", seu$KPMP_celltype), "TAL",
         ifelse(grepl("EC-", seu$KPMP_celltype), "EC", seu$KPMP_celltype2)))
seu$DCT_celltype <- ifelse(
  seu$KPMP_celltype %in% c("DCT", "dDCT"), "DCT", "Non-DCT")

# ── Filter to T2D, remove CRC-55 ─────────────────────────────────────────────
seu <- subset(seu, subset = record_id != "CRC-55")
seu <- subset(seu, subset = group == "Type_2_Diabetes")
dir.base <- "C:/Users/netio/Documents/UofW/Projects/drug_combinations_EC/"
cat("Cells after filtering:", ncol(seu), "\n")
cat("Donors:", length(unique(seu$record_id)), "\n")

# ── Compute pooled size factors ───────────────────────────────────────────────
sce_full <- as.SingleCellExperiment(seu)
sce_full <- computeSumFactors(sce_full, clusters = seu$celltype2)
seu$pooled_offset <- log(sizeFactors(sce_full))

# ══════════════════════════════════════════════════════════════════════════════
# MEDICATIONS — load, strip _1, map columns, merge into Seurat metadata
# ══════════════════════════════════════════════════════════════════════════════
meds_raw <- readxl::read_xlsx("Biopsies_w_mrn_Oct3.xlsx")

# Check available columns
med_cols_available <- names(meds_raw) %>% str_subset("tzd|glp1|sglt|mfm|metf")
cat("Available medication columns:\n"); print(med_cols_available)

# Select _1 columns only, exclude ever_ prefix, then strip _1 suffix
meds <- meds_raw %>%
  dplyr::select(mrn, ends_with("_1"), -starts_with("ever_"))
names(meds) <- gsub("_1$", "", names(meds))
cat("Stripped med columns:\n")
print(names(meds) %>% str_subset("tzd|glp1|sglt|mfm|metf"))

# Drug column mapping (no epic_ prefix in this file)
DRUG_COLS <- list(
  TZD   = "tzd",
  GLP1  = "glp1ra",
  SGLT2 = "sglti2",
  MFM   = "mfm"
)
DRUG_LABELS <- list(
  TZD   = "Pioglitazone (TZD)",
  GLP1  = "GLP-1RA",
  SGLT2 = "SGLT2i",
  MFM   = "Metformin"
)

# Safety check
missing_cols <- setdiff(unlist(DRUG_COLS), names(meds))
if (length(missing_cols) > 0) {
  stop("Missing medication columns: ", paste(missing_cols, collapse = ", "),
       "\nAvailable: ", paste(names(meds), collapse = ", "))
}

# Recode to Yes/No
drug_col_names <- unlist(DRUG_COLS, use.names = FALSE)
meds_slim <- meds %>%
  dplyr::select(mrn, all_of(drug_col_names)) %>%
  mutate(across(
    all_of(drug_col_names),
    ~ ifelse(. %in% c("Yes", 1, "1", TRUE), "Yes", "No")
  ))
cat("meds_slim created. Rows:", nrow(meds_slim), "\n")

# Merge into Seurat metadata
meds_slim <- meds_slim %>%
  mutate(mrn = as.character(mrn)) %>%
  filter(!is.na(mrn)) %>%
  distinct(mrn, .keep_all = TRUE)
meta <- seu@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  mutate(mrn = as.character(mrn)) %>%
  left_join(meds_slim, by = "mrn")
rownames(meta) <- meta$cell_barcode
seu@meta.data <- meta[colnames(seu), ]
cat("Meds merged into Seurat. Columns added:", paste(unlist(DRUG_COLS), collapse = ", "), "\n")

# ── Subset to Endothelial Cells only ─────────────────────────────────────────
cat("celltype2 values in seu:\n")
print(table(seu$celltype2))
cat("Drug column distributions in full seu:\n")
print(table(seu@meta.data[["tzd"]], useNA = "always"))

seu_ec <- subset(seu, celltype2 == "EC")
cat("EC cells:", ncol(seu_ec), "| donors:", length(unique(seu_ec$record_id)), "\n")

cat("Drug column distributions in EC subset:\n")
for (col in unlist(DRUG_COLS)) {
  cat(col, ":\n"); print(table(seu_ec@meta.data[[col]], useNA = "always"))
}

# ══════════════════════════════════════════════════════════════════════════════
# GENE SETS
# ══════════════════════════════════════════════════════════════════════════════
get_genesets <- function() {
  h  <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(.$gs_name) %>% lapply(function(x) x$gene_symbol)
  bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(.$gs_name) %>% lapply(function(x) x$gene_symbol)
  list(Hallmark = h, GOBP = bp)
}
cat("Loading gene sets...\n")
genesets <- get_genesets()

# ══════════════════════════════════════════════════════════════════════════════
# NEBULA HELPER
# ══════════════════════════════════════════════════════════════════════════════
run_nebula_binary <- function(seu_sub, drug_col, label,
                              covars = c("sex", "age", "hba1c")) {
  cat("  Running NEBULA for:", label, "\n")
  keep <- complete.cases(seu_sub@meta.data[, c(drug_col, covars)])
  seu_sub <- seu_sub[, keep]
  if (sum(seu_sub@meta.data[[drug_col]] == "Yes") < 3) {
    warning("  Fewer than 3 Yes-donors for ", label, " — skipping")
    return(NULL)
  }
  seu_sub@meta.data[[drug_col]] <- factor(seu_sub@meta.data[[drug_col]],
                                          levels = c("No", "Yes"))
  sce_sub <- as.SingleCellExperiment(seu_sub)
  df_neb <- seu_sub@meta.data %>%
    dplyr::select(record_id, all_of(c(drug_col, covars)), pooled_offset) %>%
    distinct()
  grp <- group_cell(count  = assay(sce_sub, "counts"),
                    id     = seu_sub$record_id,
                    pred   = df_neb,
                    offset = df_neb$pooled_offset)
  res <- nebula(grp$count, grp$id, pred = grp$pred, offset = grp$offset,
                model = "NBLMM", ncore = 4)$summary %>%
    as_tibble() %>%
    rename(gene = 1) %>%
    rename_with(~ str_replace(., paste0("logFC_", drug_col, "Yes"), "logFC"),
                starts_with("logFC_")) %>%
    rename_with(~ str_replace(., paste0("p_", drug_col, "Yes"), "pval"),
                starts_with("p_")) %>%
    filter(!is.na(logFC), !is.na(pval)) %>%
    mutate(fdr = p.adjust(pval, "fdr"), drug = label, drug_col = drug_col)
  res
}

# ══════════════════════════════════════════════════════════════════════════════
# GSEA HELPER
# ══════════════════════════════════════════════════════════════════════════════
run_gsea <- function(de_res, label, gs_list, gs_name) {
  ranked <- de_res %>%
    arrange(desc(logFC)) %>%
    dplyr::select(gene, logFC) %>%
    deframe()
  ranked <- ranked[!duplicated(names(ranked))]
  set.seed(42)
  fgsea(pathways = gs_list, stats = ranked,
        minSize = 15, maxSize = 500, nPerm = 10000) %>%
    as_tibble() %>%
    mutate(drug    = label,
           gs_type = gs_name,
           leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
}

# ══════════════════════════════════════════════════════════════════════════════
# DEFINE ALL CONDITIONS: single drugs + pairwise combinations
# ══════════════════════════════════════════════════════════════════════════════
drug_names <- names(DRUG_COLS)

single_conditions <- lapply(drug_names, function(d) {
  list(label = d, drugs = d, type = "single")
})
combo_conditions <- combn(drug_names, 2, simplify = FALSE) %>%
  lapply(function(pair) {
    list(label = paste(pair, collapse = "+"), drugs = pair, type = "combo")
  })
all_conditions <- c(single_conditions, combo_conditions)
cat("Conditions to test:", length(all_conditions), "\n")
cat(sapply(all_conditions, `[[`, "label"), sep = "\n")

# ══════════════════════════════════════════════════════════════════════════════
# RUN NEBULA FOR EACH CONDITION
# ══════════════════════════════════════════════════════════════════════════════
all_de <- list()

for (cond in all_conditions) {
  lbl  <- cond$label
  drgs <- cond$drugs
  
  if (cond$type == "single") {
    dcol <- DRUG_COLS[[drgs]]
    de   <- run_nebula_binary(seu_ec, dcol, lbl)
    if (!is.null(de)) all_de[[lbl]] <- de
    
  } else {
    cols    <- unlist(DRUG_COLS[drgs])
    tmp_col <- paste0("combo_", paste(drgs, collapse = "_"))
    seu_ec@meta.data[[tmp_col]] <- apply(
      seu_ec@meta.data[, cols, drop = FALSE], 1,
      function(row) {
        if (all(row == "Yes"))  "Yes"
        else if (all(row == "No")) "No"
        else NA_character_
      }
    )
    de <- run_nebula_binary(seu_ec, tmp_col, lbl)
    if (!is.null(de)) all_de[[lbl]] <- de
    seu_ec@meta.data[[tmp_col]] <- NULL
  }
  
  if (!is.null(all_de[[lbl]])) {
    write.csv(all_de[[lbl]],
              paste0(dir.base, "NEBULA/DE_EC_", lbl, ".csv"), row.names = FALSE)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# GSEA HELPER
# ══════════════════════════════════════════════════════════════════════════════
run_gsea <- function(de_res, label, gs_list, gs_name) {
  ranked <- de_res %>%
    arrange(desc(logFC)) %>%
    dplyr::select(gene, logFC) %>%
    deframe()
  ranked <- ranked[!duplicated(names(ranked))]
  set.seed(42)
  fgsea(pathways = gs_list, stats = ranked,
        minSize = 15, maxSize = 500, nPermSimple = 10000) %>%
    as_tibble() %>%
    mutate(drug    = label,
           gs_type = gs_name,
           leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN GSEA FOR EACH CONDITION x CELL TYPE x GENE SET
# ══════════════════════════════════════════════════════════════════════════════
all_gsea <- list()

for (key in names(all_de)) {
  ct  <- str_split(key, "__")[[1]][1]
  lbl <- str_split(key, "__")[[1]][2]
  de  <- all_de[[key]]
  cat("GSEA for:", key, "\n")
  for (gs_name in names(genesets)) {
    gkey <- paste(key, gs_name, sep = "|")
    gsea_res <- run_gsea(de, lbl, genesets[[gs_name]], gs_name)
    gsea_res$celltype <- ct
    all_gsea[[gkey]] <- gsea_res
  }
}

gsea_df <- bind_rows(all_gsea) %>%
  mutate(condition_type = if_else(str_detect(drug, "\\+"), "combo", "single"))

write.csv(gsea_df, paste0(dir.base, "GSEA/all_conditions_gsea.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# SYNERGY SCORE
# Synergy = NES(combo) - max(NES(drug_A), NES(drug_B))
# ══════════════════════════════════════════════════════════════════════════════
compute_synergy <- function(gsea_df, gs_type_sel = "Hallmark") {
  dat <- gsea_df %>% filter(gs_type == gs_type_sel)
  single_nes <- dat %>%
    filter(condition_type == "single") %>%
    dplyr::select(pathway, drug, NES) %>%
    rename(single_drug = drug, single_NES = NES)
  combo_nes <- dat %>%
    filter(condition_type == "combo") %>%
    dplyr::select(pathway, drug, NES, pval) %>%
    rename(combo_label = drug, combo_NES = NES, combo_pval = pval)
  
  combo_nes %>%
    mutate(d1 = str_split(combo_label, "\\+", simplify = TRUE)[,1],
           d2 = str_split(combo_label, "\\+", simplify = TRUE)[,2]) %>%
    left_join(single_nes %>% rename(NES_d1 = single_NES),
              by = c("pathway", "d1" = "single_drug")) %>%
    left_join(single_nes %>% rename(NES_d2 = single_NES),
              by = c("pathway", "d2" = "single_drug")) %>%
    mutate(
      max_single_NES = pmax(NES_d1, NES_d2, na.rm = TRUE),
      synergy_score  = combo_NES - max_single_NES,
      synergy_class  = case_when(
        synergy_score >  0.5 ~ "Synergistic",
        synergy_score < -0.5 ~ "Antagonistic",
        TRUE                 ~ "Additive"
      )
    )
}

synergy_h   <- compute_synergy(gsea_df, "Hallmark")
synergy_bp  <- compute_synergy(gsea_df, "GOBP")
synergy_all <- bind_rows(synergy_h, synergy_bp)
write.csv(synergy_all, paste0(dir.base, "GSEA/synergy_scores_EC.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATIONS
# ══════════════════════════════════════════════════════════════════════════════
clean_path <- function(x) {
  x %>% str_remove("^HALLMARK_|^GOBP_|^GO_") %>%
    str_replace_all("_", " ") %>%
    tools::toTitleCase() %>%
    str_trunc(55)
}

drug_colors <- c(
  TZD          = "#E41A1C",
  GLP1         = "#377EB8",
  SGLT2        = "#4DAF4A",
  MFM          = "#984EA3",
  "TZD+GLP1"   = "#FF7F00",
  "TZD+SGLT2"  = "#A65628",
  "TZD+MFM"    = "#F781BF",
  "GLP1+SGLT2" = "#999999",
  "GLP1+MFM"   = "#66C2A5",
  "SGLT2+MFM"  = "#FC8D62"
)

# ── 5A. NES Heatmap ───────────────────────────────────────────────────────────
plot_nes_heatmap <- function(gsea_df, gs_sel = "Hallmark", n_paths = 30,
                             fdr_thresh = 0.25) {
  dat <- gsea_df %>%
    filter(gs_type == gs_sel, padj < fdr_thresh) %>%
    mutate(pathway_clean = clean_path(pathway))
  top_paths <- dat %>%
    group_by(pathway) %>%
    summarise(var_NES = var(NES, na.rm = TRUE), .groups = "drop") %>%
    slice_max(var_NES, n = n_paths) %>% pull(pathway)
  mat <- dat %>%
    filter(pathway %in% top_paths) %>%
    dplyr::select(pathway_clean, drug, NES) %>%
    pivot_wider(names_from = drug, values_from = NES, values_fill = 0) %>%
    tibble::column_to_rownames("pathway_clean") %>%
    as.matrix()
  singles <- intersect(drug_names, colnames(mat))
  combos  <- setdiff(colnames(mat), singles)
  mat     <- mat[, c(singles, combos), drop = FALSE]
  col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#D6604D"))
  col_annot <- HeatmapAnnotation(
    Type = if_else(colnames(mat) %in% drug_names, "Single", "Combination"),
    col  = list(Type = c(Single = "#555555", Combination = "#DDAA33"))
  )
  ht <- Heatmap(mat, name = "NES", col = col_fun,
                top_annotation  = col_annot,
                cluster_rows    = TRUE, cluster_columns = FALSE,
                row_names_gp    = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_title    = paste0(gs_sel, " – NES across Drug Conditions (EC)"),
                column_title_gp = gpar(fontsize = 13, fontface = "bold")
  )
  pdf(paste0(dir.base, "Figures/Heatmaps/NES_heatmap_", gs_sel, "_EC.pdf"),
      width = 14, height = 12)
  draw(ht); dev.off()
  cat("  Saved NES heatmap:", gs_sel, "\n")
}

plot_nes_heatmap(gsea_df, "Hallmark")
plot_nes_heatmap(gsea_df, "GOBP", n_paths = 40)

# ── 5B. Bubble Plot ───────────────────────────────────────────────────────────
plot_bubble_grid <- function(gsea_df, gs_sel = "Hallmark", n_top = 15,
                             fdr_thresh = 0.25) {
  dat <- gsea_df %>%
    filter(gs_type == gs_sel) %>%
    mutate(pathway_clean = clean_path(pathway))
  top_paths <- dat %>%
    filter(padj < fdr_thresh) %>%
    group_by(drug) %>%
    slice_max(abs(NES), n = n_top) %>%
    pull(pathway) %>% unique()
  plot_dat <- dat %>%
    filter(pathway %in% top_paths) %>%
    mutate(
      drug = factor(drug, levels = c(drug_names,
                                     combn(drug_names, 2, paste, collapse = "+"))),
      sig  = padj < fdr_thresh
    )
  p <- ggplot(plot_dat, aes(x = drug, y = reorder(pathway_clean, NES),
                            color = NES, size = -log10(pmax(pval, 1e-10)))) +
    geom_point(aes(alpha = sig)) +
    scale_alpha_manual(values = c("FALSE" = 0.2, "TRUE" = 0.9), guide = "none") +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                          midpoint = 0, name = "NES") +
    scale_size_continuous(range = c(1, 8), name = "-log10(p)") +
    geom_vline(xintercept = length(drug_names) + 0.5,
               linetype = "dashed", color = "grey40") +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 40, hjust = 1, face = "bold"),
          axis.text.y  = element_text(size = 8),
          panel.grid.major.x = element_blank(),
          plot.title   = element_text(face = "bold", hjust = 0.5)) +
    labs(x = NULL, y = NULL,
         title = paste0(gs_sel, ": Drug Effects in Endothelial Cells"))
  ggsave(paste0(dir.base, "Figures/Bubbles/bubble_", gs_sel, "_EC.pdf"),
         p, width = 16, height = 14)
  ggsave(paste0(dir.base, "Figures/Bubbles/bubble_", gs_sel, "_EC.png"),
         p, width = 16, height = 14, dpi = 300)
  cat("  Saved bubble plot:", gs_sel, "\n")
}

plot_bubble_grid(gsea_df, "Hallmark")
plot_bubble_grid(gsea_df, "GOBP", n_top = 20)

# ── 5C. Synergy Heatmap ───────────────────────────────────────────────────────
plot_synergy_heatmap <- function(synergy_df, gs_sel = "Hallmark",
                                 n_paths = 30, fdr_thresh = 0.25) {
  dat <- synergy_df %>%
    filter(gs_type == gs_sel, combo_pval < fdr_thresh) %>%
    mutate(pathway_clean = clean_path(pathway))
  top_paths <- dat %>%
    group_by(pathway) %>%
    summarise(max_syn = max(abs(synergy_score), na.rm = TRUE), .groups = "drop") %>%
    slice_max(max_syn, n = n_paths) %>% pull(pathway)
  mat <- dat %>%
    filter(pathway %in% top_paths) %>%
    dplyr::select(pathway_clean, combo_label, synergy_score) %>%
    pivot_wider(names_from = combo_label, values_from = synergy_score,
                values_fill = 0) %>%
    tibble::column_to_rownames("pathway_clean") %>%
    as.matrix()
  lim     <- max(abs(mat), na.rm = TRUE)
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#4393C3", "white", "#D6604D"))
  ht <- Heatmap(mat, name = "Synergy\n(ΔNES)", col = col_fun,
                cluster_rows = TRUE, cluster_columns = TRUE,
                row_names_gp    = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 11, fontface = "bold"),
                column_title    = paste0(gs_sel, " – Drug Combination Synergy Score (EC)\n",
                                         "ΔNES = Combo NES − max(Single Drug NES)"),
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  v <- mat[i, j]
                  if (!is.na(v) && abs(v) > 0.5)
                    grid.text(sprintf("%.2f", v), x, y, gp = gpar(fontsize = 8))
                }
  )
  pdf(paste0(dir.base, "Figures/Synergy/synergy_heatmap_", gs_sel, "_EC.pdf"),
      width = 12, height = 14)
  draw(ht); dev.off()
  cat("  Saved synergy heatmap:", gs_sel, "\n")
}

plot_synergy_heatmap(synergy_h,  "Hallmark")
plot_synergy_heatmap(synergy_bp, "GOBP", n_paths = 40)

# ── 5D. Synergy Bar Charts ────────────────────────────────────────────────────
plot_synergy_bars <- function(synergy_df, gs_sel = "Hallmark", n_top = 10) {
  dat <- synergy_df %>%
    filter(gs_type == gs_sel) %>%
    mutate(pathway_clean = clean_path(pathway))
  plots <- lapply(unique(dat$combo_label), function(combo) {
    sub <- dat %>%
      filter(combo_label == combo) %>%
      slice_max(abs(synergy_score), n = n_top * 2) %>%
      arrange(synergy_score) %>%
      mutate(pathway_clean = factor(pathway_clean, levels = pathway_clean),
             direction = if_else(synergy_score > 0, "Synergistic", "Antagonistic"))
    ggplot(sub, aes(x = synergy_score, y = pathway_clean, fill = direction)) +
      geom_col() +
      scale_fill_manual(values = c(Synergistic = "#D6604D", Antagonistic = "#4393C3")) +
      geom_vline(xintercept = 0, color = "black") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 10)) +
      labs(x = "Synergy Score (ΔNES)", y = NULL,
           title = paste0(combo, " vs. best single drug"))
  })
  combined <- wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title    = paste0(gs_sel, ": Drug Combination Synergy in EC"),
      subtitle = "Red = combo stronger than best single drug; Blue = antagonism",
      theme    = theme(plot.title = element_text(face = "bold", size = 13))
    )
  ggsave(paste0(dir.base, "Figures/Synergy/synergy_bars_", gs_sel, "_EC.pdf"),
         combined, width = 18, height = 14)
  ggsave(paste0(dir.base, "Figures/Synergy/synergy_bars_", gs_sel, "_EC.png"),
         combined, width = 18, height = 14, dpi = 300)
  cat("  Saved synergy bar chart:", gs_sel, "\n")
}

plot_synergy_bars(synergy_h,  "Hallmark")
plot_synergy_bars(synergy_bp, "GOBP")

# ── 5E. Radar Chart by Biological Theme ──────────────────────────────────────
plot_radar_by_theme <- function(gsea_df, fdr_thresh = 0.25) {
  hallmark_themes <- list(
    Inflammation  = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                      "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
    Metabolism    = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM",
                      "HALLMARK_GLYCOLYSIS", "HALLMARK_ADIPOGENESIS"),
    Angiogenesis  = c("HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_VEGFA"),
    Fibrosis      = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING",
                      "HALLMARK_NOTCH_SIGNALING"),
    Proliferation = c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"),
    Apoptosis     = c("HALLMARK_APOPTOSIS", "HALLMARK_P53_PATHWAY",
                      "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
  )
  theme_nes <- lapply(names(hallmark_themes), function(theme) {
    gsea_df %>%
      filter(gs_type == "Hallmark", pathway %in% hallmark_themes[[theme]]) %>%
      group_by(drug, condition_type) %>%
      summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
      mutate(theme = theme)
  }) %>% bind_rows()
  
  p <- ggplot(theme_nes, aes(x = theme, y = mean_NES, group = drug,
                             color = drug, fill = drug)) +
    geom_polygon(alpha = 0.1, linewidth = 0.8) +
    geom_point(size = 2) +
    coord_polar() +
    scale_color_manual(values = drug_colors) +
    scale_fill_manual(values  = drug_colors) +
    facet_wrap(~ condition_type, labeller = labeller(
      condition_type = c(single = "Single Drugs", combo = "Combinations"))) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(size = 9),
          legend.title = element_blank(),
          plot.title   = element_text(face = "bold", hjust = 0.5)) +
    labs(title = "Biological Theme Profiles – Endothelial Cells", x = NULL, y = "Mean NES")
  ggsave(paste0(dir.base, "Figures/NES_Matrices/radar_themes_EC.pdf"), p, width = 14, height = 7)
  ggsave(paste0(dir.base, "Figures/NES_Matrices/radar_themes_EC.png"), p, width = 14, height = 7, dpi = 300)
  cat("  Saved radar chart\n")
}

plot_radar_by_theme(gsea_df)

# ── 5F. UpSet Plot: shared significant pathways ───────────────────────────────
plot_path_overlap <- function(gsea_df, gs_sel = "Hallmark", fdr_thresh = 0.25) {
  sig_paths <- gsea_df %>%
    filter(gs_type == gs_sel, padj < fdr_thresh) %>%
    dplyr::select(drug, pathway)
  drug_list <- split(sig_paths$pathway, sig_paths$drug)
  mat_up    <- make_comb_mat(drug_list)
  pdf(paste0(dir.base, "Figures/NES_Matrices/upset_", gs_sel, "_EC.pdf"), width = 14, height = 7)
  draw(UpSet(mat_up,
             set_order  = c(drug_names, names(drug_list)[!(names(drug_list) %in% drug_names)]),
             comb_order = order(comb_size(mat_up), decreasing = TRUE),
             top_annotation  = upset_top_annotation(mat_up,  add_numbers = TRUE),
             left_annotation = upset_left_annotation(mat_up, add_numbers = TRUE)
  ))
  dev.off()
  cat("  Saved UpSet plot:", gs_sel, "\n")
}

plot_path_overlap(gsea_df, "Hallmark")

# ── 5G. Summary Table ─────────────────────────────────────────────────────────
summary_tbl <- gsea_df %>%
  filter(padj < 0.25) %>%
  group_by(drug, condition_type, gs_type) %>%
  summarise(
    n_sig     = n(),
    n_pos_NES = sum(NES > 0),
    n_neg_NES = sum(NES < 0),
    top_up    = paste(head(pathway[NES > 0][order(-NES[NES > 0])], 3), collapse = "; "),
    top_down  = paste(head(pathway[NES < 0][order(NES[NES < 0])],  3), collapse = "; "),
    .groups   = "drop"
  )

write.csv(summary_tbl,
          paste0(dir.base, "GSEA/summary_significant_pathways_EC.csv"), row.names = FALSE)

cat("\n════════════════════════════════════════════════\n")
cat("  Analysis complete! All outputs in:\n  ", dir.base, "\n")
cat("════════════════════════════════════════════════\n")
print(summary_tbl %>% dplyr::select(drug, condition_type, gs_type, n_sig))



