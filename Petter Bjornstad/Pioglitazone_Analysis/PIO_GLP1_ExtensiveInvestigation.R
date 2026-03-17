# ============================================================================
# GLP-1 x PIOGLITAZONE EXTENDED ANALYSIS
# 1. Four-group demographics (PIO-/GLP1-, PIO+/GLP1-, PIO-/GLP1+, PIO+/GLP1+)
# 2. GLP-1RA main effect analysis (mirroring the pioglitazone sensitivity script)
# 3. PT subtype deep-dive
# 4. Cross-comparison with interaction results from sensitivity script
# ============================================================================

library(scran);      library(tidyverse);   library(patchwork)
library(cowplot);    library(ggpubr);      library(rstatix)
library(data.table); library(pheatmap);    library(readxl)
library(SingleCellExperiment); library(scater); library(Seurat)
library(dplyr);      library(nebula);      library(Matrix)
library(ggplot2);    library(gt);          library(gtsummary)
library(fgsea);      library(msigdbr);     library(stringr)
library(gridExtra);  library(grid);        library(ggrepel)
library(RColorBrewer)

# ============================================================================
# DIRECTORY SETUP
# ============================================================================

# *** Point these at your existing objects & outputs ***
dir.sensitivity <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/Sensitivity/'
dir.interaction <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/Sensitivity/Interaction_PIO_GLP1/'
dir.root        <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/GLP1_Extended/'

for (d in c('Demographics/', 'GLP1_Main/NEBULA/', 'GLP1_Main/GSEA/Raw/',
            'GLP1_Main/GSEA/DotPlots/', 'GLP1_Main/GSEA/Heatmaps/',
            'GLP1_Main/VolcanoPlots/', 'PT_Subtypes/NEBULA/',
            'PT_Subtypes/GSEA/Raw/', 'PT_Subtypes/VolcanoPlots/',
            'Comparison/')) {
  dir.create(paste0(dir.root, d), recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# SHARED CONSTANTS
# ============================================================================

celltypes_vec <- c('All', 'PT', 'TAL', 'EC', 'POD', 'DCT', 'IC')

# PT subtypes from KPMP_celltype
pt_subtypes_vec <- c('aPT', 'PT-S1/S2', 'PT-S3')

GSEA_PARAMS <- list(minSize = 15, maxSize = 500, nPermSimple = 10000)
P_THRESH    <- 0.05

hallmark_list <- split(msigdbr(species = "Homo sapiens", category = "H")$gene_symbol,
                       msigdbr(species = "Homo sapiens", category = "H")$gs_name)
go_bp_list    <- split(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")$gene_symbol,
                       msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")$gs_name)
geneset_types <- list(Hallmark = hallmark_list, GO_BP = go_bp_list)

# ============================================================================
# ASSUMES scrna_small IS ALREADY LOADED WITH COLUMNS:
#   group_labels  (PIO: Yes/No)
#   epic_glp1ra_1 (GLP1RA: Yes/No)
#   pooled_offset, celltype1, celltype2, DCT_celltype, celltype_rpca, KPMP_celltype
# If starting fresh, load your .RData and re-run the setup from the
# sensitivity script through the "POOLED OFFSET" section first.
# ============================================================================

# Assign four-group label
scrna_small$four_group <- case_when(
  scrna_small$group_labels  == "Yes" & scrna_small$epic_glp1ra_1 == "Yes" ~ "PIO+GLP1+",
  scrna_small$group_labels  == "Yes" & scrna_small$epic_glp1ra_1 == "No"  ~ "PIO+GLP1-",
  scrna_small$group_labels  == "No"  & scrna_small$epic_glp1ra_1 == "Yes" ~ "PIO-GLP1+",
  scrna_small$group_labels  == "No"  & scrna_small$epic_glp1ra_1 == "No"  ~ "PIO-GLP1-",
  TRUE ~ NA_character_
)
scrna_small$four_group <- factor(scrna_small$four_group,
                                 levels = c("PIO-GLP1-", "PIO+GLP1-", "PIO-GLP1+", "PIO+GLP1+"))

# ============================================================================
# SECTION 1 — FOUR-GROUP DEMOGRAPHICS
# ============================================================================

cat("\n=== SECTION 1: Four-Group Demographics ===\n")

patient_meta <- scrna_small@meta.data %>%
  filter(!duplicated(record_id)) %>%
  filter(!is.na(four_group))

# ---- 1a. Patient-level N per group ------------------------------------------

group_counts <- patient_meta %>%
  count(four_group, name = "n_patients") %>%
  mutate(pct = round(100 * n_patients / sum(n_patients), 1))

cat("\nPatient counts per four-group:\n"); print(group_counts)
write.csv(group_counts, paste0(dir.root, 'Demographics/PatientCounts_FourGroup.csv'), row.names = FALSE)

# ---- 1b. Cell counts per group x cell type ----------------------------------

cell_counts <- scrna_small@meta.data %>%
  filter(!is.na(four_group)) %>%
  count(four_group, celltype2, name = "n_cells")

write.csv(cell_counts, paste0(dir.root, 'Demographics/CellCounts_FourGroup_byCelltype.csv'), row.names = FALSE)

p_cells <- ggplot(cell_counts, aes(x = celltype2, y = n_cells, fill = four_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("PIO-GLP1-" = "#999999",
                               "PIO+GLP1-" = "#E69F00",
                               "PIO-GLP1+" = "#56B4E9",
                               "PIO+GLP1+" = "#009E73")) +
  theme_bw() +
  labs(x = "Cell Type", y = "Number of Cells", fill = "Group",
       title = "Cell counts by four-group (PIO x GLP-1RA)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title  = element_text(face = "bold"))

ggsave(paste0(dir.root, 'Demographics/CellCounts_FourGroup_barplot.pdf'), p_cells, width = 10, height = 5)
ggsave(paste0(dir.root, 'Demographics/CellCounts_FourGroup_barplot.png'), p_cells, width = 10, height = 5, dpi = 200)

# ---- 1c. Continuous covariate summary table (BMI, age if available) ---------

cont_vars <- intersect(c("bmi", "age"), colnames(patient_meta))

if (length(cont_vars) > 0) {
  demo_summary <- patient_meta %>%
    dplyr::select(four_group, all_of(cont_vars)) %>%
    pivot_longer(-four_group) %>%
    group_by(four_group, name) %>%
    summarise(
      mean = round(mean(value, na.rm = TRUE), 2),
      sd   = round(sd(value,   na.rm = TRUE), 2),
      n    = sum(!is.na(value)),
      .groups = "drop"
    )
  write.csv(demo_summary, paste0(dir.root, 'Demographics/ContinuousVars_FourGroup.csv'), row.names = FALSE)
  cat("\nContinuous variable summary:\n"); print(demo_summary)
}

# ---- 1d. Sex distribution per group -----------------------------------------

if ("sex" %in% colnames(patient_meta)) {
  sex_table <- patient_meta %>%
    filter(!is.na(sex)) %>%
    count(four_group, sex) %>%
    group_by(four_group) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()
  write.csv(sex_table, paste0(dir.root, 'Demographics/Sex_FourGroup.csv'), row.names = FALSE)
  cat("\nSex distribution:\n"); print(sex_table)
}

# ---- 1e. Metformin co-use per group -----------------------------------------

if ("epic_mfm_1" %in% colnames(patient_meta)) {
  mfm_table <- patient_meta %>%
    filter(!is.na(epic_mfm_1)) %>%
    count(four_group, epic_mfm_1) %>%
    group_by(four_group) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()
  write.csv(mfm_table, paste0(dir.root, 'Demographics/Metformin_FourGroup.csv'), row.names = FALSE)
  cat("\nMetformin co-use:\n"); print(mfm_table)
}

# ---- 1f. UMAP colored by four-group -----------------------------------------

if ("UMAP_1" %in% colnames(scrna_small@meta.data)) {
  umap_df <- scrna_small@meta.data %>%
    filter(!is.na(four_group)) %>%
    dplyr::select(UMAP_1, UMAP_2, four_group, celltype2)
  
  p_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = four_group)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = c("PIO-GLP1-" = "#999999",
                                  "PIO+GLP1-" = "#E69F00",
                                  "PIO-GLP1+" = "#56B4E9",
                                  "PIO+GLP1+" = "#009E73"),
                       guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_classic() +
    labs(color = "Group", title = "UMAP — PIO x GLP-1RA Groups") +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(paste0(dir.root, 'Demographics/UMAP_FourGroup.pdf'), p_umap, width = 8, height = 6)
  ggsave(paste0(dir.root, 'Demographics/UMAP_FourGroup.png'), p_umap, width = 8, height = 6, dpi = 200)
}

cat("Section 1 complete.\n")

# ============================================================================
# SECTION 2 — GLP-1RA MAIN EFFECT NEBULA
# (Mirrors the PIO sensitivity analysis but with GLP-1RA as exposure)
# ============================================================================

cat("\n=== SECTION 2: GLP-1RA Main Effect (NEBULA) ===\n")

# GLP-1RA covariate models to test
GLP1_COVARIATES <- list(
  base      = NULL,
  plus_pio  = "group_labels",   # adjust for PIO use
  plus_bmi  = "bmi",
  plus_sex  = "sex",
  plus_mfm  = "epic_mfm_1"
)

run_nebula_glp1_main <- function(so_obj, dir.out, celltype,
                                 covariate_name = NULL, model_label = "base") {
  # Cell type subsetting
  if (celltype == 'All') {
    so_c <- so_obj
  } else if (celltype %in% c('TAL', 'EC', 'PT')) {
    so_c <- subset(so_obj, celltype2 == celltype)
  } else if (celltype == 'IC') {
    so_c <- subset(so_obj, celltype1 %in% c("B", "T", "MON", "MAC"))
  } else if (celltype == 'DCT') {
    so_c <- subset(so_obj, DCT_celltype == "DCT")
  } else if (celltype == 'POD') {
    so_c <- subset(so_obj, celltype1 == "POD")
  } else {
    if      (celltype %in% unique(so_obj$celltype2)) so_c <- subset(so_obj, celltype2 == celltype)
    else if (celltype %in% unique(so_obj$celltype1)) so_c <- subset(so_obj, celltype1 == celltype)
    else { cat(paste0("Cell type '", celltype, "' not found. Skipping.\n")); return(invisible(NULL)) }
  }
  
  DefaultAssay(so_c) <- "RNA"
  ct_clean <- str_replace_all(str_replace_all(celltype, "/", "_"), "-", "_")
  cat(paste0("\n=== GLP1 main | ", ct_clean, " | model: ", model_label, " ===\n"))
  
  counts_mat <- round(GetAssayData(so_c, layer = "counts"))
  meta_g     <- so_c@meta.data
  
  required_cols <- c("epic_glp1ra_1", if (!is.null(covariate_name)) covariate_name)
  missing_cols  <- setdiff(required_cols, colnames(meta_g))
  if (length(missing_cols) > 0) {
    cat("Missing columns:", paste(missing_cols, collapse = ", "), "— skipping.\n")
    return(invisible(NULL))
  }
  
  complete_idx <- complete.cases(meta_g[, required_cols, drop = FALSE])
  meta_g    <- meta_g[complete_idx, ]
  counts_mat <- counts_mat[, complete_idx]
  
  tmp_df  <- meta_g %>% distinct(record_id, .keep_all = TRUE)
  num_yes <- sum(tmp_df$epic_glp1ra_1 == "Yes", na.rm = TRUE)
  num_no  <- sum(tmp_df$epic_glp1ra_1 == "No",  na.rm = TRUE)
  cat("Participants — GLP1+:", num_yes, "GLP1-:", num_no, "\n")
  
  formula_str <- if (is.null(covariate_name)) "~epic_glp1ra_1" else
    paste0("~epic_glp1ra_1 + ", covariate_name)
  pred_gene <- model.matrix(as.formula(formula_str), data = meta_g)
  
  lib    <- meta_g$pooled_offset
  data_g <- group_cell(count = counts_mat, id = meta_g$record_id,
                       pred = pred_gene, offset = lib)
  if (is.null(data_g))
    data_g <- list(count = counts_mat, id = meta_g$record_id,
                   pred = pred_gene, offset = lib)
  
  result <- nebula(count = data_g$count, id = data_g$id,
                   pred = data_g$pred, ncore = 1, reml = TRUE,
                   model = "NBLMM", output_re = TRUE, covariance = TRUE,
                   offset = data_g$offset)
  
  out_df <- as.data.frame(result)
  out_df$num_cells    <- nrow(meta_g)
  out_df$num_glp1_yes <- num_yes
  out_df$num_glp1_no  <- num_no
  out_df$model_label  <- model_label
  out_df$covariate    <- ifelse(is.null(covariate_name), "none", covariate_name)
  
  out_file <- paste0(dir.out, 'GLP1_Main/NEBULA/GLP1_NEBULA_', ct_clean, '_', model_label, '.csv')
  write.csv(out_df, out_file, row.names = FALSE)
  cat("Saved:", out_file, "\n")
  return(invisible(out_df))
}

# Run GLP-1 main effect across all cell types x models
for (ct in celltypes_vec) {
  for (m in names(GLP1_COVARIATES)) {
    run_nebula_glp1_main(so_obj = scrna_small, dir.out = dir.root,
                         celltype = ct, covariate_name = GLP1_COVARIATES[[m]],
                         model_label = m)
  }
}

# Helper: load GLP-1 effect
load_glp1_effect <- function(dir.root, ct_clean, model_label) {
  fp <- paste0(dir.root, 'GLP1_Main/NEBULA/GLP1_NEBULA_', ct_clean, '_', model_label, '.csv')
  if (!file.exists(fp)) return(NULL)
  df <- read.csv(fp)
  # GLP1RA Yes coefficient column
  glp1_logfc <- grep("logFC.*glp1ra.*Yes|logFC.*epic_glp1ra_1Yes", names(df), value = TRUE, ignore.case = TRUE)
  glp1_p     <- grep("^summary\\.p_.*glp1ra.*Yes|^summary\\.p_.*epic_glp1ra_1Yes", names(df), value = TRUE, ignore.case = TRUE)
  if (length(glp1_logfc) == 0) {
    # fallback: second logFC column (first is intercept)
    glp1_logfc <- grep("logFC", names(df), value = TRUE)[2]
    glp1_p     <- grep("^summary\\.p_", names(df), value = TRUE)[2]
  }
  df %>%
    transmute(
      gene   = summary.gene,
      logFC  = .data[[glp1_logfc[1]]],
      pvalue = .data[[glp1_p[1]]],
      fdr    = p.adjust(pvalue, method = "BH")
    )
}

# ---- 2a. GLP-1 Volcano Plots ------------------------------------------------

cat("\n--- GLP-1 volcano plots ---\n")

for (ct in celltypes_vec) {
  ct_clean  <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  vol_list  <- list()
  
  for (m in names(GLP1_COVARIATES)) {
    df <- load_glp1_effect(dir.root, ct_clean, m)
    if (is.null(df)) next
    
    df <- df %>%
      mutate(diffexp = case_when(
        pvalue < P_THRESH & logFC > 0 ~ "Up",
        pvalue < P_THRESH & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )) %>%
      arrange(pvalue) %>%
      mutate(label = ifelse(row_number() <= 8, gene, NA)) %>%
      filter(abs(logFC) < 10)
    
    title_str <- if (m == "base") "GLP-1 base" else paste0("GLP-1 +", m)
    
    vol_list[[m]] <- ggplot(df, aes(x = logFC, y = -log10(pvalue), color = diffexp, label = label)) +
      geom_point(size = 1.2, alpha = 0.7) +
      ggrepel::geom_text_repel(size = 2.2, color = "black", max.overlaps = 12) +
      scale_color_manual(values = c(Down = "steelblue", NS = "grey70", Up = "firebrick"),
                         labels = c(Down = "Down in GLP-1+", NS = "NS", Up = "Up in GLP-1+")) +
      geom_hline(yintercept = -log10(P_THRESH), color = "blue", linetype = "dashed", linewidth = 0.4) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
      theme_classic() +
      labs(x = "LogFC", y = "-log10(p)", color = NULL,
           title = paste0(ct, ": ", title_str)) +
      theme(plot.title = element_text(face = "bold", size = 10),
            legend.position = "bottom", aspect.ratio = 1)
  }
  
  if (length(vol_list) > 1) {
    combined_vol <- wrap_plots(vol_list, ncol = 3) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
    ggsave(paste0(dir.root, 'GLP1_Main/VolcanoPlots/GLP1_Volcano_', ct_clean, '.pdf'),
           combined_vol, width = 18, height = 6)
    ggsave(paste0(dir.root, 'GLP1_Main/VolcanoPlots/GLP1_Volcano_', ct_clean, '.png'),
           combined_vol, width = 18, height = 6, dpi = 200)
    cat("Saved GLP-1 volcano for", ct, "\n")
  }
}

# ---- 2b. GLP-1 GSEA ---------------------------------------------------------

cat("\n--- GLP-1 GSEA ---\n")

all_glp1_gsea <- list()

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  for (m in names(GLP1_COVARIATES)) {
    df <- load_glp1_effect(dir.root, ct_clean, m)
    if (is.null(df)) next
    
    ranked <- df %>%
      filter(!is.na(logFC)) %>%
      arrange(desc(logFC)) %>%
      { setNames(.$logFC, .$gene) }
    
    for (gs_name in names(geneset_types)) {
      key <- paste(ct, m, gs_name, sep = "__")
      cat("  GLP1 GSEA:", key, "\n")
      set.seed(42)
      res <- fgsea(pathways = geneset_types[[gs_name]], stats = ranked,
                   minSize = GSEA_PARAMS$minSize, maxSize = GSEA_PARAMS$maxSize,
                   nPermSimple = GSEA_PARAMS$nPermSimple)
      res$celltype      <- ct
      res$model         <- m
      res$geneset_type  <- gs_name
      res$pathway_clean <- gsub("HALLMARK_|GOBP_", "", res$pathway) %>%
        gsub("_", " ", .) %>% tools::toTitleCase()
      all_glp1_gsea[[key]] <- res
    }
  }
}

# Save raw GSEA
combined_glp1_gsea <- bind_rows(all_glp1_gsea) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(combined_glp1_gsea, paste0(dir.root, 'GLP1_Main/GSEA/Raw/GLP1_all_gsea.csv'), row.names = FALSE)

# Dot plots (base model, top pathways)
for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  for (gs_name in names(geneset_types)) {
    key     <- paste(ct, "base", gs_name, sep = "__")
    base_df <- all_glp1_gsea[[key]]
    if (is.null(base_df)) next
    
    top_paths <- base_df %>%
      filter(pval < P_THRESH) %>%
      arrange(pval) %>%
      slice_head(n = 20) %>%
      pull(pathway)
    if (length(top_paths) == 0) top_paths <- base_df %>% arrange(pval) %>% slice_head(n = 20) %>% pull(pathway)
    if (length(top_paths) == 0) next
    
    plot_df <- lapply(names(GLP1_COVARIATES), function(m) {
      df <- all_glp1_gsea[[paste(ct, m, gs_name, sep = "__")]]
      if (is.null(df)) return(NULL)
      df %>% filter(pathway %in% top_paths) %>%
        dplyr::select(pathway, pathway_clean, NES, pval) %>%
        mutate(model = m)
    }) %>% bind_rows() %>%
      mutate(model = factor(model, levels = names(GLP1_COVARIATES)),
             pathway_clean = factor(pathway_clean,
                                    levels = base_df %>% filter(pathway %in% top_paths) %>%
                                      arrange(NES) %>% pull(pathway_clean)))
    
    p <- ggplot(plot_df, aes(x = model, y = pathway_clean, color = NES, size = -log10(pval))) +
      geom_point() +
      geom_point(data = plot_df %>% filter(pval < P_THRESH),
                 shape = 8, color = "black", size = 1.5) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_size_continuous(range = c(1, 8)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title  = element_text(face = "bold", size = 11)) +
      labs(x = "Model", y = NULL,
           title = paste0(ct, " — ", gs_name, "\nGLP-1 effect: top pathways across covariate models"),
           caption = "★ = p < 0.05")
    
    h <- max(5, length(top_paths) * 0.35 + 2)
    ggsave(paste0(dir.root, 'GLP1_Main/GSEA/DotPlots/GLP1_DotPlot_', ct_clean, '_', gs_name, '.pdf'),
           p, width = 9, height = h)
    ggsave(paste0(dir.root, 'GLP1_Main/GSEA/DotPlots/GLP1_DotPlot_', ct_clean, '_', gs_name, '.png'),
           p, width = 9, height = h, dpi = 200)
    cat("  Saved GLP-1 dot plot:", ct, gs_name, "\n")
  }
}

cat("Section 2 complete.\n")

# ============================================================================
# SECTION 3 — PT SUBTYPE DEEP-DIVE
# ============================================================================

cat("\n=== SECTION 3: PT Subtype Analysis ===\n")

# Assign PT subtype label from KPMP_celltype
scrna_small$pt_subtype <- case_when(
  scrna_small$KPMP_celltype == "aPT"     ~ "aPT",
  scrna_small$KPMP_celltype == "PT-S1/S2" ~ "PT-S1/S2",
  scrna_small$KPMP_celltype == "PT-S3"   ~ "PT-S3",
  TRUE ~ NA_character_
)

# Observed PT subtypes in this dataset
observed_pt <- scrna_small@meta.data %>%
  filter(!is.na(pt_subtype)) %>%
  count(pt_subtype) %>%
  filter(n > 50)   # require at least 50 cells
cat("PT subtypes with >50 cells:\n"); print(observed_pt)

run_nebula_pt_subtype <- function(so_obj, dir.out, pt_sub, exposure = "group_labels",
                                  covariate_name = NULL, model_label = "pio_base") {
  so_c <- subset(so_obj, pt_subtype == pt_sub)
  n_cells <- ncol(so_c)
  if (n_cells < 50) {
    cat(paste0("  Skipping ", pt_sub, " — only ", n_cells, " cells.\n"))
    return(invisible(NULL))
  }
  
  DefaultAssay(so_c) <- "RNA"
  sub_clean <- str_replace_all(str_replace_all(pt_sub, "/", "_"), "-", "_")
  cat(paste0("\n=== PT subtype: ", sub_clean, " | exposure: ", exposure,
             " | model: ", model_label, " ===\n"))
  
  counts_mat <- round(GetAssayData(so_c, layer = "counts"))
  meta_g     <- so_c@meta.data
  
  required_cols <- c(exposure, if (!is.null(covariate_name)) covariate_name)
  complete_idx  <- complete.cases(meta_g[, required_cols, drop = FALSE])
  meta_g     <- meta_g[complete_idx, ]
  counts_mat  <- counts_mat[, complete_idx]
  
  cat("Cells with complete data:", nrow(meta_g), "\n")
  tmp_df <- meta_g %>% distinct(record_id, .keep_all = TRUE)
  cat("Participants:", nrow(tmp_df), "| Exposed Yes:",
      sum(tmp_df[[exposure]] == "Yes", na.rm = TRUE), "\n")
  
  formula_str <- if (is.null(covariate_name)) paste0("~", exposure) else
    paste0("~", exposure, " + ", covariate_name)
  pred_gene <- model.matrix(as.formula(formula_str), data = meta_g)
  
  lib    <- meta_g$pooled_offset
  data_g <- group_cell(count = counts_mat, id = meta_g$record_id,
                       pred = pred_gene, offset = lib)
  if (is.null(data_g))
    data_g <- list(count = counts_mat, id = meta_g$record_id,
                   pred = pred_gene, offset = lib)
  
  result <- nebula(count = data_g$count, id = data_g$id,
                   pred = data_g$pred, ncore = 1, reml = TRUE,
                   model = "NBLMM", output_re = TRUE, covariance = TRUE,
                   offset = data_g$offset)
  
  out_df <- as.data.frame(result)
  out_df$pt_subtype   <- pt_sub
  out_df$exposure     <- exposure
  out_df$model_label  <- model_label
  out_df$num_cells    <- nrow(meta_g)
  
  out_file <- paste0(dir.out, 'PT_Subtypes/NEBULA/PT_NEBULA_', sub_clean, '_', model_label, '.csv')
  write.csv(out_df, out_file, row.names = FALSE)
  cat("Saved:", out_file, "\n")
  return(invisible(out_df))
}

# Run for each PT subtype: PIO main effect, GLP-1 main effect, PIO+GLP-1 adjusted
pt_models <- list(
  pio_base    = list(exposure = "group_labels",  covariate = NULL),
  pio_glp1adj = list(exposure = "group_labels",  covariate = "epic_glp1ra_1"),
  glp1_base   = list(exposure = "epic_glp1ra_1", covariate = NULL),
  glp1_pioadj = list(exposure = "epic_glp1ra_1", covariate = "group_labels")
)

for (pt_sub in observed_pt$pt_subtype) {
  for (m in names(pt_models)) {
    run_nebula_pt_subtype(
      so_obj       = scrna_small,
      dir.out      = dir.root,
      pt_sub       = pt_sub,
      exposure     = pt_models[[m]]$exposure,
      covariate_name = pt_models[[m]]$covariate,
      model_label  = m
    )
  }
}

# Helper: load PT subtype effect
load_pt_effect <- function(dir.root, sub_clean, model_label, exposure) {
  fp <- paste0(dir.root, 'PT_Subtypes/NEBULA/PT_NEBULA_', sub_clean, '_', model_label, '.csv')
  if (!file.exists(fp)) return(NULL)
  df <- read.csv(fp)
  # Find the "Yes" coefficient for the given exposure
  logfc_col <- grep(paste0("logFC.*", exposure, ".*Yes|logFC.*", gsub("_", "", exposure), ".*Yes"),
                    names(df), value = TRUE, ignore.case = TRUE)[1]
  p_col     <- grep(paste0("^summary\\.p_.*", exposure, ".*Yes|^summary\\.p_.*", gsub("_", "", exposure), ".*Yes"),
                    names(df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(logfc_col)) {
    logfc_col <- grep("logFC", names(df), value = TRUE)[2]
    p_col     <- grep("^summary\\.p_", names(df), value = TRUE)[2]
  }
  df %>%
    transmute(
      gene   = summary.gene,
      logFC  = .data[[logfc_col]],
      pvalue = .data[[p_col]],
      fdr    = p.adjust(pvalue, method = "BH")
    )
}

# ---- 3a. PT subtype volcano plots (PIO base) --------------------------------

cat("\n--- PT subtype volcano plots ---\n")

for (pt_sub in observed_pt$pt_subtype) {
  sub_clean <- str_replace_all(str_replace_all(pt_sub, "/", "_"), "-", "_")
  vol_list  <- list()
  
  for (m in names(pt_models)) {
    exp_name <- pt_models[[m]]$exposure
    df <- load_pt_effect(dir.root, sub_clean, m, exp_name)
    if (is.null(df)) next
    
    df <- df %>%
      mutate(diffexp = case_when(
        pvalue < P_THRESH & logFC > 0 ~ "Up",
        pvalue < P_THRESH & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )) %>%
      arrange(pvalue) %>%
      mutate(label = ifelse(row_number() <= 8, gene, NA)) %>%
      filter(abs(logFC) < 10)
    
    exposure_label <- if (grepl("group_labels", exp_name)) "PIO" else "GLP-1"
    color_up   <- if (grepl("group_labels", exp_name)) "purple" else "firebrick"
    color_down <- if (grepl("group_labels", exp_name)) "orange" else "steelblue"
    
    vol_list[[m]] <- ggplot(df, aes(x = logFC, y = -log10(pvalue), color = diffexp, label = label)) +
      geom_point(size = 1.2, alpha = 0.7) +
      ggrepel::geom_text_repel(size = 2.2, color = "black", max.overlaps = 12) +
      scale_color_manual(values = c(Down = color_down, NS = "grey70", Up = color_up)) +
      geom_hline(yintercept = -log10(P_THRESH), color = "blue", linetype = "dashed", linewidth = 0.4) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
      theme_classic() +
      labs(x = "LogFC", y = "-log10(p)", color = NULL,
           title = paste0(pt_sub, ": ", exposure_label, " (", m, ")")) +
      theme(plot.title = element_text(face = "bold", size = 9),
            legend.position = "bottom", aspect.ratio = 1)
  }
  
  if (length(vol_list) > 1) {
    combined_vol <- wrap_plots(vol_list, ncol = 2) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
    ggsave(paste0(dir.root, 'PT_Subtypes/VolcanoPlots/PTsub_Volcano_', sub_clean, '.pdf'),
           combined_vol, width = 14, height = 10)
    ggsave(paste0(dir.root, 'PT_Subtypes/VolcanoPlots/PTsub_Volcano_', sub_clean, '.png'),
           combined_vol, width = 14, height = 10, dpi = 200)
    cat("Saved PT subtype volcano:", pt_sub, "\n")
  }
}

# ---- 3b. PT subtype GSEA (PIO base) -----------------------------------------

cat("\n--- PT subtype GSEA ---\n")

all_pt_gsea <- list()

for (pt_sub in observed_pt$pt_subtype) {
  sub_clean <- str_replace_all(str_replace_all(pt_sub, "/", "_"), "-", "_")
  
  for (m in names(pt_models)) {
    exp_name <- pt_models[[m]]$exposure
    df <- load_pt_effect(dir.root, sub_clean, m, exp_name)
    if (is.null(df)) next
    
    ranked <- df %>% filter(!is.na(logFC)) %>%
      arrange(desc(logFC)) %>% { setNames(.$logFC, .$gene) }
    
    for (gs_name in names(geneset_types)) {
      key <- paste(pt_sub, m, gs_name, sep = "__")
      cat("  PT GSEA:", key, "\n")
      set.seed(42)
      res <- fgsea(pathways = geneset_types[[gs_name]], stats = ranked,
                   minSize = GSEA_PARAMS$minSize, maxSize = GSEA_PARAMS$maxSize,
                   nPermSimple = GSEA_PARAMS$nPermSimple)
      res$pt_subtype    <- pt_sub
      res$model         <- m
      res$exposure      <- exp_name
      res$geneset_type  <- gs_name
      res$pathway_clean <- gsub("HALLMARK_|GOBP_", "", res$pathway) %>%
        gsub("_", " ", .) %>% tools::toTitleCase()
      all_pt_gsea[[key]] <- res
    }
  }
}

combined_pt_gsea <- bind_rows(all_pt_gsea) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(combined_pt_gsea, paste0(dir.root, 'PT_Subtypes/GSEA/Raw/PT_all_gsea.csv'), row.names = FALSE)

# ---- 3c. PT subtype heatmap: PIO effect NES across subtypes -----------------

cat("\n--- PT subtype NES comparison heatmap ---\n")

for (gs_name in names(geneset_types)) {
  # Gather base PIO NES for all PT subtypes
  pt_nes_list <- lapply(observed_pt$pt_subtype, function(pt_sub) {
    key <- paste(pt_sub, "pio_base", gs_name, sep = "__")
    df  <- all_pt_gsea[[key]]
    if (is.null(df)) return(NULL)
    df %>% dplyr::select(pathway, pathway_clean, NES, pval) %>%
      rename_with(~ paste0(c("NES", "pval"), "_", str_replace_all(str_replace_all(pt_sub, "/", "_"), "-", "_")),
                  c("NES", "pval"))
  })
  pt_nes_list <- Filter(Negate(is.null), pt_nes_list)
  if (length(pt_nes_list) < 2) next
  
  wide <- Reduce(function(a, b) full_join(a, b, by = c("pathway", "pathway_clean")), pt_nes_list)
  
  # Filter to pathways sig in any PT subtype (base PIO model)
  pval_cols <- grep("^pval_", names(wide), value = TRUE)
  nes_cols  <- grep("^NES_",  names(wide), value = TRUE)
  
  sig_any  <- rowSums(wide[, pval_cols] < P_THRESH, na.rm = TRUE) > 0
  wide_sig <- wide[sig_any, ]
  if (nrow(wide_sig) == 0) next
  if (nrow(wide_sig) > 50) wide_sig <- wide_sig[order(rowMeans(wide_sig[, pval_cols], na.rm = TRUE))[1:50], ]
  
  mat <- wide_sig %>% column_to_rownames("pathway_clean") %>%
    dplyr::select(all_of(nes_cols)) %>% as.matrix()
  colnames(mat) <- str_remove(colnames(mat), "^NES_")
  
  pval_mat <- wide_sig %>% column_to_rownames("pathway_clean") %>%
    dplyr::select(all_of(pval_cols)) %>% as.matrix()
  
  cell_labels <- matrix(
    ifelse(pval_mat < P_THRESH, paste0(round(mat, 2), "*"), round(mat, 2)),
    nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)
  )
  col_lim <- max(abs(mat), na.rm = TRUE)
  
  png(paste0(dir.root, 'PT_Subtypes/GSEA/PT_NES_Heatmap_PIOeffect_', gs_name, '.png'),
      width = 900, height = max(400, nrow(mat) * 22 + 200), res = 120)
  pheatmap::pheatmap(mat,
                     display_numbers = cell_labels, number_fontsize = 7,
                     color  = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                     breaks = seq(-col_lim, col_lim, length.out = 101),
                     cluster_cols = FALSE, cluster_rows = nrow(mat) > 1,
                     border_color = "white",
                     main = paste0("PT Subtypes — ", gs_name, "\nPIO NES (* = p<0.05)"),
                     fontsize_row = 8, fontsize_col = 10)
  dev.off()
  cat("Saved PT subtype NES heatmap:", gs_name, "\n")
}

cat("Section 3 complete.\n")

# ============================================================================
# SECTION 4 — CROSS-COMPARISON WITH SENSITIVITY SCRIPT RESULTS
# ============================================================================

cat("\n=== SECTION 4: Cross-Comparison ===\n")

# Helper: load PIO sensitivity results (from existing script outputs)
load_pio_sensitivity <- function(dir.sensitivity, ct_clean, model_label) {
  fp <- paste0(dir.sensitivity, 'NEBULA/NEBULA_', ct_clean, '_', model_label, '.csv')
  if (!file.exists(fp)) return(NULL)
  df <- read.csv(fp)
  df %>%
    transmute(
      gene   = summary.gene,
      logFC  = summary.logFC_group_labelsYes,
      pvalue = summary.p_group_labelsYes,
      fdr    = p.adjust(pvalue, method = "BH")
    )
}

# ---- 4a. PIO vs GLP-1 DEG overlap (base models, each cell type) -------------

cat("\n--- PIO vs GLP-1 DEG overlap ---\n")

overlap_rows <- list()

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  pio_df  <- load_pio_sensitivity(dir.sensitivity, ct_clean, "base")
  glp1_df <- load_glp1_effect(dir.root, ct_clean, "base")
  if (is.null(pio_df) | is.null(glp1_df)) next
  
  merged <- inner_join(
    pio_df  %>% dplyr::select(gene, logFC_pio = logFC, pval_pio = pvalue),
    glp1_df %>% dplyr::select(gene, logFC_glp1 = logFC, pval_glp1 = pvalue),
    by = "gene"
  ) %>%
    mutate(
      sig_pio  = pval_pio  < P_THRESH,
      sig_glp1 = pval_glp1 < P_THRESH,
      status   = case_when(
        sig_pio & sig_glp1  ~ "Sig in both",
        sig_pio & !sig_glp1 ~ "PIO only",
        !sig_pio & sig_glp1 ~ "GLP-1 only",
        TRUE                ~ "NS in both"
      ),
      concordant = ifelse(sig_pio & sig_glp1, sign(logFC_pio) == sign(logFC_glp1), NA)
    )
  
  r_val <- round(cor(merged$logFC_pio, merged$logFC_glp1, use = "complete.obs"), 3)
  
  overlap_rows[[ct]] <- data.frame(
    celltype       = ct,
    n_genes_tested = nrow(merged),
    n_sig_pio      = sum(merged$sig_pio),
    n_sig_glp1     = sum(merged$sig_glp1),
    n_sig_both     = sum(merged$status == "Sig in both"),
    n_pio_only     = sum(merged$status == "PIO only"),
    n_glp1_only    = sum(merged$status == "GLP-1 only"),
    n_concordant   = sum(merged$concordant, na.rm = TRUE),
    n_discordant   = sum(!merged$concordant, na.rm = TRUE),
    r_logFC        = r_val
  )
  
  # Scatter: PIO logFC vs GLP-1 logFC
  delta   <- abs(merged$logFC_pio - merged$logFC_glp1)
  label_v <- ifelse(rank(-delta) <= 8 & merged$status != "NS in both", merged$gene, NA)
  
  p <- ggplot(merged, aes(x = logFC_pio, y = logFC_glp1, color = status, label = label_v)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    ggrepel::geom_text_repel(size = 2.5, color = "black", max.overlaps = 15) +
    scale_color_manual(values = c("Sig in both" = "#d73027", "PIO only" = "#E69F00",
                                  "GLP-1 only" = "#56B4E9", "NS in both" = "grey80")) +
    theme_bw() +
    labs(x = "PIO LogFC (base model)", y = "GLP-1 LogFC (base model)", color = NULL,
         title = paste0(ct, ": PIO vs GLP-1 gene-level effect comparison"),
         subtitle = paste0("r = ", r_val)) +
    theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")
  
  ggsave(paste0(dir.root, 'Comparison/Scatter_PIOvsGLP1_', ct_clean, '.pdf'), p, width = 8, height = 7)
  ggsave(paste0(dir.root, 'Comparison/Scatter_PIOvsGLP1_', ct_clean, '.png'), p, width = 8, height = 7, dpi = 200)
  cat("Saved PIO vs GLP-1 scatter for", ct, "\n")
}

overlap_summary <- bind_rows(overlap_rows)
write.csv(overlap_summary, paste0(dir.root, 'Comparison/PIO_vs_GLP1_DEG_Overlap.csv'), row.names = FALSE)
cat("\nDEG overlap summary:\n"); print(overlap_summary)

# ---- 4b. GSEA NES: PIO vs GLP-1 concordance heatmap per cell type -----------

cat("\n--- PIO vs GLP-1 NES concordance heatmaps ---\n")

# Load sensitivity GSEA (must already exist)
load_sensitivity_gsea <- function(dir.sensitivity, ct, model, gs_name) {
  fp <- paste0(dir.sensitivity, 'GSEA/Raw/all_gsea_sensitivity.csv')
  if (!file.exists(fp)) return(NULL)
  df <- data.table::fread(fp)
  df %>% filter(celltype == ct, model == model, geneset_type == gs_name) %>%
    as.data.frame()
}

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  for (gs_name in names(geneset_types)) {
    pio_gsea  <- load_sensitivity_gsea(dir.sensitivity, ct, "base", gs_name)
    glp1_key  <- paste(ct, "base", gs_name, sep = "__")
    glp1_gsea <- all_glp1_gsea[[glp1_key]]
    if (is.null(pio_gsea) | is.null(glp1_gsea)) next
    
    # Interaction gsea from existing script
    int_key   <- paste0(ct, "__", gs_name)
    int_gsea_fp <- paste0(dir.interaction, 'GSEA/Raw/all_interaction_gsea_results.csv')
    int_gsea  <- if (file.exists(int_gsea_fp)) {
      data.table::fread(int_gsea_fp) %>%
        filter(celltype == ct, geneset_type == gs_name) %>% as.data.frame()
    } else NULL
    
    sig_paths <- unique(c(
      pio_gsea  %>% filter(pval < P_THRESH) %>% pull(pathway),
      glp1_gsea %>% filter(pval < P_THRESH) %>% pull(pathway),
      if (!is.null(int_gsea)) int_gsea %>% filter(pval < P_THRESH) %>% pull(pathway)
    ))
    if (length(sig_paths) == 0) next
    if (length(sig_paths) > 50) {
      # Trim by mean NES across the two primary comparisons
      top_idx <- union(
        pio_gsea  %>% filter(pathway %in% sig_paths) %>% arrange(pval) %>% slice_head(n = 25) %>% pull(pathway),
        glp1_gsea %>% filter(pathway %in% sig_paths) %>% arrange(pval) %>% slice_head(n = 25) %>% pull(pathway)
      )
      sig_paths <- top_idx
    }
    
    build_col <- function(df, label, paths) {
      df %>% filter(pathway %in% paths) %>%
        dplyr::select(pathway, pathway_clean, NES, pval) %>%
        rename_with(~ paste0(c("NES", "pval"), "_", label), c("NES", "pval"))
    }
    
    wide <- build_col(pio_gsea, "PIO", sig_paths)
    wide <- full_join(wide, build_col(glp1_gsea, "GLP1", sig_paths), by = "pathway")
    
    if (!is.null(int_gsea) && nrow(int_gsea) > 0) {
      # Reconstruct pathway_clean if missing
      if (!"pathway_clean" %in% names(int_gsea))
        int_gsea$pathway_clean <- gsub("HALLMARK_|GOBP_", "", int_gsea$pathway) %>%
          gsub("_", " ", .) %>% tools::toTitleCase()
      wide <- full_join(wide, build_col(int_gsea, "Interaction", sig_paths), by = "pathway")
    }
    
    # Build NES matrix
    nes_cols <- grep("^NES_", names(wide), value = TRUE)
    mat <- wide %>%
      filter(pathway %in% sig_paths) %>%
      column_to_rownames("pathway_clean") %>%
      dplyr::select(all_of(nes_cols)) %>% as.matrix()
    colnames(mat) <- str_remove(colnames(mat), "^NES_")
    
    pval_cols <- grep("^pval_", names(wide), value = TRUE)
    pval_mat <- wide %>%
      filter(pathway %in% sig_paths) %>%
      column_to_rownames("pathway_clean") %>%
      dplyr::select(all_of(pval_cols)) %>% as.matrix()
    
    cell_labels <- matrix(
      ifelse(pval_mat < P_THRESH, paste0(round(mat, 2), "*"), round(mat, 2)),
      nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)
    )
    col_lim <- max(abs(mat), na.rm = TRUE)
    
    png(paste0(dir.root, 'Comparison/NES_Compare_PIO_GLP1_Interaction_', ct_clean, '_', gs_name, '.png'),
        width = 1000, height = max(400, nrow(mat) * 22 + 200), res = 120)
    pheatmap::pheatmap(mat,
                       display_numbers = cell_labels, number_fontsize = 7,
                       color  = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                       breaks = seq(-col_lim, col_lim, length.out = 101),
                       cluster_cols = FALSE, cluster_rows = nrow(mat) > 1,
                       border_color = "white",
                       main = paste0(ct, " — ", gs_name,
                                     "\nNES: PIO vs GLP-1 vs Interaction (* = p<0.05)"),
                       fontsize_row = 8, fontsize_col = 10)
    dev.off()
    cat("  Saved concordance heatmap:", ct, gs_name, "\n")
  }
}

# ---- 4c. Robust shared DEGs: sig in PIO base AND GLP-1 base, same direction -

cat("\n--- Shared robust DEGs (PIO + GLP-1, concordant direction) ---\n")

shared_degs <- list()

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  pio_df   <- load_pio_sensitivity(dir.sensitivity, ct_clean, "base")
  glp1_df  <- load_glp1_effect(dir.root, ct_clean, "base")
  if (is.null(pio_df) | is.null(glp1_df)) next
  
  shared_degs[[ct]] <- inner_join(
    pio_df  %>% filter(pvalue < P_THRESH) %>% dplyr::select(gene, logFC_pio  = logFC, pval_pio  = pvalue, fdr_pio  = fdr),
    glp1_df %>% filter(pvalue < P_THRESH) %>% dplyr::select(gene, logFC_glp1 = logFC, pval_glp1 = pvalue, fdr_glp1 = fdr),
    by = "gene"
  ) %>%
    filter(sign(logFC_pio) == sign(logFC_glp1)) %>%
    mutate(direction = ifelse(logFC_pio > 0, "Up in both", "Down in both"),
           celltype  = ct) %>%
    arrange(pval_pio)
}

shared_combined <- bind_rows(shared_degs)
write.csv(shared_combined, paste0(dir.root, 'Comparison/SharedConcordantDEGs_PIO_GLP1.csv'), row.names = FALSE)
cat("Concordant shared DEGs:\n")
print(shared_combined %>% count(celltype, direction))

# ---- 4d. Summary count table: PIO / GLP-1 / Interaction sig genes -----------

cat("\n--- Summary: sig gene counts across all three analyses ---\n")

sig_count_rows <- list()

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  pio_df   <- load_pio_sensitivity(dir.sensitivity, ct_clean, "base")
  glp1_df  <- load_glp1_effect(dir.root, ct_clean, "base")
  
  # Load interaction results from prior script
  int_fp <- paste0(dir.interaction, 'NEBULA/Interaction_NEBULA_', ct_clean, '.csv')
  int_df <- if (file.exists(int_fp)) {
    df <- read.csv(int_fp)
    logfc_c <- grep("logFC.*group_labels.*glp1ra|logFC.*glp1ra.*group_labels",
                    names(df), value = TRUE, ignore.case = TRUE)[1]
    p_c     <- grep("summary\\.p_.*group_labels.*glp1ra|summary\\.p_.*glp1ra.*group_labels",
                    names(df), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(logfc_c) && !is.na(p_c)) {
      df %>% transmute(gene = summary.gene, logFC = .data[[logfc_c]],
                       pvalue = .data[[p_c]], fdr = p.adjust(pvalue, "BH"))
    } else NULL
  } else NULL
  
  sig_count_rows[[ct]] <- data.frame(
    celltype      = ct,
    n_sig_pio     = if (!is.null(pio_df))  sum(pio_df$pvalue  < P_THRESH, na.rm = TRUE) else NA,
    n_sig_glp1    = if (!is.null(glp1_df)) sum(glp1_df$pvalue < P_THRESH, na.rm = TRUE) else NA,
    n_sig_int     = if (!is.null(int_df))  sum(int_df$pvalue  < P_THRESH, na.rm = TRUE) else NA,
    n_shared_conc = if (ct %in% names(shared_degs)) nrow(shared_degs[[ct]]) else 0
  )
}

sig_count_df <- bind_rows(sig_count_rows)
write.csv(sig_count_df, paste0(dir.root, 'Comparison/SigGeneCounts_PIO_GLP1_Interaction.csv'), row.names = FALSE)
cat("\nSignificant gene count summary:\n"); print(sig_count_df)

# ============================================================================
# DONE
# ============================================================================

cat("\n=== GLP-1 Extended Analysis Complete ===\n")
cat("Root output directory:", dir.root, "\n\n")
cat("Key outputs:\n")
cat("  Demographics/                        — four-group patient/cell counts, UMAP\n")
cat("  GLP1_Main/NEBULA/                    — GLP-1 NEBULA results (all models)\n")
cat("  GLP1_Main/VolcanoPlots/              — GLP-1 volcano grids\n")
cat("  GLP1_Main/GSEA/                      — GLP-1 GSEA dot plots + raw results\n")
cat("  PT_Subtypes/NEBULA/                  — PT subtype NEBULA (PIO + GLP-1, adj)\n")
cat("  PT_Subtypes/VolcanoPlots/            — PT subtype volcano grids\n")
cat("  PT_Subtypes/GSEA/                    — PT subtype GSEA + NES heatmap\n")
cat("  Comparison/Scatter_PIOvsGLP1_*       — gene-level PIO vs GLP-1 logFC scatter\n")
cat("  Comparison/NES_Compare_*             — PIO / GLP-1 / Interaction NES heatmap\n")
cat("  Comparison/SharedConcordantDEGs_*    — genes sig & concordant in both drugs\n")
cat("  Comparison/SigGeneCounts_*.csv       — counts summary across all three analyses\n")