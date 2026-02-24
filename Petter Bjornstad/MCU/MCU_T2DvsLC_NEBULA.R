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

# ── 0. GENE LISTS ─────────────────────────────────────────────────────────────
mcu_complex_genes <- c("MCU", "MCUR1", "MICU1", "MICU2", "SMDT1",
                       "VDAC1", "VDAC2", "VDAC3")

calcium_signaling_genes <- c("ITPR1", "ITPR2", "ITPR3",
                             "RYR1", "RYR2",
                             "SLC8A1", "SLC8A3",
                             "ATP2A1", "ATP2A2", "ATP2A3",
                             "CALM1", "CALM2", "CALM3",
                             "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
                             "SLC25A23", "SLC25A24", "SLC25A25")

tca_cycle_genes <- c("CS", "ACO1", "ACO2",
                     "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G",
                     "OGDH", "OGDHL", "DLST", "DLD",
                     "SUCLA2", "SUCLG1", "SUCLG2",
                     "SDHA", "SDHB", "SDHC", "SDHD",
                     "FH", "MDH1", "MDH2")

pyruvate_genes <- c("PDHA1", "PDHA2", "PDHB",
                    "PDK1", "PDK2", "PDK3", "PDK4",
                    "PDP1", "PDP2",
                    "PC",
                    "MPC1", "MPC2",
                    "LDHA", "LDHB",
                    "ME1", "ME2", "ME3")

all_genes_of_interest <- unique(c(mcu_complex_genes, calcium_signaling_genes,
                                  tca_cycle_genes, pyruvate_genes))

# ── 1. CELL TYPE SUBSETTING HELPER ───────────────────────────────────────────
# Mirrors the established logic used across your sex / pioglitazone / SGLT2i analyses:
#   PT, TAL, EC, POD        → celltype2 column
#   DCTall                  → DCT_celltype column  (value == 'DCT')
#   intercalated            → KPMP_celltype column (IC-A, IC-B, tPC-IC, aIC)
#   All                     → no subsetting

subset_by_celltype <- function(so_obj, celltype) {
  if (celltype == "All") {
    return(so_obj)
  } else if (celltype %in% c("PT", "TAL", "EC", "POD")) {
    return(subset(so_obj, celltype2 == celltype))
  } else if (celltype == "DCTall") {
    return(subset(so_obj, DCT_celltype == "DCT"))
  } else if (celltype == "intercalated") {
    return(subset(so_obj, KPMP_celltype %in% c("IC-A", "IC-B", "tPC-IC", "aIC")))
  } else {
    stop("Unknown celltype: ", celltype,
         ". Add a new branch to subset_by_celltype() if needed.")
  }
}

# ── 2. SEURAT OBJECT & POOLED OFFSET ─────────────────────────────────────────
# Offset MUST be computed on the full object BEFORE any cell-type subsetting.

so_full <- readRDS("")   # ← adjust path

counts_full <- round(GetAssayData(so_full, layer = "counts"))
sce_full    <- SingleCellExperiment(assays = list(counts = counts_full))
sce_full    <- computeSumFactors(sce_full)
so_full@meta.data$pooled_offset <- sizeFactors(sce_full)
rm(sce_full, counts_full); gc()

# ── 3. MAIN ANALYSIS FUNCTION (sequential, no parallelisation) ───────────────
run_nebula_celltype <- function(so_obj,
                                celltype,
                                genes_to_test     = all_genes_of_interest,
                                predictor_formula = ~diabetes,   # ← adjust
                                id_col            = "kit_id",
                                dir_out           = "results/NEBULA_MCU_TCA/") {
  
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  message("\n========== Cell type: ", celltype, " ==========")
  
  # --- subset using the correct column for this cell type
  so_sub <- subset_by_celltype(so_obj, celltype)
  
  # --- keep only genes present in this subset
  genes_available <- intersect(genes_to_test,
                               rownames(GetAssayData(so_sub, layer = "counts")))
  message("Testing ", length(genes_available), " / ", length(genes_to_test),
          " requested genes present in data.")
  
  if (length(genes_available) == 0) {
    message("No genes available for ", celltype, " – skipping.")
    return(NULL)
  }
  
  counts_path <- round(GetAssayData(so_sub, layer = "counts"))
  
  results_list <- vector("list", length(genes_available))
  names(results_list) <- genes_available
  
  for (g in genes_available) {
    message("  Gene: ", g)
    
    result_entry <- tryCatch({
      
      count_gene <- counts_path[g, , drop = FALSE]
      meta_gene  <- so_sub@meta.data
      
      # Drop cells with NA in any predictor variable
      pred_vars <- all.vars(predictor_formula)
      keep      <- complete.cases(meta_gene[, pred_vars, drop = FALSE])
      if (sum(keep) == 0) stop("No complete cases after NA removal.")
      meta_gene  <- meta_gene[keep, ]
      count_gene <- count_gene[, keep, drop = FALSE]
      
      pred_gene <- model.matrix(predictor_formula, data = meta_gene)
      library   <- meta_gene$pooled_offset
      
      data_g <- group_cell(count  = count_gene,
                           id     = meta_gene[[id_col]],
                           pred   = pred_gene,
                           offset = library)
      
      if (is.null(data_g)) {
        data_g <- list(count   = count_gene,
                       id      = meta_gene[[id_col]],
                       pred    = pred_gene,
                       library = library)
      }
      
      result <- nebula(count      = data_g$count,
                       id         = data_g$id,
                       pred       = data_g$pred,
                       ncore      = 1,
                       reml       = TRUE,
                       model      = "NBLMM",
                       output_re  = TRUE,
                       covariance = TRUE,
                       offset     = data_g$library)
      
      list(gene = g, result = result)
      
    }, error = function(e) {
      message("    ERROR on gene ", g, ": ", e$message)
      NULL
    })
    
    results_list[[g]] <- result_entry
  }
  
  # --- tidy results into a data frame
  results_df <- bind_rows(lapply(results_list, function(x) {
    if (is.null(x)) return(NULL)
    df      <- as.data.frame(x$result$summary)
    df$gene <- x$gene
    df
  }))
  
  if (nrow(results_df) == 0) {
    message("No results returned for ", celltype)
    return(NULL)
  }
  
  # BH-adjust every p-value column
  p_cols <- grep("^p_", colnames(results_df), value = TRUE)
  for (pc in p_cols) {
    results_df[[paste0("padj_", sub("^p_", "", pc))]] <-
      p.adjust(results_df[[pc]], method = "BH")
  }
  
  results_df$celltype <- celltype
  
  outfile <- file.path(dir_out,
                       paste0("NEBULA_", gsub("[/ ]", "_", celltype), ".csv"))
  write.csv(results_df, outfile, row.names = FALSE)
  message("Saved: ", outfile)
  
  return(results_df)
}

# ── 4. RUN ACROSS ALL CELL TYPES ─────────────────────────────────────────────
# Adjust predictor_formula as needed:
#   ~diabetes                   → main effect only
#   ~diabetes + sex + age       → covariate-adjusted
#   ~diabetes * sex             → interaction model
#   ~medication_class           → medication comparison

cell_types_to_run <- c("All", "PT", "TAL", "EC", "POD", "DCTall", "intercalated")

all_results <- lapply(cell_types_to_run, function(ct) {
  run_nebula_celltype(
    so_obj            = so_full,
    celltype          = ct,
    genes_to_test     = all_genes_of_interest,
    predictor_formula = ~diabetes,    # ← YOUR predictor here
    id_col            = "kit_id",
    dir_out           = "results/NEBULA_MCU_TCA/"
  )
})

all_results_df <- bind_rows(all_results)
write.csv(all_results_df,
          "results/NEBULA_MCU_TCA/NEBULA_all_celltypes_combined.csv",
          row.names = FALSE)

# ── 5. VOLCANO PLOT ───────────────────────────────────────────────────────────
# After running, check colnames(all_results_df) to confirm exact column names.
# NEBULA typically outputs e.g. "logFC_diabetesTRUE" and "p_diabetesTRUE".

plot_volcano_celltype <- function(df, celltype, predictor_term,
                                  fc_cutoff = 0.25, p_cutoff = 0.05) {
  sub_df  <- df %>% filter(celltype == !!celltype)
  lfc_col  <- grep(paste0("logFC.*", predictor_term), colnames(sub_df), value = TRUE)[1]
  padj_col <- grep(paste0("padj.*", predictor_term), colnames(sub_df), value = TRUE)[1]
  
  if (is.na(lfc_col) | is.na(padj_col)) {
    message("Could not find logFC or padj columns. Check colnames(df).")
    return(NULL)
  }
  
  sub_df <- sub_df %>%
    mutate(
      module = case_when(
        gene %in% mcu_complex_genes       ~ "MCU complex",
        gene %in% calcium_signaling_genes ~ "Ca2+ signaling",
        gene %in% tca_cycle_genes         ~ "TCA cycle",
        gene %in% pyruvate_genes          ~ "Pyruvate",
        TRUE                              ~ "Other"
      ),
      neg_log10p = -log10(.data[[padj_col]]),
      logFC      = .data[[lfc_col]],
      label      = ifelse(.data[[padj_col]] < p_cutoff & abs(logFC) > fc_cutoff, gene, NA)
    )
  
  module_colors <- c("MCU complex"    = "#E63946",
                     "Ca2+ signaling" = "#F4A261",
                     "TCA cycle"      = "#2A9D8F",
                     "Pyruvate"       = "#457B9D",
                     "Other"          = "grey70")
  
  ggplot(sub_df, aes(x = logFC, y = neg_log10p, color = module, label = label)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 20, na.rm = TRUE) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(p_cutoff),          linetype = "dashed", color = "grey50") +
    scale_color_manual(values = module_colors) +
    labs(title    = paste0(celltype, " – Diabetes effect"),
         subtitle = paste0("NEBULA-NBLMM | predictor: ", predictor_term),
         x        = "log fold change",
         y        = expression(-log[10](p[adj])),
         color    = "Gene module") +
    theme_bw(base_size = 13)
}

# Example – adjust predictor_term to match your actual column names
plot_volcano_celltype(all_results_df, celltype = "PT", predictor_term = "diabetes")

# ── 6. SUMMARY HEATMAP (genes × cell types) ──────────────────────────────────
# Check and adjust column names as needed after inspecting all_results_df
lfc_col  <- grep("logFC.*diabetes", colnames(all_results_df), value = TRUE)[1]
padj_col <- grep("padj.*diabetes",  colnames(all_results_df), value = TRUE)[1]

heat_df <- all_results_df %>%
  filter(gene %in% all_genes_of_interest) %>%
  dplyr::select(gene, celltype, all_of(c(lfc_col, padj_col))) %>%
  rename(logFC = all_of(lfc_col), padj = all_of(padj_col))

lfc_mat <- heat_df %>%
  pivot_wider(id_cols = gene, names_from = celltype, values_from = logFC) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

sig_mat <- heat_df %>%
  pivot_wider(id_cols = gene, names_from = celltype, values_from = padj) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()
sig_labels <- ifelse(sig_mat < 0.05, "*", "")

row_annot <- data.frame(
  Module = case_when(
    rownames(lfc_mat) %in% mcu_complex_genes       ~ "MCU complex",
    rownames(lfc_mat) %in% calcium_signaling_genes ~ "Ca2+ signaling",
    rownames(lfc_mat) %in% tca_cycle_genes         ~ "TCA cycle",
    rownames(lfc_mat) %in% pyruvate_genes          ~ "Pyruvate",
    TRUE                                            ~ "Other"
  ),
  row.names = rownames(lfc_mat)
)

annot_colors <- list(Module = c("MCU complex"    = "#E63946",
                                "Ca2+ signaling" = "#F4A261",
                                "TCA cycle"      = "#2A9D8F",
                                "Pyruvate"       = "#457B9D"))

pheatmap(lfc_mat,
         display_numbers   = sig_labels,
         fontsize_number   = 10,
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         annotation_row    = row_annot,
         annotation_colors = annot_colors,
         color             = colorRampPalette(c("#457B9D", "white", "#E63946"))(100),
         breaks            = seq(-2, 2, length.out = 101),
         na_col            = "grey90",
         main              = "NEBULA logFC (diabetes) – MCU / Ca2+ / TCA / Pyruvate",
         filename          = "results/NEBULA_MCU_TCA/heatmap_all_celltypes.pdf",
         width = 12, height = 14)

